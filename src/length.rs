use crate::classify::{
    Classification, ClassificationIndex, Classifier, Thresholds, apply_discriminatory_filter,
    build_classification_index, load_classification_index,
};
use crate::{
    IndexKind, Layout, ProcessingStats, StdinTargets, TargetSource, check_index_complexity,
    create_spinner, format_bp, format_bp_per_sec, process_input, resolve_targets, sample_inputs,
    sample_limit_reached_io_error,
};
use anyhow::Result;
use indicatif::ProgressBar;
use paraseq::Record;
use paraseq::parallel::{PairedParallelProcessor, ParallelProcessor};
use parking_lot::Mutex;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

const AMBIGUOUS_LABEL: &str = "ambiguous";
const UNCLASSIFIED_LABEL: &str = "unclassified";
const ALL_LABEL: &str = "all";

/// Per-sample length histogram output. Buckets are indexed as
/// `0..group_names.len()` for groups, then `ambiguous`, then `unclassified`.
#[derive(Debug)]
struct LengthHistogramResult {
    sample_name: String,
    total_seqs_processed: u64,
    total_bp_processed: u64,
    bucket_histograms: Vec<Vec<(usize, usize)>>,
    bucket_seqs: Vec<u64>,
    bucket_bases: Vec<u64>,
}

/// `lenhist` run settings
pub struct LengthHistogramConfig {
    /// Fastx file, directory of groups, `.sk` index, or `-` to disable group filtering
    pub targets_path: PathBuf,
    pub individual: bool,
    pub sample_paths: Vec<Vec<PathBuf>>,
    pub sample_names: Vec<String>,
    pub layout: Layout,
    pub kmer_length: u8,
    pub smer_length: u8,
    pub complexity: f32,
    pub abs_threshold: u64,
    pub rel_threshold: f64,
    pub discriminatory: bool,
    pub threads: usize,
    pub output_path: Option<PathBuf>,
    pub quiet: bool,
    pub limit_bp: Option<u64>,
    /// Whether targets_path disables filtering, placing everything in one bucket
    pub no_filter: bool,
}

impl LengthHistogramConfig {
    pub fn execute(&self) -> Result<()> {
        run_lenhist(self)
    }
}

/// Per-bucket local state (cleared on each batch flush)
#[derive(Clone, Default)]
struct BucketState {
    histogram: HashMap<usize, usize>,
    seqs: u64,
    bases: u64,
}

/// Worker that bins lengths per classification bucket
#[derive(Clone)]
struct LengthHistogramProcessor {
    classifier: Classifier,
    no_filter: bool,

    local_stats: ProcessingStats,
    local_buckets: Vec<BucketState>,

    // Global state
    global_stats: Arc<Mutex<ProcessingStats>>,
    global_buckets: Arc<Vec<Mutex<BucketState>>>,
    spinner: Option<Arc<Mutex<ProgressBar>>>,
    start_time: Instant,
    limit_bp: Option<u64>,
}

impl LengthHistogramProcessor {
    #[allow(clippy::too_many_arguments)]
    fn new(
        classifier: Classifier,
        no_filter: bool,
        global_buckets: Arc<Vec<Mutex<BucketState>>>,
        global_stats: Arc<Mutex<ProcessingStats>>,
        spinner: Option<Arc<Mutex<ProgressBar>>>,
        start_time: Instant,
        limit_bp: Option<u64>,
    ) -> Self {
        let bucket_count = classifier.num_groups + 2;
        let local_buckets = (0..bucket_count).map(|_| BucketState::default()).collect();

        Self {
            classifier,
            no_filter,
            local_stats: ProcessingStats::default(),
            local_buckets,
            global_stats,
            global_buckets,
            spinner,
            start_time,
            limit_bp,
        }
    }

    fn update_spinner(&self) {
        if let Some(ref spinner) = self.spinner {
            let stats = self.global_stats.lock();
            let elapsed = self.start_time.elapsed();
            let seqs_per_sec = stats.total_seqs as f64 / elapsed.as_secs_f64();
            let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();

            spinner.lock().set_message(format!(
                "Processing sample: {} seqs ({}). {:.0} seqs/s ({})",
                stats.total_seqs,
                format_bp(stats.total_bp as usize),
                seqs_per_sec,
                format_bp_per_sec(bp_per_sec)
            ));
        }
    }
}

impl LengthHistogramProcessor {
    /// Bin mates separately under their pooled classification
    fn process(&mut self, seqs: &[&[u8]]) -> paraseq::Result<()> {
        if self
            .limit_bp
            .is_some_and(|limit| self.global_stats.lock().total_bp >= limit)
        {
            self.flush();
            return Err(paraseq::Error::Io(sample_limit_reached_io_error()));
        }

        let bucket_idx = if self.no_filter {
            0
        } else {
            let num_groups = self.classifier.num_groups;
            match self.classifier.classify(seqs) {
                Classification::Classified(g) => g,
                Classification::Ambiguous(_) => num_groups,
                Classification::Unclassified => num_groups + 1,
            }
        };

        let bucket = &mut self.local_buckets[bucket_idx];
        for seq in seqs {
            *bucket.histogram.entry(seq.len()).or_insert(0) += 1;
            bucket.seqs += 1;
            bucket.bases += seq.len() as u64;
            self.local_stats.total_seqs += 1;
            self.local_stats.total_bp += seq.len() as u64;
        }
        Ok(())
    }

    fn flush(&mut self) {
        // Merge local buckets into global
        for (i, local) in self.local_buckets.iter_mut().enumerate() {
            if local.seqs == 0 && local.histogram.is_empty() {
                continue;
            }
            let mut global = self.global_buckets[i].lock();
            for (&length, &count) in &local.histogram {
                *global.histogram.entry(length).or_insert(0) += count;
            }
            global.seqs += local.seqs;
            global.bases += local.bases;
            local.histogram.clear();
            local.seqs = 0;
            local.bases = 0;
        }

        // Update global stats
        let mut stats = self.global_stats.lock();
        stats.total_seqs += self.local_stats.total_seqs;
        stats.total_bp += self.local_stats.total_bp;

        // Update spinner every 0.1 Gbp
        let current_progress = stats.total_bp / 100_000_000;
        if current_progress > stats.last_reported {
            drop(stats);
            self.update_spinner();
            self.global_stats.lock().last_reported = current_progress;
        }

        self.local_stats = ProcessingStats::default();
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for LengthHistogramProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        self.process(&[&record.seq()])
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        self.flush();
        Ok(())
    }
}

impl<Rf: Record> PairedParallelProcessor<Rf> for LengthHistogramProcessor {
    fn process_record_pair(&mut self, record1: Rf, record2: Rf) -> paraseq::Result<()> {
        self.process(&[&record1.seq(), &record2.seq()])
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        self.flush();
        Ok(())
    }
}

fn process_seqs_input(
    input: &[PathBuf],
    layout: Layout,
    classifier: &Classifier,
    threads: usize,
    quiet: bool,
    no_filter: bool,
    limit_bp: Option<u64>,
) -> Result<(Vec<BucketState>, u64, u64)> {
    let spinner = create_spinner(quiet)?;
    if let Some(ref pb) = spinner {
        pb.lock().set_message("Processing sample: 0 seqs (0bp)");
    }

    let num_groups = classifier.num_groups;
    let bucket_count = num_groups + 2;
    let global_buckets: Arc<Vec<Mutex<BucketState>>> = Arc::new(
        (0..bucket_count)
            .map(|_| Mutex::new(BucketState::default()))
            .collect(),
    );
    let global_stats = Arc::new(Mutex::new(ProcessingStats::default()));

    let start_time = Instant::now();
    let mut processor = LengthHistogramProcessor::new(
        classifier.clone(),
        no_filter,
        Arc::clone(&global_buckets),
        Arc::clone(&global_stats),
        spinner.clone(),
        start_time,
        limit_bp,
    );

    process_input(input, layout, &mut processor, threads)?;

    if let Some(ref pb) = spinner {
        pb.lock().finish_with_message("");
    }

    let stats = global_stats.lock().clone();
    let buckets: Vec<BucketState> = global_buckets
        .iter()
        .map(|m| std::mem::take(&mut *m.lock()))
        .collect();

    if !quiet {
        let elapsed = start_time.elapsed();
        let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();
        let classified: u64 = buckets[..num_groups].iter().map(|b| b.seqs).sum();
        let ambiguous = buckets[num_groups].seqs;
        let unclassified = buckets[num_groups + 1].seqs;
        eprintln!(
            "Sample: {} records ({}), {} classified, {} ambiguous, {} unclassified ({})",
            stats.total_seqs,
            format_bp(stats.total_bp as usize),
            classified,
            ambiguous,
            unclassified,
            format_bp_per_sec(bp_per_sec)
        );
    }

    Ok((buckets, stats.total_seqs, stats.total_bp))
}

/// Process a single sample's sequences and calculate per-bucket length histograms
fn process_single_sample(
    sample_paths: &[PathBuf],
    sample_name: &str,
    classifier: &Classifier,
    config: &LengthHistogramConfig,
) -> Result<LengthHistogramResult> {
    // Silence per-sample progress for >1 sample
    let quiet_sample = config.quiet || config.sample_paths.len() > 1;

    let bucket_count = classifier.num_groups + 2;
    let mut combined_buckets: Vec<BucketState> = vec![BucketState::default(); bucket_count];
    let mut total_seqs = 0u64;
    let mut total_bp = 0u64;

    for input in sample_inputs(sample_paths, config.layout) {
        let (file_buckets, file_seqs, file_bp) = process_seqs_input(
            input,
            config.layout,
            classifier,
            config.threads,
            quiet_sample,
            config.no_filter,
            config.limit_bp.map(|limit| limit.saturating_sub(total_bp)),
        )?;

        for (i, b) in file_buckets.into_iter().enumerate() {
            for (length, count) in b.histogram {
                *combined_buckets[i].histogram.entry(length).or_insert(0) += count;
            }
            combined_buckets[i].seqs += b.seqs;
            combined_buckets[i].bases += b.bases;
        }

        total_seqs += file_seqs;
        total_bp += file_bp;

        if let Some(limit) = config.limit_bp
            && total_bp >= limit
        {
            break;
        }
    }

    let mut bucket_histograms = Vec::with_capacity(bucket_count);
    let mut bucket_seqs = Vec::with_capacity(bucket_count);
    let mut bucket_bases = Vec::with_capacity(bucket_count);
    for b in combined_buckets {
        let mut hist: Vec<(usize, usize)> = b.histogram.into_iter().collect();
        hist.sort_by_key(|(length, _)| *length);
        bucket_histograms.push(hist);
        bucket_seqs.push(b.seqs);
        bucket_bases.push(b.bases);
    }

    Ok(LengthHistogramResult {
        sample_name: sample_name.to_string(),
        total_seqs_processed: total_seqs,
        total_bp_processed: total_bp,
        bucket_histograms,
        bucket_seqs,
        bucket_bases,
    })
}

pub fn run_lenhist(config: &LengthHistogramConfig) -> Result<()> {
    let version = env!("CARGO_PKG_VERSION").to_string();

    let mut options = format!(
        "k={}, s={}, threads={}",
        config.kmer_length, config.smer_length, config.threads
    );

    if !config.no_filter {
        options.push_str(&format!(
            ", abs_threshold={}, rel_threshold={}",
            config.abs_threshold, config.rel_threshold
        ));
        if config.discriminatory {
            options.push_str(", discriminatory");
        }
    }
    options.push_str(config.layout.option_label());

    if config.sample_paths.len() > 1 {
        options.push_str(&format!(", samples={}", config.sample_paths.len()));
    }

    if config.complexity > 0.0 {
        options.push_str(&format!(", complexity={}", config.complexity));
    }

    if let Some(limit) = config.limit_bp {
        options.push_str(&format!(", limit={}", format_bp(limit as usize)));
    }

    eprintln!("Skope v{}; mode: lenhist; options: {}", version, options);

    // Load or build the classification index
    let (index, group_names, kmer_length, smer_length) = if config.no_filter {
        if !config.quiet {
            eprintln!("Groups: none (target filtering disabled, single \"all\" bucket)");
        }
        // Degenerate 1-group "all" setup; the index is never queried in this mode
        let empty_index = if config.kmer_length <= 32 {
            ClassificationIndex::U64(HashMap::with_hasher(crate::FixedRapidHasher))
        } else {
            ClassificationIndex::U128(HashMap::with_hasher(crate::FixedRapidHasher))
        };
        (
            empty_index,
            vec![ALL_LABEL.to_string()],
            config.kmer_length,
            config.smer_length,
        )
    } else {
        let src = config.targets_path.to_string_lossy().to_string();
        match resolve_targets(
            &config.targets_path,
            IndexKind::Classify,
            StdinTargets::Reject,
        )? {
            TargetSource::Index(path) => {
                let (index, group_names, k, s, index_complexity) =
                    load_classification_index(&path)?;
                check_index_complexity(config.complexity, index_complexity)?;
                if !config.quiet {
                    eprintln!(
                        "Index: {} k-mers, {} groups, k={}, s={}",
                        index.len(),
                        group_names.len(),
                        k,
                        s,
                    );
                    for (i, name) in group_names.iter().enumerate() {
                        eprintln!("  [{}] {}", i, name);
                    }
                }
                (index, group_names, k, s)
            }
            source => {
                let (index, group_names) = build_classification_index(
                    source.groups(),
                    &src,
                    source.splits_records(config.individual),
                    config.kmer_length,
                    config.smer_length,
                    config.complexity,
                    config.threads,
                    config.quiet,
                )?;
                (index, group_names, config.kmer_length, config.smer_length)
            }
        }
    };

    let mut index = index;
    if !config.no_filter && config.discriminatory {
        let removed = apply_discriminatory_filter(&mut index);
        if !config.quiet {
            eprintln!(
                "Discriminatory mode: removed {} shared k-mers, {} unique k-mers remain",
                removed,
                index.len()
            );
        }
    }

    let classifier = Classifier::new(
        Arc::new(index),
        group_names.len(),
        kmer_length,
        smer_length,
        Thresholds {
            abs: config.abs_threshold,
            rel: config.rel_threshold,
        },
        false,
    );

    // Process each sample in parallel
    use rayon::prelude::*;

    let is_multisample = config.sample_paths.len() > 1;
    let completed = if is_multisample && !config.quiet {
        eprint!(
            "\x1B[2K\rSamples: processed 0 of {}…",
            config.sample_paths.len()
        );
        Some(Arc::new(Mutex::new(0usize)))
    } else {
        if !config.quiet {
            eprint!("\x1B[2K\r");
        }
        None
    };

    let sample_results: Vec<LengthHistogramResult> = config
        .sample_paths
        .par_iter()
        .zip(&config.sample_names)
        .map(|(sample_paths, sample_name)| {
            let result = process_single_sample(sample_paths, sample_name, &classifier, config);

            if let Some(ref counter) = completed {
                let mut count = counter.lock();
                *count += 1;
                eprint!(
                    "\rSamples: processed {} of {}…",
                    *count,
                    config.sample_paths.len()
                );
            }

            result
        })
        .collect::<Result<Vec<_>>>()?;

    if is_multisample && !config.quiet {
        eprintln!();
    }

    // Output TSV
    let writer: Box<dyn Write> = if let Some(path) = &config.output_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };

    let mut csv_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);

    csv_writer.write_record([
        "sample",
        "group",
        "length",
        "count",
        "total_seqs_processed",
        "total_bp_processed",
        "group_seqs",
        "group_bases",
    ])?;

    for sample in &sample_results {
        let n = group_names.len();
        for (bucket_idx, hist) in sample.bucket_histograms.iter().enumerate() {
            if hist.is_empty() {
                continue;
            }
            let group_label: &str = if bucket_idx < n {
                &group_names[bucket_idx]
            } else if bucket_idx == n {
                AMBIGUOUS_LABEL
            } else {
                UNCLASSIFIED_LABEL
            };
            let group_seqs = sample.bucket_seqs[bucket_idx];
            let group_bases = sample.bucket_bases[bucket_idx];
            for (length, count) in hist {
                csv_writer.write_record([
                    &sample.sample_name,
                    group_label,
                    &length.to_string(),
                    &count.to_string(),
                    &sample.total_seqs_processed.to_string(),
                    &sample.total_bp_processed.to_string(),
                    &group_seqs.to_string(),
                    &group_bases.to_string(),
                ])?;
            }
        }
    }

    csv_writer.flush()?;

    Ok(())
}

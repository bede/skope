use crate::classify::{
    Classification, ClassificationIndex, apply_discriminatory_filter, build_index_in_memory,
    classify_seq_kmers, load_index,
};
use crate::query::TimingStats;
use crate::syncmers::{Buffers, KmerHasher};
use crate::{
    ProcessingStats, create_spinner, format_bp, format_bp_per_sec, handle_process_result,
    reader_with_inferred_batch_size, sample_limit_reached_io_error,
};
use anyhow::Result;
use indicatif::ProgressBar;
use paraseq::Record;
use paraseq::parallel::{ParallelProcessor, ParallelReader};
use parking_lot::Mutex;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

const AMBIGUOUS_LABEL: &str = "ambiguous";
const UNCLASSIFIED_LABEL: &str = "unclassified";
const ALL_LABEL: &str = "all";

/// Per-sample length histogram output. Buckets are indexed as
/// `0..group_names.len()` for groups, then `ambiguous`, then `unclassified`.
#[derive(Debug, Clone)]
pub struct LengthHistogramResult {
    pub sample_name: String,
    pub seq_files: Vec<String>,
    pub total_seqs_processed: u64,
    pub total_bp_processed: u64,
    pub group_names: Vec<String>,
    pub bucket_histograms: Vec<Vec<(usize, usize)>>,
    pub bucket_seqs: Vec<u64>,
    pub bucket_bases: Vec<u64>,
}

/// Combined `lenhist` output
#[derive(Debug, Clone)]
pub struct LengthHistogramReport {
    pub version: String,
    pub index_source: String,
    pub parameters: LengthHistogramParameters,
    pub samples: Vec<LengthHistogramResult>,
    pub timing: TimingStats,
}

/// `lenhist` parameters
#[derive(Debug, Clone)]
pub struct LengthHistogramParameters {
    pub kmer_length: u8,
    pub smer_length: u8,
    pub threads: usize,
    pub min_hits: u64,
    pub min_fraction: f64,
    pub discriminatory: bool,
}

/// `lenhist` run settings
pub struct LengthHistogramConfig {
    /// Path to .sk index file, directory of groups, or "-" to disable group filtering
    pub index_path: PathBuf,
    pub sample_paths: Vec<Vec<PathBuf>>,
    pub sample_names: Vec<String>,
    pub kmer_length: u8,
    pub smer_length: u8,
    pub min_hits: u64,
    pub min_fraction: f64,
    pub discriminatory: bool,
    pub threads: usize,
    pub output_path: Option<PathBuf>,
    pub quiet: bool,
    pub limit_bp: Option<u64>,
    /// True when index_path == "-": no filtering, all reads go to a single "all" bucket
    pub include_all_seqs: bool,
}

impl LengthHistogramConfig {
    pub fn execute(&self) -> Result<()> {
        run_length_histogram_analysis(self)
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
    kmer_length: u8,
    smer_length: u8,
    hasher: KmerHasher,
    index: Arc<ClassificationIndex>,
    num_groups: usize,
    min_hits: u64,
    min_fraction: f64,
    include_all_seqs: bool,

    // Local buffers
    buffers: Buffers,
    hits: [u64; 128],
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
        kmer_length: u8,
        smer_length: u8,
        index: Arc<ClassificationIndex>,
        num_groups: usize,
        min_hits: u64,
        min_fraction: f64,
        include_all_seqs: bool,
        global_buckets: Arc<Vec<Mutex<BucketState>>>,
        global_stats: Arc<Mutex<ProcessingStats>>,
        spinner: Option<Arc<Mutex<ProgressBar>>>,
        start_time: Instant,
        limit_bp: Option<u64>,
    ) -> Self {
        let buffers = if kmer_length <= 32 {
            Buffers::new_u64()
        } else {
            Buffers::new_u128()
        };

        let bucket_count = num_groups + 2;
        let local_buckets = (0..bucket_count).map(|_| BucketState::default()).collect();

        Self {
            kmer_length,
            smer_length,
            hasher: KmerHasher::new(smer_length as usize),
            index,
            num_groups,
            min_hits,
            min_fraction,
            include_all_seqs,
            buffers,
            hits: [0u64; 128],
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

impl<Rf: Record> ParallelProcessor<Rf> for LengthHistogramProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        if let Some(limit) = self.limit_bp {
            let global_bp = self.global_stats.lock().total_bp;
            if global_bp >= limit {
                return Err(paraseq::parallel::ProcessError::IoError(
                    sample_limit_reached_io_error(),
                ));
            }
        }

        let seq = record.seq();
        let seq_length = seq.len();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq_length as u64;

        let bucket_idx = if self.include_all_seqs {
            0
        } else {
            let (_total_kmers, classification) = classify_seq_kmers(
                &seq,
                &self.hasher,
                self.kmer_length,
                self.smer_length,
                &mut self.buffers,
                &mut self.hits,
                self.num_groups,
                &self.index,
                self.min_hits,
                self.min_fraction,
            );
            match classification {
                Classification::Classified(g) => g,
                Classification::Ambiguous(_) => self.num_groups,
                Classification::Unclassified => self.num_groups + 1,
            }
        };

        let bucket = &mut self.local_buckets[bucket_idx];
        *bucket.histogram.entry(seq_length).or_insert(0) += 1;
        bucket.seqs += 1;
        bucket.bases += seq_length as u64;

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
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
        {
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

        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn process_seqs_file(
    seq_path: &Path,
    index: Arc<ClassificationIndex>,
    num_groups: usize,
    kmer_length: u8,
    smer_length: u8,
    min_hits: u64,
    min_fraction: f64,
    threads: usize,
    quiet: bool,
    include_all_seqs: bool,
    limit_bp: Option<u64>,
) -> Result<(Vec<BucketState>, u64, u64)> {
    let in_path = if seq_path.to_string_lossy() == "-" {
        None
    } else {
        Some(seq_path)
    };
    let reader = reader_with_inferred_batch_size(in_path)?;

    let spinner = create_spinner(quiet)?;
    if let Some(ref pb) = spinner {
        pb.lock().set_message("Processing sample: 0 seqs (0bp)");
    }

    let bucket_count = num_groups + 2;
    let global_buckets: Arc<Vec<Mutex<BucketState>>> = Arc::new(
        (0..bucket_count)
            .map(|_| Mutex::new(BucketState::default()))
            .collect(),
    );
    let global_stats = Arc::new(Mutex::new(ProcessingStats::default()));

    let start_time = Instant::now();
    let mut processor = LengthHistogramProcessor::new(
        kmer_length,
        smer_length,
        index,
        num_groups,
        min_hits,
        min_fraction,
        include_all_seqs,
        Arc::clone(&global_buckets),
        Arc::clone(&global_stats),
        spinner.clone(),
        start_time,
        limit_bp,
    );

    let process_result = reader.process_parallel(&mut processor, threads);
    handle_process_result(process_result)?;

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
    index: Arc<ClassificationIndex>,
    group_names: &[String],
    kmer_length: u8,
    smer_length: u8,
    config: &LengthHistogramConfig,
) -> Result<LengthHistogramResult> {
    // Silence per-sample progress for >1 sample
    let quiet_sample = config.quiet || config.sample_paths.len() > 1;

    let num_groups = group_names.len();
    let bucket_count = num_groups + 2;
    let mut combined_buckets: Vec<BucketState> = vec![BucketState::default(); bucket_count];
    let mut total_seqs = 0u64;
    let mut total_bp = 0u64;

    for seq_path in sample_paths {
        let (file_buckets, file_seqs, file_bp) = process_seqs_file(
            seq_path,
            Arc::clone(&index),
            num_groups,
            kmer_length,
            smer_length,
            config.min_hits,
            config.min_fraction,
            config.threads,
            quiet_sample,
            config.include_all_seqs,
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
        seq_files: sample_paths
            .iter()
            .map(|p| {
                if p.to_string_lossy() == "-" {
                    "stdin".to_string()
                } else {
                    p.to_string_lossy().to_string()
                }
            })
            .collect(),
        total_seqs_processed: total_seqs,
        total_bp_processed: total_bp,
        group_names: group_names.to_vec(),
        bucket_histograms,
        bucket_seqs,
        bucket_bases,
    })
}

pub fn run_length_histogram_analysis(config: &LengthHistogramConfig) -> Result<()> {
    let start_time = Instant::now();
    let version = env!("CARGO_PKG_VERSION").to_string();

    let mut options = format!(
        "k={}, s={}, threads={}",
        config.kmer_length, config.smer_length, config.threads
    );

    if !config.include_all_seqs {
        options.push_str(&format!(
            ", min_hits={}, min_fraction={:.2}",
            config.min_hits, config.min_fraction
        ));
        if config.discriminatory {
            options.push_str(", discriminatory");
        }
    }

    if config.sample_paths.len() > 1 {
        options.push_str(&format!(", samples={}", config.sample_paths.len()));
    }

    if let Some(limit) = config.limit_bp {
        options.push_str(&format!(", limit={}", format_bp(limit as usize)));
    }

    eprintln!("Skope v{}; mode: lenhist; options: {}", version, options);

    // Load or build the classification index
    let index_start = Instant::now();
    let (index, group_names, kmer_length, smer_length, index_source) = if config.include_all_seqs {
        if !config.quiet {
            eprintln!("Groups: none (target filtering disabled, single \"all\" bucket)");
        }
        // Degenerate 1-group "all" setup; the index is never queried in this mode.
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
            "none (target filtering disabled)".to_string(),
        )
    } else if config.index_path.is_dir() {
        let (index, group_names) = build_index_in_memory(
            &config.index_path,
            config.kmer_length,
            config.smer_length,
            config.threads,
            config.quiet,
        )?;
        let src = config.index_path.to_string_lossy().to_string();
        (
            index,
            group_names,
            config.kmer_length,
            config.smer_length,
            src,
        )
    } else {
        let (index, group_names, k, s) = load_index(&config.index_path)?;
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
        (
            index,
            group_names,
            k,
            s,
            config.index_path.to_string_lossy().to_string(),
        )
    };

    if !config.include_all_seqs
        && (kmer_length != config.kmer_length || smer_length != config.smer_length)
        && !config.quiet
    {
        eprintln!(
            "Note: using k={}, s={} from index (CLI k/s ignored)",
            kmer_length, smer_length
        );
    }

    let mut index = index;
    if !config.include_all_seqs && config.discriminatory {
        let removed = apply_discriminatory_filter(&mut index);
        if !config.quiet {
            eprintln!(
                "Discriminatory mode: removed {} shared syncmers, {} unique syncmers remain",
                removed,
                index.len()
            );
        }
    }

    let index = Arc::new(index);
    let index_time = index_start.elapsed();

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
            let result = process_single_sample(
                sample_paths,
                sample_name,
                Arc::clone(&index),
                &group_names,
                kmer_length,
                smer_length,
                config,
            );

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

    let seqs_time = start_time.elapsed().as_secs_f64() - index_time.as_secs_f64();
    let total_time = start_time.elapsed();

    let report = LengthHistogramReport {
        version: format!("skope {}", version),
        index_source,
        parameters: LengthHistogramParameters {
            kmer_length,
            smer_length,
            threads: config.threads,
            min_hits: config.min_hits,
            min_fraction: config.min_fraction,
            discriminatory: config.discriminatory,
        },
        samples: sample_results,
        timing: TimingStats {
            reference_processing_time: index_time.as_secs_f64(),
            seqs_processing_time: seqs_time,
            analysis_time: 0.0,
            total_time: total_time.as_secs_f64(),
            seqs_per_second: 0.0,
            bp_per_second: 0.0,
        },
    };

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

    for sample in &report.samples {
        let n = sample.group_names.len();
        for (bucket_idx, hist) in sample.bucket_histograms.iter().enumerate() {
            if hist.is_empty() {
                continue;
            }
            let group_label: &str = if bucket_idx < n {
                &sample.group_names[bucket_idx]
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

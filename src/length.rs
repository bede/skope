use crate::{
    ProcessingStats, RapidHashSet, create_spinner, format_bp, format_bp_per_sec,
    handle_process_result, reader_with_inferred_batch_size,
};
use crate::containment::{MinimizerSet, TimingStats, process_targets_file};
use crate::minimizers::{Buffers, KmerHasher, MinimizerVec, fill_syncmers};
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

/// Result for a single sample's length histogram
#[derive(Debug, Clone)]
pub struct LengthHistogramResult {
    pub sample_name: String,
    pub reads_files: Vec<String>,
    pub total_reads_processed: u64,
    pub total_bp_processed: u64,
    pub reads_with_hits: u64,
    pub reads_without_hits: u64,
    pub length_histogram: Vec<(usize, usize)>, // (length, count) pairs
}

/// Complete report containing all samples
#[derive(Debug, Clone)]
pub struct LengthHistogramReport {
    pub version: String,
    pub targets_file: String,
    pub parameters: LengthHistogramParameters,
    pub samples: Vec<LengthHistogramResult>,
    pub timing: TimingStats,
}

/// Parameters used for analysis
#[derive(Debug, Clone)]
pub struct LengthHistogramParameters {
    pub kmer_length: u8,
    pub smer_length: u8,
    pub threads: usize,
}

/// Configuration for length histogram analysis
pub struct LengthHistogramConfig {
    pub targets_path: PathBuf,
    pub reads_paths: Vec<Vec<PathBuf>>, // Each sample is a Vec of file paths
    pub sample_names: Vec<String>,
    pub kmer_length: u8,
    pub smer_length: u8,
    pub threads: usize,
    pub output_path: Option<PathBuf>,
    pub quiet: bool,
    pub limit_bp: Option<u64>,
    pub include_all_reads: bool, // true when targets_path == "-" (no filtering)
}

impl LengthHistogramConfig {
    pub fn execute(&self) -> Result<()> {
        run_length_histogram_analysis(self)
    }
}

/// Processor for collecting read lengths with syncmer hits
#[derive(Clone)]
struct LengthHistogramProcessor {
    kmer_length: u8,
    smer_length: u8,
    hasher: KmerHasher,
    targets_minimizers: Arc<MinimizerSet>,
    include_all_reads: bool,

    // Local buffers
    buffers: Buffers,
    local_stats: ProcessingStats,
    local_histogram: HashMap<usize, usize>, // length -> count
    local_reads_with_hits: u64,
    local_reads_without_hits: u64,

    // Global state
    global_stats: Arc<Mutex<ProcessingStats>>,
    global_histogram: Arc<Mutex<HashMap<usize, usize>>>,
    global_reads_with_hits: Arc<Mutex<u64>>,
    global_reads_without_hits: Arc<Mutex<u64>>,
    spinner: Option<Arc<Mutex<ProgressBar>>>,
    start_time: Instant,
    limit_bp: Option<u64>,
}

impl LengthHistogramProcessor {
    fn new(
        kmer_length: u8,
        smer_length: u8,
        targets_minimizers: Arc<MinimizerSet>,
        include_all_reads: bool,
        spinner: Option<Arc<Mutex<ProgressBar>>>,
        start_time: Instant,
        limit_bp: Option<u64>,
    ) -> Self {
        let buffers = if kmer_length <= 32 {
            Buffers::new_u64()
        } else {
            Buffers::new_u128()
        };

        Self {
            kmer_length,
            smer_length,
            hasher: KmerHasher::new(smer_length as usize),
            targets_minimizers,
            include_all_reads,
            buffers,
            local_stats: ProcessingStats::default(),
            local_histogram: HashMap::new(),
            local_reads_with_hits: 0,
            local_reads_without_hits: 0,
            global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
            global_histogram: Arc::new(Mutex::new(HashMap::new())),
            global_reads_with_hits: Arc::new(Mutex::new(0)),
            global_reads_without_hits: Arc::new(Mutex::new(0)),
            spinner,
            start_time,
            limit_bp,
        }
    }

    fn update_spinner(&self) {
        if let Some(ref spinner) = self.spinner {
            let stats = self.global_stats.lock();
            let elapsed = self.start_time.elapsed();
            let reads_per_sec = stats.total_seqs as f64 / elapsed.as_secs_f64();
            let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();

            spinner.lock().set_message(format!(
                "Processing sample: {} reads ({}). {:.0} reads/s ({})",
                stats.total_seqs,
                format_bp(stats.total_bp as usize),
                reads_per_sec,
                format_bp_per_sec(bp_per_sec)
            ));
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for LengthHistogramProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        // Check if we've reached sample limit
        if let Some(limit) = self.limit_bp {
            let global_bp = self.global_stats.lock().total_bp;
            if global_bp >= limit {
                return Err(paraseq::parallel::ProcessError::IoError(
                    std::io::Error::new(std::io::ErrorKind::Interrupted, "Sample limit reached"),
                ));
            }
        }

        let seq = record.seq();
        let read_length = seq.len();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += read_length as u64;

        let has_hit = if self.include_all_reads {
            // Include all reads when no target filtering
            true
        } else {
            fill_syncmers(
                &seq,
                &self.hasher,
                self.kmer_length,
                self.smer_length,
                &mut self.buffers,
            );

            // Check if ANY syncmer hits the reference set (early termination)
            match (&self.buffers.minimizers, &*self.targets_minimizers) {
                (MinimizerVec::U64(vec), MinimizerSet::U64(targets_set)) => {
                    vec.iter().any(|minimizer| targets_set.contains(minimizer))
                }
                (MinimizerVec::U128(vec), MinimizerSet::U128(targets_set)) => {
                    vec.iter().any(|minimizer| targets_set.contains(minimizer))
                }
                _ => panic!("Mismatch between MinimizerVec and MinimizerSet types"),
            }
        };

        if has_hit {
            *self.local_histogram.entry(read_length).or_insert(0) += 1;
            self.local_reads_with_hits += 1;
        } else {
            self.local_reads_without_hits += 1;
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Merge local histogram into global
        {
            let mut global_hist = self.global_histogram.lock();
            for (&length, &count) in &self.local_histogram {
                *global_hist.entry(length).or_insert(0) += count;
            }
            self.local_histogram.clear();
        }

        // Merge hit counts
        {
            *self.global_reads_with_hits.lock() += self.local_reads_with_hits;
            *self.global_reads_without_hits.lock() += self.local_reads_without_hits;
            self.local_reads_with_hits = 0;
            self.local_reads_without_hits = 0;
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

fn process_reads_file(
    reads_path: &Path,
    targets_minimizers: Arc<MinimizerSet>,
    kmer_length: u8,
    smer_length: u8,
    threads: usize,
    quiet: bool,
    include_all_reads: bool,
    limit_bp: Option<u64>,
) -> Result<(HashMap<usize, usize>, u64, u64, u64, u64)> {
    let in_path = if reads_path.to_string_lossy() == "-" {
        None
    } else {
        Some(reads_path)
    };
    let reader = reader_with_inferred_batch_size(in_path)?;

    let spinner = create_spinner(quiet)?;
    if let Some(ref pb) = spinner {
        pb.lock().set_message("Processing sample: 0 reads (0bp)");
    }

    let start_time = Instant::now();
    let mut processor = LengthHistogramProcessor::new(
        kmer_length,
        smer_length,
        targets_minimizers,
        include_all_reads,
        spinner.clone(),
        start_time,
        limit_bp,
    );

    let process_result = reader.process_parallel(&mut processor, threads);
    handle_process_result(process_result)?;

    if let Some(ref pb) = spinner {
        pb.lock().finish_with_message("");
    }

    let stats = processor.global_stats.lock().clone();
    let histogram = Arc::try_unwrap(processor.global_histogram)
        .unwrap()
        .into_inner();
    let reads_with_hits = *processor.global_reads_with_hits.lock();
    let reads_without_hits = *processor.global_reads_without_hits.lock();

    if !quiet {
        let elapsed = start_time.elapsed();
        let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();
        eprintln!(
            "Sample: {} records ({}), {} with hits, {} without hits ({})",
            stats.total_seqs,
            format_bp(stats.total_bp as usize),
            reads_with_hits,
            reads_without_hits,
            format_bp_per_sec(bp_per_sec)
        );
    }

    Ok((
        histogram,
        stats.total_seqs,
        stats.total_bp,
        reads_with_hits,
        reads_without_hits,
    ))
}

/// Process a single sample's reads and calculate length histogram
fn process_single_sample(
    _idx: usize,
    reads_paths: &[PathBuf],
    sample_name: &str,
    targets_minimizers: Arc<MinimizerSet>,
    config: &LengthHistogramConfig,
) -> Result<LengthHistogramResult> {
    // Silence per-sample progress for >1 sample
    let quiet_sample = config.quiet || config.reads_paths.len() > 1;

    // Initialize empty histogram
    let mut combined_histogram: HashMap<usize, usize> = HashMap::new();
    let mut total_reads = 0u64;
    let mut total_bp = 0u64;
    let mut total_reads_with_hits = 0u64;
    let mut total_reads_without_hits = 0u64;

    // Process each file and accumulate results
    for reads_path in reads_paths {
        let (file_histogram, file_reads, file_bp, file_with_hits, file_without_hits) =
            process_reads_file(
                reads_path,
                Arc::clone(&targets_minimizers),
                config.kmer_length,
                config.smer_length,
                config.threads,
                quiet_sample,
                config.include_all_reads,
                config.limit_bp.map(|limit| limit.saturating_sub(total_bp)),
            )?;

        // Merge histograms
        for (length, count) in file_histogram {
            *combined_histogram.entry(length).or_insert(0) += count;
        }

        total_reads += file_reads;
        total_bp += file_bp;
        total_reads_with_hits += file_with_hits;
        total_reads_without_hits += file_without_hits;

        // Check if limit reached
        if let Some(limit) = config.limit_bp
            && total_bp >= limit {
                break;
            }
    }

    // Convert histogram to sorted vec
    let mut length_histogram: Vec<(usize, usize)> = combined_histogram.into_iter().collect();
    length_histogram.sort_by_key(|(length, _)| *length);

    Ok(LengthHistogramResult {
        sample_name: sample_name.to_string(),
        reads_files: reads_paths
            .iter()
            .map(|p| {
                if p.to_string_lossy() == "-" {
                    "stdin".to_string()
                } else {
                    p.to_string_lossy().to_string()
                }
            })
            .collect(),
        total_reads_processed: total_reads,
        total_bp_processed: total_bp,
        reads_with_hits: total_reads_with_hits,
        reads_without_hits: total_reads_without_hits,
        length_histogram,
    })
}

pub fn run_length_histogram_analysis(config: &LengthHistogramConfig) -> Result<()> {
    let start_time = Instant::now();
    let version = env!("CARGO_PKG_VERSION").to_string();

    if !config.quiet {
        let mut options = format!(
            "k={}, s={}, threads={}",
            config.kmer_length, config.smer_length, config.threads
        );

        if config.reads_paths.len() > 1 {
            options.push_str(&format!(", samples={}", config.reads_paths.len()));
        }

        if let Some(limit) = config.limit_bp {
            options.push_str(&format!(", limit={}", format_bp(limit as usize)));
        }

        eprintln!("Grate v{}; mode: len; options: {}", version, options);
    }

    // Process targets file OR skip if including all reads
    let (targets_minimizers, targets_time) = if config.include_all_reads {
        if !config.quiet {
            eprintln!("Targets: none (target filtering disabled)");
        }
        // Create empty set - won't be used

        let empty_set = if config.kmer_length <= 32 {
            MinimizerSet::U64(RapidHashSet::default())
        } else {
            MinimizerSet::U128(RapidHashSet::default())
        };
        (Arc::new(empty_set), std::time::Duration::from_secs(0))
    } else {
        // Process targets file
        let targets_start = Instant::now();
        let targets = process_targets_file(
            &config.targets_path,
            config.kmer_length,
            config.smer_length,
            config.quiet,
        )?;
        let targets_time = targets_start.elapsed();

        // Build set of all unique syncmers across targets
        if !config.quiet {
            eprint!("Building syncmer set…\r");
        }
        let targets_minimizers = if config.kmer_length <= 32 {
    
            let mut set = RapidHashSet::default();
            for target in &targets {
                if let MinimizerSet::U64(target_set) = &target.minimizers {
                    set.extend(target_set.iter());
                }
            }
            MinimizerSet::U64(set)
        } else {
    
            let mut set = RapidHashSet::default();
            for target in &targets {
                if let MinimizerSet::U128(target_set) = &target.minimizers {
                    set.extend(target_set.iter());
                }
            }
            MinimizerSet::U128(set)
        };

        let total_target_syncmers = targets_minimizers.len();
        let total_bp: usize = targets.iter().map(|t| t.length).sum();

        if !config.quiet {
            eprint!("\x1B[2K\r"); // Clear line
            eprintln!(
                "Targets: {} records ({}), {} unique syncmers",
                targets.len(),
                format_bp(total_bp),
                total_target_syncmers
            );
        }

        (Arc::new(targets_minimizers), targets_time)
    };

    // Process each sample in parallel
    use rayon::prelude::*;

    let is_multisample = config.reads_paths.len() > 1;
    let completed = if is_multisample && !config.quiet {
        eprint!(
            "\x1B[2K\rSamples: processed 0 of {}…",
            config.reads_paths.len()
        );
        Some(Arc::new(Mutex::new(0usize)))
    } else {
        if !config.quiet {
            eprint!("\x1B[2K\r");
        }
        None
    };

    let sample_results: Vec<LengthHistogramResult> = config
        .reads_paths
        .par_iter()
        .zip(&config.sample_names)
        .enumerate()
        .map(|(idx, (reads_paths, sample_name))| {
            let result = process_single_sample(
                idx,
                reads_paths,
                sample_name,
                Arc::clone(&targets_minimizers),
                config,
            );

            // Report completion for multisample runs
            if let Some(ref counter) = completed {
                let mut count = counter.lock();
                *count += 1;
                eprint!(
                    "\rSamples: processed {} of {}…",
                    *count,
                    config.reads_paths.len()
                );
            }

            result
        })
        .collect::<Result<Vec<_>>>()?;

    if is_multisample && !config.quiet {
        eprintln!();
    }

    let reads_time = start_time.elapsed().as_secs_f64() - targets_time.as_secs_f64();
    let total_time = start_time.elapsed();

    // Create report
    let report = LengthHistogramReport {
        version: format!("grate {}", version),
        targets_file: if config.include_all_reads {
            "none (target filtering disabled)".to_string()
        } else {
            config.targets_path.to_string_lossy().to_string()
        },
        parameters: LengthHistogramParameters {
            kmer_length: config.kmer_length,
            smer_length: config.smer_length,
            threads: config.threads,
        },
        samples: sample_results,
        timing: TimingStats {
            reference_processing_time: targets_time.as_secs_f64(),
            reads_processing_time: reads_time,
            analysis_time: 0.0,
            total_time: total_time.as_secs_f64(),
            reads_per_second: 0.0,
            bp_per_second: 0.0,
        },
    };

    // Output CSV
    let writer: Box<dyn Write> = if let Some(path) = &config.output_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };

    let mut csv_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);

    // Write header
    csv_writer.write_record([
        "sample",
        "length",
        "count",
        "total_reads_processed",
        "total_bp_processed",
        "reads_with_hits",
        "reads_without_hits",
    ])?;

    // Write one row per (sample, length, count) tuple
    for sample in &report.samples {
        for (length, count) in &sample.length_histogram {
            csv_writer.write_record([
                &sample.sample_name,
                &length.to_string(),
                &count.to_string(),
                &sample.total_reads_processed.to_string(),
                &sample.total_bp_processed.to_string(),
                &sample.reads_with_hits.to_string(),
                &sample.reads_without_hits.to_string(),
            ])?;
        }
    }

    csv_writer.flush()?;

    Ok(())
}

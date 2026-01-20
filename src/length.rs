use crate::containment::{
    format_bp, format_bp_per_sec, process_targets_file, MinimizerSet, TargetInfo, TimingStats,
};
use crate::minimizers::{fill_minimizers, Buffers, KmerHasher, MinimizerVec};
use anyhow::Result;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use paraseq::fastx::Reader;
use paraseq::parallel::{ParallelProcessor, ParallelReader};
use paraseq::Record;
use parking_lot::Mutex;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

/// Per-target length histogram result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TargetLengthResult {
    pub target: String,
    pub length: usize,
    pub reads_with_hits: u64,
    pub bp_from_hits: u64,
    pub length_histogram: Vec<(usize, usize)>, // (length, count) pairs
}

/// Aggregate statistics across all targets
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TotalLengthStats {
    pub total_targets: usize,
    pub reads_with_hits: u64,      // Deduplicated across targets
    pub reads_without_hits: u64,
    pub bp_from_hits: u64,
    pub length_histogram: Vec<(usize, usize)>, // (length, count) pairs
}

/// Result for a single sample's length histogram (updated structure)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleLengthResults {
    pub sample_name: String,
    pub reads_files: Vec<String>,
    pub total_reads_processed: u64,
    pub total_bp_processed: u64,
    pub targets: Vec<TargetLengthResult>,
    pub total_stats: TotalLengthStats,
    pub timing: TimingStats,
}

/// Complete report containing all samples
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LengthHistogramReport {
    pub version: String,
    pub targets_file: String,
    pub parameters: LengthHistogramParameters,
    pub samples: Vec<SampleLengthResults>,
    pub timing: TimingStats,
}

/// Parameters used for analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LengthHistogramParameters {
    pub kmer_length: u8,
    pub window_size: u8,
    pub threads: usize,
    pub discriminatory: bool,
}

/// Configuration for length histogram analysis
pub struct LengthHistogramConfig {
    pub targets_path: PathBuf,
    pub reads_paths: Vec<Vec<PathBuf>>, // Each sample is a Vec of file paths
    pub sample_names: Vec<String>,
    pub kmer_length: u8,
    pub window_size: u8,
    pub threads: usize,
    pub output_path: Option<PathBuf>,
    pub quiet: bool,
    pub limit_bp: Option<u64>,
    pub include_all_reads: bool, // true when targets_path == "-" (no filtering)
    pub discriminatory: bool,    // Skip reads that hit multiple targets
}

impl LengthHistogramConfig {
    pub fn execute(&self) -> Result<()> {
        run_length_histogram_analysis(self)
    }
}

fn reader_with_inferred_batch_size(
    in_path: Option<&Path>,
) -> Result<Reader<Box<dyn std::io::Read + Send>>> {
    let mut reader = paraseq::fastx::Reader::from_optional_path(in_path)?;
    reader.update_batch_size_in_bp(256 * 1024)?;
    Ok(reader)
}

#[derive(Clone, Default, Debug)]
struct ProcessingStats {
    total_seqs: u64,
    total_bp: u64,
    last_reported: u64,
}

/// Reverse index mapping minimizers to their target indices
#[derive(Clone)]
enum MinimizerToTargets {
    U64(HashMap<u64, SmallVec<[usize; 2]>>),
    U128(HashMap<u128, SmallVec<[usize; 2]>>),
}

impl MinimizerToTargets {
    /// Build the reverse index from targets
    fn build(targets: &[TargetInfo]) -> Self {
        // Determine type from first target
        let is_u64 = targets
            .first()
            .map(|t| t.minimizers.is_u64())
            .unwrap_or(true);

        if is_u64 {
            let mut map: HashMap<u64, SmallVec<[usize; 2]>> = HashMap::new();
            for (idx, target) in targets.iter().enumerate() {
                if let MinimizerSet::U64(set) = &target.minimizers {
                    for &minimizer in set {
                        map.entry(minimizer).or_default().push(idx);
                    }
                }
            }
            MinimizerToTargets::U64(map)
        } else {
            let mut map: HashMap<u128, SmallVec<[usize; 2]>> = HashMap::new();
            for (idx, target) in targets.iter().enumerate() {
                if let MinimizerSet::U128(set) = &target.minimizers {
                    for &minimizer in set {
                        map.entry(minimizer).or_default().push(idx);
                    }
                }
            }
            MinimizerToTargets::U128(map)
        }
    }
}

/// Processor for collecting read lengths with minimizer hits (per-target tracking)
#[derive(Clone)]
struct LengthHistogramProcessor {
    kmer_length: u8,
    window_size: u8,
    hasher: KmerHasher,
    minimizer_to_targets: Arc<MinimizerToTargets>,
    _num_targets: usize,
    include_all_reads: bool,
    discriminatory: bool,

    // Local buffers
    buffers: Buffers,
    local_stats: ProcessingStats,
    local_per_target_histograms: Vec<HashMap<usize, usize>>,  // Per-target: length -> count
    local_per_target_bp: Vec<u64>,                            // Per-target: total bp from hits
    local_total_histogram: HashMap<usize, usize>,             // Total (deduplicated): length -> count
    local_total_bp: u64,                                       // Total bp from hits
    local_reads_with_hits: u64,
    local_reads_without_hits: u64,

    // Global state
    global_stats: Arc<Mutex<ProcessingStats>>,
    global_per_target_histograms: Arc<Mutex<Vec<HashMap<usize, usize>>>>,
    global_per_target_bp: Arc<Mutex<Vec<u64>>>,
    global_total_histogram: Arc<Mutex<HashMap<usize, usize>>>,
    global_total_bp: Arc<Mutex<u64>>,
    global_reads_with_hits: Arc<Mutex<u64>>,
    global_reads_without_hits: Arc<Mutex<u64>>,
    spinner: Option<Arc<Mutex<ProgressBar>>>,
    start_time: Instant,
    limit_bp: Option<u64>,
}

impl LengthHistogramProcessor {
    fn new(
        kmer_length: u8,
        window_size: u8,
        minimizer_to_targets: Arc<MinimizerToTargets>,
        num_targets: usize,
        include_all_reads: bool,
        discriminatory: bool,
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
            window_size,
            hasher: KmerHasher::new(kmer_length as usize),
            minimizer_to_targets,
            _num_targets: num_targets,
            include_all_reads,
            discriminatory,
            buffers,
            local_stats: ProcessingStats::default(),
            local_per_target_histograms: (0..num_targets).map(|_| HashMap::new()).collect(),
            local_per_target_bp: vec![0u64; num_targets],
            local_total_histogram: HashMap::new(),
            local_total_bp: 0,
            local_reads_with_hits: 0,
            local_reads_without_hits: 0,
            global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
            global_per_target_histograms: Arc::new(Mutex::new(
                (0..num_targets).map(|_| HashMap::new()).collect(),
            )),
            global_per_target_bp: Arc::new(Mutex::new(vec![0u64; num_targets])),
            global_total_histogram: Arc::new(Mutex::new(HashMap::new())),
            global_total_bp: Arc::new(Mutex::new(0)),
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

        if self.include_all_reads {
            // Include all reads when no target filtering - count in TOTAL only
            *self.local_total_histogram.entry(read_length).or_insert(0) += 1;
            self.local_total_bp += read_length as u64;
            self.local_reads_with_hits += 1;
            return Ok(());
        }

        // Extract minimizers from read
        fill_minimizers(
            &seq,
            &self.hasher,
            self.kmer_length,
            self.window_size,
            &mut self.buffers,
        );

        // Find ALL matching targets via reverse index
        let mut matching_targets: SmallVec<[usize; 4]> = SmallVec::new();

        match (&self.buffers.minimizers, &*self.minimizer_to_targets) {
            (MinimizerVec::U64(vec), MinimizerToTargets::U64(index)) => {
                for minimizer in vec {
                    if let Some(targets) = index.get(minimizer) {
                        for &target_idx in targets {
                            if !matching_targets.contains(&target_idx) {
                                matching_targets.push(target_idx);
                            }
                        }
                    }
                }
            }
            (MinimizerVec::U128(vec), MinimizerToTargets::U128(index)) => {
                for minimizer in vec {
                    if let Some(targets) = index.get(minimizer) {
                        for &target_idx in targets {
                            if !matching_targets.contains(&target_idx) {
                                matching_targets.push(target_idx);
                            }
                        }
                    }
                }
            }
            _ => panic!("Mismatch between MinimizerVec and MinimizerToTargets types"),
        }

        // If discriminatory mode and read hits multiple targets, skip it
        if self.discriminatory && matching_targets.len() > 1 {
            self.local_reads_without_hits += 1;
            return Ok(());
        }

        if matching_targets.is_empty() {
            self.local_reads_without_hits += 1;
        } else {
            // Update per-target histograms for each matching target
            for &target_idx in &matching_targets {
                *self.local_per_target_histograms[target_idx]
                    .entry(read_length)
                    .or_insert(0) += 1;
                self.local_per_target_bp[target_idx] += read_length as u64;
            }

            // Update TOTAL histogram once (automatic deduplication)
            *self.local_total_histogram.entry(read_length).or_insert(0) += 1;
            self.local_total_bp += read_length as u64;
            self.local_reads_with_hits += 1;
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Merge local per-target histograms into global
        {
            let mut global_histograms = self.global_per_target_histograms.lock();
            let mut global_bp = self.global_per_target_bp.lock();
            for (idx, local_hist) in self.local_per_target_histograms.iter_mut().enumerate() {
                for (&length, &count) in local_hist.iter() {
                    *global_histograms[idx].entry(length).or_insert(0) += count;
                }
                local_hist.clear();
                global_bp[idx] += self.local_per_target_bp[idx];
                self.local_per_target_bp[idx] = 0;
            }
        }

        // Merge local total histogram into global
        {
            let mut global_hist = self.global_total_histogram.lock();
            for (&length, &count) in &self.local_total_histogram {
                *global_hist.entry(length).or_insert(0) += count;
            }
            self.local_total_histogram.clear();

            *self.global_total_bp.lock() += self.local_total_bp;
            self.local_total_bp = 0;
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

/// Result from processing a single reads file
struct FileProcessingResult {
    per_target_histograms: Vec<HashMap<usize, usize>>,
    per_target_bp: Vec<u64>,
    total_histogram: HashMap<usize, usize>,
    total_bp: u64,
    total_reads: u64,
    total_bp_processed: u64,
    reads_with_hits: u64,
    reads_without_hits: u64,
}

fn process_reads_file(
    reads_path: &Path,
    minimizer_to_targets: Arc<MinimizerToTargets>,
    num_targets: usize,
    kmer_length: u8,
    window_size: u8,
    threads: usize,
    quiet: bool,
    include_all_reads: bool,
    discriminatory: bool,
    limit_bp: Option<u64>,
) -> Result<FileProcessingResult> {
    let in_path = if reads_path.to_string_lossy() == "-" {
        None
    } else {
        Some(reads_path)
    };
    let reader = reader_with_inferred_batch_size(in_path)?;

    // Progress bar
    let spinner = if !quiet {
        let pb = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr());
        pb.set_style(
            ProgressStyle::default_spinner()
                .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                .template("{msg}")?,
        );
        pb.set_message("Processing sample: 0 reads (0bp)");
        Some(Arc::new(Mutex::new(pb)))
    } else {
        None
    };

    let start_time = Instant::now();
    let mut processor = LengthHistogramProcessor::new(
        kmer_length,
        window_size,
        minimizer_to_targets,
        num_targets,
        include_all_reads,
        discriminatory,
        spinner.clone(),
        start_time,
        limit_bp,
    );

    // Process reads - may terminate early if sample limit reached
    let process_result = reader.process_parallel(&mut processor, threads);

    // Check if we stopped due to sampling
    let stopped_early = if let Err(ref e) = process_result {
        e.to_string().contains("Sample limit reached")
    } else {
        false
    };

    // If it's not a sample stop, propagate the error
    if !stopped_early {
        process_result?;
    }

    // Finish spinner
    if let Some(ref pb) = spinner {
        pb.lock().finish_with_message("");
    }

    let stats = processor.global_stats.lock().clone();
    let per_target_histograms = Arc::try_unwrap(processor.global_per_target_histograms)
        .unwrap()
        .into_inner();
    let per_target_bp = Arc::try_unwrap(processor.global_per_target_bp)
        .unwrap()
        .into_inner();
    let total_histogram = Arc::try_unwrap(processor.global_total_histogram)
        .unwrap()
        .into_inner();
    let total_bp = *processor.global_total_bp.lock();
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

    Ok(FileProcessingResult {
        per_target_histograms,
        per_target_bp,
        total_histogram,
        total_bp,
        total_reads: stats.total_seqs,
        total_bp_processed: stats.total_bp,
        reads_with_hits,
        reads_without_hits,
    })
}

/// Process a single sample's reads and calculate length histogram
fn process_single_sample(
    _idx: usize,
    reads_paths: &[PathBuf],
    sample_name: &str,
    targets: &[TargetInfo],
    minimizer_to_targets: Arc<MinimizerToTargets>,
    config: &LengthHistogramConfig,
) -> Result<SampleLengthResults> {
    let reads_start = Instant::now();

    // Silence per-sample progress for >1 sample
    let quiet_sample = config.quiet || config.reads_paths.len() > 1;

    let num_targets = targets.len();

    // Initialize accumulators
    let mut combined_per_target_histograms: Vec<HashMap<usize, usize>> =
        (0..num_targets).map(|_| HashMap::new()).collect();
    let mut combined_per_target_bp: Vec<u64> = vec![0u64; num_targets];
    let mut combined_total_histogram: HashMap<usize, usize> = HashMap::new();
    let mut combined_total_bp = 0u64;
    let mut total_reads = 0u64;
    let mut total_bp_processed = 0u64;
    let mut total_reads_with_hits = 0u64;
    let mut total_reads_without_hits = 0u64;

    // Process each file and accumulate results
    for reads_path in reads_paths {
        let result = process_reads_file(
            reads_path,
            Arc::clone(&minimizer_to_targets),
            num_targets,
            config.kmer_length,
            config.window_size,
            config.threads,
            quiet_sample,
            config.include_all_reads,
            config.discriminatory,
            config.limit_bp.map(|limit| limit.saturating_sub(total_bp_processed)),
        )?;

        // Merge per-target histograms
        for (idx, hist) in result.per_target_histograms.into_iter().enumerate() {
            for (length, count) in hist {
                *combined_per_target_histograms[idx].entry(length).or_insert(0) += count;
            }
            combined_per_target_bp[idx] += result.per_target_bp[idx];
        }

        // Merge total histogram
        for (length, count) in result.total_histogram {
            *combined_total_histogram.entry(length).or_insert(0) += count;
        }
        combined_total_bp += result.total_bp;

        total_reads += result.total_reads;
        total_bp_processed += result.total_bp_processed;
        total_reads_with_hits += result.reads_with_hits;
        total_reads_without_hits += result.reads_without_hits;

        // Check if limit reached
        if let Some(limit) = config.limit_bp {
            if total_bp_processed >= limit {
                break;
            }
        }
    }

    let reads_time = reads_start.elapsed();

    // Build per-target results
    let target_results: Vec<TargetLengthResult> = targets
        .iter()
        .enumerate()
        .map(|(idx, target)| {
            let hist = &combined_per_target_histograms[idx];
            let reads_with_hits: u64 = hist.values().map(|&c| c as u64).sum();

            let mut length_histogram: Vec<(usize, usize)> = hist.iter().map(|(&l, &c)| (l, c)).collect();
            length_histogram.sort_by_key(|(length, _)| *length);

            TargetLengthResult {
                target: target.name.clone(),
                length: target.length,
                reads_with_hits,
                bp_from_hits: combined_per_target_bp[idx],
                length_histogram,
            }
        })
        .collect();

    // Build total stats
    let mut total_histogram_vec: Vec<(usize, usize)> =
        combined_total_histogram.into_iter().collect();
    total_histogram_vec.sort_by_key(|(length, _)| *length);

    let total_stats = TotalLengthStats {
        total_targets: num_targets,
        reads_with_hits: total_reads_with_hits,
        reads_without_hits: total_reads_without_hits,
        bp_from_hits: combined_total_bp,
        length_histogram: total_histogram_vec,
    };

    let reads_per_second = total_reads as f64 / reads_time.as_secs_f64();
    let bp_per_second = total_bp_processed as f64 / reads_time.as_secs_f64();

    Ok(SampleLengthResults {
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
        total_bp_processed,
        targets: target_results,
        total_stats,
        timing: TimingStats {
            reference_processing_time: 0.0, // Not per-sample
            reads_processing_time: reads_time.as_secs_f64(),
            analysis_time: 0.0,
            total_time: reads_time.as_secs_f64(),
            reads_per_second,
            bp_per_second,
        },
    })
}

pub fn run_length_histogram_analysis(config: &LengthHistogramConfig) -> Result<()> {
    let start_time = Instant::now();
    let version = env!("CARGO_PKG_VERSION").to_string();

    if !config.quiet {
        let mut options = format!(
            "k={}, w={}, threads={}",
            config.kmer_length, config.window_size, config.threads
        );

        if config.reads_paths.len() > 1 {
            options.push_str(&format!(", samples={}", config.reads_paths.len()));
        }

        if config.discriminatory {
            options.push_str(", discriminatory");
        }

        if let Some(limit) = config.limit_bp {
            options.push_str(&format!(", limit={}", format_bp(limit as usize)));
        }

        eprintln!("Grate v{}; mode: len; options: {}", version, options);
    }

    // Process targets file OR skip if including all reads
    let (targets, minimizer_to_targets, targets_time) = if config.include_all_reads {
        if !config.quiet {
            eprintln!("Targets: none (target filtering disabled)");
        }
        // Create empty structures - won't be used
        let empty_index = if config.kmer_length <= 32 {
            MinimizerToTargets::U64(HashMap::new())
        } else {
            MinimizerToTargets::U128(HashMap::new())
        };
        (Vec::new(), Arc::new(empty_index), std::time::Duration::from_secs(0))
    } else {
        // Process targets file
        let targets_start = Instant::now();
        let targets = process_targets_file(
            &config.targets_path,
            config.kmer_length,
            config.window_size,
            config.quiet,
        )?;
        let targets_time = targets_start.elapsed();

        // Build reverse index: minimizer -> target indices
        if !config.quiet {
            eprint!("Building minimizer index…\r");
        }
        let minimizer_to_targets = MinimizerToTargets::build(&targets);

        let total_unique_minimizers = match &minimizer_to_targets {
            MinimizerToTargets::U64(map) => map.len(),
            MinimizerToTargets::U128(map) => map.len(),
        };
        let total_bp: usize = targets.iter().map(|t| t.length).sum();

        if !config.quiet {
            eprint!("\x1B[2K\r"); // Clear line
            eprintln!(
                "Targets: {} records ({}), {} unique minimizers",
                targets.len(),
                format_bp(total_bp),
                total_unique_minimizers
            );
        }

        (targets, Arc::new(minimizer_to_targets), targets_time)
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

    let sample_results: Vec<SampleLengthResults> = config
        .reads_paths
        .par_iter()
        .zip(&config.sample_names)
        .enumerate()
        .map(|(idx, (reads_paths, sample_name))| {
            let result = process_single_sample(
                idx,
                reads_paths,
                sample_name,
                &targets,
                Arc::clone(&minimizer_to_targets),
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
            window_size: config.window_size,
            threads: config.threads,
            discriminatory: config.discriminatory,
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

    // Output JSON
    let writer: Box<dyn Write> = if let Some(path) = &config.output_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };

    let mut writer = writer;
    serde_json::to_writer_pretty(&mut writer, &report)?;
    writeln!(writer)?;

    Ok(())
}

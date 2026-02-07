use crate::{
    ProcessingStats, RapidHashSet, create_spinner, format_bp, format_bp_per_sec,
    handle_process_result, reader_with_inferred_batch_size, sample_limit_reached_io_error,
};
use crate::minimizers::{
    Buffers, KmerHasher, MinimizerVec, fill_syncmers, fill_syncmers_with_positions,
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

/// Alias for abund counts
type CountDepth = u16;

/// Sort order for results
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    Original,    // Original order from reference file (default)
    Target,      // Alphabetical by target name
    Sample,      // Alphabetical by sample name
    Containment, // Descending by containment1 (highest first)
}

/// Zero-cost (hopefully?) abstraction over u64 and u128 minimizer sets
#[derive(Debug, Clone)]
pub enum MinimizerSet {
    U64(RapidHashSet<u64>),
    U128(RapidHashSet<u128>),
}

impl MinimizerSet {
    pub fn len(&self) -> usize {
        match self {
            MinimizerSet::U64(set) => set.len(),
            MinimizerSet::U128(set) => set.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn is_u64(&self) -> bool {
        matches!(self, MinimizerSet::U64(_))
    }
}

#[derive(Debug, Clone)]
pub struct TargetInfo {
    pub name: String,
    pub length: usize,
    pub minimizers: MinimizerSet,
    pub minimizer_positions: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct ContainmentResult {
    pub target: String,
    pub length: usize,
    pub target_kmers: usize,
    pub contained_minimizers: usize,
    pub containment1: f64,
    pub median_nz_abundance: f64,
    pub abundance_histogram: Vec<(CountDepth, usize)>, // (abundance, count)
    pub containment_at_threshold: HashMap<usize, f64>, // threshold -> containment
    pub hits_at_threshold: HashMap<usize, usize>,      // threshold -> hit count
    pub sample_name: Option<String>,                   // Only used in multi-sample mode
}

#[derive(Debug, Clone)]
pub struct ContainmentParameters {
    pub kmer_length: u8,
    pub smer_length: u8,
    pub threads: usize,
    pub abundance_thresholds: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct TotalStats {
    pub total_targets: usize,
    pub target_kmers: usize,
    pub total_contained_minimizers: usize,
    pub total_containment1: f64,
    pub total_median_nz_abundance: f64,
    pub total_seqs_processed: u64,
    pub total_bp_processed: u64,
    pub total_containment_at_threshold: HashMap<usize, f64>, // threshold -> overall containment
}

#[derive(Debug, Clone)]
pub struct TimingStats {
    pub reference_processing_time: f64,
    pub seqs_processing_time: f64,
    pub analysis_time: f64,
    pub total_time: f64,
    pub seqs_per_second: f64,
    pub bp_per_second: f64,
}

/// Results for a single sample in multi-sample mode
#[derive(Debug, Clone)]
pub struct SampleResults {
    pub sample_name: String,
    pub seq_files: Vec<String>, // Multiple files per sample
    pub targets: Vec<ContainmentResult>,
    pub total_stats: TotalStats,
    pub timing: TimingStats,
}

/// Report containing results for one or more samples
#[derive(Debug, Clone)]
pub struct Report {
    pub version: String,
    pub targets_file: String,
    pub parameters: ContainmentParameters,
    pub samples: Vec<SampleResults>,
    pub total_timing: TimingStats,
}

/// Output format options
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Table,
    Tsv,
}

pub struct ContainmentConfig {
    pub targets_path: PathBuf,
    pub sample_paths: Vec<Vec<PathBuf>>, // Each sample is a Vec of file paths
    pub sample_names: Vec<String>,
    pub kmer_length: u8,
    pub smer_length: u8,
    pub threads: usize,
    pub output_path: Option<PathBuf>,
    pub quiet: bool,
    pub output_format: OutputFormat,
    pub abundance_thresholds: Vec<usize>,
    pub discriminatory: bool,
    pub limit_bp: Option<u64>,
    pub sort_order: SortOrder,
    pub dump_positions_path: Option<PathBuf>,
    pub no_total: bool,
}

impl ContainmentConfig {
    pub fn execute(&self) -> Result<()> {
        run_containment_analysis(self)
    }
}

/// Processor for collecting target sequence records with syncmers
#[derive(Clone)]
struct TargetsProcessor {
    kmer_length: u8,
    smer_length: u8,
    hasher: KmerHasher,
    buffers: Buffers,
    positions: Vec<usize>,
    targets: Arc<Mutex<Vec<TargetInfo>>>,

    // Progress tracking
    local_stats: ProcessingStats,
    global_stats: Arc<Mutex<ProcessingStats>>,
    spinner: Option<Arc<Mutex<ProgressBar>>>,
    start_time: std::time::Instant,
}

impl TargetsProcessor {
    fn new(
        kmer_length: u8,
        smer_length: u8,
        spinner: Option<Arc<Mutex<ProgressBar>>>,
        start_time: std::time::Instant,
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
            buffers,
            positions: Vec::new(),
            targets: Arc::new(Mutex::new(Vec::new())),
            local_stats: ProcessingStats::default(),
            global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
            spinner,
            start_time,
        }
    }

    fn update_spinner(&self) {
        if let Some(ref spinner) = self.spinner {
            let stats = self.global_stats.lock();
            let elapsed = self.start_time.elapsed();
            let seqs_per_sec = stats.total_seqs as f64 / elapsed.as_secs_f64();
            let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();

            spinner.lock().set_message(format!(
                "Collecting target syncmers: {} seqs ({}). {:.0} seqs/s ({})",
                stats.total_seqs,
                format_bp(stats.total_bp as usize),
                seqs_per_sec,
                format_bp_per_sec(bp_per_sec)
            ));
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for TargetsProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        let sequence = record.seq();
        let target_name = String::from_utf8_lossy(record.id()).to_string();

        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += sequence.len() as u64;

        fill_syncmers_with_positions(
            &sequence,
            &self.hasher,
            self.kmer_length,
            self.smer_length,
            &mut self.buffers,
            &mut self.positions,
        );

        // Build unique syncmer set for this target
        let minimizers = match &self.buffers.minimizers {
            MinimizerVec::U64(vec) => {
                let set: RapidHashSet<u64> = vec.iter().copied().collect();
                MinimizerSet::U64(set)
            }
            MinimizerVec::U128(vec) => {
                let set: RapidHashSet<u128> = vec.iter().copied().collect();
                MinimizerSet::U128(set)
            }
        };

        let minimizer_positions = self.positions.clone();

        self.targets.lock().push(TargetInfo {
            name: target_name,
            length: sequence.len(),
            minimizers,
            minimizer_positions,
        });

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;

            // Update spinner every 0.1 Gbp
            let current_progress = stats.total_bp / 100_000_000; // 0.1 Gbp increments
            if current_progress > stats.last_reported {
                drop(stats); // Release lock before updating spinner
                self.update_spinner();
                self.global_stats.lock().last_reported = current_progress;
            }

            self.local_stats = ProcessingStats::default();
        }

        Ok(())
    }
}

pub fn process_targets_file(
    targets_path: &Path,
    kmer_length: u8,
    smer_length: u8,
    quiet: bool,
) -> Result<Vec<TargetInfo>> {
    let in_path = if targets_path.to_string_lossy() == "-" {
        None
    } else {
        Some(targets_path)
    };

    let reader = reader_with_inferred_batch_size(in_path)?;

    let spinner = create_spinner(quiet)?;

    let start_time = std::time::Instant::now();
    let mut processor =
        TargetsProcessor::new(kmer_length, smer_length, spinner.clone(), start_time);

    // Single thread to preserve order
    reader.process_parallel(&mut processor, 1)?;

    // Finish spinner and clear
    if let Some(ref pb) = spinner {
        pb.lock().finish_and_clear();
    }

    let targets = Arc::try_unwrap(processor.targets).unwrap().into_inner();

    Ok(targets)
}

/// Processor for counting syncmer depths from sequences
#[derive(Clone)]
struct SeqsProcessor {
    kmer_length: u8,
    smer_length: u8,
    hasher: KmerHasher,
    targets_minimizers: Arc<MinimizerSet>,

    // Local buffers
    buffers: Buffers,
    local_stats: ProcessingStats,
    local_counts_u64: Option<HashMap<u64, CountDepth>>,
    local_counts_u128: Option<HashMap<u128, CountDepth>>,

    // Global state
    global_stats: Arc<Mutex<ProcessingStats>>,
    global_counts_u64: Arc<Mutex<Option<HashMap<u64, CountDepth>>>>,
    global_counts_u128: Arc<Mutex<Option<HashMap<u128, CountDepth>>>>,
    spinner: Option<Arc<Mutex<ProgressBar>>>,
    start_time: Instant,
    limit_bp: Option<u64>,
}

impl SeqsProcessor {
    fn new(
        kmer_length: u8,
        smer_length: u8,
        targets_minimizers: Arc<MinimizerSet>,
        spinner: Option<Arc<Mutex<ProgressBar>>>,
        start_time: Instant,
        limit_bp: Option<u64>,
    ) -> Self {
        let buffers = if kmer_length <= 32 {
            Buffers::new_u64()
        } else {
            Buffers::new_u128()
        };

        let (local_counts_u64, local_counts_u128) = if kmer_length <= 32 {
            (Some(HashMap::new()), None)
        } else {
            (None, Some(HashMap::new()))
        };

        let (global_counts_u64, global_counts_u128) = if kmer_length <= 32 {
            (
                Arc::new(Mutex::new(Some(HashMap::new()))),
                Arc::new(Mutex::new(None)),
            )
        } else {
            (
                Arc::new(Mutex::new(None)),
                Arc::new(Mutex::new(Some(HashMap::new()))),
            )
        };

        Self {
            kmer_length,
            smer_length,
            hasher: KmerHasher::new(smer_length as usize),
            targets_minimizers,
            buffers,
            local_stats: ProcessingStats::default(),
            local_counts_u64,
            local_counts_u128,
            global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
            global_counts_u64,
            global_counts_u128,
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

impl<Rf: Record> ParallelProcessor<Rf> for SeqsProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
if let Some(limit) = self.limit_bp {            let global_bp = self.global_stats.lock().total_bp;
            if global_bp >= limit {
                return Err(paraseq::parallel::ProcessError::IoError(
                    sample_limit_reached_io_error(),
                ));
            }
        }

        let seq = record.seq();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq.len() as u64;

        fill_syncmers(
            &seq,
            &self.hasher,
            self.kmer_length,
            self.smer_length,
            &mut self.buffers,
        );

        // Count syncmers present in targets
        match (&self.buffers.minimizers, &*self.targets_minimizers) {
            (MinimizerVec::U64(vec), MinimizerSet::U64(targets_set)) => {
                let local_counts = self.local_counts_u64.as_mut().unwrap();
                for &minimizer in vec {
                    if targets_set.contains(&minimizer) {
                        local_counts
                            .entry(minimizer)
                            .and_modify(|e| *e = e.saturating_add(1))
                            .or_insert(1);
                    }
                }
            }
            (MinimizerVec::U128(vec), MinimizerSet::U128(targets_set)) => {
                let local_counts = self.local_counts_u128.as_mut().unwrap();
                for &minimizer in vec {
                    if targets_set.contains(&minimizer) {
                        local_counts
                            .entry(minimizer)
                            .and_modify(|e| *e = e.saturating_add(1))
                            .or_insert(1);
                    }
                }
            }
            _ => panic!("Mismatch between MinimizerVec and MinimizerSet types"),
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Merge local into global counts
        if let Some(local) = &mut self.local_counts_u64 {
            let mut global = self.global_counts_u64.lock();
            let global_map = global.as_mut().unwrap();
            for (&minimizer, &count) in local.iter() {
                global_map
                    .entry(minimizer)
                    .and_modify(|e| *e = e.saturating_add(count))
                    .or_insert(count);
            }
            local.clear();
        } else {
            let mut global = self.global_counts_u128.lock();
            let global_map = global.as_mut().unwrap();
            let local = self.local_counts_u128.as_mut().unwrap();
            for (&minimizer, &count) in local.iter() {
                global_map
                    .entry(minimizer)
                    .and_modify(|e| *e = e.saturating_add(count))
                    .or_insert(count);
            }
            local.clear();
        }

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;

            // Update spinner every 0.1 Gbp
            let current_progress = stats.total_bp / 100_000_000; // 0.1 Gbp increments
            if current_progress > stats.last_reported {
                drop(stats); // Release lock before updating spinner
                self.update_spinner();
                self.global_stats.lock().last_reported = current_progress;
            }

            self.local_stats = ProcessingStats::default();
        }

        Ok(())
    }
}

/// Enum to abstract over u64 and u128 abundance maps
enum AbundanceMap {
    U64(HashMap<u64, CountDepth>),
    U128(HashMap<u128, CountDepth>),
}

fn process_seqs_file(
    seq_path: &Path,
    targets_minimizers: Arc<MinimizerSet>,
    kmer_length: u8,
    smer_length: u8,
    threads: usize,
    quiet: bool,
    limit_bp: Option<u64>,
) -> Result<(AbundanceMap, u64, u64)> {
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

    let total_target_syncmers = targets_minimizers.len();

    let start_time = Instant::now();
    let mut processor = SeqsProcessor::new(
        kmer_length,
        smer_length,
        targets_minimizers,
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
    let abundance_map = if kmer_length <= 32 {
        let map = Arc::try_unwrap(processor.global_counts_u64)
            .unwrap()
            .into_inner()
            .unwrap();
        AbundanceMap::U64(map)
    } else {
        let map = Arc::try_unwrap(processor.global_counts_u128)
            .unwrap()
            .into_inner()
            .unwrap();
        AbundanceMap::U128(map)
    };

    if !quiet {
        let elapsed = start_time.elapsed();
        let unique_syncmers = match &abundance_map {
            AbundanceMap::U64(m) => m.len(),
            AbundanceMap::U128(m) => m.len(),
        };
        let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();
        eprintln!(
            "Sample: {} records ({}), found {} of {} distinct target syncmers ({})",
            stats.total_seqs,
            format_bp(stats.total_bp as usize),
            unique_syncmers,
            total_target_syncmers,
            format_bp_per_sec(bp_per_sec)
        );
    }

    Ok((abundance_map, stats.total_seqs, stats.total_bp))
}

fn calculate_containment_statistics(
    targets: &[TargetInfo],
    abundance_map: &AbundanceMap,
    abundance_thresholds: &[usize],
    sample_name: Option<String>,
) -> Vec<ContainmentResult> {
    targets
        .iter()
        .map(|target| {
            let target_kmers = target.minimizers.len();
            let mut abundances: Vec<CountDepth> = Vec::new();
            let mut contained_count = 0;

            // Collect abundances for all unique minimizers in this target
            match (&target.minimizers, abundance_map) {
                (MinimizerSet::U64(set), AbundanceMap::U64(map)) => {
                    for &minimizer in set {
                        let abundance = map.get(&minimizer).copied().unwrap_or(0);
                        abundances.push(abundance);
                        if abundance > 0 {
                            contained_count += 1;
                        }
                    }
                }
                (MinimizerSet::U128(set), AbundanceMap::U128(map)) => {
                    for &minimizer in set {
                        let abundance = map.get(&minimizer).copied().unwrap_or(0);
                        abundances.push(abundance);
                        if abundance > 0 {
                            contained_count += 1;
                        }
                    }
                }
                _ => panic!("Mismatch between MinimizerSet and AbundanceMap types"),
            }

            let containment1 = if target_kmers > 0 {
                contained_count as f64 / target_kmers as f64
            } else {
                0.0
            };

            // Ignore zero-abundance minimizers for median calc
            let non_zero_abundances: Vec<CountDepth> =
                abundances.iter().copied().filter(|&a| a > 0).collect();

            // Calculate median w/o zeros
            let mut sorted_abundances = non_zero_abundances;
            sorted_abundances.sort_unstable();
            let median_nz_abundance = if sorted_abundances.is_empty() {
                0.0
            } else if sorted_abundances.len().is_multiple_of(2) {
                let mid = sorted_abundances.len() / 2;
                (sorted_abundances[mid - 1] + sorted_abundances[mid]) as f64 / 2.0
            } else {
                sorted_abundances[sorted_abundances.len() / 2] as f64
            };

            // Abundance hist
            let mut abundance_counts: HashMap<CountDepth, usize> = HashMap::new();
            for abundance in &abundances {
                *abundance_counts.entry(*abundance).or_insert(0) += 1;
            }
            let mut abundance_histogram: Vec<(CountDepth, usize)> =
                abundance_counts.into_iter().collect();
            abundance_histogram.sort_by_key(|(abundance, _)| *abundance);

            // Calculate containment at threshold
            let mut containment_at_threshold = HashMap::new();
            let mut hits_at_threshold = HashMap::new();
            for &threshold in abundance_thresholds {
                let count_at_threshold = abundances
                    .iter()
                    .filter(|&&a| a >= threshold as CountDepth)
                    .count();

                // Store the hit count
                hits_at_threshold.insert(threshold, count_at_threshold);

                let containment_value = if target_kmers > 0 {
                    count_at_threshold as f64 / target_kmers as f64
                } else {
                    0.0
                };
                containment_at_threshold.insert(threshold, containment_value);
            }

            ContainmentResult {
                target: target.name.clone(),
                length: target.length,
                target_kmers,
                contained_minimizers: contained_count,
                containment1,
                median_nz_abundance,
                abundance_histogram,
                containment_at_threshold,
                hits_at_threshold,
                sample_name: sample_name.clone(),
            }
        })
        .collect()
}

fn truncate_string(s: &str, max_len: usize) -> String {
    if s.len() <= max_len {
        s.to_string()
    } else {
        format!("{}...", &s[..max_len - 3])
    }
}


/// Process a single sample's sequences and calculate statistics
fn process_single_sample(
    _idx: usize,
    sample_paths: &[PathBuf], // Multiple files per sample
    sample_name: &str,
    targets: &[TargetInfo],
    targets_minimizers: Arc<MinimizerSet>,
    config: &ContainmentConfig,
) -> Result<SampleResults> {
    // Silence per-sample progress for >1 sample
    let quiet_sample = config.quiet || config.sample_paths.len() > 1;

    let seqs_start = Instant::now();

    // Initialise empty abundance map based on k-mer length
    let mut combined_abundance_map = if config.kmer_length <= 32 {
        AbundanceMap::U64(HashMap::new())
    } else {
        AbundanceMap::U128(HashMap::new())
    };
    let mut total_seqs = 0u64;
    let mut total_bp = 0u64;

    // Process each file and accumulate results
    for seq_path in sample_paths {
        let (file_abundance_map, file_seqs, file_bp) = process_seqs_file(
            seq_path,
            Arc::clone(&targets_minimizers),
            config.kmer_length,
            config.smer_length,
            config.threads,
            quiet_sample,
            config.limit_bp.map(|limit| limit.saturating_sub(total_bp)),
        )?;

        // Merge abundance maps
        match (&mut combined_abundance_map, file_abundance_map) {
            (AbundanceMap::U64(combined), AbundanceMap::U64(new)) => {
                for (minimizer, count) in new {
                    combined
                        .entry(minimizer)
                        .and_modify(|e| *e = e.saturating_add(count))
                        .or_insert(count);
                }
            }
            (AbundanceMap::U128(combined), AbundanceMap::U128(new)) => {
                for (minimizer, count) in new {
                    combined
                        .entry(minimizer)
                        .and_modify(|e| *e = e.saturating_add(count))
                        .or_insert(count);
                }
            }
            _ => panic!("Mismatch between AbundanceMap types"),
        }

        total_seqs += file_seqs;
        total_bp += file_bp;

        // Check if limit reached
        if let Some(limit) = config.limit_bp
            && total_bp >= limit {
                break;
            }
    }

    let seqs_time = seqs_start.elapsed();
    let abundance_map = combined_abundance_map;

    // Calculate containment statistics for this sample
    let analysis_start = Instant::now();
    let containment_results = calculate_containment_statistics(
        targets,
        &abundance_map,
        &config.abundance_thresholds,
        Some(sample_name.to_string()),
    );
    let analysis_time = analysis_start.elapsed();

    // Overall stats for this sample
    let target_kmers: usize = containment_results.iter().map(|r| r.target_kmers).sum();
    let total_contained_minimizers: usize = containment_results
        .iter()
        .map(|r| r.contained_minimizers)
        .sum();
    let total_containment1 = if target_kmers > 0 {
        total_contained_minimizers as f64 / target_kmers as f64
    } else {
        0.0
    };

    let total_median_nz_abundance = if target_kmers > 0 {
        containment_results
            .iter()
            .map(|r| r.median_nz_abundance * r.target_kmers as f64)
            .sum::<f64>()
            / target_kmers as f64
    } else {
        0.0
    };

    let seqs_per_second = total_seqs as f64 / seqs_time.as_secs_f64();
    let bp_per_second = total_bp as f64 / seqs_time.as_secs_f64();

    // Calculate overall containment at each threshold
    let mut total_containment_at_threshold = HashMap::new();
    for &threshold in &config.abundance_thresholds {
        let total_at_threshold: usize = containment_results
            .iter()
            .map(|r| {
                let containment = r.containment_at_threshold.get(&threshold).unwrap_or(&0.0);
                (containment * r.target_kmers as f64) as usize
            })
            .sum();
        let overall_containment_value = if target_kmers > 0 {
            total_at_threshold as f64 / target_kmers as f64
        } else {
            0.0
        };
        total_containment_at_threshold.insert(threshold, overall_containment_value);
    }

    Ok(SampleResults {
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
        targets: containment_results,
        total_stats: TotalStats {
            total_targets: targets.len(),
            target_kmers,
            total_contained_minimizers,
            total_containment1,
            total_median_nz_abundance,
            total_seqs_processed: total_seqs,
            total_bp_processed: total_bp,
            total_containment_at_threshold,
        },
        timing: TimingStats {
            reference_processing_time: 0.0, // Not per-sample
            seqs_processing_time: seqs_time.as_secs_f64(),
            analysis_time: analysis_time.as_secs_f64(),
            total_time: seqs_time.as_secs_f64() + analysis_time.as_secs_f64(),
            seqs_per_second,
            bp_per_second,
        },
    })
}

pub fn run_containment_analysis(config: &ContainmentConfig) -> Result<()> {
    let start_time = Instant::now();
    let version = env!("CARGO_PKG_VERSION").to_string();

    if !config.quiet {
        let mut options = format!(
            "k={}, s={}, threads={}",
            config.kmer_length, config.smer_length, config.threads
        );

        if config.sample_paths.len() > 1 {
            options.push_str(&format!(", samples={}", config.sample_paths.len()));
        }

        if !config.abundance_thresholds.is_empty() {
            options.push_str(&format!(
                ", abundance-thresholds={}",
                config
                    .abundance_thresholds
                    .iter()
                    .map(|t| t.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            ));
        }

        if config.discriminatory {
            options.push_str(", discriminatory");
        }

        if let Some(limit) = config.limit_bp {
            options.push_str(&format!(", limit={}", format_bp(limit as usize)));
        }

        options.push_str(&format!(
            ", format={}",
            match config.output_format {
                OutputFormat::Table => "table",
                OutputFormat::Tsv => "tsv",
            }
        ));

        eprintln!("Skope v{}; mode: con; options: {}", version, options);
    }

    // Process targets file
    let targets_start = Instant::now();
    let mut targets = process_targets_file(
        &config.targets_path,
        config.kmer_length,
        config.smer_length,
        config.quiet,
    )?;
    let targets_time = targets_start.elapsed();

    // Count syncmers shared between targets
    if !config.quiet {
        eprint!("Counting shared syncmers…\r");
    }
    let mut minimizer_target_counts: HashMap<u64, usize> = HashMap::new();
    let mut minimizer_target_counts_u128: HashMap<u128, usize> = HashMap::new();

    for target in &targets {
        match &target.minimizers {
            MinimizerSet::U64(set) => {
                for &minimizer in set {
                    *minimizer_target_counts.entry(minimizer).or_insert(0) += 1;
                }
            }
            MinimizerSet::U128(set) => {
                for &minimizer in set {
                    *minimizer_target_counts_u128.entry(minimizer).or_insert(0) += 1;
                }
            }
        }
    }

    let shared_minimizers = if config.kmer_length <= 32 {
        minimizer_target_counts
            .values()
            .filter(|&&count| count > 1)
            .count()
    } else {
        minimizer_target_counts_u128
            .values()
            .filter(|&&count| count > 1)
            .count()
    };

    let unique_across_all = if config.kmer_length <= 32 {
        minimizer_target_counts.len()
    } else {
        minimizer_target_counts_u128.len()
    };

    // Apply discriminatory filtering if enabled
    if config.discriminatory {
        for target in &mut targets {
            match &mut target.minimizers {
                MinimizerSet::U64(set) => {
                    set.retain(|minimizer| {
                        minimizer_target_counts
                            .get(minimizer)
                            .is_none_or(|&count| count == 1)
                    });
                }
                MinimizerSet::U128(set) => {
                    set.retain(|minimizer| {
                        minimizer_target_counts_u128
                            .get(minimizer)
                            .is_none_or(|&count| count == 1)
                    });
                }
            }
        }
    }

    if !config.quiet {
        eprint!("\r"); // Clear space
        let total_unique_syncmers: usize = targets.iter().map(|t| t.minimizers.len()).sum();
        let total_bp: usize = targets.iter().map(|t| t.length).sum();

        if config.discriminatory {
            eprintln!(
                "Targets: {} records ({}), {} discriminatory syncmers ({} shared syncmers dropped)",
                targets.len(),
                format_bp(total_bp),
                total_unique_syncmers,
                shared_minimizers
            );
        } else {
            let shared_pct = if unique_across_all > 0 {
                shared_minimizers as f64 / unique_across_all as f64 * 100.0
            } else {
                0.0
            };
            eprintln!(
                "Targets: {} records ({}), {} syncmers, of which {} ({:.1}%) shared by multiple targets",
                targets.len(),
                format_bp(total_bp),
                total_unique_syncmers,
                shared_minimizers,
                shared_pct
            );
        }
    }

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

    let targets_minimizers = Arc::new(targets_minimizers);

    // Dump syncmer positions if requested
    if let Some(ref path) = config.dump_positions_path {
        let mut file = BufWriter::new(File::create(path)?);
        for target in &targets {
            for &pos in &target.minimizer_positions {
                writeln!(file, "{}\t{}", target.name, pos)?;
            }
        }
        file.flush()?;
        if !config.quiet {
            let total_positions: usize = targets.iter().map(|t| t.minimizer_positions.len()).sum();
            eprintln!(
                "Dumped {} syncmer positions to {}",
                total_positions,
                path.display()
            );
        }
    }

    // Process each sample in parallel
    use rayon::prelude::*;

    let is_multisample = config.sample_paths.len() > 1;
    let completed = if is_multisample && !config.quiet {
        // Give us a blank line to overwrite
        eprint!(
            "\x1B[2K\rSamples: processed 0 of {}…",
            config.sample_paths.len()
        );
        Some(Arc::new(Mutex::new(0usize)))
    } else {
        // Give us a blank line to overwrite
        if !config.quiet {
            eprint!("\x1B[2K\r");
        }
        None
    };

    let mut sample_results_with_idx: Vec<(usize, SampleResults)> = config
        .sample_paths
        .par_iter()
        .zip(&config.sample_names)
        .enumerate()
        .map(|(idx, (sample_paths, sample_name))| {
            let result = process_single_sample(
                idx,
                sample_paths, // Now a &Vec<PathBuf>
                sample_name,
                &targets,
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
                    config.sample_paths.len()
                );
            }

            result.map(|r| (idx, r))
        })
        .collect::<Result<Vec<_>>>()?;

    // Sort by original index to maintain CLI argument order
    sample_results_with_idx.sort_by_key(|(idx, _)| *idx);
    let sample_results: Vec<SampleResults> = sample_results_with_idx
        .into_iter()
        .map(|(_, result)| result)
        .collect();

    if is_multisample && !config.quiet {
        eprintln!();
    }

    let total_seqs_time: f64 = sample_results
        .iter()
        .map(|s| s.timing.seqs_processing_time)
        .sum();
    let total_analysis_time: f64 = sample_results.iter().map(|s| s.timing.analysis_time).sum();

    let total_time = start_time.elapsed();

    // Create report
    let report = Report {
        version: format!("skope {}", version),
        targets_file: config.targets_path.to_string_lossy().to_string(),
        parameters: ContainmentParameters {
            kmer_length: config.kmer_length,
            smer_length: config.smer_length,
            threads: config.threads,
            abundance_thresholds: config.abundance_thresholds.clone(),
        },
        samples: sample_results,
        total_timing: TimingStats {
            reference_processing_time: targets_time.as_secs_f64(),
            seqs_processing_time: total_seqs_time,
            analysis_time: total_analysis_time,
            total_time: total_time.as_secs_f64(),
            seqs_per_second: 0.0, // Not meaningful across samples
            bp_per_second: 0.0,    // Not meaningful across samples
        },
    };

    // Output results
    output_results(
        &report,
        config.output_path.as_ref(),
        config.output_format,
        config.sort_order,
        config.no_total,
    )?;

    Ok(())
}

/// Sort containment results based on the specified sort order (per-sample sorting)
fn sort_results(results: &mut [ContainmentResult], sort_order: SortOrder) {
    match sort_order {
        SortOrder::Original => {
            // No sorting needed - keep original order
        }
        SortOrder::Target => {
            results.sort_by(|a, b| a.target.cmp(&b.target));
        }
        SortOrder::Sample => {
            // For per-sample sorting (CSV/JSON), this doesn't make sense
            // but we'll keep it for consistency - sorts by target name
            results.sort_by(|a, b| a.target.cmp(&b.target));
        }
        SortOrder::Containment => {
            // Sort by containment descending (highest first)
            results.sort_by(|a, b| {
                b.containment1
                    .partial_cmp(&a.containment1)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        }
    }
}

fn output_results(
    report: &Report,
    output_path: Option<&PathBuf>,
    output_format: OutputFormat,
    sort_order: SortOrder,
    no_total: bool,
) -> Result<()> {
    let writer: Box<dyn Write> = if let Some(path) = output_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };

    let mut writer = writer;

    match output_format {
        OutputFormat::Tsv => {
            // TSV: per-sample sorting (keep existing approach)
            let mut sorted_report = report.clone();
            if sort_order != SortOrder::Original {
                for sample in &mut sorted_report.samples {
                    sort_results(&mut sample.targets, sort_order);
                }
            }
            output_tsv(&mut writer, &sorted_report, no_total)?;
        }
        OutputFormat::Table => {
            // Table: full cross-sample sorting
            output_table_sorted(&mut writer, report, sort_order, no_total)?;
        }
    }

    Ok(())
}

/// Flattened row for cross-sample sorting in table output
#[derive(Clone)]
struct TableRow<'a> {
    target: &'a str,
    sample_name: &'a str,
    sample_display: String,
    result: &'a ContainmentResult,
}

/// Output table with cross-sample sorting
fn output_table_sorted(
    writer: &mut dyn Write,
    report: &Report,
    sort_order: SortOrder,
    no_total: bool,
) -> Result<()> {
    let mut thresholds = report.parameters.abundance_thresholds.clone();
    thresholds.sort_unstable();

    // Build header with target column first
    let mut header = format!(
        "{:<50} | {:<20} | {:>13}",
        "target", "sample", "containment1"
    );
    for threshold in &thresholds {
        header.push_str(&format!(" | {:>15}", format!("containment{}", threshold)));
    }
    header.push_str(&format!(" | {:>18}", "median_nz_abundance"));

    let separator = "-".repeat(header.len());

    // Output header
    writeln!(writer)?;
    writeln!(writer, "{}", header)?;
    writeln!(writer, "{}", separator)?;

    // Flatten all rows (including TOTAL rows)
    let mut rows: Vec<TableRow> = Vec::new();

    for sample in &report.samples {
        let sample_display = if sample.seq_files.len() > 1 {
            format!(
                "{} ({} files)",
                sample.sample_name,
                sample.seq_files.len()
            )
        } else {
            sample.sample_name.clone()
        };

        // Add regular target rows
        for result in &sample.targets {
            rows.push(TableRow {
                target: &result.target,
                sample_name: &sample.sample_name,
                sample_display: sample_display.clone(),
                result,
            });
        }
    }

    // Sort the rows based on sort_order
    match sort_order {
        SortOrder::Original => {
            // Keep original order (already in order from nested loops)
        }
        SortOrder::Target => {
            // Sort by target name (alphabetically), then by sample order
            rows.sort_by(|a, b| {
                a.target
                    .cmp(b.target)
                    .then_with(|| a.sample_name.cmp(b.sample_name))
            });
        }
        SortOrder::Sample => {
            // Sort by sample name (alphabetically), then by target order
            rows.sort_by(|a, b| {
                a.sample_name
                    .cmp(b.sample_name)
                    .then_with(|| a.target.cmp(b.target))
            });
        }
        SortOrder::Containment => {
            // Sort by containment1 (descending)
            rows.sort_by(|a, b| {
                b.result
                    .containment1
                    .partial_cmp(&a.result.containment1)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        }
    }

    // Output sorted rows
    for row in &rows {
        let target_with_info = format!("{} ({})", row.result.target, format_bp(row.result.length));

        let mut output_row = format!(
            "{:<50} | {:<20} | {:>12.2}%",
            truncate_string(&target_with_info, 50),
            truncate_string(&row.sample_display, 20),
            row.result.containment1 * 100.0
        );
        for threshold in &thresholds {
            let containment = row
                .result
                .containment_at_threshold
                .get(threshold)
                .unwrap_or(&0.0);
            output_row.push_str(&format!(" | {:>14.2}%", containment * 100.0));
        }
        output_row.push_str(&format!(" | {:>18.0}", row.result.median_nz_abundance,));
        writeln!(writer, "{}", output_row)?;
    }

    // Output TOTAL rows (grouped by sample, after regular rows)
    if !no_total {
        for sample in &report.samples {
            let sample_display = if sample.seq_files.len() > 1 {
                format!(
                    "{} ({} files)",
                    sample.sample_name,
                    sample.seq_files.len()
                )
            } else {
                sample.sample_name.clone()
            };

            let mut total_row = format!(
                "{:<50} | {:<20} | {:>12.2}%",
                "TOTAL",
                truncate_string(&sample_display, 20),
                sample.total_stats.total_containment1 * 100.0
            );
            for threshold in &thresholds {
                let containment = sample
                    .total_stats
                    .total_containment_at_threshold
                    .get(threshold)
                    .unwrap_or(&0.0);
                total_row.push_str(&format!(" | {:>14.2}%", containment * 100.0));
            }
            total_row.push_str(&format!(
                " | {:>18.0}",
                sample.total_stats.total_median_nz_abundance,
            ));
            writeln!(writer, "{}", total_row)?;
        }
    }

    Ok(())
}

fn output_tsv(writer: &mut dyn Write, report: &Report, no_total: bool) -> Result<()> {
    let mut thresholds = report.parameters.abundance_thresholds.clone();
    thresholds.sort_unstable();

    // Build header with target column first
    let mut header = "target\tsample\tcontainment1\tcontainment1_hits".to_string();
    for threshold in &thresholds {
        header.push_str(&format!(
            "\tcontainment{}\tcontainment{}_hits",
            threshold, threshold
        ));
    }
    header.push_str("\ttarget_length\ttarget_kmers\tmedian_nz_abundance");
    writeln!(writer, "{}", header)?;

    // Output data rows for all samples
    for sample in &report.samples {
        for result in &sample.targets {
            let mut row = format!(
                "{}\t{}\t{:.5}\t{}",
                result.target,
                sample.sample_name,
                result.containment1,
                result.contained_minimizers // Use contained_minimizers for threshold 1 hits
            );
            for threshold in &thresholds {
                let containment = result
                    .containment_at_threshold
                    .get(threshold)
                    .unwrap_or(&0.0);
                let hits = result.hits_at_threshold.get(threshold).unwrap_or(&0);
                row.push_str(&format!("\t{:.5}\t{}", containment, hits));
            }
            row.push_str(&format!(
                "\t{}\t{}\t{:.0}",
                result.length, result.target_kmers, result.median_nz_abundance,
            ));
            writeln!(writer, "{}", row)?;
        }

        // Add TOTAL row for this sample
        if !no_total {
            let mut total_row = format!(
                "{}\t{}\t{:.5}\t{}",
                "TOTAL",
                sample.sample_name,
                sample.total_stats.total_containment1,
                sample.total_stats.total_contained_minimizers
            );
            for threshold in &thresholds {
                let containment = sample
                    .total_stats
                    .total_containment_at_threshold
                    .get(threshold)
                    .unwrap_or(&0.0);
                let hits = (containment * sample.total_stats.target_kmers as f64).round() as usize;
                total_row.push_str(&format!("\t{:.5}\t{}", containment, hits));
            }
            total_row.push_str(&format!(
                "\t0\t{}\t{:.0}",
                sample.total_stats.target_kmers, sample.total_stats.total_median_nz_abundance,
            ));
            writeln!(writer, "{}", total_row)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_truncate_string() {
        assert_eq!(truncate_string("short", 10), "short");
        assert_eq!(truncate_string("verylongstring", 8), "veryl...");
    }
}

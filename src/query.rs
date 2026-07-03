use crate::stats::{WILSON_Z_95, normal_survival, wilson_interval};
use crate::syncmers::{
    Buffers, KmerHasher, SyncmerVec, decode_u64, decode_u128, fill_syncmers,
    fill_syncmers_with_positions,
};
use crate::{
    ProcessingStats, RapidHashSet, create_spinner, discover_sequence_groups, format_bp,
    format_bp_per_sec, handle_process_result, reader_with_inferred_batch_size,
    sample_limit_reached_io_error,
};
use anyhow::{Context, Result};
use indicatif::ProgressBar;
use paraseq::Record;
use paraseq::parallel::{ParallelProcessor, ParallelReader};
use parking_lot::Mutex;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, BufWriter, Read, Write};
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

/// Zero-cost (hopefully?) abstraction over u64 and u128 syncmer sets
#[derive(Debug, Clone)]
enum SyncmerSet {
    U64(RapidHashSet<u64>),
    U128(RapidHashSet<u128>),
}

impl SyncmerSet {
    fn len(&self) -> usize {
        match self {
            SyncmerSet::U64(set) => set.len(),
            SyncmerSet::U128(set) => set.len(),
        }
    }

    fn extend(&mut self, other: &SyncmerSet) {
        match (self, other) {
            (SyncmerSet::U64(a), SyncmerSet::U64(b)) => a.extend(b.iter().copied()),
            (SyncmerSet::U128(a), SyncmerSet::U128(b)) => a.extend(b.iter().copied()),
            _ => panic!("Cannot extend SyncmerSet: mismatched variants"),
        }
    }

    /// Remove all syncmers from --background stream, returning count
    fn retain_not_in(&mut self, other: &SyncmerSet) -> usize {
        match (self, other) {
            (SyncmerSet::U64(a), SyncmerSet::U64(b)) => {
                let before = a.len();
                a.retain(|s| !b.contains(s));
                before - a.len()
            }
            (SyncmerSet::U128(a), SyncmerSet::U128(b)) => {
                let before = a.len();
                a.retain(|s| !b.contains(s));
                before - a.len()
            }
            _ => panic!("Cannot mask SyncmerSet: mismatched variants"),
        }
    }
}

#[derive(Debug, Clone)]
enum PositionedSyncmers {
    U64(Vec<(u64, usize)>),
    U128(Vec<(u128, usize)>),
    Empty,
}

impl PositionedSyncmers {
    fn offset_positions(&mut self, offset: usize) {
        match self {
            PositionedSyncmers::U64(entries) => {
                for (_, position) in entries {
                    *position += offset;
                }
            }
            PositionedSyncmers::U128(entries) => {
                for (_, position) in entries {
                    *position += offset;
                }
            }
            PositionedSyncmers::Empty => {}
        }
    }

    fn extend(&mut self, other: PositionedSyncmers) {
        match self {
            PositionedSyncmers::U64(a) => match other {
                PositionedSyncmers::U64(b) => a.extend(b),
                PositionedSyncmers::Empty => {}
                _ => panic!("Cannot extend PositionedSyncmers: mismatched variants"),
            },
            PositionedSyncmers::U128(a) => match other {
                PositionedSyncmers::U128(b) => a.extend(b),
                PositionedSyncmers::Empty => {}
                _ => panic!("Cannot extend PositionedSyncmers: mismatched variants"),
            },
            PositionedSyncmers::Empty => *self = other,
        }
    }

    fn retain_in_set(&mut self, syncmers: &SyncmerSet) {
        match (self, syncmers) {
            (PositionedSyncmers::U64(entries), SyncmerSet::U64(set)) => {
                entries.retain(|(syncmer, _)| set.contains(syncmer));
            }
            (PositionedSyncmers::U128(entries), SyncmerSet::U128(set)) => {
                entries.retain(|(syncmer, _)| set.contains(syncmer));
            }
            (PositionedSyncmers::Empty, _) => {}
            _ => panic!("Cannot filter PositionedSyncmers: mismatched variants"),
        }
    }

    fn positions(&self) -> Vec<usize> {
        match self {
            PositionedSyncmers::U64(entries) => {
                entries.iter().map(|(_, position)| *position).collect()
            }
            PositionedSyncmers::U128(entries) => {
                entries.iter().map(|(_, position)| *position).collect()
            }
            PositionedSyncmers::Empty => Vec::new(),
        }
    }

    fn len(&self) -> usize {
        match self {
            PositionedSyncmers::U64(entries) => entries.len(),
            PositionedSyncmers::U128(entries) => entries.len(),
            PositionedSyncmers::Empty => 0,
        }
    }
}

#[derive(Debug, Clone)]
struct TargetInfo {
    name: String,
    length: usize,
    syncmers: SyncmerSet,
    syncmer_positions: Vec<usize>,
    positioned_syncmers: PositionedSyncmers,
}

#[derive(Debug, Clone, Copy)]
pub struct PatchinessResult {
    pub z: f64,
    pub p: f64,
}

#[derive(Debug, Clone)]
pub struct ContainmentResult {
    pub target: String,
    pub length: usize,
    pub target_kmers: usize,
    pub containment1_hits: usize,
    pub containment1: f64,
    pub ani_est: Option<f64>,
    pub median_nz_abundance: f64,
    pub containment_at_threshold: HashMap<usize, f64>, // threshold -> containment
    pub hits_at_threshold: HashMap<usize, usize>,      // threshold -> hit count
    pub patchiness: Option<PatchinessResult>,
}

#[derive(Debug, Clone)]
pub struct ContainmentParameters {
    pub kmer_length: u8,
    pub smer_length: u8,
    pub threads: usize,
    pub abundance_thresholds: Vec<usize>,
    pub confidence: bool,
}

#[derive(Debug, Clone)]
pub struct TotalStats {
    pub total_targets: usize,
    pub target_kmers: usize,
    pub total_containment1_hits: usize,
    pub total_containment1: f64,
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

pub struct ContainmentConfig {
    pub targets_path: PathBuf,
    pub background_paths: Vec<PathBuf>, // Off-target sequences to mask (empty = none)
    pub sample_paths: Vec<Vec<PathBuf>>, // Each sample is a Vec of file paths
    pub sample_names: Vec<String>,
    pub kmer_length: u8,
    pub smer_length: u8,
    pub threads: usize,
    pub output_path: Option<PathBuf>,
    pub quiet: bool,
    pub abundance_thresholds: Vec<usize>,
    pub discriminatory: bool,
    pub positions: bool,
    pub spacing: u16,
    pub individual: bool,
    pub limit_bp: Option<u64>,
    pub sort_order: SortOrder,
    pub dump_syncmers_path: Option<PathBuf>,
    pub no_total: bool,
    pub confidence: bool,
}

impl ContainmentConfig {
    pub fn execute(&self) -> Result<()> {
        run_containment_analysis(self)
    }
}

fn normalize_abundance_thresholds(thresholds: &[usize]) -> Vec<usize> {
    let mut thresholds: Vec<usize> = thresholds
        .iter()
        .copied()
        .filter(|&threshold| threshold > 1)
        .collect();
    thresholds.sort_unstable();
    thresholds.dedup();
    thresholds
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
    spacing: u16,
    collect_positions: bool,

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
        targets: Arc<Mutex<Vec<TargetInfo>>>,
        global_stats: Arc<Mutex<ProcessingStats>>,
        spinner: Option<Arc<Mutex<ProgressBar>>>,
        start_time: std::time::Instant,
        spacing: u16,
        collect_positions: bool,
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
            targets,
            spacing,
            collect_positions,
            local_stats: ProcessingStats::default(),
            global_stats,
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

        if self.collect_positions {
            fill_syncmers_with_positions(
                &sequence,
                &self.hasher,
                self.kmer_length,
                self.smer_length,
                &mut self.buffers,
                &mut self.positions,
                self.spacing,
            );
        } else {
            fill_syncmers(
                &sequence,
                &self.hasher,
                self.kmer_length,
                self.smer_length,
                &mut self.buffers,
                self.spacing,
            );
        }

        // Build unique syncmer set for this target
        let syncmers = match &self.buffers.syncmers {
            SyncmerVec::U64(vec) => {
                let set: RapidHashSet<u64> = vec.iter().copied().collect();
                SyncmerSet::U64(set)
            }
            SyncmerVec::U128(vec) => {
                let set: RapidHashSet<u128> = vec.iter().copied().collect();
                SyncmerSet::U128(set)
            }
        };

        let positioned_syncmers = if self.collect_positions {
            match &self.buffers.syncmers {
                SyncmerVec::U64(vec) => PositionedSyncmers::U64(
                    vec.iter()
                        .copied()
                        .zip(self.positions.iter().copied())
                        .collect(),
                ),
                SyncmerVec::U128(vec) => PositionedSyncmers::U128(
                    vec.iter()
                        .copied()
                        .zip(self.positions.iter().copied())
                        .collect(),
                ),
            }
        } else {
            PositionedSyncmers::Empty
        };

        let syncmer_positions = if self.collect_positions {
            self.positions.clone()
        } else {
            Vec::new()
        };

        self.targets.lock().push(TargetInfo {
            name: target_name,
            length: sequence.len(),
            syncmers,
            syncmer_positions,
            positioned_syncmers,
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

fn process_targets_file(
    targets_path: &Path,
    kmer_length: u8,
    smer_length: u8,
    quiet: bool,
    spacing: u16,
    collect_positions: bool,
) -> Result<Vec<TargetInfo>> {
    let in_path = if targets_path.to_string_lossy() == "-" {
        None
    } else {
        Some(targets_path)
    };

    let reader = reader_with_inferred_batch_size(in_path)?;

    let spinner = create_spinner(quiet)?;

    let start_time = std::time::Instant::now();
    let targets: Arc<Mutex<Vec<TargetInfo>>> = Arc::new(Mutex::new(Vec::new()));
    let global_stats = Arc::new(Mutex::new(ProcessingStats::default()));

    let mut processor = TargetsProcessor::new(
        kmer_length,
        smer_length,
        Arc::clone(&targets),
        Arc::clone(&global_stats),
        spinner.clone(),
        start_time,
        spacing,
        collect_positions,
    );

    // Single thread to preserve order
    reader.process_parallel(&mut processor, 1)?;

    // Finish spinner and clear
    if let Some(ref pb) = spinner {
        pb.lock().finish_and_clear();
    }

    // Drop processor so its Arc clones are released
    drop(processor);
    let targets = Arc::try_unwrap(targets).unwrap().into_inner();

    Ok(targets)
}

/// Merge multiple TargetInfos into one with the given name
fn merge_targets(targets: Vec<TargetInfo>, name: String) -> Result<TargetInfo> {
    if targets.is_empty() {
        return Err(anyhow::anyhow!("No targets to merge for '{}'", name));
    }

    let mut iter = targets.into_iter();
    let mut merged = iter.next().unwrap();
    merged.name = name;

    for t in iter {
        let offset = merged.length;
        let mut positioned_syncmers = t.positioned_syncmers;
        positioned_syncmers.offset_positions(offset);
        merged.syncmers.extend(&t.syncmers);
        merged.syncmer_positions.extend(
            t.syncmer_positions
                .into_iter()
                .map(|position| position + offset),
        );
        merged.positioned_syncmers.extend(positioned_syncmers);
        merged.length += t.length;
    }

    Ok(merged)
}

/// Process a directory of fastx files/subdirectories as targets
/// Each top-level fastx file becomes one target (all records merged).
/// Each subdirectory becomes one target (all fastx files within merged).
fn process_targets_dir(
    dir_path: &Path,
    kmer_length: u8,
    smer_length: u8,
    quiet: bool,
    spacing: u16,
    collect_positions: bool,
) -> Result<Vec<TargetInfo>> {
    let groups = discover_sequence_groups(dir_path)?;

    let mut results: Vec<TargetInfo> = Vec::with_capacity(groups.len());
    for group in groups {
        let mut all_targets: Vec<TargetInfo> = Vec::new();
        for file_path in &group.files {
            let targets = process_targets_file(
                file_path,
                kmer_length,
                smer_length,
                quiet,
                spacing,
                collect_positions,
            )?;
            all_targets.extend(targets);
        }
        results.push(merge_targets(all_targets, group.name)?);
    }

    Ok(results)
}

/// Processor for counting syncmer depths from sequences
#[derive(Clone)]
struct SeqsProcessor {
    kmer_length: u8,
    smer_length: u8,
    hasher: KmerHasher,
    targets_syncmers: Arc<SyncmerSet>,

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
    spinner_label: &'static str,
    start_time: Instant,
    limit_bp: Option<u64>,
}

impl SeqsProcessor {
    #[allow(clippy::too_many_arguments)]
    fn new(
        kmer_length: u8,
        smer_length: u8,
        targets_syncmers: Arc<SyncmerSet>,
        global_counts_u64: Arc<Mutex<Option<HashMap<u64, CountDepth>>>>,
        global_counts_u128: Arc<Mutex<Option<HashMap<u128, CountDepth>>>>,
        global_stats: Arc<Mutex<ProcessingStats>>,
        spinner: Option<Arc<Mutex<ProgressBar>>>,
        spinner_label: &'static str,
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

        Self {
            kmer_length,
            smer_length,
            hasher: KmerHasher::new(smer_length as usize),
            targets_syncmers,
            buffers,
            local_stats: ProcessingStats::default(),
            local_counts_u64,
            local_counts_u128,
            global_stats,
            global_counts_u64,
            global_counts_u128,
            spinner,
            spinner_label,
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
                "{}: {} seqs ({}). {:.0} seqs/s ({})",
                self.spinner_label,
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
        if let Some(limit) = self.limit_bp {
            let global_bp = self.global_stats.lock().total_bp;
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
            1, // samples are always dense; thinning the sample side biases containment
        );

        // Count syncmers present in targets
        match (&self.buffers.syncmers, &*self.targets_syncmers) {
            (SyncmerVec::U64(vec), SyncmerSet::U64(targets_set)) => {
                let local_counts = self.local_counts_u64.as_mut().unwrap();
                for &syncmer in vec {
                    if targets_set.contains(&syncmer) {
                        local_counts
                            .entry(syncmer)
                            .and_modify(|e| *e = e.saturating_add(1))
                            .or_insert(1);
                    }
                }
            }
            (SyncmerVec::U128(vec), SyncmerSet::U128(targets_set)) => {
                let local_counts = self.local_counts_u128.as_mut().unwrap();
                for &syncmer in vec {
                    if targets_set.contains(&syncmer) {
                        local_counts
                            .entry(syncmer)
                            .and_modify(|e| *e = e.saturating_add(1))
                            .or_insert(1);
                    }
                }
            }
            _ => panic!("Mismatch between SyncmerVec and SyncmerSet types"),
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Merge local into global counts
        if let Some(local) = &mut self.local_counts_u64 {
            let mut global = self.global_counts_u64.lock();
            let global_map = global.as_mut().unwrap();
            for (&syncmer, &count) in local.iter() {
                global_map
                    .entry(syncmer)
                    .and_modify(|e| *e = e.saturating_add(count))
                    .or_insert(count);
            }
            local.clear();
        } else {
            let mut global = self.global_counts_u128.lock();
            let global_map = global.as_mut().unwrap();
            let local = self.local_counts_u128.as_mut().unwrap();
            for (&syncmer, &count) in local.iter() {
                global_map
                    .entry(syncmer)
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

const MIN_TARGET_KMERS_FOR_ANI_ADJUSTMENT: usize = 50;
const MIN_NONZERO_KMERS_FOR_ANI_ADJUSTMENT: usize = 25;
const MIN_DEPTH_BIN_FOR_ANI_ADJUSTMENT: usize = 3;
const MAX_MEDIAN_DEPTH_FOR_ANI_ADJUSTMENT: f64 = 2.0;
// Gates for displaying an ANI estimate: minimum target syncmers and minimum ANI.
const MIN_TARGET_KMERS_FOR_ANI_DISPLAY: usize = 50;
const MIN_ANI_FOR_DISPLAY: f64 = 0.90;
const MIN_TARGET_KMERS_FOR_PATCHINESS: usize = 50;
const MIN_HITS_FOR_PATCHINESS: usize = 10;
const MIN_MISSES_FOR_PATCHINESS: usize = 10;

fn estimate_lambda(
    abundance_histogram: &[(CountDepth, usize)],
    target_kmers: usize,
    median_nz_abundance: f64,
) -> Option<f64> {
    if target_kmers < MIN_TARGET_KMERS_FOR_ANI_ADJUSTMENT
        || median_nz_abundance > MAX_MEDIAN_DEPTH_FOR_ANI_ADJUSTMENT
    {
        return None;
    }

    let nonzero_count: usize = abundance_histogram
        .iter()
        .filter(|(abundance, _)| *abundance > 0)
        .map(|(_, count)| *count)
        .sum();
    if nonzero_count < MIN_NONZERO_KMERS_FOR_ANI_ADJUSTMENT {
        return None;
    }

    let (mode_abundance, mode_count) = abundance_histogram
        .iter()
        .filter(|(abundance, _)| *abundance > 0)
        .max_by_key(|(abundance, count)| (*count, *abundance))?;

    let next_abundance = mode_abundance.checked_add(1)?;
    let next_count = abundance_histogram
        .iter()
        .find(|(abundance, _)| *abundance == next_abundance)
        .map(|(_, count)| *count)
        .unwrap_or(0);

    if *mode_count < MIN_DEPTH_BIN_FOR_ANI_ADJUSTMENT
        || next_count < MIN_DEPTH_BIN_FOR_ANI_ADJUSTMENT
    {
        return None;
    }

    let lambda = next_count as f64 / *mode_count as f64 * next_abundance as f64;
    lambda
        .is_finite()
        .then_some(lambda)
        .filter(|lambda| *lambda > 0.0)
}

fn poisson_tail_at_least(lambda: f64, threshold: usize) -> Option<f64> {
    if !lambda.is_finite() || lambda <= 0.0 || threshold == 0 {
        return None;
    }

    if threshold == 1 {
        return Some(1.0 - (-lambda).exp());
    }

    let mut term = (-lambda).exp();
    let mut cdf = term;
    for k in 1..threshold {
        term *= lambda / k as f64;
        cdf += term;
    }

    let tail = 1.0 - cdf;
    tail.is_finite().then_some(tail.clamp(0.0, 1.0))
}

fn coverage_adjusted_kmer_identity(containment: f64, lambda: f64, threshold: usize) -> Option<f64> {
    let detection_probability = poisson_tail_at_least(lambda, threshold)?;
    if detection_probability <= 0.0 {
        return None;
    }
    Some((containment / detection_probability).clamp(0.0, 1.0))
}

fn ani_estimate(
    containment: f64,
    kmer_length: u8,
    lambda: Option<f64>,
    target_kmers: usize,
) -> Option<f64> {
    // Need enough syncmers and at least one contained syncmer.
    if kmer_length == 0 || target_kmers < MIN_TARGET_KMERS_FOR_ANI_DISPLAY || containment <= 0.0 {
        return None;
    }

    let kmer_identity = lambda
        .and_then(|lambda| coverage_adjusted_kmer_identity(containment, lambda, 1))
        .unwrap_or_else(|| containment.clamp(0.0, 1.0));

    let ani = kmer_identity.powf(1.0 / kmer_length as f64);

    // Suppress low-confidence estimates.
    (ani >= MIN_ANI_FOR_DISPLAY).then_some(ani)
}

fn calculate_patchiness_from_labels(labels: &[bool]) -> Option<PatchinessResult> {
    let n = labels.len();
    if n < MIN_TARGET_KMERS_FOR_PATCHINESS {
        return None;
    }

    let hits = labels.iter().filter(|&&hit| hit).count();
    let misses = n - hits;
    if hits < MIN_HITS_FOR_PATCHINESS || misses < MIN_MISSES_FOR_PATCHINESS {
        return None;
    }

    let runs = labels
        .windows(2)
        .filter(|window| window[0] != window[1])
        .count()
        + 1;

    let n = n as f64;
    let hits = hits as f64;
    let misses = misses as f64;
    let expected = 1.0 + (2.0 * hits * misses) / n;
    let variance = (2.0 * hits * misses * (2.0 * hits * misses - n)) / (n * n * (n - 1.0));

    if !variance.is_finite() || variance <= 0.0 {
        return None;
    }

    let runs_z = (runs as f64 - expected) / variance.sqrt();
    let patchiness_z = -runs_z;
    Some(PatchinessResult {
        z: patchiness_z,
        p: normal_survival(patchiness_z),
    })
}

fn calculate_patchiness_u64(
    positioned_syncmers: &[(u64, usize)],
    abundances: &HashMap<u64, CountDepth>,
    threshold: usize,
) -> Option<PatchinessResult> {
    if threshold == 0 {
        return None;
    }

    let mut occurrence_counts: HashMap<u64, usize> = HashMap::new();
    for &(syncmer, _) in positioned_syncmers {
        *occurrence_counts.entry(syncmer).or_insert(0) += 1;
    }

    let labels: Vec<bool> = positioned_syncmers
        .iter()
        .filter(|(syncmer, _)| {
            occurrence_counts
                .get(syncmer)
                .is_some_and(|&count| count == 1)
        })
        .map(|(syncmer, _)| {
            abundances.get(syncmer).copied().unwrap_or(0) >= threshold as CountDepth
        })
        .collect();

    calculate_patchiness_from_labels(&labels)
}

fn calculate_patchiness_u128(
    positioned_syncmers: &[(u128, usize)],
    abundances: &HashMap<u128, CountDepth>,
    threshold: usize,
) -> Option<PatchinessResult> {
    if threshold == 0 {
        return None;
    }

    let mut occurrence_counts: HashMap<u128, usize> = HashMap::new();
    for &(syncmer, _) in positioned_syncmers {
        *occurrence_counts.entry(syncmer).or_insert(0) += 1;
    }

    let labels: Vec<bool> = positioned_syncmers
        .iter()
        .filter(|(syncmer, _)| {
            occurrence_counts
                .get(syncmer)
                .is_some_and(|&count| count == 1)
        })
        .map(|(syncmer, _)| {
            abundances.get(syncmer).copied().unwrap_or(0) >= threshold as CountDepth
        })
        .collect();

    calculate_patchiness_from_labels(&labels)
}

fn calculate_patchiness(
    positioned_syncmers: &PositionedSyncmers,
    abundance_map: &AbundanceMap,
    threshold: usize,
) -> Option<PatchinessResult> {
    match (positioned_syncmers, abundance_map) {
        (PositionedSyncmers::U64(positioned), AbundanceMap::U64(map)) => {
            calculate_patchiness_u64(positioned, map, threshold)
        }
        (PositionedSyncmers::U128(positioned), AbundanceMap::U128(map)) => {
            calculate_patchiness_u128(positioned, map, threshold)
        }
        (PositionedSyncmers::Empty, _) => None,
        _ => panic!("Mismatch between PositionedSyncmers and AbundanceMap types"),
    }
}

fn process_seqs_file(
    seq_path: &Path,
    targets_syncmers: Arc<SyncmerSet>,
    kmer_length: u8,
    smer_length: u8,
    threads: usize,
    quiet: bool,
    limit_bp: Option<u64>,
    label: &'static str,
) -> Result<(AbundanceMap, u64, u64)> {
    let in_path = if seq_path.to_string_lossy() == "-" {
        None
    } else {
        Some(seq_path)
    };
    let reader = reader_with_inferred_batch_size(in_path)?;

    let spinner = create_spinner(quiet)?;
    if let Some(ref pb) = spinner {
        pb.lock().set_message(format!("{label}: 0 seqs (0bp)"));
    }

    let total_target_syncmers = targets_syncmers.len();

    let start_time = Instant::now();
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
    let global_stats = Arc::new(Mutex::new(ProcessingStats::default()));

    let mut processor = SeqsProcessor::new(
        kmer_length,
        smer_length,
        targets_syncmers,
        Arc::clone(&global_counts_u64),
        Arc::clone(&global_counts_u128),
        Arc::clone(&global_stats),
        spinner.clone(),
        label,
        start_time,
        limit_bp,
    );

    let process_result = reader.process_parallel(&mut processor, threads);
    handle_process_result(process_result)?;

    if let Some(ref pb) = spinner {
        pb.lock().finish_with_message("");
    }

    // Drop processor so its Arc clones are released
    drop(processor);
    let stats = global_stats.lock().clone();
    let abundance_map = if kmer_length <= 32 {
        let map = Arc::try_unwrap(global_counts_u64)
            .unwrap()
            .into_inner()
            .unwrap();
        AbundanceMap::U64(map)
    } else {
        let map = Arc::try_unwrap(global_counts_u128)
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
    kmer_length: u8,
    calculate_patchiness_metrics: bool,
) -> Vec<ContainmentResult> {
    targets
        .iter()
        .map(|target| {
            let target_kmers = target.syncmers.len();
            let mut abundances: Vec<CountDepth> = Vec::new();
            let mut contained_count = 0;

            // Collect abundances for all unique syncmers in this target
            match (&target.syncmers, abundance_map) {
                (SyncmerSet::U64(set), AbundanceMap::U64(map)) => {
                    for &syncmer in set {
                        let abundance = map.get(&syncmer).copied().unwrap_or(0);
                        abundances.push(abundance);
                        if abundance > 0 {
                            contained_count += 1;
                        }
                    }
                }
                (SyncmerSet::U128(set), AbundanceMap::U128(map)) => {
                    for &syncmer in set {
                        let abundance = map.get(&syncmer).copied().unwrap_or(0);
                        abundances.push(abundance);
                        if abundance > 0 {
                            contained_count += 1;
                        }
                    }
                }
                _ => panic!("Mismatch between SyncmerSet and AbundanceMap types"),
            }

            let containment1 = if target_kmers > 0 {
                contained_count as f64 / target_kmers as f64
            } else {
                0.0
            };

            // Ignore zero-abundance syncmers for median calc
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

            let lambda = estimate_lambda(&abundance_histogram, target_kmers, median_nz_abundance);
            let ani_est = ani_estimate(containment1, kmer_length, lambda, target_kmers);
            let patchiness = if calculate_patchiness_metrics {
                calculate_patchiness(&target.positioned_syncmers, abundance_map, 1)
            } else {
                None
            };

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
                containment1_hits: contained_count,
                containment1,
                ani_est,
                median_nz_abundance,
                containment_at_threshold,
                hits_at_threshold,
                patchiness,
            }
        })
        .collect()
}

/// Process a single sample's sequences and calculate statistics
fn process_single_sample(
    _idx: usize,
    sample_paths: &[PathBuf], // Multiple files per sample
    sample_name: &str,
    targets: &[TargetInfo],
    targets_syncmers: Arc<SyncmerSet>,
    abundance_thresholds: &[usize],
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
            Arc::clone(&targets_syncmers),
            config.kmer_length,
            config.smer_length,
            config.threads,
            quiet_sample,
            config.limit_bp.map(|limit| limit.saturating_sub(total_bp)),
            "Processing sample",
        )?;

        // Merge abundance maps
        match (&mut combined_abundance_map, file_abundance_map) {
            (AbundanceMap::U64(combined), AbundanceMap::U64(new)) => {
                for (syncmer, count) in new {
                    combined
                        .entry(syncmer)
                        .and_modify(|e| *e = e.saturating_add(count))
                        .or_insert(count);
                }
            }
            (AbundanceMap::U128(combined), AbundanceMap::U128(new)) => {
                for (syncmer, count) in new {
                    combined
                        .entry(syncmer)
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
            && total_bp >= limit
        {
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
        abundance_thresholds,
        config.kmer_length,
        config.confidence,
    );
    let analysis_time = analysis_start.elapsed();

    // Overall stats for this sample
    let target_kmers: usize = containment_results.iter().map(|r| r.target_kmers).sum();
    let total_containment1_hits: usize = containment_results
        .iter()
        .map(|r| r.containment1_hits)
        .sum();
    let total_containment1 = if target_kmers > 0 {
        total_containment1_hits as f64 / target_kmers as f64
    } else {
        0.0
    };

    let seqs_per_second = total_seqs as f64 / seqs_time.as_secs_f64();
    let bp_per_second = total_bp as f64 / seqs_time.as_secs_f64();

    // Calculate overall containment at each threshold
    let mut total_containment_at_threshold = HashMap::new();
    for &threshold in abundance_thresholds {
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
            total_containment1_hits,
            total_containment1,
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

/// Combined set of every target's syncmers
fn build_union(targets: &[TargetInfo], kmer_length: u8) -> SyncmerSet {
    if kmer_length <= 32 {
        let mut set = RapidHashSet::default();
        for t in targets {
            if let SyncmerSet::U64(s) = &t.syncmers {
                set.extend(s.iter());
            }
        }
        SyncmerSet::U64(set)
    } else {
        let mut set = RapidHashSet::default();
        for t in targets {
            if let SyncmerSet::U128(s) = &t.syncmers {
                set.extend(s.iter());
            }
        }
        SyncmerSet::U128(set)
    }
}

/// Drop target syncmers shared with background sequences. Streaming via `process_seqs_file`
/// keeps peak memory bounded by the targets, not the background. Returns count removed.
fn mask_background(
    targets: &mut [TargetInfo],
    union: &SyncmerSet,
    background_paths: &[PathBuf],
    kmer_length: u8,
    smer_length: u8,
    threads: usize,
    quiet: bool,
) -> Result<usize> {
    let union = Arc::new(union.clone());
    let mut rm_u64: RapidHashSet<u64> = RapidHashSet::default();
    let mut rm_u128: RapidHashSet<u128> = RapidHashSet::default();
    for path in background_paths {
        let (map, _, _) = process_seqs_file(
            path,
            Arc::clone(&union),
            kmer_length,
            smer_length,
            threads,
            quiet,
            None,
            "Masking background",
        )?;
        match map {
            AbundanceMap::U64(m) => rm_u64.extend(m.into_keys()),
            AbundanceMap::U128(m) => rm_u128.extend(m.into_keys()),
        }
    }

    let to_remove = if kmer_length <= 32 {
        SyncmerSet::U64(rm_u64)
    } else {
        SyncmerSet::U128(rm_u128)
    };
    let removed = to_remove.len();

    // Apply removal to each target; re-sync positions like the discriminatory filter
    for target in targets.iter_mut() {
        target.syncmers.retain_not_in(&to_remove);
        target.positioned_syncmers.retain_in_set(&target.syncmers);
        target.syncmer_positions = target.positioned_syncmers.positions();
    }
    Ok(removed)
}

// ── Query index (.sk, kind = query) ─────────────────────────────────────────
// Header + wincode metadata, then per target `count` raw entries of `ceil(k/4)` value
// bytes, each with a trailing u32 position when positions are stored.
use crate::{INDEX_MAGIC, IndexKind, read_index_kind};

const QUERY_INDEX_VERSION: u8 = 1;
const QUERY_INDEX_FLAG_POSITIONS: u8 = 0b001;
const QUERY_INDEX_FLAG_BACKGROUND: u8 = 0b010;

// magic, kind, version, k, s, spacing, flags, n_targets
type QueryIndexHeader = ([u8; 4], u8, u8, u8, u8, u16, u8, u32);
type QueryIndexMeta = Vec<(String, u64, u64)>; // (name, length, entry_count)

/// Config for `skope index build-query`
pub struct BuildQueryConfig {
    pub targets_path: PathBuf,
    pub background_paths: Vec<PathBuf>,
    pub kmer_length: u8,
    pub smer_length: u8,
    pub spacing: u16,
    pub individual: bool,
    pub positions: bool,
    pub threads: usize,
    pub output_path: Option<PathBuf>,
    pub quiet: bool,
}

struct QueryIndex {
    targets: Vec<TargetInfo>,
    has_positions: bool,
}

/// True if `path` is a skope query index
pub fn is_query_index(path: &Path) -> bool {
    read_index_kind(path) == Some(IndexKind::Query)
}

/// Read k, s, spacing from a query index header without loading entries
pub fn read_query_index_meta(path: &Path) -> Result<(u8, u8, u16)> {
    let mut buf = [0u8; 64];
    let n = File::open(path)?.read(&mut buf)?;
    let mut cursor = wincode::io::Cursor::new(&buf[..n]);
    let (magic, kind, version, k, s, spacing, _f, _n): QueryIndexHeader =
        wincode::deserialize_from(&mut cursor).context("Failed to read query index header")?;
    validate_query_index_header(&magic, kind, version)?;
    Ok((k, s, spacing))
}

fn validate_query_index_header(magic: &[u8; 4], kind: u8, version: u8) -> Result<()> {
    if magic != INDEX_MAGIC || IndexKind::from_byte(kind) != Some(IndexKind::Query) {
        return Err(anyhow::anyhow!("Not a skope query index"));
    }
    if version != QUERY_INDEX_VERSION {
        return Err(anyhow::anyhow!(
            "Unsupported query index version: {version} (expected {QUERY_INDEX_VERSION})"
        ));
    }
    Ok(())
}

fn save_query_index(
    targets: &[TargetInfo],
    kmer_length: u8,
    smer_length: u8,
    spacing: u16,
    flags: u8,
    output_path: Option<&Path>,
) -> Result<()> {
    let has_positions = flags & QUERY_INDEX_FLAG_POSITIONS != 0;
    let kmer_bytes = (kmer_length as usize).div_ceil(4);

    if has_positions && targets.iter().any(|t| t.length > u32::MAX as usize) {
        return Err(anyhow::anyhow!(
            "Target exceeds {} bp; positions cannot be stored (build without --positions)",
            u32::MAX
        ));
    }

    let mut writer: Box<dyn Write> = match output_path {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };

    let meta: QueryIndexMeta = targets
        .iter()
        .map(|t| {
            let count = if has_positions {
                t.positioned_syncmers.len()
            } else {
                t.syncmers.len()
            } as u64;
            (t.name.clone(), t.length as u64, count)
        })
        .collect();

    let header: QueryIndexHeader = (
        *INDEX_MAGIC,
        IndexKind::Query as u8,
        QUERY_INDEX_VERSION,
        kmer_length,
        smer_length,
        spacing,
        flags,
        targets.len() as u32,
    );
    writer
        .write_all(&wincode::serialize(&header).context("Failed to encode query index header")?)?;
    writer
        .write_all(&wincode::serialize(&meta).context("Failed to encode query index metadata")?)?;

    for t in targets {
        if has_positions {
            match &t.positioned_syncmers {
                PositionedSyncmers::U64(entries) => {
                    for &(v, pos) in entries {
                        writer.write_all(&v.to_le_bytes()[..kmer_bytes])?;
                        writer.write_all(&(pos as u32).to_le_bytes())?;
                    }
                }
                PositionedSyncmers::U128(entries) => {
                    for &(v, pos) in entries {
                        writer.write_all(&v.to_le_bytes()[..kmer_bytes])?;
                        writer.write_all(&(pos as u32).to_le_bytes())?;
                    }
                }
                PositionedSyncmers::Empty => {}
            }
        } else {
            match &t.syncmers {
                SyncmerSet::U64(set) => {
                    for &v in set {
                        writer.write_all(&v.to_le_bytes()[..kmer_bytes])?;
                    }
                }
                SyncmerSet::U128(set) => {
                    for &v in set {
                        writer.write_all(&v.to_le_bytes()[..kmer_bytes])?;
                    }
                }
            }
        }
    }
    writer.flush()?;
    Ok(())
}

fn load_query_index(path: &Path) -> Result<QueryIndex> {
    let bytes = fs::read(path)
        .with_context(|| format!("Failed to open query index: {}", path.display()))?;
    let mut cursor = wincode::io::Cursor::new(bytes.as_slice());

    let (magic, kind, version, kmer_length, _smer_length, _spacing, flags, n_targets): QueryIndexHeader =
        wincode::deserialize_from(&mut cursor).context("Failed to decode query index header")?;
    validate_query_index_header(&magic, kind, version)?;

    let meta: QueryIndexMeta =
        wincode::deserialize_from(&mut cursor).context("Failed to decode query index metadata")?;
    if meta.len() != n_targets as usize {
        return Err(anyhow::anyhow!("Query index target count mismatch"));
    }

    let has_positions = flags & QUERY_INDEX_FLAG_POSITIONS != 0;
    let kmer_bytes = (kmer_length as usize).div_ceil(4);
    let entry = if has_positions {
        kmer_bytes + 4
    } else {
        kmer_bytes
    };
    let raw = &bytes[cursor.position()..];

    let read_u64 = |b: &[u8]| {
        let mut v = [0u8; 8];
        v[..kmer_bytes].copy_from_slice(&b[..kmer_bytes]);
        u64::from_le_bytes(v)
    };
    let read_u128 = |b: &[u8]| {
        let mut v = [0u8; 16];
        v[..kmer_bytes].copy_from_slice(&b[..kmer_bytes]);
        u128::from_le_bytes(v)
    };

    let mut off = 0usize;
    let mut targets = Vec::with_capacity(n_targets as usize);
    for (name, length, count) in meta {
        let count = count as usize;
        let end = off + count * entry;
        if end > raw.len() {
            return Err(anyhow::anyhow!("Query index truncated"));
        }
        let block = &raw[off..end];
        off = end;

        let mut target = TargetInfo {
            name,
            length: length as usize,
            syncmers: if kmer_length <= 32 {
                SyncmerSet::U64(RapidHashSet::default())
            } else {
                SyncmerSet::U128(RapidHashSet::default())
            },
            syncmer_positions: Vec::new(),
            positioned_syncmers: PositionedSyncmers::Empty,
        };

        if has_positions {
            let mut positions = Vec::with_capacity(count);
            if kmer_length <= 32 {
                let mut entries = Vec::with_capacity(count);
                let mut set = RapidHashSet::default();
                for i in 0..count {
                    let e = &block[i * entry..];
                    let v = read_u64(e);
                    let pos = u32::from_le_bytes(e[kmer_bytes..kmer_bytes + 4].try_into().unwrap())
                        as usize;
                    entries.push((v, pos));
                    positions.push(pos);
                    set.insert(v);
                }
                target.syncmers = SyncmerSet::U64(set);
                target.positioned_syncmers = PositionedSyncmers::U64(entries);
            } else {
                let mut entries = Vec::with_capacity(count);
                let mut set = RapidHashSet::default();
                for i in 0..count {
                    let e = &block[i * entry..];
                    let v = read_u128(e);
                    let pos = u32::from_le_bytes(e[kmer_bytes..kmer_bytes + 4].try_into().unwrap())
                        as usize;
                    entries.push((v, pos));
                    positions.push(pos);
                    set.insert(v);
                }
                target.syncmers = SyncmerSet::U128(set);
                target.positioned_syncmers = PositionedSyncmers::U128(entries);
            }
            target.syncmer_positions = positions;
        } else if kmer_length <= 32 {
            let mut set = RapidHashSet::default();
            for i in 0..count {
                set.insert(read_u64(&block[i * entry..]));
            }
            target.syncmers = SyncmerSet::U64(set);
        } else {
            let mut set = RapidHashSet::default();
            for i in 0..count {
                set.insert(read_u128(&block[i * entry..]));
            }
            target.syncmers = SyncmerSet::U128(set);
        }

        targets.push(target);
    }

    Ok(QueryIndex {
        targets,
        has_positions,
    })
}

/// Build and serialize a query index (with optional background masking)
pub fn build_query_index(config: &BuildQueryConfig) -> Result<()> {
    let start = Instant::now();
    let version = env!("CARGO_PKG_VERSION");
    eprintln!(
        "Skope v{version}; mode: index build-query; options: k={}, s={}, spacing={}, threads={}",
        config.kmer_length, config.smer_length, config.spacing, config.threads
    );

    let collect_positions = config.positions;
    let mut targets = if config.targets_path.is_dir() {
        process_targets_dir(
            &config.targets_path,
            config.kmer_length,
            config.smer_length,
            config.quiet,
            config.spacing,
            collect_positions,
        )?
    } else {
        let per_record = process_targets_file(
            &config.targets_path,
            config.kmer_length,
            config.smer_length,
            config.quiet,
            config.spacing,
            collect_positions,
        )?;
        if config.individual || per_record.len() <= 1 {
            per_record
        } else {
            let name = crate::derive_sample_name(&config.targets_path, false);
            vec![merge_targets(per_record, name)?]
        }
    };

    let mut flags = 0u8;
    if config.positions {
        flags |= QUERY_INDEX_FLAG_POSITIONS;
    }

    if !config.background_paths.is_empty() {
        let union = build_union(&targets, config.kmer_length);
        let removed = mask_background(
            &mut targets,
            &union,
            &config.background_paths,
            config.kmer_length,
            config.smer_length,
            config.threads,
            config.quiet,
        )?;
        flags |= QUERY_INDEX_FLAG_BACKGROUND;
        if !config.quiet {
            eprintln!("Masked {removed} background syncmers");
        }
    }

    save_query_index(
        &targets,
        config.kmer_length,
        config.smer_length,
        config.spacing,
        flags,
        config.output_path.as_deref(),
    )?;

    if !config.quiet {
        let total: usize = targets.iter().map(|t| t.syncmers.len()).sum();
        eprintln!(
            "Wrote query index: {} targets, {total} syncmers, positions={} in {:.1}s",
            targets.len(),
            config.positions,
            start.elapsed().as_secs_f64()
        );
    }
    Ok(())
}

pub fn run_containment_analysis(config: &ContainmentConfig) -> Result<()> {
    let start_time = Instant::now();
    let version = env!("CARGO_PKG_VERSION").to_string();
    let abundance_thresholds = normalize_abundance_thresholds(&config.abundance_thresholds);

    let mut options = format!(
        "k={}, s={}, threads={}",
        config.kmer_length, config.smer_length, config.threads
    );

    if config.sample_paths.len() > 1 {
        options.push_str(&format!(", samples={}", config.sample_paths.len()));
    }

    if !abundance_thresholds.is_empty() {
        options.push_str(&format!(
            ", abundance-thresholds={}",
            abundance_thresholds
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

    let is_targets_dir = config.targets_path.is_dir();
    let index_kind = if is_targets_dir {
        None
    } else {
        read_index_kind(&config.targets_path)
    };
    if let Some(kind) = index_kind
        && kind != IndexKind::Query
    {
        return Err(anyhow::anyhow!(
            "{} is a skope {kind:?} index, not a query index",
            config.targets_path.display()
        ));
    }
    let from_index = index_kind == Some(IndexKind::Query);
    eprintln!(
        "Skope v{}; mode: query{}; options: {}",
        version,
        if from_index {
            " (from index)"
        } else if is_targets_dir {
            " (from directory)"
        } else {
            ""
        },
        options
    );

    // Load a prebuilt query index else extract targets from fastx
    let targets_start = Instant::now();
    let need_positions = config.dump_syncmers_path.is_some() || config.confidence;
    let collect_positions = config.positions || need_positions;
    let mut targets = if from_index {
        let index = load_query_index(&config.targets_path)?;
        if need_positions && !index.has_positions {
            return Err(anyhow::anyhow!(
                "Query index built without --positions; cannot use --confidence or --dump-syncmers"
            ));
        }
        index.targets
    } else if is_targets_dir {
        process_targets_dir(
            &config.targets_path,
            config.kmer_length,
            config.smer_length,
            config.quiet,
            config.spacing,
            collect_positions,
        )?
    } else {
        let per_record = process_targets_file(
            &config.targets_path,
            config.kmer_length,
            config.smer_length,
            config.quiet,
            config.spacing,
            collect_positions,
        )?;
        if config.individual || per_record.len() <= 1 {
            per_record
        } else {
            let name = crate::derive_sample_name(&config.targets_path, false);
            vec![merge_targets(per_record, name)?]
        }
    };
    let targets_time = targets_start.elapsed();

    // Count syncmers shared between targets (only meaningful with >1 target)
    let (shared_syncmers, unique_across_all) = if targets.len() > 1 {
        if !config.quiet {
            eprint!("Counting shared syncmers…\r");
        }
        let mut syncmer_target_counts: HashMap<u64, usize> = HashMap::new();
        let mut syncmer_target_counts_u128: HashMap<u128, usize> = HashMap::new();

        for target in &targets {
            match &target.syncmers {
                SyncmerSet::U64(set) => {
                    for &syncmer in set {
                        *syncmer_target_counts.entry(syncmer).or_insert(0) += 1;
                    }
                }
                SyncmerSet::U128(set) => {
                    for &syncmer in set {
                        *syncmer_target_counts_u128.entry(syncmer).or_insert(0) += 1;
                    }
                }
            }
        }

        let shared = if config.kmer_length <= 32 {
            syncmer_target_counts
                .values()
                .filter(|&&count| count > 1)
                .count()
        } else {
            syncmer_target_counts_u128
                .values()
                .filter(|&&count| count > 1)
                .count()
        };

        let unique = if config.kmer_length <= 32 {
            syncmer_target_counts.len()
        } else {
            syncmer_target_counts_u128.len()
        };

        // Apply discriminatory filtering if enabled
        if config.discriminatory {
            for target in &mut targets {
                match &mut target.syncmers {
                    SyncmerSet::U64(set) => {
                        set.retain(|syncmer| {
                            syncmer_target_counts
                                .get(syncmer)
                                .is_none_or(|&count| count == 1)
                        });
                    }
                    SyncmerSet::U128(set) => {
                        set.retain(|syncmer| {
                            syncmer_target_counts_u128
                                .get(syncmer)
                                .is_none_or(|&count| count == 1)
                        });
                    }
                }
                target.positioned_syncmers.retain_in_set(&target.syncmers);
                target.syncmer_positions = target.positioned_syncmers.positions();
            }
        }

        (shared, unique)
    } else {
        (0, targets.first().map(|t| t.syncmers.len()).unwrap_or(0))
    };

    // Runtime masking
    if !config.background_paths.is_empty() {
        let union = build_union(&targets, config.kmer_length);
        let removed = mask_background(
            &mut targets,
            &union,
            &config.background_paths,
            config.kmer_length,
            config.smer_length,
            config.threads,
            config.quiet,
        )?;
        if !config.quiet {
            eprintln!("Masked {removed} background syncmers");
        }
    }

    if !config.quiet {
        eprint!("\r"); // Clear space
        let total_unique_syncmers: usize = targets.iter().map(|t| t.syncmers.len()).sum();
        let total_bp: usize = targets.iter().map(|t| t.length).sum();

        if config.discriminatory {
            eprintln!(
                "Targets: {} records ({}), {} discriminatory syncmers ({} shared syncmers dropped)",
                targets.len(),
                format_bp(total_bp),
                total_unique_syncmers,
                shared_syncmers
            );
        } else {
            let shared_pct = if unique_across_all > 0 {
                shared_syncmers as f64 / unique_across_all as f64 * 100.0
            } else {
                0.0
            };
            eprintln!(
                "Targets: {} records ({}), {} syncmers, of which {} ({:.1}%) shared by multiple targets",
                targets.len(),
                format_bp(total_bp),
                total_unique_syncmers,
                shared_syncmers,
                shared_pct
            );
        }
    }

    // Build set of all unique syncmers across targets
    if !config.quiet {
        eprint!("Building syncmer set…\r");
    }
    let targets_syncmers = Arc::new(build_union(&targets, config.kmer_length));

    // Dump syncmers (position + k-mer sequence) if requested
    if let Some(ref path) = config.dump_syncmers_path {
        let mut file = BufWriter::new(File::create(path)?);
        let k = config.kmer_length;
        let mut total = 0usize;
        for target in &targets {
            match &target.positioned_syncmers {
                PositionedSyncmers::U64(entries) => {
                    for &(value, pos) in entries {
                        let kmer = decode_u64(value, k);
                        writeln!(
                            file,
                            "{}\t{}\t{}",
                            target.name,
                            pos,
                            String::from_utf8_lossy(&kmer)
                        )?;
                        total += 1;
                    }
                }
                PositionedSyncmers::U128(entries) => {
                    for &(value, pos) in entries {
                        let kmer = decode_u128(value, k);
                        writeln!(
                            file,
                            "{}\t{}\t{}",
                            target.name,
                            pos,
                            String::from_utf8_lossy(&kmer)
                        )?;
                        total += 1;
                    }
                }
                PositionedSyncmers::Empty => {}
            }
        }
        file.flush()?;
        if !config.quiet {
            eprintln!("Dumped {} syncmers to {}", total, path.display());
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
                Arc::clone(&targets_syncmers),
                &abundance_thresholds,
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
            abundance_thresholds: abundance_thresholds.clone(),
            confidence: config.confidence,
        },
        samples: sample_results,
        total_timing: TimingStats {
            reference_processing_time: targets_time.as_secs_f64(),
            seqs_processing_time: total_seqs_time,
            analysis_time: total_analysis_time,
            total_time: total_time.as_secs_f64(),
            seqs_per_second: 0.0, // Not meaningful across samples
            bp_per_second: 0.0,   // Not meaningful across samples
        },
    };

    // Output results
    output_results(
        &report,
        config.output_path.as_ref(),
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
    sort_order: SortOrder,
    no_total: bool,
) -> Result<()> {
    let writer: Box<dyn Write> = if let Some(path) = output_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };

    let mut writer = writer;

    let mut sorted_report = report.clone();
    if sort_order != SortOrder::Original {
        for sample in &mut sorted_report.samples {
            sort_results(&mut sample.targets, sort_order);
        }
    }
    output_tsv(&mut writer, &sorted_report, no_total)?;

    Ok(())
}

/// Format a Wilson 95% confidence interval for a containment proportion as
/// `low-high` bounds (e.g. `0.49010-0.94330`).
fn format_containment_ci(hits: usize, target_kmers: usize) -> String {
    let (lo, hi) = wilson_interval(hits, target_kmers, WILSON_Z_95);
    format!("{:.3}-{:.3}", lo, hi)
}

fn format_ani_est(value: Option<f64>) -> String {
    value
        .map(|value| format!("{:.3}", value))
        .unwrap_or_else(|| "-".to_string())
}

fn format_patchiness(value: Option<PatchinessResult>) -> String {
    const PATCHINESS_P_THRESHOLD: f64 = 0.05;
    const MIN_DISPLAY_P: f64 = 1.0e-16;

    let Some(value) = value else {
        return "-".to_string();
    };

    if value.z <= 0.0 || value.p > PATCHINESS_P_THRESHOLD {
        return "-".to_string();
    }

    if value.p < MIN_DISPLAY_P {
        format!("{:.1}|<{:.0e}", value.z, MIN_DISPLAY_P)
    } else {
        format!("{:.1}|{:.0e}", value.z, value.p)
    }
}

fn output_tsv(writer: &mut dyn Write, report: &Report, no_total: bool) -> Result<()> {
    let mut thresholds = report.parameters.abundance_thresholds.clone();
    thresholds.sort_unstable();
    let confidence = report.parameters.confidence;

    // Build header with target column first
    let mut header = "target\tsample\tcontainment1\tcontainment1_hits".to_string();
    if confidence {
        header.push_str("\tcontainment1_ci");
        header.push_str("\tpatchiness");
        header.push_str("\tani_est");
    }
    for threshold in &thresholds {
        header.push_str(&format!(
            "\tcontainment{}\tcontainment{}_hits",
            threshold, threshold
        ));
        if confidence {
            header.push_str(&format!("\tcontainment{}_ci", threshold));
        }
    }
    header
        .push_str("\tmedian_nz_abundance\ttarget_kmers\ttarget_length\tsample_seqs\tsample_bases");
    writeln!(writer, "{}", header)?;

    // Output data rows for all samples
    for sample in &report.samples {
        for result in &sample.targets {
            let mut row = format!(
                "{}\t{}\t{:.3}\t{}",
                result.target,
                sample.sample_name,
                result.containment1,
                result.containment1_hits // Use containment1_hits for threshold 1 hits
            );
            if confidence {
                row.push_str(&format!(
                    "\t{}",
                    format_containment_ci(result.containment1_hits, result.target_kmers)
                ));
                row.push_str(&format!("\t{}", format_patchiness(result.patchiness)));
                row.push_str(&format!("\t{}", format_ani_est(result.ani_est)));
            }
            for threshold in &thresholds {
                let containment = result
                    .containment_at_threshold
                    .get(threshold)
                    .unwrap_or(&0.0);
                let hits = result.hits_at_threshold.get(threshold).unwrap_or(&0);
                row.push_str(&format!("\t{:.3}\t{}", containment, hits));
                if confidence {
                    row.push_str(&format!(
                        "\t{}",
                        format_containment_ci(*hits, result.target_kmers)
                    ));
                }
            }
            row.push_str(&format!(
                "\t{:.0}\t{}\t{}\t{}\t{}",
                result.median_nz_abundance,
                result.target_kmers,
                result.length,
                sample.total_stats.total_seqs_processed,
                sample.total_stats.total_bp_processed,
            ));
            writeln!(writer, "{}", row)?;
        }

        // Add TOTAL row for this sample
        if !no_total {
            let mut total_row = format!(
                "{}\t{}\t{:.3}\t{}",
                "TOTAL",
                sample.sample_name,
                sample.total_stats.total_containment1,
                sample.total_stats.total_containment1_hits
            );
            if confidence {
                total_row.push_str(&format!(
                    "\t{}",
                    format_containment_ci(
                        sample.total_stats.total_containment1_hits,
                        sample.total_stats.target_kmers
                    )
                ));
                total_row.push_str("\t-");
                total_row.push_str("\t-");
            }
            for threshold in &thresholds {
                let containment = sample
                    .total_stats
                    .total_containment_at_threshold
                    .get(threshold)
                    .unwrap_or(&0.0);
                let hits = (containment * sample.total_stats.target_kmers as f64).round() as usize;
                total_row.push_str(&format!("\t{:.3}\t{}", containment, hits));
                if confidence {
                    total_row.push_str(&format!(
                        "\t{}",
                        format_containment_ci(hits, sample.total_stats.target_kmers)
                    ));
                }
            }
            total_row.push_str(&format!(
                "\t-\t{}\t0\t{}\t{}",
                sample.total_stats.target_kmers,
                sample.total_stats.total_seqs_processed,
                sample.total_stats.total_bp_processed,
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
    fn test_patchiness_positive_for_clustered_hits() {
        let labels: Vec<bool> = (0..100).map(|idx| idx < 50).collect();
        let patchiness = calculate_patchiness_from_labels(&labels).unwrap();

        assert!(patchiness.z > 0.0);
        assert!(patchiness.p < 0.001);
    }

    #[test]
    fn test_patchiness_negative_for_alternating_hits() {
        let labels: Vec<bool> = (0..100).map(|idx| idx % 2 == 0).collect();
        let patchiness = calculate_patchiness_from_labels(&labels).unwrap();

        assert!(patchiness.z < 0.0);
        assert!(patchiness.p > 0.999);
    }

    #[test]
    fn test_patchiness_skips_degenerate_counts() {
        let labels: Vec<bool> = (0..100).map(|idx| idx < 5).collect();

        assert!(calculate_patchiness_from_labels(&labels).is_none());
    }

    #[test]
    fn test_format_patchiness_is_warning_only() {
        assert_eq!(
            format_patchiness(Some(PatchinessResult { z: -1.0, p: 0.9 })),
            "-"
        );
        assert_eq!(
            format_patchiness(Some(PatchinessResult { z: 1.0, p: 0.2 })),
            "-"
        );
        assert_eq!(
            format_patchiness(Some(PatchinessResult { z: 3.0, p: 0.001 })),
            "3.0|1e-3"
        );
    }

    #[test]
    fn test_format_patchiness_uses_p_floor() {
        assert_eq!(
            format_patchiness(Some(PatchinessResult { z: 9.0, p: 0.0 })),
            "9.0|<1e-16"
        );
    }
}

use crate::kmers::{FracMinHash, Kdust, KmerVec, Kmers};
use crate::stats::{WILSON_Z_95, normal_survival, wilson_interval};
use crate::{
    DEFAULT_KMER_LENGTH, DEFAULT_SMER_LENGTH, FixedRapidHasher, INDEX_MAGIC, IndexKind, Input,
    ProcessingStats, RapidHashSet, SeqProcessor, check_index_complexity, format_bp, process_stream,
    validate_complexity, validate_fraction, validate_k_s,
};
use anyhow::{Context, Result};
use parking_lot::Mutex;
use std::collections::HashMap;
use std::hash::Hash;
use std::io::{Read, Write};
use std::sync::Arc;

#[cfg(feature = "cli")]
mod io;
#[cfg(feature = "cli")]
pub use io::{
    BuildQueryConfig, ContainmentConfig, print_query_index_info, read_query_index_meta,
    run_build_query, run_query,
};

/// Alias for abund counts
type CountDepth = u16;

type CountMap<T> = HashMap<T, CountDepth, FixedRapidHasher>;

/// Sort order for results
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    Original,    // Original order from reference file (default)
    Target,      // Alphabetical by target name
    Containment, // Descending by containment1 (highest first)
}

impl std::str::FromStr for SortOrder {
    type Err = anyhow::Error;

    fn from_str(sort: &str) -> Result<Self> {
        match sort {
            "input" => Ok(Self::Original),
            "target" => Ok(Self::Target),
            "containment" => Ok(Self::Containment),
            _ => Err(anyhow::anyhow!(
                "Invalid sort order {sort}: expected containment, target or input"
            )),
        }
    }
}

/// Zero-cost (hopefully?) abstraction over u64 and u128 k-mer sets
#[derive(Debug, Clone)]
enum KmerSet {
    U64(RapidHashSet<u64>),
    U128(RapidHashSet<u128>),
}

impl KmerSet {
    /// Empty set as wide as `kmer_length` needs
    fn new(kmer_length: u8) -> Self {
        if kmer_length <= 32 {
            KmerSet::U64(RapidHashSet::default())
        } else {
            KmerSet::U128(RapidHashSet::default())
        }
    }

    fn len(&self) -> usize {
        match self {
            KmerSet::U64(set) => set.len(),
            KmerSet::U128(set) => set.len(),
        }
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn extend(&mut self, other: &KmerSet) {
        match (self, other) {
            (KmerSet::U64(a), KmerSet::U64(b)) => a.extend(b.iter().copied()),
            (KmerSet::U128(a), KmerSet::U128(b)) => a.extend(b.iter().copied()),
            _ => panic!("Cannot extend KmerSet: mismatched variants"),
        }
    }

    /// Remove all k-mers from --background stream, returning count
    fn retain_not_in(&mut self, other: &KmerSet) -> usize {
        match (self, other) {
            (KmerSet::U64(a), KmerSet::U64(b)) => {
                let before = a.len();
                a.retain(|s| !b.contains(s));
                before - a.len()
            }
            (KmerSet::U128(a), KmerSet::U128(b)) => {
                let before = a.len();
                a.retain(|s| !b.contains(s));
                before - a.len()
            }
            _ => panic!("Cannot mask KmerSet: mismatched variants"),
        }
    }
}

#[derive(Debug, Clone)]
enum PositionedKmers {
    U64(Vec<(u64, usize)>),
    U128(Vec<(u128, usize)>),
    Empty,
}

impl PositionedKmers {
    fn offset_positions(&mut self, offset: usize) {
        match self {
            PositionedKmers::U64(entries) => {
                for (_, position) in entries {
                    *position += offset;
                }
            }
            PositionedKmers::U128(entries) => {
                for (_, position) in entries {
                    *position += offset;
                }
            }
            PositionedKmers::Empty => {}
        }
    }

    fn extend(&mut self, other: PositionedKmers) {
        match self {
            PositionedKmers::U64(a) => match other {
                PositionedKmers::U64(b) => a.extend(b),
                PositionedKmers::Empty => {}
                _ => panic!("Cannot extend PositionedKmers: mismatched variants"),
            },
            PositionedKmers::U128(a) => match other {
                PositionedKmers::U128(b) => a.extend(b),
                PositionedKmers::Empty => {}
                _ => panic!("Cannot extend PositionedKmers: mismatched variants"),
            },
            PositionedKmers::Empty => *self = other,
        }
    }

    fn retain_in_set(&mut self, kmers: &KmerSet) {
        match (self, kmers) {
            (PositionedKmers::U64(entries), KmerSet::U64(set)) => {
                entries.retain(|(kmer, _)| set.contains(kmer));
            }
            (PositionedKmers::U128(entries), KmerSet::U128(set)) => {
                entries.retain(|(kmer, _)| set.contains(kmer));
            }
            (PositionedKmers::Empty, _) => {}
            _ => panic!("Cannot filter PositionedKmers: mismatched variants"),
        }
    }

    fn positions(&self) -> Vec<usize> {
        match self {
            PositionedKmers::U64(entries) => {
                entries.iter().map(|(_, position)| *position).collect()
            }
            PositionedKmers::U128(entries) => {
                entries.iter().map(|(_, position)| *position).collect()
            }
            PositionedKmers::Empty => Vec::new(),
        }
    }

    fn len(&self) -> usize {
        match self {
            PositionedKmers::U64(entries) => entries.len(),
            PositionedKmers::U128(entries) => entries.len(),
            PositionedKmers::Empty => 0,
        }
    }
}

#[derive(Debug, Clone)]
struct TargetInfo {
    name: String,
    length: usize,
    kmers: KmerSet,
    kmer_positions: Vec<usize>,
    positioned_kmers: PositionedKmers,
}

#[derive(Debug, Clone, Copy)]
struct PatchinessResult {
    z: f64,
    p: f64,
}

#[derive(Debug)]
struct ContainmentResult {
    target: String,
    length: usize,
    target_kmers: usize,
    containment1_hits: usize,
    containment1: f64,
    ani_est: Option<f64>,
    median_nz_abundance: f64,
    containment_at_threshold: HashMap<usize, f64>, // threshold -> containment
    hits_at_threshold: HashMap<usize, usize>,      // threshold -> hit count
    patchiness: Option<PatchinessResult>,
}

#[derive(Debug)]
struct TotalStats {
    target_kmers: usize,
    total_containment1_hits: usize,
    total_containment1: f64,
    total_seqs_processed: u64,
    total_bp_processed: u64,
    total_hits_at_threshold: HashMap<usize, usize>,
}

/// Results for a single sample in multi-sample mode
#[derive(Debug)]
struct SampleResults {
    sample_name: String,
    targets: Vec<ContainmentResult>,
    total_stats: TotalStats,
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

/// Processor for collecting target sequence records with k-mers
#[derive(Clone)]
struct TargetsProcessor {
    kmers: Kmers,
    fmh: FracMinHash,
    kdust: Kdust,
    positions: Vec<usize>,
    targets: Arc<Mutex<Vec<TargetInfo>>>,
    collect_positions: bool,
}

impl TargetsProcessor {
    fn new(
        kmer_length: u8,
        smer_length: u8,
        fmh: FracMinHash,
        kdust: Kdust,
        collect_positions: bool,
    ) -> Self {
        Self {
            kmers: Kmers::new(kmer_length, smer_length),
            fmh,
            kdust,
            positions: Vec::new(),
            targets: Arc::default(),
            collect_positions,
        }
    }
}

impl SeqProcessor for TargetsProcessor {
    fn process(&mut self, id: &[u8], seqs: &[&[u8]]) -> paraseq::Result<()> {
        let mut positions = self.collect_positions.then_some(&mut self.positions);
        let values = match positions.as_deref_mut() {
            Some(positions) => self.kmers.pool_with_positions(seqs, positions),
            None => self.kmers.pool(seqs),
        };
        self.fmh.retain(values, positions.as_deref_mut());
        self.kdust.retain(values, positions);

        // Build unique k-mer set for this target
        let kmers = match &*values {
            KmerVec::U64(vec) => {
                let set: RapidHashSet<u64> = vec.iter().copied().collect();
                KmerSet::U64(set)
            }
            KmerVec::U128(vec) => {
                let set: RapidHashSet<u128> = vec.iter().copied().collect();
                KmerSet::U128(set)
            }
        };

        let positioned_kmers = if self.collect_positions {
            match &*values {
                KmerVec::U64(vec) => PositionedKmers::U64(
                    vec.iter()
                        .copied()
                        .zip(self.positions.iter().copied())
                        .collect(),
                ),
                KmerVec::U128(vec) => PositionedKmers::U128(
                    vec.iter()
                        .copied()
                        .zip(self.positions.iter().copied())
                        .collect(),
                ),
            }
        } else {
            PositionedKmers::Empty
        };

        let kmer_positions = if self.collect_positions {
            self.positions.clone()
        } else {
            Vec::new()
        };

        self.targets.lock().push(TargetInfo {
            name: String::from_utf8_lossy(id).to_string(),
            length: seqs.iter().map(|seq| seq.len()).sum(),
            kmers,
            kmer_positions,
            positioned_kmers,
        });

        Ok(())
    }

    fn flush(&mut self) -> paraseq::Result<()> {
        Ok(())
    }
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
        let mut positioned_kmers = t.positioned_kmers;
        positioned_kmers.offset_positions(offset);
        merged.kmers.extend(&t.kmers);
        merged.kmer_positions.extend(
            t.kmer_positions
                .into_iter()
                .map(|position| position + offset),
        );
        merged.positioned_kmers.extend(positioned_kmers);
        merged.length += t.length;
    }

    Ok(merged)
}

/// Targets of one group. A lone file's records stay apart under `--individual`, or when there is
/// only one, keeping their names. Otherwise they merge under the group's name
fn group_targets(
    targets: Vec<TargetInfo>,
    name: &str,
    lone_file: bool,
    individual: bool,
) -> Result<Vec<TargetInfo>> {
    if lone_file && (individual || targets.len() <= 1) {
        Ok(targets)
    } else {
        Ok(vec![merge_targets(targets, name.to_string())?])
    }
}

/// Warning naming any targets without k-mers
fn empty_targets_warning(targets: &[TargetInfo]) -> Option<String> {
    let empty: Vec<&str> = targets
        .iter()
        .filter(|t| t.kmers.is_empty())
        .map(|t| t.name.as_str())
        .collect();
    if empty.is_empty() {
        return None;
    }

    const SHOWN: usize = 3;
    let named = empty
        .iter()
        .take(SHOWN)
        .copied()
        .collect::<Vec<_>>()
        .join(", ");
    let more = match empty.len().saturating_sub(SHOWN) {
        0 => String::new(),
        n => format!(" and {n} more"),
    };
    Some(format!(
        "Warning: {} of {} targets yielded no k-mers: {named}{more}",
        empty.len(),
        targets.len()
    ))
}

/// Processor for counting k-mer depths from sequences
#[derive(Clone)]
struct CountProcessor {
    kmers: Kmers,
    targets_kmers: Arc<KmerSet>,
    local_counts: AbundanceMap,
    global_counts: Arc<Mutex<AbundanceMap>>,
}

impl CountProcessor {
    fn new(
        kmer_length: u8,
        smer_length: u8,
        targets_kmers: Arc<KmerSet>,
        global_counts: Arc<Mutex<AbundanceMap>>,
    ) -> Self {
        Self {
            kmers: Kmers::new(kmer_length, smer_length),
            targets_kmers,
            local_counts: AbundanceMap::new(kmer_length),
            global_counts,
        }
    }
}

impl SeqProcessor for CountProcessor {
    fn process(&mut self, _id: &[u8], seqs: &[&[u8]]) -> paraseq::Result<()> {
        let distinct = seqs.len() > 1;
        match (
            self.kmers.pool(seqs),
            &*self.targets_kmers,
            &mut self.local_counts,
        ) {
            (KmerVec::U64(kmers), KmerSet::U64(targets), AbundanceMap::U64(counts)) => {
                count_target_kmers(kmers, targets, counts, distinct)
            }
            (KmerVec::U128(kmers), KmerSet::U128(targets), AbundanceMap::U128(counts)) => {
                count_target_kmers(kmers, targets, counts, distinct)
            }
            _ => unreachable!("k-mer width does not match targets"),
        }
        Ok(())
    }

    fn flush(&mut self) -> paraseq::Result<()> {
        match (&mut self.local_counts, &mut *self.global_counts.lock()) {
            (AbundanceMap::U64(local), AbundanceMap::U64(global)) => {
                merge_counts(local.drain(), global)
            }
            (AbundanceMap::U128(local), AbundanceMap::U128(global)) => {
                merge_counts(local.drain(), global)
            }
            _ => unreachable!("k-mer width does not match targets"),
        }
        Ok(())
    }
}

/// Merge counts with saturation
fn merge_counts<T: Eq + Hash>(
    counts: impl IntoIterator<Item = (T, CountDepth)>,
    global: &mut CountMap<T>,
) {
    for (kmer, count) in counts {
        global
            .entry(kmer)
            .and_modify(|e| *e = e.saturating_add(count))
            .or_insert(count);
    }
}

/// Count occurrences, or distinct hits for pooled mates
fn count_target_kmers<T: Copy + Ord + Hash>(
    kmers: &mut Vec<T>,
    targets: &RapidHashSet<T>,
    counts: &mut CountMap<T>,
    distinct: bool,
) {
    kmers.retain(|kmer| targets.contains(kmer));
    if distinct {
        kmers.sort_unstable();
        kmers.dedup();
    }
    for &kmer in kmers.iter() {
        counts
            .entry(kmer)
            .and_modify(|e| *e = e.saturating_add(1))
            .or_insert(1);
    }
}

/// Enum to abstract over u64 and u128 abundance maps
#[derive(Clone)]
enum AbundanceMap {
    U64(CountMap<u64>),
    U128(CountMap<u128>),
}

impl AbundanceMap {
    fn new(kmer_length: u8) -> Self {
        if kmer_length <= 32 {
            Self::U64(CountMap::default())
        } else {
            Self::U128(CountMap::default())
        }
    }

    fn len(&self) -> usize {
        match self {
            Self::U64(map) => map.len(),
            Self::U128(map) => map.len(),
        }
    }

    /// The counted k-mers
    fn into_kmers(self) -> KmerSet {
        match self {
            Self::U64(map) => KmerSet::U64(map.into_keys().collect()),
            Self::U128(map) => KmerSet::U128(map.into_keys().collect()),
        }
    }
}

const MIN_TARGET_KMERS_FOR_ANI_ADJUSTMENT: usize = 50;
const MIN_NONZERO_KMERS_FOR_ANI_ADJUSTMENT: usize = 25;
const MIN_DEPTH_BIN_FOR_ANI_ADJUSTMENT: usize = 3;
const MAX_MEDIAN_DEPTH_FOR_ANI_ADJUSTMENT: f64 = 2.0;
// Gates for displaying an ANI estimate: minimum target k-mers and minimum ANI
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
    // Need enough k-mers and at least one contained k-mer
    if kmer_length == 0 || target_kmers < MIN_TARGET_KMERS_FOR_ANI_DISPLAY || containment <= 0.0 {
        return None;
    }

    let kmer_identity = lambda
        .and_then(|lambda| coverage_adjusted_kmer_identity(containment, lambda, 1))
        .unwrap_or_else(|| containment.clamp(0.0, 1.0));

    let ani = kmer_identity.powf(1.0 / kmer_length as f64);

    // Suppress low-confidence estimates
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
    positioned_kmers: &[(u64, usize)],
    abundances: &CountMap<u64>,
    threshold: usize,
) -> Option<PatchinessResult> {
    if threshold == 0 {
        return None;
    }

    let mut occurrence_counts: HashMap<u64, usize> = HashMap::new();
    for &(kmer, _) in positioned_kmers {
        *occurrence_counts.entry(kmer).or_insert(0) += 1;
    }

    let labels: Vec<bool> = positioned_kmers
        .iter()
        .filter(|(kmer, _)| occurrence_counts.get(kmer).is_some_and(|&count| count == 1))
        .map(|(kmer, _)| abundances.get(kmer).copied().unwrap_or(0) >= threshold as CountDepth)
        .collect();

    calculate_patchiness_from_labels(&labels)
}

fn calculate_patchiness_u128(
    positioned_kmers: &[(u128, usize)],
    abundances: &CountMap<u128>,
    threshold: usize,
) -> Option<PatchinessResult> {
    if threshold == 0 {
        return None;
    }

    let mut occurrence_counts: HashMap<u128, usize> = HashMap::new();
    for &(kmer, _) in positioned_kmers {
        *occurrence_counts.entry(kmer).or_insert(0) += 1;
    }

    let labels: Vec<bool> = positioned_kmers
        .iter()
        .filter(|(kmer, _)| occurrence_counts.get(kmer).is_some_and(|&count| count == 1))
        .map(|(kmer, _)| abundances.get(kmer).copied().unwrap_or(0) >= threshold as CountDepth)
        .collect();

    calculate_patchiness_from_labels(&labels)
}

fn calculate_patchiness(
    positioned_kmers: &PositionedKmers,
    abundance_map: &AbundanceMap,
    threshold: usize,
) -> Option<PatchinessResult> {
    match (positioned_kmers, abundance_map) {
        (PositionedKmers::U64(positioned), AbundanceMap::U64(map)) => {
            calculate_patchiness_u64(positioned, map, threshold)
        }
        (PositionedKmers::U128(positioned), AbundanceMap::U128(map)) => {
            calculate_patchiness_u128(positioned, map, threshold)
        }
        (PositionedKmers::Empty, _) => None,
        _ => panic!("Mismatch between PositionedKmers and AbundanceMap types"),
    }
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
            let target_kmers = target.kmers.len();
            let mut abundances: Vec<CountDepth> = Vec::new();
            let mut contained_count = 0;

            // Collect abundances for all unique k-mers in this target
            match (&target.kmers, abundance_map) {
                (KmerSet::U64(set), AbundanceMap::U64(map)) => {
                    for &kmer in set {
                        let abundance = map.get(&kmer).copied().unwrap_or(0);
                        abundances.push(abundance);
                        if abundance > 0 {
                            contained_count += 1;
                        }
                    }
                }
                (KmerSet::U128(set), AbundanceMap::U128(map)) => {
                    for &kmer in set {
                        let abundance = map.get(&kmer).copied().unwrap_or(0);
                        abundances.push(abundance);
                        if abundance > 0 {
                            contained_count += 1;
                        }
                    }
                }
                _ => panic!("Mismatch between KmerSet and AbundanceMap types"),
            }

            let containment1 = if target_kmers > 0 {
                contained_count as f64 / target_kmers as f64
            } else {
                0.0
            };

            // Ignore zero-abundance k-mers for median calc
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
                calculate_patchiness(&target.positioned_kmers, abundance_map, 1)
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

impl SampleResults {
    /// Results of one sample, totalled over its targets
    fn new(
        sample_name: String,
        targets: Vec<ContainmentResult>,
        abundance_thresholds: &[usize],
        total_seqs_processed: u64,
        total_bp_processed: u64,
    ) -> Self {
        let target_kmers: usize = targets.iter().map(|r| r.target_kmers).sum();
        let total_containment1_hits: usize = targets.iter().map(|r| r.containment1_hits).sum();
        let total_containment1 = if target_kmers > 0 {
            total_containment1_hits as f64 / target_kmers as f64
        } else {
            0.0
        };

        // Sum integer hits to avoid float truncation
        let total_hits_at_threshold = abundance_thresholds
            .iter()
            .map(|&threshold| {
                let hits = targets
                    .iter()
                    .map(|r| r.hits_at_threshold[&threshold])
                    .sum();
                (threshold, hits)
            })
            .collect();

        Self {
            sample_name,
            targets,
            total_stats: TotalStats {
                target_kmers,
                total_containment1_hits,
                total_containment1,
                total_seqs_processed,
                total_bp_processed,
                total_hits_at_threshold,
            },
        }
    }
}

/// Summary of a counted stream, less its rate
fn counted_summary(label: &str, stats: &ProcessingStats, found: usize, total: usize) -> String {
    format!(
        "{label}: {} records ({}), found {found} of {total} distinct target k-mers",
        stats.total_seqs,
        format_bp(stats.total_bp as usize)
    )
}

/// Combined set of every target's k-mers
fn build_union(targets: &[TargetInfo], kmer_length: u8) -> KmerSet {
    if kmer_length <= 32 {
        let mut set = RapidHashSet::default();
        for t in targets {
            if let KmerSet::U64(s) = &t.kmers {
                set.extend(s.iter());
            }
        }
        KmerSet::U64(set)
    } else {
        let mut set = RapidHashSet::default();
        for t in targets {
            if let KmerSet::U128(s) = &t.kmers {
                set.extend(s.iter());
            }
        }
        KmerSet::U128(set)
    }
}

/// Drop background k-mers from every target, re-syncing positions like the discriminatory filter
fn mask_targets(targets: &mut [TargetInfo], background: &KmerSet) {
    for target in targets {
        target.kmers.retain_not_in(background);
        target.positioned_kmers.retain_in_set(&target.kmers);
        target.kmer_positions = target.positioned_kmers.positions();
    }
}

/// Count k-mers shared by multiple targets, keeping only unshared k-mers when `discriminatory`.
/// Returns the shared and distinct k-mer counts
fn share_kmers(
    targets: &mut [TargetInfo],
    kmer_length: u8,
    discriminatory: bool,
) -> (usize, usize) {
    if targets.len() <= 1 {
        return (0, targets.first().map(|t| t.kmers.len()).unwrap_or(0));
    }
    let mut kmer_target_counts: HashMap<u64, usize> = HashMap::new();
    let mut kmer_target_counts_u128: HashMap<u128, usize> = HashMap::new();

    for target in targets.iter() {
        match &target.kmers {
            KmerSet::U64(set) => {
                for &kmer in set {
                    *kmer_target_counts.entry(kmer).or_insert(0) += 1;
                }
            }
            KmerSet::U128(set) => {
                for &kmer in set {
                    *kmer_target_counts_u128.entry(kmer).or_insert(0) += 1;
                }
            }
        }
    }

    let shared = if kmer_length <= 32 {
        kmer_target_counts
            .values()
            .filter(|&&count| count > 1)
            .count()
    } else {
        kmer_target_counts_u128
            .values()
            .filter(|&&count| count > 1)
            .count()
    };

    let unique = if kmer_length <= 32 {
        kmer_target_counts.len()
    } else {
        kmer_target_counts_u128.len()
    };

    // Apply discriminatory filtering if enabled
    if discriminatory {
        for target in targets.iter_mut() {
            match &mut target.kmers {
                KmerSet::U64(set) => {
                    set.retain(|kmer| kmer_target_counts.get(kmer).is_none_or(|&count| count == 1));
                }
                KmerSet::U128(set) => {
                    set.retain(|kmer| {
                        kmer_target_counts_u128
                            .get(kmer)
                            .is_none_or(|&count| count == 1)
                    });
                }
            }
            target.positioned_kmers.retain_in_set(&target.kmers);
            target.kmer_positions = target.positioned_kmers.positions();
        }
    }

    (shared, unique)
}

/// The targets summary line
fn targets_summary(
    targets: &[TargetInfo],
    discriminatory: bool,
    shared_kmers: usize,
    unique_across_all: usize,
) -> String {
    let total_unique_kmers: usize = targets.iter().map(|t| t.kmers.len()).sum();
    let total_bp: usize = targets.iter().map(|t| t.length).sum();
    if discriminatory {
        format!(
            "Targets: {} records ({}), {} discriminatory k-mers ({} shared k-mers dropped)",
            targets.len(),
            format_bp(total_bp),
            total_unique_kmers,
            shared_kmers
        )
    } else {
        let shared_pct = if unique_across_all > 0 {
            shared_kmers as f64 / unique_across_all as f64 * 100.0
        } else {
            0.0
        };
        format!(
            "Targets: {} records ({}), {} k-mers, of which {} ({:.1}%) shared by multiple targets",
            targets.len(),
            format_bp(total_bp),
            total_unique_kmers,
            shared_kmers,
            shared_pct
        )
    }
}

// ── Query index (.sk, kind = query) ─────────────────────────────────────────
// Header + wincode metadata, then per target `count` raw entries of `ceil(k/4)` value
// bytes, each with a trailing u32 position when positions are stored.
const QUERY_INDEX_VERSION: u8 = 1;
const QUERY_INDEX_FLAG_POSITIONS: u8 = 0b001;
const QUERY_INDEX_FLAG_BACKGROUND: u8 = 0b010;

// magic, kind, version, k, s, flags, n_targets, fraction (FracMinHash), complexity (kdust)
type QueryIndexHeader = ([u8; 4], u8, u8, u8, u8, u8, u32, f64, f32);
type QueryIndexMeta = Vec<(String, u64, u64)>; // (name, length, entry_count)

struct QueryIndex {
    kmer_length: u8,
    smer_length: u8,
    flags: u8,
    fraction: f64,
    complexity: f32,
    targets: Vec<TargetInfo>,
}

impl QueryIndex {
    /// Reject an index without positions when they are needed
    fn require_positions(&self, needed: bool) -> Result<()> {
        if needed && self.flags & QUERY_INDEX_FLAG_POSITIONS == 0 {
            return Err(anyhow::anyhow!(
                "Query index built without --positions; cannot use --confidence or --dump-kmers"
            ));
        }
        Ok(())
    }
}

/// Reject a fraction that disagrees with a prebuilt index's.
/// Exact comparison: a CLI-parsed float round-trips through the header bit-identically
fn check_index_fraction(fraction: f64, index_fraction: f64) -> Result<()> {
    if fraction < 1.0 && fraction != index_fraction {
        return Err(anyhow::anyhow!(
            "--fraction {fraction} conflicts with the query index's fraction {index_fraction}; omit --fraction to use the index's"
        ));
    }
    Ok(())
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

/// Serialise targets as a query index
fn write_query_index(
    writer: &mut dyn Write,
    targets: &[TargetInfo],
    kmer_length: u8,
    smer_length: u8,
    fraction: f64,
    complexity: f32,
    flags: u8,
) -> Result<()> {
    let has_positions = flags & QUERY_INDEX_FLAG_POSITIONS != 0;
    let kmer_bytes = (kmer_length as usize).div_ceil(4);

    if has_positions && targets.iter().any(|t| t.length > u32::MAX as usize) {
        return Err(anyhow::anyhow!(
            "Target exceeds {} bp; positions cannot be stored (build without --positions)",
            u32::MAX
        ));
    }

    let meta: QueryIndexMeta = targets
        .iter()
        .map(|t| {
            let count = if has_positions {
                t.positioned_kmers.len()
            } else {
                t.kmers.len()
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
        flags,
        targets.len() as u32,
        fraction,
        complexity,
    );
    writer
        .write_all(&wincode::serialize(&header).context("Failed to encode query index header")?)?;
    writer
        .write_all(&wincode::serialize(&meta).context("Failed to encode query index metadata")?)?;

    for t in targets {
        if has_positions {
            match &t.positioned_kmers {
                PositionedKmers::U64(entries) => {
                    for &(v, pos) in entries {
                        writer.write_all(&v.to_le_bytes()[..kmer_bytes])?;
                        writer.write_all(&(pos as u32).to_le_bytes())?;
                    }
                }
                PositionedKmers::U128(entries) => {
                    for &(v, pos) in entries {
                        writer.write_all(&v.to_le_bytes()[..kmer_bytes])?;
                        writer.write_all(&(pos as u32).to_le_bytes())?;
                    }
                }
                PositionedKmers::Empty => {}
            }
        } else {
            match &t.kmers {
                KmerSet::U64(set) => {
                    for &v in set {
                        writer.write_all(&v.to_le_bytes()[..kmer_bytes])?;
                    }
                }
                KmerSet::U128(set) => {
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

/// Deserialise a query index
fn read_query_index(bytes: &[u8]) -> Result<QueryIndex> {
    let mut cursor = wincode::io::Cursor::new(bytes);

    let (magic, kind, version, kmer_length, smer_length, flags, n_targets, fraction, complexity): QueryIndexHeader =
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
            kmers: if kmer_length <= 32 {
                KmerSet::U64(RapidHashSet::default())
            } else {
                KmerSet::U128(RapidHashSet::default())
            },
            kmer_positions: Vec::new(),
            positioned_kmers: PositionedKmers::Empty,
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
                target.kmers = KmerSet::U64(set);
                target.positioned_kmers = PositionedKmers::U64(entries);
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
                target.kmers = KmerSet::U128(set);
                target.positioned_kmers = PositionedKmers::U128(entries);
            }
            target.kmer_positions = positions;
        } else if kmer_length <= 32 {
            let mut set = RapidHashSet::default();
            for i in 0..count {
                set.insert(read_u64(&block[i * entry..]));
            }
            target.kmers = KmerSet::U64(set);
        } else {
            let mut set = RapidHashSet::default();
            for i in 0..count {
                set.insert(read_u128(&block[i * entry..]));
            }
            target.kmers = KmerSet::U128(set);
        }

        targets.push(target);
    }

    Ok(QueryIndex {
        kmer_length,
        smer_length,
        flags,
        fraction,
        complexity,
        targets,
    })
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

fn output_tsv(
    writer: &mut dyn Write,
    samples: &[SampleResults],
    thresholds: &[usize],
    confidence: bool,
    no_total: bool,
) -> Result<()> {
    // Build header with target column first
    let mut header = "target\tsample\tcontainment1\tcontainment1_hits".to_string();
    if confidence {
        header.push_str("\tcontainment1_ci");
        header.push_str("\tpatchiness");
        header.push_str("\tani_est");
    }
    for threshold in thresholds {
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
    for sample in samples {
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
            for threshold in thresholds {
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
            for threshold in thresholds {
                let target_kmers = sample.total_stats.target_kmers;
                let hits = sample.total_stats.total_hits_at_threshold[threshold];
                let containment = if target_kmers > 0 {
                    hits as f64 / target_kmers as f64
                } else {
                    0.0
                };
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

// ── Streaming session ───────────────────────────────────────────────────────

/// Parameters of a [`QuerySession`], as for `skope query`
#[derive(Debug, Clone)]
pub struct QueryOptions {
    pub kmer_length: u8,
    /// 0 takes every k-mer
    pub smer_length: u8,
    pub individual: bool,
    pub discriminatory: bool,
    pub confidence: bool,
    pub abundance_thresholds: Vec<usize>,
    pub sort_order: SortOrder,
    pub no_total: bool,
    pub fraction: f64,
    pub complexity: f32,
    pub limit_bp: Option<u64>,
}

impl Default for QueryOptions {
    fn default() -> Self {
        Self {
            kmer_length: DEFAULT_KMER_LENGTH,
            smer_length: DEFAULT_SMER_LENGTH,
            individual: false,
            discriminatory: false,
            confidence: false,
            abundance_thresholds: vec![10],
            sort_order: SortOrder::Containment,
            no_total: false,
            fraction: 1.0,
            complexity: 0.0,
            limit_bp: None,
        }
    }
}

/// `skope index build-query` and `skope query` over fastx streams on one thread, e.g. in wasm.
/// Add target streams or load an index, optionally add background streams, then add samples.
/// One target stream acts as a file, several as a directory of files
pub struct QuerySession {
    options: QueryOptions,
    thresholds: Vec<usize>,
    /// Name and targets of each target stream, until grouped
    groups: Vec<(String, Vec<TargetInfo>)>,
    targets: Vec<TargetInfo>,
    flags: u8,
    /// Union counted against and the background k-mers found in it, until masked
    background: Option<(Arc<KmerSet>, KmerSet)>,
    /// Union counted against and the targets summary, once prepared
    prepared: Option<(Arc<KmerSet>, String)>,
    samples: Vec<SampleResults>,
}

impl QuerySession {
    pub fn new(options: QueryOptions) -> Result<Self> {
        validate_k_s(options.kmer_length, options.smer_length)?;
        validate_fraction(options.fraction)?;
        validate_complexity(options.complexity)?;
        Ok(Self {
            thresholds: normalize_abundance_thresholds(&options.abundance_thresholds),
            options,
            groups: Vec::new(),
            targets: Vec::new(),
            flags: 0,
            background: None,
            prepared: None,
            samples: Vec::new(),
        })
    }

    /// K-mer length, an index's once loaded
    pub fn kmer_length(&self) -> u8 {
        self.options.kmer_length
    }

    /// S-mer length, an index's once loaded
    pub fn smer_length(&self) -> u8 {
        self.options.smer_length
    }

    /// Add a target stream named after its file
    pub fn add_targets(
        &mut self,
        name: &str,
        stream: impl Read,
        on_batch: impl FnMut(&ProcessingStats),
    ) -> Result<()> {
        anyhow::ensure!(
            self.targets.is_empty() && self.background.is_none() && self.prepared.is_none(),
            "Targets must precede background and samples, and cannot join an index"
        );
        let QueryOptions {
            kmer_length: k,
            smer_length: s,
            confidence,
            ..
        } = self.options;
        let mut processor = TargetsProcessor::new(
            k,
            s,
            FracMinHash::from_fraction(self.options.fraction),
            Kdust::from_threshold(self.options.complexity, k),
            confidence,
        );
        process_stream(Input::File(stream), &mut processor, None, on_batch)?;
        if confidence {
            self.flags |= QUERY_INDEX_FLAG_POSITIONS;
        }
        let targets = std::mem::take(&mut *processor.targets.lock());
        self.groups.push((name.to_string(), targets));
        Ok(())
    }

    /// Load a prebuilt query index in place of targets. Its k, s, fraction and complexity prevail
    pub fn load_index(&mut self, bytes: &[u8]) -> Result<()> {
        anyhow::ensure!(
            self.groups.is_empty() && self.targets.is_empty() && self.background.is_none(),
            "An index replaces targets, so must come first"
        );
        let index = read_query_index(bytes)?;
        check_index_fraction(self.options.fraction, index.fraction)?;
        check_index_complexity(self.options.complexity, index.complexity)?;
        index.require_positions(self.options.confidence)?;
        self.options.kmer_length = index.kmer_length;
        self.options.smer_length = index.smer_length;
        self.options.fraction = index.fraction;
        self.options.complexity = index.complexity;
        self.flags = index.flags;
        self.targets = index.targets;
        Ok(())
    }

    /// Group target streams, rejecting duplicate names as for a directory
    fn group(&mut self) -> Result<()> {
        let lone_file = self.groups.len() == 1;
        for (i, (name, _)) in self.groups.iter().enumerate() {
            anyhow::ensure!(
                self.groups[..i].iter().all(|(other, _)| other != name),
                "Duplicate target name '{name}'"
            );
        }
        for (name, targets) in std::mem::take(&mut self.groups) {
            let grouped = group_targets(targets, &name, lone_file, self.options.individual)?;
            self.targets.extend(grouped);
        }
        Ok(())
    }

    /// Find target k-mers in a background stream, to drop them from targets. Returns its summary
    pub fn add_background(
        &mut self,
        stream: impl Read,
        on_batch: impl FnMut(&ProcessingStats),
    ) -> Result<String> {
        anyhow::ensure!(self.prepared.is_none(), "Background must precede samples");
        self.group()?;
        let (k, s) = (self.options.kmer_length, self.options.smer_length);
        let (union, found) = self
            .background
            .get_or_insert_with(|| (Arc::new(build_union(&self.targets, k)), KmerSet::new(k)));
        let (counts, stats) = count_stream(Input::File(stream), union, k, s, None, on_batch)?;
        let summary = counted_summary("Background", &stats, counts.len(), union.len());
        found.extend(&counts.into_kmers());
        self.flags |= QUERY_INDEX_FLAG_BACKGROUND;
        Ok(summary)
    }

    /// Group targets and drop background k-mers
    fn finish_targets(&mut self) -> Result<()> {
        self.group()?;
        if let Some((_, found)) = self.background.take() {
            mask_targets(&mut self.targets, &found);
        }
        Ok(())
    }

    /// Query index of the targets, as `skope index build-query` writes. Precedes samples
    pub fn index(&mut self) -> Result<Vec<u8>> {
        anyhow::ensure!(
            self.prepared.is_none(),
            "The index must be written before samples are added"
        );
        self.finish_targets()?;
        let mut bytes = Vec::new();
        write_query_index(
            &mut bytes,
            &self.targets,
            self.options.kmer_length,
            self.options.smer_length,
            self.options.fraction,
            self.options.complexity,
            self.flags,
        )?;
        Ok(bytes)
    }

    /// Ready targets for counting, applying `--discriminatory`. Returns the targets summary
    pub fn prepare(&mut self) -> Result<String> {
        if self.prepared.is_none() {
            self.finish_targets()?;
            let k = self.options.kmer_length;
            let discriminatory = self.options.discriminatory;
            let (shared, distinct) = share_kmers(&mut self.targets, k, discriminatory);
            let mut summary = targets_summary(&self.targets, discriminatory, shared, distinct);
            if let Some(warning) = empty_targets_warning(&self.targets) {
                summary = format!("{summary}\n{warning}");
            }
            self.prepared = Some((Arc::new(build_union(&self.targets, k)), summary));
        }
        Ok(self
            .prepared
            .as_ref()
            .map(|(_, summary)| summary.clone())
            .unwrap_or_default())
    }

    /// Count a sample's target k-mers, preparing targets first. Returns its summary
    pub fn add_sample<R: Read>(
        &mut self,
        name: &str,
        input: Input<R>,
        on_batch: impl FnMut(&ProcessingStats),
    ) -> Result<String> {
        anyhow::ensure!(
            self.samples.iter().all(|sample| sample.sample_name != name),
            "Duplicate sample name '{name}'"
        );
        self.prepare()?;
        let union = Arc::clone(&self.prepared.as_ref().unwrap().0);
        let QueryOptions {
            kmer_length: k,
            smer_length: s,
            confidence,
            sort_order,
            limit_bp,
            ..
        } = self.options;
        let (counts, stats) = count_stream(input, &union, k, s, limit_bp, on_batch)?;
        let summary = counted_summary("Sample", &stats, counts.len(), union.len());
        let mut results = calculate_containment_statistics(
            &self.targets,
            &counts,
            &self.thresholds,
            k,
            confidence,
        );
        sort_results(&mut results, sort_order);
        self.samples.push(SampleResults::new(
            name.to_string(),
            results,
            &self.thresholds,
            stats.total_seqs,
            stats.total_bp,
        ));
        Ok(summary)
    }

    /// TSV of every sample so far, as `skope query` writes
    pub fn tsv(&self) -> Result<String> {
        let mut tsv = Vec::new();
        output_tsv(
            &mut tsv,
            &self.samples,
            &self.thresholds,
            self.options.confidence,
            self.options.no_total,
        )?;
        Ok(String::from_utf8(tsv)?)
    }
}

/// Count a stream's target k-mers on one thread
fn count_stream(
    input: Input<impl Read>,
    union: &Arc<KmerSet>,
    kmer_length: u8,
    smer_length: u8,
    limit_bp: Option<u64>,
    on_batch: impl FnMut(&ProcessingStats),
) -> Result<(AbundanceMap, ProcessingStats)> {
    let counts = Arc::new(Mutex::new(AbundanceMap::new(kmer_length)));
    let mut processor = CountProcessor::new(
        kmer_length,
        smer_length,
        Arc::clone(union),
        Arc::clone(&counts),
    );
    let stats = process_stream(input, &mut processor, limit_bp, on_batch)?;
    drop(processor);
    Ok((Arc::into_inner(counts).unwrap().into_inner(), stats))
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

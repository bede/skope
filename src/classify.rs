use crate::kmers::{Kdust, KmerVec, Kmers};
use crate::{
    FixedRapidHasher, IndexKind, Layout, ProcessingStats, Progress, RapidHashSet, SeqProcessor,
    StdinTargets, TargetGroup, TargetSource, check_index_complexity, complexity_info_line,
    format_bp, format_bp_per_sec, index_writer, output_writer, process_input, reader_for_path,
    resolve_targets, sample_inputs, validate_index_output,
};
use anyhow::{Context, Result};
use paraseq::Record;
use paraseq::parallel::{ParallelProcessor, ParallelReader};
use parking_lot::Mutex;
use std::collections::HashMap;
use std::fs::{self, File};
use std::hash::Hash;
use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

// Index header constants and metadata
use crate::INDEX_MAGIC;
const CLASSIFICATION_INDEX_VERSION: u8 = 1;
/// Maximum groups representable by a u128 bitmask
pub const MAX_GROUPS: usize = 128;
/// Marker for `--individual` group-cap errors crossing the paraseq boundary
const TOO_MANY_RECORDS_MSG: &str = "Too many records for --individual";
// magic, kind, version, k, s, num_groups, complexity (kdust)
type ClassificationIndexHeader = ([u8; 4], u8, u8, u8, u8, u8, f32);

/// Classification index mapping k-mers to group bitmasks (up to 128 groups)
#[derive(Clone)]
pub enum ClassificationIndex {
    U64(HashMap<u64, u128, FixedRapidHasher>),
    U128(HashMap<u128, u128, FixedRapidHasher>),
}

impl ClassificationIndex {
    fn new(kmer_length: u8) -> Self {
        if kmer_length <= 32 {
            Self::U64(HashMap::with_hasher(FixedRapidHasher))
        } else {
            Self::U128(HashMap::with_hasher(FixedRapidHasher))
        }
    }

    pub fn len(&self) -> usize {
        match self {
            ClassificationIndex::U64(m) => m.len(),
            ClassificationIndex::U128(m) => m.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Remove k-mers shared across groups, keeping only group-unique k-mers
/// Returns how many shared k-mers were removed
pub(crate) fn apply_discriminatory_filter(index: &mut ClassificationIndex) -> usize {
    match index {
        ClassificationIndex::U64(map) => {
            let before = map.len();
            map.retain(|_, bitmask| bitmask.count_ones() == 1);
            before - map.len()
        }
        ClassificationIndex::U128(map) => {
            let before = map.len();
            map.retain(|_, bitmask| bitmask.count_ones() == 1);
            before - map.len()
        }
    }
}

// Classification result types
#[derive(Debug, Clone, Copy)]
pub(crate) enum Classification {
    Unclassified,
    Classified(usize),
    Ambiguous(u128),
}

#[derive(Debug, Clone, Default)]
struct GroupCounts {
    seqs: u64,
    bases: u64,
}

#[derive(Debug, Clone, Default)]
struct SampleClassificationResult {
    counts: ClassifyCounts,
    total_seqs: u64,
    total_bases: u64,
}

// Configuration structs
pub struct BuildClassifyConfig {
    pub targets_path: PathBuf,
    pub individual: bool,
    pub kmer_length: u8,
    pub smer_length: u8,
    pub complexity: f32,
    pub threads: usize,
    pub output_path: Option<PathBuf>,
    pub quiet: bool,
}

pub struct ClassifyConfig {
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
    pub threads: usize,
    pub limit_bp: Option<u64>,
    pub output_path: Option<PathBuf>,
    pub per_seq: bool,
    pub discriminatory: bool,
    pub quiet: bool,
}

/// Collect k-mers from one group FASTA file
#[derive(Clone)]
struct GroupKmerProcessor {
    kmers: Kmers,
    kdust: Kdust,
    group_bit: u128,

    local_map: ClassificationIndex,
    global_map: Arc<Mutex<ClassificationIndex>>,
    progress: Progress,

    /// Source name for `--individual` group-cap errors
    source: String,

    /// `--individual`: one group per record. Each record's bit is its index in this vec,
    /// so the run must be single-threaded for the order to be deterministic.
    individual_names: Option<Arc<Mutex<Vec<String>>>>,
}

impl GroupKmerProcessor {
    fn new(
        kmer_length: u8,
        smer_length: u8,
        kdust: Kdust,
        group_bit: u128,
        global_map: Arc<Mutex<ClassificationIndex>>,
        progress: Progress,
        individual_names: Option<Arc<Mutex<Vec<String>>>>,
        source: String,
    ) -> Self {
        Self {
            kmers: Kmers::new(kmer_length, smer_length),
            kdust,
            group_bit,
            local_map: ClassificationIndex::new(kmer_length),
            global_map,
            progress,
            source,
            individual_names,
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for GroupKmerProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        let group_bit = match &self.individual_names {
            Some(names) => {
                let mut names = names.lock();
                if names.len() >= MAX_GROUPS {
                    return Err(paraseq::Error::Io(std::io::Error::other(format!(
                        "{TOO_MANY_RECORDS_MSG}: {} has more than {MAX_GROUPS} records, \
                             the maximum number of groups",
                        self.source
                    ))));
                }
                names.push(String::from_utf8_lossy(record.id()).to_string());
                1u128 << (names.len() - 1)
            }
            None => self.group_bit,
        };

        let seq = record.seq();
        self.progress.add(1, seq.len() as u64)?;
        let kmers = self.kmers.fill(&seq);
        self.kdust.retain(kmers, None);

        match (&*kmers, &mut self.local_map) {
            (KmerVec::U64(kmers), ClassificationIndex::U64(map)) => {
                merge_group_bits(kmers.iter().map(|&kmer| (kmer, group_bit)), map)
            }
            (KmerVec::U128(kmers), ClassificationIndex::U128(map)) => {
                merge_group_bits(kmers.iter().map(|&kmer| (kmer, group_bit)), map)
            }
            _ => unreachable!("k-mer width does not match index"),
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        match (&mut self.local_map, &mut *self.global_map.lock()) {
            (ClassificationIndex::U64(local), ClassificationIndex::U64(global)) => {
                merge_group_bits(local.drain(), global)
            }
            (ClassificationIndex::U128(local), ClassificationIndex::U128(global)) => {
                merge_group_bits(local.drain(), global)
            }
            _ => unreachable!("k-mer width does not match index"),
        }
        self.progress.flush();
        Ok(())
    }
}

fn merge_group_bits<K: Eq + Hash>(
    entries: impl IntoIterator<Item = (K, u128)>,
    map: &mut HashMap<K, u128, FixedRapidHasher>,
) {
    for (kmer, bits) in entries {
        *map.entry(kmer).or_insert(0) |= bits;
    }
}

fn group_counts<K>(map: &HashMap<K, u128, FixedRapidHasher>, group_bit: u128) -> (usize, usize) {
    map.values().fold((0, 0), |(total, unique), &bits| {
        (
            total + usize::from(bits & group_bit != 0),
            unique + usize::from(bits == group_bit),
        )
    })
}

/// Build an in-memory classification index from target groups
pub(crate) fn build_classification_index(
    groups: &[TargetGroup],
    source: &str,
    individual: bool,
    kmer_length: u8,
    smer_length: u8,
    complexity: f32,
    threads: usize,
    quiet: bool,
) -> Result<(ClassificationIndex, Vec<String>)> {
    if groups.len() > MAX_GROUPS {
        return Err(anyhow::anyhow!(
            "Too many groups: {} (max {MAX_GROUPS}). Each top-level fastx file or subdirectory is one group.",
            groups.len()
        ));
    }

    if !quiet {
        eprintln!(
            "Groups: {} (from {})",
            if individual {
                "one per record".to_string()
            } else {
                groups.len().to_string()
            },
            source
        );
    }

    let global_map = Arc::new(Mutex::new(ClassificationIndex::new(kmer_length)));

    let kdust = Kdust::from_threshold(complexity, kmer_length);

    if individual {
        let names = build_individual_groups(
            groups,
            kmer_length,
            smer_length,
            kdust,
            Arc::clone(&global_map),
            quiet,
        )?;
        return Ok((finish_index(global_map, quiet), names));
    }

    let group_names: Vec<String> = groups.iter().map(|g| g.name.clone()).collect();

    for (group_idx, group) in groups.iter().enumerate() {
        let group_bit = 1u128 << group_idx;
        let progress = Progress::new("Collecting group k-mers", quiet, None)?;

        for group_file in &group.files {
            let mut processor = GroupKmerProcessor::new(
                kmer_length,
                smer_length,
                kdust,
                group_bit,
                Arc::clone(&global_map),
                progress.clone(),
                None,
                String::new(),
            );

            let reader = reader_for_path(group_file)?;
            reader.process_parallel(&mut processor, threads)?;
        }

        let stats = progress.finish();

        if !quiet {
            let (group_kmers, unique_kmers) = match &*global_map.lock() {
                ClassificationIndex::U64(map) => group_counts(map, group_bit),
                ClassificationIndex::U128(map) => group_counts(map, group_bit),
            };
            eprintln!(
                "  [{}] {} ({} file{}): {} seqs ({}), {} k-mers ({} unique)",
                group_idx,
                group.name,
                group.files.len(),
                if group.files.len() == 1 { "" } else { "s" },
                stats.total_seqs,
                format_bp(stats.total_bp as usize),
                group_kmers,
                unique_kmers,
            );
        }
    }

    Ok((finish_index(global_map, quiet), group_names))
}

/// Unwrap accumulated k-mer maps into a classification index
fn finish_index(global_map: Arc<Mutex<ClassificationIndex>>, quiet: bool) -> ClassificationIndex {
    let index = Arc::into_inner(global_map).unwrap().into_inner();

    if !quiet {
        let shared = match &index {
            ClassificationIndex::U64(m) => m.values().filter(|v| v.count_ones() > 1).count(),
            ClassificationIndex::U128(m) => m.values().filter(|v| v.count_ones() > 1).count(),
        };
        eprintln!(
            "Index: {} total k-mers, {} shared across groups",
            index.len(),
            shared
        );
    }

    index
}

/// Index a single fastx file one group per record (`--individual`), returning the record
/// names in file order. Single-threaded so those names match the bit indices.
fn build_individual_groups(
    groups: &[TargetGroup],
    kmer_length: u8,
    smer_length: u8,
    kdust: Kdust,
    global_map: Arc<Mutex<ClassificationIndex>>,
    quiet: bool,
) -> Result<Vec<String>> {
    let [group] = groups else {
        return Err(anyhow::anyhow!(
            "--individual expects a single fastx file, got {} groups",
            groups.len()
        ));
    };
    let [file] = group.files.as_slice() else {
        return Err(anyhow::anyhow!(
            "--individual expects a single fastx file, got {} files",
            group.files.len()
        ));
    };

    let names = Arc::new(Mutex::new(Vec::new()));
    let mut processor = GroupKmerProcessor::new(
        kmer_length,
        smer_length,
        kdust,
        0,
        global_map,
        Progress::new("Collecting record k-mers", quiet, None)?,
        Some(Arc::clone(&names)),
        group.name.clone(),
    );

    let reader = reader_for_path(file)?;
    // Surface the group-cap error on its own rather than wrapped as an I/O failure
    if let Err(err) = reader.process_parallel(&mut processor, 1) {
        return Err(match &err {
            paraseq::Error::Io(io_err) if io_err.to_string().starts_with(TOO_MANY_RECORDS_MSG) => {
                anyhow::anyhow!("{io_err}")
            }
            _ => err.into(),
        });
    }
    let stats = processor.progress.finish();
    drop(processor);

    let names = Arc::try_unwrap(names).unwrap().into_inner();
    if !quiet {
        eprintln!(
            "  {} records from {} ({})",
            names.len(),
            group.name,
            format_bp(stats.total_bp as usize),
        );
    }

    Ok(names)
}

pub fn run_build_classify(config: &BuildClassifyConfig) -> Result<()> {
    validate_index_output(config.output_path.as_deref())?;
    let start_time = Instant::now();
    let version = env!("CARGO_PKG_VERSION");

    let mut options = String::new();
    if config.complexity > 0.0 {
        options.push_str(&format!(", complexity={}", config.complexity));
    }
    eprintln!(
        "Skope v{}; mode: index build; options: k={}, s={}, threads={}{}",
        version, config.kmer_length, config.smer_length, config.threads, options
    );

    let source = resolve_targets(
        &config.targets_path,
        IndexKind::Classify,
        StdinTargets::Accept,
    )?;
    if let TargetSource::Index(path) = &source {
        return Err(anyhow::anyhow!(
            "{} is already a skope classification index",
            path.display()
        ));
    }

    let (index, group_names) = build_classification_index(
        source.groups(),
        &config.targets_path.to_string_lossy(),
        source.splits_records(config.individual),
        config.kmer_length,
        config.smer_length,
        config.complexity,
        config.threads,
        config.quiet,
    )?;

    save_index(
        &index,
        &group_names,
        config.kmer_length,
        config.smer_length,
        config.complexity,
        config.output_path.as_deref(),
    )?;

    if !config.quiet {
        let elapsed = start_time.elapsed();
        eprintln!("Done in {:.1}s", elapsed.as_secs_f64());
    }

    Ok(())
}

// Index serialization
fn save_index(
    index: &ClassificationIndex,
    group_names: &[String],
    kmer_length: u8,
    smer_length: u8,
    complexity: f32,
    output_path: Option<&Path>,
) -> Result<()> {
    let mut writer = index_writer(output_path)?;

    let header: ClassificationIndexHeader = (
        *INDEX_MAGIC,
        IndexKind::Classify as u8,
        CLASSIFICATION_INDEX_VERSION,
        kmer_length,
        smer_length,
        group_names.len() as u8,
        complexity,
    );

    let header_bytes = wincode::serialize(&header).context("Failed to encode index header")?;
    writer.write_all(&header_bytes)?;

    let names_bytes = wincode::serialize(&group_names).context("Failed to encode group names")?;
    writer.write_all(&names_bytes)?;

    let count = index.len() as u64;
    let count_bytes = wincode::serialize(&count).context("Failed to encode entry count")?;
    writer.write_all(&count_bytes)?;

    let kmer_bytes = (kmer_length as usize).div_ceil(4); // ceil(k / 4)
    match index {
        ClassificationIndex::U64(map) => {
            for (&kmer, &bitmask) in map {
                let kmer_le = kmer.to_le_bytes();
                writer.write_all(&kmer_le[..kmer_bytes])?;
                writer.write_all(&bitmask.to_le_bytes())?;
            }
        }
        ClassificationIndex::U128(map) => {
            for (&kmer, &bitmask) in map {
                let kmer_le = kmer.to_le_bytes();
                writer.write_all(&kmer_le[..kmer_bytes])?;
                writer.write_all(&bitmask.to_le_bytes())?;
            }
        }
    }

    writer.flush()?;
    Ok(())
}

/// Print human-readable metadata for a classification index (`skope index info`)
pub fn print_classification_index_info(path: &Path) -> Result<()> {
    let mut reader = BufReader::new(File::open(path)?);
    let (magic, kind, version, kmer_length, smer_length, num_groups, complexity): ClassificationIndexHeader =
        wincode::deserialize_from(&mut reader).context("Failed to decode index header")?;
    if &magic != INDEX_MAGIC || IndexKind::from_byte(kind) != Some(IndexKind::Classify) {
        return Err(anyhow::anyhow!(
            "{} is not a skope classification index",
            path.display()
        ));
    }
    if version != CLASSIFICATION_INDEX_VERSION {
        return Err(anyhow::anyhow!(
            "Unsupported index format version: {version} (expected {CLASSIFICATION_INDEX_VERSION})"
        ));
    }
    let group_names: Vec<String> =
        wincode::deserialize_from(&mut reader).context("Failed to decode group names")?;
    if group_names.len() != num_groups as usize {
        return Err(anyhow::anyhow!(
            "Group count mismatch: header says {num_groups} but found {} names",
            group_names.len()
        ));
    }
    let count: u64 =
        wincode::deserialize_from(&mut reader).context("Failed to decode entry count")?;

    println!("Index information:");
    println!("  Format: classify");
    println!("  Format version: {version}");
    println!("  K-mer length (k): {kmer_length}");
    println!("  S-mer length (s): {smer_length}");
    println!("  Groups: {num_groups}");
    println!("  Distinct k-mers: {count}");
    println!("{}", complexity_info_line(complexity));
    for name in &group_names {
        println!("    - {name}");
    }
    Ok(())
}

pub fn load_classification_index(
    path: &Path,
) -> Result<(ClassificationIndex, Vec<String>, u8, u8, f32)> {
    let file_bytes =
        fs::read(path).with_context(|| format!("Failed to open index file: {}", path.display()))?;
    let mut cursor = wincode::io::Cursor::new(file_bytes.as_slice());

    let header: ClassificationIndexHeader =
        wincode::deserialize_from(&mut cursor).context("Failed to decode index header")?;
    let (magic, kind, format_version, kmer_length, smer_length, num_groups, complexity) = header;

    if &magic != INDEX_MAGIC || IndexKind::from_byte(kind) != Some(IndexKind::Classify) {
        return Err(anyhow::anyhow!(
            "{} is not a skope classification index",
            path.display()
        ));
    }

    if format_version != CLASSIFICATION_INDEX_VERSION {
        return Err(anyhow::anyhow!(
            "Unsupported index format version: {} (expected {})",
            format_version,
            CLASSIFICATION_INDEX_VERSION
        ));
    }

    let group_names: Vec<String> =
        wincode::deserialize_from(&mut cursor).context("Failed to decode group names")?;

    if group_names.len() != num_groups as usize {
        return Err(anyhow::anyhow!(
            "Group count mismatch: header says {} but found {} names",
            num_groups,
            group_names.len()
        ));
    }

    let count_u64: u64 =
        wincode::deserialize_from(&mut cursor).context("Failed to decode entry count")?;
    let count = usize::try_from(count_u64)
        .with_context(|| format!("Entry count is too large for this platform: {count_u64}"))?;

    let kmer_bytes = (kmer_length as usize).div_ceil(4);
    let entry_size = kmer_bytes + 16; // k-mer bytes + group bitmask (u128)

    let raw_data = &file_bytes[cursor.position()..];

    let expected_size = count * entry_size;
    if raw_data.len() < expected_size {
        return Err(anyhow::anyhow!(
            "Index file truncated: expected {} bytes of entries, got {}",
            expected_size,
            raw_data.len()
        ));
    }

    let index = if kmer_length <= 32 {
        let mut map: HashMap<u64, u128, FixedRapidHasher> =
            HashMap::with_capacity_and_hasher(count, FixedRapidHasher);
        for i in 0..count {
            let offset = i * entry_size;
            let mut kmer_buf = [0u8; 8];
            kmer_buf[..kmer_bytes].copy_from_slice(&raw_data[offset..offset + kmer_bytes]);
            let kmer = u64::from_le_bytes(kmer_buf);

            let bitmask_offset = offset + kmer_bytes;
            let bitmask = u128::from_le_bytes(
                raw_data[bitmask_offset..bitmask_offset + 16]
                    .try_into()
                    .unwrap(),
            );

            map.insert(kmer, bitmask);
        }
        ClassificationIndex::U64(map)
    } else {
        let mut map: HashMap<u128, u128, FixedRapidHasher> =
            HashMap::with_capacity_and_hasher(count, FixedRapidHasher);
        for i in 0..count {
            let offset = i * entry_size;
            let mut kmer_buf = [0u8; 16];
            kmer_buf[..kmer_bytes].copy_from_slice(&raw_data[offset..offset + kmer_bytes]);
            let kmer = u128::from_le_bytes(kmer_buf);

            let bitmask_offset = offset + kmer_bytes;
            let bitmask = u128::from_le_bytes(
                raw_data[bitmask_offset..bitmask_offset + 16]
                    .try_into()
                    .unwrap(),
            );

            map.insert(kmer, bitmask);
        }
        ClassificationIndex::U128(map)
    };

    Ok((index, group_names, kmer_length, smer_length, complexity))
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Thresholds {
    pub abs: u64,
    pub rel: f64,
}

impl Thresholds {
    /// Minimum hits, nondecreasing with distinct k-mer count
    pub(crate) fn required_hits(&self, kmers: usize) -> u64 {
        // ceil(0.1 * 30) overshoots to 4
        let kmers = kmers.max(1) as f64;
        let lower = (self.rel * kmers) as u64;
        let rel = lower + u64::from((lower as f64) / kmers < self.rel);
        self.abs.max(rel).max(1)
    }
}

/// Classify pooled k-mers by distinct hits
#[derive(Clone)]
pub(crate) struct Classifier {
    index: Arc<ClassificationIndex>,
    pub num_groups: usize,
    thresholds: Thresholds,
    /// Count all distinct hits for --per-seq
    exact: bool,
    kmers: Kmers,
    seen_u64: RapidHashSet<u64>,
    seen_u128: RapidHashSet<u128>,
    /// Complete and distinct only in exact mode
    pub hits: [u64; MAX_GROUPS],
    /// Distinct count in exact mode, otherwise an upper bound
    pub kmer_count: usize,
}

impl Classifier {
    pub(crate) fn new(
        index: Arc<ClassificationIndex>,
        num_groups: usize,
        kmer_length: u8,
        smer_length: u8,
        thresholds: Thresholds,
        exact: bool,
    ) -> Self {
        Self {
            index,
            num_groups,
            thresholds,
            exact,
            kmers: Kmers::new(kmer_length, smer_length),
            seen_u64: RapidHashSet::default(),
            seen_u128: RapidHashSet::default(),
            hits: [0; MAX_GROUPS],
            kmer_count: 0,
        }
    }

    pub(crate) fn classify(&mut self, seqs: &[&[u8]]) -> Classification {
        let kmers = self.kmers.pool(seqs);
        if self.exact {
            kmers.sort_dedup();
        }
        let positions = kmers.len();
        let hits = &mut self.hits[..self.num_groups];
        hits.fill(0);

        // Positions bound the distinct count from above
        let pass = self.thresholds.required_hits(positions);
        // Stop at one match for a single group, two for ambiguity
        let stop = if self.exact {
            usize::MAX
        } else {
            self.num_groups.min(2)
        };
        let dedup = !self.exact && pass > 1;
        let hit_kmers = match (&*kmers, &*self.index) {
            (KmerVec::U64(kmers), ClassificationIndex::U64(map)) => {
                count_hits(kmers, map, &mut self.seen_u64, dedup, hits, pass, stop)
            }
            (KmerVec::U128(kmers), ClassificationIndex::U128(map)) => {
                count_hits(kmers, map, &mut self.seen_u128, dedup, hits, pass, stop)
            }
            _ => unreachable!("k-mer width does not match index"),
        };

        // Hit k-mers bound the distinct count from below
        // Resolve groups between the lower and upper requirements by deduping
        let undecided = dedup
            && hit_kmers.is_some_and(|hit_kmers| {
                let floor = self.thresholds.required_hits(hit_kmers);
                hits.iter().any(|h| (floor..pass).contains(h))
            });
        self.kmer_count = if undecided {
            kmers.sort_dedup()
        } else {
            positions
        };

        let required = self.thresholds.required_hits(self.kmer_count);
        let groups = hits
            .iter()
            .enumerate()
            .filter(|&(_, &h)| h >= required)
            .fold(0u128, |groups, (group_idx, _)| groups | 1 << group_idx);
        match groups.count_ones() {
            0 => Classification::Unclassified,
            1 => Classification::Classified(groups.trailing_zeros() as usize),
            _ => Classification::Ambiguous(groups),
        }
    }
}

/// Count hit k-mers, returning None after `stop` groups reach `pass`
fn count_hits<T: Copy + Eq + Hash>(
    kmers: &[T],
    index: &HashMap<T, u128, FixedRapidHasher>,
    seen: &mut RapidHashSet<T>,
    dedup: bool,
    hits: &mut [u64],
    pass: u64,
    stop: usize,
) -> Option<usize> {
    seen.clear();
    let (mut hit_kmers, mut passed) = (0, 0);
    for kmer in kmers {
        if let Some(&mask) = index.get(kmer)
            && (!dedup || seen.insert(*kmer))
        {
            hit_kmers += 1;
            for group_idx in set_bits(mask) {
                hits[group_idx] += 1;
                passed += usize::from(hits[group_idx] == pass);
            }
            if passed >= stop {
                return None;
            }
        }
    }
    Some(hit_kmers)
}

/// Group indices set in a bitmask, ascending
fn set_bits(mut mask: u128) -> impl Iterator<Item = usize> {
    std::iter::from_fn(move || {
        (mask != 0).then(|| {
            let group_idx = mask.trailing_zeros() as usize;
            mask &= mask - 1;
            group_idx
        })
    })
}

#[derive(Debug, Clone, Default)]
struct ClassifyCounts {
    groups: Vec<GroupCounts>,
    ambiguous_seqs: u64,
    ambiguous_bases: u64,
    unclassified_seqs: u64,
    unclassified_bases: u64,
}

impl ClassifyCounts {
    fn new(num_groups: usize) -> Self {
        Self {
            groups: vec![GroupCounts::default(); num_groups],
            ambiguous_seqs: 0,
            ambiguous_bases: 0,
            unclassified_seqs: 0,
            unclassified_bases: 0,
        }
    }
    fn add(&mut self, classification: Classification, seqs: u64, bases: u64) {
        match classification {
            Classification::Classified(group_idx) => {
                self.groups[group_idx].seqs += seqs;
                self.groups[group_idx].bases += bases;
            }
            Classification::Ambiguous(_) => {
                self.ambiguous_seqs += seqs;
                self.ambiguous_bases += bases;
            }
            Classification::Unclassified => {
                self.unclassified_seqs += seqs;
                self.unclassified_bases += bases;
            }
        }
    }

    fn merge(&mut self, other: &mut Self) {
        for (total, local) in self.groups.iter_mut().zip(&mut other.groups) {
            total.seqs += local.seqs;
            total.bases += local.bases;
            *local = GroupCounts::default();
        }
        self.ambiguous_seqs += std::mem::take(&mut other.ambiguous_seqs);
        self.ambiguous_bases += std::mem::take(&mut other.ambiguous_bases);
        self.unclassified_seqs += std::mem::take(&mut other.unclassified_seqs);
        self.unclassified_bases += std::mem::take(&mut other.unclassified_bases);
    }
}

#[derive(Clone)]
enum ClassifyOutput {
    Summary {
        local: ClassifyCounts,
        global: Arc<Mutex<ClassifyCounts>>,
    },
    PerSeq {
        group_names: Arc<Vec<String>>,
        sample_name: String,
        local: Vec<u8>,
        writer: Arc<Mutex<Box<dyn Write + Send>>>,
    },
}

impl ClassifyOutput {
    fn add_record(
        &mut self,
        seq_id: &[u8],
        seqs: u64,
        seq_len: u64,
        total_kmers: usize,
        classification: Classification,
        hits: &[u64; 128],
    ) {
        let (group_names, sample_name, local) = match self {
            Self::Summary { local, .. } => {
                local.add(classification, seqs, seq_len);
                return;
            }
            Self::PerSeq {
                group_names,
                sample_name,
                local,
                ..
            } => (group_names, sample_name, local),
        };

        let seq_id = String::from_utf8_lossy(seq_id);
        match classification {
            Classification::Classified(group_idx) => {
                let _ = writeln!(
                    local,
                    "{sample_name}\t{seq_id}\tclassified\t{}\t{}\t{total_kmers}\t{seq_len}",
                    group_names[group_idx], hits[group_idx],
                );
            }
            Classification::Unclassified => {
                let _ = writeln!(
                    local,
                    "{sample_name}\t{seq_id}\tunclassified\t.\t0\t{total_kmers}\t{seq_len}",
                );
            }
            Classification::Ambiguous(mask) => {
                let groups = set_bits(mask)
                    .map(|i| group_names[i].as_str())
                    .collect::<Vec<_>>()
                    .join(",");
                let hits_str = set_bits(mask)
                    .map(|i| hits[i].to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                let _ = writeln!(
                    local,
                    "{sample_name}\t{seq_id}\tambiguous\t{groups}\t{hits_str}\t{total_kmers}\t{seq_len}",
                );
            }
        }
    }

    fn flush(&mut self) -> paraseq::Result<()> {
        match self {
            Self::Summary { local, global } => global.lock().merge(local),
            Self::PerSeq { local, writer, .. } if !local.is_empty() => {
                writer.lock().write_all(local).map_err(paraseq::Error::Io)?;
                local.clear();
            }
            Self::PerSeq { .. } => {}
        }
        Ok(())
    }
}

#[derive(Clone)]
struct ClassifyProcessor {
    classifier: Classifier,
    output: ClassifyOutput,
}

impl SeqProcessor for ClassifyProcessor {
    fn process(&mut self, seq_id: &[u8], seqs: &[&[u8]]) -> paraseq::Result<()> {
        let bases = seqs.iter().map(|seq| seq.len() as u64).sum();
        let classification = self.classifier.classify(seqs);
        self.output.add_record(
            seq_id,
            seqs.len() as u64,
            bases,
            self.classifier.kmer_count,
            classification,
            &self.classifier.hits,
        );
        Ok(())
    }

    fn flush(&mut self) -> paraseq::Result<()> {
        self.output.flush()
    }
}

pub fn run_classification(config: &ClassifyConfig) -> Result<()> {
    let start_time = Instant::now();
    let version = env!("CARGO_PKG_VERSION");

    let source = resolve_targets(
        &config.targets_path,
        IndexKind::Classify,
        StdinTargets::Reject,
    )?;

    let (index, group_names, kmer_length, smer_length) = if !matches!(
        source,
        TargetSource::Index(_)
    ) {
        let mut options = config.layout.option_label().to_string();
        if config.complexity > 0.0 {
            options.push_str(&format!(", complexity={}", config.complexity));
        }
        if let Some(limit) = config.limit_bp {
            options.push_str(&format!(", limit_bp={}", limit));
        }
        eprintln!(
            "Skope v{}; mode: classify (from {}); options: k={}, s={}, threads={}, abs_threshold={}, rel_threshold={}{}",
            version,
            if matches!(source, TargetSource::Directory(_)) {
                "directory"
            } else {
                "fastx"
            },
            config.kmer_length,
            config.smer_length,
            config.threads,
            config.abs_threshold,
            config.rel_threshold,
            options
        );

        let (index, group_names) = build_classification_index(
            source.groups(),
            &config.targets_path.to_string_lossy(),
            source.splits_records(config.individual),
            config.kmer_length,
            config.smer_length,
            config.complexity,
            config.threads,
            config.quiet,
        )?;

        (index, group_names, config.kmer_length, config.smer_length)
    } else {
        let limit_str = config
            .limit_bp
            .map_or(String::new(), |v| format!(", limit_bp={}", v));
        eprintln!(
            "Skope v{}; mode: classify (from index); options: threads={}, abs_threshold={}, rel_threshold={}{}{}",
            version,
            config.threads,
            config.abs_threshold,
            config.rel_threshold,
            config.layout.option_label(),
            limit_str
        );

        let load_start = Instant::now();
        let (index, group_names, k, s, index_complexity) =
            load_classification_index(&config.targets_path)?;
        check_index_complexity(config.complexity, index_complexity)?;

        if !config.quiet {
            let elapsed = load_start.elapsed();
            eprintln!(
                "Index: {} k-mers, {} groups, k={}, s={} (loaded in {:.1}s)",
                index.len(),
                group_names.len(),
                k,
                s,
                elapsed.as_secs_f64()
            );
            for (i, name) in group_names.iter().enumerate() {
                eprintln!("  [{}] {}", i, name);
            }
        }

        (index, group_names, k, s)
    };

    let mut index = index;
    if config.discriminatory {
        let removed = apply_discriminatory_filter(&mut index);
        if !config.quiet {
            eprintln!(
                "Discriminatory mode: removed {} shared k-mers, {} unique k-mers remain",
                removed,
                index.len()
            );
        }
    }

    let group_names = Arc::new(group_names);
    let num_groups = group_names.len();
    let classifier = Classifier::new(
        Arc::new(index),
        num_groups,
        kmer_length,
        smer_length,
        Thresholds {
            abs: config.abs_threshold,
            rel: config.rel_threshold,
        },
        config.per_seq,
    );

    use rayon::prelude::*;
    if config.per_seq {
        let writer = Arc::new(Mutex::new(output_writer(config.output_path.as_deref())?));

        {
            let mut w = writer.lock();
            writeln!(
                w,
                "sample\tseq_id\tclassification\tgroup\thits\tseq_kmers\tseq_length"
            )?;
        }

        for (sample_paths, sample_name) in config.sample_paths.iter().zip(&config.sample_names) {
            process_sample_files(
                sample_paths,
                config.layout,
                sample_name,
                &classifier,
                config.threads,
                config.quiet,
                config.limit_bp,
                || ClassifyOutput::PerSeq {
                    group_names: Arc::clone(&group_names),
                    sample_name: sample_name.clone(),
                    local: Vec::new(),
                    writer: Arc::clone(&writer),
                },
            )?;
        }

        writer.lock().flush()?;
    } else {
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

        let sample_results: Vec<(String, SampleClassificationResult)> = config
            .sample_paths
            .par_iter()
            .zip(&config.sample_names)
            .map(|(sample_paths, sample_name)| {
                let counts = Arc::new(Mutex::new(ClassifyCounts::new(num_groups)));
                let result = process_sample_files(
                    sample_paths,
                    config.layout,
                    sample_name,
                    &classifier,
                    config.threads,
                    config.quiet || is_multisample,
                    config.limit_bp,
                    || ClassifyOutput::Summary {
                        local: ClassifyCounts::new(num_groups),
                        global: Arc::clone(&counts),
                    },
                )
                .map(|stats| SampleClassificationResult {
                    counts: counts.lock().clone(),
                    total_seqs: stats.total_seqs,
                    total_bases: stats.total_bp,
                });

                if let Some(ref counter) = completed {
                    let mut count = counter.lock();
                    *count += 1;
                    eprint!(
                        "\rSamples: processed {} of {}…",
                        *count,
                        config.sample_paths.len()
                    );
                }

                result.map(|r| (sample_name.clone(), r))
            })
            .collect::<Result<Vec<_>>>()?;

        if is_multisample && !config.quiet {
            eprintln!();
        }

        let mut writer = output_writer(config.output_path.as_deref())?;

        writeln!(writer, "sample\tgroup\tseqs_pct\tseqs\tbases_pct\tbases")?;

        for (sample_name, result) in &sample_results {
            let total_seqs = result.total_seqs as f64;
            let total_bases = result.total_bases as f64;

            let mut rows: Vec<(&str, u64, f64, u64, f64)> = Vec::new();
            for (group_idx, counts) in result.counts.groups.iter().enumerate() {
                let pct_seqs = if total_seqs > 0.0 {
                    counts.seqs as f64 / total_seqs * 100.0
                } else {
                    0.0
                };
                let pct_bases = if total_bases > 0.0 {
                    counts.bases as f64 / total_bases * 100.0
                } else {
                    0.0
                };
                rows.push((
                    &group_names[group_idx],
                    counts.seqs,
                    pct_seqs,
                    counts.bases,
                    pct_bases,
                ));
            }

            {
                let pct_seqs = if total_seqs > 0.0 {
                    result.counts.ambiguous_seqs as f64 / total_seqs * 100.0
                } else {
                    0.0
                };
                let pct_bases = if total_bases > 0.0 {
                    result.counts.ambiguous_bases as f64 / total_bases * 100.0
                } else {
                    0.0
                };
                rows.push((
                    "ambiguous",
                    result.counts.ambiguous_seqs,
                    pct_seqs,
                    result.counts.ambiguous_bases,
                    pct_bases,
                ));
            }

            {
                let pct_seqs = if total_seqs > 0.0 {
                    result.counts.unclassified_seqs as f64 / total_seqs * 100.0
                } else {
                    0.0
                };
                let pct_bases = if total_bases > 0.0 {
                    result.counts.unclassified_bases as f64 / total_bases * 100.0
                } else {
                    0.0
                };
                rows.push((
                    "unclassified",
                    result.counts.unclassified_seqs,
                    pct_seqs,
                    result.counts.unclassified_bases,
                    pct_bases,
                ));
            }

            rows.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

            for (group_name, seqs, pct_seqs, bases, pct_bases) in &rows {
                writeln!(
                    writer,
                    "{}\t{}\t{:.3}\t{}\t{:.3}\t{}",
                    sample_name, group_name, pct_seqs, seqs, pct_bases, bases,
                )?;
            }
        }

        writer.flush()?;
    }

    if !config.quiet {
        let elapsed = start_time.elapsed();
        eprintln!("Done in {:.1}s", elapsed.as_secs_f64());
    }

    Ok(())
}

/// Stream each file of one sample through a `ClassifyProcessor`, honouring `limit_bp` across files
#[allow(clippy::too_many_arguments)]
fn process_sample_files(
    sample_paths: &[PathBuf],
    layout: Layout,
    sample_name: &str,
    classifier: &Classifier,
    threads: usize,
    quiet: bool,
    limit_bp: Option<u64>,
    mut make_output: impl FnMut() -> ClassifyOutput,
) -> Result<ProcessingStats> {
    let mut totals = ProcessingStats::default();
    for input in sample_inputs(sample_paths, layout)? {
        if limit_bp.is_some_and(|limit| totals.total_bp >= limit) {
            break;
        }
        let file_start = Instant::now();
        let processor = ClassifyProcessor {
            classifier: classifier.clone(),
            output: make_output(),
        };
        let file_limit = limit_bp.map(|limit| limit.saturating_sub(totals.total_bp));
        let progress = Progress::new("Classifying", quiet, file_limit)?;
        let stats = process_input(input, layout, processor, progress, threads)?;
        totals.total_seqs += stats.total_seqs;
        totals.total_bp += stats.total_bp;
        if !quiet {
            let bp_per_sec = stats.total_bp as f64 / file_start.elapsed().as_secs_f64();
            eprintln!(
                "Sample {}: {} seqs ({}) ({})",
                sample_name,
                stats.total_seqs,
                format_bp(stats.total_bp as usize),
                format_bp_per_sec(bp_per_sec)
            );
        }
    }
    Ok(totals)
}

#[cfg(test)]
mod tests {
    use super::*;

    const K: u8 = 15;
    const S: u8 = 7;
    const ANY_HIT: Thresholds = Thresholds { abs: 1, rel: 0.0 };

    fn pseudo_dna(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed | 1;
        (0..n)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[((x >> 33) & 0b11) as usize]
            })
            .collect()
    }

    fn revcomp(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect()
    }

    fn index_of(groups: &[&[u8]], smer_length: u8) -> Arc<ClassificationIndex> {
        let mut extractor = Kmers::new(K, smer_length);
        let mut map = HashMap::with_hasher(FixedRapidHasher);
        for (group_idx, seq) in groups.iter().enumerate() {
            let KmerVec::U64(kmers) = extractor.fill(seq) else {
                unreachable!()
            };
            for &kmer in kmers.iter() {
                *map.entry(kmer).or_insert(0) |= 1u128 << group_idx;
            }
        }
        Arc::new(ClassificationIndex::U64(map))
    }

    #[test]
    fn test_required_hits_is_fewest_passing_hits() {
        for abs in [0, 1, 3] {
            for rel in [0.0, 0.01, 0.1, 0.25, 0.49, 0.5, 0.9, 0.99, 1.0] {
                let thresholds = Thresholds { abs, rel };
                let mut previous = 0;
                for kmers in 0..2000usize {
                    let n = kmers.max(1) as f64;
                    let fewest = (0..).find(|&h| h as f64 / n >= rel).unwrap();
                    let required = thresholds.required_hits(kmers);
                    assert_eq!(required, abs.max(fewest).max(1), "{abs} {rel} {kmers}");
                    assert!(required >= previous, "decreased at {abs} {rel} {kmers}");
                    previous = required;
                }
            }
        }
    }

    #[test]
    fn test_rel_threshold_is_not_rounded() {
        let required = |rel, kmers| Thresholds { abs: 1, rel }.required_hits(kmers);
        assert_eq!(required(0.49, 5), 3);
        assert_eq!(required(0.1, 30), 3);
    }

    #[test]
    fn test_fast_matches_exact() {
        let (a, b) = (pseudo_dna(400, 1), pseudo_dna(400, 2));
        let reads = [
            a.clone(),
            b.clone(),
            pseudo_dna(400, 3),
            [&a[..], &b[..]].concat(),
            a[..100].repeat(6),
            [&a[..60], &pseudo_dna(340, 4)[..]].concat(),
            a[..10].to_vec(),
            Vec::new(),
        ];
        let (mut substituted, mut stopped) = (false, false);
        for groups in [vec![&a[..]], vec![&a[..], &b[..]]] {
            let index = index_of(&groups, S);
            for abs in [0, 1, 2, 5, 40] {
                for rel in [0.0, 0.05, 0.2, 0.5, 0.9, 1.0] {
                    let thresholds = Thresholds { abs, rel };
                    let classifier = |exact| {
                        Classifier::new(Arc::clone(&index), groups.len(), K, S, thresholds, exact)
                    };
                    let (mut fast, mut exact) = (classifier(false), classifier(true));
                    for r1 in &reads {
                        for r2 in &reads {
                            for seqs in [&[&r1[..]][..], &[&r1[..], &r2[..]]] {
                                let verdicts = (fast.classify(seqs), exact.classify(seqs));
                                match verdicts {
                                    (
                                        Classification::Classified(f),
                                        Classification::Classified(e),
                                    ) if f == e => {}
                                    (
                                        Classification::Ambiguous(_),
                                        Classification::Ambiguous(_),
                                    )
                                    | (
                                        Classification::Unclassified,
                                        Classification::Unclassified,
                                    ) => {}
                                    _ => panic!("{verdicts:?} at abs={abs} rel={rel}"),
                                }
                                substituted |= fast.kmer_count > exact.kmer_count;
                                stopped |= (0..groups.len()).any(|g| fast.hits[g] < exact.hits[g]);
                            }
                        }
                    }
                }
            }
        }
        assert!(substituted, "positions never stood in for distinct k-mers");
        assert!(stopped, "hit counting never stopped early");
    }

    #[test]
    fn test_rel_threshold_uses_distinct_kmers() {
        let a = pseudo_dna(400, 1);
        let repeat = b"ACGT".repeat(100);
        let thresholds = Thresholds { abs: 1, rel: 0.9 };
        for exact in [false, true] {
            let mut classifier = Classifier::new(index_of(&[&a], 0), 1, K, 0, thresholds, exact);
            assert!(matches!(
                classifier.classify(&[&a, &repeat]),
                Classification::Classified(0)
            ));
            assert!(classifier.kmer_count < a.len());
        }
    }

    #[test]
    fn test_zero_hits_never_match() {
        let a = pseudo_dna(400, 1);
        let thresholds = Thresholds { abs: 0, rel: 0.0 };
        let mut classifier = Classifier::new(index_of(&[&a], S), 1, K, S, thresholds, false);
        for seq in [&pseudo_dna(400, 3)[..], &a[..10], b""] {
            assert!(matches!(
                classifier.classify(&[seq]),
                Classification::Unclassified
            ));
        }
    }

    #[test]
    fn test_pooled_mates_count_distinct_kmers() {
        let a = pseudo_dna(400, 1);
        let mut classifier = Classifier::new(index_of(&[&a], S), 1, K, S, ANY_HIT, true);
        classifier.classify(&[&a]);
        let counts = (classifier.hits[0], classifier.kmer_count);
        assert!(counts.0 > 0);
        classifier.classify(&[&a, &revcomp(&a)]);
        assert_eq!((classifier.hits[0], classifier.kmer_count), counts);
    }
}

use crate::syncmers::{Buffers, Kdust, KmerHasher, SyncmerVec, fill_syncmers};
use crate::{
    FixedRapidHasher, IndexKind, ProcessingStats, StdinTargets, TargetGroup, TargetSource,
    check_index_complexity, complexity_info_line, create_spinner, format_bp, format_bp_per_sec,
    handle_process_result, reader_for_path,
    reader_with_inferred_batch_size, resolve_targets, sample_limit_reached_io_error,
};
use anyhow::{Context, Result};
use indicatif::ProgressBar;
use paraseq::Record;
use paraseq::parallel::{ParallelProcessor, ParallelReader};
use parking_lot::Mutex;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Write};
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

/// Classification index mapping syncmers to group bitmasks (up to 128 groups)
#[derive(Clone)]
pub enum ClassificationIndex {
    U64(HashMap<u64, u128, FixedRapidHasher>),
    U128(HashMap<u128, u128, FixedRapidHasher>),
}

impl ClassificationIndex {
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

/// Remove syncmers shared across groups, keeping only group-unique syncmers
/// Returns how many shared syncmers were removed
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

/// Collect syncmers from one group FASTA file
#[derive(Clone)]
struct GroupKmerProcessor {
    kmer_length: u8,
    smer_length: u8,
    hasher: KmerHasher,
    kdust: Kdust,
    buffers: Buffers,
    group_bit: u128,

    // Thread-local k-mer map
    local_map_u64: Option<HashMap<u64, u128, FixedRapidHasher>>,
    local_map_u128: Option<HashMap<u128, u128, FixedRapidHasher>>,

    // Shared global state
    global_map_u64: Arc<Mutex<Option<HashMap<u64, u128, FixedRapidHasher>>>>,
    global_map_u128: Arc<Mutex<Option<HashMap<u128, u128, FixedRapidHasher>>>>,
    local_stats: ProcessingStats,
    global_stats: Arc<Mutex<ProcessingStats>>,

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
        global_map_u64: Arc<Mutex<Option<HashMap<u64, u128, FixedRapidHasher>>>>,
        global_map_u128: Arc<Mutex<Option<HashMap<u128, u128, FixedRapidHasher>>>>,
        global_stats: Arc<Mutex<ProcessingStats>>,
        individual_names: Option<Arc<Mutex<Vec<String>>>>,
        source: String,
    ) -> Self {
        let buffers = if kmer_length <= 32 {
            Buffers::new_u64()
        } else {
            Buffers::new_u128()
        };

        let (local_map_u64, local_map_u128) = if kmer_length <= 32 {
            (Some(HashMap::with_hasher(FixedRapidHasher)), None)
        } else {
            (None, Some(HashMap::with_hasher(FixedRapidHasher)))
        };

        Self {
            kmer_length,
            smer_length,
            hasher: KmerHasher::new(smer_length as usize),
            kdust,
            buffers,
            group_bit,
            local_map_u64,
            local_map_u128,
            global_map_u64,
            global_map_u128,
            local_stats: ProcessingStats::default(),
            global_stats,
            source,
            individual_names,
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for GroupKmerProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        let group_bit = match &self.individual_names {
            Some(names) => {
                let mut names = names.lock();
                if names.len() >= MAX_GROUPS {
                    return Err(paraseq::parallel::ProcessError::IoError(
                        std::io::Error::other(format!(
                            "{TOO_MANY_RECORDS_MSG}: {} has more than {MAX_GROUPS} records, \
                             the maximum number of groups",
                            self.source
                        )),
                    ));
                }
                names.push(String::from_utf8_lossy(record.id()).to_string());
                1u128 << (names.len() - 1)
            }
            None => self.group_bit,
        };

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
        self.kdust.retain(&mut self.buffers.syncmers, None);

        match &self.buffers.syncmers {
            SyncmerVec::U64(vec) => {
                let local = self.local_map_u64.as_mut().unwrap();
                for &kmer in vec {
                    *local.entry(kmer).or_insert(0) |= group_bit;
                }
            }
            SyncmerVec::U128(vec) => {
                let local = self.local_map_u128.as_mut().unwrap();
                for &kmer in vec {
                    *local.entry(kmer).or_insert(0) |= group_bit;
                }
            }
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        if let Some(local) = &mut self.local_map_u64 {
            let mut global = self.global_map_u64.lock();
            let global_map = global.as_mut().unwrap();
            for (&kmer, &bits) in local.iter() {
                *global_map.entry(kmer).or_insert(0) |= bits;
            }
            local.clear();
        } else if let Some(local) = &mut self.local_map_u128 {
            let mut global = self.global_map_u128.lock();
            let global_map = global.as_mut().unwrap();
            for (&kmer, &bits) in local.iter() {
                *global_map.entry(kmer).or_insert(0) |= bits;
            }
            local.clear();
        }

        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;
            self.local_stats = ProcessingStats::default();
        }

        Ok(())
    }
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

    let global_map_u64: Arc<Mutex<Option<HashMap<u64, u128, FixedRapidHasher>>>> =
        if kmer_length <= 32 {
            Arc::new(Mutex::new(Some(HashMap::with_hasher(FixedRapidHasher))))
        } else {
            Arc::new(Mutex::new(None))
        };

    let global_map_u128: Arc<Mutex<Option<HashMap<u128, u128, FixedRapidHasher>>>> =
        if kmer_length > 32 {
            Arc::new(Mutex::new(Some(HashMap::with_hasher(FixedRapidHasher))))
        } else {
            Arc::new(Mutex::new(None))
        };

    let kdust = Kdust::from_threshold(complexity, kmer_length);

    if individual {
        let names = build_individual_groups(
            groups,
            kmer_length,
            smer_length,
            kdust,
            Arc::clone(&global_map_u64),
            Arc::clone(&global_map_u128),
            quiet,
        )?;
        return Ok((
            finish_index(global_map_u64, global_map_u128, kmer_length, quiet),
            names,
        ));
    }

    let group_names: Vec<String> = groups.iter().map(|g| g.name.clone()).collect();

    for (group_idx, group) in groups.iter().enumerate() {
        let group_bit = 1u128 << group_idx;
        let global_stats = Arc::new(Mutex::new(ProcessingStats::default()));

        for group_file in &group.files {
            let mut processor = GroupKmerProcessor::new(
                kmer_length,
                smer_length,
                kdust,
                group_bit,
                Arc::clone(&global_map_u64),
                Arc::clone(&global_map_u128),
                Arc::clone(&global_stats),
                None,
                String::new(),
            );

            let reader = reader_for_path(group_file)?;
            reader.process_parallel(&mut processor, threads)?;
        }

        let stats = global_stats.lock().clone();

        let (group_kmers, unique_kmers) = if kmer_length <= 32 {
            let map = global_map_u64.lock();
            let map = map.as_ref().unwrap();
            let group_kmers = map.values().filter(|&&v| v & group_bit != 0).count();
            let unique = map
                .values()
                .filter(|&&v| v & group_bit != 0 && v.count_ones() == 1)
                .count();
            (group_kmers, unique)
        } else {
            let map = global_map_u128.lock();
            let map = map.as_ref().unwrap();
            let group_kmers = map.values().filter(|&&v| v & group_bit != 0).count();
            let unique = map
                .values()
                .filter(|&&v| v & group_bit != 0 && v.count_ones() == 1)
                .count();
            (group_kmers, unique)
        };

        if !quiet {
            eprintln!(
                "  [{}] {} ({} file{}): {} seqs ({}), {} syncmers ({} unique)",
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

    Ok((
        finish_index(global_map_u64, global_map_u128, kmer_length, quiet),
        group_names,
    ))
}

type GlobalMap<K> = Arc<Mutex<Option<HashMap<K, u128, FixedRapidHasher>>>>;

/// Unwrap accumulated k-mer maps into a classification index
fn finish_index(
    global_map_u64: GlobalMap<u64>,
    global_map_u128: GlobalMap<u128>,
    kmer_length: u8,
    quiet: bool,
) -> ClassificationIndex {
    let index = if kmer_length <= 32 {
        let map = Arc::try_unwrap(global_map_u64)
            .unwrap()
            .into_inner()
            .unwrap();
        ClassificationIndex::U64(map)
    } else {
        let map = Arc::try_unwrap(global_map_u128)
            .unwrap()
            .into_inner()
            .unwrap();
        ClassificationIndex::U128(map)
    };

    if !quiet {
        let shared = match &index {
            ClassificationIndex::U64(m) => m.values().filter(|v| v.count_ones() > 1).count(),
            ClassificationIndex::U128(m) => m.values().filter(|v| v.count_ones() > 1).count(),
        };
        eprintln!(
            "Index: {} total syncmers, {} shared across groups",
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
    global_map_u64: GlobalMap<u64>,
    global_map_u128: GlobalMap<u128>,
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
    let global_stats = Arc::new(Mutex::new(ProcessingStats::default()));
    let mut processor = GroupKmerProcessor::new(
        kmer_length,
        smer_length,
        kdust,
        0,
        global_map_u64,
        global_map_u128,
        Arc::clone(&global_stats),
        Some(Arc::clone(&names)),
        group.name.clone(),
    );

    let reader = reader_for_path(file)?;
    // Surface the group-cap error on its own rather than wrapped as an I/O failure
    if let Err(err) = reader.process_parallel(&mut processor, 1) {
        return Err(match &err {
            paraseq::parallel::ProcessError::IoError(io_err)
                if io_err.to_string().starts_with(TOO_MANY_RECORDS_MSG) =>
            {
                anyhow::anyhow!("{io_err}")
            }
            _ => err.into(),
        });
    }
    drop(processor);

    let names = Arc::try_unwrap(names).unwrap().into_inner();
    if !quiet {
        let stats = global_stats.lock();
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
        config.output_path.as_ref(),
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
    output_path: Option<&PathBuf>,
) -> Result<()> {
    let writer: Box<dyn Write> = if let Some(path) = output_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };
    let mut writer = writer;

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
    println!("  Format: classify (open syncmer set)");
    println!("  Format version: {version}");
    println!("  K-mer length (k): {kmer_length}");
    println!("  S-mer length (s): {smer_length}");
    println!("  Groups: {num_groups}");
    println!("  Distinct syncmers: {count}");
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

/// Classify one sequecne using per-group hit counts
fn classify_seq(
    hits: &mut [u64; 128],
    num_groups: usize,
    total_kmers: usize,
    abs_threshold: u64,
    rel_threshold: f64,
) -> Classification {
    if total_kmers == 0 {
        return Classification::Unclassified;
    }

    let mut matching_mask = 0u128;
    let mut match_count = 0u32;
    let mut single_match = 0usize;

    for (group_idx, &group_hits) in hits[..num_groups].iter().enumerate() {
        if group_hits >= abs_threshold && (group_hits as f64 / total_kmers as f64) >= rel_threshold
        {
            matching_mask |= 1u128 << group_idx;
            match_count += 1;
            single_match = group_idx;
        }
    }

    match match_count {
        0 => Classification::Unclassified,
        1 => Classification::Classified(single_match),
        _ => Classification::Ambiguous(matching_mask),
    }
}

/// Update the classification spinner with current throughput stats
fn update_classify_spinner(
    spinner: &Option<Arc<Mutex<ProgressBar>>>,
    global_stats: &Mutex<ProcessingStats>,
    start_time: Instant,
) {
    if let Some(spinner) = spinner {
        let stats = global_stats.lock();
        let elapsed = start_time.elapsed();
        let seqs_per_sec = stats.total_seqs as f64 / elapsed.as_secs_f64();
        let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();

        spinner.lock().set_message(format!(
            "Classifying: {} seqs ({}). {:.0} seqs/s ({})",
            stats.total_seqs,
            format_bp(stats.total_bp as usize),
            seqs_per_sec,
            format_bp_per_sec(bp_per_sec)
        ));
    }
}

/// Classify one sequence and populate per-group hit counts
pub(crate) fn classify_seq_kmers(
    seq: &[u8],
    hasher: &KmerHasher,
    kmer_length: u8,
    smer_length: u8,
    buffers: &mut Buffers,
    hits: &mut [u64; 128],
    num_groups: usize,
    index: &ClassificationIndex,
    abs_threshold: u64,
    rel_threshold: f64,
) -> (usize, Classification) {
    fill_syncmers(seq, hasher, kmer_length, smer_length, buffers);

    for h in hits[..num_groups].iter_mut() {
        *h = 0;
    }

    let total_kmers = buffers.syncmers.len();

    match (&buffers.syncmers, index) {
        (SyncmerVec::U64(vec), ClassificationIndex::U64(map)) => {
            for &kmer in vec {
                if let Some(&bitmask) = map.get(&kmer) {
                    for group_idx in set_bits(bitmask) {
                        hits[group_idx] += 1;
                    }
                }
            }
        }
        (SyncmerVec::U128(vec), ClassificationIndex::U128(map)) => {
            for &kmer in vec {
                if let Some(&bitmask) = map.get(&kmer) {
                    for group_idx in set_bits(bitmask) {
                        hits[group_idx] += 1;
                    }
                }
            }
        }
        _ => panic!("Mismatch between SyncmerVec and ClassificationIndex types"),
    }

    let classification = classify_seq(hits, num_groups, total_kmers, abs_threshold, rel_threshold);

    (total_kmers, classification)
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
    fn add(&mut self, classification: Classification, seq_len: u64) {
        match classification {
            Classification::Classified(group_idx) => {
                self.groups[group_idx].seqs += 1;
                self.groups[group_idx].bases += seq_len;
            }
            Classification::Ambiguous(_) => {
                self.ambiguous_seqs += 1;
                self.ambiguous_bases += seq_len;
            }
            Classification::Unclassified => {
                self.unclassified_seqs += 1;
                self.unclassified_bases += seq_len;
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
        writer: Arc<Mutex<BufWriter<Box<dyn Write + Send>>>>,
    },
}

impl ClassifyOutput {
    fn add_record(
        &mut self,
        seq_id: &[u8],
        seq_len: usize,
        total_kmers: usize,
        classification: Classification,
        hits: &[u64; 128],
    ) {
        let (group_names, sample_name, local) = match self {
            Self::Summary { local, .. } => {
                local.add(classification, seq_len as u64);
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

    fn flush(&mut self) -> paraseq::parallel::Result<()> {
        match self {
            Self::Summary { local, global } => global.lock().merge(local),
            Self::PerSeq { local, writer, .. } if !local.is_empty() => {
                writer
                    .lock()
                    .write_all(local)
                    .map_err(paraseq::parallel::ProcessError::IoError)?;
                local.clear();
            }
            Self::PerSeq { .. } => {}
        }
        Ok(())
    }
}

#[derive(Clone)]
struct ClassifyProcessor {
    kmer_length: u8,
    smer_length: u8,
    hasher: KmerHasher,
    index: Arc<ClassificationIndex>,
    num_groups: usize,
    abs_threshold: u64,
    rel_threshold: f64,
    buffers: Buffers,
    hits: [u64; 128],
    output: ClassifyOutput,
    local_stats: ProcessingStats,
    global_stats: Arc<Mutex<ProcessingStats>>,
    spinner: Option<Arc<Mutex<ProgressBar>>>,
    start_time: Instant,
    limit_bp: Option<u64>,
}

impl ClassifyProcessor {
    #[allow(clippy::too_many_arguments)]
    fn new(
        kmer_length: u8,
        smer_length: u8,
        index: Arc<ClassificationIndex>,
        num_groups: usize,
        abs_threshold: u64,
        rel_threshold: f64,
        output: ClassifyOutput,
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

        Self {
            kmer_length,
            smer_length,
            hasher: KmerHasher::new(smer_length as usize),
            index,
            num_groups,
            abs_threshold,
            rel_threshold,
            buffers,
            hits: [0u64; 128],
            output,
            local_stats: ProcessingStats::default(),
            global_stats,
            spinner,
            start_time,
            limit_bp,
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for ClassifyProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        if let Some(limit) = self.limit_bp {
            let global_bp = self.global_stats.lock().total_bp;
            if global_bp >= limit {
                ParallelProcessor::<Rf>::on_batch_complete(self)?;
                return Err(paraseq::parallel::ProcessError::IoError(
                    sample_limit_reached_io_error(),
                ));
            }
        }

        let seq = record.seq();
        let seq_len = seq.len();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq_len as u64;

        let (total_kmers, classification) = classify_seq_kmers(
            &seq,
            &self.hasher,
            self.kmer_length,
            self.smer_length,
            &mut self.buffers,
            &mut self.hits,
            self.num_groups,
            &self.index,
            self.abs_threshold,
            self.rel_threshold,
        );

        self.output.add_record(
            record.id(),
            seq_len,
            total_kmers,
            classification,
            &self.hits,
        );

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        self.output.flush()?;
        let update_progress = {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;
            let current_progress = stats.total_bp / 100_000_000;
            if current_progress > stats.last_reported {
                stats.last_reported = current_progress;
                true
            } else {
                false
            }
        };
        if update_progress {
            update_classify_spinner(&self.spinner, &self.global_stats, self.start_time);
        }
        self.local_stats = ProcessingStats::default();

        Ok(())
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
        let mut options = String::new();
        if config.complexity > 0.0 {
            options.push_str(&format!(", complexity={}", config.complexity));
        }
        if let Some(limit) = config.limit_bp {
            options.push_str(&format!(", limit_bp={}", limit));
        }
        eprintln!(
            "Skope v{}; mode: classify (from {}); options: k={}, s={}, threads={}, abs_threshold={}, rel_threshold={:.2}{}",
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
            "Skope v{}; mode: classify (from index); options: threads={}, abs_threshold={}, rel_threshold={:.2}{}",
            version, config.threads, config.abs_threshold, config.rel_threshold, limit_str
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
                "Discriminatory mode: removed {} shared syncmers, {} unique syncmers remain",
                removed,
                index.len()
            );
        }
    }

    let index = Arc::new(index);
    let group_names = Arc::new(group_names);
    let num_groups = group_names.len();

    use rayon::prelude::*;
    if config.per_seq {
        let writer: Box<dyn Write + Send> = if let Some(path) = &config.output_path {
            Box::new(BufWriter::new(File::create(path)?))
        } else {
            Box::new(BufWriter::new(io::stdout()))
        };
        let writer = Arc::new(Mutex::new(BufWriter::new(writer)));

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
                sample_name,
                kmer_length,
                smer_length,
                &index,
                num_groups,
                config.abs_threshold,
                config.rel_threshold,
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
                    sample_name,
                    kmer_length,
                    smer_length,
                    &index,
                    num_groups,
                    config.abs_threshold,
                    config.rel_threshold,
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

        let writer: Box<dyn Write> = if let Some(path) = &config.output_path {
            Box::new(BufWriter::new(File::create(path)?))
        } else {
            Box::new(BufWriter::new(io::stdout()))
        };
        let mut writer = writer;

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
    sample_name: &str,
    kmer_length: u8,
    smer_length: u8,
    index: &Arc<ClassificationIndex>,
    num_groups: usize,
    abs_threshold: u64,
    rel_threshold: f64,
    threads: usize,
    quiet: bool,
    limit_bp: Option<u64>,
    mut make_output: impl FnMut() -> ClassifyOutput,
) -> Result<ProcessingStats> {
    let mut totals = ProcessingStats::default();
    for seq_path in sample_paths {
        if limit_bp.is_some_and(|limit| totals.total_bp >= limit) {
            break;
        }
        let in_path = (seq_path.to_string_lossy() != "-").then_some(seq_path.as_path());
        let spinner = create_spinner(quiet)?;
        let file_start = Instant::now();
        let global_stats = Arc::new(Mutex::new(ProcessingStats::default()));
        let mut processor = ClassifyProcessor::new(
            kmer_length,
            smer_length,
            Arc::clone(index),
            num_groups,
            abs_threshold,
            rel_threshold,
            make_output(),
            Arc::clone(&global_stats),
            spinner.clone(),
            file_start,
            limit_bp.map(|limit| limit.saturating_sub(totals.total_bp)),
        );

        let reader = reader_with_inferred_batch_size(in_path)?;
        handle_process_result(reader.process_parallel(&mut processor, threads))?;
        if let Some(ref pb) = spinner {
            pb.lock().finish_and_clear();
        }

        let stats = global_stats.lock().clone();
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

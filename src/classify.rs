use crate::syncmers::{Buffers, KmerHasher, SyncmerVec, fill_syncmers};
use crate::{
    FixedRapidHasher, ProcessingStats, create_spinner, format_bp, format_bp_per_sec,
    handle_process_result, reader_with_inferred_batch_size, sample_limit_reached_io_error,
};
use anyhow::{Context, Result};
use indicatif::ProgressBar;
use paraseq::Record;
use paraseq::parallel::{ParallelProcessor, ParallelReader};
use parking_lot::Mutex;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

// Index header constants and metadata
const INDEX_MAGIC: &[u8; 4] = b"GRAT";
const INDEX_FORMAT_VERSION: u8 = 1;
type ClassificationIndexHeader = ([u8; 4], u8, u8, u8, u8);

/// Classification index mapping k-mers to group bitmasks (up to 64 groups)
#[derive(Clone)]
pub enum ClassificationIndex {
    U64(HashMap<u64, u64, FixedRapidHasher>),
    U128(HashMap<u128, u64, FixedRapidHasher>),
}

impl ClassificationIndex {
    fn len(&self) -> usize {
        match self {
            ClassificationIndex::U64(m) => m.len(),
            ClassificationIndex::U128(m) => m.len(),
        }
    }
}

/// Remove k-mers shared across groups, keeping only group-unique k-mers
/// Returns how many shared k-mers were removed
fn apply_discriminatory_filter(index: &mut ClassificationIndex) -> usize {
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
enum Classification {
    Unclassified,
    Classified(usize),
    Ambiguous(u64),
}

#[derive(Debug, Clone, Default)]
struct GroupCounts {
    seqs: u64,
    bases: u64,
}

#[derive(Debug, Clone, Default)]
struct SampleClassificationResult {
    group_counts: Vec<GroupCounts>,
    ambiguous_seqs: u64,
    ambiguous_bases: u64,
    unclassified_seqs: u64,
    unclassified_bases: u64,
    total_seqs: u64,
    total_bases: u64,
}

// Configuration structs
pub struct BuildConfig {
    pub groups_dir: PathBuf,
    pub kmer_length: u8,
    pub smer_length: u8,
    pub threads: usize,
    pub output_path: Option<PathBuf>,
    pub quiet: bool,
}

pub struct ClassifyConfig {
    pub index_path: PathBuf,
    pub sample_paths: Vec<Vec<PathBuf>>,
    pub sample_names: Vec<String>,
    pub kmer_length: u8,
    pub smer_length: u8,
    pub min_hits: u64,
    pub min_fraction: f64,
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
    kmer_length: u8,
    smer_length: u8,
    hasher: KmerHasher,
    buffers: Buffers,
    group_bit: u64,

    // Thread-local k-mer map
    local_map_u64: Option<HashMap<u64, u64, FixedRapidHasher>>,
    local_map_u128: Option<HashMap<u128, u64, FixedRapidHasher>>,

    // Shared global state
    global_map_u64: Arc<Mutex<Option<HashMap<u64, u64, FixedRapidHasher>>>>,
    global_map_u128: Arc<Mutex<Option<HashMap<u128, u64, FixedRapidHasher>>>>,
    local_stats: ProcessingStats,
    global_stats: Arc<Mutex<ProcessingStats>>,
}

impl GroupKmerProcessor {
    fn new(
        kmer_length: u8,
        smer_length: u8,
        group_bit: u64,
        global_map_u64: Arc<Mutex<Option<HashMap<u64, u64, FixedRapidHasher>>>>,
        global_map_u128: Arc<Mutex<Option<HashMap<u128, u64, FixedRapidHasher>>>>,
        global_stats: Arc<Mutex<ProcessingStats>>,
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
            buffers,
            group_bit,
            local_map_u64,
            local_map_u128,
            global_map_u64,
            global_map_u128,
            local_stats: ProcessingStats::default(),
            global_stats,
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for GroupKmerProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
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

        match &self.buffers.syncmers {
            SyncmerVec::U64(vec) => {
                let local = self.local_map_u64.as_mut().unwrap();
                for &kmer in vec {
                    *local.entry(kmer).or_insert(0) |= self.group_bit;
                }
            }
            SyncmerVec::U128(vec) => {
                let local = self.local_map_u128.as_mut().unwrap();
                for &kmer in vec {
                    *local.entry(kmer).or_insert(0) |= self.group_bit;
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

use crate::{derive_sample_name, find_fastx_files};

/// Build a classification index in memory from group FASTA files in a directory
fn build_index_in_memory(
    groups_dir: &Path,
    kmer_length: u8,
    smer_length: u8,
    threads: usize,
    quiet: bool,
) -> Result<(ClassificationIndex, Vec<String>)> {
    let group_files = find_fastx_files(groups_dir)?;

    if group_files.len() > 64 {
        return Err(anyhow::anyhow!(
            "Too many groups: {} (max 64). Each FASTA file in the directory is one group.",
            group_files.len()
        ));
    }

    let group_names: Vec<String> = group_files
        .iter()
        .map(|p| derive_sample_name(p, false))
        .collect();

    {
        let mut seen = std::collections::HashSet::new();
        for name in &group_names {
            if !seen.insert(name) {
                return Err(anyhow::anyhow!(
                    "Duplicate group name derived from filenames: '{}'. Rename the FASTA files to have unique names.",
                    name
                ));
            }
        }
    }

    if !quiet {
        eprintln!(
            "Groups: {} files in {}",
            group_files.len(),
            groups_dir.display()
        );
    }

    let global_map_u64: Arc<Mutex<Option<HashMap<u64, u64, FixedRapidHasher>>>> =
        if kmer_length <= 32 {
            Arc::new(Mutex::new(Some(HashMap::with_hasher(FixedRapidHasher))))
        } else {
            Arc::new(Mutex::new(None))
        };

    let global_map_u128: Arc<Mutex<Option<HashMap<u128, u64, FixedRapidHasher>>>> =
        if kmer_length > 32 {
            Arc::new(Mutex::new(Some(HashMap::with_hasher(FixedRapidHasher))))
        } else {
            Arc::new(Mutex::new(None))
        };

    for (group_idx, (group_file, group_name)) in group_files.iter().zip(&group_names).enumerate() {
        let group_bit = 1u64 << group_idx;
        let global_stats = Arc::new(Mutex::new(ProcessingStats::default()));

        let mut processor = GroupKmerProcessor::new(
            kmer_length,
            smer_length,
            group_bit,
            Arc::clone(&global_map_u64),
            Arc::clone(&global_map_u128),
            Arc::clone(&global_stats),
        );

        let reader = reader_with_inferred_batch_size(Some(group_file))?;
        reader.process_parallel(&mut processor, threads)?;

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
                "  [{}] {}: {} seqs ({}), {} syncmers ({} unique)",
                group_idx,
                group_name,
                stats.total_seqs,
                format_bp(stats.total_bp as usize),
                group_kmers,
                unique_kmers,
            );
        }
    }

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

    Ok((index, group_names))
}

pub fn build_classification_index(config: &BuildConfig) -> Result<()> {
    let start_time = Instant::now();
    let version = env!("CARGO_PKG_VERSION");

    eprintln!(
        "Skope v{}; mode: index build; options: k={}, s={}, threads={}",
        version, config.kmer_length, config.smer_length, config.threads
    );

    let (index, group_names) = build_index_in_memory(
        &config.groups_dir,
        config.kmer_length,
        config.smer_length,
        config.threads,
        config.quiet,
    )?;

    save_index(
        &index,
        &group_names,
        config.kmer_length,
        config.smer_length,
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
        INDEX_FORMAT_VERSION,
        kmer_length,
        smer_length,
        group_names.len() as u8,
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

pub fn load_index(path: &Path) -> Result<(ClassificationIndex, Vec<String>, u8, u8)> {
    let file_bytes =
        fs::read(path).with_context(|| format!("Failed to open index file: {}", path.display()))?;
    let mut cursor = wincode::io::Cursor::new(file_bytes.as_slice());

    let header: ClassificationIndexHeader =
        wincode::deserialize_from(&mut cursor).context("Failed to decode index header")?;
    let (magic, format_version, kmer_length, smer_length, num_groups) = header;

    if &magic != INDEX_MAGIC {
        return Err(anyhow::anyhow!(
            "Not a skope classification index (invalid magic bytes)"
        ));
    }

    if format_version != INDEX_FORMAT_VERSION {
        return Err(anyhow::anyhow!(
            "Unsupported index format version: {} (expected {})",
            format_version,
            INDEX_FORMAT_VERSION
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
    let entry_size = kmer_bytes + 8; // k-mer bytes + group bitmask

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
        let mut map: HashMap<u64, u64, FixedRapidHasher> =
            HashMap::with_capacity_and_hasher(count, FixedRapidHasher);
        for i in 0..count {
            let offset = i * entry_size;
            let mut kmer_buf = [0u8; 8];
            kmer_buf[..kmer_bytes].copy_from_slice(&raw_data[offset..offset + kmer_bytes]);
            let kmer = u64::from_le_bytes(kmer_buf);

            let bitmask_offset = offset + kmer_bytes;
            let bitmask = u64::from_le_bytes(
                raw_data[bitmask_offset..bitmask_offset + 8]
                    .try_into()
                    .unwrap(),
            );

            map.insert(kmer, bitmask);
        }
        ClassificationIndex::U64(map)
    } else {
        let mut map: HashMap<u128, u64, FixedRapidHasher> =
            HashMap::with_capacity_and_hasher(count, FixedRapidHasher);
        for i in 0..count {
            let offset = i * entry_size;
            let mut kmer_buf = [0u8; 16];
            kmer_buf[..kmer_bytes].copy_from_slice(&raw_data[offset..offset + kmer_bytes]);
            let kmer = u128::from_le_bytes(kmer_buf);

            let bitmask_offset = offset + kmer_bytes;
            let bitmask = u64::from_le_bytes(
                raw_data[bitmask_offset..bitmask_offset + 8]
                    .try_into()
                    .unwrap(),
            );

            map.insert(kmer, bitmask);
        }
        ClassificationIndex::U128(map)
    };

    Ok((index, group_names, kmer_length, smer_length))
}

/// Classify one sequecne using per-group hit counts
fn classify_seq(
    hits: &mut [u64; 64],
    num_groups: usize,
    total_kmers: usize,
    min_hits: u64,
    min_fraction: f64,
) -> Classification {
    if total_kmers == 0 {
        return Classification::Unclassified;
    }

    let mut matching_mask = 0u64;
    let mut match_count = 0u32;
    let mut single_match = 0usize;

    for (group_idx, &group_hits) in hits[..num_groups].iter().enumerate() {
        if group_hits >= min_hits && (group_hits as f64 / total_kmers as f64) >= min_fraction {
            matching_mask |= 1u64 << group_idx;
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
fn classify_seq_kmers(
    seq: &[u8],
    hasher: &KmerHasher,
    kmer_length: u8,
    smer_length: u8,
    buffers: &mut Buffers,
    hits: &mut [u64; 64],
    num_groups: usize,
    index: &ClassificationIndex,
    min_hits: u64,
    min_fraction: f64,
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
                    let mut bits = bitmask;
                    while bits != 0 {
                        let group_idx = bits.trailing_zeros() as usize;
                        hits[group_idx] += 1;
                        bits &= bits - 1;
                    }
                }
            }
        }
        (SyncmerVec::U128(vec), ClassificationIndex::U128(map)) => {
            for &kmer in vec {
                if let Some(&bitmask) = map.get(&kmer) {
                    let mut bits = bitmask;
                    while bits != 0 {
                        let group_idx = bits.trailing_zeros() as usize;
                        hits[group_idx] += 1;
                        bits &= bits - 1;
                    }
                }
            }
        }
        _ => panic!("Mismatch between SyncmerVec and ClassificationIndex types"),
    }

    let classification = classify_seq(hits, num_groups, total_kmers, min_hits, min_fraction);

    (total_kmers, classification)
}

/// Shared state for summary classification mode
#[derive(Clone)]
struct GlobalClassifyState {
    group_seqs: Vec<u64>,
    group_bases: Vec<u64>,
    ambiguous_seqs: u64,
    ambiguous_bases: u64,
    unclassified_seqs: u64,
    unclassified_bases: u64,
    stats: ProcessingStats,
}

impl GlobalClassifyState {
    fn new(num_groups: usize) -> Self {
        Self {
            group_seqs: vec![0; num_groups],
            group_bases: vec![0; num_groups],
            ambiguous_seqs: 0,
            ambiguous_bases: 0,
            unclassified_seqs: 0,
            unclassified_bases: 0,
            stats: ProcessingStats::default(),
        }
    }
}

/// Processor for summary classification mode
#[derive(Clone)]
struct ClassifySummaryProcessor {
    kmer_length: u8,
    smer_length: u8,
    hasher: KmerHasher,
    index: Arc<ClassificationIndex>,
    num_groups: usize,
    min_hits: u64,
    min_fraction: f64,

    buffers: Buffers,
    hits: [u64; 64],

    local_group_seqs: Vec<u64>,
    local_group_bases: Vec<u64>,
    local_ambiguous_seqs: u64,
    local_ambiguous_bases: u64,
    local_unclassified_seqs: u64,
    local_unclassified_bases: u64,
    local_stats: ProcessingStats,

    global: Arc<Mutex<GlobalClassifyState>>,
    spinner: Option<Arc<Mutex<ProgressBar>>>,
    start_time: Instant,
    limit_bp: Option<u64>,
}

impl ClassifySummaryProcessor {
    fn new(
        kmer_length: u8,
        smer_length: u8,
        index: Arc<ClassificationIndex>,
        num_groups: usize,
        min_hits: u64,
        min_fraction: f64,
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
            min_hits,
            min_fraction,
            buffers,
            hits: [0u64; 64],
            local_group_seqs: vec![0; num_groups],
            local_group_bases: vec![0; num_groups],
            local_ambiguous_seqs: 0,
            local_ambiguous_bases: 0,
            local_unclassified_seqs: 0,
            local_unclassified_bases: 0,
            local_stats: ProcessingStats::default(),
            global: Arc::new(Mutex::new(GlobalClassifyState::new(num_groups))),
            spinner,
            start_time,
            limit_bp,
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for ClassifySummaryProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        if let Some(limit) = self.limit_bp {
            let global_bp = self.global.lock().stats.total_bp;
            if global_bp >= limit {
                return Err(paraseq::parallel::ProcessError::IoError(
                    sample_limit_reached_io_error(),
                ));
            }
        }

        let seq = record.seq();
        let seq_len = seq.len();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq_len as u64;

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
            Classification::Classified(group_idx) => {
                self.local_group_seqs[group_idx] += 1;
                self.local_group_bases[group_idx] += seq_len as u64;
            }
            Classification::Ambiguous(_) => {
                self.local_ambiguous_seqs += 1;
                self.local_ambiguous_bases += seq_len as u64;
            }
            Classification::Unclassified => {
                self.local_unclassified_seqs += 1;
                self.local_unclassified_bases += seq_len as u64;
            }
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        let mut g = self.global.lock();

        for i in 0..self.num_groups {
            g.group_seqs[i] += self.local_group_seqs[i];
            g.group_bases[i] += self.local_group_bases[i];
            self.local_group_seqs[i] = 0;
            self.local_group_bases[i] = 0;
        }

        g.ambiguous_seqs += self.local_ambiguous_seqs;
        g.ambiguous_bases += self.local_ambiguous_bases;
        g.unclassified_seqs += self.local_unclassified_seqs;
        g.unclassified_bases += self.local_unclassified_bases;
        self.local_ambiguous_seqs = 0;
        self.local_ambiguous_bases = 0;
        self.local_unclassified_seqs = 0;
        self.local_unclassified_bases = 0;

        g.stats.total_seqs += self.local_stats.total_seqs;
        g.stats.total_bp += self.local_stats.total_bp;

        let current_progress = g.stats.total_bp / 100_000_000;
        if current_progress > g.stats.last_reported {
            g.stats.last_reported = current_progress;
            drop(g);
            if let Some(ref spinner) = self.spinner {
                let g = self.global.lock();
                let elapsed = self.start_time.elapsed();
                let seqs_per_sec = g.stats.total_seqs as f64 / elapsed.as_secs_f64();
                let bp_per_sec = g.stats.total_bp as f64 / elapsed.as_secs_f64();
                spinner.lock().set_message(format!(
                    "Classifying: {} seqs ({}). {:.0} seqs/s ({})",
                    g.stats.total_seqs,
                    format_bp(g.stats.total_bp as usize),
                    seqs_per_sec,
                    format_bp_per_sec(bp_per_sec)
                ));
            }
        }

        self.local_stats = ProcessingStats::default();

        Ok(())
    }
}

/// Processor for per-sequence classification output mode
#[derive(Clone)]
struct ClassifyPerSeqProcessor {
    kmer_length: u8,
    smer_length: u8,
    hasher: KmerHasher,
    index: Arc<ClassificationIndex>,
    num_groups: usize,
    group_names: Arc<Vec<String>>,
    min_hits: u64,
    min_fraction: f64,
    sample_name: String,

    buffers: Buffers,
    hits: [u64; 64],

    local_output: Vec<u8>,
    output_writer: Arc<Mutex<BufWriter<Box<dyn Write + Send>>>>,

    local_stats: ProcessingStats,
    global_stats: Arc<Mutex<ProcessingStats>>,
    spinner: Option<Arc<Mutex<ProgressBar>>>,
    start_time: Instant,
    limit_bp: Option<u64>,
}

impl ClassifyPerSeqProcessor {
    fn new(
        kmer_length: u8,
        smer_length: u8,
        index: Arc<ClassificationIndex>,
        num_groups: usize,
        group_names: Arc<Vec<String>>,
        min_hits: u64,
        min_fraction: f64,
        sample_name: String,
        output_writer: Arc<Mutex<BufWriter<Box<dyn Write + Send>>>>,
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
            group_names,
            min_hits,
            min_fraction,
            sample_name,
            buffers,
            hits: [0u64; 64],
            local_output: Vec::with_capacity(64 * 1024),
            output_writer,
            local_stats: ProcessingStats::default(),
            global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
            spinner,
            start_time,
            limit_bp,
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for ClassifyPerSeqProcessor {
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
        let seq_id = String::from_utf8_lossy(record.id());
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
            self.min_hits,
            self.min_fraction,
        );

        use std::fmt::Write as FmtWrite;

        match classification {
            Classification::Classified(group_idx) => {
                let _ = writeln!(
                    self.local_output,
                    "{}\t{}\tclassified\t{}\t{}\t{}\t{}",
                    self.sample_name,
                    seq_id,
                    self.group_names[group_idx],
                    self.hits[group_idx],
                    total_kmers,
                    seq_len,
                );
            }
            Classification::Unclassified => {
                let _ = writeln!(
                    self.local_output,
                    "{}\t{}\tunclassified\t.\t0\t{}\t{}",
                    self.sample_name, seq_id, total_kmers, seq_len,
                );
            }
            Classification::Ambiguous(mask) => {
                let mut groups = String::new();
                let mut hits_str = String::new();
                let mut bits = mask;
                let mut first = true;
                while bits != 0 {
                    let group_idx = bits.trailing_zeros() as usize;
                    if !first {
                        groups.push(',');
                        hits_str.push(',');
                    }
                    groups.push_str(&self.group_names[group_idx]);
                    let _ = write!(hits_str, "{}", self.hits[group_idx]);
                    first = false;
                    bits &= bits - 1;
                }

                let _ = writeln!(
                    self.local_output,
                    "{}\t{}\tambiguous\t{}\t{}\t{}\t{}",
                    self.sample_name, seq_id, groups, hits_str, total_kmers, seq_len,
                );
            }
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        if !self.local_output.is_empty() {
            let mut writer = self.output_writer.lock();
            writer
                .write_all(&self.local_output)
                .map_err(paraseq::parallel::ProcessError::IoError)?;
            self.local_output.clear();
        }

        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;

            let current_progress = stats.total_bp / 100_000_000;
            if current_progress > stats.last_reported {
                drop(stats);
                update_classify_spinner(&self.spinner, &self.global_stats, self.start_time);
                self.global_stats.lock().last_reported = current_progress;
            }

            self.local_stats = ProcessingStats::default();
        }

        Ok(())
    }
}

pub fn run_classification(config: &ClassifyConfig) -> Result<()> {
    let start_time = Instant::now();
    let version = env!("CARGO_PKG_VERSION");

    let (index, group_names, kmer_length, smer_length) = if config.index_path.is_dir() {
        let limit_str = config
            .limit_bp
            .map_or(String::new(), |v| format!(", limit_bp={}", v));
        eprintln!(
            "Skope v{}; mode: classify (from directory); options: k={}, s={}, threads={}, min_hits={}, min_fraction={:.2}{}",
            version,
            config.kmer_length,
            config.smer_length,
            config.threads,
            config.min_hits,
            config.min_fraction,
            limit_str
        );

        let (index, group_names) = build_index_in_memory(
            &config.index_path,
            config.kmer_length,
            config.smer_length,
            config.threads,
            config.quiet,
        )?;

        (index, group_names, config.kmer_length, config.smer_length)
    } else {
        let limit_str = config
            .limit_bp
            .map_or(String::new(), |v| format!(", limit_bp={}", v));
        eprintln!(
            "Skope v{}; mode: classify (from index); options: threads={}, min_hits={}, min_fraction={:.2}{}",
            version, config.threads, config.min_hits, config.min_fraction, limit_str
        );

        let load_start = Instant::now();
        let (index, group_names, k, s) = load_index(&config.index_path)?;

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
            for seq_path in sample_paths {
                let in_path = if seq_path.to_string_lossy() == "-" {
                    None
                } else {
                    Some(seq_path.as_path())
                };

                let spinner = create_spinner(config.quiet)?;

                let pr_start = Instant::now();

                let mut processor = ClassifyPerSeqProcessor::new(
                    kmer_length,
                    smer_length,
                    Arc::clone(&index),
                    num_groups,
                    Arc::clone(&group_names),
                    config.min_hits,
                    config.min_fraction,
                    sample_name.clone(),
                    Arc::clone(&writer),
                    spinner.clone(),
                    pr_start,
                    config.limit_bp,
                );

                let reader = reader_with_inferred_batch_size(in_path)?;
                let result = reader.process_parallel(&mut processor, config.threads);
                handle_process_result(result)?;

                if let Some(ref pb) = spinner {
                    pb.lock().finish_and_clear();
                }

                let stats = processor.global_stats.lock().clone();
                if !config.quiet {
                    let elapsed = pr_start.elapsed();
                    let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();
                    eprintln!(
                        "Sample {}: {} seqs ({}) ({})",
                        sample_name,
                        stats.total_seqs,
                        format_bp(stats.total_bp as usize),
                        format_bp_per_sec(bp_per_sec)
                    );
                }
            }
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
                let result = process_sample_summary(
                    sample_paths,
                    sample_name,
                    kmer_length,
                    smer_length,
                    &index,
                    num_groups,
                    config.min_hits,
                    config.min_fraction,
                    config.threads,
                    config.quiet || is_multisample,
                    config.limit_bp,
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
            for (group_idx, counts) in result.group_counts.iter().enumerate() {
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
                    result.ambiguous_seqs as f64 / total_seqs * 100.0
                } else {
                    0.0
                };
                let pct_bases = if total_bases > 0.0 {
                    result.ambiguous_bases as f64 / total_bases * 100.0
                } else {
                    0.0
                };
                rows.push((
                    "ambiguous",
                    result.ambiguous_seqs,
                    pct_seqs,
                    result.ambiguous_bases,
                    pct_bases,
                ));
            }

            {
                let pct_seqs = if total_seqs > 0.0 {
                    result.unclassified_seqs as f64 / total_seqs * 100.0
                } else {
                    0.0
                };
                let pct_bases = if total_bases > 0.0 {
                    result.unclassified_bases as f64 / total_bases * 100.0
                } else {
                    0.0
                };
                rows.push((
                    "unclassified",
                    result.unclassified_seqs,
                    pct_seqs,
                    result.unclassified_bases,
                    pct_bases,
                ));
            }

            rows.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

            for (group_name, seqs, pct_seqs, bases, pct_bases) in &rows {
                writeln!(
                    writer,
                    "{}\t{}\t{:.2}\t{}\t{:.2}\t{}",
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

/// Process one sample in summary mode
fn process_sample_summary(
    sample_paths: &[PathBuf],
    sample_name: &str,
    kmer_length: u8,
    smer_length: u8,
    index: &Arc<ClassificationIndex>,
    num_groups: usize,
    min_hits: u64,
    min_fraction: f64,
    threads: usize,
    quiet: bool,
    limit_bp: Option<u64>,
) -> Result<SampleClassificationResult> {
    let mut combined = SampleClassificationResult {
        group_counts: vec![GroupCounts::default(); num_groups],
        ..Default::default()
    };

    for seq_path in sample_paths {
        let in_path = if seq_path.to_string_lossy() == "-" {
            None
        } else {
            Some(seq_path.as_path())
        };

        let spinner = create_spinner(quiet)?;

        let file_start = Instant::now();

        let mut processor = ClassifySummaryProcessor::new(
            kmer_length,
            smer_length,
            Arc::clone(index),
            num_groups,
            min_hits,
            min_fraction,
            spinner.clone(),
            file_start,
            limit_bp.map(|l| l.saturating_sub(combined.total_bases)),
        );

        let reader = reader_with_inferred_batch_size(in_path)?;
        let result = reader.process_parallel(&mut processor, threads);
        handle_process_result(result)?;

        if let Some(ref pb) = spinner {
            pb.lock().finish_and_clear();
        }

        let g = processor.global.lock();
        for i in 0..num_groups {
            combined.group_counts[i].seqs += g.group_seqs[i];
            combined.group_counts[i].bases += g.group_bases[i];
        }

        combined.ambiguous_seqs += g.ambiguous_seqs;
        combined.ambiguous_bases += g.ambiguous_bases;
        combined.unclassified_seqs += g.unclassified_seqs;
        combined.unclassified_bases += g.unclassified_bases;

        combined.total_seqs += g.stats.total_seqs;
        combined.total_bases += g.stats.total_bp;
        let stats = g.stats.clone();
        drop(g);

        if !quiet {
            let elapsed = file_start.elapsed();
            let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();
            eprintln!(
                "Sample {}: {} seqs ({}) ({})",
                sample_name,
                stats.total_seqs,
                format_bp(stats.total_bp as usize),
                format_bp_per_sec(bp_per_sec)
            );
        }

        if let Some(limit) = limit_bp
            && combined.total_bases >= limit
        {
            break;
        }
    }

    Ok(combined)
}

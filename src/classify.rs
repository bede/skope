use crate::containment::{FixedRapidHasher, format_bp, format_bp_per_sec};
use crate::minimizers::{Buffers, KmerHasher, MinimizerVec, fill_syncmers};
use anyhow::{Context, Result};
use bincode::{Decode, Encode};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use paraseq::Record;
use paraseq::fastx::Reader;
use paraseq::parallel::{ParallelProcessor, ParallelReader};
use parking_lot::Mutex;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read as IoRead, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

// ── Index header ──────────────────────────────────────────────────────────────

const INDEX_MAGIC: &[u8; 4] = b"GRAT";
const INDEX_FORMAT_VERSION: u8 = 1;

#[derive(Debug, Clone, Encode, Decode)]
struct ClassificationIndexHeader {
    magic: [u8; 4],
    format_version: u8,
    kmer_length: u8,
    smer_length: u8,
    num_groups: u8,
}

// ── Index ─────────────────────────────────────────────────────────────────────

/// Classification index: k-mer → group bitmask (u64, max 64 groups)
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

// ── Classification result types ───────────────────────────────────────────────

#[derive(Debug, Clone, Copy)]
enum Classification {
    Unclassified,
    Classified(usize),
    Ambiguous(u64),
}

#[derive(Debug, Clone, Default)]
struct GroupCounts {
    reads: u64,
    bases: u64,
    group_kmers_seen: u64,
    unique_kmers_seen: u64,
}

#[derive(Debug, Clone, Default)]
struct SampleClassificationResult {
    group_counts: Vec<GroupCounts>,
    ambiguous_reads: u64,
    ambiguous_bases: u64,
    unclassified_reads: u64,
    unclassified_bases: u64,
    total_reads: u64,
    total_bases: u64,
}

// ── Configs ───────────────────────────────────────────────────────────────────

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
    pub reads_paths: Vec<Vec<PathBuf>>,
    pub sample_names: Vec<String>,
    pub kmer_length: u8,
    pub smer_length: u8,
    pub min_hits: u64,
    pub min_fraction: f64,
    pub threads: usize,
    pub limit_bp: Option<u64>,
    pub output_path: Option<PathBuf>,
    pub per_read: bool,
    pub quiet: bool,
}

// ── Index building ────────────────────────────────────────────────────────────

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

/// Processor for collecting k-mers from a single group's FASTA files
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

    // Global state
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
            (
                Some(HashMap::with_hasher(FixedRapidHasher)),
                None,
            )
        } else {
            (
                None,
                Some(HashMap::with_hasher(FixedRapidHasher)),
            )
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

        match &self.buffers.minimizers {
            MinimizerVec::U64(vec) => {
                let local = self.local_map_u64.as_mut().unwrap();
                for &kmer in vec {
                    *local.entry(kmer).or_insert(0) |= self.group_bit;
                }
            }
            MinimizerVec::U128(vec) => {
                let local = self.local_map_u128.as_mut().unwrap();
                for &kmer in vec {
                    *local.entry(kmer).or_insert(0) |= self.group_bit;
                }
            }
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Merge local into global
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

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;
            self.local_stats = ProcessingStats::default();
        }

        Ok(())
    }
}

/// Find FASTA files in the groups directory (each file = one group)
fn find_group_files(groups_dir: &Path) -> Result<Vec<PathBuf>> {
    const FASTA_EXTENSIONS: &[&str] = &[
        ".fasta", ".fa", ".fastq", ".fq",
        ".fasta.gz", ".fa.gz", ".fastq.gz", ".fq.gz",
        ".fasta.xz", ".fa.xz", ".fastq.xz", ".fq.xz",
        ".fasta.zst", ".fa.zst", ".fastq.zst", ".fq.zst",
    ];

    let mut files = Vec::new();

    let entries = std::fs::read_dir(groups_dir)
        .with_context(|| format!("Failed to read directory: {}", groups_dir.display()))?;

    for entry in entries {
        let entry = entry
            .with_context(|| format!("Failed to read directory entry in {}", groups_dir.display()))?;

        let path = entry.path();

        // Skip hidden files
        if let Some(name) = path.file_name().and_then(|n| n.to_str())
            && name.starts_with('.')
        {
            continue;
        }

        // Follow symlinks
        let metadata = std::fs::metadata(&path)
            .with_context(|| format!("Failed to access: {}", path.display()))?;

        if !metadata.is_file() {
            continue;
        }

        let path_str = path.to_string_lossy().to_lowercase();
        if FASTA_EXTENSIONS.iter().any(|ext| path_str.ends_with(ext)) {
            files.push(path);
        }
    }

    if files.is_empty() {
        return Err(anyhow::anyhow!(
            "Directory contains no FASTA files: {}",
            groups_dir.display()
        ));
    }

    files.sort();
    Ok(files)
}

/// Derive group name from FASTA filename (strip extensions)
fn derive_group_name(path: &Path) -> String {
    let filename = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("unknown");

    let mut name = filename.to_string();
    let extensions = [".xz", ".gz", ".zst", ".fasta", ".fa"];

    loop {
        let original_len = name.len();
        for ext in &extensions {
            if name.ends_with(ext) {
                name = name[..name.len() - ext.len()].to_string();
                break;
            }
        }
        if name.len() == original_len {
            break;
        }
    }

    if name.is_empty() { filename.to_string() } else { name }
}

pub fn build_classification_index(config: &BuildConfig) -> Result<()> {
    let start_time = Instant::now();
    let version = env!("CARGO_PKG_VERSION");

    if !config.quiet {
        eprintln!(
            "Grate v{}; mode: index build; options: k={}, s={}, threads={}",
            version, config.kmer_length, config.smer_length, config.threads
        );
    }

    // Find group files
    let group_files = find_group_files(&config.groups_dir)?;

    if group_files.len() > 64 {
        return Err(anyhow::anyhow!(
            "Too many groups: {} (max 64). Each FASTA file in the directory is one group.",
            group_files.len()
        ));
    }

    let group_names: Vec<String> = group_files.iter().map(|p| derive_group_name(p)).collect();

    // Check for duplicate names
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

    if !config.quiet {
        eprintln!(
            "Groups: {} files in {}",
            group_files.len(),
            config.groups_dir.display()
        );
        for (i, name) in group_names.iter().enumerate() {
            eprintln!("  [{}] {}", i, name);
        }
    }

    // Build index: process each group's FASTA file(s)
    let global_map_u64: Arc<Mutex<Option<HashMap<u64, u64, FixedRapidHasher>>>> =
        if config.kmer_length <= 32 {
            Arc::new(Mutex::new(Some(HashMap::with_hasher(FixedRapidHasher))))
        } else {
            Arc::new(Mutex::new(None))
        };

    let global_map_u128: Arc<Mutex<Option<HashMap<u128, u64, FixedRapidHasher>>>> =
        if config.kmer_length > 32 {
            Arc::new(Mutex::new(Some(HashMap::with_hasher(FixedRapidHasher))))
        } else {
            Arc::new(Mutex::new(None))
        };

    for (group_idx, (group_file, group_name)) in
        group_files.iter().zip(&group_names).enumerate()
    {
        let group_bit = 1u64 << group_idx;
        let global_stats = Arc::new(Mutex::new(ProcessingStats::default()));

        let mut processor = GroupKmerProcessor::new(
            config.kmer_length,
            config.smer_length,
            group_bit,
            Arc::clone(&global_map_u64),
            Arc::clone(&global_map_u128),
            Arc::clone(&global_stats),
        );

        let reader = reader_with_inferred_batch_size(Some(group_file))?;
        reader.process_parallel(&mut processor, config.threads)?;

        let stats = global_stats.lock().clone();

        // Count k-mers belonging to this group
        let (group_kmers, unique_kmers) = if config.kmer_length <= 32 {
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

        if !config.quiet {
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

    // Build the index enum
    let index = if config.kmer_length <= 32 {
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

    if !config.quiet {
        let shared = match &index {
            ClassificationIndex::U64(m) => {
                m.values().filter(|v| v.count_ones() > 1).count()
            }
            ClassificationIndex::U128(m) => {
                m.values().filter(|v| v.count_ones() > 1).count()
            }
        };
        eprintln!(
            "Index: {} total syncmers, {} shared across groups",
            index.len(),
            shared
        );
    }

    // Save index
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

// ── Index serialization ───────────────────────────────────────────────────────

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

    let header = ClassificationIndexHeader {
        magic: *INDEX_MAGIC,
        format_version: INDEX_FORMAT_VERSION,
        kmer_length,
        smer_length,
        num_groups: group_names.len() as u8,
    };

    let bincode_config = bincode::config::standard().with_fixed_int_encoding();

    // Write header
    let header_bytes =
        bincode::encode_to_vec(&header, bincode_config)
            .context("Failed to encode index header")?;
    writer.write_all(&header_bytes)?;

    // Write group names
    let names_bytes = bincode::encode_to_vec(group_names, bincode_config)
        .context("Failed to encode group names")?;
    writer.write_all(&names_bytes)?;

    // Write entry count
    let count = index.len();
    let count_bytes = bincode::encode_to_vec(count, bincode_config)
        .context("Failed to encode entry count")?;
    writer.write_all(&count_bytes)?;

    // Write entries as raw bytes
    let kmer_bytes = (kmer_length as usize).div_ceil(4); // ceil(k/4)

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
    let file = File::open(path)
        .with_context(|| format!("Failed to open index file: {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let bincode_config = bincode::config::standard().with_fixed_int_encoding();

    // Read header
    let header: ClassificationIndexHeader =
        bincode::decode_from_std_read(&mut reader, bincode_config)
            .context("Failed to decode index header")?;

    if &header.magic != INDEX_MAGIC {
        return Err(anyhow::anyhow!(
            "Not a grate classification index (invalid magic bytes)"
        ));
    }

    if header.format_version != INDEX_FORMAT_VERSION {
        return Err(anyhow::anyhow!(
            "Unsupported index format version: {} (expected {})",
            header.format_version,
            INDEX_FORMAT_VERSION
        ));
    }

    // Read group names
    let group_names: Vec<String> =
        bincode::decode_from_std_read(&mut reader, bincode_config)
            .context("Failed to decode group names")?;

    if group_names.len() != header.num_groups as usize {
        return Err(anyhow::anyhow!(
            "Group count mismatch: header says {} but found {} names",
            header.num_groups,
            group_names.len()
        ));
    }

    // Read entry count
    let count: usize = bincode::decode_from_std_read(&mut reader, bincode_config)
        .context("Failed to decode entry count")?;

    // Read entries
    let kmer_bytes = (header.kmer_length as usize).div_ceil(4);
    let entry_size = kmer_bytes + 8; // kmer + bitmask

    // Read all remaining data
    let mut raw_data = Vec::new();
    reader.read_to_end(&mut raw_data)?;

    let expected_size = count * entry_size;
    if raw_data.len() < expected_size {
        return Err(anyhow::anyhow!(
            "Index file truncated: expected {} bytes of entries, got {}",
            expected_size,
            raw_data.len()
        ));
    }

    let index = if header.kmer_length <= 32 {
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

    Ok((index, group_names, header.kmer_length, header.smer_length))
}

// ── Classification pipeline ───────────────────────────────────────────────────

/// Classify a single read's k-mers against the index
fn classify_read(
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
        if group_hits >= min_hits
            && (group_hits as f64 / total_kmers as f64) >= min_fraction
        {
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

/// Processor for classifying reads (summary mode)
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

    // Thread-local accumulators
    local_group_reads: Vec<u64>,
    local_group_bases: Vec<u64>,
    local_group_kmer_hits: Vec<u64>,
    local_ambiguous_reads: u64,
    local_ambiguous_bases: u64,
    local_unclassified_reads: u64,
    local_unclassified_bases: u64,
    local_stats: ProcessingStats,

    // Global state
    global_group_reads: Arc<Mutex<Vec<u64>>>,
    global_group_bases: Arc<Mutex<Vec<u64>>>,
    global_group_kmer_hits: Arc<Mutex<Vec<u64>>>,
    global_ambiguous_reads: Arc<Mutex<u64>>,
    global_ambiguous_bases: Arc<Mutex<u64>>,
    global_unclassified_reads: Arc<Mutex<u64>>,
    global_unclassified_bases: Arc<Mutex<u64>>,
    global_stats: Arc<Mutex<ProcessingStats>>,
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
            local_group_reads: vec![0; num_groups],
            local_group_bases: vec![0; num_groups],
            local_group_kmer_hits: vec![0; num_groups],
            local_ambiguous_reads: 0,
            local_ambiguous_bases: 0,
            local_unclassified_reads: 0,
            local_unclassified_bases: 0,
            local_stats: ProcessingStats::default(),
            global_group_reads: Arc::new(Mutex::new(vec![0; num_groups])),
            global_group_bases: Arc::new(Mutex::new(vec![0; num_groups])),
            global_group_kmer_hits: Arc::new(Mutex::new(vec![0; num_groups])),
            global_ambiguous_reads: Arc::new(Mutex::new(0)),
            global_ambiguous_bases: Arc::new(Mutex::new(0)),
            global_unclassified_reads: Arc::new(Mutex::new(0)),
            global_unclassified_bases: Arc::new(Mutex::new(0)),
            global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
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
                "Classifying: {} reads ({}). {:.0} reads/s ({})",
                stats.total_seqs,
                format_bp(stats.total_bp as usize),
                reads_per_sec,
                format_bp_per_sec(bp_per_sec)
            ));
        }
    }

    fn process_read_kmers(&mut self, seq: &[u8]) -> (usize, Classification) {
        fill_syncmers(
            seq,
            &self.hasher,
            self.kmer_length,
            self.smer_length,
            &mut self.buffers,
        );

        // Reset hits
        for h in self.hits[..self.num_groups].iter_mut() {
            *h = 0;
        }

        let total_kmers = self.buffers.minimizers.len();

        match (&self.buffers.minimizers, &*self.index) {
            (MinimizerVec::U64(vec), ClassificationIndex::U64(map)) => {
                for &kmer in vec {
                    if let Some(&bitmask) = map.get(&kmer) {
                        // Increment hit count for each group in the bitmask
                        let mut bits = bitmask;
                        while bits != 0 {
                            let group_idx = bits.trailing_zeros() as usize;
                            self.hits[group_idx] += 1;
                            bits &= bits - 1; // Clear lowest set bit
                        }
                    }
                }
            }
            (MinimizerVec::U128(vec), ClassificationIndex::U128(map)) => {
                for &kmer in vec {
                    if let Some(&bitmask) = map.get(&kmer) {
                        let mut bits = bitmask;
                        while bits != 0 {
                            let group_idx = bits.trailing_zeros() as usize;
                            self.hits[group_idx] += 1;
                            bits &= bits - 1;
                        }
                    }
                }
            }
            _ => panic!("Mismatch between MinimizerVec and ClassificationIndex types"),
        }

        let classification = classify_read(
            &mut self.hits,
            self.num_groups,
            total_kmers,
            self.min_hits,
            self.min_fraction,
        );

        (total_kmers, classification)
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for ClassifySummaryProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        // Check limit
        if let Some(limit) = self.limit_bp {
            let global_bp = self.global_stats.lock().total_bp;
            if global_bp >= limit {
                return Err(paraseq::parallel::ProcessError::IoError(
                    std::io::Error::new(std::io::ErrorKind::Interrupted, "Sample limit reached"),
                ));
            }
        }

        let seq = record.seq();
        let read_len = seq.len();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += read_len as u64;

        let (_total_kmers, classification) = self.process_read_kmers(&seq);

        match classification {
            Classification::Classified(group_idx) => {
                self.local_group_reads[group_idx] += 1;
                self.local_group_bases[group_idx] += read_len as u64;
                self.local_group_kmer_hits[group_idx] += self.hits[group_idx];
            }
            Classification::Ambiguous(_) => {
                self.local_ambiguous_reads += 1;
                self.local_ambiguous_bases += read_len as u64;
            }
            Classification::Unclassified => {
                self.local_unclassified_reads += 1;
                self.local_unclassified_bases += read_len as u64;
            }
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Merge group counts
        {
            let mut gr = self.global_group_reads.lock();
            let mut gb = self.global_group_bases.lock();
            let mut gk = self.global_group_kmer_hits.lock();
            for i in 0..self.num_groups {
                gr[i] += self.local_group_reads[i];
                gb[i] += self.local_group_bases[i];
                gk[i] += self.local_group_kmer_hits[i];
                self.local_group_reads[i] = 0;
                self.local_group_bases[i] = 0;
                self.local_group_kmer_hits[i] = 0;
            }
        }

        *self.global_ambiguous_reads.lock() += self.local_ambiguous_reads;
        *self.global_ambiguous_bases.lock() += self.local_ambiguous_bases;
        *self.global_unclassified_reads.lock() += self.local_unclassified_reads;
        *self.global_unclassified_bases.lock() += self.local_unclassified_bases;
        self.local_ambiguous_reads = 0;
        self.local_ambiguous_bases = 0;
        self.local_unclassified_reads = 0;
        self.local_unclassified_bases = 0;

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;

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

/// Processor for per-read classification output
#[derive(Clone)]
struct ClassifyPerReadProcessor {
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

    // Thread-local output buffer
    local_output: Vec<u8>,
    output_writer: Arc<Mutex<BufWriter<Box<dyn Write + Send>>>>,

    local_stats: ProcessingStats,
    global_stats: Arc<Mutex<ProcessingStats>>,
    spinner: Option<Arc<Mutex<ProgressBar>>>,
    start_time: Instant,
    limit_bp: Option<u64>,
}

impl ClassifyPerReadProcessor {
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

    fn update_spinner(&self) {
        if let Some(ref spinner) = self.spinner {
            let stats = self.global_stats.lock();
            let elapsed = self.start_time.elapsed();
            let reads_per_sec = stats.total_seqs as f64 / elapsed.as_secs_f64();
            let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();

            spinner.lock().set_message(format!(
                "Classifying: {} reads ({}). {:.0} reads/s ({})",
                stats.total_seqs,
                format_bp(stats.total_bp as usize),
                reads_per_sec,
                format_bp_per_sec(bp_per_sec)
            ));
        }
    }

    fn process_read_kmers(&mut self, seq: &[u8]) -> (usize, Classification) {
        fill_syncmers(
            seq,
            &self.hasher,
            self.kmer_length,
            self.smer_length,
            &mut self.buffers,
        );

        for h in self.hits[..self.num_groups].iter_mut() {
            *h = 0;
        }

        let total_kmers = self.buffers.minimizers.len();

        match (&self.buffers.minimizers, &*self.index) {
            (MinimizerVec::U64(vec), ClassificationIndex::U64(map)) => {
                for &kmer in vec {
                    if let Some(&bitmask) = map.get(&kmer) {
                        let mut bits = bitmask;
                        while bits != 0 {
                            let group_idx = bits.trailing_zeros() as usize;
                            self.hits[group_idx] += 1;
                            bits &= bits - 1;
                        }
                    }
                }
            }
            (MinimizerVec::U128(vec), ClassificationIndex::U128(map)) => {
                for &kmer in vec {
                    if let Some(&bitmask) = map.get(&kmer) {
                        let mut bits = bitmask;
                        while bits != 0 {
                            let group_idx = bits.trailing_zeros() as usize;
                            self.hits[group_idx] += 1;
                            bits &= bits - 1;
                        }
                    }
                }
            }
            _ => panic!("Mismatch between MinimizerVec and ClassificationIndex types"),
        }

        let classification = classify_read(
            &mut self.hits,
            self.num_groups,
            total_kmers,
            self.min_hits,
            self.min_fraction,
        );

        (total_kmers, classification)
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for ClassifyPerReadProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        if let Some(limit) = self.limit_bp {
            let global_bp = self.global_stats.lock().total_bp;
            if global_bp >= limit {
                return Err(paraseq::parallel::ProcessError::IoError(
                    std::io::Error::new(std::io::ErrorKind::Interrupted, "Sample limit reached"),
                ));
            }
        }

        let seq = record.seq();
        let read_id = String::from_utf8_lossy(record.id());
        let read_len = seq.len();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += read_len as u64;

        let (total_kmers, classification) = self.process_read_kmers(&seq);

        use std::fmt::Write as FmtWrite;

        match classification {
            Classification::Classified(group_idx) => {
                let _ = writeln!(
                    self.local_output,
                    "{}\t{}\tclassified\t{}\t{}\t{}\t{}",
                    self.sample_name,
                    read_id,
                    self.group_names[group_idx],
                    self.hits[group_idx],
                    total_kmers,
                    read_len,
                );
            }
            Classification::Unclassified => {
                let _ = writeln!(
                    self.local_output,
                    "{}\t{}\tunclassified\t.\t0\t{}\t{}",
                    self.sample_name, read_id, total_kmers, read_len,
                );
            }
            Classification::Ambiguous(mask) => {
                // Build comma-separated group names and hits
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
                    self.sample_name, read_id, groups, hits_str, total_kmers, read_len,
                );
            }
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Flush local output to global writer
        if !self.local_output.is_empty() {
            let mut writer = self.output_writer.lock();
            writer
                .write_all(&self.local_output)
                .map_err(paraseq::parallel::ProcessError::IoError)?;
            self.local_output.clear();
        }

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;

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

// ── Main classification entry point ───────────────────────────────────────────

pub fn run_classification(config: &ClassifyConfig) -> Result<()> {
    let start_time = Instant::now();
    let version = env!("CARGO_PKG_VERSION");

    // Load or build index
    let (index, group_names, kmer_length, smer_length) = if config.index_path.is_dir() {
        if !config.quiet {
            eprintln!(
                "Grate v{}; mode: classify (from directory); options: k={}, s={}, threads={}, min_hits={}, min_fraction={:.3}",
                version, config.kmer_length, config.smer_length, config.threads, config.min_hits, config.min_fraction
            );
        }

        // Build index on the fly
        let build_config = BuildConfig {
            groups_dir: config.index_path.clone(),
            kmer_length: config.kmer_length,
            smer_length: config.smer_length,
            threads: config.threads,
            output_path: None, // Don't save
            quiet: config.quiet,
        };

        // We need to build the index but not save it - inline the building logic
        let group_files = find_group_files(&config.index_path)?;
        if group_files.len() > 64 {
            return Err(anyhow::anyhow!(
                "Too many groups: {} (max 64)",
                group_files.len()
            ));
        }

        let group_names: Vec<String> = group_files.iter().map(|p| derive_group_name(p)).collect();
        {
            let mut seen = std::collections::HashSet::new();
            for name in &group_names {
                if !seen.insert(name) {
                    return Err(anyhow::anyhow!(
                        "Duplicate group name: '{}'",
                        name
                    ));
                }
            }
        }

        if !config.quiet {
            eprintln!(
                "Groups: {} files in {}",
                group_files.len(),
                config.index_path.display()
            );
        }

        let global_map_u64: Arc<Mutex<Option<HashMap<u64, u64, FixedRapidHasher>>>> =
            if build_config.kmer_length <= 32 {
                Arc::new(Mutex::new(Some(HashMap::with_hasher(FixedRapidHasher))))
            } else {
                Arc::new(Mutex::new(None))
            };

        let global_map_u128: Arc<Mutex<Option<HashMap<u128, u64, FixedRapidHasher>>>> =
            if build_config.kmer_length > 32 {
                Arc::new(Mutex::new(Some(HashMap::with_hasher(FixedRapidHasher))))
            } else {
                Arc::new(Mutex::new(None))
            };

        for (group_idx, (group_file, group_name)) in
            group_files.iter().zip(&group_names).enumerate()
        {
            let group_bit = 1u64 << group_idx;
            let global_stats = Arc::new(Mutex::new(ProcessingStats::default()));

            let mut processor = GroupKmerProcessor::new(
                build_config.kmer_length,
                build_config.smer_length,
                group_bit,
                Arc::clone(&global_map_u64),
                Arc::clone(&global_map_u128),
                Arc::clone(&global_stats),
            );

            let reader = reader_with_inferred_batch_size(Some(group_file))?;
            reader.process_parallel(&mut processor, build_config.threads)?;

            let stats = global_stats.lock();
            if !config.quiet {
                eprintln!(
                    "  [{}] {}: {} seqs ({})",
                    group_idx,
                    group_name,
                    stats.total_seqs,
                    format_bp(stats.total_bp as usize),
                );
            }
        }

        let index = if build_config.kmer_length <= 32 {
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

        if !config.quiet {
            eprintln!("Index: {} k-mers (built on-the-fly)", index.len());
        }

        (index, group_names, build_config.kmer_length, build_config.smer_length)
    } else {
        // Load pre-built index
        if !config.quiet {
            eprintln!(
                "Grate v{}; mode: classify (from index); options: threads={}, min_hits={}, min_fraction={:.3}",
                version, config.threads, config.min_hits, config.min_fraction
            );
        }

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

    let index = Arc::new(index);
    let group_names = Arc::new(group_names);
    let num_groups = group_names.len();

    // Process samples
    use rayon::prelude::*;

    if config.per_read {
        // Per-read output mode
        let writer: Box<dyn Write + Send> = if let Some(path) = &config.output_path {
            Box::new(BufWriter::new(File::create(path)?))
        } else {
            Box::new(BufWriter::new(io::stdout()))
        };
        let writer = Arc::new(Mutex::new(BufWriter::new(writer)));

        // Write header
        {
            let mut w = writer.lock();
            writeln!(
                w,
                "sample\tread_id\tclassification\tgroup\thits\tread_kmers\tread_length"
            )?;
        }

        for (sample_paths, sample_name) in
            config.reads_paths.iter().zip(&config.sample_names)
        {
            for reads_path in sample_paths {
                let in_path = if reads_path.to_string_lossy() == "-" {
                    None
                } else {
                    Some(reads_path.as_path())
                };

                let spinner = if !config.quiet {
                    let pb = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr());
                    pb.set_style(
                        ProgressStyle::default_spinner()
                            .tick_strings(&[
                                "⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏",
                            ])
                            .template("{msg}")?,
                    );
                    Some(Arc::new(Mutex::new(pb)))
                } else {
                    None
                };

                let pr_start = Instant::now();

                let mut processor = ClassifyPerReadProcessor::new(
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

                let stopped_early = if let Err(ref e) = result {
                    e.to_string().contains("Sample limit reached")
                } else {
                    false
                };
                if !stopped_early {
                    result?;
                }

                if let Some(ref pb) = spinner {
                    pb.lock().finish_and_clear();
                }

                let stats = processor.global_stats.lock().clone();
                if !config.quiet {
                    let elapsed = pr_start.elapsed();
                    let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();
                    eprintln!(
                        "Sample {}: {} reads ({}) ({})",
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
        // Summary mode
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

        let sample_results: Vec<(String, SampleClassificationResult)> = config
            .reads_paths
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
                        config.reads_paths.len()
                    );
                }

                result.map(|r| (sample_name.clone(), r))
            })
            .collect::<Result<Vec<_>>>()?;

        if is_multisample && !config.quiet {
            eprintln!();
        }

        // Output summary TSV
        let writer: Box<dyn Write> = if let Some(path) = &config.output_path {
            Box::new(BufWriter::new(File::create(path)?))
        } else {
            Box::new(BufWriter::new(io::stdout()))
        };
        let mut writer = writer;

        writeln!(
            writer,
            "sample\tgroup\treads\tbases\tpct_reads\tpct_bases\tgroup_kmers\tunique_kmers"
        )?;

        for (sample_name, result) in &sample_results {
            let total_reads = result.total_reads as f64;
            let total_bases = result.total_bases as f64;

            for (group_idx, counts) in result.group_counts.iter().enumerate() {
                let pct_reads = if total_reads > 0.0 {
                    counts.reads as f64 / total_reads * 100.0
                } else {
                    0.0
                };
                let pct_bases = if total_bases > 0.0 {
                    counts.bases as f64 / total_bases * 100.0
                } else {
                    0.0
                };

                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}\t{}",
                    sample_name,
                    group_names[group_idx],
                    counts.reads,
                    counts.bases,
                    pct_reads,
                    pct_bases,
                    counts.group_kmers_seen,
                    counts.unique_kmers_seen,
                )?;
            }

            // Ambiguous row
            let pct_reads = if total_reads > 0.0 {
                result.ambiguous_reads as f64 / total_reads * 100.0
            } else {
                0.0
            };
            let pct_bases = if total_bases > 0.0 {
                result.ambiguous_bases as f64 / total_bases * 100.0
            } else {
                0.0
            };
            writeln!(
                writer,
                "{}\tambiguous\t{}\t{}\t{:.2}\t{:.2}\t.\t.",
                sample_name, result.ambiguous_reads, result.ambiguous_bases, pct_reads, pct_bases,
            )?;

            // Unclassified row
            let pct_reads = if total_reads > 0.0 {
                result.unclassified_reads as f64 / total_reads * 100.0
            } else {
                0.0
            };
            let pct_bases = if total_bases > 0.0 {
                result.unclassified_bases as f64 / total_bases * 100.0
            } else {
                0.0
            };
            writeln!(
                writer,
                "{}\tunclassified\t{}\t{}\t{:.2}\t{:.2}\t.\t.",
                sample_name,
                result.unclassified_reads,
                result.unclassified_bases,
                pct_reads,
                pct_bases,
            )?;
        }

        writer.flush()?;
    }

    if !config.quiet {
        let elapsed = start_time.elapsed();
        eprintln!("Done in {:.1}s", elapsed.as_secs_f64());
    }

    Ok(())
}

/// Process a single sample in summary mode
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

    for reads_path in sample_paths {
        let in_path = if reads_path.to_string_lossy() == "-" {
            None
        } else {
            Some(reads_path.as_path())
        };

        let spinner = if !quiet {
            let pb = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr());
            pb.set_style(
                ProgressStyle::default_spinner()
                    .tick_strings(&[
                        "⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏",
                    ])
                    .template("{msg}")?,
            );
            Some(Arc::new(Mutex::new(pb)))
        } else {
            None
        };

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

        let stopped_early = if let Err(ref e) = result {
            e.to_string().contains("Sample limit reached")
        } else {
            false
        };
        if !stopped_early {
            result?;
        }

        if let Some(ref pb) = spinner {
            pb.lock().finish_and_clear();
        }

        // Merge results
        let gr = processor.global_group_reads.lock();
        let gb = processor.global_group_bases.lock();
        let gk = processor.global_group_kmer_hits.lock();
        for i in 0..num_groups {
            combined.group_counts[i].reads += gr[i];
            combined.group_counts[i].bases += gb[i];
            combined.group_counts[i].group_kmers_seen += gk[i];
        }

        combined.ambiguous_reads += *processor.global_ambiguous_reads.lock();
        combined.ambiguous_bases += *processor.global_ambiguous_bases.lock();
        combined.unclassified_reads += *processor.global_unclassified_reads.lock();
        combined.unclassified_bases += *processor.global_unclassified_bases.lock();

        let stats = processor.global_stats.lock().clone();
        combined.total_reads += stats.total_seqs;
        combined.total_bases += stats.total_bp;

        if !quiet {
            let elapsed = file_start.elapsed();
            let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();
            eprintln!(
                "Sample {}: {} reads ({}) ({})",
                sample_name,
                stats.total_seqs,
                format_bp(stats.total_bp as usize),
                format_bp_per_sec(bp_per_sec)
            );
        }

        // Check limit
        if let Some(limit) = limit_bp
            && combined.total_bases >= limit
        {
            break;
        }
    }

    Ok(combined)
}

use crate::{
    FixedRapidHasher, ProcessingStats, create_spinner, format_bp, format_bp_per_sec,
    handle_process_result, reader_with_inferred_batch_size,
};
use crate::minimizers::{Buffers, KmerHasher, MinimizerVec, fill_syncmers};
use anyhow::{Context, Result};
use bincode::{Decode, Encode};
use indicatif::ProgressBar;
use paraseq::Record;
use paraseq::parallel::{ParallelProcessor, ParallelReader};
use parking_lot::Mutex;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read as IoRead, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

// Index header constants and metadata.

const INDEX_MAGIC: &[u8; 4] = b"GRAT";
const INDEX_FORMAT_VERSION: u8 = 1;

#[derive(Debug, Clone, Encode, Decode)]
struct ClassificationIndexHeader {
    magic: [u8; 4],
    format_version: u8,
    kmer_length: u8,
    smer_length: u8,
    num_classes: u8,
}

/// Classification index mapping k-mers to class bitmasks (up to 64 classes).
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

/// Remove k-mers shared across classes, keeping only class-unique k-mers.
/// Returns how many shared k-mers were removed.
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

// Classification result types.

#[derive(Debug, Clone, Copy)]
enum Classification {
    Unclassified,
    Classified(usize),
    Ambiguous(u64),
}

#[derive(Debug, Clone, Default)]
struct ClassCounts {
    seqs: u64,
    bases: u64,
}

#[derive(Debug, Clone, Default)]
struct SampleClassificationResult {
    class_counts: Vec<ClassCounts>,
    ambiguous_seqs: u64,
    ambiguous_bases: u64,
    unclassified_seqs: u64,
    unclassified_bases: u64,
    total_seqs: u64,
    total_bases: u64,
}

// Configuration structs.

pub struct BuildConfig {
    pub classes_dir: PathBuf,
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

/// Collect k-mers from one class FASTA file.
#[derive(Clone)]
struct GroupKmerProcessor {
    kmer_length: u8,
    smer_length: u8,
    hasher: KmerHasher,
    buffers: Buffers,
    class_bit: u64,

    // Thread-local k-mer map.
    local_map_u64: Option<HashMap<u64, u64, FixedRapidHasher>>,
    local_map_u128: Option<HashMap<u128, u64, FixedRapidHasher>>,

    // Shared global state.
    global_map_u64: Arc<Mutex<Option<HashMap<u64, u64, FixedRapidHasher>>>>,
    global_map_u128: Arc<Mutex<Option<HashMap<u128, u64, FixedRapidHasher>>>>,
    local_stats: ProcessingStats,
    global_stats: Arc<Mutex<ProcessingStats>>,
}

impl GroupKmerProcessor {
    fn new(
        kmer_length: u8,
        smer_length: u8,
        class_bit: u64,
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
            class_bit,
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
                    *local.entry(kmer).or_insert(0) |= self.class_bit;
                }
            }
            MinimizerVec::U128(vec) => {
                let local = self.local_map_u128.as_mut().unwrap();
                for &kmer in vec {
                    *local.entry(kmer).or_insert(0) |= self.class_bit;
                }
            }
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Merge local maps into the shared global map.
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

        // Accumulate per-thread stats into shared stats.
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

/// Build a classification index in memory from class FASTA files in a directory.
fn build_index_in_memory(
    classes_dir: &Path,
    kmer_length: u8,
    smer_length: u8,
    threads: usize,
    quiet: bool,
) -> Result<(ClassificationIndex, Vec<String>)> {
    let class_files = find_fastx_files(classes_dir)?;

    if class_files.len() > 64 {
        return Err(anyhow::anyhow!(
            "Too many classes: {} (max 64). Each FASTA file in the directory is one class.",
            class_files.len()
        ));
    }

    let class_names: Vec<String> = class_files.iter().map(|p| derive_sample_name(p, false)).collect();

    // Check for duplicate class names.
    {
        let mut seen = std::collections::HashSet::new();
        for name in &class_names {
            if !seen.insert(name) {
                return Err(anyhow::anyhow!(
                    "Duplicate class name derived from filenames: '{}'. Rename the FASTA files to have unique names.",
                    name
                ));
            }
        }
    }

    if !quiet {
        eprintln!(
            "Classes: {} files in {}",
            class_files.len(),
            classes_dir.display()
        );
    }

    // Build the index by processing each class FASTA file.
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

    for (class_idx, (class_file, class_name)) in
        class_files.iter().zip(&class_names).enumerate()
    {
        let class_bit = 1u64 << class_idx;
        let global_stats = Arc::new(Mutex::new(ProcessingStats::default()));

        let mut processor = GroupKmerProcessor::new(
            kmer_length,
            smer_length,
            class_bit,
            Arc::clone(&global_map_u64),
            Arc::clone(&global_map_u128),
            Arc::clone(&global_stats),
        );

        let reader = reader_with_inferred_batch_size(Some(class_file))?;
        reader.process_parallel(&mut processor, threads)?;

        let stats = global_stats.lock().clone();

        // Count k-mers observed for this class.
        let (class_kmers, unique_kmers) = if kmer_length <= 32 {
            let map = global_map_u64.lock();
            let map = map.as_ref().unwrap();
            let class_kmers = map.values().filter(|&&v| v & class_bit != 0).count();
            let unique = map
                .values()
                .filter(|&&v| v & class_bit != 0 && v.count_ones() == 1)
                .count();
            (class_kmers, unique)
        } else {
            let map = global_map_u128.lock();
            let map = map.as_ref().unwrap();
            let class_kmers = map.values().filter(|&&v| v & class_bit != 0).count();
            let unique = map
                .values()
                .filter(|&&v| v & class_bit != 0 && v.count_ones() == 1)
                .count();
            (class_kmers, unique)
        };

        if !quiet {
            eprintln!(
                "  [{}] {}: {} seqs ({}), {} syncmers ({} unique)",
                class_idx,
                class_name,
                stats.total_seqs,
                format_bp(stats.total_bp as usize),
                class_kmers,
                unique_kmers,
            );
        }
    }

    // Materialize the final index variant.
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
            ClassificationIndex::U64(m) => {
                m.values().filter(|v| v.count_ones() > 1).count()
            }
            ClassificationIndex::U128(m) => {
                m.values().filter(|v| v.count_ones() > 1).count()
            }
        };
        eprintln!(
            "Index: {} total syncmers, {} shared across classes",
            index.len(),
            shared
        );
    }

    Ok((index, class_names))
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

    let (index, class_names) = build_index_in_memory(
        &config.classes_dir,
        config.kmer_length,
        config.smer_length,
        config.threads,
        config.quiet,
    )?;

    save_index(
        &index,
        &class_names,
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

// Index serialization.

fn save_index(
    index: &ClassificationIndex,
    class_names: &[String],
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
        num_classes: class_names.len() as u8,
    };

    let bincode_config = bincode::config::standard().with_fixed_int_encoding();

    // Write header.
    let header_bytes =
        bincode::encode_to_vec(&header, bincode_config)
            .context("Failed to encode index header")?;
    writer.write_all(&header_bytes)?;

    // Write class names.
    let names_bytes = bincode::encode_to_vec(class_names, bincode_config)
        .context("Failed to encode class names")?;
    writer.write_all(&names_bytes)?;

    // Write entry count.
    let count = index.len();
    let count_bytes = bincode::encode_to_vec(count, bincode_config)
        .context("Failed to encode entry count")?;
    writer.write_all(&count_bytes)?;

    // Write entries as raw bytes.
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
    let file = File::open(path)
        .with_context(|| format!("Failed to open index file: {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let bincode_config = bincode::config::standard().with_fixed_int_encoding();

    // Read header.
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

    // Read class names.
    let class_names: Vec<String> =
        bincode::decode_from_std_read(&mut reader, bincode_config)
            .context("Failed to decode class names")?;

    if class_names.len() != header.num_classes as usize {
        return Err(anyhow::anyhow!(
            "Class count mismatch: header says {} but found {} names",
            header.num_classes,
            class_names.len()
        ));
    }

    // Read entry count.
    let count: usize = bincode::decode_from_std_read(&mut reader, bincode_config)
        .context("Failed to decode entry count")?;

    // Read entries.
    let kmer_bytes = (header.kmer_length as usize).div_ceil(4);
    let entry_size = kmer_bytes + 8; // k-mer bytes + class bitmask

    // Read the remaining entry payload.
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

    Ok((index, class_names, header.kmer_length, header.smer_length))
}

/// Classify one sequecne using per-class hit counts.
fn classify_seq(
    hits: &mut [u64; 64],
    num_classes: usize,
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

    for (class_idx, &class_hits) in hits[..num_classes].iter().enumerate() {
        if class_hits >= min_hits
            && (class_hits as f64 / total_kmers as f64) >= min_fraction
        {
            matching_mask |= 1u64 << class_idx;
            match_count += 1;
            single_match = class_idx;
        }
    }

    match match_count {
        0 => Classification::Unclassified,
        1 => Classification::Classified(single_match),
        _ => Classification::Ambiguous(matching_mask),
    }
}

/// Update the classification spinner with current throughput stats.
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

/// Classify one sequence and populate per-class hit counts.
fn classify_seq_kmers(
    seq: &[u8],
    hasher: &KmerHasher,
    kmer_length: u8,
    smer_length: u8,
    buffers: &mut Buffers,
    hits: &mut [u64; 64],
    num_classes: usize,
    index: &ClassificationIndex,
    min_hits: u64,
    min_fraction: f64,
) -> (usize, Classification) {
    fill_syncmers(seq, hasher, kmer_length, smer_length, buffers);

    for h in hits[..num_classes].iter_mut() {
        *h = 0;
    }

    let total_kmers = buffers.minimizers.len();

    match (&buffers.minimizers, index) {
        (MinimizerVec::U64(vec), ClassificationIndex::U64(map)) => {
            for &kmer in vec {
                if let Some(&bitmask) = map.get(&kmer) {
                    let mut bits = bitmask;
                    while bits != 0 {
                        let class_idx = bits.trailing_zeros() as usize;
                        hits[class_idx] += 1;
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
                        let class_idx = bits.trailing_zeros() as usize;
                        hits[class_idx] += 1;
                        bits &= bits - 1;
                    }
                }
            }
        }
        _ => panic!("Mismatch between MinimizerVec and ClassificationIndex types"),
    }

    let classification = classify_seq(
        hits,
        num_classes,
        total_kmers,
        min_hits,
        min_fraction,
    );

    (total_kmers, classification)
}

/// Shared state for summary classification mode.
#[derive(Clone)]
struct GlobalClassifyState {
    class_seqs: Vec<u64>,
    class_bases: Vec<u64>,
    ambiguous_seqs: u64,
    ambiguous_bases: u64,
    unclassified_seqs: u64,
    unclassified_bases: u64,
    stats: ProcessingStats,
}

impl GlobalClassifyState {
    fn new(num_classes: usize) -> Self {
        Self {
            class_seqs: vec![0; num_classes],
            class_bases: vec![0; num_classes],
            ambiguous_seqs: 0,
            ambiguous_bases: 0,
            unclassified_seqs: 0,
            unclassified_bases: 0,
            stats: ProcessingStats::default(),
        }
    }
}

/// Processor for summary classification mode.
#[derive(Clone)]
struct ClassifySummaryProcessor {
    kmer_length: u8,
    smer_length: u8,
    hasher: KmerHasher,
    index: Arc<ClassificationIndex>,
    num_classes: usize,
    min_hits: u64,
    min_fraction: f64,

    buffers: Buffers,
    hits: [u64; 64],

    // Thread-local accumulators.
    local_class_seqs: Vec<u64>,
    local_class_bases: Vec<u64>,
    local_ambiguous_seqs: u64,
    local_ambiguous_bases: u64,
    local_unclassified_seqs: u64,
    local_unclassified_bases: u64,
    local_stats: ProcessingStats,

    // Shared global state protected by a single lock.
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
        num_classes: usize,
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
            num_classes,
            min_hits,
            min_fraction,
            buffers,
            hits: [0u64; 64],
            local_class_seqs: vec![0; num_classes],
            local_class_bases: vec![0; num_classes],
            local_ambiguous_seqs: 0,
            local_ambiguous_bases: 0,
            local_unclassified_seqs: 0,
            local_unclassified_bases: 0,
            local_stats: ProcessingStats::default(),
            global: Arc::new(Mutex::new(GlobalClassifyState::new(num_classes))),
            spinner,
            start_time,
            limit_bp,
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for ClassifySummaryProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        // Stop once this sample reaches the base limit.
        if let Some(limit) = self.limit_bp {
            let global_bp = self.global.lock().stats.total_bp;
            if global_bp >= limit {
                return Err(paraseq::parallel::ProcessError::IoError(
                    std::io::Error::new(std::io::ErrorKind::Interrupted, "Sample limit reached"),
                ));
            }
        }

        let seq = record.seq();
        let seq_len = seq.len();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq_len as u64;

        let (_total_kmers, classification) = classify_seq_kmers(
            &seq, &self.hasher, self.kmer_length, self.smer_length,
            &mut self.buffers, &mut self.hits, self.num_classes,
            &self.index, self.min_hits, self.min_fraction,
        );

        match classification {
            Classification::Classified(class_idx) => {
                self.local_class_seqs[class_idx] += 1;
                self.local_class_bases[class_idx] += seq_len as u64;
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

        for i in 0..self.num_classes {
            g.class_seqs[i] += self.local_class_seqs[i];
            g.class_bases[i] += self.local_class_bases[i];
            self.local_class_seqs[i] = 0;
            self.local_class_bases[i] = 0;
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
            // Refresh spinner progress.
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

/// Processor for per-sequence classification output mode.
#[derive(Clone)]
struct ClassifyPerSeqProcessor {
    kmer_length: u8,
    smer_length: u8,
    hasher: KmerHasher,
    index: Arc<ClassificationIndex>,
    num_classes: usize,
    class_names: Arc<Vec<String>>,
    min_hits: u64,
    min_fraction: f64,
    sample_name: String,

    buffers: Buffers,
    hits: [u64; 64],

    // Thread-local output buffer.
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
        num_classes: usize,
        class_names: Arc<Vec<String>>,
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
            num_classes,
            class_names,
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
                    std::io::Error::new(std::io::ErrorKind::Interrupted, "Sample limit reached"),
                ));
            }
        }

        let seq = record.seq();
        let seq_id = String::from_utf8_lossy(record.id());
        let seq_len = seq.len();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq_len as u64;

        let (total_kmers, classification) = classify_seq_kmers(
            &seq, &self.hasher, self.kmer_length, self.smer_length,
            &mut self.buffers, &mut self.hits, self.num_classes,
            &self.index, self.min_hits, self.min_fraction,
        );

        use std::fmt::Write as FmtWrite;

        match classification {
            Classification::Classified(class_idx) => {
                let _ = writeln!(
                    self.local_output,
                    "{}\t{}\tclassified\t{}\t{}\t{}\t{}",
                    self.sample_name,
                    seq_id,
                    self.class_names[class_idx],
                    self.hits[class_idx],
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
                // Build comma-separated class names and hit counts.
                let mut classes = String::new();
                let mut hits_str = String::new();
                let mut bits = mask;
                let mut first = true;
                while bits != 0 {
                    let class_idx = bits.trailing_zeros() as usize;
                    if !first {
                        classes.push(',');
                        hits_str.push(',');
                    }
                    classes.push_str(&self.class_names[class_idx]);
                    let _ = write!(hits_str, "{}", self.hits[class_idx]);
                    first = false;
                    bits &= bits - 1;
                }

                let _ = writeln!(
                    self.local_output,
                    "{}\t{}\tambiguous\t{}\t{}\t{}\t{}",
                    self.sample_name, seq_id, classes, hits_str, total_kmers, seq_len,
                );
            }
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Flush buffered output into the shared writer.
        if !self.local_output.is_empty() {
            let mut writer = self.output_writer.lock();
            writer
                .write_all(&self.local_output)
                .map_err(paraseq::parallel::ProcessError::IoError)?;
            self.local_output.clear();
        }

        // Update shared progress stats.
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

    // Load a prebuilt index or build one from a class directory.
    let (index, class_names, kmer_length, smer_length) = if config.index_path.is_dir() {
        if !config.quiet {
            eprintln!(
                "Grate v{}; mode: classify (from directory); options: k={}, s={}, threads={}, min_hits={}, min_fraction={:.2}",
                version, config.kmer_length, config.smer_length, config.threads, config.min_hits, config.min_fraction
            );
        }

        let (index, class_names) = build_index_in_memory(
            &config.index_path,
            config.kmer_length,
            config.smer_length,
            config.threads,
            config.quiet,
        )?;

        (index, class_names, config.kmer_length, config.smer_length)
    } else {
        // Load prebuilt index.
        if !config.quiet {
            eprintln!(
                "Grate v{}; mode: classify (from index); options: threads={}, min_hits={}, min_fraction={:.2}",
                version, config.threads, config.min_hits, config.min_fraction
            );
        }

        let load_start = Instant::now();
        let (index, class_names, k, s) = load_index(&config.index_path)?;

        if !config.quiet {
            let elapsed = load_start.elapsed();
            eprintln!(
                "Index: {} k-mers, {} classes, k={}, s={} (loaded in {:.1}s)",
                index.len(),
                class_names.len(),
                k,
                s,
                elapsed.as_secs_f64()
            );
            for (i, name) in class_names.iter().enumerate() {
                eprintln!("  [{}] {}", i, name);
            }
        }

        (index, class_names, k, s)
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
    let class_names = Arc::new(class_names);
    let num_classes = class_names.len();

    // Process samples.
    use rayon::prelude::*;

    if config.per_seq {
        // Per-sequence output mode.
        let writer: Box<dyn Write + Send> = if let Some(path) = &config.output_path {
            Box::new(BufWriter::new(File::create(path)?))
        } else {
            Box::new(BufWriter::new(io::stdout()))
        };
        let writer = Arc::new(Mutex::new(BufWriter::new(writer)));

        // Write output header.
        {
            let mut w = writer.lock();
            writeln!(
                w,
                "sample\tseq_id\tclassification\tclass\thits\tseq_kmers\tseq_length"
            )?;
        }

        for (sample_paths, sample_name) in
            config.sample_paths.iter().zip(&config.sample_names)
        {
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
                    num_classes,
                    Arc::clone(&class_names),
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
        // Summary mode.
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
                    num_classes,
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

        // Write summary TSV.
        let writer: Box<dyn Write> = if let Some(path) = &config.output_path {
            Box::new(BufWriter::new(File::create(path)?))
        } else {
            Box::new(BufWriter::new(io::stdout()))
        };
        let mut writer = writer;

        writeln!(
            writer,
            "sample\tclass\tseqs_pct\tseqs\tbases_pct\tbases"
        )?;

        for (sample_name, result) in &sample_results {
            let total_seqs = result.total_seqs as f64;
            let total_bases = result.total_bases as f64;

            // Collect rows as (class_name, seqs, seq_pct, bases, base_pct).
            let mut rows: Vec<(&str, u64, f64, u64, f64)> = Vec::new();

            for (class_idx, counts) in result.class_counts.iter().enumerate() {
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
                rows.push((&class_names[class_idx], counts.seqs, pct_seqs, counts.bases, pct_bases));
            }

            // Add ambiguous row.
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
                rows.push(("ambiguous", result.ambiguous_seqs, pct_seqs, result.ambiguous_bases, pct_bases));
            }

            // Add unclassified row.
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
                rows.push(("unclassified", result.unclassified_seqs, pct_seqs, result.unclassified_bases, pct_bases));
            }

            // Sort by descending sequnce percentage.
            rows.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

            for (class_name, seqs, pct_seqs, bases, pct_bases) in &rows {
                writeln!(
                    writer,
                    "{}\t{}\t{:.2}\t{}\t{:.2}\t{}",
                    sample_name, class_name, pct_seqs, seqs, pct_bases, bases,
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

/// Process one sample in summary mode.
fn process_sample_summary(
    sample_paths: &[PathBuf],
    sample_name: &str,
    kmer_length: u8,
    smer_length: u8,
    index: &Arc<ClassificationIndex>,
    num_classes: usize,
    min_hits: u64,
    min_fraction: f64,
    threads: usize,
    quiet: bool,
    limit_bp: Option<u64>,
) -> Result<SampleClassificationResult> {
    let mut combined = SampleClassificationResult {
        class_counts: vec![ClassCounts::default(); num_classes],
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
            num_classes,
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

        // Merge this file's results into the sample totals.
        let g = processor.global.lock();
        for i in 0..num_classes {
            combined.class_counts[i].seqs += g.class_seqs[i];
            combined.class_counts[i].bases += g.class_bases[i];
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

        // Stop when the sample-level base limit is reached.
        if let Some(limit) = limit_bp
            && combined.total_bases >= limit
        {
            break;
        }
    }

    Ok(combined)
}

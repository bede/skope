//! Accelerated genome containment estimation using open syncmers or all k-mers
#[cfg(feature = "cli")]
pub mod classify;
#[cfg(feature = "cli")]
mod io;
pub mod kmers;
#[cfg(feature = "cli")]
pub mod length;
pub mod query;
pub mod stats;

use anyhow::Result;
use paraseq::Record;
use std::collections::HashSet;
use std::hash::BuildHasher;
use std::path::{Path, PathBuf};

// Re-export the main functionality
#[cfg(feature = "cli")]
pub use io::*;

#[cfg(feature = "cli")]
pub use query::{
    BuildQueryConfig, ContainmentConfig, read_query_index_meta, run_build_query, run_query,
};
pub use query::{QueryOptions, QuerySession, SortOrder};

#[cfg(feature = "cli")]
pub use length::{LengthHistogramConfig, run_lenhist};

#[cfg(feature = "cli")]
pub use classify::{BuildClassifyConfig, ClassifyConfig, run_build_classify, run_classification};

pub use kmers::{
    DEFAULT_KMER_LENGTH, DEFAULT_SMER_LENGTH, FracMinHash, Kdust, KmerVec, Kmers, calculate_kdust,
    decode_u64, decode_u128,
};

// ── Shared types ──────────────────────────────────────────────────────────────

/// BuildHasher using rapidhash with fixed seed for fast init
#[derive(Clone, Default)]
pub struct FixedRapidHasher;

impl BuildHasher for FixedRapidHasher {
    type Hasher = rapidhash::fast::RapidHasher<'static>;

    fn build_hasher(&self) -> Self::Hasher {
        rapidhash::fast::SeedableState::fixed().build_hasher()
    }
}

/// RapidHashSet using rapidhash with fixed seed for fast init
pub type RapidHashSet<T> = HashSet<T, FixedRapidHasher>;

// ── On-disk index header ────────────────────────────────────────────────────
// Every skope index begins with a shared prefix: 4-byte magic, 1-byte kind, 1-byte
// version. Detection reads the magic and kind; each kind interprets the rest by version.

/// Magic identifying any on-disk skope index
pub const INDEX_MAGIC: &[u8; 4] = b"SKPE";

/// Index kind, stored on disk as the single byte after the magic. Encode with `kind as u8`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum IndexKind {
    Classify = 0,
    Query = 1,
}

impl IndexKind {
    /// Parse the on-disk kind byte
    pub fn from_byte(byte: u8) -> Option<Self> {
        match byte {
            0 => Some(Self::Classify),
            1 => Some(Self::Query),
            _ => None,
        }
    }
}

impl std::fmt::Display for IndexKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Classify => "classification",
            Self::Query => "query",
        })
    }
}

/// Reject a CLI complexity threshold that disagrees with a prebuilt index's.
/// Exact comparison: a CLI-parsed f32 round-trips through the header bit-identically
pub fn check_index_complexity(cli: f32, index: f32) -> Result<()> {
    if cli > 0.0 && cli != index {
        return Err(anyhow::anyhow!(
            "--complexity {cli} conflicts with the index's complexity {index}; omit --complexity to use the index's"
        ));
    }
    Ok(())
}

/// Render the `index info` complexity line for either index kind
pub fn complexity_info_line(complexity: f32) -> String {
    if complexity > 0.0 {
        format!("  Complexity filter: {complexity} (kdust)")
    } else {
        "  Complexity filter: 0 (retain all)".to_string()
    }
}

// ── Shared utilities ──────────────────────────────────────────────────────────

#[derive(Clone, Default, Debug)]
pub struct ProcessingStats {
    pub total_seqs: u64,
    pub total_bp: u64,
    pub last_reported: u64,
}

pub fn format_bp(bp: u64) -> String {
    if bp >= 1_000_000_000 {
        format!("{:.1}Gbp", bp as f64 / 1_000_000_000.0)
    } else if bp >= 1_000_000 {
        format!("{:.1}Mbp", bp as f64 / 1_000_000.0)
    } else if bp >= 1_000 {
        format!("{:.1}Kbp", bp as f64 / 1_000.0)
    } else {
        format!("{}bp", bp)
    }
}

pub fn format_bp_per_sec(bp_per_sec: f64) -> String {
    if bp_per_sec >= 1_000_000_000.0 {
        format!("{:.1} Gbp/s", bp_per_sec / 1_000_000_000.0)
    } else if bp_per_sec >= 1_000_000.0 {
        format!("{:.1} Mbp/s", bp_per_sec / 1_000_000.0)
    } else if bp_per_sec >= 1_000.0 {
        format!("{:.1} Kbp/s", bp_per_sec / 1_000.0)
    } else {
        format!("{:.0} bp/s", bp_per_sec)
    }
}

/// How to read one input, from paths or streams. A path's format, fastx or CBQ, is detected on
/// opening
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Input<T = PathBuf> {
    /// Fastx, or CBQ with its own pairing
    File(T),
    /// Fastx of consecutive mates
    Interleaved(T),
    /// Separate fastx mate files
    Paired(T, T),
}

/// A processor fed one record, or a pair's mates together
pub(crate) trait SeqProcessor: Clone + Send {
    fn process(&mut self, id: &[u8], seqs: &[&[u8]]) -> paraseq::Result<()>;

    /// Merge thread-local state, once per batch
    fn flush(&mut self) -> paraseq::Result<()>;
}

/// Unwrap a processor's own error from paraseq's
pub(crate) fn processor_error(err: paraseq::Error) -> anyhow::Error {
    match err {
        paraseq::Error::Process(e) => anyhow::Error::from_boxed(e),
        e => e.into(),
    }
}

type StreamReader<'a> = paraseq::fastx::Reader<Box<dyn std::io::Read + 'a>>;

/// Fastx reader over a stream, decompressed as its magic dictates
fn stream_reader<'a>(stream: impl std::io::Read + 'a) -> Result<StreamReader<'a>> {
    let (stream, _format) = niffler::get_reader(Box::new(stream))?;
    Ok(paraseq::fastx::Reader::new(stream)?)
}

fn feed<P: SeqProcessor>(
    processor: &mut P,
    stats: &mut ProcessingStats,
    id: &[u8],
    seqs: &[&[u8]],
) -> Result<()> {
    stats.total_seqs += seqs.len() as u64;
    stats.total_bp += seqs.iter().map(|seq| seq.len() as u64).sum::<u64>();
    processor.process(id, seqs).map_err(processor_error)
}

/// Feed `processor` from fastx streams on one thread, flushing and reporting running totals once
/// per batch. Batches and pairing errors match paraseq's, as does `process_input`'s base limit
pub(crate) fn process_stream<'a, P: SeqProcessor>(
    input: Input<impl std::io::Read + 'a>,
    processor: &mut P,
    limit_bp: Option<u64>,
    mut on_batch: impl FnMut(&ProcessingStats),
) -> Result<ProcessingStats> {
    let mut stats = ProcessingStats::default();
    if limit_bp == Some(0) {
        return Ok(stats);
    }
    // Whether the batch reached the limit
    let mut flush = |processor: &mut P, stats: &ProcessingStats| -> Result<bool> {
        processor.flush().map_err(processor_error)?;
        on_batch(stats);
        Ok(limit_bp.is_some_and(|limit| stats.total_bp >= limit))
    };
    match input {
        Input::File(stream) => {
            let mut reader = stream_reader(stream)?;
            reader.update_batch_size_in_bp(1024 * 1024)?;
            let mut set = reader.new_record_set();
            while set.fill(&mut reader)? {
                for record in set.iter() {
                    let record = record?;
                    feed(processor, &mut stats, record.id(), &[&record.seq()])?;
                }
                if flush(processor, &stats)? {
                    break;
                }
            }
        }
        Input::Interleaved(stream) => {
            let mut reader = stream_reader(stream)?;
            let mut set = reader.new_record_set();
            while set.fill(&mut reader)? {
                let records = set.iter().collect::<paraseq::Result<Vec<_>>>()?;
                if records.len() % 2 != 0 {
                    Err(paraseq::Error::IncompatibleInterleavedSetSize(
                        records.len(),
                    ))?;
                }
                for [r1, r2] in records.as_chunks::<2>().0 {
                    feed(processor, &mut stats, r1.id(), &[&r1.seq(), &r2.seq()])?;
                }
                if flush(processor, &stats)? {
                    break;
                }
            }
        }
        Input::Paired(stream1, stream2) => {
            let (mut reader1, mut reader2) = (stream_reader(stream1)?, stream_reader(stream2)?);
            let (mut set1, mut set2) = (reader1.new_record_set(), reader2.new_record_set());
            loop {
                match (set1.fill(&mut reader1)?, set2.fill(&mut reader2)?) {
                    (false, false) => break,
                    (false, true) => Err(paraseq::Error::PairedRecordMismatch("R1"))?,
                    (true, false) => Err(paraseq::Error::PairedRecordMismatch("R2"))?,
                    (true, true) => {}
                }
                let records1 = set1.iter().collect::<paraseq::Result<Vec<_>>>()?;
                let records2 = set2.iter().collect::<paraseq::Result<Vec<_>>>()?;
                if records1.len() != records2.len() {
                    Err(paraseq::Error::IncompatibleRecordSetSizes(
                        records1.len(),
                        records2.len(),
                    ))?;
                }
                for (r1, r2) in records1.iter().zip(&records2) {
                    feed(processor, &mut stats, r1.id(), &[&r1.seq(), &r2.seq()])?;
                }
                if flush(processor, &stats)? {
                    break;
                }
            }
        }
    }
    Ok(stats)
}

/// Derive sample/group name from file path by stripping directory and extensions
pub fn derive_sample_name(path: &Path, is_directory: bool) -> String {
    // Handle stdin
    if path.to_string_lossy() == "-" {
        return "stdin".to_string();
    }

    // Get filename without directory
    let filename = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("unknown");

    // For directories, use the directory name directly
    if is_directory {
        return filename.to_string();
    }

    // For files, strip known extensions
    let mut name = filename.to_string();
    let extensions = [
        ".xz", ".gz", ".zst", ".fasta", ".fa", ".fastq", ".fq", ".cbq", ".cba",
    ];

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

    if name.is_empty() {
        filename.to_string()
    } else {
        name
    }
}

/// Strip R1 mate tags (_1, _R1, .r1, _R1_001) from the sample name
pub fn derive_mates_name(r1: &Path) -> String {
    let name = derive_sample_name(r1, false);
    let tagless = name.strip_suffix("_001").unwrap_or(&name);
    tagless
        .strip_suffix('1')
        .map(|s| s.strip_suffix(['R', 'r']).unwrap_or(s))
        .and_then(|s| s.strip_suffix(['_', '.']))
        .filter(|s| !s.is_empty())
        .map_or_else(|| name.clone(), str::to_string)
}

/// Validate the FracMinHash fraction is in (0, 1]
pub fn validate_fraction(fraction: f64) -> Result<()> {
    if !(fraction > 0.0 && fraction <= 1.0) {
        return Err(anyhow::anyhow!(
            "Invalid --fraction {fraction}: must be in (0, 1] (1 = retain all)"
        ));
    }
    Ok(())
}

/// Check the kdust threshold is in [0, 1], rejecting NaN and inf
pub fn validate_complexity(complexity: f32) -> Result<()> {
    if !(0.0..=1.0).contains(&complexity) {
        return Err(anyhow::anyhow!(
            "Invalid --complexity {complexity}: must be in [0, 1] (0 = retain all)"
        ));
    }
    Ok(())
}

/// Validate k-mer and s-mer size constraints (s = 0 selects every k-mer)
pub fn validate_k_s(kmer_length: u8, smer_length: u8) -> Result<()> {
    let k = kmer_length as usize;
    let s = smer_length as usize;

    // s = 0 bypasses syncmer selection, so the s constraints do not apply
    let s_ok = s == 0 || (s < k && s <= 32 && !s.is_multiple_of(2));
    if k > 61 || k.is_multiple_of(2) || !s_ok {
        return Err(anyhow::anyhow!(
            "Invalid k-s combination: k={}, s={} (constraints: k<=61, k odd, s odd, 1<=s<k, s<=32; s=0 takes every k-mer)",
            k,
            s
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_bp() {
        assert_eq!(format_bp(500), "500bp");
        assert_eq!(format_bp(1500), "1.5Kbp");
        assert_eq!(format_bp(1500000), "1.5Mbp");
        assert_eq!(format_bp(1500000000), "1.5Gbp");
    }

    #[test]
    fn test_derive_mates_name_drops_mate_tag() {
        for (r1, name) in [
            ("dir/sample_1.fq.gz", "sample"),
            ("sample_R1.fastq", "sample"),
            ("sample_R1_001.fastq.gz", "sample"),
            ("rsviruses17900.r1.fastq", "rsviruses17900"),
            ("liver1.fq", "liver1"),
            ("reads_001.fq", "reads_001"),
            ("_1.fq", "_1"),
        ] {
            assert_eq!(derive_mates_name(Path::new(r1)), name, "{r1}");
        }
    }

    #[test]
    fn test_format_bp_per_sec() {
        assert_eq!(format_bp_per_sec(500.0), "500 bp/s");
        assert_eq!(format_bp_per_sec(1500.0), "1.5 Kbp/s");
        assert_eq!(format_bp_per_sec(1500000.0), "1.5 Mbp/s");
        assert_eq!(format_bp_per_sec(1500000000.0), "1.5 Gbp/s");
    }
}

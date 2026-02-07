//! # Skope
//!
//! A fast open syncmer-based containment analysis tool for genomic sequences
//!
//! Skope analyses the containment of open syncmers from reference sequences in sequence datasets,
//! providing detailed statistics on containment and abundance per target sequence
pub mod classify;
pub mod query;
pub mod length;
pub mod minimizers;

use anyhow::Result;
use parking_lot::Mutex;
use std::collections::HashSet;
use std::hash::BuildHasher;
use std::path::{Path, PathBuf};
use std::sync::Arc;

// Re-export the main functionality
pub use query::{
    ContainmentConfig, ContainmentParameters, ContainmentResult, OutputFormat, Report,
    SampleResults, SortOrder, TargetInfo, TimingStats, TotalStats, run_containment_analysis,
};

pub use length::{
    LengthHistogramConfig, LengthHistogramParameters, LengthHistogramReport, LengthHistogramResult,
    run_length_histogram_analysis,
};

pub use classify::{BuildConfig, ClassifyConfig, build_classification_index, run_classification};

pub use minimizers::{
    Buffers, DEFAULT_KMER_LENGTH, DEFAULT_SMER_LENGTH, KmerHasher, MinimizerVec, decode_u64,
    decode_u128, fill_syncmers, fill_syncmers_with_positions,
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

// ── Shared utilities ──────────────────────────────────────────────────────────

#[derive(Clone, Default, Debug)]
pub struct ProcessingStats {
    pub total_seqs: u64,
    pub total_bp: u64,
    pub last_reported: u64,
}

pub const SAMPLE_LIMIT_REACHED_MSG: &str = "sample base limit reached";

pub fn sample_limit_reached_io_error() -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::Interrupted, SAMPLE_LIMIT_REACHED_MSG)
}

fn is_sample_limit_error(err: &paraseq::parallel::ProcessError) -> bool {
    match err {
        paraseq::parallel::ProcessError::IoError(io_err) => {
            io_err.kind() == std::io::ErrorKind::Interrupted
                && io_err.to_string() == SAMPLE_LIMIT_REACHED_MSG
        }
        _ => false,
    }
}

pub fn reader_with_inferred_batch_size(
    in_path: Option<&Path>,
) -> Result<paraseq::fastx::Reader<Box<dyn std::io::Read + Send>>> {
    let mut reader = paraseq::fastx::Reader::from_optional_path(in_path)?;
    reader.update_batch_size_in_bp(256 * 1024)?;
    Ok(reader)
}

pub fn format_bp(bp: usize) -> String {
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

/// Create a spinner progress bar for status display, or None if quiet
pub fn create_spinner(quiet: bool) -> Result<Option<Arc<Mutex<indicatif::ProgressBar>>>> {
    use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};

    if quiet {
        Ok(None)
    } else {
        let pb = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr());
        pb.set_style(
            ProgressStyle::default_spinner()
                .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                .template("{msg}")?,
        );
        Ok(Some(Arc::new(Mutex::new(pb))))
    }
}

/// Treat sample-limit interruptions as normal completion
pub fn handle_process_result(result: std::result::Result<(), paraseq::parallel::ProcessError>) -> Result<()> {
    match result {
        Ok(()) => Ok(()),
        Err(e) if is_sample_limit_error(&e) => Ok(()),
        Err(e) => Err(anyhow::anyhow!(e)),
    }
}

/// Find all fastx files in a directory (non-recursive, following symlinks)
pub fn find_fastx_files(dir_path: &Path) -> Result<Vec<PathBuf>> {
    const FASTX_EXTENSIONS: &[&str] = &[
        ".fasta", ".fa", ".fastq", ".fq",
        ".fasta.gz", ".fa.gz", ".fastq.gz", ".fq.gz",
        ".fasta.xz", ".fa.xz", ".fastq.xz", ".fq.xz",
        ".fasta.zst", ".fa.zst", ".fastq.zst", ".fq.zst",
    ];

    let mut files = Vec::new();

    let entries = std::fs::read_dir(dir_path)
        .map_err(|e| anyhow::anyhow!("Failed to read directory: {}: {}", dir_path.display(), e))?;

    for entry in entries {
        let entry = entry
            .map_err(|e| anyhow::anyhow!("Failed to read directory entry in {}: {}", dir_path.display(), e))?;

        let path = entry.path();

        // Skip hidden files
        if let Some(name) = path.file_name().and_then(|n| n.to_str())
            && name.starts_with('.')
        {
            continue;
        }

        // Follow symlinks
        let metadata = std::fs::metadata(&path)
            .map_err(|e| anyhow::anyhow!("Failed to access: {}: {}", path.display(), e))?;

        if !metadata.is_file() {
            continue;
        }

        let path_str = path.to_string_lossy().to_lowercase();
        if FASTX_EXTENSIONS.iter().any(|ext| path_str.ends_with(ext)) {
            files.push(path);
        }
    }

    if files.is_empty() {
        return Err(anyhow::anyhow!(
            "Directory contains no fastx files: {}. Expected files with extensions: {}",
            dir_path.display(),
            FASTX_EXTENSIONS.join(", ")
        ));
    }

    files.sort();
    Ok(files)
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
    let extensions = [".xz", ".gz", ".zst", ".fasta", ".fa", ".fastq", ".fq"];

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

/// Validate k-mer and s-mer size constraints for open syncmers
pub fn validate_k_s(kmer_length: u8, smer_length: u8) -> Result<()> {
    let k = kmer_length as usize;
    let s = smer_length as usize;

    if k > 61 || s >= k || !(1..=32).contains(&s) || k.is_multiple_of(2) || s.is_multiple_of(2) {
        return Err(anyhow::anyhow!(
            "Invalid k-s combination: k={}, s={} (constraints: k<=61, k odd, s odd, 1<=s<k, s<=32)",
            k, s
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
    fn test_format_bp_per_sec() {
        assert_eq!(format_bp_per_sec(500.0), "500 bp/s");
        assert_eq!(format_bp_per_sec(1500.0), "1.5 Kbp/s");
        assert_eq!(format_bp_per_sec(1500000.0), "1.5 Mbp/s");
        assert_eq!(format_bp_per_sec(1500000000.0), "1.5 Gbp/s");
    }
}

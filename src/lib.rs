//! Accelerated genome containment estimation using syncmers
pub mod classify;
pub mod length;
pub mod query;
pub mod stats;
pub mod syncmers;

use anyhow::Result;
use parking_lot::Mutex;
use std::collections::HashSet;
use std::hash::BuildHasher;
use std::path::{Path, PathBuf};
use std::sync::Arc;

// Re-export the main functionality
pub use query::{
    BuildQueryConfig, ContainmentConfig, SortOrder, read_query_index_meta, run_build_query,
    run_query,
};

pub use length::{LengthHistogramConfig, run_lenhist};

pub use classify::{BuildClassifyConfig, ClassifyConfig, run_build_classify, run_classification};

pub use syncmers::{
    Buffers, DEFAULT_KMER_LENGTH, DEFAULT_SMER_LENGTH, FracMinHash, Kdust, KmerHasher, SyncmerVec,
    calculate_kdust, decode_u64, decode_u128, fill_syncmers, fill_syncmers_with_positions,
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

/// Index kind if `path` is a skope index, else `None`
pub fn read_index_kind(path: &Path) -> Option<IndexKind> {
    use std::io::Read;
    let mut buf = [0u8; 5];
    std::fs::File::open(path).ok()?.read_exact(&mut buf).ok()?;
    if &buf[..4] == INDEX_MAGIC {
        IndexKind::from_byte(buf[4])
    } else {
        None
    }
}

/// Effective k/s for `targets`: a prebuilt index supplies its own, noting any conflicting CLI
/// values, otherwise CLI values fall back to defaults and are validated
pub fn resolve_k_s(targets: &Path, k: Option<u8>, s: Option<u8>, quiet: bool) -> Result<(u8, u8)> {
    if read_index_kind(targets).is_some() {
        // Both index kinds share the header prefix (magic, kind, version, k, s)
        use std::io::Read;
        let mut buf = [0u8; 8];
        std::fs::File::open(targets)?.read_exact(&mut buf)?;
        let (index_k, index_s) = (buf[6], buf[7]);
        let conflicts = k.is_some_and(|k| k != index_k) || s.is_some_and(|s| s != index_s);
        if conflicts && !quiet {
            eprintln!("Note: using k={index_k}, s={index_s} from index (CLI k/s ignored)");
        }
        return Ok((index_k, index_s));
    }
    let (k, s) = (
        k.unwrap_or(DEFAULT_KMER_LENGTH),
        s.unwrap_or(DEFAULT_SMER_LENGTH),
    );
    validate_k_s(k, s)?;
    Ok((k, s))
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

/// Print metadata for a skope index of either kind (`skope index info`)
pub fn run_index_info(path: &Path) -> Result<()> {
    let start = std::time::Instant::now();
    match read_index_kind(path) {
        Some(IndexKind::Query) => query::print_query_index_info(path)?,
        Some(IndexKind::Classify) => classify::print_classification_index_info(path)?,
        None => {
            return Err(anyhow::anyhow!("{} is not a skope index", path.display()));
        }
    }
    eprintln!("Loaded index info in {:.2?}", start.elapsed());
    Ok(())
}

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

/// Read one target path, treating `-` as stdin
pub fn reader_for_path(
    path: &Path,
) -> Result<paraseq::fastx::Reader<Box<dyn std::io::Read + Send>>> {
    let in_path = if path.to_string_lossy() == "-" {
        None
    } else {
        Some(path)
    };
    reader_with_inferred_batch_size(in_path)
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
pub fn handle_process_result(
    result: std::result::Result<(), paraseq::parallel::ProcessError>,
) -> Result<()> {
    match result {
        Ok(()) => Ok(()),
        Err(e) if is_sample_limit_error(&e) => Ok(()),
        Err(e) => Err(anyhow::anyhow!(e)),
    }
}

const FASTX_EXTENSIONS: &[&str] = &[
    ".fasta",
    ".fa",
    ".fastq",
    ".fq",
    ".fasta.gz",
    ".fa.gz",
    ".fastq.gz",
    ".fq.gz",
    ".fasta.xz",
    ".fa.xz",
    ".fastq.xz",
    ".fq.xz",
    ".fasta.zst",
    ".fa.zst",
    ".fastq.zst",
    ".fq.zst",
];

/// Check whether a path has a recognised fastx extension
pub fn is_fastx_file(path: &Path) -> bool {
    let path_str = path.to_string_lossy().to_lowercase();
    FASTX_EXTENSIONS.iter().any(|ext| path_str.ends_with(ext))
}

/// One sequence group discovered under a target/class directory:
/// either a single top-level FASTX file or all FASTX files directly inside a subdirectory.
#[derive(Debug, Clone)]
pub struct TargetGroup {
    pub name: String,
    pub files: Vec<PathBuf>,
}

/// Discover sequence groups under a directory, one per top-level FASTX file or subdirectory.
///
/// - Top-level FASTX files are returned first (sorted), then top-level subdirectories (sorted).
/// - Each subdirectory must contain at least one FASTX file directly inside it; nested
///   sub-subdirectories are not allowed and are rejected with an error.
/// - Hidden top-level entries and hidden files inside subdirectories are skipped.
/// - Duplicate derived group names (including `foo.fa` colliding with `foo/`) are rejected.
pub fn discover_target_groups(dir_path: &Path) -> Result<Vec<TargetGroup>> {
    let mut top_files: Vec<PathBuf> = Vec::new();
    let mut top_subdirs: Vec<PathBuf> = Vec::new();

    let entries = std::fs::read_dir(dir_path)
        .map_err(|e| anyhow::anyhow!("Failed to read directory: {}: {}", dir_path.display(), e))?;

    for entry in entries {
        let entry = entry.map_err(|e| {
            anyhow::anyhow!(
                "Failed to read directory entry in {}: {}",
                dir_path.display(),
                e
            )
        })?;
        let path = entry.path();

        if let Some(name) = path.file_name().and_then(|n| n.to_str())
            && name.starts_with('.')
        {
            continue;
        }

        let metadata = std::fs::metadata(&path)
            .map_err(|e| anyhow::anyhow!("Failed to access: {}: {}", path.display(), e))?;

        if metadata.is_dir() {
            top_subdirs.push(path);
        } else if metadata.is_file() && is_fastx_file(&path) {
            top_files.push(path);
        }
    }

    top_files.sort();
    top_subdirs.sort();

    if top_files.is_empty() && top_subdirs.is_empty() {
        return Err(anyhow::anyhow!(
            "Directory contains no fastx files or subdirectories: {}",
            dir_path.display()
        ));
    }

    let mut groups: Vec<TargetGroup> = Vec::with_capacity(top_files.len() + top_subdirs.len());

    for file in top_files {
        groups.push(TargetGroup {
            name: derive_sample_name(&file, false),
            files: vec![file],
        });
    }

    for subdir in top_subdirs {
        let mut subdir_files: Vec<PathBuf> = Vec::new();

        let entries = std::fs::read_dir(&subdir).map_err(|e| {
            anyhow::anyhow!("Failed to read directory: {}: {}", subdir.display(), e)
        })?;

        for entry in entries {
            let entry = entry.map_err(|e| {
                anyhow::anyhow!(
                    "Failed to read directory entry in {}: {}",
                    subdir.display(),
                    e
                )
            })?;
            let path = entry.path();

            if let Some(name) = path.file_name().and_then(|n| n.to_str())
                && name.starts_with('.')
            {
                continue;
            }

            let metadata = std::fs::metadata(&path)
                .map_err(|e| anyhow::anyhow!("Failed to access: {}: {}", path.display(), e))?;

            if metadata.is_dir() {
                return Err(anyhow::anyhow!(
                    "Nested subdirectory not allowed in group directory '{}': {}. Group directories must contain only fastx files.",
                    subdir.display(),
                    path.display()
                ));
            }

            if metadata.is_file() && is_fastx_file(&path) {
                subdir_files.push(path);
            }
        }

        if subdir_files.is_empty() {
            return Err(anyhow::anyhow!(
                "Group directory contains no fastx files: {}",
                subdir.display()
            ));
        }

        subdir_files.sort();
        groups.push(TargetGroup {
            name: derive_sample_name(&subdir, true),
            files: subdir_files,
        });
    }

    let mut seen = HashSet::new();
    for group in &groups {
        if !seen.insert(group.name.clone()) {
            return Err(anyhow::anyhow!(
                "Duplicate group name '{}' derived from directory entries. Top-level files and subdirectories must produce unique names (e.g. avoid 'foo.fa' next to 'foo/').",
                group.name
            ));
        }
    }

    Ok(groups)
}

/// Whether a path is a FIFO or process substitution
#[cfg(unix)]
pub fn is_special_input_path(path: &Path) -> bool {
    use std::os::unix::fs::FileTypeExt;
    let path_str = path.to_string_lossy();
    if path_str.starts_with("/dev/fd/") || path_str.starts_with("/proc/self/fd/") {
        return true;
    }
    if let Ok(metadata) = std::fs::metadata(path) {
        return metadata.file_type().is_fifo();
    }
    false
}

#[cfg(not(unix))]
pub fn is_special_input_path(_path: &Path) -> bool {
    false
}

/// Cheap sniff for fastx content: a plain `>`/`@` header, or a compression magic. Extensions
/// are too varied to check (`.fna` and friends are legitimate), so this reads the first bytes,
/// purely to give a clear error for non-sequence files.
fn looks_like_fastx(path: &Path) -> bool {
    use std::io::Read;
    let mut buf = [0u8; 8];
    let Ok(mut file) = std::fs::File::open(path) else {
        return true; // let the reader report the real I/O error
    };
    let Ok(n) = file.read(&mut buf) else {
        return true;
    };
    let buf = &buf[..n];

    const MAGICS: [&[u8]; 4] = [
        &[0x1f, 0x8b],                   // gzip
        &[0x28, 0xb5, 0x2f, 0xfd],       // zstd
        &[0xfd, b'7', b'z', b'X', b'Z'], // xz
        b"BZh",                          // bzip2
    ];
    if MAGICS.iter().any(|magic| buf.starts_with(magic)) {
        return true;
    }

    matches!(
        buf.iter().find(|b| !b.is_ascii_whitespace()),
        Some(b'>') | Some(b'@')
    )
}

/// What a `<TARGETS>` positional resolved to. `File` and `Directory` differ because
/// `--individual` splits the records of a file, but a directory's children already group it.
#[derive(Debug, Clone)]
pub enum TargetSource {
    /// Prebuilt `.sk` index of the expected kind
    Index(PathBuf),
    /// Fastx file or stdin, merged unless `--individual`
    File(TargetGroup),
    /// Directory grouped by top-level file or subdirectory
    Directory(Vec<TargetGroup>),
}

impl TargetSource {
    /// Groups to build, empty for a prebuilt index
    pub fn groups(&self) -> &[TargetGroup] {
        match self {
            Self::Index(_) => &[],
            Self::File(group) => std::slice::from_ref(group),
            Self::Directory(groups) => groups,
        }
    }

    /// Whether `--individual` should split records into groups
    pub fn splits_records(&self, individual: bool) -> bool {
        individual && matches!(self, Self::File(_))
    }
}

/// Name of the group formed from stdin
pub const STDIN_TARGET_NAME: &str = "stdin";

/// Whether a subcommand can take its targets from stdin. Only `index build-*` can:
/// elsewhere the samples already claim `-`, and there is only one stdin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StdinTargets {
    Accept,
    Reject,
}

/// Resolve `<TARGETS>` for a subcommand expecting index kind `expect`
///
/// Accepts, in order: `-` (only where `stdin` is [`StdinTargets::Accept`]), a prebuilt `.sk`
/// index, a directory (one group per top-level fastx file or subdirectory), or a single fastx
/// file (one group). Pipes are always rejected: paraseq cannot sniff compression on a
/// non-seekable path, so they fail later with a confusing parse error.
pub fn resolve_targets(
    path: &Path,
    expect: IndexKind,
    stdin: StdinTargets,
) -> Result<TargetSource> {
    if path.to_string_lossy() == "-" {
        if stdin == StdinTargets::Reject {
            return Err(anyhow::anyhow!(
                "Targets cannot be read from stdin; `-` is accepted only for samples. \
                 Give a fastx file, a directory of fastx files, or a skope {expect} index."
            ));
        }
        return Ok(TargetSource::File(TargetGroup {
            name: STDIN_TARGET_NAME.to_string(),
            files: vec![path.to_path_buf()],
        }));
    }

    if is_special_input_path(path) {
        return Err(anyhow::anyhow!(
            "Targets cannot be read from a pipe ({}). \
             Give a fastx file, a directory of fastx files, or a skope {expect} index.",
            path.display()
        ));
    }

    if !path.exists() {
        return Err(anyhow::anyhow!("Path does not exist: {}", path.display()));
    }

    if let Some(kind) = read_index_kind(path) {
        if kind != expect {
            return Err(anyhow::anyhow!(
                "{} is a skope {kind} index, not a {expect} index",
                path.display()
            ));
        }
        return Ok(TargetSource::Index(path.to_path_buf()));
    }

    if path.is_dir() {
        return Ok(TargetSource::Directory(discover_target_groups(path)?));
    }

    if path.is_file() {
        if !looks_like_fastx(path) {
            return Err(anyhow::anyhow!(
                "{} is not a fastx file, a directory of fastx files, or a skope {expect} index",
                path.display()
            ));
        }
        return Ok(TargetSource::File(TargetGroup {
            name: derive_sample_name(path, false),
            files: vec![path.to_path_buf()],
        }));
    }

    Err(anyhow::anyhow!(
        "{} is neither a fastx file, a directory of fastx files, nor a skope {expect} index",
        path.display()
    ))
}

/// Find all fastx files in a directory (non-recursive, following symlinks)
pub fn find_fastx_files(dir_path: &Path) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();

    let entries = std::fs::read_dir(dir_path)
        .map_err(|e| anyhow::anyhow!("Failed to read directory: {}: {}", dir_path.display(), e))?;

    for entry in entries {
        let entry = entry.map_err(|e| {
            anyhow::anyhow!(
                "Failed to read directory entry in {}: {}",
                dir_path.display(),
                e
            )
        })?;

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

        if is_fastx_file(&path) {
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

/// Recursive fastq hunting under a directory. Descends real subdirectories
/// only (symlinked dirs skipped, so cycles can't loop); hidden entries skipped.
pub fn find_fastx_files_recursive(dir_path: &Path) -> Result<Vec<PathBuf>> {
    fn walk(dir: &Path, out: &mut Vec<PathBuf>) -> Result<()> {
        for entry in std::fs::read_dir(dir)? {
            let entry = entry?;
            let path = entry.path();
            if path
                .file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with('.'))
            {
                continue;
            }
            if entry.file_type()?.is_dir() {
                walk(&path, out)?;
            } else if is_fastx_file(&path) {
                out.push(path);
            }
        }
        Ok(())
    }

    let mut files = Vec::new();
    walk(dir_path, &mut files)?;
    files.sort();
    if files.is_empty() {
        return Err(anyhow::anyhow!(
            "Directory contains no fastx files (recursively): {}",
            dir_path.display()
        ));
    }
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

    if name.is_empty() {
        filename.to_string()
    } else {
        name
    }
}

/// Validate k-mer and s-mer size constraints for open syncmers
pub fn validate_k_s(kmer_length: u8, smer_length: u8) -> Result<()> {
    let k = kmer_length as usize;
    let s = smer_length as usize;

    if k > 61 || s >= k || !(1..=32).contains(&s) || k.is_multiple_of(2) || s.is_multiple_of(2) {
        return Err(anyhow::anyhow!(
            "Invalid k-s combination: k={}, s={} (constraints: k<=61, k odd, s odd, 1<=s<k, s<=32)",
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
    fn test_format_bp_per_sec() {
        assert_eq!(format_bp_per_sec(500.0), "500 bp/s");
        assert_eq!(format_bp_per_sec(1500.0), "1.5 Kbp/s");
        assert_eq!(format_bp_per_sec(1500000.0), "1.5 Mbp/s");
        assert_eq!(format_bp_per_sec(1500000000.0), "1.5 Gbp/s");
    }
}

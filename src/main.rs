use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use std::collections::HashSet;
use std::path::{Path, PathBuf};

const DEFAULT_KMER_LENGTH: u8 = 31;
const DEFAULT_SMER_SIZE: u8 = 15;

/// Derive sample name from file path by stripping directory and extensions
fn derive_sample_name(path: &Path, is_directory: bool) -> String {
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

    // Strip known extensions in order (handles chained extensions like .fastq.gz)
    let extensions = [".xz", ".gz", ".zst", ".fasta", ".fa", ".fastq", ".fq"];

    loop {
        let original_len = name.len();
        for ext in &extensions {
            if name.ends_with(ext) {
                name = name[..name.len() - ext.len()].to_string();
                break;
            }
        }
        // If no extension was removed, we're done
        if name.len() == original_len {
            break;
        }
    }

    // Fallback if everything was stripped
    if name.is_empty() {
        filename.to_string()
    } else {
        name
    }
}

/// Validate that sample names are unique
fn validate_sample_names(names: &[String]) -> Result<()> {
    let mut seen = HashSet::new();
    let mut duplicates = Vec::new();

    for name in names {
        if !seen.insert(name) {
            if !duplicates.contains(name) {
                duplicates.push(name.clone());
            }
        }
    }

    if !duplicates.is_empty() {
        return Err(anyhow::anyhow!(
            "Duplicate sample names detected: {}. Please rename files or use --sample-names to specify unique names.",
            duplicates.join(", ")
        ));
    }

    Ok(())
}

/// Find all fastx files in a directory (non-recursive, following symlinks)
fn find_fastx_files_in_dir(dir_path: &Path) -> Result<Vec<PathBuf>> {
    // Extensions to match
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

    let mut fastx_files = Vec::new();

    // Read directory entries
    let entries = std::fs::read_dir(dir_path)
        .with_context(|| format!("Failed to read directory: {}", dir_path.display()))?;

    for entry in entries {
        let entry = entry
            .with_context(|| format!("Failed to read directory entry in {}", dir_path.display()))?;

        let path = entry.path();

        // Skip hidden files (starting with '.')
        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            if name.starts_with('.') {
                continue;
            }
        }

        // Follow symlinks via metadata()
        let metadata = std::fs::metadata(&path)
            .with_context(|| format!("Failed to access: {}", path.display()))?;

        // Only process files (not subdirectories)
        if !metadata.is_file() {
            continue;
        }

        // Check if file has fastx extension (case-insensitive)
        let path_str = path.to_string_lossy().to_lowercase();
        if FASTX_EXTENSIONS.iter().any(|ext| path_str.ends_with(ext)) {
            fastx_files.push(path);
        }
    }

    // Error if no fastx files found
    if fastx_files.is_empty() {
        return Err(anyhow::anyhow!(
            "Directory contains no fastx files: {}. Expected files with extensions: {}",
            dir_path.display(),
            FASTX_EXTENSIONS.join(", ")
        ));
    }

    // Sort for deterministic ordering
    fastx_files.sort();

    Ok(fastx_files)
}

/// Expand sample inputs (files and directories) into lists of files per sample
fn expand_sample_inputs(inputs: &[PathBuf]) -> Result<(Vec<Vec<PathBuf>>, Vec<bool>)> {
    let mut expanded_samples = Vec::new();
    let mut is_directory = Vec::new();

    for input in inputs {
        // Check stdin special case
        if input.to_string_lossy() == "-" {
            expanded_samples.push(vec![input.clone()]);
            is_directory.push(false);
            continue;
        }

        // Check path exists
        if !input.exists() {
            return Err(anyhow::anyhow!("Path does not exist: {}", input.display()));
        }

        let metadata = std::fs::metadata(input)
            .with_context(|| format!("Failed to access path: {}", input.display()))?;

        if metadata.is_file() {
            expanded_samples.push(vec![input.clone()]);
            is_directory.push(false);
        } else if metadata.is_dir() {
            let files = find_fastx_files_in_dir(input)?;
            expanded_samples.push(files);
            is_directory.push(true);
        } else {
            return Err(anyhow::anyhow!(
                "Path is neither a regular file nor directory: {}",
                input.display()
            ));
        }
    }

    Ok((expanded_samples, is_directory))
}

/// Parse sample string with K/M/G/T suffix into bp count
fn parse_sample(s: &str) -> Result<u64> {
    let s = s.trim().to_uppercase();
    let (num_str, multiplier) = if s.ends_with('T') {
        (&s[..s.len() - 1], 1_000_000_000_000u64)
    } else if s.ends_with('G') {
        (&s[..s.len() - 1], 1_000_000_000u64)
    } else if s.ends_with('M') {
        (&s[..s.len() - 1], 1_000_000u64)
    } else if s.ends_with('K') {
        (&s[..s.len() - 1], 1_000u64)
    } else {
        (s.as_str(), 1u64)
    };

    let num: u64 = num_str
        .parse()
        .with_context(|| format!("Invalid sample value: {}", s))?;

    Ok(num * multiplier)
}

#[derive(Parser)]
#[command(author, version, about = "Streaming containment and abundance estimation using open syncmers", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Calculate open syncmer containment & abundance in fastx files or directories thereof
    Con {
        /// Path to fasta file containing target sequence record(s)
        targets: PathBuf,

        /// Path(s) to fastx files/dirs (- for stdin). Each file/dir is treated as a separate sample.
        #[arg(required = true)]
        samples: Vec<PathBuf>,

        // Algorithm parameters
        /// Syncmer length (1-61)
        #[arg(short = 'k', long = "kmer-length", default_value_t = DEFAULT_KMER_LENGTH, value_parser = clap::value_parser!(u8).range(1..=61))]
        kmer_length: u8,

        /// S-mer size for open syncmer selection (s < k, s must be odd)
        #[arg(short = 's', long = "smer-size", default_value_t = DEFAULT_SMER_SIZE)]
        smer_size: u8,

        /// Comma-separated abundance thresholds for containment calculation
        #[arg(
            short = 'a',
            long = "abundance-thresholds",
            value_delimiter = ',',
            default_value = "10"
        )]
        abundance_thresholds: Vec<usize>,

        /// Consider only syncmers unique to each target
        #[arg(short = 'd', long = "discriminatory", default_value_t = false)]
        discriminatory: bool,

        // Processing options
        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Terminate read processing after approximately this many bases (e.g. 50M, 10G)
        #[arg(short = 'l', long = "limit")]
        limit: Option<String>,

        // Output options
        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Output format
        #[arg(short = 'f', long = "format", default_value = "table", value_parser = ["table", "csv", "json"])]
        format: String,

        /// Comma-separated sample names (default is file/dir name without extension)
        #[arg(short = 'n', long = "names", value_delimiter = ',')]
        sample_names: Option<Vec<String>>,

        /// Sort displayed results: o=original, t=target, S=sample, c=containment (descending)
        #[arg(short = 'S', long = "sort", default_value = "o", value_parser = ["o", "t", "s", "c"])]
        sort: String,

        /// Suppress progress reporting
        #[arg(short = 'q', long = "quiet", default_value_t = false)]
        quiet: bool,
    },
    /// Generate length histogram for reads with one or more syncmer hits to target sequences
    Len {
        /// Path to fasta file containing target sequence record(s) (- to disable target filtering)
        targets: PathBuf,

        /// Path(s) to fastx files/dirs (- for stdin). Each file/dir is treated as a separate sample.
        #[arg(required = true)]
        samples: Vec<PathBuf>,

        // Algorithm parameters
        /// Syncmer length (1-61)
        #[arg(short = 'k', long = "kmer-length", default_value_t = DEFAULT_KMER_LENGTH, value_parser = clap::value_parser!(u8).range(1..=61))]
        kmer_length: u8,

        /// S-mer size for open syncmer selection (s < k, s must be odd)
        #[arg(short = 's', long = "smer-size", default_value_t = DEFAULT_SMER_SIZE)]
        smer_size: u8,

        // Processing options
        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Terminate read processing after approximately this many bases (e.g. 50M, 10G)
        #[arg(short = 'l', long = "limit")]
        limit: Option<String>,

        // Output options
        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Comma-separated sample names (default is file/dir name without extension)
        #[arg(short = 'n', long = "names", value_delimiter = ',')]
        sample_names: Option<Vec<String>>,

        /// Suppress progress reporting
        #[arg(short = 'q', long = "quiet", default_value_t = false)]
        quiet: bool,
    },
}

fn main() -> Result<()> {
    // Check we have either AVX2 or NEON for SIMD acceleration
    #[cfg(not(any(target_feature = "avx2", target_feature = "neon")))]
    {
        eprintln!(
            "Warning: SIMD acceleration is unavailable. For best performance, compile with `cargo build --release -C target-cpu=native`"
        );
    }

    let cli = Cli::parse();

    match &cli.command {
        Commands::Con {
            targets,
            samples,
            sample_names,
            kmer_length,
            smer_size,
            threads,
            output,
            quiet,
            format,
            abundance_thresholds,
            discriminatory,
            limit,
            sort,
        } => {
            // Expand directories to lists of files
            let (expanded_reads, is_directory) = expand_sample_inputs(&samples)?;

            // Derive or validate sample names
            let derived_sample_names: Vec<String> = if let Some(names) = sample_names {
                // User-provided names
                if names.len() != expanded_reads.len() {
                    return Err(anyhow::anyhow!(
                        "Number of sample names ({}) must match number of samples ({})",
                        names.len(),
                        expanded_reads.len()
                    ));
                }
                names.clone()
            } else {
                // Derive from filenames/directory names
                samples
                    .iter()
                    .zip(&is_directory)
                    .map(|(p, &is_dir)| derive_sample_name(p, is_dir))
                    .collect()
            };

            // Validate uniqueness
            validate_sample_names(&derived_sample_names)?;
            // Validate k-mer and s-mer size constraints for open syncmers
            let k = *kmer_length as usize;
            let s = *smer_size as usize;

            // Check constraints:
            // - k <= 61 (fits in packed representation)
            // - 1 <= s < k (valid s-mer within k-mer)
            // - s <= 32 (s-mer must fit in hasher's u64 representation)
            // - k must be odd (for canonical strand determination)
            // - s must be odd (for open syncmers, w = k - s + 1 must be odd)
            if k > 61 || s >= k || s < 1 || s > 32 || k % 2 == 0 || s % 2 == 0 {
                return Err(anyhow::anyhow!(
                    "Invalid k-s combination: k={}, s={} (constraints: k<=61, k odd, s odd, 1<=s<k, s<=32)",
                    k,
                    s
                ));
            }

            // Configure thread pool if specified (non-zero)
            if *threads > 0 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(*threads)
                    .build_global()
                    .context("Failed to initialize thread pool")?;
            }

            // Parse outfmt
            let output_format = match format.as_str() {
                "table" => grate::OutputFormat::Table,
                "csv" => grate::OutputFormat::Csv,
                "json" => grate::OutputFormat::Json,
                _ => unreachable!("clap should have validated the format"),
            };

            // Parse limit if provided
            let limit_bp = if let Some(s) = limit {
                Some(parse_sample(s)?)
            } else {
                None
            };

            // Parse sort order
            let sort_order = match sort.as_str() {
                "o" => grate::SortOrder::Original,
                "t" => grate::SortOrder::Target,
                "s" => grate::SortOrder::Sample,
                "c" => grate::SortOrder::Containment,
                _ => unreachable!("clap should have validated the sort order"),
            };

            let config = grate::ContainmentConfig {
                targets_path: targets.clone(),
                reads_paths: expanded_reads,
                sample_names: derived_sample_names,
                kmer_length: *kmer_length,
                smer_size: *smer_size,
                threads: *threads,
                output_path: if output == "-" {
                    None
                } else {
                    Some(PathBuf::from(output))
                },
                quiet: *quiet,
                output_format,
                abundance_thresholds: abundance_thresholds.clone(),
                discriminatory: *discriminatory,
                limit_bp,
                sort_order,
            };

            config
                .execute()
                .context("Failed to run containment analysis")?;
        }
        Commands::Len {
            targets,
            samples,
            sample_names,
            kmer_length,
            smer_size,
            threads,
            output,
            quiet,
            limit,
        } => {
            // Expand directories to lists of files
            let (expanded_reads, is_directory) = expand_sample_inputs(&samples)?;

            // Derive or validate sample names
            let derived_sample_names: Vec<String> = if let Some(names) = sample_names {
                // User-provided names
                if names.len() != expanded_reads.len() {
                    return Err(anyhow::anyhow!(
                        "Number of sample names ({}) must match number of samples ({})",
                        names.len(),
                        expanded_reads.len()
                    ));
                }
                names.clone()
            } else {
                // Derive from filenames/directory names
                samples
                    .iter()
                    .zip(&is_directory)
                    .map(|(p, &is_dir)| derive_sample_name(p, is_dir))
                    .collect()
            };

            // Validate uniqueness
            validate_sample_names(&derived_sample_names)?;

            // Validate k-mer and s-mer size constraints for open syncmers
            let k = *kmer_length as usize;
            let s = *smer_size as usize;

            // Check constraints:
            // - k <= 61 (fits in packed representation)
            // - 1 <= s < k (valid s-mer within k-mer)
            // - s <= 32 (s-mer must fit in hasher's u64 representation)
            // - k must be odd (for canonical strand determination)
            // - s must be odd (for open syncmers, w = k - s + 1 must be odd)
            if k > 61 || s >= k || s < 1 || s > 32 || k % 2 == 0 || s % 2 == 0 {
                return Err(anyhow::anyhow!(
                    "Invalid k-s combination: k={}, s={} (constraints: k<=61, k odd, s odd, 1<=s<k, s<=32)",
                    k,
                    s
                ));
            }

            // Configure thread pool if specified (non-zero)
            if *threads > 0 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(*threads)
                    .build_global()
                    .context("Failed to initialize thread pool")?;
            }

            // Parse limit if provided
            let limit_bp = if let Some(s) = limit {
                Some(parse_sample(s)?)
            } else {
                None
            };

            // Detect if user wants all reads (no target filtering)
            let include_all_reads = targets.to_string_lossy() == "-";

            let config = grate::LengthHistogramConfig {
                targets_path: targets.clone(),
                reads_paths: expanded_reads,
                sample_names: derived_sample_names,
                kmer_length: *kmer_length,
                smer_size: *smer_size,
                threads: *threads,
                output_path: if output == "-" {
                    None
                } else {
                    Some(PathBuf::from(output))
                },
                quiet: *quiet,
                limit_bp,
                include_all_reads,
            };

            config
                .execute()
                .context("Failed to run length histogram analysis")?;
        }
    }

    Ok(())
}

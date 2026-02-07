use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use std::collections::HashSet;
use std::path::PathBuf;

use skope::{derive_sample_name, find_fastx_files, validate_k_s};

const DEFAULT_KMER_LENGTH: u8 = 31;
const DEFAULT_SMER_LENGTH: u8 = 9;

/// Validate that sample names are unique
fn validate_sample_names(names: &[String]) -> Result<()> {
    let mut seen = HashSet::new();
    let mut duplicates = Vec::new();

    for name in names {
        if !seen.insert(name)
            && !duplicates.contains(name) {
                duplicates.push(name.clone());
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
            let files = find_fastx_files(input)?;
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
#[command(author, version, about = "Containment and abundance estimation using open syncmers", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum IndexCommands {
    /// Build a classification index from a directory of FASTA files (one per group)
    Build {
        /// Directory containing FASTA files (one per group)
        groups: PathBuf,

        /// K-mer length (1-61, must be odd)
        #[arg(short = 'k', long = "kmer-length", default_value_t = DEFAULT_KMER_LENGTH, value_parser = clap::value_parser!(u8).range(1..=61))]
        kmer_length: u8,

        /// S-mer length used for open syncmer selection (s < k, s must be odd)
        #[arg(short = 's', long = "smer-length", default_value_t = DEFAULT_SMER_LENGTH)]
        smer_length: u8,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Path to output index file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Suppress progress reporting
        #[arg(short = 'q', long = "quiet", default_value_t = false)]
        quiet: bool,
    },
}

#[derive(Subcommand)]
enum Commands {
    /// Estimate k-mer containment & abundance in fastx file(s) or directories thereof
    Query {
        /// Path to fastx file containing target sequence record(s)
        targets: PathBuf,

        /// Path(s) to fastx files/dirs (- for stdin). Each file/dir is treated as a separate sample
        #[arg(required = true)]
        samples: Vec<PathBuf>,

        // Algorithm parameters
        /// K-mer length (1-61)
        #[arg(short = 'k', long = "kmer-length", default_value_t = DEFAULT_KMER_LENGTH, value_parser = clap::value_parser!(u8).range(1..=61))]
        kmer_length: u8,

        /// S-mer length used for open syncmer selection (s < k, s must be odd)
        #[arg(short = 's', long = "smer-length", default_value_t = DEFAULT_SMER_LENGTH)]
        smer_length: u8,

        /// Comma-separated additional abundance thresholds for containment estimation
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

        /// Terminate processing after approximately this many bases (e.g. 50M, 10G)
        #[arg(short = 'l', long = "limit")]
        limit: Option<String>,

        // Output options
        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Output format
        #[arg(short = 'f', long = "format", default_value = "tsv", value_parser = ["tsv", "table"])]
        format: String,

        /// Comma-separated sample names (default is file/dir name without extension)
        #[arg(short = 'n', long = "names", value_delimiter = ',')]
        sample_names: Option<Vec<String>>,

        /// Sort displayed results: c=containment (descending), t=target, s=sample, o=original
        #[arg(short = 'S', long = "sort", default_value = "c", value_parser = ["c", "t", "o", "s"])]
        sort: String,

        /// Suppress progress reporting
        #[arg(short = 'q', long = "quiet", default_value_t = false)]
        quiet: bool,

        /// Suppress TOTAL summary rows in output
        #[arg(long = "no-total", default_value_t = false)]
        no_total: bool,

        /// Dump open syncmer positions to TSV file (target\tposition)
        #[arg(long = "dump-positions")]
        dump_positions: Option<String>,
    },

    /// Classify sequences into groups based on k-mer membership
    Classify {
        /// Path to .idx index file or directory of FASTA files (one per group)
        index: PathBuf,

        /// Path(s) to fastx files/dirs (- for stdin)
        #[arg(required = true)]
        samples: Vec<PathBuf>,

        /// K-mer length (only used when index is a directory) (1-61, must be odd)
        #[arg(short = 'k', long = "kmer-length", default_value_t = DEFAULT_KMER_LENGTH, value_parser = clap::value_parser!(u8).range(1..=61))]
        kmer_length: u8,

        /// S-mer length (only used when index is a directory)
        #[arg(short = 's', long = "smer-length", default_value_t = DEFAULT_SMER_LENGTH)]
        smer_length: u8,

        /// Minimum k-mer hits to classify a sequence to a group
        #[arg(short = 'm', long = "min-hits", default_value_t = 1)]
        min_hits: u64,

        /// Minimum fraction of sequence k-mers hitting a group
        #[arg(short = 'r', long = "min-fraction", default_value_t = 0.0)]
        min_fraction: f64,

        /// Consider only k-mers unique to each group
        #[arg(short = 'd', long = "discriminatory", default_value_t = false)]
        discriminatory: bool,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Terminate processing after approximately this many bases (e.g. 50M, 10G)
        #[arg(short = 'l', long = "limit")]
        limit: Option<String>,

        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Output per-sequence classifications instead of summary
        #[arg(long = "per-seq", default_value_t = false)]
        per_seq: bool,

        /// Comma-separated sample names (default is file/dir name without extension)
        #[arg(short = 'n', long = "names", value_delimiter = ',')]
        sample_names: Option<Vec<String>>,

        /// Suppress progress reporting
        #[arg(short = 'q', long = "quiet", default_value_t = false)]
        quiet: bool,
    },

    /// Generate length histogram for sequences with k-mer hits to target sequence(s)
    Lenhist {
        /// Path to fastx file containing target sequence record(s) (- to disable target filtering)
        targets: PathBuf,

        /// Path(s) to fastx files/dirs (- for stdin). Each file/dir is treated as a separate sample
        #[arg(required = true)]
        samples: Vec<PathBuf>,

        // Algorithm parameters
        /// k-mer length (1-61)
        #[arg(short = 'k', long = "kmer-length", default_value_t = DEFAULT_KMER_LENGTH, value_parser = clap::value_parser!(u8).range(1..=61))]
        kmer_length: u8,

        /// S-mer length used for open syncmer selection (s < k, s must be odd)
        #[arg(short = 's', long = "smer-length", default_value_t = DEFAULT_SMER_LENGTH)]
        smer_length: u8,

        // Processing options
        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Terminate processing after approximately this many bases (e.g. 50M, 10G)
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

    /// Build and manage classification indexes
    #[command(subcommand)]
    Index(IndexCommands),
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
        Commands::Index(index_cmd) => match index_cmd {
            IndexCommands::Build {
                groups,
                kmer_length,
                smer_length,
                threads,
                output,
                quiet,
            } => {
                validate_k_s(*kmer_length, *smer_length)?;

                if !groups.is_dir() {
                    return Err(anyhow::anyhow!(
                        "Groups path must be a directory: {}",
                        groups.display()
                    ));
                }

                if *threads > 0 {
                    rayon::ThreadPoolBuilder::new()
                        .num_threads(*threads)
                        .build_global()
                        .context("Failed to initialise thread pool")?;
                }

                let config = skope::BuildConfig {
                    groups_dir: groups.clone(),
                    kmer_length: *kmer_length,
                    smer_length: *smer_length,
                    threads: *threads,
                    output_path: if output == "-" {
                        None
                    } else {
                        Some(PathBuf::from(output))
                    },
                    quiet: *quiet,
                };

                skope::build_classification_index(&config)
                    .context("Failed to build classification index")?;
            }
        },

        Commands::Classify {
            index,
            samples,
            sample_names,
            kmer_length,
            smer_length,
            min_hits,
            min_fraction,
            discriminatory,
            threads,
            limit,
            output,
            per_seq,
            quiet,
        } => {
            let (expanded_samples, is_directory) = expand_sample_inputs(samples)?;

            let derived_sample_names: Vec<String> = if let Some(names) = sample_names {
                if names.len() != expanded_samples.len() {
                    return Err(anyhow::anyhow!(
                        "Number of sample names ({}) must match number of samples ({})",
                        names.len(),
                        expanded_samples.len()
                    ));
                }
                names.clone()
            } else {
                samples
                    .iter()
                    .zip(&is_directory)
                    .map(|(p, &is_dir)| derive_sample_name(p, is_dir))
                    .collect()
            };

            validate_sample_names(&derived_sample_names)?;

            // Only validate k/s when index is a directory
            if index.is_dir() {
                validate_k_s(*kmer_length, *smer_length)?;
            }

            if *threads > 0 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(*threads)
                    .build_global()
                    .context("Failed to initialise thread pool")?;
            }

            let limit_bp = if let Some(s) = limit {
                Some(parse_sample(s)?)
            } else {
                None
            };

            let config = skope::ClassifyConfig {
                index_path: index.clone(),
                sample_paths: expanded_samples,
                sample_names: derived_sample_names,
                kmer_length: *kmer_length,
                smer_length: *smer_length,
                min_hits: *min_hits,
                min_fraction: *min_fraction,
                threads: *threads,
                limit_bp,
                output_path: if output == "-" {
                    None
                } else {
                    Some(PathBuf::from(output))
                },
                per_seq: *per_seq,
                discriminatory: *discriminatory,
                quiet: *quiet,
            };

            skope::run_classification(&config)
                .context("Failed to run classification")?;
        }

        Commands::Query {
            targets,
            samples,
            sample_names,
            kmer_length,
            smer_length,
            threads,
            output,
            quiet,
            format,
            abundance_thresholds,
            discriminatory,
            limit,
            sort,
            dump_positions,
            no_total,
        } => {
            // Expand directories to lists of files
            let (expanded_samples, is_directory) = expand_sample_inputs(samples)?;

            // Derive or validate sample names
            let derived_sample_names: Vec<String> = if let Some(names) = sample_names {
                // User-provided names
                if names.len() != expanded_samples.len() {
                    return Err(anyhow::anyhow!(
                        "Number of sample names ({}) must match number of samples ({})",
                        names.len(),
                        expanded_samples.len()
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
            validate_k_s(*kmer_length, *smer_length)?;

            // Configure thread pool if specified (non-zero)
            if *threads > 0 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(*threads)
                    .build_global()
                    .context("Failed to initialise thread pool")?;
            }

            // Parse outfmt
            let output_format = match format.as_str() {
                "table" => skope::OutputFormat::Table,
                "tsv" => skope::OutputFormat::Tsv,
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
                "o" => skope::SortOrder::Original,
                "t" => skope::SortOrder::Target,
                "s" => skope::SortOrder::Sample,
                "c" => skope::SortOrder::Containment,
                _ => unreachable!("clap should have validated the sort order"),
            };

            let config = skope::ContainmentConfig {
                targets_path: targets.clone(),
                sample_paths: expanded_samples,
                sample_names: derived_sample_names,
                kmer_length: *kmer_length,
                smer_length: *smer_length,
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
                dump_positions_path: dump_positions.as_ref().map(PathBuf::from),
                no_total: *no_total,
            };

            config
                .execute()
                .context("Failed to run containment analysis")?;
        }
        Commands::Lenhist {
            targets,
            samples,
            sample_names,
            kmer_length,
            smer_length,
            threads,
            output,
            quiet,
            limit,
        } => {
            // Expand directories to lists of files
            let (expanded_samples, is_directory) = expand_sample_inputs(samples)?;

            // Derive or validate sample names
            let derived_sample_names: Vec<String> = if let Some(names) = sample_names {
                // User-provided names
                if names.len() != expanded_samples.len() {
                    return Err(anyhow::anyhow!(
                        "Number of sample names ({}) must match number of samples ({})",
                        names.len(),
                        expanded_samples.len()
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

            validate_k_s(*kmer_length, *smer_length)?;

            // Configure thread pool if specified (non-zero)
            if *threads > 0 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(*threads)
                    .build_global()
                    .context("Failed to initialise thread pool")?;
            }

            // Parse limit if provided
            let limit_bp = if let Some(s) = limit {
                Some(parse_sample(s)?)
            } else {
                None
            };

            // Detect if user wants all seqs (no target filtering)
            let include_all_seqs = targets.to_string_lossy() == "-";

            let config = skope::LengthHistogramConfig {
                targets_path: targets.clone(),
                sample_paths: expanded_samples,
                sample_names: derived_sample_names,
                kmer_length: *kmer_length,
                smer_length: *smer_length,
                threads: *threads,
                output_path: if output == "-" {
                    None
                } else {
                    Some(PathBuf::from(output))
                },
                quiet: *quiet,
                limit_bp,
                include_all_seqs,
            };

            config
                .execute()
                .context("Failed to run length histogram analysis")?;
        }
    }

    Ok(())
}

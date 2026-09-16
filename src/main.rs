use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use std::collections::HashSet;
use std::path::PathBuf;

use skope::{
    DEFAULT_KMER_LENGTH, DEFAULT_SMER_LENGTH, derive_sample_name, find_fastx_files,
    find_fastx_files_recursive, is_special_input_path, resolve_k_s, validate_k_s,
};

/// Check the kdust threshold is in [0, 1], rejecting NaN and inf
fn validate_complexity(complexity: f32) -> Result<()> {
    if !(0.0..=1.0).contains(&complexity) {
        return Err(anyhow::anyhow!(
            "Invalid --complexity {complexity}: must be in [0, 1] (0 = retain all)"
        ));
    }
    Ok(())
}

/// Validate the FracMinHash fraction is in (0, 1]
fn validate_fraction(fraction: f64) -> Result<()> {
    if !(fraction > 0.0 && fraction <= 1.0) {
        return Err(anyhow::anyhow!(
            "Invalid --fraction {fraction}: must be in (0, 1] (1 = retain all)"
        ));
    }
    Ok(())
}

fn k_help() -> String {
    format!("K-mer length (1-61, default {DEFAULT_KMER_LENGTH}), read from a prebuilt index")
}

fn s_help() -> String {
    format!("S-mer length (odd, s < k, default {DEFAULT_SMER_LENGTH}), read from a prebuilt index")
}

fn validate_sample_names(names: &[String]) -> Result<()> {
    let mut seen = HashSet::new();
    let mut duplicates = Vec::new();

    for name in names {
        if !seen.insert(name) && !duplicates.contains(name) {
            duplicates.push(name.clone());
        }
    }

    if !duplicates.is_empty() {
        return Err(anyhow::anyhow!(
            "Duplicate sample names detected: {}. Please rename files or use --names to specify unique names.",
            duplicates.join(", ")
        ));
    }

    Ok(())
}

#[derive(Debug)]
struct PreparedSamples {
    paths: Vec<Vec<PathBuf>>,
    names: Vec<String>,
}

fn prepare_samples(inputs: &[PathBuf], names: Option<&[String]>) -> Result<PreparedSamples> {
    if let Some(names) = names
        && names.len() != inputs.len()
    {
        return Err(anyhow::anyhow!(
            "Number of sample names ({}) must match number of samples ({})",
            names.len(),
            inputs.len()
        ));
    }

    let mut paths = Vec::with_capacity(inputs.len());
    let mut prepared_names = Vec::with_capacity(inputs.len());

    for (i, input) in inputs.iter().enumerate() {
        let (sample_paths, is_directory) =
            if input.to_string_lossy() == "-" || is_special_input_path(input) {
                (vec![input.clone()], false)
            } else {
                if !input.exists() {
                    return Err(anyhow::anyhow!("Path does not exist: {}", input.display()));
                }
                let metadata = std::fs::metadata(input)
                    .with_context(|| format!("Failed to access path: {}", input.display()))?;
                if metadata.is_file() {
                    (vec![input.clone()], false)
                } else if metadata.is_dir() {
                    (find_fastx_files(input)?, true)
                } else {
                    return Err(anyhow::anyhow!(
                        "Path is neither a regular file nor directory: {}",
                        input.display()
                    ));
                }
            };
        paths.push(sample_paths);
        prepared_names.push(names.map_or_else(
            || derive_sample_name(input, is_directory),
            |names| names[i].clone(),
        ));
    }

    validate_sample_names(&prepared_names)?;
    Ok(PreparedSamples {
        paths,
        names: prepared_names,
    })
}

fn initialise_thread_pool(threads: usize) -> Result<()> {
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .context("Failed to initialise thread pool")?;
    }
    Ok(())
}

fn output_path(output: &str) -> Option<PathBuf> {
    (output != "-").then(|| PathBuf::from(output))
}

fn parse_limit(limit: Option<&str>) -> Result<Option<u64>> {
    limit.map(parse_bases).transpose()
}

/// Expand background inputs into a flat list of fastx files (directories searched recursively)
fn expand_background_inputs(inputs: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for input in inputs {
        if input.to_string_lossy() == "-" || is_special_input_path(input) {
            files.push(input.clone());
            continue;
        }
        if !input.exists() {
            return Err(anyhow::anyhow!(
                "Background path does not exist: {}",
                input.display()
            ));
        }
        if std::fs::metadata(input)?.is_dir() {
            files.extend(find_fastx_files_recursive(input)?);
        } else {
            files.push(input.clone());
        }
    }
    Ok(files)
}

/// Parse a base-count string with K/M/G/T suffix into a bp count
fn parse_bases(s: &str) -> Result<u64> {
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
        .with_context(|| format!("Invalid base count: {}", s))?;

    Ok(num * multiplier)
}

#[derive(Parser)]
#[command(author, version, about = "Containment and abundance estimation using open syncmers or all k-mers", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum IndexCommands {
    /// Build a classification index (.sk) from fastx file(s) or a directory of groups (alpha)
    #[command(alias = "build")]
    BuildClassify {
        /// Path to fastx file (single group unless -i), directory of fastx files/subdirs (one group per child file/subdir), or - for stdin
        targets: PathBuf,

        /// Treat each fastx record as a separate group (single fastx file only)
        #[arg(short = 'i', long = "individual", default_value_t = false)]
        individual: bool,

        /// K-mer length (1-61, must be odd)
        #[arg(short = 'k', long = "kmer", value_name = "K", default_value_t = DEFAULT_KMER_LENGTH, value_parser = clap::value_parser!(u8).range(1..=61))]
        kmer_length: u8,

        /// S-mer length used for open syncmer selection (s < k, s must be odd)
        #[arg(short = 's', long = "smer", value_name = "S", default_value_t = DEFAULT_SMER_LENGTH)]
        smer_length: u8,

        /// Evaluate all k-mers, bypassing syncmer selection (alias for s=0)
        #[arg(
            long = "all-kmers",
            default_value_t = false,
            conflicts_with = "smer_length"
        )]
        all_kmers: bool,

        /// Discard target k-mers below this kdust complexity [0, 1] (0 = retain all)
        #[arg(long = "complexity", value_name = "FLOAT", default_value_t = 0.0)]
        complexity: f32,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Path to output index file (.sk) (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Suppress progress reporting
        #[arg(short = 'q', long = "quiet", default_value_t = false)]
        quiet: bool,
    },

    /// Build a query index (.sk) from target fastx file(s), optionally masking background k-mers (alpha)
    BuildQuery {
        /// Path to fastx file (single target unless -i), directory of fastx files/subdirs (one target per child file/subdir), or - for stdin
        targets: PathBuf,

        /// Path to fastx file(s) whose k-mers we wish to drop from our targets
        #[arg(short = 'b', long = "background")]
        background: Vec<PathBuf>,

        /// K-mer length (1-61)
        #[arg(short = 'k', long = "kmer", value_name = "K", default_value_t = DEFAULT_KMER_LENGTH, value_parser = clap::value_parser!(u8).range(1..=61))]
        kmer_length: u8,

        /// S-mer length used for open syncmer selection (s < k, s must be odd)
        #[arg(short = 's', long = "smer", value_name = "S", default_value_t = DEFAULT_SMER_LENGTH)]
        smer_length: u8,

        /// Treat each fastx record as a separate target (default: merge records into one target)
        #[arg(short = 'i', long = "individual", default_value_t = false)]
        individual: bool,

        /// Store k-mer positions (needed for --confidence/--dump-kmers at query time)
        #[arg(short = 'p', long = "positions", default_value_t = false)]
        positions: bool,

        /// Evaluate all k-mers, bypassing syncmer selection (alias for s=0)
        #[arg(
            long = "all-kmers",
            default_value_t = false,
            conflicts_with = "smer_length"
        )]
        all_kmers: bool,

        /// Discard target k-mers below this kdust complexity [0, 1] (0 = retain all)
        #[arg(long = "complexity", value_name = "FLOAT", default_value_t = 0.0)]
        complexity: f32,

        /// FracMinHash fraction of target k-mers to keep [0, 1]
        #[arg(
            short = 'f',
            long = "fraction",
            value_name = "FLOAT",
            default_value_t = 1.0
        )]
        fraction: f64,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Path to output index file (.sk) (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Suppress progress reporting
        #[arg(short = 'q', long = "quiet", default_value_t = false)]
        quiet: bool,
    },

    /// Show metadata for an index (.sk) (alpha)
    Info {
        /// Path to index file (.sk)
        index: PathBuf,
    },
}

#[derive(Subcommand)]
enum Commands {
    /// Estimate target containment & abundance in fastx file(s) or directories thereof using open syncmers or all k-mers
    Query {
        /// Path to fastx file (single target unless -i), directory of fastx files/subdirs (one target per child file/subdir) or query index (.sk)
        targets: PathBuf,

        /// Path(s) to fastx files/dirs (- for stdin). Each file/dir is treated as a separate sample
        #[arg(required = true)]
        samples: Vec<PathBuf>,

        #[arg(short = 'k', long = "kmer", value_name = "K", value_parser = clap::value_parser!(u8).range(1..=61), help = k_help())]
        kmer_length: Option<u8>,

        #[arg(short = 's', long = "smer", value_name = "S", help = s_help())]
        smer_length: Option<u8>,

        /// Treat each fastx record as separate target (default: merge records into one target)
        #[arg(short = 'i', long = "individual", default_value_t = false)]
        individual: bool,

        /// Report confidence intervals, ANI estimates, and patchiness columns
        #[arg(short = 'c', long = "confidence", default_value_t = false)]
        confidence: bool,

        /// Consider only k-mers unique to each target
        #[arg(short = 'd', long = "discriminatory", default_value_t = false)]
        discriminatory: bool,

        /// Evaluate all k-mers, bypassing syncmer selection (alias for s=0)
        #[arg(
            long = "all-kmers",
            default_value_t = false,
            conflicts_with = "smer_length"
        )]
        all_kmers: bool,

        /// Discard target k-mers below this kdust complexity [0, 1] (0 = retain all)
        #[arg(long = "complexity", value_name = "FLOAT", default_value_t = 0.0)]
        complexity: f32,

        /// FracMinHash fraction of target k-mers to keep [0, 1]
        #[arg(
            short = 'f',
            long = "fraction",
            value_name = "FLOAT",
            default_value_t = 1.0
        )]
        fraction: f64,

        /// Comma-separated additional abundance thresholds for containment estimation
        #[arg(
            short = 'a',
            long = "abundance-thresholds",
            value_name = "INT,...",
            value_delimiter = ',',
            default_value = "10"
        )]
        abundance_thresholds: Vec<usize>,

        /// Path to fastx file(s) whose k-mers we wish to drop from our targets
        #[arg(short = 'b', long = "background")]
        background: Vec<PathBuf>,

        /// Terminate processing after approximately this many bases (e.g. 50M, 10G)
        #[arg(short = 'l', long = "limit", value_name = "BASES")]
        limit: Option<String>,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Comma-separated sample names (default is file/dir name without extension)
        #[arg(
            short = 'n',
            long = "names",
            value_name = "NAME,...",
            value_delimiter = ','
        )]
        sample_names: Option<Vec<String>>,

        /// Sort results
        #[arg(long = "sort", default_value = "containment", value_parser = ["containment", "target", "input"])]
        sort: String,

        /// Dump selected target k-mers to TSV file (target, position, kmer)
        #[arg(long = "dump-kmers", value_name = "FILE")]
        dump_kmers: Option<PathBuf>,

        /// Suppress TOTAL summary rows in output
        #[arg(long = "no-total", default_value_t = false)]
        no_total: bool,

        /// Suppress progress reporting
        #[arg(short = 'q', long = "quiet", default_value_t = false)]
        quiet: bool,
    },

    /// Classify sequences into groups by k-mer content (alpha)
    Classify {
        /// Path to fastx file (single group unless -i), directory of fastx files/subdirs (one group per child file/subdir) or classification index (.sk)
        targets: PathBuf,

        /// Treat each fastx record as a separate group (single fastx file only)
        #[arg(short = 'i', long = "individual", default_value_t = false)]
        individual: bool,

        /// Path(s) to fastx files/dirs (- for stdin)
        #[arg(required = true)]
        samples: Vec<PathBuf>,

        #[arg(short = 'k', long = "kmer", value_name = "K", value_parser = clap::value_parser!(u8).range(1..=61), help = k_help())]
        kmer_length: Option<u8>,

        #[arg(short = 's', long = "smer", value_name = "S", help = s_help())]
        smer_length: Option<u8>,

        /// Consider only k-mers unique to each group
        #[arg(short = 'd', long = "discriminatory", default_value_t = false)]
        discriminatory: bool,

        /// Evaluate all k-mers, bypassing syncmer selection (alias for s=0)
        #[arg(
            long = "all-kmers",
            default_value_t = false,
            conflicts_with = "smer_length"
        )]
        all_kmers: bool,

        /// Discard target k-mers below this kdust complexity [0, 1] (0 = retain all)
        #[arg(long = "complexity", value_name = "FLOAT", default_value_t = 0.0)]
        complexity: f32,

        /// Minimum absolute number of k-mer hits for a match
        #[arg(
            short = 'a',
            long = "abs-threshold",
            value_name = "ABS_THRESHOLD",
            default_value_t = 1
        )]
        abs_threshold: u64,

        /// Minimum relative proportion (0.0-1.0) of k-mer hits for a match
        #[arg(
            short = 'r',
            long = "rel-threshold",
            value_name = "REL_THRESHOLD",
            default_value_t = 0.0
        )]
        rel_threshold: f64,

        /// Terminate processing after approximately this many bases (e.g. 50M, 10G)
        #[arg(short = 'l', long = "limit", value_name = "BASES")]
        limit: Option<String>,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Comma-separated sample names (default is file/dir name without extension)
        #[arg(
            short = 'n',
            long = "names",
            value_name = "NAME,...",
            value_delimiter = ','
        )]
        sample_names: Option<Vec<String>>,

        /// Output per-sequence classifications instead of summary
        #[arg(long = "per-seq", default_value_t = false)]
        per_seq: bool,

        /// Suppress progress reporting
        #[arg(short = 'q', long = "quiet", default_value_t = false)]
        quiet: bool,
    },

    /// Generate per-group length histograms based on k-mer classification (alpha)
    Lenhist {
        /// Path to fastx file (single group unless -i), directory of fastx files/subdirs (one group per child file/subdir), classification index (.sk), or - to disable group filtering (single "all" bucket)
        targets: PathBuf,

        /// Treat each fastx record as a separate group (single fastx file only)
        #[arg(short = 'i', long = "individual", default_value_t = false)]
        individual: bool,

        /// Path(s) to fastx files/dirs (- for stdin). Each file/dir is treated as a separate sample
        #[arg(required = true)]
        samples: Vec<PathBuf>,

        // Algorithm parameters
        #[arg(short = 'k', long = "kmer", value_name = "K", value_parser = clap::value_parser!(u8).range(1..=61), help = k_help())]
        kmer_length: Option<u8>,

        #[arg(short = 's', long = "smer", value_name = "S", help = s_help())]
        smer_length: Option<u8>,

        /// Consider only k-mers unique to each group
        #[arg(short = 'd', long = "discriminatory", default_value_t = false)]
        discriminatory: bool,

        /// Evaluate all k-mers, bypassing syncmer selection (alias for s=0)
        #[arg(
            long = "all-kmers",
            default_value_t = false,
            conflicts_with = "smer_length"
        )]
        all_kmers: bool,

        /// Discard target k-mers below this kdust complexity [0, 1] (0 = retain all)
        #[arg(long = "complexity", value_name = "FLOAT", default_value_t = 0.0)]
        complexity: f32,

        /// Minimum absolute number of k-mer hits for a match
        #[arg(
            short = 'a',
            long = "abs-threshold",
            value_name = "ABS_THRESHOLD",
            default_value_t = 1
        )]
        abs_threshold: u64,

        /// Minimum relative proportion (0.0-1.0) of k-mer hits for a match
        #[arg(
            short = 'r',
            long = "rel-threshold",
            value_name = "REL_THRESHOLD",
            default_value_t = 0.0
        )]
        rel_threshold: f64,

        // Processing options
        /// Terminate processing after approximately this many bases (e.g. 50M, 10G)
        #[arg(short = 'l', long = "limit", value_name = "BASES")]
        limit: Option<String>,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        // Output options
        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Comma-separated sample names (default is file/dir name without extension)
        #[arg(
            short = 'n',
            long = "names",
            value_name = "NAME,...",
            value_delimiter = ','
        )]
        sample_names: Option<Vec<String>>,

        /// Suppress progress reporting
        #[arg(short = 'q', long = "quiet", default_value_t = false)]
        quiet: bool,
    },

    /// Build and manage query and classification indexes (alpha)
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
            IndexCommands::BuildClassify {
                targets,
                individual,
                kmer_length,
                smer_length,
                all_kmers,
                complexity,
                threads,
                output,
                quiet,
            } => {
                let smer_length = if *all_kmers { 0 } else { *smer_length };
                validate_k_s(*kmer_length, smer_length)?;
                validate_complexity(*complexity)?;

                initialise_thread_pool(*threads)?;

                let config = skope::BuildClassifyConfig {
                    targets_path: targets.clone(),
                    individual: *individual,
                    kmer_length: *kmer_length,
                    smer_length,
                    complexity: *complexity,
                    threads: *threads,
                    output_path: output_path(output),
                    quiet: *quiet,
                };

                skope::run_build_classify(&config)
                    .context("Failed to build classification index")?;
            }

            IndexCommands::BuildQuery {
                targets,
                background,
                kmer_length,
                smer_length,
                all_kmers,
                individual,
                positions,
                fraction,
                complexity,
                threads,
                output,
                quiet,
            } => {
                let smer_length = if *all_kmers { 0 } else { *smer_length };
                validate_k_s(*kmer_length, smer_length)?;
                validate_fraction(*fraction)?;
                validate_complexity(*complexity)?;

                initialise_thread_pool(*threads)?;

                let config = skope::BuildQueryConfig {
                    targets_path: targets.clone(),
                    background_paths: expand_background_inputs(background)?,
                    kmer_length: *kmer_length,
                    smer_length,
                    individual: *individual,
                    positions: *positions,
                    threads: *threads,
                    output_path: output_path(output),
                    quiet: *quiet,
                    fraction: *fraction,
                    complexity: *complexity,
                };

                skope::run_build_query(&config).context("Failed to build query index")?;
            }

            IndexCommands::Info { index } => {
                skope::run_index_info(index).context("Failed to run index info command")?;
            }
        },

        Commands::Classify {
            targets,
            individual,
            samples,
            sample_names,
            kmer_length,
            smer_length,
            all_kmers,
            abs_threshold,
            rel_threshold,
            discriminatory,
            complexity,
            threads,
            limit,
            output,
            per_seq,
            quiet,
        } => {
            let prepared = prepare_samples(samples, sample_names.as_deref())?;
            validate_complexity(*complexity)?;
            let (kmer_length, smer_length) = resolve_k_s(
                targets,
                *kmer_length,
                if *all_kmers { Some(0) } else { *smer_length },
                *quiet,
            )?;
            initialise_thread_pool(*threads)?;

            let config = skope::ClassifyConfig {
                targets_path: targets.clone(),
                individual: *individual,
                sample_paths: prepared.paths,
                sample_names: prepared.names,
                kmer_length,
                smer_length,
                complexity: *complexity,
                abs_threshold: *abs_threshold,
                rel_threshold: *rel_threshold,
                threads: *threads,
                limit_bp: parse_limit(limit.as_deref())?,
                output_path: output_path(output),
                per_seq: *per_seq,
                discriminatory: *discriminatory,
                quiet: *quiet,
            };

            skope::run_classification(&config).context("Failed to run classification")?;
        }

        Commands::Query {
            targets,
            samples,
            sample_names,
            kmer_length,
            smer_length,
            all_kmers,
            threads,
            output,
            quiet,
            abundance_thresholds,
            discriminatory,
            individual,
            limit,
            sort,
            dump_kmers,
            no_total,
            confidence,
            background,
            fraction,
            complexity,
        } => {
            let prepared = prepare_samples(samples, sample_names.as_deref())?;
            let background_paths = expand_background_inputs(background)?;
            validate_fraction(*fraction)?;
            validate_complexity(*complexity)?;
            let (kmer_length, smer_length) = resolve_k_s(
                targets,
                *kmer_length,
                if *all_kmers { Some(0) } else { *smer_length },
                *quiet,
            )?;
            initialise_thread_pool(*threads)?;

            // Parse sort order
            let sort_order = match sort.as_str() {
                "input" => skope::SortOrder::Original,
                "target" => skope::SortOrder::Target,
                "containment" => skope::SortOrder::Containment,
                _ => unreachable!("clap should have validated the sort order"),
            };

            let config = skope::ContainmentConfig {
                targets_path: targets.clone(),
                background_paths,
                sample_paths: prepared.paths,
                sample_names: prepared.names,
                kmer_length,
                smer_length,
                threads: *threads,
                output_path: output_path(output),
                quiet: *quiet,
                abundance_thresholds: abundance_thresholds.clone(),
                discriminatory: *discriminatory,
                individual: *individual,
                limit_bp: parse_limit(limit.as_deref())?,
                sort_order,
                dump_kmers_path: dump_kmers.clone(),
                no_total: *no_total,
                confidence: *confidence,
                fraction: *fraction,
                complexity: *complexity,
            };

            config
                .execute()
                .context("Failed to run containment analysis")?;
        }
        Commands::Lenhist {
            targets,
            individual,
            samples,
            sample_names,
            kmer_length,
            smer_length,
            all_kmers,
            abs_threshold,
            rel_threshold,
            discriminatory,
            complexity,
            threads,
            output,
            quiet,
            limit,
        } => {
            let prepared = prepare_samples(samples, sample_names.as_deref())?;
            validate_complexity(*complexity)?;
            let (kmer_length, smer_length) = resolve_k_s(
                targets,
                *kmer_length,
                if *all_kmers { Some(0) } else { *smer_length },
                *quiet,
            )?;
            initialise_thread_pool(*threads)?;

            // `-` keeps its lenhist-specific meaning: no filtering, single "all" bucket
            let no_filter = targets.to_string_lossy() == "-";

            let config = skope::LengthHistogramConfig {
                targets_path: targets.clone(),
                individual: *individual,
                sample_paths: prepared.paths,
                sample_names: prepared.names,
                kmer_length,
                smer_length,
                complexity: *complexity,
                abs_threshold: *abs_threshold,
                rel_threshold: *rel_threshold,
                discriminatory: *discriminatory,
                threads: *threads,
                output_path: output_path(output),
                quiet: *quiet,
                limit_bp: parse_limit(limit.as_deref())?,
                no_filter,
            };

            config
                .execute()
                .context("Failed to run length histogram analysis")?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn prepare_samples_expands_inputs_and_derives_names() {
        let temp = TempDir::new().unwrap();
        let file = temp.path().join("single.fastq");
        let directory = temp.path().join("reads");
        std::fs::write(&file, b"@r\nACGT\n+\nIIII\n").unwrap();
        std::fs::create_dir(&directory).unwrap();
        std::fs::write(directory.join("part.fa"), b">r\nACGT\n").unwrap();

        let prepared = prepare_samples(&[file.clone(), directory.clone()], None).unwrap();
        assert_eq!(prepared.names, ["single", "reads"]);
        assert_eq!(
            prepared.paths,
            [vec![file], vec![directory.join("part.fa")]]
        );
    }

    #[test]
    fn prepare_samples_validates_supplied_names() {
        let error = prepare_samples(&[PathBuf::from("-")], Some(&[])).unwrap_err();
        assert!(error.to_string().contains("must match number of samples"));

        let names = ["same".to_string(), "same".to_string()];
        let error =
            prepare_samples(&[PathBuf::from("-"), PathBuf::from("-")], Some(&names)).unwrap_err();
        assert!(error.to_string().contains("Duplicate sample names"));
    }
}

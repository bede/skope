use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand};
use std::collections::HashSet;
use std::path::{Path, PathBuf};

use skope::{
    DEFAULT_KMER_LENGTH, DEFAULT_SMER_LENGTH, Input, derive_mates_name, derive_sample_name,
    find_seq_files, find_seq_files_recursive, is_special_input_path, resolve_k_s, validate_k_s,
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

fn validate_rel_threshold(rel_threshold: f64) -> Result<()> {
    if !(0.0..=1.0).contains(&rel_threshold) {
        return Err(anyhow::anyhow!(
            "Invalid --rel-threshold {rel_threshold}: must be in [0, 1]"
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

/// Sample arguments shared by query, classify and lenhist
#[derive(Args)]
struct SampleArgs {
    /// Samples as fastx/CBQ files/dirs (- for stdin), or comma-separated mate files (R1,R2)
    #[arg(required = true)]
    samples: Vec<PathBuf>,

    /// Comma-separated sample names (default is file/dir name without extension)
    #[arg(
        short = 'n',
        long = "names",
        value_name = "NAME,...",
        value_delimiter = ','
    )]
    names: Option<Vec<String>>,

    /// Sample files outside R1,R2 pairs are interleaved pairs
    #[arg(long = "interleaved", default_value_t = false)]
    interleaved: bool,

    /// Terminate processing after approximately this many bases (e.g. 50M, 10G)
    #[arg(short = 'l', long = "limit", value_name = "BASES", value_parser = parse_bases)]
    limit: Option<u64>,
}

impl SampleArgs {
    fn prepare(&self) -> Result<PreparedSamples> {
        prepare_samples(&self.samples, self.names.as_deref(), self.interleaved)
    }
}

#[derive(Debug)]
struct PreparedSamples {
    inputs: Vec<Vec<Input>>,
    names: Vec<String>,
}

/// Split an `R1,R2` argument into its two mate files
fn split_mates(arg: &Path) -> Result<(PathBuf, PathBuf)> {
    let arg = arg.to_string_lossy();
    let Some((r1, r2)) = arg
        .split_once(',')
        .filter(|(r1, r2)| !r1.is_empty() && !r2.is_empty() && !r2.contains(','))
    else {
        anyhow::bail!("Mate pairs must be two comma-separated files (R1,R2), got {arg}");
    };
    for mate in [r1, r2].map(Path::new) {
        anyhow::ensure!(
            mate.is_file() || is_special_input_path(mate),
            "Mate is not a file: {}",
            mate.display()
        );
    }
    Ok((r1.into(), r2.into()))
}

/// Resolve each sample argument: an existing path wins, then `R1,R2` mates.
/// `--interleaved` applies to every file outside a mate pair
fn prepare_samples(
    args: &[PathBuf],
    names: Option<&[String]>,
    interleaved: bool,
) -> Result<PreparedSamples> {
    if let Some(names) = names
        && names.len() != args.len()
    {
        return Err(anyhow::anyhow!(
            "Number of sample names ({}) must match number of samples ({})",
            names.len(),
            args.len()
        ));
    }

    let file = |path| match interleaved {
        true => Input::Interleaved(path),
        false => Input::File(path),
    };
    let mut inputs = Vec::with_capacity(args.len());
    let mut prepared_names = Vec::with_capacity(args.len());

    for (i, arg) in args.iter().enumerate() {
        let (sample, name) =
            if arg.as_os_str() == "-" || is_special_input_path(arg) || arg.is_file() {
                (vec![file(arg.clone())], derive_sample_name(arg, false))
            } else if arg.is_dir() {
                let files = find_seq_files(arg)?.into_iter().map(file).collect();
                (files, derive_sample_name(arg, true))
            } else if arg.exists() {
                anyhow::bail!(
                    "Path is neither a regular file nor directory: {}",
                    arg.display()
                );
            } else if arg.to_string_lossy().contains(',') {
                let (r1, r2) = split_mates(arg)?;
                let name = derive_mates_name(&r1);
                (vec![Input::Paired(r1, r2)], name)
            } else {
                anyhow::bail!("Path does not exist: {}", arg.display());
            };
        inputs.push(sample);
        prepared_names.push(names.map_or(name, |names| names[i].clone()));
    }

    validate_sample_names(&prepared_names)?;
    Ok(PreparedSamples {
        inputs,
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

/// Expand background paths into a flat list of sequence files (directories searched recursively)
fn expand_background_paths(inputs: &[PathBuf]) -> Result<Vec<PathBuf>> {
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
            files.extend(find_seq_files_recursive(input)?);
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
    /// Build a classification index (.sk) from fastx/CBQ file(s) or a directory of groups (alpha)
    #[command(alias = "build")]
    BuildClassify {
        /// Path to fastx/CBQ file (single group unless -i), directory of such files/subdirs (one group per child file/subdir), or - for stdin
        targets: PathBuf,

        /// Treat each record as a separate group (single file only)
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

    /// Build a query index (.sk) from target fastx/CBQ file(s), optionally masking background k-mers (alpha)
    BuildQuery {
        /// Path to fastx/CBQ file (single target unless -i), directory of such files/subdirs (one target per child file/subdir), or - for stdin
        targets: PathBuf,

        /// Path to fastx/CBQ file(s) whose k-mers we wish to drop from our targets
        #[arg(short = 'b', long = "background")]
        background: Vec<PathBuf>,

        /// K-mer length (1-61)
        #[arg(short = 'k', long = "kmer", value_name = "K", default_value_t = DEFAULT_KMER_LENGTH, value_parser = clap::value_parser!(u8).range(1..=61))]
        kmer_length: u8,

        /// S-mer length used for open syncmer selection (s < k, s must be odd)
        #[arg(short = 's', long = "smer", value_name = "S", default_value_t = DEFAULT_SMER_LENGTH)]
        smer_length: u8,

        /// Treat each record as a separate target (default: merge records into one target)
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
    /// Estimate target containment & abundance in sequence collections using open syncmers or all k-mers
    Query {
        /// Path to fastx/CBQ file (single target unless -i), directory of such files/subdirs (one target per child file/subdir) or query index (.sk)
        targets: PathBuf,

        #[arg(short = 'k', long = "kmer", value_name = "K", value_parser = clap::value_parser!(u8).range(1..=61), help = k_help())]
        kmer_length: Option<u8>,

        #[arg(short = 's', long = "smer", value_name = "S", help = s_help())]
        smer_length: Option<u8>,

        /// Treat each record as a separate target (default: merge records into one target)
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

        /// Path to fastx/CBQ file(s) whose k-mers we wish to drop from our targets
        #[arg(short = 'b', long = "background")]
        background: Vec<PathBuf>,

        #[command(flatten)]
        samples: SampleArgs,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

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
        /// Path to fastx/CBQ file (single group unless -i), directory of such files/subdirs (one group per child file/subdir) or classification index (.sk)
        targets: PathBuf,

        /// Treat each record as a separate group (single file only)
        #[arg(short = 'i', long = "individual", default_value_t = false)]
        individual: bool,

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

        /// Minimum distinct k-mer hits for a match
        #[arg(
            short = 'a',
            long = "abs-threshold",
            value_name = "ABS_THRESHOLD",
            default_value_t = 1
        )]
        abs_threshold: u64,

        /// Minimum proportion [0, 1] of distinct k-mers hit for a match
        #[arg(
            short = 'r',
            long = "rel-threshold",
            value_name = "REL_THRESHOLD",
            default_value_t = 0.0
        )]
        rel_threshold: f64,

        #[command(flatten)]
        samples: SampleArgs,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Output per-sequence classifications instead of summary
        #[arg(long = "per-seq", default_value_t = false)]
        per_seq: bool,

        /// Suppress progress reporting
        #[arg(short = 'q', long = "quiet", default_value_t = false)]
        quiet: bool,
    },

    /// Generate per-group length histograms based on k-mer classification (alpha)
    Lenhist {
        /// Path to fastx/CBQ file (single group unless -i), directory of such files/subdirs (one group per child file/subdir), classification index (.sk), or - to disable group filtering (single "all" bucket)
        targets: PathBuf,

        /// Treat each record as a separate group (single file only)
        #[arg(short = 'i', long = "individual", default_value_t = false)]
        individual: bool,

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

        /// Minimum distinct k-mer hits for a match
        #[arg(
            short = 'a',
            long = "abs-threshold",
            value_name = "ABS_THRESHOLD",
            default_value_t = 1
        )]
        abs_threshold: u64,

        /// Minimum proportion [0, 1] of distinct k-mers hit for a match
        #[arg(
            short = 'r',
            long = "rel-threshold",
            value_name = "REL_THRESHOLD",
            default_value_t = 0.0
        )]
        rel_threshold: f64,

        // Processing options
        #[command(flatten)]
        samples: SampleArgs,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        // Output options
        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

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
                    background_paths: expand_background_paths(background)?,
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
            kmer_length,
            smer_length,
            all_kmers,
            discriminatory,
            complexity,
            abs_threshold,
            rel_threshold,
            threads,
            output,
            per_seq,
            quiet,
        } => {
            let prepared = samples.prepare()?;
            validate_complexity(*complexity)?;
            validate_rel_threshold(*rel_threshold)?;
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
                sample_inputs: prepared.inputs,
                sample_names: prepared.names,
                kmer_length,
                smer_length,
                complexity: *complexity,
                abs_threshold: *abs_threshold,
                rel_threshold: *rel_threshold,
                threads: *threads,
                limit_bp: samples.limit,
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
            kmer_length,
            smer_length,
            all_kmers,
            threads,
            output,
            quiet,
            abundance_thresholds,
            discriminatory,
            individual,
            sort,
            dump_kmers,
            no_total,
            confidence,
            background,
            fraction,
            complexity,
        } => {
            let prepared = samples.prepare()?;
            let background_paths = expand_background_paths(background)?;
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
                sample_inputs: prepared.inputs,
                sample_names: prepared.names,
                kmer_length,
                smer_length,
                threads: *threads,
                output_path: output_path(output),
                quiet: *quiet,
                abundance_thresholds: abundance_thresholds.clone(),
                discriminatory: *discriminatory,
                individual: *individual,
                limit_bp: samples.limit,
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
            kmer_length,
            smer_length,
            all_kmers,
            discriminatory,
            complexity,
            abs_threshold,
            rel_threshold,
            threads,
            output,
            quiet,
        } => {
            let prepared = samples.prepare()?;
            validate_complexity(*complexity)?;
            validate_rel_threshold(*rel_threshold)?;
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
                sample_inputs: prepared.inputs,
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
                limit_bp: samples.limit,
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
        let directory = temp.path().join("collection");
        std::fs::write(&file, b"@r\nACGT\n+\nIIII\n").unwrap();
        std::fs::create_dir(&directory).unwrap();
        std::fs::write(directory.join("part.fa"), b">r\nACGT\n").unwrap();

        let args = [file.clone(), directory.clone()];
        let prepared = prepare_samples(&args, None, false).unwrap();
        assert_eq!(prepared.names, ["single", "collection"]);
        assert_eq!(
            prepared.inputs,
            [
                vec![Input::File(file.clone())],
                vec![Input::File(directory.join("part.fa"))]
            ]
        );
        let prepared = prepare_samples(&args[..1], None, true).unwrap();
        assert_eq!(prepared.inputs, [vec![Input::Interleaved(file)]]);
    }

    #[test]
    fn prepare_samples_splits_comma_mates() {
        let temp = TempDir::new().unwrap();
        let (r1, r2) = (temp.path().join("s_R1.fq"), temp.path().join("s_R2.fq"));
        std::fs::write(&r1, b"@r\nACGT\n+\nIIII\n").unwrap();
        std::fs::write(&r2, b"@r\nACGT\n+\nIIII\n").unwrap();
        let comma = |a: &Path, b: &Path| PathBuf::from(format!("{},{}", a.display(), b.display()));

        // Mates stay paired under --interleaved and mix with single files
        let args = [comma(&r1, &r2), r1.clone()];
        let prepared = prepare_samples(&args, None, true).unwrap();
        assert_eq!(prepared.names, ["s", "s_R1"]);
        assert_eq!(
            prepared.inputs,
            [
                vec![Input::Paired(r1.clone(), r2.clone())],
                vec![Input::Interleaved(r1.clone())]
            ]
        );

        // A literal path wins over comma splitting
        let literal = comma(&r1, Path::new("x"));
        std::fs::write(&literal, b">r\nACGT\n").unwrap();
        let prepared = prepare_samples(std::slice::from_ref(&literal), None, false).unwrap();
        assert_eq!(prepared.inputs, [vec![Input::File(literal)]]);

        let missing = comma(&r1, &temp.path().join("y"));
        let error = prepare_samples(&[missing], None, false).unwrap_err();
        assert!(error.to_string().contains("Mate is not a file"));
        for bad in [",", "a,", ",b", "a,b,c"].map(PathBuf::from) {
            let error = prepare_samples(&[bad], None, false).unwrap_err();
            assert!(error.to_string().contains("R1,R2"), "{error}");
        }
    }

    #[test]
    fn prepare_samples_validates_supplied_names() {
        let error = prepare_samples(&[PathBuf::from("-")], Some(&[]), false).unwrap_err();
        assert!(error.to_string().contains("must match number of samples"));

        let names = ["same".to_string(), "same".to_string()];
        let error = prepare_samples(
            &[PathBuf::from("-"), PathBuf::from("-")],
            Some(&names),
            false,
        )
        .unwrap_err();
        assert!(error.to_string().contains("Duplicate sample names"));
    }
}

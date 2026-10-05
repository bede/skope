//! Native query commands over paths, multithreaded
use super::*;
use crate::io::{
    Progress, StdinTargets, TargetSource, index_writer, output_writer, process_input,
    resolve_targets, validate_index_output,
};
use crate::kmers::{decode_u64, decode_u128};
use crate::{complexity_info_line, format_bp_per_sec};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::time::Instant;

pub struct ContainmentConfig {
    pub targets_path: PathBuf,
    pub background_paths: Vec<PathBuf>, // Off-target sequences to mask (empty = none)
    pub sample_inputs: Vec<Vec<Input>>,
    pub sample_names: Vec<String>,
    pub kmer_length: u8,
    pub smer_length: u8,
    pub threads: usize,
    pub output_path: Option<PathBuf>,
    pub quiet: bool,
    pub abundance_thresholds: Vec<usize>,
    pub discriminatory: bool,
    pub individual: bool,
    pub limit_bp: Option<u64>,
    pub sort_order: SortOrder,
    pub dump_kmers_path: Option<PathBuf>,
    pub no_total: bool,
    pub confidence: bool,
    pub fraction: f64,
    pub complexity: f32,
}

impl ContainmentConfig {
    pub fn execute(&self) -> Result<()> {
        run_query(self)
    }
}

/// Targets of one file, one per record or mate pair
fn process_targets_file(
    targets_path: &Path,
    kmer_length: u8,
    smer_length: u8,
    fmh: FracMinHash,
    kdust: Kdust,
    quiet: bool,
    collect_positions: bool,
) -> Result<Vec<TargetInfo>> {
    let processor = TargetsProcessor::new(kmer_length, smer_length, fmh, kdust, collect_positions);
    let targets = Arc::clone(&processor.targets);
    // Single thread to preserve order
    process_input(
        &Input::File(targets_path.to_path_buf()),
        processor,
        Progress::new("Collecting target k-mers", quiet, None)?,
        1,
    )?;
    Ok(Arc::into_inner(targets).unwrap().into_inner())
}

/// Process resolved targets into `TargetInfo`s
///
/// A directory always yields one target per child, named after that child. A single sequence
/// file yields one target named after the file, unless `--individual` is set or it holds a
/// single record, in which case records keep their own names.
fn process_target_groups(
    source: &TargetSource,
    individual: bool,
    kmer_length: u8,
    smer_length: u8,
    fmh: FracMinHash,
    kdust: Kdust,
    quiet: bool,
    collect_positions: bool,
) -> Result<Vec<TargetInfo>> {
    let groups = source.groups();
    let per_record = matches!(source, TargetSource::File(_));

    let mut results: Vec<TargetInfo> = Vec::with_capacity(groups.len());
    for group in groups {
        let mut all_targets: Vec<TargetInfo> = Vec::new();
        for file_path in &group.files {
            let targets = process_targets_file(
                file_path,
                kmer_length,
                smer_length,
                fmh,
                kdust,
                quiet,
                collect_positions,
            )?;
            all_targets.extend(targets);
        }
        results.extend(group_targets(
            all_targets,
            &group.name,
            per_record,
            individual,
        )?);
    }

    Ok(results)
}

#[allow(clippy::too_many_arguments)]
fn process_seqs_input(
    input: &Input,
    targets_kmers: Arc<KmerSet>,
    kmer_length: u8,
    smer_length: u8,
    threads: usize,
    quiet: bool,
    limit_bp: Option<u64>,
    label: &'static str,
    summary_label: &'static str,
) -> Result<(AbundanceMap, u64, u64)> {
    let start_time = Instant::now();
    let total_target_kmers = targets_kmers.len();
    let global_counts = Arc::new(Mutex::new(AbundanceMap::new(kmer_length)));
    let processor = CountProcessor::new(
        kmer_length,
        smer_length,
        targets_kmers,
        Arc::clone(&global_counts),
    );
    let progress = Progress::new(label, quiet, limit_bp)?;
    let stats = process_input(input, processor, progress, threads)?;
    let abundance_map = Arc::into_inner(global_counts).unwrap().into_inner();

    if !quiet {
        let bp_per_sec = stats.total_bp as f64 / start_time.elapsed().as_secs_f64();
        eprintln!(
            "{} ({})",
            counted_summary(
                summary_label,
                &stats,
                abundance_map.len(),
                total_target_kmers
            ),
            format_bp_per_sec(bp_per_sec)
        );
    }

    Ok((abundance_map, stats.total_seqs, stats.total_bp))
}

/// Process a single sample's sequences and calculate statistics
fn process_sample(
    sample: &[Input],
    sample_name: &str,
    targets: &[TargetInfo],
    targets_kmers: Arc<KmerSet>,
    abundance_thresholds: &[usize],
    config: &ContainmentConfig,
) -> Result<SampleResults> {
    // Silence per-sample progress for >1 sample
    let quiet_sample = config.quiet || config.sample_inputs.len() > 1;

    let mut combined_abundance_map = AbundanceMap::new(config.kmer_length);
    let mut total_seqs = 0u64;
    let mut total_bp = 0u64;

    // Process each file and accumulate results
    for input in sample {
        let (file_abundance_map, file_seqs, file_bp) = process_seqs_input(
            input,
            Arc::clone(&targets_kmers),
            config.kmer_length,
            config.smer_length,
            config.threads,
            quiet_sample,
            config.limit_bp.map(|limit| limit.saturating_sub(total_bp)),
            "Processing sample",
            "Sample",
        )?;

        // Merge abundance maps
        match (&mut combined_abundance_map, file_abundance_map) {
            (AbundanceMap::U64(combined), AbundanceMap::U64(new)) => {
                merge_counts(new, combined);
            }
            (AbundanceMap::U128(combined), AbundanceMap::U128(new)) => {
                merge_counts(new, combined);
            }
            _ => panic!("Mismatch between AbundanceMap types"),
        }

        total_seqs += file_seqs;
        total_bp += file_bp;

        // Check if limit reached
        if let Some(limit) = config.limit_bp
            && total_bp >= limit
        {
            break;
        }
    }

    let abundance_map = combined_abundance_map;

    // Calculate containment statistics for this sample
    let mut containment_results = calculate_containment_statistics(
        targets,
        &abundance_map,
        abundance_thresholds,
        config.kmer_length,
        config.confidence,
    );
    sort_results(&mut containment_results, config.sort_order);

    Ok(SampleResults::new(
        sample_name.to_string(),
        containment_results,
        abundance_thresholds,
        total_seqs,
        total_bp,
    ))
}

/// Stream background masking and return the number of k-mers removed
fn mask_background(
    targets: &mut [TargetInfo],
    union: &KmerSet,
    background_paths: &[PathBuf],
    kmer_length: u8,
    smer_length: u8,
    threads: usize,
    quiet: bool,
) -> Result<usize> {
    let union = Arc::new(union.clone());
    let mut removed = KmerSet::new(kmer_length);
    for path in background_paths {
        let (map, _, _) = process_seqs_input(
            &Input::File(path.clone()),
            Arc::clone(&union),
            kmer_length,
            smer_length,
            threads,
            quiet,
            None,
            "Masking background",
            "Background",
        )?;
        removed.extend(&map.into_kmers());
    }
    mask_targets(targets, &removed);
    Ok(removed.len())
}

/// Config for `skope index build-query`
pub struct BuildQueryConfig {
    pub targets_path: PathBuf,
    pub background_paths: Vec<PathBuf>,
    pub kmer_length: u8,
    pub smer_length: u8,
    pub individual: bool,
    pub positions: bool,
    pub threads: usize,
    pub output_path: Option<PathBuf>,
    pub quiet: bool,
    /// Fraction for stable FracMinHash selection in `(0, 1]`
    pub fraction: f64,
    /// Minimum kdust in `[0, 1]`, 0 retains all
    pub complexity: f32,
}

/// Read a query index header without loading entries
pub fn read_query_index_meta(path: &Path) -> Result<(u8, u8, f64, f32)> {
    let mut buf = [0u8; 64];
    let n = File::open(path)?.read(&mut buf)?;
    let mut cursor = wincode::io::Cursor::new(&buf[..n]);
    let (magic, kind, version, k, s, _flags, _n, fraction, complexity): QueryIndexHeader =
        wincode::deserialize_from(&mut cursor).context("Failed to read query index header")?;
    validate_query_index_header(&magic, kind, version)?;
    Ok((k, s, fraction, complexity))
}

/// Print human-readable metadata for a query index (`skope index info`)
pub fn print_query_index_info(path: &Path) -> Result<()> {
    let mut reader = BufReader::new(File::open(path)?);
    let (magic, kind, version, kmer_length, smer_length, flags, n_targets, fraction, complexity): QueryIndexHeader =
        wincode::deserialize_from(&mut reader).context("Failed to read query index header")?;
    validate_query_index_header(&magic, kind, version)?;
    let meta: QueryIndexMeta =
        wincode::deserialize_from(&mut reader).context("Failed to decode query index metadata")?;
    if meta.len() != n_targets as usize {
        return Err(anyhow::anyhow!("Query index target count mismatch"));
    }

    let total_kmers: u64 = meta.iter().map(|(_, _, count)| count).sum();
    let total_bp: u64 = meta.iter().map(|(_, length, _)| length).sum();

    println!("Index information:");
    println!("  Format: query");
    println!("  Format version: {version}");
    println!("  K-mer length (k): {kmer_length}");
    println!("  S-mer length (s): {smer_length}");
    println!("  Targets: {n_targets}");
    println!("  Total k-mers: {total_kmers}");
    println!("  Total length: {}", format_bp(total_bp as usize));
    println!(
        "  FracMinHash fraction: {}",
        if fraction >= 1.0 {
            "1 (retain all)".to_string()
        } else {
            format!("{fraction} (keeps ~{:.0}% of k-mers)", fraction * 100.0)
        }
    );
    println!("{}", complexity_info_line(complexity));
    println!(
        "  K-mer positions: {}",
        if flags & QUERY_INDEX_FLAG_POSITIONS != 0 {
            "stored"
        } else {
            "not stored"
        }
    );
    println!(
        "  Background masking: {}",
        if flags & QUERY_INDEX_FLAG_BACKGROUND != 0 {
            "applied"
        } else {
            "none"
        }
    );
    Ok(())
}

fn save_query_index(
    targets: &[TargetInfo],
    kmer_length: u8,
    smer_length: u8,
    fraction: f64,
    complexity: f32,
    flags: u8,
    output_path: Option<&Path>,
) -> Result<()> {
    write_query_index(
        &mut index_writer(output_path)?,
        targets,
        kmer_length,
        smer_length,
        fraction,
        complexity,
        flags,
    )
}

fn load_query_index(path: &Path) -> Result<QueryIndex> {
    let bytes = fs::read(path)
        .with_context(|| format!("Failed to open query index: {}", path.display()))?;
    read_query_index(&bytes)
}

/// Build and serialize a query index (with optional background masking)
pub fn run_build_query(config: &BuildQueryConfig) -> Result<()> {
    validate_index_output(config.output_path.as_deref())?;
    let start = Instant::now();
    let version = env!("CARGO_PKG_VERSION");
    let mut options = String::new();
    if config.fraction < 1.0 {
        options.push_str(&format!(", fraction={}", config.fraction));
    }
    if config.complexity > 0.0 {
        options.push_str(&format!(", complexity={}", config.complexity));
    }
    eprintln!(
        "Skope v{version}; mode: index build-query; options: k={}, s={}, threads={}{}",
        config.kmer_length, config.smer_length, config.threads, options
    );

    let fmh = FracMinHash::from_fraction(config.fraction);
    let kdust = Kdust::from_threshold(config.complexity, config.kmer_length);
    let collect_positions = config.positions;
    let source = resolve_targets(&config.targets_path, IndexKind::Query, StdinTargets::Accept)?;
    if let TargetSource::Index(path) = &source {
        return Err(anyhow::anyhow!(
            "{} is already a skope query index",
            path.display()
        ));
    }
    let mut targets = process_target_groups(
        &source,
        config.individual,
        config.kmer_length,
        config.smer_length,
        fmh,
        kdust,
        config.quiet,
        collect_positions,
    )?;

    let mut flags = 0u8;
    if config.positions {
        flags |= QUERY_INDEX_FLAG_POSITIONS;
    }

    if !config.background_paths.is_empty() {
        let union = build_union(&targets, config.kmer_length);
        let removed = mask_background(
            &mut targets,
            &union,
            &config.background_paths,
            config.kmer_length,
            config.smer_length,
            config.threads,
            config.quiet,
        )?;
        flags |= QUERY_INDEX_FLAG_BACKGROUND;
        // Per-file "found N of M" lines already report single-file totals;
        // only the multi-file aggregate adds information.
        if !config.quiet && config.background_paths.len() > 1 {
            eprintln!(
                "Masked {removed} background k-mers across {} files",
                config.background_paths.len()
            );
        }
    }

    save_query_index(
        &targets,
        config.kmer_length,
        config.smer_length,
        config.fraction,
        config.complexity,
        flags,
        config.output_path.as_deref(),
    )?;

    if !config.quiet {
        let total: usize = targets.iter().map(|t| t.kmers.len()).sum();
        eprintln!(
            "Wrote query index ({} targets, {total} k-mers, positions={}) in {:.1}s",
            targets.len(),
            config.positions,
            start.elapsed().as_secs_f64()
        );
    }
    Ok(())
}

pub fn run_query(config: &ContainmentConfig) -> Result<()> {
    let version = env!("CARGO_PKG_VERSION").to_string();
    let abundance_thresholds = normalize_abundance_thresholds(&config.abundance_thresholds);

    let mut options = format!(
        "k={}, s={}, threads={}",
        config.kmer_length, config.smer_length, config.threads
    );

    if config.sample_inputs.len() > 1 {
        options.push_str(&format!(", samples={}", config.sample_inputs.len()));
    }

    if !abundance_thresholds.is_empty() {
        options.push_str(&format!(
            ", abundance-thresholds={}",
            abundance_thresholds
                .iter()
                .map(|t| t.to_string())
                .collect::<Vec<_>>()
                .join(",")
        ));
    }

    if config.discriminatory {
        options.push_str(", discriminatory");
    }

    if let Some(limit) = config.limit_bp {
        options.push_str(&format!(", limit={}", format_bp(limit as usize)));
    }

    let source = resolve_targets(&config.targets_path, IndexKind::Query, StdinTargets::Reject)?;
    let from_index = matches!(source, TargetSource::Index(_));
    // Stored fraction and complexity trump
    // A CLI-parsed float round-trips through the header bit-identically, so compare exactly
    let (effective_fraction, effective_complexity) = if from_index {
        let (_k, _s, index_fraction, index_complexity) =
            read_query_index_meta(&config.targets_path)?;
        check_index_fraction(config.fraction, index_fraction)?;
        check_index_complexity(config.complexity, index_complexity)?;
        (index_fraction, index_complexity)
    } else {
        (config.fraction, config.complexity)
    };
    if effective_fraction < 1.0 {
        options.push_str(&format!(", fraction={effective_fraction}"));
    }
    if effective_complexity > 0.0 {
        options.push_str(&format!(", complexity={effective_complexity}"));
    }
    eprintln!(
        "Skope v{}; mode: query{}; options: {}",
        version,
        match source {
            TargetSource::Index(_) => " (from index)",
            TargetSource::Directory(_) => " (from directory)",
            TargetSource::File(_) => "",
        },
        options
    );

    // Load a prebuilt query index else extract targets from sequence files
    let need_positions = config.dump_kmers_path.is_some() || config.confidence;
    let mut targets = if let TargetSource::Index(path) = &source {
        let index = load_query_index(path)?;
        index.require_positions(need_positions)?;
        index.targets
    } else {
        process_target_groups(
            &source,
            config.individual,
            config.kmer_length,
            config.smer_length,
            FracMinHash::from_fraction(effective_fraction),
            Kdust::from_threshold(effective_complexity, config.kmer_length),
            config.quiet,
            need_positions,
        )?
    };

    // Runtime masking
    if !config.background_paths.is_empty() {
        let union = build_union(&targets, config.kmer_length);
        let removed = mask_background(
            &mut targets,
            &union,
            &config.background_paths,
            config.kmer_length,
            config.smer_length,
            config.threads,
            config.quiet,
        )?;
        // Per-file "found N of M" lines already report single-file totals;
        // only the multi-file aggregate adds information.
        if !config.quiet && config.background_paths.len() > 1 {
            eprintln!(
                "Masked {removed} background k-mers across {} files",
                config.background_paths.len()
            );
        }
    }

    // Count k-mers shared between targets (only meaningful with >1 target)
    if !config.quiet && targets.len() > 1 {
        eprint!("Counting shared k-mers…\r");
    }
    let (shared_kmers, unique_across_all) =
        share_kmers(&mut targets, config.kmer_length, config.discriminatory);

    if !config.quiet {
        eprint!("\r"); // Clear space
        eprintln!(
            "{}",
            targets_summary(
                &targets,
                config.discriminatory,
                shared_kmers,
                unique_across_all
            )
        );
        if let Some(warning) = empty_targets_warning(&targets) {
            eprintln!("{warning}");
        }
    }

    // Build set of all unique k-mers across targets
    if !config.quiet {
        eprint!("Building k-mer set…\r");
    }
    let targets_kmers = Arc::new(build_union(&targets, config.kmer_length));

    // Dump k-mers (position + k-mer sequence) if requested
    if let Some(ref path) = config.dump_kmers_path {
        let mut file = BufWriter::new(File::create(path)?);
        let k = config.kmer_length;
        let mut total = 0usize;
        for target in &targets {
            match &target.positioned_kmers {
                PositionedKmers::U64(entries) => {
                    for &(value, pos) in entries {
                        let kmer = decode_u64(value, k);
                        writeln!(
                            file,
                            "{}\t{}\t{}",
                            target.name,
                            pos,
                            String::from_utf8_lossy(&kmer)
                        )?;
                        total += 1;
                    }
                }
                PositionedKmers::U128(entries) => {
                    for &(value, pos) in entries {
                        let kmer = decode_u128(value, k);
                        writeln!(
                            file,
                            "{}\t{}\t{}",
                            target.name,
                            pos,
                            String::from_utf8_lossy(&kmer)
                        )?;
                        total += 1;
                    }
                }
                PositionedKmers::Empty => {}
            }
        }
        file.flush()?;
        if !config.quiet {
            eprintln!("Dumped {} k-mers to {}", total, path.display());
        }
    }

    // Process each sample in parallel
    use rayon::prelude::*;

    let is_multisample = config.sample_inputs.len() > 1;
    let completed = if is_multisample && !config.quiet {
        // Give us a blank line to overwrite
        eprint!(
            "\x1B[2K\rSamples: processed 0 of {}…",
            config.sample_inputs.len()
        );
        Some(Arc::new(Mutex::new(0usize)))
    } else {
        // Give us a blank line to overwrite
        if !config.quiet {
            eprint!("\x1B[2K\r");
        }
        None
    };

    let sample_results: Vec<SampleResults> = config
        .sample_inputs
        .par_iter()
        .zip(&config.sample_names)
        .map(|(sample, sample_name)| {
            let result = process_sample(
                sample,
                sample_name,
                &targets,
                Arc::clone(&targets_kmers),
                &abundance_thresholds,
                config,
            );

            // Report completion for multisample runs
            if let Some(ref counter) = completed {
                let mut count = counter.lock();
                *count += 1;
                eprint!(
                    "\rSamples: processed {} of {}…",
                    *count,
                    config.sample_inputs.len()
                );
            }

            result
        })
        .collect::<Result<Vec<_>>>()?;

    if is_multisample && !config.quiet {
        eprintln!();
    }

    // Output results
    output_results(
        &sample_results,
        config.output_path.as_deref(),
        &abundance_thresholds,
        config.confidence,
        config.no_total,
    )?;

    Ok(())
}

fn output_results(
    samples: &[SampleResults],
    output_path: Option<&Path>,
    abundance_thresholds: &[usize],
    confidence: bool,
    no_total: bool,
) -> Result<()> {
    let mut writer = output_writer(output_path)?;
    output_tsv(
        &mut writer,
        samples,
        abundance_thresholds,
        confidence,
        no_total,
    )?;

    Ok(())
}

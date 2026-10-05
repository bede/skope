//! `QuerySession` over streams matches `skope query` and `skope index build-query` over files
use skope::{
    BuildQueryConfig, ContainmentConfig, Input, QueryOptions, QuerySession, SortOrder,
    derive_sample_name,
};
use std::fs::File;
use std::path::{Path, PathBuf};
use tempfile::{NamedTempFile, TempDir};

const TARGETS: &str = "data/zmrp21.viruses.fa";
const TARGET_DIR: &str = "data/target-dir";
const SAMPLE: &str = "data/rsviruses17900.1k.fastq.zst";

fn open(path: impl AsRef<Path>) -> File {
    File::open(path).unwrap()
}

fn cli_config(targets: &Path, options: &QueryOptions) -> ContainmentConfig {
    ContainmentConfig {
        targets_path: targets.to_path_buf(),
        background_paths: Vec::new(),
        sample_inputs: Vec::new(),
        sample_names: Vec::new(),
        kmer_length: options.kmer_length,
        smer_length: options.smer_length,
        threads: 1,
        output_path: None,
        quiet: true,
        abundance_thresholds: options.abundance_thresholds.clone(),
        discriminatory: options.discriminatory,
        individual: options.individual,
        limit_bp: options.limit_bp,
        sort_order: options.sort_order,
        dump_kmers_path: None,
        no_total: options.no_total,
        confidence: options.confidence,
        fraction: options.fraction,
        complexity: options.complexity,
    }
}

fn cli_query(
    targets: &Path,
    background: &[PathBuf],
    samples: Vec<(String, Input)>,
    options: &QueryOptions,
) -> String {
    let out = NamedTempFile::new().unwrap();
    let (sample_names, sample_inputs) = samples
        .into_iter()
        .map(|(name, input)| (name, vec![input]))
        .unzip();
    skope::run_query(&ContainmentConfig {
        background_paths: background.to_vec(),
        sample_inputs,
        sample_names,
        output_path: Some(out.path().to_path_buf()),
        ..cli_config(targets, options)
    })
    .unwrap();
    std::fs::read_to_string(out.path()).unwrap()
}

fn cli_build_query(targets: &Path, background: &[PathBuf], options: &QueryOptions) -> Vec<u8> {
    let out = NamedTempFile::new().unwrap();
    skope::run_build_query(&BuildQueryConfig {
        targets_path: targets.to_path_buf(),
        background_paths: background.to_vec(),
        kmer_length: options.kmer_length,
        smer_length: options.smer_length,
        individual: options.individual,
        positions: options.confidence,
        threads: 1,
        output_path: Some(out.path().to_path_buf()),
        quiet: true,
        fraction: options.fraction,
        complexity: options.complexity,
    })
    .unwrap();
    std::fs::read(out.path()).unwrap()
}

/// Session with each target file added under its derived name
fn session(targets: &[PathBuf], background: &[PathBuf], options: &QueryOptions) -> QuerySession {
    let mut session = QuerySession::new(options.clone()).unwrap();
    for path in targets {
        let name = derive_sample_name(path, false);
        session.add_targets(&name, open(path), |_| {}).unwrap();
    }
    for path in background {
        session.add_background(open(path), |_| {}).unwrap();
    }
    session
}

fn session_query(
    mut session: QuerySession,
    samples: Vec<(String, Input<File>)>,
) -> (QuerySession, String) {
    for (name, input) in samples {
        session.add_sample(&name, input, |_| {}).unwrap();
    }
    let tsv = session.tsv().unwrap();
    (session, tsv)
}

fn sample() -> Vec<(String, Input)> {
    vec![("sample".to_string(), Input::File(SAMPLE.into()))]
}

fn sample_stream() -> Vec<(String, Input<File>)> {
    vec![("sample".to_string(), Input::File(open(SAMPLE)))]
}

fn target_dir_files() -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = std::fs::read_dir(TARGET_DIR)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "fa"))
        .collect();
    files.sort();
    files
}

#[test]
fn test_session_matches_query_across_options() {
    let cases = [
        QueryOptions::default(),
        QueryOptions {
            individual: true,
            discriminatory: true,
            confidence: true,
            abundance_thresholds: vec![10, 2],
            sort_order: SortOrder::Target,
            ..Default::default()
        },
        QueryOptions {
            individual: true,
            kmer_length: 35,
            smer_length: 0,
            fraction: 0.5,
            complexity: 0.3,
            no_total: true,
            sort_order: SortOrder::Original,
            ..Default::default()
        },
        QueryOptions {
            individual: true,
            limit_bp: Some(40_000),
            ..Default::default()
        },
    ];
    for options in cases {
        let (_, tsv) = session_query(session(&[TARGETS.into()], &[], &options), sample_stream());
        assert_eq!(
            tsv,
            cli_query(Path::new(TARGETS), &[], sample(), &options),
            "{options:?}"
        );
        assert!(tsv.lines().count() > 2, "{tsv}");
    }
}

#[test]
fn test_session_targets_act_as_directory() {
    let options = QueryOptions {
        confidence: true,
        ..Default::default()
    };
    let (_, tsv) = session_query(session(&target_dir_files(), &[], &options), sample_stream());
    assert_eq!(
        tsv,
        cli_query(Path::new(TARGET_DIR), &[], sample(), &options)
    );
}

#[test]
fn test_session_masks_background() {
    let options = QueryOptions {
        individual: true,
        ..Default::default()
    };
    let background = vec![target_dir_files()[0].clone()];
    let (_, tsv) = session_query(
        session(&[TARGETS.into()], &background, &options),
        sample_stream(),
    );
    assert_eq!(
        tsv,
        cli_query(Path::new(TARGETS), &background, sample(), &options)
    );
}

#[test]
fn test_session_index_matches_build_query() {
    let dir = TempDir::new().unwrap();
    let options = QueryOptions {
        individual: true,
        confidence: true,
        fraction: 0.5,
        complexity: 0.2,
        ..Default::default()
    };
    let background = vec![target_dir_files()[0].clone()];
    let cli_index = cli_build_query(Path::new(TARGETS), &background, &options);
    let index = session(&[TARGETS.into()], &background, &options)
        .index()
        .unwrap();
    assert_eq!(index, cli_index);

    // The index serves both, and fraction and complexity can be left to it
    let path = dir.path().join("t.sk");
    std::fs::write(&path, &index).unwrap();
    let from_index = QueryOptions {
        fraction: 1.0,
        complexity: 0.0,
        ..options.clone()
    };
    let mut loaded = QuerySession::new(from_index.clone()).unwrap();
    loaded.load_index(&index).unwrap();
    let (_, tsv) = session_query(loaded, sample_stream());
    assert_eq!(tsv, cli_query(&path, &[], sample(), &from_index));

    let mut conflicting = QuerySession::new(QueryOptions {
        fraction: 0.25,
        ..options
    })
    .unwrap();
    assert!(conflicting.load_index(&index).is_err());
}

#[test]
fn test_session_pairs_mates() {
    let dir = TempDir::new().unwrap();
    let seq = std::fs::read_to_string(TARGETS)
        .unwrap()
        .lines()
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .collect::<String>();
    let qual = "I".repeat(150);
    let record = |id: String, seq: &str| format!("@{id}\n{seq}\n+\n{qual}\n");
    let (r1, r2): (Vec<String>, Vec<String>) = (0..400)
        .map(|i| {
            let (a, b) = (&seq[i * 7..i * 7 + 150], &seq[i * 5 + 300..i * 5 + 450]);
            (record(format!("p{i}/1"), a), record(format!("p{i}/2"), b))
        })
        .unzip();
    let il: String = r1.iter().zip(&r2).map(|(a, b)| format!("{a}{b}")).collect();
    let (r1, r1_short, r2) = (r1.concat(), r1[..100].concat(), r2.concat());
    let paths = ["r1.fq", "r2.fq", "il.fq"].map(|name| dir.path().join(name));
    for (path, data) in paths.iter().zip([&r1, &r2, &il]) {
        std::fs::write(path, data).unwrap();
    }
    let short = (r1_short.into_bytes(), r2.into_bytes());
    let [r1, r2, il] = paths;
    let options = QueryOptions::default();
    let samples = vec![
        ("paired".to_string(), Input::Paired(r1.clone(), r2.clone())),
        ("interleaved".to_string(), Input::Interleaved(il.clone())),
        ("single".to_string(), Input::File(il.clone())),
    ];
    let streams = vec![
        ("paired".to_string(), Input::Paired(open(&r1), open(&r2))),
        ("interleaved".to_string(), Input::Interleaved(open(&il))),
        ("single".to_string(), Input::File(open(&il))),
    ];
    let (_, tsv) = session_query(session(&[TARGETS.into()], &[], &options), streams);
    assert_eq!(tsv, cli_query(Path::new(TARGETS), &[], samples, &options));

    // Unequal mates fail as in the CLI
    let r1_short_path = dir.path().join("r1_short.fq");
    std::fs::write(&r1_short_path, &short.0).unwrap();
    let cli_err = skope::run_query(&ContainmentConfig {
        sample_inputs: vec![vec![Input::Paired(r1_short_path, r2.clone())]],
        sample_names: vec!["short".to_string()],
        ..cli_config(Path::new(TARGETS), &options)
    })
    .unwrap_err();
    let mut session = session(&[TARGETS.into()], &[], &options);
    let mates = Input::Paired(&short.0[..], &short.1[..]);
    let err = session.add_sample("short", mates, |_| {}).unwrap_err();
    assert_eq!(err.to_string(), cli_err.to_string());
}

#[test]
fn test_session_rejects_out_of_order_calls() {
    let options = QueryOptions::default();
    let mut session = session(&[TARGETS.into()], &[], &options);
    session
        .add_sample("a", Input::File(open(SAMPLE)), |_| {})
        .unwrap();
    assert!(
        session
            .add_sample("a", Input::File(open(SAMPLE)), |_| {})
            .is_err()
    );
    assert!(session.add_targets("t", open(TARGETS), |_| {}).is_err());
    assert!(session.add_background(open(TARGETS), |_| {}).is_err());
    assert!(session.index().is_err());
    assert!(
        QuerySession::new(QueryOptions {
            kmer_length: 30,
            ..Default::default()
        })
        .is_err()
    );
}

use skope::{
    BuildClassifyConfig, ClassifyConfig, ContainmentConfig, LengthHistogramConfig, SortOrder,
    discover_target_groups,
};
use std::path::PathBuf;
use tempfile::{NamedTempFile, TempDir};

#[test]
fn test_multisample_processing() {
    let config = ContainmentConfig {
        background_paths: Vec::new(),
        targets_path: PathBuf::from("data/zmrp21.viruses.fa"),
        sample_paths: vec![
            vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")],
            vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")],
        ],
        sample_names: vec!["sample1".to_string(), "sample2".to_string()],
        kmer_length: 31,
        smer_length: 15,
        complexity: 0.0,
        threads: 2,
        output_path: None,
        quiet: true,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_kmers_path: None,
        confidence: false,
        fraction: 1.0,
        no_total: false,
        individual: true,
    };

    assert!(skope::run_query(&config).is_ok());
}

#[test]
fn test_multisample_tsv_structure() {
    let temp_output = NamedTempFile::new().unwrap();
    let output_path = temp_output.path().to_path_buf();

    let config = ContainmentConfig {
        background_paths: Vec::new(),
        targets_path: PathBuf::from("data/zmrp21.viruses.fa"),
        sample_paths: vec![
            vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")],
            vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")],
        ],
        sample_names: vec!["sample_a".to_string(), "sample_b".to_string()],
        kmer_length: 31,
        smer_length: 15,
        complexity: 0.0,
        threads: 2,
        output_path: Some(output_path.clone()),
        quiet: true,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_kmers_path: None,
        confidence: false,
        fraction: 1.0,
        no_total: false,
        individual: true,
    };

    skope::run_query(&config).unwrap();

    let tsv_str = std::fs::read_to_string(&output_path).unwrap();
    let lines: Vec<&str> = tsv_str.lines().collect();

    // Check header
    assert!(lines[0].starts_with("target\tsample\t"));
    assert!(lines[0].contains("\tsample_seqs\tsample_bases"));

    // Count data rows per sample (excluding header and TOTAL rows)
    let sample_a_rows = lines
        .iter()
        .skip(1)
        .filter(|line| line.contains("\tsample_a\t") && !line.starts_with("TOTAL"))
        .count();
    let sample_b_rows = lines
        .iter()
        .skip(1)
        .filter(|line| line.contains("\tsample_b\t") && !line.starts_with("TOTAL"))
        .count();

    assert_eq!(sample_a_rows, 16, "Expected 16 target rows for sample_a");
    assert_eq!(sample_b_rows, 16, "Expected 16 target rows for sample_b");

    // Check TOTAL rows exist
    let total_rows = lines
        .iter()
        .filter(|line| line.starts_with("TOTAL\t"))
        .count();
    assert_eq!(total_rows, 2, "Expected 2 TOTAL rows (one per sample)");

    let last_sample_a = lines
        .iter()
        .rposition(|line| line.contains("\tsample_a\t"))
        .unwrap();
    let first_sample_b = lines
        .iter()
        .position(|line| line.contains("\tsample_b\t"))
        .unwrap();
    assert!(
        last_sample_a < first_sample_b,
        "samples should remain grouped in input order"
    );

    let header_cols: Vec<&str> = lines[0].split('\t').collect();
    let median_idx = header_cols
        .iter()
        .position(|&col| col == "median_nz_abundance")
        .expect("median_nz_abundance column not found");
    for total_row in lines.iter().filter(|line| line.starts_with("TOTAL\t")) {
        assert_eq!(total_row.split('\t').nth(median_idx), Some("-"));
    }

    // Check sample_seqs and sample_bases are non-zero and equal across both samples
    // (both samples use the same input file)
    let get_last_two_cols = |line: &str| -> (u64, u64) {
        let cols: Vec<&str> = line.split('\t').collect();
        let n = cols.len();
        (cols[n - 2].parse().unwrap(), cols[n - 1].parse().unwrap())
    };
    let first_a = lines
        .iter()
        .skip(1)
        .find(|l| l.contains("\tsample_a\t") && !l.starts_with("TOTAL"))
        .unwrap();
    let first_b = lines
        .iter()
        .skip(1)
        .find(|l| l.contains("\tsample_b\t") && !l.starts_with("TOTAL"))
        .unwrap();
    let (seqs_a, bases_a) = get_last_two_cols(first_a);
    let (seqs_b, bases_b) = get_last_two_cols(first_b);
    assert!(seqs_a > 0, "sample_seqs should be non-zero");
    assert!(bases_a > 0, "sample_bases should be non-zero");
    assert_eq!(seqs_a, seqs_b, "sample_seqs should be equal for same input");
    assert_eq!(
        bases_a, bases_b,
        "sample_bases should be equal for same input"
    );
}

#[test]
fn test_confidence_outputs_ani_and_patchiness_columns() {
    let temp_output = NamedTempFile::new().unwrap();

    let config = ContainmentConfig {
        background_paths: Vec::new(),
        targets_path: PathBuf::from("data/zmrp21.viruses.fa"),
        sample_paths: vec![vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")]],
        sample_names: vec!["sample".to_string()],
        kmer_length: 31,
        smer_length: 15,
        complexity: 0.0,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_kmers_path: None,
        confidence: true,
        fraction: 1.0,
        no_total: false,
        individual: true,
    };

    skope::run_query(&config).unwrap();
    let tsv_str = std::fs::read_to_string(temp_output.path()).unwrap();
    let lines: Vec<&str> = tsv_str.lines().collect();
    let header_cols: Vec<&str> = lines[0].split('\t').collect();
    let ani_est_idx = header_cols
        .iter()
        .position(|&col| col == "ani_est")
        .expect("ani_est column not found");
    let patchiness_idx = header_cols
        .iter()
        .position(|&col| col == "patchiness")
        .expect("patchiness column not found");

    let first_target = lines
        .iter()
        .skip(1)
        .find(|line| !line.starts_with("TOTAL\t"))
        .unwrap();
    let ani_est = first_target.split('\t').nth(ani_est_idx).unwrap();
    // ani_est is "-" when suppressed, otherwise a containment ANI in [0, 1]
    assert!(
        ani_est == "-" || (0.0..=1.0).contains(&ani_est.parse::<f64>().unwrap()),
        "unexpected ani_est value: {ani_est}"
    );

    let patchiness = first_target.split('\t').nth(patchiness_idx).unwrap();
    assert!(patchiness == "-" || patchiness.contains('|'));
}

#[test]
fn test_sort_target() {
    let temp_output = NamedTempFile::new().unwrap();
    let config = ContainmentConfig {
        background_paths: Vec::new(),
        targets_path: PathBuf::from("data/zmrp21.viruses.fa"),
        sample_paths: vec![vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")]],
        sample_names: vec!["test".to_string()],
        kmer_length: 31,
        smer_length: 15,
        complexity: 0.0,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Target,
        dump_kmers_path: None,
        confidence: false,
        fraction: 1.0,
        no_total: false,
        individual: true,
    };

    skope::run_query(&config).unwrap();
    let tsv_str = std::fs::read_to_string(temp_output.path()).unwrap();
    let lines: Vec<&str> = tsv_str.lines().collect();

    // Extract target names (first column), skip header and TOTAL rows
    let targets: Vec<&str> = lines
        .iter()
        .skip(1)
        .filter(|line| !line.starts_with("TOTAL"))
        .map(|line| line.split('\t').next().unwrap())
        .collect();

    // Verify sorted alphabetically by target
    for i in 0..targets.len() - 1 {
        assert!(
            targets[i] <= targets[i + 1],
            "Targets not sorted: {} > {}",
            targets[i],
            targets[i + 1]
        );
    }
}

#[test]
fn test_sort_containment() {
    let temp_output = NamedTempFile::new().unwrap();
    let config = ContainmentConfig {
        background_paths: Vec::new(),
        targets_path: PathBuf::from("data/zmrp21.viruses.fa"),
        sample_paths: vec![vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")]],
        sample_names: vec!["test".to_string()],
        kmer_length: 31,
        smer_length: 15,
        complexity: 0.0,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Containment,
        dump_kmers_path: None,
        confidence: false,
        fraction: 1.0,
        no_total: false,
        individual: true,
    };

    skope::run_query(&config).unwrap();
    let tsv_str = std::fs::read_to_string(temp_output.path()).unwrap();
    let lines: Vec<&str> = tsv_str.lines().collect();

    // Find containment1 column index from header
    let header_cols: Vec<&str> = lines[0].split('\t').collect();
    let containment_idx = header_cols
        .iter()
        .position(|&c| c == "containment1")
        .expect("containment1 column not found");

    // Extract containment values, skip header and TOTAL rows
    let containments: Vec<f64> = lines
        .iter()
        .skip(1)
        .filter(|line| !line.starts_with("TOTAL"))
        .map(|line| {
            line.split('\t')
                .nth(containment_idx)
                .unwrap()
                .parse::<f64>()
                .unwrap()
        })
        .collect();

    // Verify sorted by containment descending
    for i in 0..containments.len() - 1 {
        assert!(
            containments[i] >= containments[i + 1],
            "Containments not sorted descending: {} < {}",
            containments[i],
            containments[i + 1]
        );
    }

    // With non-overlapping kmers, absolute containment is lower than before,
    // but the top hit should still be non-zero.
    assert!(
        containments[0] > 0.0,
        "Highest containment {} should be > 0.0",
        containments[0]
    );
}

#[test]
fn test_length_histogram() {
    // Wrap the viruses fasta in a temp dir so it's treated as a single group "viruses"
    let targets_path = TempDir::new().unwrap();
    std::fs::copy(
        "data/zmrp21.viruses.fa",
        targets_path.path().join("viruses.fa"),
    )
    .unwrap();

    let temp_output = NamedTempFile::new().unwrap();

    let config = LengthHistogramConfig {
        individual: false,
        targets_path: targets_path.path().to_path_buf(),
        sample_paths: vec![vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")]],
        sample_names: vec!["test".to_string()],
        kmer_length: 31,
        smer_length: 15,
        complexity: 0.0,
        abs_threshold: 1,
        rel_threshold: 0.0,
        discriminatory: false,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        limit_bp: None,
        no_filter: false,
    };

    skope::run_lenhist(&config).unwrap();

    // Verify TSV output has new per-group header and data rows
    let content = std::fs::read_to_string(temp_output.path()).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(
        lines[0],
        "sample\tgroup\tlength\tcount\ttotal_seqs_processed\ttotal_bp_processed\tgroup_seqs\tgroup_bases",
        "Header should include group column"
    );
    assert!(lines.len() > 1, "Should have data rows");
    assert!(
        lines[1].starts_with("test\tviruses\t"),
        "Data rows should start with sample\\tgroup, got: {}",
        lines[1]
    );
}

#[test]
fn test_length_histogram_all_seqs() {
    let temp_output = NamedTempFile::new().unwrap();

    let config = LengthHistogramConfig {
        individual: false,
        targets_path: PathBuf::from("-"),
        sample_paths: vec![vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")]],
        sample_names: vec!["test".to_string()],
        kmer_length: 31,
        smer_length: 15,
        complexity: 0.0,
        abs_threshold: 1,
        rel_threshold: 0.0,
        discriminatory: false,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        limit_bp: None,
        no_filter: true,
    };

    skope::run_lenhist(&config).unwrap();

    let content = std::fs::read_to_string(temp_output.path()).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert!(lines.len() > 1, "Should have data rows");
    assert!(
        lines[1].starts_with("test\tall\t"),
        "All reads should go to the 'all' bucket, got: {}",
        lines[1]
    );
}

fn write_fasta(path: &std::path::Path, header: &str, seq: &str) {
    use std::io::Write;
    let mut f = std::fs::File::create(path).unwrap();
    writeln!(f, ">{}", header).unwrap();
    writeln!(f, "{}", seq).unwrap();
}

fn classify_to(
    targets: &std::path::Path,
    sample_paths: Vec<Vec<PathBuf>>,
    sample_names: &[&str],
    output: &std::path::Path,
    per_seq: bool,
    limit_bp: Option<u64>,
) {
    skope::run_classification(&ClassifyConfig {
        targets_path: targets.to_path_buf(),
        individual: false,
        sample_paths,
        sample_names: sample_names.iter().map(|name| name.to_string()).collect(),
        kmer_length: 15,
        smer_length: 7,
        complexity: 0.0,
        abs_threshold: 1,
        rel_threshold: 0.0,
        threads: 1,
        limit_bp,
        output_path: Some(output.to_path_buf()),
        per_seq,
        discriminatory: false,
        quiet: true,
    })
    .unwrap();
}

const SEQ_A: &str = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAA";
const SEQ_B: &str = "GTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCCCTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGAATGTGTACATACGCTCTTACTGCGGTCGCGTCTAATAATATACATTTGCTTCGTTGACTAGCAACCCAGGGCTATAGCTA";

#[test]
fn test_discover_target_groups_mixed_layout() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();

    // Top-level file
    write_fasta(&root.join("class_a.fa"), "a1", SEQ_A);
    // Subdirectory with two fastx files
    std::fs::create_dir(root.join("class_b")).unwrap();
    write_fasta(&root.join("class_b/part1.fa"), "b1", SEQ_A);
    write_fasta(&root.join("class_b/part2.fa"), "b2", SEQ_B);
    // Hidden top-level entries (skipped)
    write_fasta(&root.join(".hidden.fa"), "h", SEQ_A);
    std::fs::create_dir(root.join(".hidden_dir")).unwrap();
    // Non-fastx file (skipped)
    std::fs::write(root.join("README.txt"), b"ignore me").unwrap();
    // Hidden file inside subdir (skipped)
    write_fasta(&root.join("class_b/.skip.fa"), "x", SEQ_A);

    let groups = discover_target_groups(root).unwrap();
    assert_eq!(groups.len(), 2);
    assert_eq!(groups[0].name, "class_a");
    assert_eq!(groups[0].files.len(), 1);
    assert!(groups[0].files[0].ends_with("class_a.fa"));
    assert_eq!(groups[1].name, "class_b");
    assert_eq!(groups[1].files.len(), 2);
    assert!(groups[1].files[0].ends_with("part1.fa"));
    assert!(groups[1].files[1].ends_with("part2.fa"));
}

#[test]
fn test_discover_target_groups_nested_subdir_errors() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    std::fs::create_dir_all(root.join("class_a/nested")).unwrap();
    write_fasta(&root.join("class_a/part1.fa"), "a", SEQ_A);

    let err = discover_target_groups(root).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("Nested subdirectory"), "got: {}", msg);
}

#[test]
fn test_discover_target_groups_empty_subdir_errors() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    std::fs::create_dir(root.join("class_a")).unwrap();
    std::fs::write(root.join("class_a/notes.txt"), b"no fastx here").unwrap();

    let err = discover_target_groups(root).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("no fastx files"), "got: {}", msg);
    assert!(msg.contains("class_a"), "got: {}", msg);
}

#[test]
fn test_discover_target_groups_duplicate_name_errors() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    write_fasta(&root.join("foo.fa"), "f", SEQ_A);
    std::fs::create_dir(root.join("foo")).unwrap();
    write_fasta(&root.join("foo/inner.fa"), "i", SEQ_B);

    let err = discover_target_groups(root).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("Duplicate group name"), "got: {}", msg);
    assert!(msg.contains("foo"), "got: {}", msg);
}

#[test]
fn test_query_directory_mixed_layout() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    write_fasta(&root.join("class_a.fa"), "a1", SEQ_A);
    std::fs::create_dir(root.join("class_b")).unwrap();
    write_fasta(&root.join("class_b/part1.fa"), "b1", SEQ_A);
    write_fasta(&root.join("class_b/part2.fa"), "b2", SEQ_B);

    let sample = NamedTempFile::new().unwrap();
    write_fasta(sample.path(), "s1", SEQ_A);

    let output = NamedTempFile::new().unwrap();
    let config = ContainmentConfig {
        background_paths: Vec::new(),
        targets_path: root.to_path_buf(),
        sample_paths: vec![vec![sample.path().to_path_buf()]],
        sample_names: vec!["s".to_string()],
        kmer_length: 15,
        smer_length: 7,
        complexity: 0.0,
        threads: 1,
        output_path: Some(output.path().to_path_buf()),
        quiet: true,
        abundance_thresholds: vec![1],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_kmers_path: None,
        confidence: false,
        fraction: 1.0,
        no_total: true,
        individual: false,
    };
    skope::run_query(&config).unwrap();

    let content = std::fs::read_to_string(output.path()).unwrap();
    let target_names: Vec<&str> = content
        .lines()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    assert!(
        target_names.contains(&"class_a"),
        "missing class_a: {:?}",
        target_names
    );
    assert!(
        target_names.contains(&"class_b"),
        "missing class_b: {:?}",
        target_names
    );
    assert_eq!(target_names.len(), 2);
}

// Two targets sharing a common region but with distinct unique tails. The shared
// region yields cross-target kmers; the unique tails do not.
const DISC_COMMON: &str = "GATTACAGGCATCCTAGCTAGGACTTGCAACATGCTTAGCCATGGAACTGTCCAGTTACGGATCCTAGGCATTAGCCAGTTCATGGACTTAGCGGATCCTA";
const DISC_UNIQUE_A: &str = "TTGCAACGGTACCATTAGCGGATCCTTAGCAACATGCTTAGCCATGGAACTGTCCAGTTACGGATCCTAGGCATTAGCCAGTTCATGGACTTAGCGGATCC";
const DISC_UNIQUE_B: &str = "CCAGTTACGGATCCTAGGCATTAGCCAGTTCATGGACTTAGCGGATCCTAGCTAGGACTTGCAACATGCTTAGCCATGGAACTGTCCAGTTACGGATCCTA";

// Parse a --dump-kmers TSV into target -> set of k-mers (col 1 -> col 3)
fn parse_dump(
    path: &std::path::Path,
) -> std::collections::HashMap<String, std::collections::HashSet<String>> {
    let content = std::fs::read_to_string(path).unwrap();
    let mut map: std::collections::HashMap<String, std::collections::HashSet<String>> =
        std::collections::HashMap::new();
    for line in content.lines() {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 3, "expected 3 columns, got: {:?}", cols);
        map.entry(cols[0].to_string())
            .or_default()
            .insert(cols[2].to_string());
    }
    map
}

#[test]
fn test_dump_kmers_respects_discriminatory() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    let seq_a = format!("{}{}", DISC_COMMON, DISC_UNIQUE_A);
    let seq_b = format!("{}{}", DISC_COMMON, DISC_UNIQUE_B);
    write_fasta(&root.join("tA.fa"), "tA", &seq_a);
    write_fasta(&root.join("tB.fa"), "tB", &seq_b);

    let sample = NamedTempFile::new().unwrap();
    {
        use std::io::Write;
        let mut f = std::fs::File::create(sample.path()).unwrap();
        writeln!(f, ">r1\n{}\n>r2\n{}", seq_a, seq_b).unwrap();
    }

    let make_config = |discriminatory: bool, dump: PathBuf| ContainmentConfig {
        background_paths: Vec::new(),
        targets_path: root.to_path_buf(),
        sample_paths: vec![vec![sample.path().to_path_buf()]],
        sample_names: vec!["s".to_string()],
        kmer_length: 15,
        smer_length: 7,
        complexity: 0.0,
        threads: 1,
        output_path: None,
        quiet: true,
        abundance_thresholds: vec![1],
        discriminatory,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_kmers_path: Some(dump),
        confidence: false,
        fraction: 1.0,
        no_total: true,
        individual: false,
    };

    let plain_dump = NamedTempFile::new().unwrap();
    skope::run_query(&make_config(false, plain_dump.path().to_path_buf())).unwrap();
    let plain = parse_dump(plain_dump.path());

    let disc_dump = NamedTempFile::new().unwrap();
    skope::run_query(&make_config(true, disc_dump.path().to_path_buf())).unwrap();
    let disc = parse_dump(disc_dump.path());

    let plain_a = &plain["tA"];
    let plain_b = &plain["tB"];
    let shared: Vec<_> = plain_a.intersection(plain_b).collect();
    assert!(
        !shared.is_empty(),
        "expected cross-target shared kmers in plain dump"
    );

    let disc_a = &disc["tA"];
    let disc_b = &disc["tB"];
    let disc_shared: Vec<_> = disc_a.intersection(disc_b).collect();
    assert!(
        disc_shared.is_empty(),
        "discriminatory dump must contain no cross-target shared kmers, found: {:?}",
        disc_shared
    );

    let plain_total = plain_a.len() + plain_b.len();
    let disc_total = disc_a.len() + disc_b.len();
    assert!(
        disc_total < plain_total,
        "discriminatory dump ({}) should have fewer kmers than plain ({})",
        disc_total,
        plain_total
    );
}

#[test]
fn test_classify_build_mixed_layout() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    write_fasta(&root.join("class_a.fa"), "a1", SEQ_A);
    std::fs::create_dir(root.join("class_b")).unwrap();
    write_fasta(&root.join("class_b/part1.fa"), "b1", SEQ_A);
    write_fasta(&root.join("class_b/part2.fa"), "b2", SEQ_B);

    let idx_out = NamedTempFile::new().unwrap();
    let config = BuildClassifyConfig {
        individual: false,
        targets_path: root.to_path_buf(),
        kmer_length: 15,
        smer_length: 7,
        complexity: 0.0,
        threads: 1,
        output_path: Some(idx_out.path().to_path_buf()),
        quiet: true,
    };
    skope::run_build_classify(&config).unwrap();

    // Index file should exist and be non-empty
    let meta = std::fs::metadata(idx_out.path()).unwrap();
    assert!(meta.len() > 0);

    // Run classify against it with a sample whose seq matches class_b
    let sample = NamedTempFile::new().unwrap();
    write_fasta(sample.path(), "s1", SEQ_B);

    let classify_out = NamedTempFile::new().unwrap();
    let cfg = ClassifyConfig {
        individual: false,
        targets_path: idx_out.path().to_path_buf(),
        sample_paths: vec![vec![sample.path().to_path_buf()]],
        sample_names: vec!["s".to_string()],
        kmer_length: 15,
        smer_length: 7,
        complexity: 0.0,
        abs_threshold: 1,
        rel_threshold: 0.0,
        threads: 1,
        limit_bp: None,
        output_path: Some(classify_out.path().to_path_buf()),
        per_seq: false,
        discriminatory: false,
        quiet: true,
    };
    skope::run_classification(&cfg).unwrap();

    let content = std::fs::read_to_string(classify_out.path()).unwrap();
    // class_b should appear in the output (the sample matches it)
    assert!(content.contains("class_b"), "classify output: {}", content);
}

#[test]
fn test_classify_too_many_groups_errors() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    for i in 0..129 {
        write_fasta(&root.join(format!("g{:03}.fa", i)), "x", SEQ_A);
    }

    let idx_out = NamedTempFile::new().unwrap();
    let config = BuildClassifyConfig {
        individual: false,
        targets_path: root.to_path_buf(),
        kmer_length: 15,
        smer_length: 7,
        complexity: 0.0,
        threads: 1,
        output_path: Some(idx_out.path().to_path_buf()),
        quiet: true,
    };
    let err = skope::run_build_classify(&config).unwrap_err();
    let msg = format!("{:#}", err);
    assert!(msg.contains("Too many groups"), "got: {}", msg);
}

#[test]
fn test_query_cli_omits_noop_options() {
    let query_help = std::process::Command::new(env!("CARGO_BIN_EXE_skope"))
        .args(["query", "--help"])
        .output()
        .unwrap();
    assert!(query_help.status.success());
    let query_help = String::from_utf8(query_help.stdout).unwrap();
    assert!(!query_help.contains("--positions"));
    assert!(query_help.contains("possible values: containment, target, input"));

    let build_help = std::process::Command::new(env!("CARGO_BIN_EXE_skope"))
        .args(["index", "build-query", "--help"])
        .output()
        .unwrap();
    assert!(build_help.status.success());
    assert!(
        String::from_utf8(build_help.stdout)
            .unwrap()
            .contains("--positions")
    );
}

#[test]
fn test_query_cli_k_s_defaults_and_index_precedence() {
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("target.fa");
    let sample = dir.path().join("sample.fa");
    let index = dir.path().join("target.sk");
    write_fasta(&target, "target", SEQ_A);
    write_fasta(&sample, "sample", SEQ_A);
    build_index(target.clone(), vec![], index.clone());

    let raw_out = dir.path().join("raw.tsv");
    let raw = std::process::Command::new(env!("CARGO_BIN_EXE_skope"))
        .args(["query", "-q", "-o"])
        .arg(raw_out)
        .arg(&target)
        .arg(&sample)
        .output()
        .unwrap();
    assert!(raw.status.success());
    assert!(String::from_utf8_lossy(&raw.stderr).contains("options: k=31, s=9"));

    let matching_out = dir.path().join("matching.tsv");
    let matching = std::process::Command::new(env!("CARGO_BIN_EXE_skope"))
        .args(["query", "-q", "-k", "15", "-o"])
        .arg(matching_out)
        .arg(&index)
        .arg(&sample)
        .output()
        .unwrap();
    assert!(matching.status.success());
    assert!(!String::from_utf8_lossy(&matching.stderr).contains("CLI k/s ignored"));

    let conflicting_out = dir.path().join("conflicting.tsv");
    let conflicting = std::process::Command::new(env!("CARGO_BIN_EXE_skope"))
        .args(["query", "-k", "31", "-o"])
        .arg(conflicting_out)
        .arg(&index)
        .arg(&sample)
        .output()
        .unwrap();
    assert!(conflicting.status.success());
    assert!(String::from_utf8_lossy(&conflicting.stderr).contains("CLI k/s ignored"));
}

#[cfg(unix)]
#[test]
fn test_fifo_sample_input() {
    use nix::sys::stat;
    use nix::unistd;
    use std::io::Write;

    let tmp_dir = TempDir::new().unwrap();
    let fifo_path = tmp_dir.path().join("test.fastq");

    // Create named pipe
    unistd::mkfifo(&fifo_path, stat::Mode::S_IRWXU).unwrap();

    let fifo_path_clone = fifo_path.clone();
    let writer = std::thread::spawn(move || {
        // Write minimal FASTQ data
        let fastq = b"@read1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
        let mut f = std::fs::File::create(&fifo_path_clone).unwrap();
        f.write_all(fastq).unwrap();
    });

    // Run the binary with the FIFO as sample input
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_skope"))
        .args([
            "query",
            "data/zmrp21.viruses.fa",
            fifo_path.to_str().unwrap(),
            "-q",
        ])
        .output()
        .expect("failed to execute skope");

    writer.join().unwrap();

    assert!(
        output.status.success(),
        "skope query with FIFO failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<&str> = stdout.lines().collect();
    assert!(
        lines.len() > 1,
        "Expected header + data rows from FIFO input"
    );
    assert!(
        lines[0].starts_with("target\tsample\t"),
        "Expected TSV header"
    );
}

fn build_index_from_stdin(subcommand: &str, expected_kind: skope::IndexKind) {
    use std::io::Write;
    use std::process::{Command, Stdio};

    let dir = TempDir::new().unwrap();
    let index = dir.path().join("stdin.sk");
    let mut child = Command::new(env!("CARGO_BIN_EXE_skope"))
        .args([
            "index", subcommand, "-", "-k", "15", "-s", "7", "-t", "1", "-q", "-o",
        ])
        .arg(&index)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to execute skope");

    child
        .stdin
        .take()
        .unwrap()
        .write_all(format!(">target\n{SEQ_A}\n").as_bytes())
        .unwrap();

    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "skope index {subcommand} with stdin failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(skope::read_index_kind(&index), Some(expected_kind));
}

#[test]
fn test_build_query_accepts_stdin_targets() {
    build_index_from_stdin("build-query", skope::IndexKind::Query);
}

#[test]
fn test_build_classify_accepts_stdin_targets() {
    build_index_from_stdin("build-classify", skope::IndexKind::Classify);
}

fn build_index(targets: PathBuf, background: Vec<PathBuf>, out: PathBuf) {
    skope::run_build_query(&skope::BuildQueryConfig {
        targets_path: targets,
        background_paths: background,
        kmer_length: 15,
        smer_length: 7,
        complexity: 0.0,
        individual: false,
        positions: false,
        threads: 1,
        output_path: Some(out),
        quiet: true,
        fraction: 1.0,
    })
    .unwrap();
}

fn query_to_tsv(targets_path: PathBuf, sample: &std::path::Path, out: &std::path::Path) {
    skope::run_query(&ContainmentConfig {
        background_paths: Vec::new(),
        targets_path,
        sample_paths: vec![vec![sample.to_path_buf()]],
        sample_names: vec!["s".to_string()],
        kmer_length: 15,
        smer_length: 7,
        complexity: 0.0,
        threads: 1,
        output_path: Some(out.to_path_buf()),
        quiet: true,
        abundance_thresholds: vec![1],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_kmers_path: None,
        confidence: false,
        fraction: 1.0,
        no_total: true,
        individual: false,
    })
    .unwrap();
}

#[test]
fn test_query_index_matches_fastx() {
    let dir = TempDir::new().unwrap();
    let (target, sample, index) = (
        dir.path().join("t.fa"),
        dir.path().join("s.fa"),
        dir.path().join("t.sk"),
    );
    write_fasta(&target, "t1", SEQ_A);
    write_fasta(&sample, "s1", SEQ_A);
    build_index(target.clone(), vec![], index.clone());

    let bytes = std::fs::read(&index).unwrap();
    assert_eq!(&bytes[..4], skope::INDEX_MAGIC);
    assert_eq!(bytes[4], skope::IndexKind::Query as u8);

    let (i, f) = (dir.path().join("i.tsv"), dir.path().join("f.tsv"));
    query_to_tsv(index, &sample, &i);
    query_to_tsv(target, &sample, &f);
    let from_index = std::fs::read_to_string(i).unwrap();
    assert_eq!(from_index, std::fs::read_to_string(f).unwrap());
    // Equality is vacuous unless the target yielded k-mers
    let tk = from_index.lines().next().unwrap();
    let tk = tk.split('\t').position(|h| h == "target_kmers").unwrap();
    assert_ne!(
        from_index.lines().nth(1).unwrap().split('\t').nth(tk),
        Some("0")
    );
}

#[test]
fn test_query_index_positions_required_only_when_consumed() {
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("t.fa");
    let sample = dir.path().join("s.fa");
    let positioned_index = dir.path().join("positioned.sk");
    let plain_index = dir.path().join("plain.sk");
    write_fasta(&target, "t1", SEQ_A);
    write_fasta(&sample, "s1", SEQ_A);

    let build = |output_path: PathBuf, positions: bool| {
        skope::run_build_query(&skope::BuildQueryConfig {
            targets_path: target.clone(),
            background_paths: vec![],
            kmer_length: 15,
            smer_length: 7,
            complexity: 0.0,
            individual: false,
            positions,
            threads: 1,
            output_path: Some(output_path),
            quiet: true,
            fraction: 1.0,
        })
        .unwrap();
    };
    build(positioned_index.clone(), true);
    build(plain_index.clone(), false);

    let query = |targets_path: PathBuf, output_path: PathBuf| {
        skope::run_query(&ContainmentConfig {
            targets_path,
            background_paths: vec![],
            sample_paths: vec![vec![sample.clone()]],
            sample_names: vec!["s".to_string()],
            kmer_length: 15,
            smer_length: 7,
            complexity: 0.0,
            threads: 1,
            output_path: Some(output_path),
            quiet: true,
            abundance_thresholds: vec![10],
            discriminatory: false,
            individual: false,
            limit_bp: None,
            sort_order: SortOrder::Original,
            dump_kmers_path: None,
            no_total: true,
            confidence: true,
            fraction: 1.0,
        })
    };

    assert!(query(positioned_index, dir.path().join("positioned.tsv")).is_ok());
    let error = query(plain_index, dir.path().join("plain.tsv")).unwrap_err();
    assert!(error.to_string().contains("built without --positions"));
}

// Deterministic pseudo-random DNA for thinning tests
fn pseudo_dna_string(n: usize, seed: u64) -> String {
    let mut x = seed | 1;
    let bases = *b"ACGT";
    (0..n)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            bases[((x >> 33) & 0b11) as usize] as char
        })
        .collect()
}

fn build_index_frac(targets: PathBuf, out: PathBuf, fraction: f64) {
    skope::run_build_query(&skope::BuildQueryConfig {
        targets_path: targets,
        background_paths: vec![],
        kmer_length: 15,
        smer_length: 7,
        complexity: 0.0,
        individual: false,
        positions: false,
        threads: 1,
        output_path: Some(out),
        quiet: true,
        fraction,
    })
    .unwrap();
}

fn query_frac(
    targets_path: PathBuf,
    sample: &std::path::Path,
    out: &std::path::Path,
    fraction: f64,
) {
    skope::run_query(&ContainmentConfig {
        background_paths: Vec::new(),
        targets_path,
        sample_paths: vec![vec![sample.to_path_buf()]],
        sample_names: vec!["s".to_string()],
        kmer_length: 15,
        smer_length: 7,
        complexity: 0.0,
        threads: 1,
        output_path: Some(out.to_path_buf()),
        quiet: true,
        abundance_thresholds: vec![1],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_kmers_path: None,
        confidence: false,
        fraction,
        no_total: true,
        individual: false,
    })
    .unwrap();
}

// Index-time and query-time thinning must yield identical containment
#[test]
fn test_thinning_index_matches_fastx() {
    let dir = TempDir::new().unwrap();
    let seq = pseudo_dna_string(20_000, 3);
    let (target, sample, index) = (
        dir.path().join("t.fa"),
        dir.path().join("s.fa"),
        dir.path().join("t01.sk"),
    );
    write_fasta(&target, "t1", &seq);
    write_fasta(&sample, "s1", &seq);
    build_index_frac(target.clone(), index.clone(), 0.1);

    // Header stores the fraction
    let (_k, _s, frac, _c) = skope::read_query_index_meta(&index).unwrap();
    assert!((frac - 0.1).abs() < 1e-4, "stored fraction {frac} != 0.1");

    let (i, f) = (dir.path().join("i.tsv"), dir.path().join("f.tsv"));
    query_frac(index, &sample, &i, 1.0); // Index fraction wins, CLI default ignored
    query_frac(target, &sample, &f, 0.1);
    assert_eq!(
        std::fs::read_to_string(i).unwrap(),
        std::fs::read_to_string(f).unwrap()
    );
}

#[test]
fn test_query_index_background_masks_everything() {
    let dir = TempDir::new().unwrap();
    let (target, sample, index) = (
        dir.path().join("t.fa"),
        dir.path().join("s.fa"),
        dir.path().join("m.sk"),
    );
    write_fasta(&target, "t1", SEQ_A);
    write_fasta(&sample, "s1", SEQ_A);
    build_index(target.clone(), vec![target.clone()], index.clone());

    let out = dir.path().join("o.tsv");
    query_to_tsv(index, &sample, &out);
    let content = std::fs::read_to_string(out).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    let tk = lines[0]
        .split('\t')
        .position(|h| h == "target_kmers")
        .unwrap();
    assert_eq!(lines[1].split('\t').nth(tk).unwrap(), "0");
}

// ── Uniform <TARGETS> resolution ──────────────────────────────────────────────

fn lenhist_groups(
    targets_path: PathBuf,
    individual: bool,
    sample: &std::path::Path,
) -> Vec<String> {
    let out = NamedTempFile::new().unwrap();
    skope::run_lenhist(&LengthHistogramConfig {
        targets_path,
        individual,
        sample_paths: vec![vec![sample.to_path_buf()]],
        sample_names: vec!["s".to_string()],
        kmer_length: 15,
        smer_length: 7,
        complexity: 0.0,
        abs_threshold: 1,
        rel_threshold: 0.0,
        discriminatory: false,
        threads: 1,
        output_path: Some(out.path().to_path_buf()),
        quiet: true,
        limit_bp: None,
        no_filter: false,
    })
    .unwrap();

    let content = std::fs::read_to_string(out.path()).unwrap();
    let mut groups: Vec<String> = content
        .lines()
        .skip(1)
        .map(|line| line.split('\t').nth(1).unwrap().to_string())
        .collect();
    groups.sort();
    groups.dedup();
    groups
}

#[test]
fn test_lenhist_accepts_bare_fastx_file() {
    // Regression for fastx files misread as classification indexes
    let dir = TempDir::new().unwrap();
    let (targets, sample) = (dir.path().join("segments.fa"), dir.path().join("s.fa"));
    let seq = pseudo_dna_string(500, 11);
    write_fasta(&targets, "seg_L", &seq);
    write_fasta(&sample, "r1", &seq);

    // Records merge into one group named after the file
    assert_eq!(
        lenhist_groups(targets, false, &sample),
        vec!["segments".to_string()]
    );
}

#[test]
fn test_lenhist_individual_splits_records() {
    let dir = TempDir::new().unwrap();
    let (targets, sample) = (dir.path().join("segments.fa"), dir.path().join("s.fa"));
    let (seg_l, seg_m) = (pseudo_dna_string(500, 11), pseudo_dna_string(500, 29));
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&targets).unwrap();
        writeln!(f, ">seg_L\n{seg_l}\n>seg_M\n{seg_m}").unwrap();
    }
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&sample).unwrap();
        writeln!(f, ">r1\n{seg_l}\n>r2\n{seg_m}").unwrap();
    }

    assert_eq!(
        lenhist_groups(targets.clone(), false, &sample),
        vec!["segments".to_string()],
        "without --individual the records merge into one group"
    );
    assert_eq!(
        lenhist_groups(targets, true, &sample),
        vec!["seg_L".to_string(), "seg_M".to_string()],
        "with --individual each record becomes its own group"
    );
}

#[test]
fn test_classify_accepts_bare_fastx_file() {
    let dir = TempDir::new().unwrap();
    let (targets, sample, out) = (
        dir.path().join("refs.fa"),
        dir.path().join("s.fa"),
        dir.path().join("o.tsv"),
    );
    let seq = pseudo_dna_string(500, 11);
    write_fasta(&targets, "r1", &seq);
    write_fasta(&sample, "s1", &seq);

    skope::run_classification(&ClassifyConfig {
        targets_path: targets,
        individual: false,
        sample_paths: vec![vec![sample]],
        sample_names: vec!["s".to_string()],
        kmer_length: 15,
        smer_length: 7,
        complexity: 0.0,
        abs_threshold: 1,
        rel_threshold: 0.0,
        threads: 1,
        limit_bp: None,
        output_path: Some(out.clone()),
        per_seq: false,
        discriminatory: false,
        quiet: true,
    })
    .unwrap();

    // The matching read lands in the file-named group
    let content = std::fs::read_to_string(out).unwrap();
    let row = content
        .lines()
        .find(|line| line.split('\t').nth(1) == Some("refs"))
        .unwrap_or_else(|| panic!("no 'refs' group row in: {content}"));
    assert_eq!(
        row.split('\t').nth(3),
        Some("1"),
        "expected the read classified into 'refs': {content}"
    );
}

#[test]
fn test_classify_limit_spans_files_and_resets_per_sample() {
    let dir = TempDir::new().unwrap();
    let targets = dir.path().join("targets.fa");
    let first = dir.path().join("first.fa");
    let second = dir.path().join("second.fa");
    write_fasta(&targets, "target", SEQ_A);
    write_fasta(&first, "first", SEQ_A);
    write_fasta(&second, "second", SEQ_A);

    let summary_out = dir.path().join("summary.tsv");
    classify_to(
        &targets,
        vec![vec![first.clone(), second.clone()]],
        &["sample"],
        &summary_out,
        false,
        Some(1),
    );

    let summary_seqs: u64 = std::fs::read_to_string(summary_out)
        .unwrap()
        .lines()
        .skip(1)
        .map(|line| line.split('\t').nth(3).unwrap().parse::<u64>().unwrap())
        .sum();
    assert_eq!(summary_seqs, 1, "summary mode should skip the second file");

    let per_seq_out = dir.path().join("per-seq.tsv");
    classify_to(
        &targets,
        vec![vec![first.clone(), second.clone()], vec![first, second]],
        &["sample_a", "sample_b"],
        &per_seq_out,
        true,
        Some(1),
    );

    let output = std::fs::read_to_string(per_seq_out).unwrap();
    let rows: Vec<&str> = output.lines().skip(1).collect();
    assert_eq!(rows.len(), 2, "each sample should receive its own limit");
    assert!(rows.iter().any(|row| row.starts_with("sample_a\tfirst\t")));
    assert!(rows.iter().any(|row| row.starts_with("sample_b\tfirst\t")));
    assert!(!output.contains("\tsecond\t"));
}

#[test]
fn test_classify_modes_report_all_outcomes() {
    let dir = TempDir::new().unwrap();
    let targets = dir.path().join("targets");
    let sample = dir.path().join("sample.fa");
    std::fs::create_dir(&targets).unwrap();

    let seq_a = pseudo_dna_string(500, 11);
    let seq_b = pseudo_dna_string(500, 29);
    let unrelated = pseudo_dna_string(500, 47);
    write_fasta(&targets.join("a.fa"), "a", &seq_a);
    write_fasta(&targets.join("b.fa"), "b", &seq_b);
    std::fs::write(
        &sample,
        format!(
            ">classified\n{seq_a}\n>ambiguous\n{seq_a}NNNNNNNNNNNNNNN{seq_b}\n>unclassified\n{unrelated}\n"
        ),
    )
    .unwrap();

    let summary = dir.path().join("summary.tsv");
    classify_to(
        &targets,
        vec![vec![sample.clone()]],
        &["sample"],
        &summary,
        false,
        None,
    );
    let summary = std::fs::read_to_string(summary).unwrap();
    for (group, seqs) in [
        ("a", "1"),
        ("b", "0"),
        ("ambiguous", "1"),
        ("unclassified", "1"),
    ] {
        let row = summary
            .lines()
            .find(|line| line.split('\t').nth(1) == Some(group))
            .unwrap();
        assert_eq!(row.split('\t').nth(3), Some(seqs));
    }

    let per_seq = dir.path().join("per-seq.tsv");
    classify_to(
        &targets,
        vec![vec![sample]],
        &["sample"],
        &per_seq,
        true,
        None,
    );
    let per_seq = std::fs::read_to_string(per_seq).unwrap();
    for (id, outcome) in [
        ("classified", "classified"),
        ("ambiguous", "ambiguous"),
        ("unclassified", "unclassified"),
    ] {
        let row = per_seq
            .lines()
            .find(|line| line.split('\t').nth(1) == Some(id))
            .unwrap();
        assert_eq!(row.split('\t').nth(2), Some(outcome));
    }
}

#[test]
fn test_build_classify_accepts_bare_fastx_file() {
    let dir = TempDir::new().unwrap();
    let (targets, index) = (dir.path().join("refs.fa"), dir.path().join("i.sk"));
    write_fasta(&targets, "r1", SEQ_A);

    skope::run_build_classify(&BuildClassifyConfig {
        targets_path: targets,
        individual: false,
        kmer_length: 15,
        smer_length: 7,
        complexity: 0.0,
        threads: 1,
        output_path: Some(index.clone()),
        quiet: true,
    })
    .unwrap();

    assert_eq!(
        skope::read_index_kind(&index),
        Some(skope::IndexKind::Classify)
    );
}

#[test]
fn test_wrong_index_kind_names_both_kinds() {
    let dir = TempDir::new().unwrap();
    let (targets, query_index, classify_index) = (
        dir.path().join("t.fa"),
        dir.path().join("q.sk"),
        dir.path().join("c.sk"),
    );
    write_fasta(&targets, "t1", SEQ_A);
    build_index(targets.clone(), vec![], query_index.clone());
    skope::run_build_classify(&BuildClassifyConfig {
        targets_path: targets,
        individual: false,
        kmer_length: 15,
        smer_length: 7,
        complexity: 0.0,
        threads: 1,
        output_path: Some(classify_index.clone()),
        quiet: true,
    })
    .unwrap();

    let err = skope::resolve_targets(
        &query_index,
        skope::IndexKind::Classify,
        skope::StdinTargets::Reject,
    )
    .unwrap_err()
    .to_string();
    assert!(
        err.contains("is a skope query index, not a classification index"),
        "got: {err}"
    );

    let err = skope::resolve_targets(
        &classify_index,
        skope::IndexKind::Query,
        skope::StdinTargets::Reject,
    )
    .unwrap_err()
    .to_string();
    assert!(
        err.contains("is a skope classification index, not a query index"),
        "got: {err}"
    );
}

#[test]
fn test_resolve_targets_input_forms() {
    let dir = TempDir::new().unwrap();
    let file = dir.path().join("t.fa");
    write_fasta(&file, "t1", SEQ_A);

    // Single fastx file grouped by filename
    match skope::resolve_targets(&file, skope::IndexKind::Query, skope::StdinTargets::Reject)
        .unwrap()
    {
        skope::TargetSource::File(group) => {
            assert_eq!(group.name, "t");
            assert_eq!(group.files, vec![file.clone()]);
        }
        other => panic!("expected File, got {other:?}"),
    }

    // Directory grouped by child
    let groups_dir = TempDir::new().unwrap();
    write_fasta(&groups_dir.path().join("a.fa"), "a1", SEQ_A);
    write_fasta(&groups_dir.path().join("b.fa"), "b1", SEQ_B);
    match skope::resolve_targets(
        groups_dir.path(),
        skope::IndexKind::Query,
        skope::StdinTargets::Reject,
    )
    .unwrap()
    {
        skope::TargetSource::Directory(groups) => {
            let names: Vec<&str> = groups.iter().map(|g| g.name.as_str()).collect();
            assert_eq!(names, vec!["a", "b"]);
        }
        other => panic!("expected Directory, got {other:?}"),
    }

    // Prebuilt index
    let index = dir.path().join("i.sk");
    build_index(file, vec![], index.clone());
    assert!(matches!(
        skope::resolve_targets(&index, skope::IndexKind::Query, skope::StdinTargets::Reject)
            .unwrap(),
        skope::TargetSource::Index(_)
    ));

    // stdin rejected when samples claim `-`, accepted by build commands
    let dash = std::path::Path::new("-");
    let err = skope::resolve_targets(dash, skope::IndexKind::Query, skope::StdinTargets::Reject)
        .unwrap_err()
        .to_string();
    assert!(err.contains("cannot be read from stdin"), "got: {err}");

    match skope::resolve_targets(dash, skope::IndexKind::Query, skope::StdinTargets::Accept)
        .unwrap()
    {
        skope::TargetSource::File(group) => assert_eq!(group.name, "stdin"),
        other => panic!("expected File, got {other:?}"),
    }

    // Neither fastx nor index
    let junk = dir.path().join("notes.txt");
    std::fs::write(&junk, "not a sequence file at all\n").unwrap();
    let err = skope::resolve_targets(
        &junk,
        skope::IndexKind::Classify,
        skope::StdinTargets::Reject,
    )
    .unwrap_err()
    .to_string();
    assert!(err.contains("is not a fastx file"), "got: {err}");

    // Missing path
    let err = skope::resolve_targets(
        &dir.path().join("nope.fa"),
        skope::IndexKind::Query,
        skope::StdinTargets::Reject,
    )
    .unwrap_err()
    .to_string();
    assert!(err.contains("does not exist"), "got: {err}");
}

fn build_index_complexity(targets: PathBuf, out: PathBuf, complexity: f32) {
    skope::run_build_query(&skope::BuildQueryConfig {
        targets_path: targets,
        background_paths: vec![],
        kmer_length: 31,
        smer_length: 9,
        complexity,
        individual: false,
        positions: false,
        threads: 1,
        output_path: Some(out),
        quiet: true,
        fraction: 1.0,
    })
    .unwrap()
}

fn query_complexity(
    targets_path: PathBuf,
    out: &std::path::Path,
    complexity: f32,
) -> anyhow::Result<()> {
    skope::run_query(&ContainmentConfig {
        background_paths: Vec::new(),
        targets_path,
        sample_paths: vec![vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")]],
        sample_names: vec!["s".to_string()],
        kmer_length: 31,
        smer_length: 9,
        complexity,
        threads: 1,
        output_path: Some(out.to_path_buf()),
        quiet: true,
        abundance_thresholds: vec![1],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_kmers_path: None,
        confidence: false,
        fraction: 1.0,
        no_total: true,
        individual: false,
    })
}

/// Sum of the target_kmers column, the kmers surviving filtering
fn target_kmers(tsv: &std::path::Path) -> u64 {
    let text = std::fs::read_to_string(tsv).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    let idx = lines[0]
        .split('\t')
        .position(|col| col == "target_kmers")
        .expect("target_kmers column not found");
    lines[1..]
        .iter()
        .map(|line| line.split('\t').nth(idx).unwrap().parse::<u64>().unwrap())
        .sum()
}

// Index-time and query-time filtering must give identical containment
#[test]
fn test_complexity_index_matches_fastx() {
    let dir = TempDir::new().unwrap();
    let (target, index) = (
        PathBuf::from("data/zmrp21.viruses.fa"),
        dir.path().join("t90.sk"),
    );
    build_index_complexity(target.clone(), index.clone(), 0.9);

    // Header stores the threshold
    // An f32 threshold round-trips through the header bit-identically
    let (_k, _s, _frac, complexity) = skope::read_query_index_meta(&index).unwrap();
    assert_eq!(complexity, 0.9, "stored complexity {complexity} != 0.9");

    let (i, f) = (dir.path().join("i.tsv"), dir.path().join("f.tsv"));
    query_complexity(index, &i, 0.0).unwrap(); // Index threshold wins, CLI default ignored
    query_complexity(target, &f, 0.9).unwrap();
    assert_eq!(
        std::fs::read_to_string(i).unwrap(),
        std::fs::read_to_string(f).unwrap()
    );
}

#[test]
fn test_complexity_drops_low_complexity_kmers() {
    let dir = TempDir::new().unwrap();
    let target = PathBuf::from("data/zmrp21.viruses.fa");
    let (off, on) = (dir.path().join("off.tsv"), dir.path().join("on.tsv"));

    query_complexity(target.clone(), &off, 0.0).unwrap();
    query_complexity(target, &on, 0.95).unwrap();

    let (unfiltered, filtered) = (target_kmers(&off), target_kmers(&on));
    assert!(unfiltered > 0);
    assert!(
        filtered < unfiltered,
        "filtering kept {filtered} of {unfiltered} kmers"
    );
}

#[test]
fn test_complexity_conflicting_with_index_is_rejected() {
    let dir = TempDir::new().unwrap();
    let index = dir.path().join("t90.sk");
    build_index_complexity(PathBuf::from("data/zmrp21.viruses.fa"), index.clone(), 0.9);

    let out = dir.path().join("out.tsv");
    let err = query_complexity(index, &out, 0.5).unwrap_err();
    assert!(
        err.to_string()
            .contains("conflicts with the index's complexity"),
        "unexpected error: {err}"
    );
}

// ── All-kmers mode (s = 0, --all-kmers) ─────────────────────────────────────

fn query_all_kmers(
    targets_path: PathBuf,
    out: &std::path::Path,
    smer_length: u8,
) -> anyhow::Result<()> {
    skope::run_query(&ContainmentConfig {
        background_paths: Vec::new(),
        targets_path,
        sample_paths: vec![vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")]],
        sample_names: vec!["s".to_string()],
        kmer_length: 31,
        smer_length,
        complexity: 0.0,
        threads: 1,
        output_path: Some(out.to_path_buf()),
        quiet: true,
        abundance_thresholds: vec![1],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_kmers_path: None,
        confidence: false,
        fraction: 1.0,
        no_total: true,
        individual: true,
    })
}

/// Write one record per k-mer, the case syncmer selection cannot serve
fn write_kmer_records(path: &std::path::Path, k: usize, n: usize) {
    use std::io::Write;
    let source = format!("{DISC_COMMON}{DISC_UNIQUE_A}{DISC_UNIQUE_B}");
    let mut f = std::fs::File::create(path).unwrap();
    for i in 0..n {
        writeln!(f, ">km{i}\n{}", &source[i..i + k]).unwrap();
    }
}

// A fastx of bare k-mers keeps every record in all-kmers mode, but syncmer
// selection retains only ~1/w of them
#[test]
fn test_all_kmers_keeps_every_kmer_record() {
    let dir = TempDir::new().unwrap();
    let targets = dir.path().join("kmers.fa");
    write_kmer_records(&targets, 31, 200);

    let all = dir.path().join("all.tsv");
    query_all_kmers(targets.clone(), &all, 0).unwrap();
    let syncmers = dir.path().join("sync.tsv");
    query_all_kmers(targets, &syncmers, 9).unwrap();

    let all_n = target_kmers(&all);
    let sync_n = target_kmers(&syncmers);
    assert_eq!(all_n, 200, "every k-mer record should be retained");
    assert!(
        sync_n * 5 < all_n,
        "syncmer selection should keep far fewer: {sync_n} vs {all_n}"
    );
}

// A query index built with --all-kmers must reproduce the direct fastx run
#[test]
fn test_all_kmers_index_matches_fastx() {
    let dir = TempDir::new().unwrap();
    let targets = dir.path().join("t.fa");
    write_fasta(&targets, "t", &format!("{DISC_COMMON}{DISC_UNIQUE_A}"));

    let index = dir.path().join("t.sk");
    skope::run_build_query(&skope::BuildQueryConfig {
        targets_path: targets.clone(),
        background_paths: vec![],
        kmer_length: 31,
        smer_length: 0,
        complexity: 0.0,
        individual: true,
        positions: false,
        threads: 1,
        output_path: Some(index.clone()),
        quiet: true,
        fraction: 1.0,
    })
    .unwrap();

    let from_fastx = dir.path().join("fastx.tsv");
    query_all_kmers(targets, &from_fastx, 0).unwrap();
    let from_index = dir.path().join("index.tsv");
    query_all_kmers(index, &from_index, 0).unwrap();

    assert_eq!(
        std::fs::read_to_string(&from_fastx).unwrap(),
        std::fs::read_to_string(&from_index).unwrap()
    );
}

#[test]
fn test_query_individual_warns_on_empty_targets() {
    let dir = TempDir::new().unwrap();
    let targets = dir.path().join("targets.fa");
    let sample = dir.path().join("sample.fa");
    std::fs::write(
        &targets,
        format!(">long\n{SEQ_A}\n>tooshort\nACGT\n>short2\nACGT\n>short3\nACGT\n>short4\nACGT\n"),
    )
    .unwrap();
    write_fasta(&sample, "sample", SEQ_A);

    let run = |args: &[&str], out: &str| {
        std::process::Command::new(env!("CARGO_BIN_EXE_skope"))
            .args(args)
            .arg("-o")
            .arg(dir.path().join(out))
            .arg(&targets)
            .arg(&sample)
            .output()
            .unwrap()
    };

    let loud = run(&["query", "-i"], "loud.tsv");
    assert!(loud.status.success());
    let stderr = String::from_utf8_lossy(&loud.stderr).to_string();
    assert!(
        stderr.contains(
            "Warning: 4 of 5 targets yielded no k-mers: tooshort, short2, short3 and 1 more"
        ),
        "{stderr}"
    );

    let quiet = run(&["query", "-i", "-q"], "quiet.tsv");
    assert!(quiet.status.success());
    assert!(!String::from_utf8_lossy(&quiet.stderr).contains("Warning"));
}

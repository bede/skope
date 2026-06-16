use skope::{
    BuildConfig, ClassifyConfig, ContainmentConfig, LengthHistogramConfig, SortOrder,
    discover_sequence_groups,
};
use std::path::PathBuf;
use tempfile::{NamedTempFile, TempDir};

#[test]
fn test_multisample_processing() {
    let config = ContainmentConfig {
        targets_path: PathBuf::from("data/zmrp21.viruses.fa"),
        sample_paths: vec![
            vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")],
            vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")],
        ],
        sample_names: vec!["sample1".to_string(), "sample2".to_string()],
        kmer_length: 31,
        smer_length: 15,
        threads: 2,
        output_path: None,
        quiet: true,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_syncmers_path: None,
        confidence: false,
        no_total: false,
        spacing: 1,
        individual: true,
    };

    assert!(skope::run_containment_analysis(&config).is_ok());
}

#[test]
fn test_multisample_report_structure() {
    let temp_output = NamedTempFile::new().unwrap();
    let output_path = temp_output.path().to_path_buf();

    let config = ContainmentConfig {
        targets_path: PathBuf::from("data/zmrp21.viruses.fa"),
        sample_paths: vec![
            vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")],
            vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")],
        ],
        sample_names: vec!["sample_a".to_string(), "sample_b".to_string()],
        kmer_length: 31,
        smer_length: 15,
        threads: 2,
        output_path: Some(output_path.clone()),
        quiet: true,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_syncmers_path: None,
        confidence: false,
        no_total: false,
        spacing: 1,
        individual: true,
    };

    skope::run_containment_analysis(&config).unwrap();

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
        targets_path: PathBuf::from("data/zmrp21.viruses.fa"),
        sample_paths: vec![vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")]],
        sample_names: vec!["sample".to_string()],
        kmer_length: 31,
        smer_length: 15,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_syncmers_path: None,
        confidence: true,
        no_total: false,
        spacing: 31,
        individual: true,
    };

    skope::run_containment_analysis(&config).unwrap();
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
    // ani_est is "-" when suppressed, otherwise a containment ANI in [0, 1].
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
        targets_path: PathBuf::from("data/zmrp21.viruses.fa"),
        sample_paths: vec![vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")]],
        sample_names: vec!["test".to_string()],
        kmer_length: 31,
        smer_length: 15,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Target,
        dump_syncmers_path: None,
        confidence: false,
        no_total: false,
        spacing: 1,
        individual: true,
    };

    skope::run_containment_analysis(&config).unwrap();
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
        targets_path: PathBuf::from("data/zmrp21.viruses.fa"),
        sample_paths: vec![vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")]],
        sample_names: vec!["test".to_string()],
        kmer_length: 31,
        smer_length: 15,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Containment,
        dump_syncmers_path: None,
        confidence: false,
        no_total: false,
        spacing: 1,
        individual: true,
    };

    skope::run_containment_analysis(&config).unwrap();
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

    // With non-overlapping syncmers, absolute containment is lower than before,
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
    let groups_dir = TempDir::new().unwrap();
    std::fs::copy(
        "data/zmrp21.viruses.fa",
        groups_dir.path().join("viruses.fa"),
    )
    .unwrap();

    let temp_output = NamedTempFile::new().unwrap();

    let config = LengthHistogramConfig {
        index_path: groups_dir.path().to_path_buf(),
        sample_paths: vec![vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")]],
        sample_names: vec!["test".to_string()],
        kmer_length: 31,
        smer_length: 15,
        min_hits: 1,
        min_fraction: 0.0,
        discriminatory: false,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        limit_bp: None,
        include_all_seqs: false,
    };

    skope::run_length_histogram_analysis(&config).unwrap();

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
        index_path: PathBuf::from("-"),
        sample_paths: vec![vec![PathBuf::from("data/rsviruses17900.1k.fastq.zst")]],
        sample_names: vec!["test".to_string()],
        kmer_length: 31,
        smer_length: 15,
        min_hits: 1,
        min_fraction: 0.0,
        discriminatory: false,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        limit_bp: None,
        include_all_seqs: true,
    };

    skope::run_length_histogram_analysis(&config).unwrap();

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

const SEQ_A: &str =
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
const SEQ_B: &str =
    "GTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCA";

#[test]
fn test_discover_sequence_groups_mixed_layout() {
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

    let groups = discover_sequence_groups(root).unwrap();
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
fn test_discover_sequence_groups_nested_subdir_errors() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    std::fs::create_dir_all(root.join("class_a/nested")).unwrap();
    write_fasta(&root.join("class_a/part1.fa"), "a", SEQ_A);

    let err = discover_sequence_groups(root).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("Nested subdirectory"), "got: {}", msg);
}

#[test]
fn test_discover_sequence_groups_empty_subdir_errors() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    std::fs::create_dir(root.join("class_a")).unwrap();
    std::fs::write(root.join("class_a/notes.txt"), b"no fastx here").unwrap();

    let err = discover_sequence_groups(root).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("no fastx files"), "got: {}", msg);
    assert!(msg.contains("class_a"), "got: {}", msg);
}

#[test]
fn test_discover_sequence_groups_duplicate_name_errors() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    write_fasta(&root.join("foo.fa"), "f", SEQ_A);
    std::fs::create_dir(root.join("foo")).unwrap();
    write_fasta(&root.join("foo/inner.fa"), "i", SEQ_B);

    let err = discover_sequence_groups(root).unwrap_err();
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
        targets_path: root.to_path_buf(),
        sample_paths: vec![vec![sample.path().to_path_buf()]],
        sample_names: vec!["s".to_string()],
        kmer_length: 15,
        smer_length: 7,
        threads: 1,
        output_path: Some(output.path().to_path_buf()),
        quiet: true,
        abundance_thresholds: vec![1],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_syncmers_path: None,
        confidence: false,
        no_total: true,
        spacing: 1,
        individual: false,
    };
    skope::run_containment_analysis(&config).unwrap();

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
// region yields cross-target syncmers; the unique tails do not.
const DISC_COMMON: &str = "GATTACAGGCATCCTAGCTAGGACTTGCAACATGCTTAGCCATGGAACTGTCCAGTTACGGATCCTAGGCATTAGCCAGTTCATGGACTTAGCGGATCCTA";
const DISC_UNIQUE_A: &str = "TTGCAACGGTACCATTAGCGGATCCTTAGCAACATGCTTAGCCATGGAACTGTCCAGTTACGGATCCTAGGCATTAGCCAGTTCATGGACTTAGCGGATCC";
const DISC_UNIQUE_B: &str = "CCAGTTACGGATCCTAGGCATTAGCCAGTTCATGGACTTAGCGGATCCTAGCTAGGACTTGCAACATGCTTAGCCATGGAACTGTCCAGTTACGGATCCTA";

// Parse a --dump-syncmers TSV into target -> set of k-mers (col 1 -> col 3).
fn parse_dump(path: &std::path::Path) -> std::collections::HashMap<String, std::collections::HashSet<String>> {
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
fn test_dump_syncmers_respects_discriminatory() {
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
        targets_path: root.to_path_buf(),
        sample_paths: vec![vec![sample.path().to_path_buf()]],
        sample_names: vec!["s".to_string()],
        kmer_length: 15,
        smer_length: 7,
        threads: 1,
        output_path: None,
        quiet: true,
        abundance_thresholds: vec![1],
        discriminatory,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_syncmers_path: Some(dump),
        confidence: false,
        no_total: true,
        spacing: 1,
        individual: false,
    };

    let plain_dump = NamedTempFile::new().unwrap();
    skope::run_containment_analysis(&make_config(false, plain_dump.path().to_path_buf())).unwrap();
    let plain = parse_dump(plain_dump.path());

    let disc_dump = NamedTempFile::new().unwrap();
    skope::run_containment_analysis(&make_config(true, disc_dump.path().to_path_buf())).unwrap();
    let disc = parse_dump(disc_dump.path());

    let plain_a = &plain["tA"];
    let plain_b = &plain["tB"];
    let shared: Vec<_> = plain_a.intersection(plain_b).collect();
    assert!(
        !shared.is_empty(),
        "expected cross-target shared syncmers in plain dump"
    );

    let disc_a = &disc["tA"];
    let disc_b = &disc["tB"];
    let disc_shared: Vec<_> = disc_a.intersection(disc_b).collect();
    assert!(
        disc_shared.is_empty(),
        "discriminatory dump must contain no cross-target shared syncmers, found: {:?}",
        disc_shared
    );

    let plain_total = plain_a.len() + plain_b.len();
    let disc_total = disc_a.len() + disc_b.len();
    assert!(
        disc_total < plain_total,
        "discriminatory dump ({}) should have fewer syncmers than plain ({})",
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
    let config = BuildConfig {
        groups_dir: root.to_path_buf(),
        kmer_length: 15,
        smer_length: 7,
        threads: 1,
        output_path: Some(idx_out.path().to_path_buf()),
        quiet: true,
    };
    skope::build_classification_index(&config).unwrap();

    // Index file should exist and be non-empty
    let meta = std::fs::metadata(idx_out.path()).unwrap();
    assert!(meta.len() > 0);

    // Run classify against it with a sample whose seq matches class_b
    let sample = NamedTempFile::new().unwrap();
    write_fasta(sample.path(), "s1", SEQ_B);

    let classify_out = NamedTempFile::new().unwrap();
    let cfg = ClassifyConfig {
        index_path: idx_out.path().to_path_buf(),
        sample_paths: vec![vec![sample.path().to_path_buf()]],
        sample_names: vec!["s".to_string()],
        kmer_length: 15,
        smer_length: 7,
        min_hits: 1,
        min_fraction: 0.0,
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
    let config = BuildConfig {
        groups_dir: root.to_path_buf(),
        kmer_length: 15,
        smer_length: 7,
        threads: 1,
        output_path: Some(idx_out.path().to_path_buf()),
        quiet: true,
    };
    let err = skope::build_classification_index(&config).unwrap_err();
    let msg = format!("{:#}", err);
    assert!(msg.contains("Too many groups"), "got: {}", msg);
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

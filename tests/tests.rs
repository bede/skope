use grate::{ContainmentConfig, LengthHistogramConfig, OutputFormat, SortOrder};
use std::path::PathBuf;
use tempfile::NamedTempFile;

#[test]
fn test_multisample_processing() {
    let config = ContainmentConfig {
        targets_path: PathBuf::from("data/zmrp21.combined-segments.fa"),
        sample_paths: vec![
            vec![PathBuf::from("data/rsviruses17900.10k.fastq.zst")],
            vec![PathBuf::from("data/rsviruses17900.10k.fastq.zst")],
        ],
        sample_names: vec!["sample1".to_string(), "sample2".to_string()],
        kmer_length: 31,
        smer_length: 15,
        threads: 2,
        output_path: None,
        quiet: true,
        output_format: OutputFormat::Tsv,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_positions_path: None,
        no_total: false,
    };

    assert!(grate::run_containment_analysis(&config).is_ok());
}

#[test]
fn test_multisample_report_structure() {
    let temp_output = NamedTempFile::new().unwrap();
    let output_path = temp_output.path().to_path_buf();

    let config = ContainmentConfig {
        targets_path: PathBuf::from("data/zmrp21.combined-segments.fa"),
        sample_paths: vec![
            vec![PathBuf::from("data/rsviruses17900.10k.fastq.zst")],
            vec![PathBuf::from("data/rsviruses17900.10k.fastq.zst")],
        ],
        sample_names: vec!["sample_a".to_string(), "sample_b".to_string()],
        kmer_length: 31,
        smer_length: 15,
        threads: 2,
        output_path: Some(output_path.clone()),
        quiet: true,
        output_format: OutputFormat::Tsv,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
        dump_positions_path: None,
        no_total: false,
    };

    grate::run_containment_analysis(&config).unwrap();

    let tsv_str = std::fs::read_to_string(&output_path).unwrap();
    let lines: Vec<&str> = tsv_str.lines().collect();

    // Check header
    assert!(lines[0].starts_with("target\tsample\t"));

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
}

#[test]
fn test_sort_target() {
    let temp_output = NamedTempFile::new().unwrap();
    let config = ContainmentConfig {
        targets_path: PathBuf::from("data/zmrp21.combined-segments.fa"),
        sample_paths: vec![vec![PathBuf::from("data/rsviruses17900.10k.fastq.zst")]],
        sample_names: vec!["test".to_string()],
        kmer_length: 31,
        smer_length: 15,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        output_format: OutputFormat::Tsv,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Target,
        dump_positions_path: None,
        no_total: false,
    };

    grate::run_containment_analysis(&config).unwrap();
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
        targets_path: PathBuf::from("data/zmrp21.combined-segments.fa"),
        sample_paths: vec![vec![PathBuf::from("data/rsviruses17900.10k.fastq.zst")]],
        sample_names: vec!["test".to_string()],
        kmer_length: 31,
        smer_length: 15,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        output_format: OutputFormat::Tsv,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Containment,
        dump_positions_path: None,
        no_total: false,
    };

    grate::run_containment_analysis(&config).unwrap();
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

    // Highest containment should be > 60%
    assert!(
        containments[0] > 0.6,
        "Highest containment {} should be > 0.6",
        containments[0]
    );
}

#[test]
fn test_length_histogram() {
    let temp_output = NamedTempFile::new().unwrap();

    let config = LengthHistogramConfig {
        targets_path: PathBuf::from("data/zmrp21.combined-segments.fa"),
        sample_paths: vec![vec![PathBuf::from("data/rsviruses17900.10k.fastq.zst")]],
        sample_names: vec!["test".to_string()],
        kmer_length: 31,
        smer_length: 15,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        limit_bp: None,
        include_all_seqs: false,
    };

    grate::run_length_histogram_analysis(&config).unwrap();

    // Verify TSV output has header and data rows
    let content = std::fs::read_to_string(temp_output.path()).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert!(
        lines[0].starts_with("sample\tlength\tcount\t"),
        "Header should be tab-separated"
    );
    assert!(lines.len() > 1, "Should have data rows");
    assert!(
        lines[1].starts_with("test\t"),
        "Data rows should start with sample name"
    );
}

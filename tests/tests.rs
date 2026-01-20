use grate::{ContainmentConfig, LengthHistogramConfig, OutputFormat, SortOrder};
use std::path::PathBuf;
use tempfile::NamedTempFile;

#[test]
fn test_multisample_processing() {
    let config = ContainmentConfig {
        targets_path: PathBuf::from("data/zmrp21.combined-segments.fa"),
        reads_paths: vec![
            vec![PathBuf::from("data/rsviruses17900.10k.fastq.zst")],
            vec![PathBuf::from("data/rsviruses17900.10k.fastq.zst")],
        ],
        sample_names: vec!["sample1".to_string(), "sample2".to_string()],
        kmer_length: 31,
        window_size: 15,
        threads: 2,
        output_path: None,
        quiet: true,
        output_format: OutputFormat::Json,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
    };

    assert!(grate::run_containment_analysis(&config).is_ok());
}

#[test]
fn test_multisample_report_structure() {
    let temp_output = NamedTempFile::new().unwrap();
    let output_path = temp_output.path().to_path_buf();

    let config = ContainmentConfig {
        targets_path: PathBuf::from("data/zmrp21.combined-segments.fa"),
        reads_paths: vec![
            vec![PathBuf::from("data/rsviruses17900.10k.fastq.zst")],
            vec![PathBuf::from("data/rsviruses17900.10k.fastq.zst")],
        ],
        sample_names: vec!["sample_a".to_string(), "sample_b".to_string()],
        kmer_length: 31,
        window_size: 15,
        threads: 2,
        output_path: Some(output_path.clone()),
        quiet: true,
        output_format: OutputFormat::Json,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Original,
    };

    grate::run_containment_analysis(&config).unwrap();

    let json_str = std::fs::read_to_string(&output_path).unwrap();
    let report: grate::Report = serde_json::from_str(&json_str).unwrap();

    assert_eq!(report.samples.len(), 2);
    assert_eq!(report.samples[0].sample_name, "sample_a");
    assert_eq!(report.samples[1].sample_name, "sample_b");
    assert_eq!(report.samples[0].targets.len(), 16);
    assert_eq!(report.samples[1].targets.len(), 16);
}

#[test]
fn test_sort_target() {
    let temp_output = NamedTempFile::new().unwrap();
    let config = ContainmentConfig {
        targets_path: PathBuf::from("data/zmrp21.combined-segments.fa"),
        reads_paths: vec![vec![PathBuf::from("data/rsviruses17900.10k.fastq.zst")]],
        sample_names: vec!["test".to_string()],
        kmer_length: 31,
        window_size: 15,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        output_format: OutputFormat::Json,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Target,
    };

    grate::run_containment_analysis(&config).unwrap();
    let report: grate::Report =
        serde_json::from_str(&std::fs::read_to_string(temp_output.path()).unwrap()).unwrap();
    let targets = &report.samples[0].targets;

    // Verify sorted alphabetically by target
    for i in 0..targets.len() - 1 {
        assert!(targets[i].target <= targets[i + 1].target);
    }
}

#[test]
fn test_sort_containment() {
    let temp_output = NamedTempFile::new().unwrap();
    let config = ContainmentConfig {
        targets_path: PathBuf::from("data/zmrp21.combined-segments.fa"),
        reads_paths: vec![vec![PathBuf::from("data/rsviruses17900.10k.fastq.zst")]],
        sample_names: vec!["test".to_string()],
        kmer_length: 31,
        window_size: 15,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        output_format: OutputFormat::Json,
        abundance_thresholds: vec![10],
        discriminatory: false,
        limit_bp: None,
        sort_order: SortOrder::Containment,
    };

    grate::run_containment_analysis(&config).unwrap();
    let report: grate::Report =
        serde_json::from_str(&std::fs::read_to_string(temp_output.path()).unwrap()).unwrap();
    let targets = &report.samples[0].targets;

    // Verify sorted by containment descending
    for i in 0..targets.len() - 1 {
        assert!(targets[i].containment1 >= targets[i + 1].containment1);
    }
    // Highest containment should be > 60%
    assert!(targets[0].containment1 > 0.6);
}

#[test]
fn test_length_histogram() {
    let temp_output = NamedTempFile::new().unwrap();

    let config = LengthHistogramConfig {
        targets_path: PathBuf::from("data/zmrp21.combined-segments.fa"),
        reads_paths: vec![vec![PathBuf::from("data/rsviruses17900.10k.fastq.zst")]],
        sample_names: vec!["test".to_string()],
        kmer_length: 31,
        window_size: 21,
        threads: 2,
        output_path: Some(temp_output.path().to_path_buf()),
        quiet: true,
        limit_bp: None,
        include_all_reads: false,
        discriminatory: false,
    };

    grate::run_length_histogram_analysis(&config).unwrap();

    let report: grate::LengthHistogramReport =
        serde_json::from_str(&std::fs::read_to_string(temp_output.path()).unwrap()).unwrap();

    assert_eq!(report.samples.len(), 1);
    assert!(report.samples[0].total_stats.reads_with_hits > 0);
    assert!(!report.samples[0].total_stats.length_histogram.is_empty());
}

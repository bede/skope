//! # Grate
//!
//! A fast minimizer-based containment analysis tool for genomic sequences.
//!
//! Grate analyzes the containment of minimizers from reference sequences in read datasets,
//! providing detailed statistics on containment and abundance per target sequence.
//!

pub mod containment;
pub mod length;
pub mod minimizers;

// Re-export the main functionality
pub use containment::{
    run_containment_analysis, ContainmentConfig, ContainmentParameters, ContainmentResult,
    OutputFormat, Report, SampleResults, SortOrder, TargetInfo, TimingStats, TotalStats,
};

pub use length::{
    run_length_histogram_analysis, LengthHistogramConfig, LengthHistogramParameters,
    LengthHistogramReport, SampleLengthResults, TargetLengthResult, TotalLengthStats,
};

pub use minimizers::{
    decode_u128, decode_u64, fill_minimizers, fill_minimizers_with_positions, Buffers, KmerHasher,
    MinimizerVec, DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE,
};

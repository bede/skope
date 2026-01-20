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
    ContainmentConfig, ContainmentParameters, ContainmentResult, OutputFormat, Report,
    SampleResults, SortOrder, TargetInfo, TimingStats, TotalStats, run_containment_analysis,
};

pub use length::{
    LengthHistogramConfig, LengthHistogramParameters, LengthHistogramReport, LengthHistogramResult,
    run_length_histogram_analysis,
};

pub use minimizers::{
    Buffers, DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE, KmerHasher, MinimizerVec, decode_u64,
    decode_u128, fill_minimizers, fill_minimizers_with_positions,
};

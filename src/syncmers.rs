use crate::FixedRapidHasher;
use packed_seq::{PackedNSeqVec, SeqVec, unpack_base};
use std::hash::{BuildHasher, Hasher};

pub const DEFAULT_KMER_LENGTH: u8 = 31;
pub const DEFAULT_SMER_LENGTH: u8 = 9;

pub type KmerHasher = simd_minimizers::seq_hash::NtHasher<true, 1>;

/// FracMinHash: keep syncmer if mix(hash(kmer)) lte [0,1] threshold
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FracMinHash {
    threshold: u64,
}

impl FracMinHash {
    /// Keep all syncmers
    pub const NONE: FracMinHash = FracMinHash {
        threshold: u64::MAX,
    };

    pub fn from_fraction(fraction: f64) -> FracMinHash {
        if fraction >= 1.0 || fraction.is_nan() {
            FracMinHash::NONE
        } else if fraction <= 0.0 {
            FracMinHash { threshold: 0 }
        } else {
            // Scale by 2^64 in f64 space since 2^64 is unrepresentable as u64
            FracMinHash {
                threshold: (fraction * 2.0f64.powi(64)) as u64,
            }
        }
    }

    /// Whether selection is a no-op (retains everything)
    #[inline(always)]
    pub fn is_none(&self) -> bool {
        self.threshold == u64::MAX
    }

    #[inline(always)]
    fn keeps_u64(&self, kmer: u64) -> bool {
        rapid_mix(&kmer.to_le_bytes()) <= self.threshold
    }

    #[inline(always)]
    fn keeps_u128(&self, kmer: u128) -> bool {
        rapid_mix(&kmer.to_le_bytes()) <= self.threshold
    }

    pub fn retain(&self, syncmers: &mut SyncmerVec, positions: Option<&mut Vec<usize>>) {
        if self.is_none() {
            return;
        }
        match syncmers {
            SyncmerVec::U64(vec) => match positions {
                Some(pos) => retain_paired(vec, pos, |&v| self.keeps_u64(v)),
                None => vec.retain(|&v| self.keeps_u64(v)),
            },
            SyncmerVec::U128(vec) => match positions {
                Some(pos) => retain_paired(vec, pos, |&v| self.keeps_u128(v)),
                None => vec.retain(|&v| self.keeps_u128(v)),
            },
        }
    }
}

/// Retain `values` (and aligned `positions`) where `keep` holds
fn retain_paired<T: Copy>(
    values: &mut Vec<T>,
    positions: &mut Vec<usize>,
    keep: impl Fn(&T) -> bool,
) {
    debug_assert_eq!(values.len(), positions.len());
    let mut w = 0;
    for r in 0..values.len() {
        if keep(&values[r]) {
            values[w] = values[r];
            positions[w] = positions[r];
            w += 1;
        }
    }
    values.truncate(w);
    positions.truncate(w);
}

#[inline(always)]
fn rapid_mix(bytes: &[u8]) -> u64 {
    let mut h = FixedRapidHasher.build_hasher();
    h.write(bytes);
    h.finish()
}

/// Zero-cost abstraction over u64 and u128 syncmer vectors
#[derive(Debug, Clone)]
pub enum SyncmerVec {
    U64(Vec<u64>),
    U128(Vec<u128>),
}

impl SyncmerVec {
    pub fn clear(&mut self) {
        match self {
            SyncmerVec::U64(v) => v.clear(),
            SyncmerVec::U128(v) => v.clear(),
        }
    }

    pub fn len(&self) -> usize {
        match self {
            SyncmerVec::U64(v) => v.len(),
            SyncmerVec::U128(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match self {
            SyncmerVec::U64(v) => v.is_empty(),
            SyncmerVec::U128(v) => v.is_empty(),
        }
    }
}

/// Decode u64 syncmer (2-bit canonical k-mer)
pub fn decode_u64(syncmer: u64, k: u8) -> Vec<u8> {
    (0..k)
        .map(|i| {
            let base_bits = ((syncmer >> (2 * i)) & 0b11) as u8;
            unpack_base(base_bits)
        })
        .rev()
        .collect()
}

/// Decode u128 syncmer (2-bit canonical k-mer)
pub fn decode_u128(syncmer: u128, k: u8) -> Vec<u8> {
    (0..k)
        .map(|i| {
            let base_bits = ((syncmer >> (2 * i)) & 0b11) as u8;
            unpack_base(base_bits)
        })
        .rev()
        .collect()
}

/// Reusable buffers for syncmer computation
#[derive(Clone)]
pub struct Buffers {
    pub packed_nseq: PackedNSeqVec,
    pub positions: Vec<u32>,
    pub syncmers: SyncmerVec,
}

impl Buffers {
    pub fn new_u64() -> Self {
        Self {
            packed_nseq: PackedNSeqVec {
                seq: Default::default(),
                ambiguous: Default::default(),
            },
            positions: Default::default(),
            syncmers: SyncmerVec::U64(Vec::new()),
        }
    }

    pub fn new_u128() -> Self {
        Self {
            packed_nseq: PackedNSeqVec {
                seq: Default::default(),
                ambiguous: Default::default(),
            },
            positions: Default::default(),
            syncmers: SyncmerVec::U128(Vec::new()),
        }
    }
}

/// Fill syncmers vector and positions vector from sequence
pub fn fill_syncmers_with_positions(
    seq: &[u8],
    hasher: &KmerHasher,
    kmer_length: u8,
    smer_length: u8,
    buffers: &mut Buffers,
    positions_out: &mut Vec<usize>,
) {
    let Buffers {
        packed_nseq,
        positions,
        syncmers,
    } = buffers;

    packed_nseq.seq.clear();
    packed_nseq.ambiguous.clear();
    syncmers.clear();
    positions.clear();
    positions_out.clear();

    if seq.len() < kmer_length as usize {
        return;
    }

    packed_nseq.seq.push_ascii(seq);
    packed_nseq.ambiguous.push_ascii(seq);

    let s = smer_length as usize;
    let w = kmer_length as usize - s + 1;
    let m = simd_minimizers::canonical_open_syncmers(s, w)
        .hasher(hasher)
        .run_skip_ambiguous_windows(packed_nseq.as_slice(), positions);

    match syncmers {
        SyncmerVec::U64(vec) => {
            for (pos, val) in m.pos_and_values_u64() {
                vec.push(val);
                positions_out.push(pos as usize);
            }
        }
        SyncmerVec::U128(vec) => {
            for (pos, val) in m.pos_and_values_u128() {
                vec.push(val);
                positions_out.push(pos as usize);
            }
        }
    }
}

/// Fill syncmers vector from sequence (without positions)
#[inline]
pub fn fill_syncmers(
    seq: &[u8],
    hasher: &KmerHasher,
    kmer_length: u8,
    smer_length: u8,
    buffers: &mut Buffers,
) {
    let Buffers {
        packed_nseq,
        positions,
        syncmers,
    } = buffers;

    packed_nseq.seq.clear();
    packed_nseq.ambiguous.clear();
    syncmers.clear();
    positions.clear();

    if seq.len() < kmer_length as usize {
        return;
    }

    packed_nseq.seq.push_ascii(seq);
    packed_nseq.ambiguous.push_ascii(seq);

    let s = smer_length as usize;
    let w = kmer_length as usize - s + 1;
    let m = simd_minimizers::canonical_open_syncmers(s, w)
        .hasher(hasher)
        .run_skip_ambiguous_windows(packed_nseq.as_slice(), positions);

    match syncmers {
        SyncmerVec::U64(vec) => {
            for (_pos, val) in m.pos_and_values_u64() {
                vec.push(val);
            }
        }
        SyncmerVec::U128(vec) => {
            for (_pos, val) in m.pos_and_values_u128() {
                vec.push(val);
            }
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fill_syncmers() {
        let seq = b"ACGTACGTACGT";
        let k = 5;
        let s = 3;
        let hasher = KmerHasher::new(s as usize);
        let mut buffers = Buffers::new_u64();

        fill_syncmers(seq, &hasher, k, s, &mut buffers);

        // We should have at least one syncmer
        assert!(!buffers.syncmers.is_empty());

        // Test with a sequence shorter than k
        let short_seq = b"ACGT";
        fill_syncmers(short_seq, &hasher, k, s, &mut buffers);
        assert!(buffers.syncmers.is_empty());
    }

    #[test]
    fn test_fill_syncmers_with_positions() {
        let seq = b"ACGTACGTACGTACGT";
        let k = 7;
        let s = 3;
        let hasher = KmerHasher::new(s as usize);
        let mut buffers = Buffers::new_u64();
        let mut positions = Vec::new();

        fill_syncmers_with_positions(seq, &hasher, k, s, &mut buffers, &mut positions);

        // Should have same number of syncmers and positions
        assert_eq!(buffers.syncmers.len(), positions.len());

        // All positions should be valid
        for &pos in &positions {
            assert!(pos + k as usize <= seq.len());
        }
    }

    #[test]
    fn test_syncmers_match_between_apis() {
        let seq = b"ACGTTGCATGTCGCATGATGCATGAGAGCTACGTTGCATGTCGCATGATGCATGAGAGCT";
        let k = 15;
        let s = 7;
        let hasher = KmerHasher::new(s as usize);

        let mut values_only_buffers = Buffers::new_u64();
        fill_syncmers(seq, &hasher, k, s, &mut values_only_buffers);
        let values_only = match &values_only_buffers.syncmers {
            SyncmerVec::U64(v) => v.clone(),
            SyncmerVec::U128(_) => panic!("Expected u64 syncmers for k <= 32"),
        };

        let mut with_pos_buffers = Buffers::new_u64();
        let mut positions = Vec::new();
        fill_syncmers_with_positions(seq, &hasher, k, s, &mut with_pos_buffers, &mut positions);
        let with_pos_values = match &with_pos_buffers.syncmers {
            SyncmerVec::U64(v) => v.clone(),
            SyncmerVec::U128(_) => panic!("Expected u64 syncmers for k <= 32"),
        };

        // Give us same syncmers and same order from both APIs
        assert_eq!(values_only, with_pos_values);
        assert_eq!(with_pos_values.len(), positions.len());
    }

    /// Deterministic pseudo-random DNA of length `n` (LCG), for stable stats tests
    fn pseudo_dna(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed | 1;
        let bases = b"ACGT";
        (0..n)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                bases[((x >> 33) & 0b11) as usize]
            })
            .collect()
    }

    fn syncmers_u64(seq: &[u8], k: u8, s: u8, fmh: FracMinHash) -> Vec<u64> {
        let hasher = KmerHasher::new(s as usize);
        let mut buffers = Buffers::new_u64();
        fill_syncmers(seq, &hasher, k, s, &mut buffers);
        fmh.retain(&mut buffers.syncmers, None);
        match &buffers.syncmers {
            SyncmerVec::U64(v) => v.clone(),
            SyncmerVec::U128(_) => panic!("expected u64"),
        }
    }

    #[test]
    fn test_fmh_none_subset_and_fraction() {
        let seq = pseudo_dna(200_000, 42);
        let all = syncmers_u64(&seq, 31, 9, FracMinHash::NONE);
        // from_fraction(1.0) is a no-op equal to NONE
        assert_eq!(
            all,
            syncmers_u64(&seq, 31, 9, FracMinHash::from_fraction(1.0))
        );
        let kept = syncmers_u64(&seq, 31, 9, FracMinHash::from_fraction(0.1));
        // Deterministic, a strict subset, and roughly the target fraction
        assert_eq!(
            kept,
            syncmers_u64(&seq, 31, 9, FracMinHash::from_fraction(0.1))
        );
        let all_set: std::collections::HashSet<u64> = all.iter().copied().collect();
        assert!(kept.iter().all(|v| all_set.contains(v)));
        let realised = kept.len() as f64 / all.len() as f64;
        assert!(
            (realised - 0.1).abs() < 0.03,
            "realised fraction {realised}"
        );
    }
}

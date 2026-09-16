use crate::FixedRapidHasher;
use packed_seq::{PackedNSeq, PackedNSeqVec, Seq, SeqVec, unpack_base};
use std::hash::BuildHasher;

pub const DEFAULT_KMER_LENGTH: u8 = 31;
pub const DEFAULT_SMER_LENGTH: u8 = 9;

pub type SmerHasher = simd_minimizers::seq_hash::NtHasher<true, 1>;

/// Hasher for syncmer selection. Unused when s is 0, but NtHasher::new(0) underflows
pub fn make_hasher(smer_length: u8) -> SmerHasher {
    SmerHasher::new(smer_length.max(1) as usize)
}

/// FracMinHash: keep k-mer if mix(hash(kmer)) lte [0,1] threshold
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FracMinHash {
    threshold: u64,
}

impl FracMinHash {
    /// Keep all k-mers
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
        rapid_mix_u64(kmer) <= self.threshold
    }

    #[inline(always)]
    fn keeps_u128(&self, kmer: u128) -> bool {
        rapid_mix_u128(kmer) <= self.threshold
    }

    pub fn retain(&self, kmers: &mut KmerVec, positions: Option<&mut Vec<usize>>) {
        if self.is_none() {
            return;
        }
        match kmers {
            KmerVec::U64(vec) => match positions {
                Some(pos) => retain_paired(vec, pos, |&v| self.keeps_u64(v)),
                None => vec.retain(|&v| self.keeps_u64(v)),
            },
            KmerVec::U128(vec) => match positions {
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
fn rapid_mix_u64(kmer: u64) -> u64 {
    FixedRapidHasher.hash_one(kmer)
}

#[inline(always)]
fn rapid_mix_u128(kmer: u128) -> u64 {
    FixedRapidHasher.hash_one(kmer)
}

/// kdust: max-normalised DUST triplet score in [0,1] of a packed canonical k-mer
#[inline]
pub fn calculate_kdust(code: u128, kmer_length: u8) -> f32 {
    // k=3 divides by zero, k<3 has no triplets
    if kmer_length < 4 {
        return 1.0;
    }
    let k = kmer_length as usize;
    let mut counts = [0u8; 64];
    let mut score = 0u32;
    let mut tri = 0usize;
    for i in 0..k {
        tri = ((tri << 2) | ((code >> (2 * i)) & 0b11) as usize) & 0b11_1111;
        if i >= 2 {
            score += counts[tri] as u32;
            counts[tri] += 1;
        }
    }
    let l = (k - 2) as f32;
    1.0 - score as f32 / (l * (l - 1.0) / 2.0)
}

/// Discard k-mers below a kdust threshold in [0, 1]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kdust {
    threshold: f32,
    kmer_length: u8,
}

impl Kdust {
    /// Keep all k-mers
    pub const NONE: Kdust = Kdust {
        threshold: 0.0,
        kmer_length: 0,
    };

    pub fn from_threshold(threshold: f32, kmer_length: u8) -> Kdust {
        if threshold <= 0.0 || threshold.is_nan() {
            Kdust::NONE
        } else {
            Kdust {
                threshold,
                kmer_length,
            }
        }
    }

    /// Whether filtering is a no-op (retains everything)
    #[inline(always)]
    pub fn is_none(&self) -> bool {
        self.threshold <= 0.0
    }

    #[inline(always)]
    fn keeps(&self, kmer: u128) -> bool {
        calculate_kdust(kmer, self.kmer_length) >= self.threshold
    }

    pub fn retain(&self, kmers: &mut KmerVec, positions: Option<&mut Vec<usize>>) {
        if self.is_none() {
            return;
        }
        match kmers {
            KmerVec::U64(vec) => match positions {
                Some(pos) => retain_paired(vec, pos, |&v| self.keeps(v as u128)),
                None => vec.retain(|&v| self.keeps(v as u128)),
            },
            KmerVec::U128(vec) => match positions {
                Some(pos) => retain_paired(vec, pos, |&v| self.keeps(v)),
                None => vec.retain(|&v| self.keeps(v)),
            },
        }
    }
}

/// Zero-cost abstraction over u64 and u128 k-mer vectors
#[derive(Debug, Clone)]
pub enum KmerVec {
    U64(Vec<u64>),
    U128(Vec<u128>),
}

impl KmerVec {
    pub fn clear(&mut self) {
        match self {
            KmerVec::U64(v) => v.clear(),
            KmerVec::U128(v) => v.clear(),
        }
    }

    pub fn len(&self) -> usize {
        match self {
            KmerVec::U64(v) => v.len(),
            KmerVec::U128(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match self {
            KmerVec::U64(v) => v.is_empty(),
            KmerVec::U128(v) => v.is_empty(),
        }
    }
}

/// Decode u64 k-mer (2-bit canonical k-mer)
pub fn decode_u64(kmer: u64, k: u8) -> Vec<u8> {
    (0..k)
        .map(|i| {
            let base_bits = ((kmer >> (2 * i)) & 0b11) as u8;
            unpack_base(base_bits)
        })
        .rev()
        .collect()
}

/// Decode u128 k-mer (2-bit canonical k-mer)
pub fn decode_u128(kmer: u128, k: u8) -> Vec<u8> {
    (0..k)
        .map(|i| {
            let base_bits = ((kmer >> (2 * i)) & 0b11) as u8;
            unpack_base(base_bits)
        })
        .rev()
        .collect()
}

/// Reusable buffers for k-mer computation
#[derive(Clone)]
pub struct Buffers {
    pub packed_nseq: PackedNSeqVec,
    pub positions: Vec<u32>,
    pub kmers: KmerVec,
}

impl Buffers {
    pub fn new_u64() -> Self {
        Self {
            packed_nseq: PackedNSeqVec {
                seq: Default::default(),
                ambiguous: Default::default(),
            },
            positions: Default::default(),
            kmers: KmerVec::U64(Vec::new()),
        }
    }

    pub fn new_u128() -> Self {
        Self {
            packed_nseq: PackedNSeqVec {
                seq: Default::default(),
                ambiguous: Default::default(),
            },
            positions: Default::default(),
            kmers: KmerVec::U128(Vec::new()),
        }
    }
}

/// Fill k-mers vector and positions vector from sequence
pub fn fill_kmers_with_positions(
    seq: &[u8],
    hasher: &SmerHasher,
    kmer_length: u8,
    smer_length: u8,
    buffers: &mut Buffers,
    positions_out: &mut Vec<usize>,
) {
    let Buffers {
        packed_nseq,
        positions,
        kmers,
    } = buffers;

    packed_nseq.seq.clear();
    packed_nseq.ambiguous.clear();
    kmers.clear();
    positions.clear();
    positions_out.clear();

    if seq.len() < kmer_length as usize {
        return;
    }

    packed_nseq.seq.push_ascii(seq);
    packed_nseq.ambiguous.push_ascii(seq);

    // s = 0 bypasses syncmer selection, taking every canonical k-mer
    if smer_length == 0 {
        let k = kmer_length as usize;
        let PackedNSeq {
            seq: packed,
            ambiguous,
        } = packed_nseq.as_slice();
        let n = packed.len() + 1 - k;
        match kmers {
            KmerVec::U64(vec) => {
                for pos in 0..n {
                    if ambiguous.read_kmer(k, pos) == 0 {
                        vec.push(
                            packed
                                .read_kmer(k, pos)
                                .min(packed.read_revcomp_kmer(k, pos)),
                        );
                        positions_out.push(pos);
                    }
                }
            }
            KmerVec::U128(vec) => {
                for pos in 0..n {
                    if ambiguous.read_kmer(k, pos) == 0 {
                        vec.push(
                            packed
                                .read_kmer_u128(k, pos)
                                .min(packed.read_revcomp_kmer_u128(k, pos)),
                        );
                        positions_out.push(pos);
                    }
                }
            }
        }
        return;
    }

    let s = smer_length as usize;
    let w = kmer_length as usize - s + 1;
    let m = simd_minimizers::canonical_open_syncmers(s, w)
        .hasher(hasher)
        .run_skip_ambiguous_windows(packed_nseq.as_slice(), positions);

    match kmers {
        KmerVec::U64(vec) => {
            for (pos, val) in m.pos_and_values_u64() {
                vec.push(val);
                positions_out.push(pos as usize);
            }
        }
        KmerVec::U128(vec) => {
            for (pos, val) in m.pos_and_values_u128() {
                vec.push(val);
                positions_out.push(pos as usize);
            }
        }
    }
}

/// Fill k-mers vector from sequence (without positions)
#[inline]
pub fn fill_kmers(
    seq: &[u8],
    hasher: &SmerHasher,
    kmer_length: u8,
    smer_length: u8,
    buffers: &mut Buffers,
) {
    let Buffers {
        packed_nseq,
        positions,
        kmers,
    } = buffers;

    packed_nseq.seq.clear();
    packed_nseq.ambiguous.clear();
    kmers.clear();
    positions.clear();

    if seq.len() < kmer_length as usize {
        return;
    }

    packed_nseq.seq.push_ascii(seq);
    packed_nseq.ambiguous.push_ascii(seq);

    // s = 0 bypasses syncmer selection, taking every canonical k-mer
    if smer_length == 0 {
        let k = kmer_length as usize;
        let PackedNSeq {
            seq: packed,
            ambiguous,
        } = packed_nseq.as_slice();
        let n = packed.len() + 1 - k;
        match kmers {
            KmerVec::U64(vec) => {
                for pos in 0..n {
                    if ambiguous.read_kmer(k, pos) == 0 {
                        vec.push(
                            packed
                                .read_kmer(k, pos)
                                .min(packed.read_revcomp_kmer(k, pos)),
                        );
                    }
                }
            }
            KmerVec::U128(vec) => {
                for pos in 0..n {
                    if ambiguous.read_kmer(k, pos) == 0 {
                        vec.push(
                            packed
                                .read_kmer_u128(k, pos)
                                .min(packed.read_revcomp_kmer_u128(k, pos)),
                        );
                    }
                }
            }
        }
        return;
    }

    let s = smer_length as usize;
    let w = kmer_length as usize - s + 1;
    let m = simd_minimizers::canonical_open_syncmers(s, w)
        .hasher(hasher)
        .run_skip_ambiguous_windows(packed_nseq.as_slice(), positions);

    match kmers {
        KmerVec::U64(vec) => {
            for (_pos, val) in m.pos_and_values_u64() {
                vec.push(val);
            }
        }
        KmerVec::U128(vec) => {
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
    fn test_fill_kmers() {
        let seq = b"ACGTACGTACGT";
        let k = 5;
        let s = 3;
        let hasher = SmerHasher::new(s as usize);
        let mut buffers = Buffers::new_u64();

        fill_kmers(seq, &hasher, k, s, &mut buffers);

        // We should have at least one k-mer
        assert!(!buffers.kmers.is_empty());

        // Test with a sequence shorter than k
        let short_seq = b"ACGT";
        fill_kmers(short_seq, &hasher, k, s, &mut buffers);
        assert!(buffers.kmers.is_empty());
    }

    #[test]
    fn test_fill_kmers_with_positions() {
        let seq = b"ACGTACGTACGTACGT";
        let k = 7;
        let s = 3;
        let hasher = SmerHasher::new(s as usize);
        let mut buffers = Buffers::new_u64();
        let mut positions = Vec::new();

        fill_kmers_with_positions(seq, &hasher, k, s, &mut buffers, &mut positions);

        // Should have same number of k-mers and positions
        assert_eq!(buffers.kmers.len(), positions.len());

        // All positions should be valid
        for &pos in &positions {
            assert!(pos + k as usize <= seq.len());
        }
    }

    #[test]
    fn test_kmers_match_between_apis() {
        let seq = b"ACGTTGCATGTCGCATGATGCATGAGAGCTACGTTGCATGTCGCATGATGCATGAGAGCT";
        let k = 15;
        let s = 7;
        let hasher = SmerHasher::new(s as usize);

        let mut values_only_buffers = Buffers::new_u64();
        fill_kmers(seq, &hasher, k, s, &mut values_only_buffers);
        let values_only = match &values_only_buffers.kmers {
            KmerVec::U64(v) => v.clone(),
            KmerVec::U128(_) => panic!("Expected u64 k-mers for k <= 32"),
        };

        let mut with_pos_buffers = Buffers::new_u64();
        let mut positions = Vec::new();
        fill_kmers_with_positions(seq, &hasher, k, s, &mut with_pos_buffers, &mut positions);
        let with_pos_values = match &with_pos_buffers.kmers {
            KmerVec::U64(v) => v.clone(),
            KmerVec::U128(_) => panic!("Expected u64 k-mers for k <= 32"),
        };

        // Give us same k-mers and same order from both APIs
        assert_eq!(values_only, with_pos_values);
        assert_eq!(with_pos_values.len(), positions.len());
    }

    /// Pack ASCII into packed-seq's 2-bit encoding (A=0 C=1 T=2 G=3, base i at bit 2i)
    fn pack_packed_seq(seq: &[u8]) -> u128 {
        seq.iter().enumerate().fold(0u128, |acc, (i, &b)| {
            acc | (((b >> 1) & 3) as u128) << (2 * i)
        })
    }

    /// Every canonical k-mer, computed naively from the ASCII sequence
    fn naive_all_kmers(seq: &[u8], k: usize) -> Vec<u128> {
        (0..=seq.len() - k)
            .map(|i| {
                let w = &seq[i..i + k];
                pack_packed_seq(w).min(pack_packed_seq(&revcomp(w)))
            })
            .collect()
    }

    #[test]
    fn test_all_kmers_matches_naive_u64() {
        let seq = pseudo_dna(500, 7);
        let k = 31usize;
        let mut buffers = Buffers::new_u64();
        let mut positions = Vec::new();
        fill_kmers_with_positions(
            &seq,
            &make_hasher(0),
            k as u8,
            0,
            &mut buffers,
            &mut positions,
        );

        let got = match &buffers.kmers {
            KmerVec::U64(v) => v.iter().map(|&x| x as u128).collect::<Vec<_>>(),
            KmerVec::U128(_) => panic!("expected u64"),
        };
        assert_eq!(got, naive_all_kmers(&seq, k));
        assert_eq!(positions, (0..=seq.len() - k).collect::<Vec<_>>());
    }

    #[test]
    fn test_all_kmers_matches_naive_u128() {
        let seq = pseudo_dna(300, 13);
        let k = 41usize;
        let mut buffers = Buffers::new_u128();
        fill_kmers(&seq, &make_hasher(0), k as u8, 0, &mut buffers);

        let got = match &buffers.kmers {
            KmerVec::U128(v) => v.clone(),
            KmerVec::U64(_) => panic!("expected u128"),
        };
        assert_eq!(got, naive_all_kmers(&seq, k));
    }

    #[test]
    fn test_all_kmers_skips_ambiguous() {
        let k = 31usize;
        let mut seq = pseudo_dna(200, 11);
        let clean = kmers_u64(&seq, k as u8, 0, FracMinHash::NONE);
        assert_eq!(clean.len(), seq.len() - k + 1);

        seq[100] = b'N';
        let with_n = kmers_u64(&seq, k as u8, 0, FracMinHash::NONE);

        // Exactly the k k-mers spanning the N are dropped, the rest are untouched
        let mut expected = clean.clone();
        expected.drain(100 + 1 - k..=100);
        assert_eq!(with_n.len(), clean.len() - k);
        assert_eq!(with_n, expected);
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

    fn kmers_u64(seq: &[u8], k: u8, s: u8, fmh: FracMinHash) -> Vec<u64> {
        let hasher = make_hasher(s);
        let mut buffers = Buffers::new_u64();
        fill_kmers(seq, &hasher, k, s, &mut buffers);
        fmh.retain(&mut buffers.kmers, None);
        match &buffers.kmers {
            KmerVec::U64(v) => v.clone(),
            KmerVec::U128(_) => panic!("expected u64"),
        }
    }

    #[test]
    fn test_fmh_none_subset_and_fraction() {
        let seq = pseudo_dna(200_000, 42);
        let all = kmers_u64(&seq, 31, 9, FracMinHash::NONE);
        // from_fraction(1.0) is a no-op equal to NONE
        assert_eq!(all, kmers_u64(&seq, 31, 9, FracMinHash::from_fraction(1.0)));
        let kept = kmers_u64(&seq, 31, 9, FracMinHash::from_fraction(0.1));
        // Deterministic, a strict subset, and roughly the target fraction
        assert_eq!(
            kept,
            kmers_u64(&seq, 31, 9, FracMinHash::from_fraction(0.1))
        );
        let all_set: std::collections::HashSet<u64> = all.iter().copied().collect();
        assert!(kept.iter().all(|v| all_set.contains(v)));
        let realised = kept.len() as f64 / all.len() as f64;
        assert!(
            (realised - 0.1).abs() < 0.03,
            "realised fraction {realised}"
        );
    }

    /// Pack ASCII bases as k-mer values are stored
    fn pack(seq: &[u8]) -> u128 {
        seq.iter().rev().fold(0u128, |acc, &b| {
            let bits = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!("non-ACGT base"),
            };
            (acc << 2) | bits
        })
    }

    fn revcomp(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => panic!("non-ACGT base"),
            })
            .collect()
    }

    #[test]
    fn test_kdust_k3_is_not_nan() {
        // k=3 has a single triplet and would divide by zero without the guard.
        // validate_k_s forces odd k > s >= 1, so k<3 is unreachable
        let score = calculate_kdust(pack(b"AAA"), 3);
        assert!(score.is_finite(), "k=3 scored {score}");
        assert_eq!(score, 1.0);
    }

    #[test]
    fn test_kdust_is_reverse_complement_invariant() {
        // revcomp only relabels triplet types, so the summed score is unchanged
        for seed in [1u64, 7, 99] {
            let seq = pseudo_dna(31, seed);
            let rc = revcomp(&seq);
            assert_eq!(
                calculate_kdust(pack(&seq), 31),
                calculate_kdust(pack(&rc), 31)
            );
        }
        let skewed = b"AAAAACCCCCAAAAAGGGGGAAAAATTTTTA";
        assert_eq!(
            calculate_kdust(pack(skewed), 31),
            calculate_kdust(pack(&revcomp(skewed)), 31)
        );
    }

    #[test]
    fn test_kdust_scores_and_ranks_by_repetitiveness() {
        // Bounds: every triplet identical scores 0, no repeated triplet scores 1
        let distinct = b"ACGTAACCGGTTACAGATCCTGGCATTGACT";
        assert_eq!(distinct.len(), 31);
        assert_eq!(calculate_kdust(pack(distinct), 31), 1.0);

        // Only a homopolymer bottoms out at 0, shorter-period repeats sit in between
        let homopolymer = calculate_kdust(pack(&[b'A'; 31]), 31);
        let dinucleotide = calculate_kdust(pack(b"ATATATATATATATATATATATATATATATA"), 31);
        let trinucleotide = calculate_kdust(pack(b"ATTATTATTATTATTATTATTATTATTATTA"), 31);
        let random = calculate_kdust(pack(&pseudo_dna(31, 5)), 31);

        assert_eq!(homopolymer, 0.0);
        assert!(random > 0.9, "random 31-mer scored {random}");
        assert!(
            homopolymer < dinucleotide && dinucleotide < trinucleotide && trinucleotide < random,
            "expected monotonic ordering, got {homopolymer} {dinucleotide} {trinucleotide} {random}"
        );
    }

    #[test]
    fn test_kdust_none_is_a_noop() {
        let seq = pseudo_dna(2_000, 11);
        let all = kmers_u64(&seq, 31, 9, FracMinHash::NONE);
        assert!(!all.is_empty());

        let kdust = Kdust::from_threshold(0.0, 31);
        assert!(kdust.is_none());
        let mut vec = KmerVec::U64(all.clone());
        kdust.retain(&mut vec, None);
        match &vec {
            KmerVec::U64(v) => assert_eq!(*v, all),
            KmerVec::U128(_) => panic!("expected u64"),
        }
    }

    #[test]
    fn test_kdust_retain_drops_below_threshold_and_keeps_positions_aligned() {
        // Hand-built values pin retain's behaviour, not k-mer selection
        let low = pack(&[b'A'; 31]) as u64;
        let mid = pack(b"ATATATATATATATATATATATATATATATA") as u64;
        let high = pack(&pseudo_dna(31, 5)) as u64;
        assert!(calculate_kdust(high as u128, 31) > 0.9);

        let mut kmers = KmerVec::U64(vec![high, low, high, mid, low]);
        let mut positions = vec![10usize, 20, 30, 40, 50];
        Kdust::from_threshold(0.9, 31).retain(&mut kmers, Some(&mut positions));

        match &kmers {
            KmerVec::U64(v) => assert_eq!(*v, vec![high, high]),
            KmerVec::U128(_) => panic!("expected u64"),
        }
        // Survivors keep their paired positions
        assert_eq!(positions, vec![10, 30]);
    }

    #[test]
    fn test_fmh_uses_fixed_rapid_hash_of_value() {
        let fmh = FracMinHash::from_fraction(0.5);
        let kmer64 = 0x0123_4567_89ab_cdef;
        let kmer128 = 0x0123_4567_89ab_cdef_fedc_ba98_7654_3210;

        assert_eq!(
            fmh.keeps_u64(kmer64),
            FixedRapidHasher.hash_one(kmer64) <= fmh.threshold
        );
        assert_eq!(
            fmh.keeps_u128(kmer128),
            FixedRapidHasher.hash_one(kmer128) <= fmh.threshold
        );
    }
}

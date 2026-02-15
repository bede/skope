use packed_seq::{PackedNSeqVec, SeqVec, unpack_base};

pub const DEFAULT_KMER_LENGTH: u8 = 31;
pub const DEFAULT_SMER_LENGTH: u8 = 15;

pub type KmerHasher = simd_minimizers::seq_hash::NtHasher<true, 1>;

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

    let k = kmer_length as u32;
    match syncmers {
        SyncmerVec::U64(vec) => {
            let mut next_allowed: u32 = 0;
            for (pos, val) in m.pos_and_values_u64() {
                if pos >= next_allowed {
                    vec.push(val);
                    positions_out.push(pos as usize);
                    next_allowed = pos.saturating_add(k);
                }
            }
        }
        SyncmerVec::U128(vec) => {
            let mut next_allowed: u32 = 0;
            for (pos, val) in m.pos_and_values_u128() {
                if pos >= next_allowed {
                    vec.push(val);
                    positions_out.push(pos as usize);
                    next_allowed = pos.saturating_add(k);
                }
            }
        }
    }
}

/// Fill syncmers vector from sequence (without positions)
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

    let k = kmer_length as u32;
    match syncmers {
        SyncmerVec::U64(vec) => {
            let mut next_allowed: u32 = 0;
            for (pos, val) in m.pos_and_values_u64() {
                if pos >= next_allowed {
                    vec.push(val);
                    next_allowed = pos.saturating_add(k);
                }
            }
        }
        SyncmerVec::U128(vec) => {
            let mut next_allowed: u32 = 0;
            for (pos, val) in m.pos_and_values_u128() {
                if pos >= next_allowed {
                    vec.push(val);
                    next_allowed = pos.saturating_add(k);
                }
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

        // Syncmers should be non-overlapping by construction
        for window in positions.windows(2) {
            assert!(window[1] >= window[0] + k as usize);
        }
    }

    #[test]
    fn test_non_overlapping_syncmers_match_between_apis() {
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

        // Both APIs should retain the same non-overlapping syncmers in the same order.
        assert_eq!(values_only, with_pos_values);
        assert_eq!(with_pos_values.len(), positions.len());

        for window in positions.windows(2) {
            assert!(window[1] >= window[0] + k as usize);
        }
    }
}

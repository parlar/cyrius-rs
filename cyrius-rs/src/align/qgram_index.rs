/// Q-gram hash index over reference sequences.
///
/// For each q-gram (substring of length q) in the reference, stores a list
/// of (seq_id, position) hits. Uses a rolling 2-bit hash for O(1) lookup.
///
/// For q <= 12, uses a direct-addressed array of size 4^q.
/// For q > 12, falls back to HashMap.
///
/// Adapted from sv_tools/sv_align.

use std::collections::HashMap;

use crate::align::sequence::encode_base;

/// A hit in the q-gram index.
#[derive(Clone, Copy, Debug)]
pub struct QgramHit {
    /// Sequence index (0 = forward, 1 = reverse complement).
    pub seq_id: u16,
    /// 0-based start position of the q-gram within the sequence.
    pub pos: u32,
}

/// Internal storage for the q-gram table.
enum QgramTable {
    /// Direct-addressed array indexed by q-gram hash (used for q <= 12).
    Direct(Vec<Vec<QgramHit>>),
    /// HashMap fallback for large q (q > 12).
    Hash(HashMap<u64, Vec<QgramHit>>),
}

/// Q-gram index: maps q-gram hash -> list of reference hits.
pub struct QgramIndex {
    pub q: usize,
    table: QgramTable,
}

/// Maximum q for direct-addressed table (4^12 = 16M entries).
const MAX_DIRECT_Q: usize = 12;

/// Maximum hits per q-gram before it's considered repetitive and skipped.
/// Prevents excessive candidate enumeration in repeat-rich regions.
const MAX_HITS_PER_QGRAM: usize = 500;

impl QgramIndex {
    /// Build a q-gram index from a set of sequences.
    ///
    /// Each entry is (sequence_bytes, is_reverse_complement).
    /// Sequences are assigned sequential IDs starting from 0.
    /// Q-grams containing N bases are skipped.
    pub fn build(sequences: &[(&[u8], bool)], q: usize) -> Self {
        assert!(q > 0 && q <= 32, "q must be in 1..=32, got {}", q);

        let mask: u64 = if q >= 32 {
            u64::MAX
        } else {
            (1u64 << (2 * q)) - 1
        };
        let use_direct = q <= MAX_DIRECT_Q;

        let mut direct: Vec<Vec<QgramHit>> = if use_direct {
            vec![Vec::new(); 1usize << (2 * q)]
        } else {
            Vec::new()
        };
        let mut hash_map: HashMap<u64, Vec<QgramHit>> = HashMap::new();

        for (seq_idx, (seq, _is_rc)) in sequences.iter().enumerate() {
            let seq_id = seq_idx as u16;
            let mut hash: u64 = 0;
            let mut valid: usize = 0;

            for (i, &base) in seq.iter().enumerate() {
                let c = encode_base(base);
                if c >= 4 {
                    valid = 0;
                    hash = 0;
                    continue;
                }
                hash = ((hash << 2) | c as u64) & mask;
                valid += 1;

                if valid >= q {
                    let pos = (i + 1 - q) as u32;
                    let hit = QgramHit { seq_id, pos };
                    if use_direct {
                        direct[hash as usize].push(hit);
                    } else {
                        hash_map.entry(hash).or_default().push(hit);
                    }
                }
            }
        }

        let table = if use_direct {
            QgramTable::Direct(direct)
        } else {
            QgramTable::Hash(hash_map)
        };

        QgramIndex { q, table }
    }

    /// Look up a q-gram hash, returning all reference hits.
    /// Returns empty slice if the q-gram is too repetitive (>MAX_HITS_PER_QGRAM).
    #[inline]
    pub fn lookup(&self, hash: u64) -> &[QgramHit] {
        let hits: &[QgramHit] = match &self.table {
            QgramTable::Direct(vec) => {
                let idx = hash as usize;
                if idx < vec.len() {
                    vec[idx].as_slice()
                } else {
                    &[]
                }
            }
            QgramTable::Hash(map) => map.get(&hash).map_or(&[], |v| v.as_slice()),
        };
        // Skip repetitive q-grams to avoid O(n^2) in repeat regions.
        if hits.len() > MAX_HITS_PER_QGRAM {
            &[]
        } else {
            hits
        }
    }
}

/// Compute the 2-bit hash for a q-gram (ASCII bytes).
/// Returns `None` if any base is N or unknown.
pub fn qgram_hash(seq: &[u8]) -> Option<u64> {
    let q = seq.len();
    assert!(q <= 32);
    let mut hash: u64 = 0;
    for &base in seq {
        let c = encode_base(base);
        if c >= 4 {
            return None;
        }
        hash = (hash << 2) | c as u64;
    }
    Some(hash)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_qgram_hash() {
        assert_eq!(qgram_hash(b"A"), Some(0b00));
        assert_eq!(qgram_hash(b"C"), Some(0b01));
        assert_eq!(qgram_hash(b"G"), Some(0b10));
        assert_eq!(qgram_hash(b"T"), Some(0b11));
        assert_eq!(qgram_hash(b"ACGT"), Some(0b00011011));
        assert_eq!(qgram_hash(b"N"), None);
    }

    #[test]
    fn test_build_simple() {
        let seq = b"ACGTACGT";
        let idx = QgramIndex::build(&[(seq.as_slice(), false)], 4);
        let hash_acgt = qgram_hash(b"ACGT").unwrap();
        let hits = idx.lookup(hash_acgt);
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].pos, 0);
        assert_eq!(hits[1].pos, 4);
    }

    #[test]
    fn test_two_sequences() {
        let fwd = b"ACGTACGT";
        let rc = b"ACGTACGT";
        let idx =
            QgramIndex::build(&[(fwd.as_slice(), false), (rc.as_slice(), true)], 4);
        let hash = qgram_hash(b"ACGT").unwrap();
        let hits = idx.lookup(hash);
        // 2 hits in fwd + 2 hits in rc = 4
        assert_eq!(hits.len(), 4);
    }

    #[test]
    fn test_n_base_reset() {
        let seq = b"ACNGT";
        let idx = QgramIndex::build(&[(seq.as_slice(), false)], 2);
        let hash_ac = qgram_hash(b"AC").unwrap();
        assert_eq!(idx.lookup(hash_ac).len(), 1);
        let hash_gt = qgram_hash(b"GT").unwrap();
        assert_eq!(idx.lookup(hash_gt).len(), 1);
    }
}

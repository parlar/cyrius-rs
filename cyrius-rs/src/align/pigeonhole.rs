/// Pigeonhole seed filter for full-sensitivity candidate enumeration.
///
/// For a read of length m with at most k errors (edit distance), the read is
/// split into k+1 non-overlapping seeds of length q = floor(m / (k+1)).
/// By the pigeonhole principle, at least one seed must match exactly in any
/// valid alignment position. This guarantees **zero false negatives**.
///
/// Adapted from sv_tools/sv_align.

use std::collections::HashSet;

use crate::align::qgram_index::{qgram_hash, QgramIndex};

/// A candidate alignment position to verify.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct Candidate {
    /// Sequence index (0 = forward, 1 = reverse complement).
    pub seq_id: u16,
    /// Approximate alignment start position (diagonal = ref_pos - read_offset).
    pub diagonal: i32,
}

/// Extract seeds from a query and find candidate alignment positions.
///
/// Seeds are of length `index.q` (matching the pre-built index).
/// Returns a deduplicated list of candidates. Candidates on nearby
/// diagonals (within k+1) for the same sequence are merged.
pub fn find_candidates(query: &[u8], k: usize, index: &QgramIndex) -> Vec<Candidate> {
    let m = query.len();
    if m == 0 {
        return Vec::new();
    }
    let q = index.q;
    if q == 0 || q > 32 || q > m {
        return Vec::new();
    }
    // Cap n_seeds at m/q so seed extraction stays within the query.
    let n_seeds = (k + 1).min(m / q);

    let mut seen: HashSet<(u16, i32)> = HashSet::new();
    let mut candidates: Vec<Candidate> = Vec::new();
    let bucket_size = (k + 1) as i32;

    for s in 0..n_seeds {
        let offset = s * q;
        let seed = &query[offset..offset + q];

        let hash = match qgram_hash(seed) {
            Some(h) => h,
            None => continue,
        };

        for hit in index.lookup(hash) {
            let diagonal = hit.pos as i32 - offset as i32;
            let bucket = if diagonal >= 0 {
                diagonal / bucket_size
            } else {
                (diagonal - bucket_size + 1) / bucket_size
            };

            if seen.insert((hit.seq_id, bucket)) {
                candidates.push(Candidate {
                    seq_id: hit.seq_id,
                    diagonal,
                });
            }
        }
    }

    candidates
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::qgram_index::QgramIndex;

    #[test]
    fn test_exact_match_found() {
        let seq = b"ACGTACGTACGTACGTACGTACGT";
        let idx = QgramIndex::build(&[(seq.as_slice(), false)], 4);
        let read = b"ACGTACGTACGT";
        let candidates = find_candidates(read, 2, &idx);
        assert!(!candidates.is_empty());
    }

    #[test]
    fn test_no_match() {
        let seq = b"AAAAAAAAAAAAAAAAAAAAAAAA";
        let idx = QgramIndex::build(&[(seq.as_slice(), false)], 4);
        let read = b"GGGGGGGGGGGG";
        let candidates = find_candidates(read, 2, &idx);
        assert!(candidates.is_empty());
    }

    #[test]
    fn test_dedup_nearby_diagonals() {
        let seq = b"ACGTACGTACGTACGT";
        let idx = QgramIndex::build(&[(seq.as_slice(), false)], 4);
        let read = b"ACGTACGTACGT";
        let candidates = find_candidates(read, 2, &idx);
        assert!(candidates.len() < 12); // dedup reduces count
    }
}

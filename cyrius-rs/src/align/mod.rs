//! Local sequence alignment for SV breakpoint analysis.
//!
//! Maps short sequences (soft clips, mismatch segments) back to a reference
//! region to find their origin. Uses q-gram indexed seeding, pigeonhole
//! filtering for zero false negatives, and Myers bit-parallel DP for
//! verification.
//!
//! Adapted from sv_tools/sv_align. Zero external dependencies.

pub mod myers;
pub mod pigeonhole;
pub mod qgram_index;
pub mod sequence;

use myers::myers_semiglobal;
use pigeonhole::find_candidates;
use qgram_index::QgramIndex;
use sequence::{reverse_complement, PeqTable};

/// A mapping hit: position where a query maps in the reference.
#[derive(Debug, Clone)]
pub struct MappingHit {
    /// 0-based position in the forward reference.
    pub ref_pos: usize,
    /// Edit distance (substitutions + insertions + deletions).
    pub edit_distance: u32,
    /// True if the hit is on the reverse complement strand.
    pub is_reverse: bool,
}

/// Local alignment index for a reference region.
///
/// Indexes both forward and reverse complement strands. Provides
/// guaranteed-complete search: all positions where a query maps
/// within the specified error tolerance will be found.
pub struct LocalIndex {
    index: QgramIndex,
    fwd_seq: Vec<u8>,
    rc_seq: Vec<u8>,
    q: usize,
}

/// Default q-gram length. q=8 balances sensitivity for short clips (20bp)
/// with specificity in genomic reference (65K possible 8-mers).
pub const DEFAULT_Q: usize = 8;

impl LocalIndex {
    /// Build an index for a reference sequence.
    ///
    /// Both forward and reverse complement strands are indexed.
    /// `q` is the q-gram length (8 recommended for SV analysis).
    pub fn build(ref_seq: &[u8], q: usize) -> Self {
        let rc_seq = reverse_complement(ref_seq);
        let index = QgramIndex::build(&[(ref_seq, false), (&rc_seq, true)], q);
        LocalIndex {
            index,
            fwd_seq: ref_seq.to_vec(),
            rc_seq,
            q,
        }
    }

    /// Map a query sequence against the indexed reference.
    ///
    /// Returns all positions where the query aligns with edit distance
    /// <= `max_edit_distance`. Uses pigeonhole seeding (zero false negatives)
    /// and Myers bit-parallel DP for verification.
    ///
    /// Both strands are searched. Hits are deduplicated (nearby positions
    /// merged) and sorted by edit distance (best first).
    pub fn map_sequence(&self, query: &[u8], max_edit_distance: u32) -> Vec<MappingHit> {
        if query.is_empty() || query.len() < self.q {
            return Vec::new();
        }

        let k = max_edit_distance as usize;
        let m = query.len();

        let candidates = find_candidates(query, k, &self.index);
        if candidates.is_empty() {
            return Vec::new();
        }

        let peq = PeqTable::build(query);
        let seqs: [&[u8]; 2] = [&self.fwd_seq, &self.rc_seq];
        let mut hits = Vec::new();

        for cand in &candidates {
            let seq_idx = cand.seq_id as usize;
            if seq_idx >= 2 {
                continue;
            }
            let seq = seqs[seq_idx];
            let seq_len = seq.len();
            let is_rc = seq_idx == 1;

            // Verification window: [diagonal - k, diagonal + m + k].
            let window_start = (cand.diagonal - k as i32).max(0) as usize;
            let window_end_i64 = cand.diagonal as i64 + m as i64 + k as i64;
            let window_end = if window_end_i64 <= 0 {
                continue;
            } else {
                (window_end_i64 as usize).min(seq_len)
            };

            if window_start >= window_end {
                continue;
            }

            let window = &seq[window_start..window_end];

            if let Some(result) = myers_semiglobal(&peq, window, max_edit_distance) {
                let aln_end = window_start + result.text_end;
                let ref_pos = if is_rc {
                    // Convert RC position to forward-strand coordinate.
                    seq_len.saturating_sub(1 + aln_end)
                } else {
                    aln_end.saturating_sub(m - 1)
                };

                hits.push(MappingHit {
                    ref_pos,
                    edit_distance: result.edit_distance,
                    is_reverse: is_rc,
                });
            }
        }

        // Deduplicate: merge hits at nearby positions, keep best.
        hits.sort_by_key(|h| (h.ref_pos, h.edit_distance));
        deduplicate_hits(&mut hits);

        // Sort by edit distance (best first).
        hits.sort_by_key(|h| h.edit_distance);

        hits
    }

    /// Reference length (forward strand).
    pub fn ref_len(&self) -> usize {
        self.fwd_seq.len()
    }
}

/// Remove near-duplicate hits (within 10bp), keeping the best edit distance.
fn deduplicate_hits(hits: &mut Vec<MappingHit>) {
    if hits.len() <= 1 {
        return;
    }
    let mut kept = Vec::with_capacity(hits.len());
    for hit in hits.iter() {
        let dominated = kept.iter().any(|k: &MappingHit| {
            k.ref_pos.abs_diff(hit.ref_pos) <= 10
                && k.is_reverse == hit.is_reverse
                && k.edit_distance <= hit.edit_distance
        });
        if !dominated {
            kept.push(hit.clone());
        }
    }
    *hits = kept;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exact_match() {
        let reference = b"AAAAAAAAAAACGTACGTACGTACGTACGTAAAAAAAAAA";
        let index = LocalIndex::build(reference, 8);

        let query = b"ACGTACGTACGTACGTACGT";
        let hits = index.map_sequence(query, 2);
        assert!(!hits.is_empty(), "should find exact match");
        assert_eq!(hits[0].edit_distance, 0);
    }

    #[test]
    fn test_match_with_errors() {
        let mut reference = vec![b'T'; 10];
        reference.extend_from_slice(b"CGATCTAGCATGCTTACGCA");
        reference.extend_from_slice(&vec![b'T'; 10]);
        let index = LocalIndex::build(&reference, 5);

        let query = b"CGATCTAGCAAGCTTACGCA";
        let hits = index.map_sequence(query, 3);
        assert!(!hits.is_empty(), "should find match with 1 error");
        assert!(hits[0].edit_distance <= 2);
    }

    #[test]
    fn test_no_match() {
        let reference = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let index = LocalIndex::build(reference, 8);

        let query = b"GGGGGGGGGGGGGGGGGGGG";
        let hits = index.map_sequence(query, 2);
        assert!(hits.is_empty());
    }

    #[test]
    fn test_reverse_complement_match() {
        let mut reference = vec![b'T'; 10];
        reference.extend_from_slice(b"CGATCTAGCATGCTTACGCA");
        reference.extend_from_slice(&vec![b'T'; 10]);
        let index = LocalIndex::build(&reference, 5);

        let segment = b"CGATCTAGCATGCTTACGCA";
        let query = sequence::reverse_complement(segment);
        assert_ne!(&query, segment.as_slice());

        let hits = index.map_sequence(&query, 2);
        assert!(!hits.is_empty(), "should find RC match");
        assert_eq!(hits[0].edit_distance, 0);
    }

    #[test]
    fn test_multiple_hits() {
        let mut reference = vec![b'A'; 100];
        let segment = b"CGTACGTACGTACGTACGT";
        reference[10..10 + segment.len()].copy_from_slice(segment);
        reference[60..60 + segment.len()].copy_from_slice(segment);

        let index = LocalIndex::build(&reference, 8);
        let hits = index.map_sequence(segment, 0);
        assert!(hits.len() >= 2, "should find 2 hits, got {}", hits.len());
    }

    #[test]
    fn test_short_query_skipped() {
        let reference = b"ACGTACGTACGTACGTACGT";
        let index = LocalIndex::build(reference, 8);

        let hits = index.map_sequence(b"ACGT", 2);
        assert!(hits.is_empty(), "query shorter than q should return empty");
    }

    #[test]
    fn test_deletion_scenario() {
        let flank_left = b"ACGTACGTACGTACGTACGT";
        let deleted = b"GGGGGGGGGGGGGGGGGGGG";
        let flank_right = b"TGATCGATCGATCGATCGAT";

        let mut reference = Vec::new();
        reference.extend_from_slice(flank_left);
        reference.extend_from_slice(deleted);
        reference.extend_from_slice(flank_right);

        let index = LocalIndex::build(&reference, 8);

        let clip = flank_right;
        let hits = index.map_sequence(clip, 2);
        assert!(!hits.is_empty());
        let pos = hits[0].ref_pos;
        assert!(
            pos >= 37 && pos <= 41,
            "expected position near 39, got {}",
            pos
        );
        assert_eq!(hits[0].edit_distance, 0);
    }
}

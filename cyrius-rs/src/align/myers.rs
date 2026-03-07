/// Multi-word Myers bit-parallel semi-global edit distance.
///
/// Implements the Myers 1999 algorithm for approximate string matching.
/// Each DP column is encoded as bit-vectors (Pv, Mv) and advanced in
/// O(ceil(m/64)) operations per text character.
///
/// Semi-global mode: the pattern (read) must be fully consumed, but can
/// start and end at any position in the text (reference window).
///
/// Supports patterns up to 320 bp via 5 x u64 words.
///
/// Adapted from sv_tools/sv_align.

use crate::align::sequence::{encode_base, PeqTable, MAX_WORDS};

/// Result of a successful semi-global alignment.
#[derive(Debug, Clone)]
pub struct MyersResult {
    /// Edit distance (substitutions + insertions + deletions).
    pub edit_distance: u32,
    /// 0-based end position in the text (inclusive).
    pub text_end: usize,
}

/// Compute the semi-global edit distance between a pattern and text.
///
/// Returns the best alignment (minimum edit distance <= max_k) if one
/// exists, along with the text end position.
///
/// Semi-global boundary: D[0][j] = 0 for all j (free start in text).
pub fn myers_semiglobal(peq: &PeqTable, text: &[u8], max_k: u32) -> Option<MyersResult> {
    let m = peq.pattern_len;
    let n = text.len();

    if m == 0 {
        return Some(MyersResult {
            edit_distance: 0,
            text_end: 0,
        });
    }
    if n == 0 {
        return None;
    }

    let n_words = peq.n_words;

    // Initialize: Pv = all 1s, Mv = all 0s. Score starts at m.
    let mut pv = [u64::MAX; MAX_WORDS];
    let mut mv = [0u64; MAX_WORDS];
    let mut score = m as i32;

    let last_word = n_words - 1;
    let last_bit = (m - 1) % 64;
    let last_mask = 1u64 << last_bit;

    let mut best_score = (max_k as i32) + 1;
    let mut best_j: usize = 0;

    for j in 0..n {
        let c = encode_base(text[j]);

        // Pass 1: Compute Xv and Xh with addition carry chain.
        let mut xh = [0u64; MAX_WORDS];
        let mut xv = [0u64; MAX_WORDS];
        let mut carry_add: u64 = 0;

        for w in 0..n_words {
            let eq = peq.words[c][w];
            xv[w] = eq | mv[w];

            let eq_and_pv = eq & pv[w];
            let (s1, c1) = eq_and_pv.overflowing_add(pv[w]);
            let (s2, c2) = s1.overflowing_add(carry_add);
            carry_add = (c1 as u64) | (c2 as u64);

            xh[w] = (s2 ^ pv[w]) | eq;
        }

        // Pass 2: Compute Ph and Mh.
        let mut ph = [0u64; MAX_WORDS];
        let mut mh = [0u64; MAX_WORDS];

        for w in 0..n_words {
            ph[w] = mv[w] | !(xh[w] | pv[w]);
            mh[w] = pv[w] & xh[w];
        }

        // Score delta from the last active bit.
        if ph[last_word] & last_mask != 0 {
            score += 1;
        } else if mh[last_word] & last_mask != 0 {
            score -= 1;
        }

        // Pass 3: Shift Ph and Mh left, update Pv/Mv.
        let mut carry_ph: u64 = 0;
        let mut carry_mh: u64 = 0;

        for w in 0..n_words {
            let ph_shifted = (ph[w] << 1) | carry_ph;
            carry_ph = ph[w] >> 63;

            let mh_shifted = (mh[w] << 1) | carry_mh;
            carry_mh = mh[w] >> 63;

            pv[w] = mh_shifted | !(xv[w] | ph_shifted);
            mv[w] = ph_shifted & xv[w];
        }

        if score <= max_k as i32 && score < best_score {
            best_score = score;
            best_j = j;
        }
    }

    if best_score <= max_k as i32 {
        Some(MyersResult {
            edit_distance: best_score as u32,
            text_end: best_j,
        })
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::sequence::PeqTable;

    fn check_ed(pattern: &[u8], text: &[u8], expected: u32) {
        let peq = PeqTable::build(pattern);
        let result = myers_semiglobal(&peq, text, expected + 5);
        assert!(result.is_some(), "expected ed={}, got None", expected);
        assert_eq!(result.unwrap().edit_distance, expected);
    }

    #[test]
    fn test_exact_match() {
        check_ed(b"ACGT", b"ACGT", 0);
    }

    #[test]
    fn test_exact_match_embedded() {
        check_ed(b"ACGT", b"XXXACGTXXX", 0);
    }

    #[test]
    fn test_single_substitution() {
        check_ed(b"ACGT", b"ACAT", 1);
    }

    #[test]
    fn test_single_insertion() {
        check_ed(b"ACGT", b"ACGGT", 1);
    }

    #[test]
    fn test_single_deletion() {
        check_ed(b"ACGT", b"ACT", 1);
    }

    #[test]
    fn test_threshold_exceeded() {
        let peq = PeqTable::build(b"ACGT");
        let result = myers_semiglobal(&peq, b"TTTT", 1);
        assert!(result.is_none());
    }

    #[test]
    fn test_multi_word_130bp() {
        let pattern: Vec<u8> = (0..130)
            .map(|i| match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();
        let mut text = vec![b'G'; 5];
        text.extend_from_slice(&pattern);
        text.extend_from_slice(&[b'G'; 5]);
        let peq = PeqTable::build(&pattern);
        let result = myers_semiglobal(&peq, &text, 10).unwrap();
        assert_eq!(result.edit_distance, 0);
    }

    #[test]
    fn test_semiglobal_best_position() {
        let peq = PeqTable::build(b"GATT");
        let result = myers_semiglobal(&peq, b"AAAGATTCCC", 3).unwrap();
        assert_eq!(result.edit_distance, 0);
        assert_eq!(result.text_end, 6);
    }
}

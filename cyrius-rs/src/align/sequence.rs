/// DNA encoding, reverse complement, and Peq bitmask tables for Myers DP.
///
/// Adapted from sv_tools/sv_align.

/// Number of distinct encoded bases (A, C, G, T, N).
pub const ALPHABET_SIZE: usize = 5;

/// Maximum number of u64 words for Peq tables (supports patterns up to 320 bp).
pub const MAX_WORDS: usize = 5;

/// Encode an ASCII DNA base to a 0-4 integer.
/// A=0, C=1, G=2, T=3, anything else (including N)=4.
#[inline]
pub fn encode_base(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// Complement a single ASCII DNA base.
#[inline]
pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'a' => b't',
        b'T' => b'A',
        b't' => b'a',
        b'C' => b'G',
        b'c' => b'g',
        b'G' => b'C',
        b'g' => b'c',
        _ => b'N',
    }
}

/// Compute the reverse complement of a DNA sequence.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Peq (pattern equality) bitmask table for Myers bit-parallel DP.
///
/// For each base `c` in {A, C, G, T, N}, `words[c][w]` has bit `i` set
/// (where `i = position - w*64`) if `pattern[position] == c`.
///
/// Supports patterns up to `MAX_WORDS * 64 = 320` bp.
pub struct PeqTable {
    pub words: [[u64; MAX_WORDS]; ALPHABET_SIZE],
    pub n_words: usize,
    pub pattern_len: usize,
}

impl PeqTable {
    /// Build a Peq table from a DNA pattern (ASCII bytes).
    pub fn build(pattern: &[u8]) -> Self {
        let m = pattern.len();
        assert!(
            m <= MAX_WORDS * 64,
            "pattern too long: {} > {}",
            m,
            MAX_WORDS * 64
        );

        let n_words = (m + 63) / 64;
        let mut words = [[0u64; MAX_WORDS]; ALPHABET_SIZE];

        for (i, &base) in pattern.iter().enumerate() {
            let c = encode_base(base);
            let word_idx = i / 64;
            let bit_idx = i % 64;
            words[c][word_idx] |= 1u64 << bit_idx;
        }

        PeqTable {
            words,
            n_words,
            pattern_len: m,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_base() {
        assert_eq!(encode_base(b'A'), 0);
        assert_eq!(encode_base(b'C'), 1);
        assert_eq!(encode_base(b'G'), 2);
        assert_eq!(encode_base(b'T'), 3);
        assert_eq!(encode_base(b'N'), 4);
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GATTACA"), b"TGTAATC");
    }

    #[test]
    fn test_peq_simple() {
        let peq = PeqTable::build(b"ACGT");
        assert_eq!(peq.n_words, 1);
        assert_eq!(peq.words[0][0], 0b0001); // A at pos 0
        assert_eq!(peq.words[1][0], 0b0010); // C at pos 1
        assert_eq!(peq.words[2][0], 0b0100); // G at pos 2
        assert_eq!(peq.words[3][0], 0b1000); // T at pos 3
    }
}

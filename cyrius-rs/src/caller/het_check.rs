//! Hemizygosity detection via heterozygous site counting.
//!
//! D6/D7 paralog het check: counts positions where both D6 and D7 signal
//! are present. Detects *5/*5 (homozygous deletion) where D6 signal is
//! absent at nearly all positions. Cannot detect single-copy deletions
//! because the remaining D6 copy still provides signal at all positions.
//!
//! Currently used for diagnostic logging only.

/// Result of hemizygosity check.
#[derive(Debug, Clone)]
pub struct HetCheckResult {
    /// Number of D6/D7 diagnostic positions with heterozygous signal
    pub n_het_sites: usize,
    /// Number of D6/D7 positions with adequate depth
    pub n_adequate_depth: usize,
    /// Whether the sample appears hemizygous (D6/D7 level)
    pub is_hemizygous: bool,
    /// Het fraction (n_het_sites / n_adequate_depth)
    pub het_fraction: f64,
}

/// Check for hemizygosity using D6/D7 paralog diagnostic positions.
///
/// Primarily detects *5/*5 (homozygous deletion): zero D6 signal at most
/// positions means both D6 copies are deleted.
pub fn check_hemizygosity(
    snp_d6: &[usize],
    snp_d7: &[usize],
    min_depth: usize,
    min_minor_count: usize,
    min_het_frac: f64,
) -> HetCheckResult {
    let n = snp_d6.len().min(snp_d7.len());
    let mut n_het = 0;
    let mut n_adequate = 0;

    for i in 0..n {
        let total = snp_d6[i] + snp_d7[i];
        if total < min_depth {
            continue;
        }
        n_adequate += 1;

        let minor = snp_d6[i].min(snp_d7[i]);
        if minor >= min_minor_count {
            let minor_frac = minor as f64 / total as f64;
            if minor_frac >= min_het_frac {
                n_het += 1;
            }
        }
    }

    let het_fraction = if n_adequate > 0 {
        n_het as f64 / n_adequate as f64
    } else {
        0.0
    };

    let is_hemizygous = n_adequate >= 20 && n_het <= 2;

    HetCheckResult {
        n_het_sites: n_het,
        n_adequate_depth: n_adequate,
        is_hemizygous,
        het_fraction,
    }
}


//! Hemizygosity detection via heterozygous site counting.
//!
//! Two complementary approaches:
//! 1. D6/D7 paralog het check: counts positions where both D6 and D7 signal
//!    are present. Detects *5/*5 (homozygous deletion) where D6 signal is
//!    absent at nearly all positions. Cannot detect single-copy deletions
//!    because the remaining D6 copy still provides signal at all positions.
//! 2. Variant-level het check: among called variant positions, count how many
//!    are heterozygous (both alt and ref reads). Zero het variants suggests
//!    hemizygosity (one copy deleted, remaining copy is uniform).

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

/// Result of variant-level heterozygosity check.
#[derive(Debug, Clone)]
pub struct VariantHetResult {
    /// Number of variant positions with heterozygous signal (both alt > 0 and ref > 0)
    pub n_het_variants: usize,
    /// Number of variant positions with any alt reads and adequate depth
    pub n_called_variants: usize,
    /// Whether the variant pattern suggests hemizygosity
    pub is_hemizygous: bool,
}

/// Check variant-level heterozygosity.
///
/// Among positions where alt variants are called (alt > 0), check how many
/// also have substantial ref reads. In a sample with two distinct D6 copies,
/// variants unique to one copy will appear heterozygous (some ref reads from
/// the other copy). If all variant positions show only alt reads (no ref),
/// it suggests only one D6 copy exists (hemizygous).
pub fn check_variant_heterozygosity(
    var_alt: &[usize],
    var_ref: &[usize],
    min_depth: usize,
    min_alt_count: usize,
) -> VariantHetResult {
    let n = var_alt.len().min(var_ref.len());
    let mut n_called = 0;
    let mut n_het = 0;

    for i in 0..n {
        let total = var_alt[i] + var_ref[i];
        if total < min_depth || var_alt[i] < min_alt_count {
            continue;
        }
        n_called += 1;

        // Position is heterozygous if ref reads are present and substantial
        let ref_frac = var_ref[i] as f64 / total as f64;
        if var_ref[i] >= 3 && ref_frac >= 0.15 {
            n_het += 1;
        }
    }

    // Hemizygous if we have called variants but none are heterozygous
    let is_hemizygous = n_called >= 3 && n_het == 0;

    VariantHetResult {
        n_het_variants: n_het,
        n_called_variants: n_called,
        is_hemizygous,
    }
}

/// Suggest a modified CNV tag if hemizygosity is detected.
pub fn suggest_cnv_with_deletion(
    het_result: &HetCheckResult,
    current_cnvtag: &str,
    _total_cn: u32,
) -> Option<String> {
    if !het_result.is_hemizygous {
        return None;
    }

    if current_cnvtag.contains("star5") {
        return None;
    }

    log::info!(
        "het_check: hemizygous sample with CNV tag '{}' — possible deletion",
        current_cnvtag
    );

    None
}

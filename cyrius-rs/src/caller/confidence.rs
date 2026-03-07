//! Confidence scoring for CYP2D6 star allele calls.
//!
//! Computes a confidence score (0.0–1.0) based on multiple independent
//! quality signals. Each signal contributes a penalty if it indicates
//! lower reliability.

use std::collections::{HashMap, HashSet};

use crate::types::StarCombinations;

/// Confidence score result with breakdown.
#[derive(Debug, Clone)]
pub struct ConfidenceScore {
    /// Overall confidence score (0.0–1.0)
    pub score: f64,
    /// Qualitative label: HIGH, MEDIUM, LOW
    pub label: String,
    /// Individual component scores for transparency
    pub components: ConfidenceComponents,
}

/// Individual components contributing to confidence.
#[derive(Debug, Clone)]
pub struct ConfidenceComponents {
    /// Depth quality (0–1): penalised by low depth or high MAD
    pub depth_quality: f64,
    /// CN rounding quality (0–1): how cleanly raw CN rounds to integer
    pub cn_quality: f64,
    /// Match quality (0–1): unique_match=1.0, more_than_one=0.6, fuzzy=0.4, no_match=0.0
    pub match_quality: f64,
    /// Variant completeness (0–1): fraction of expected variants observed
    pub variant_completeness: f64,
    /// Variant specificity (0–1): penalised by unexplained ALT observations
    pub variant_specificity: f64,
    /// SNP consistency (0–1): how consistent D6/D7 SNP ratios are within regions
    pub snp_consistency: f64,
}

/// Input data needed for confidence calculation.
pub struct ConfidenceInput<'a> {
    pub coverage_mad: f64,
    pub median_depth: f64,
    pub total_cn_raw: f64,
    pub spacer_cn_raw: f64,
    pub call_info: Option<&'a str>,
    pub filter: Option<&'a str>,
    pub genotype: Option<&'a str>,
    pub cnv_group: Option<&'a str>,
    pub d67_snp_raw: Option<&'a str>,
    pub variant_raw_count: Option<&'a indexmap::IndexMap<String, String>>,
    pub star_combinations: &'a StarCombinations,
}

/// Compute confidence score for a CYP2D6 call.
pub fn compute_confidence(input: &ConfidenceInput) -> ConfidenceScore {
    let depth_quality = score_depth(input.median_depth, input.coverage_mad);
    let cn_quality = score_cn(input.total_cn_raw, input.spacer_cn_raw);
    let match_quality = score_match(input.call_info, input.filter);
    let (variant_completeness, variant_specificity) = score_variants(
        input.genotype,
        input.variant_raw_count,
        &input.star_combinations.dstar,
    );
    let snp_consistency = score_snp_consistency(input.d67_snp_raw);

    // Weighted geometric mean — each component can independently drag down confidence
    let weights = [
        (depth_quality, 0.10),
        (cn_quality, 0.20),
        (match_quality, 0.25),
        (variant_completeness, 0.20),
        (variant_specificity, 0.10),
        (snp_consistency, 0.15),
    ];

    let weighted_log_sum: f64 = weights
        .iter()
        .map(|(score, weight)| weight * score.max(0.01).ln())
        .sum();
    let total_weight: f64 = weights.iter().map(|(_, w)| w).sum();
    let score = (weighted_log_sum / total_weight).exp().clamp(0.0, 1.0);

    let label = if score >= 0.85 {
        "HIGH"
    } else if score >= 0.5 {
        "MEDIUM"
    } else {
        "LOW"
    }
    .to_string();

    ConfidenceScore {
        score,
        label,
        components: ConfidenceComponents {
            depth_quality,
            cn_quality,
            match_quality,
            variant_completeness,
            variant_specificity,
            snp_consistency,
        },
    }
}

// ---------------------------------------------------------------------------
// Component scoring functions
// ---------------------------------------------------------------------------

/// Score depth quality: low depth or uneven coverage → lower confidence.
fn score_depth(median_depth: f64, coverage_mad: f64) -> f64 {
    // Depth: sigmoid from 0.5 at depth=10 to 1.0 at depth≥40
    let depth_score = if median_depth >= 40.0 {
        1.0
    } else if median_depth <= 5.0 {
        0.2
    } else {
        0.2 + 0.8 * ((median_depth - 5.0) / 35.0)
    };

    // MAD: penalise if > 0.15 (normal is ~0.05-0.10)
    let mad_score = if coverage_mad <= 0.10 {
        1.0
    } else if coverage_mad >= 0.30 {
        0.3
    } else {
        1.0 - 3.5 * (coverage_mad - 0.10)
    };

    depth_score * mad_score
}

/// Score CN rounding quality: how cleanly raw values round to integers.
fn score_cn(total_cn_raw: f64, spacer_cn_raw: f64) -> f64 {
    let total_frac = (total_cn_raw - total_cn_raw.round()).abs();
    let spacer_frac = (spacer_cn_raw - spacer_cn_raw.round()).abs();

    // Perfect rounding (frac < 0.1) → 1.0, poor (frac > 0.4) → 0.3
    let total_score = if total_frac <= 0.15 {
        1.0
    } else if total_frac >= 0.40 {
        0.3
    } else {
        1.0 - 2.8 * (total_frac - 0.15)
    };

    let spacer_score = if spacer_frac <= 0.15 {
        1.0
    } else if spacer_frac >= 0.40 {
        0.3
    } else {
        1.0 - 2.8 * (spacer_frac - 0.15)
    };

    (total_score + spacer_score) / 2.0
}

/// Score match quality from call_info and filter.
fn score_match(call_info: Option<&str>, filter: Option<&str>) -> f64 {
    let base = match call_info {
        Some("unique_match") => 1.0,
        Some("more_than_one_match") => 0.7,
        Some(s) if s.starts_with("fuzzy_match") => 0.4,
        Some("no_match") | None => 0.1,
        _ => 0.5,
    };

    // Penalty for non-PASS filters
    let filter_penalty = match filter {
        Some("PASS") => 1.0,
        Some("More_than_one_possible_genotype") => 0.7,
        Some("Fuzzy_match") => 0.6,
        Some("LowQ_high_CN") => 0.4,
        Some("Not_assigned_to_haplotypes") => 0.2,
        None => 0.5,
        _ => 0.5,
    };

    base * filter_penalty
}

/// Score variant completeness and specificity.
/// Returns (completeness, specificity).
fn score_variants(
    genotype: Option<&str>,
    raw_counts: Option<&indexmap::IndexMap<String, String>>,
    dstar: &HashMap<String, String>,
) -> (f64, f64) {
    let genotype = match genotype {
        Some(g) => g,
        None => return (0.1, 1.0),
    };
    let raw_counts = match raw_counts {
        Some(rc) => rc,
        None => return (0.5, 1.0),
    };

    // Get the first diplotype (before ';')
    let first_diplotype = genotype.split(';').next().unwrap_or(genotype);

    // Extract expected variants for the called alleles
    let mut expected_variants: HashSet<String> = HashSet::new();
    for hap in first_diplotype.split('/') {
        for component in hap.split('+') {
            // Strip xN suffix
            let base = if let Some(pos) = component.find('x') {
                let suffix = &component[pos + 1..];
                if !suffix.is_empty() && suffix.chars().all(|c| c.is_ascii_digit()) {
                    &component[..pos]
                } else {
                    component
                }
            } else {
                component
            };
            let lookup = if base.starts_with('*') {
                base.to_string()
            } else {
                format!("*{}", base)
            };
            if let Some(var_key) = dstar.get(&lookup) {
                if var_key != "NA" {
                    for v in var_key.split('_') {
                        if v != "exon9gc" {
                            expected_variants.insert(v.to_string());
                        }
                    }
                }
            }
        }
    }

    if expected_variants.is_empty() {
        // *1/*1 or similar — no defining variants expected
        return (1.0, 1.0);
    }

    // Check how many expected variants are actually observed (ALT count > 0)
    let mut found = 0;
    let mut total = 0;
    for var in &expected_variants {
        total += 1;
        if let Some(count_str) = raw_counts.get(var) {
            let alt_count = parse_alt_count(count_str);
            if alt_count > 0 {
                found += 1;
            }
        }
    }

    let completeness = if total > 0 {
        found as f64 / total as f64
    } else {
        1.0
    };

    // Specificity: check for unexplained ALT variants
    // Count variants with significant ALT reads that aren't in expected set
    let mut unexplained = 0;
    let mut total_variants_checked = 0;
    for (var, count_str) in raw_counts {
        let alt_count = parse_alt_count(count_str);
        let ref_count = parse_ref_count(count_str);
        let total = alt_count + ref_count;
        if total < 10 {
            continue;
        }
        total_variants_checked += 1;
        let af = alt_count as f64 / total as f64;
        // Significant ALT (> 15% AF) that's not expected
        if af > 0.15 && !expected_variants.contains(var) {
            unexplained += 1;
        }
    }

    let specificity = if total_variants_checked > 0 {
        let unexplained_frac = unexplained as f64 / total_variants_checked as f64;
        (1.0 - unexplained_frac * 5.0).max(0.2)
    } else {
        1.0
    };

    (completeness, specificity)
}

/// Score consistency of D6/D7 SNP ratios.
fn score_snp_consistency(d67_snp_raw: Option<&str>) -> f64 {
    let raw = match d67_snp_raw {
        Some(r) => r,
        None => return 0.5,
    };

    let values: Vec<f64> = raw
        .split(',')
        .filter_map(|s| s.trim().parse::<f64>().ok())
        .collect();

    if values.len() < 10 {
        return 0.5;
    }

    // Check how consistent ratios are within expected CN regions.
    // Good calls have tight clusters at integer values.
    // Score: average deviation from nearest integer across all sites.
    let total_dev: f64 = values
        .iter()
        .map(|v| (v - v.round()).abs())
        .sum();
    let avg_dev = total_dev / values.len() as f64;

    // avg_dev < 0.15 → great, > 0.35 → poor
    if avg_dev <= 0.10 {
        1.0
    } else if avg_dev >= 0.35 {
        0.3
    } else {
        1.0 - 2.8 * (avg_dev - 0.10)
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Parse ALT count from raw_count format "alt,ref" or "alt(fwd:rev),ref"
fn parse_alt_count(s: &str) -> u32 {
    let parts: Vec<&str> = s.split(',').collect();
    if parts.is_empty() {
        return 0;
    }
    let alt_str = parts[0];
    // Handle "0(0:0)" format
    let base = if let Some(paren) = alt_str.find('(') {
        &alt_str[..paren]
    } else {
        alt_str
    };
    base.parse::<u32>().unwrap_or(0)
}

/// Parse REF count from raw_count format "alt,ref"
fn parse_ref_count(s: &str) -> u32 {
    let parts: Vec<&str> = s.split(',').collect();
    if parts.len() < 2 {
        return 0;
    }
    parts[1].trim().parse::<u32>().unwrap_or(0)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_depth_scoring() {
        assert!((score_depth(60.0, 0.05) - 1.0).abs() < 0.01);
        assert!(score_depth(10.0, 0.05) < 0.5);
        assert!(score_depth(60.0, 0.25) < 0.7);
        assert!(score_depth(5.0, 0.30) < 0.2);
    }

    #[test]
    fn test_cn_scoring() {
        assert!(score_cn(4.05, 2.02) > 0.9);
        assert!(score_cn(4.45, 2.48) < 0.5);
    }

    #[test]
    fn test_match_scoring() {
        assert!((score_match(Some("unique_match"), Some("PASS")) - 1.0).abs() < 0.01);
        assert!(score_match(Some("no_match"), None) < 0.1);
        assert!(score_match(Some("unique_match"), Some("More_than_one_possible_genotype")) < 0.8);
    }

    #[test]
    fn test_snp_consistency() {
        // All perfect integers
        let perfect = "2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3";
        assert!(score_snp_consistency(Some(perfect)) > 0.9);

        // Noisy
        let noisy = "2.4,1.6,2.3,1.7,2.5,1.8,2.2,2.6,1.5,2.8,3.4,2.7,3.3,2.8,3.5,2.6,3.2,3.4,2.7,3.3";
        assert!(score_snp_consistency(Some(noisy)) < 0.7);
    }

    #[test]
    fn test_parse_alt_count() {
        assert_eq!(parse_alt_count("26,0"), 26);
        assert_eq!(parse_alt_count("0(0:0),81"), 0);
        assert_eq!(parse_alt_count("1(1:0),106"), 1);
        assert_eq!(parse_alt_count("0,163"), 0);
    }
}

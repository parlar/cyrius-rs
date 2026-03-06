//! Post-hoc consistency checks on star allele calls.
//!
//! These checks validate the final call against independent evidence:
//! 1. Gene conversion map: spatial D6/D7 pattern should match CNV group
//! 2. Allele mismatch rate: called allele's defining variants should be observed
//! 3. Allele balance: het variant VAFs should match expected for the CN

use crate::depth_calling::snp_count::CrossingStats;
use crate::types::StarCombinations;

/// Results of all consistency checks.
#[derive(Debug, Clone)]
pub struct ConsistencyResult {
    /// Gene conversion map: chimeric pattern in a non-hybrid, or uniform in a hybrid
    pub conversion_map_flag: Option<String>,
    /// Fraction of called allele defining variants that contradict observations
    pub mismatch_rate: Option<f64>,
    pub mismatch_flag: Option<String>,
    /// Mean VAF deviation from expected at defining variant positions
    pub balance_deviation: Option<f64>,
    pub balance_flag: Option<String>,
    /// Read-pair crossing fraction inconsistency
    pub crossing_flag: Option<String>,
}

// --- Gene Conversion Map ---

const WINDOW_SIZE: usize = 5;
const MIN_DEPTH: usize = 10;
const D6_THRESHOLD: f64 = 0.35;
const D7_THRESHOLD: f64 = 0.15;
const MIN_WINDOWS_FOR_CALL: usize = 3;

/// Check if the spatial D6/D7 pattern across the gene matches the CNV group.
pub fn check_conversion_map(
    snp_d6: &[usize],
    snp_d7: &[usize],
    cnvtag: &str,
) -> Option<String> {
    let n = snp_d6.len().min(snp_d7.len());
    if n < WINDOW_SIZE {
        return None;
    }

    // Compute per-position ratios
    let ratios: Vec<Option<f64>> = (0..n)
        .map(|i| {
            let total = snp_d6[i] + snp_d7[i];
            if total >= MIN_DEPTH {
                Some(snp_d6[i] as f64 / total as f64)
            } else {
                None
            }
        })
        .collect();

    // Sliding window classification
    let mut d6_windows = 0usize;
    let mut d7_windows = 0usize;
    let mut total_windows = 0usize;

    for start in 0..=n.saturating_sub(WINDOW_SIZE) {
        let window = &ratios[start..start + WINDOW_SIZE];
        let valid: Vec<f64> = window.iter().filter_map(|&r| r).collect();
        if valid.len() < 3 {
            continue;
        }
        total_windows += 1;
        let mean = valid.iter().sum::<f64>() / valid.len() as f64;
        if mean > D6_THRESHOLD {
            d6_windows += 1;
        } else if mean < D7_THRESHOLD {
            d7_windows += 1;
        }
    }

    if total_windows < MIN_WINDOWS_FOR_CALL {
        return None;
    }

    let is_chimeric = d6_windows >= MIN_WINDOWS_FOR_CALL && d7_windows >= MIN_WINDOWS_FOR_CALL;

    if cnvtag == "cn2" && is_chimeric {
        log::info!(
            "conversion_map: chimeric pattern in cn2 (D6={}, D7={}, total={})",
            d6_windows, d7_windows, total_windows,
        );
        return Some("Conversion_map_chimeric_in_cn2".to_string());
    }

    // Note: "uniform in hybrid" check disabled — with pileup-based ratios (not k-mer),
    // normal D6 copies keep the ratio above D7_THRESHOLD even in hybrid samples.
    // The changepoint module handles hybrid detection better for this data.

    None
}

// --- Allele Mismatch Rate ---

const MISMATCH_VAF_LOW: f64 = 0.05;
const MISMATCH_RATE_THRESHOLD: f64 = 0.50;

/// Check if called star alleles' defining variants match observations.
/// Returns (mismatch_rate, flag) if the called allele has contradictions.
pub fn check_allele_mismatch(
    genotype: &str,
    var_alt: &[usize],
    var_ref: &[usize],
    var_list: &[String],
    star_combinations: &StarCombinations,
) -> (Option<f64>, Option<String>) {
    // Parse alleles from genotype like "*1/*2" or "*36+*10/*41"
    let alleles = parse_alleles(genotype);
    if alleles.is_empty() {
        return (None, None);
    }

    let mut n_checked = 0usize;
    let mut n_contradictions = 0usize;

    for allele in &alleles {
        // Look up defining variants for this allele
        let variant_key = match star_combinations.dstar.get(allele.as_str()) {
            Some(k) => k,
            None => continue,
        };

        // variant_key is like "g.42128945C>T_g.42129033G>A" or "NA" for *1
        if variant_key == "NA" {
            continue;
        }

        let defining_vars: Vec<&str> = variant_key.split('_').collect();

        for dvar in &defining_vars {
            // Find this variant in var_list
            if let Some(idx) = var_list.iter().position(|v| v == dvar) {
                if idx < var_alt.len() && idx < var_ref.len() {
                    let alt = var_alt[idx];
                    let ref_count = var_ref[idx];
                    let total = alt + ref_count;
                    if total < 5 {
                        continue; // insufficient depth
                    }
                    n_checked += 1;
                    let vaf = alt as f64 / total as f64;
                    // Allele expects this variant to be present
                    if vaf < MISMATCH_VAF_LOW {
                        n_contradictions += 1;
                    }
                }
            }
        }
    }

    if n_checked == 0 {
        return (None, None);
    }

    let rate = n_contradictions as f64 / n_checked as f64;

    if rate > MISMATCH_RATE_THRESHOLD {
        log::info!(
            "allele_mismatch: {}/{} defining variants contradicted (rate={:.2}) for genotype '{}'",
            n_contradictions, n_checked, rate, genotype,
        );
        (
            Some(rate),
            Some(format!("Allele_mismatch_rate_{:.0}pct", rate * 100.0)),
        )
    } else {
        (Some(rate), None)
    }
}

/// Parse individual star allele names from a diplotype string.
/// "*1/*36+*10" → ["*1", "*36", "*10"]
/// "*2x2/*41" → ["*2", "*41"]
fn parse_alleles(genotype: &str) -> Vec<String> {
    let mut alleles = Vec::new();
    // Split on '/' for haplotypes, then '+' for tandems
    for hap in genotype.split('/') {
        for tandem_part in hap.split('+') {
            // Strip duplication markers like "x2", "x3"
            let base = if let Some(pos) = tandem_part.find('x') {
                let suffix = &tandem_part[pos + 1..];
                if !suffix.is_empty() && suffix.chars().all(|c| c.is_ascii_digit()) {
                    &tandem_part[..pos]
                } else {
                    tandem_part
                }
            } else {
                tandem_part
            };
            let trimmed = base.trim();
            if !trimmed.is_empty() {
                alleles.push(trimmed.to_string());
            }
        }
    }
    alleles
}

// --- Allele Balance ---

const BALANCE_MIN_DEPTH: usize = 10;
const BALANCE_DEVIATION_THRESHOLD: f64 = 0.30;

/// Check VAF balance at defining variant positions for cn2 samples.
/// For cn2 het sites, expected VAF ~0.5. Returns mean absolute deviation.
pub fn check_allele_balance(
    genotype: &str,
    var_alt: &[usize],
    var_ref: &[usize],
    var_list: &[String],
    star_combinations: &StarCombinations,
    cnvtag: &str,
    total_cn: u32,
) -> (Option<f64>, Option<String>) {
    // Only meaningful for cn2 (diploid, no structural variants)
    if cnvtag != "cn2" || total_cn != 4 {
        return (None, None);
    }

    // Collect defining variants per haplotype to identify het vs hom sites
    let haplotypes: Vec<&str> = genotype.split('/').collect();
    let mut hap_variants: Vec<std::collections::HashSet<String>> = Vec::new();
    for hap in &haplotypes {
        let mut vars = std::collections::HashSet::new();
        for tandem_part in hap.split('+') {
            let allele_name = if let Some(pos) = tandem_part.find('x') {
                let suffix = &tandem_part[pos + 1..];
                if !suffix.is_empty() && suffix.chars().all(|c| c.is_ascii_digit()) {
                    tandem_part[..pos].trim()
                } else {
                    tandem_part.trim()
                }
            } else {
                tandem_part.trim()
            };
            if let Some(variant_key) = star_combinations.dstar.get(allele_name) {
                if variant_key != "NA" {
                    for v in variant_key.split('_') {
                        vars.insert(v.to_string());
                    }
                }
            }
        }
        hap_variants.push(vars);
    }

    // Only check variants that are het (present in exactly one haplotype)
    let het_variants: std::collections::HashSet<String> = if hap_variants.len() == 2 {
        let only_hap0 = hap_variants[0].difference(&hap_variants[1]).cloned();
        let only_hap1 = hap_variants[1].difference(&hap_variants[0]).cloned();
        only_hap0.chain(only_hap1).collect()
    } else {
        // Can't determine het/hom with non-diploid haplotypes
        return (None, None);
    };

    // For each het defining variant, compute VAF and deviation from 0.5
    let mut deviations = Vec::new();
    let expected_vaf = 0.5;

    for dvar in &het_variants {
        if let Some(idx) = var_list.iter().position(|v| v == dvar) {
            if idx < var_alt.len() && idx < var_ref.len() {
                let alt = var_alt[idx];
                let ref_count = var_ref[idx];
                let total = alt + ref_count;
                if total < BALANCE_MIN_DEPTH {
                    continue;
                }
                let vaf = alt as f64 / total as f64;
                // Only check het sites (VAF between 0.05 and 0.95)
                if vaf > 0.05 && vaf < 0.95 {
                    deviations.push((vaf - expected_vaf).abs());
                }
            }
        }
    }

    if deviations.len() < 2 {
        return (None, None);
    }

    let mean_dev = deviations.iter().sum::<f64>() / deviations.len() as f64;

    if mean_dev > BALANCE_DEVIATION_THRESHOLD {
        log::info!(
            "allele_balance: mean VAF deviation={:.3} ({} sites) for genotype '{}' — \
             expected ~0.5 for cn2",
            mean_dev, deviations.len(), genotype,
        );
        (
            Some(mean_dev),
            Some(format!("Allele_balance_deviation_{:.0}pct", mean_dev * 100.0)),
        )
    } else {
        (Some(mean_dev), None)
    }
}

// --- Read-Pair Crossing Fraction ---

/// Expected crossing fractions (from pdx_caller calibration.py defaults)
const CROSSING_NORMAL: f64 = 0.02;
const CROSSING_HYBRID: f64 = 0.10;
const CROSSING_MIN_READS: usize = 50;

/// Check if read-pair crossing fraction is consistent with CNV group.
/// Note: With pileup-based paralog classification at sparse diagnostic SNP positions,
/// crossing fraction is ~0.001-0.005 for ALL samples (normal and hybrid alike).
/// This check is therefore informational only — it logs the data but does not flag.
/// Effective crossing detection requires dense k-mer classification (as in pdx_caller).
pub fn check_crossing_fraction(
    crossing: &CrossingStats,
    cnvtag: &str,
) -> Option<String> {
    let total = crossing.n_consistent + crossing.n_crossing;
    if total >= CROSSING_MIN_READS {
        log::info!(
            "crossing_fraction: {:.4} ({}/{}) for CNV group '{}'",
            crossing.crossing_fraction, crossing.n_crossing, total, cnvtag,
        );
    }
    // Informational only — pileup-based crossing does not distinguish hybrids from normals
    None
}

/// Run all consistency checks and return combined result.
pub fn run_consistency_checks(
    snp_d6: &[usize],
    snp_d7: &[usize],
    var_alt: &[usize],
    var_ref: &[usize],
    var_list: &[String],
    cnvtag: &str,
    total_cn: u32,
    genotype: Option<&str>,
    star_combinations: &StarCombinations,
    crossing: &CrossingStats,
) -> ConsistencyResult {
    let conversion_map_flag = check_conversion_map(snp_d6, snp_d7, cnvtag);

    let (mismatch_rate, mismatch_flag) = if let Some(geno) = genotype {
        check_allele_mismatch(geno, var_alt, var_ref, var_list, star_combinations)
    } else {
        (None, None)
    };

    let (balance_deviation, balance_flag) = if let Some(geno) = genotype {
        check_allele_balance(geno, var_alt, var_ref, var_list, star_combinations, cnvtag, total_cn)
    } else {
        (None, None)
    };

    let crossing_flag = check_crossing_fraction(crossing, cnvtag);

    ConsistencyResult {
        conversion_map_flag,
        mismatch_rate,
        mismatch_flag,
        balance_deviation,
        balance_flag,
        crossing_flag,
    }
}

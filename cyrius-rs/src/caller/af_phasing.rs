use crate::types::StarCombinations;

/// Result of AF-based duplication phasing for a single het variant.
#[derive(Debug, Clone)]
pub struct AfPhasingResult {
    /// The het variant used for phasing (e.g., "g.42127941G>A")
    pub variant: String,
    /// Allele frequency of the variant
    pub af: f64,
    /// Which star allele carries this variant (if identified)
    pub carrier_allele: Option<String>,
    /// Estimated copy number of the carrier allele
    pub carrier_copies: u32,
    /// Estimated copy number of the other allele
    pub other_copies: u32,
    /// Total CN used
    pub total_cn: u32,
}

/// Overall AF phasing summary for a sample.
#[derive(Debug, Clone)]
pub struct AfPhasingSummary {
    /// Individual phasing results for each informative het variant
    pub results: Vec<AfPhasingResult>,
    /// Consensus phasing if multiple variants agree
    pub consensus: Option<String>,
}

/// Identify het variants and estimate per-allele copy numbers.
///
/// For CN≥3 with two candidate alleles, het variant AFs reveal which
/// allele is duplicated. If a variant has AF≈0.33 at CN=3, its carrier
/// allele has 1 copy (1/3). If AF≈0.67, carrier has 2 copies.
///
/// Algorithm (from StellarPGx dup_test_cn_3_4):
/// 1. Find het variants (0.1 < AF < 0.9)
/// 2. For each het variant, compute carrier_copies = round(AF × total_CN)
/// 3. Look up which star allele carries the variant
/// 4. Report per-allele copy number estimates
pub fn estimate_allele_copies(
    var_alt: &[usize],
    var_ref: &[usize],
    var_list: &[String],
    total_cn: u32,
    genotype: Option<&str>,
    star_combinations: &StarCombinations,
) -> AfPhasingSummary {
    if total_cn < 3 {
        return AfPhasingSummary {
            results: Vec::new(),
            consensus: None,
        };
    }

    // Parse genotype into two alleles
    let alleles: Vec<&str> = match genotype {
        Some(g) if g.contains('/') => g.split('/').collect(),
        _ => {
            return AfPhasingSummary {
                results: Vec::new(),
                consensus: None,
            }
        }
    };
    if alleles.len() != 2 {
        return AfPhasingSummary {
            results: Vec::new(),
            consensus: None,
        };
    }

    // Get variant keys for each allele from dstar.
    // For tandem alleles like "*36+*10", collect variants from ALL sub-alleles.
    let vars_a = collect_allele_variants(alleles[0], star_combinations);
    let vars_b = collect_allele_variants(alleles[1], star_combinations);

    let mut results = Vec::new();

    for (i, var_name) in var_list.iter().enumerate() {
        if i >= var_alt.len() || i >= var_ref.len() {
            break;
        }
        let alt = var_alt[i] as f64;
        let reff = var_ref[i] as f64;
        let total = alt + reff;
        if total < 10.0 {
            continue; // Skip low-depth sites
        }

        let af = alt / total;
        // Only consider het variants
        if af < 0.1 || af > 0.9 {
            continue;
        }

        // Determine which allele carries this variant
        let in_a = vars_a.iter().any(|v| v == var_name);
        let in_b = vars_b.iter().any(|v| v == var_name);

        // Only useful if variant is in exactly one allele
        let carrier = if in_a && !in_b {
            Some(alleles[0].to_string())
        } else if in_b && !in_a {
            Some(alleles[1].to_string())
        } else {
            None
        };

        let carrier_copies = (af * total_cn as f64).round() as u32;
        let other_copies = total_cn.saturating_sub(carrier_copies);

        results.push(AfPhasingResult {
            variant: var_name.clone(),
            af,
            carrier_allele: carrier,
            carrier_copies,
            other_copies,
            total_cn,
        });
    }

    // Build consensus: count how many variants agree on a phasing
    let consensus = build_consensus(&results, alleles[0], alleles[1], total_cn);

    AfPhasingSummary { results, consensus }
}

/// Strip copy number suffix (e.g., "*10x2" → "*10").
fn strip_copy_suffix(allele: &str) -> String {
    if let Some(pos) = allele.rfind('x') {
        if allele[pos + 1..].parse::<u32>().is_ok() {
            return allele[..pos].to_string();
        }
    }
    allele.to_string()
}

/// Collect all variant keys for an allele, handling tandem notation.
/// For "*36+*10", collects variants from both "*36" and "*10".
/// For "*10x2", strips the "x2" suffix before lookup.
fn collect_allele_variants(allele: &str, star_combinations: &StarCombinations) -> Vec<String> {
    let mut vars = Vec::new();
    // Split tandem notation (e.g., "*36+*10" → ["*36", "*10"])
    for sub_allele in allele.split('+') {
        let key = strip_copy_suffix(sub_allele);
        if let Some(var_str) = star_combinations.dstar.get(&key) {
            for v in var_str.split('_') {
                if !vars.contains(&v.to_string()) {
                    vars.push(v.to_string());
                }
            }
        }
    }
    vars
}

/// Build consensus description from individual phasing results.
fn build_consensus(results: &[AfPhasingResult], allele_a: &str, allele_b: &str, total_cn: u32) -> Option<String> {
    // Count votes for each (carrier, copies) assignment
    let mut votes_a: Vec<u32> = Vec::new();
    let mut votes_b: Vec<u32> = Vec::new();

    for r in results {
        match r.carrier_allele.as_deref() {
            Some(a) if a == allele_a => votes_a.push(r.carrier_copies),
            Some(b) if b == allele_b => votes_b.push(r.carrier_copies),
            _ => {}
        }
    }

    if votes_a.is_empty() && votes_b.is_empty() {
        return None;
    }

    // Find modal copy number for each allele
    let mode_a = mode(&votes_a);
    let mode_b = mode(&votes_b);

    // Check consistency
    match (mode_a, mode_b) {
        (Some((copies_a, count_a)), Some((copies_b, count_b))) => {
            if copies_a + copies_b == total_cn && count_a >= 2 && count_b >= 2 {
                Some(format!(
                    "{}x{}/{}x{} (a_votes={},b_votes={})",
                    allele_a, copies_a, allele_b, copies_b, count_a, count_b,
                ))
            } else {
                Some(format!(
                    "inconsistent: {}x{}(votes={})/{}x{}(votes={}) total_cn={}",
                    allele_a, copies_a, count_a, allele_b, copies_b, count_b,
                    total_cn,
                ))
            }
        }
        (Some((copies_a, count_a)), None) => {
            let copies_b = total_cn.saturating_sub(copies_a);
            if count_a >= 2 {
                Some(format!(
                    "{}x{}/{}x{} (a_votes={})",
                    allele_a, copies_a, allele_b, copies_b, count_a,
                ))
            } else {
                None
            }
        }
        (None, Some((copies_b, count_b))) => {
            let copies_a = total_cn.saturating_sub(copies_b);
            if count_b >= 2 {
                Some(format!(
                    "{}x{}/{}x{} (b_votes={})",
                    allele_a, copies_a, allele_b, copies_b, count_b,
                ))
            } else {
                None
            }
        }
        _ => None,
    }
}

/// Find the mode (most common value) and its count.
fn mode(values: &[u32]) -> Option<(u32, usize)> {
    if values.is_empty() {
        return None;
    }
    let mut counts = std::collections::HashMap::new();
    for &v in values {
        *counts.entry(v).or_insert(0usize) += 1;
    }
    counts.into_iter().max_by_key(|&(_, count)| count)
}

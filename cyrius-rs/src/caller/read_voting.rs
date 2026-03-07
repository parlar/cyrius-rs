//! Read-level allele voting for independent QC of star allele calls.
//!
//! Aldy-style approach: for each read overlapping target variant positions,
//! score how well it matches each star allele's expected variant profile.
//! Aggregate into diplotype-level likelihoods and compare with the pileup-based call.

use std::collections::{HashMap, HashSet};
use rust_htslib::bam::{self, Read as BamRead};
use crate::types::StarCombinations;

/// A single variant position with expected ref and alt bases.
#[derive(Debug, Clone)]
pub struct VariantSite {
    pub genomic_pos: i64,
    pub ref_base: u8,
    /// For SNPs: the alt base. For deletions: 0 (meaning "base absent").
    pub alt_base: u8,
    pub name: String,
    /// Index in the original var_list (for cross-referencing with pileup counts)
    pub var_list_index: usize,
    /// True if this is a single-base deletion (alt = base absent)
    pub is_deletion: bool,
}

/// Per-read observations at variant positions.
#[derive(Debug)]
pub struct ReadProfile {
    pub name: String,
    /// (variant_site_index, observed_base, base_quality)
    pub observations: Vec<(usize, u8, u8)>,
}

/// Result of graph-based validation.
#[derive(Debug, Clone)]
pub struct VotingResult {
    /// Top diplotypes with scores: (allele1, allele2, log10_likelihood)
    pub top_diplotypes: Vec<(String, String, f64)>,
    /// Number of informative reads (covering at least one variant site)
    pub n_informative_reads: usize,
    /// Whether the pileup-based call appears in the top N
    pub call_rank: Option<usize>,
    /// Score gap between best diplotype and the called diplotype (if found)
    pub call_score_gap: Option<f64>,
    /// Flag if call is not in top N
    pub flag: Option<String>,
}

/// Parse variant sites from the var_list.
/// Handles single-base substitutions and single-base deletions.
pub fn parse_variant_sites(var_list: &[String]) -> Vec<VariantSite> {
    let mut sites = Vec::new();
    for (idx, name) in var_list.iter().enumerate() {
        if let Some(mut site) = parse_variant(name) {
            site.var_list_index = idx;
            sites.push(site);
        }
    }
    sites
}

fn parse_variant(name: &str) -> Option<VariantSite> {
    let s = name.strip_prefix("g.")?;

    // Find where digits end
    let digit_end = s.find(|c: char| !c.is_ascii_digit())?;
    let pos: i64 = s[..digit_end].parse().ok()?;
    let rest = &s[digit_end..];

    // Single-base SNP: "G>A"
    if rest.len() == 3 && rest.as_bytes()[1] == b'>' {
        let ref_base = rest.as_bytes()[0];
        let alt_base = rest.as_bytes()[2];
        if ref_base.is_ascii_alphabetic() && alt_base.is_ascii_alphabetic() {
            return Some(VariantSite {
                genomic_pos: pos,
                ref_base: ref_base.to_ascii_uppercase(),
                alt_base: alt_base.to_ascii_uppercase(),
                name: name.to_string(),
                var_list_index: 0,
                is_deletion: false,
            });
        }
    }

    // Deletion: "delT", "delCTT", etc.
    // For multi-base deletions, use the first deleted base position
    if rest.starts_with("del") {
        let del_seq = &rest[3..];
        if !del_seq.is_empty() && del_seq.chars().all(|c| c.is_ascii_alphabetic()) {
            let first_base = del_seq.as_bytes()[0];
            return Some(VariantSite {
                genomic_pos: pos,
                ref_base: first_base.to_ascii_uppercase(),
                alt_base: 0, // absent
                name: name.to_string(),
                var_list_index: 0,
                is_deletion: true,
            });
        }
    }

    None
}

/// Build allele definitions: for each star allele, the set of variant site indices it carries.
pub fn build_allele_defs(
    dstar: &HashMap<String, String>,
    sites: &[VariantSite],
) -> Vec<(String, HashSet<usize>)> {
    let name_to_idx: HashMap<&str, usize> = sites
        .iter()
        .enumerate()
        .map(|(i, s)| (s.name.as_str(), i))
        .collect();

    let mut defs: Vec<(String, HashSet<usize>)> = dstar
        .iter()
        .map(|(allele_name, variant_key)| {
            let mut indices = HashSet::new();
            if variant_key != "NA" {
                for var_name in variant_key.split('_') {
                    if let Some(&idx) = name_to_idx.get(var_name) {
                        indices.insert(idx);
                    }
                }
            }
            (allele_name.clone(), indices)
        })
        .collect();

    defs.sort_by(|a, b| a.0.cmp(&b.0));
    defs
}

/// Collect read profiles from BAM across the CYP2D6 variant region.
pub fn collect_read_profiles(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    region_start: i64,
    region_end: i64,
    sites: &[VariantSite],
) -> Vec<ReadProfile> {
    let tid = match reader.header().tid(nchr.as_bytes()) {
        Some(t) => t,
        None => {
            log::warn!("read_voting: chromosome '{}' not found in BAM header", nchr);
            return Vec::new();
        }
    };

    // Genomic position (1-based) → site index
    let pos_to_site: HashMap<i64, usize> = sites
        .iter()
        .enumerate()
        .map(|(i, s)| (s.genomic_pos, i))
        .collect();

    log::info!(
        "read_voting: fetching {}:{}-{} (tid={}, {} variant positions)",
        nchr, region_start, region_end, tid, pos_to_site.len(),
    );

    let mut profiles: HashMap<String, Vec<(usize, u8, u8)>> = HashMap::new();

    if reader
        .fetch(bam::FetchDefinition::Region(
            tid as i32,
            region_start,
            region_end,
        ))
        .is_err()
    {
        log::warn!("read_voting: fetch failed for region");
        return Vec::new();
    }

    let mut n_reads_total = 0usize;
    let mut n_reads_passing = 0usize;

    for result in reader.records() {
        let mut record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };
        n_reads_total += 1;

        if record.is_secondary()
            || record.is_supplementary()
            || record.is_duplicate()
            || record.is_unmapped()
        {
            continue;
        }
        if record.mapq() < 10 {
            continue;
        }
        n_reads_passing += 1;

        let ref_start = record.pos(); // 0-based inclusive
        record.cache_cigar();
        let cigar = match record.cigar_cached() {
            Some(c) => c,
            None => continue,
        };
        let ref_end = cigar.end_pos(); // 0-based exclusive

        let read_name = String::from_utf8_lossy(record.qname()).to_string();
        let seq = record.seq().as_bytes();
        let qual = record.qual().to_vec();

        for (&gpos_1based, &site_idx) in &pos_to_site {
            let ref_pos_0based = gpos_1based - 1;

            if ref_pos_0based < ref_start || ref_pos_0based >= ref_end {
                continue;
            }

            // Map reference position to read position via CIGAR
            let read_pos_result = cigar.read_pos(ref_pos_0based as u32, false, false).ok().flatten();

            if sites[site_idx].is_deletion {
                // For deletion sites: None means base IS deleted (carrier),
                // Some means base is present (non-carrier)
                let (base, q) = match read_pos_result {
                    None => (0u8, 30u8), // deleted = matches the deletion allele
                    Some(rp) => {
                        let rp = rp as usize;
                        if rp < seq.len() {
                            (seq[rp], if rp < qual.len() { qual[rp] } else { 30 })
                        } else {
                            continue;
                        }
                    }
                };
                profiles
                    .entry(read_name.clone())
                    .or_default()
                    .push((site_idx, base, q));
            } else {
                // For SNP sites: need a valid read position
                if let Some(read_pos) = read_pos_result {
                    let rp = read_pos as usize;
                    if rp < seq.len() {
                        let base = seq[rp];
                        let q = if rp < qual.len() { qual[rp] } else { 30 };
                        profiles
                            .entry(read_name.clone())
                            .or_default()
                            .push((site_idx, base, q));
                    }
                }
            }
        }
    }

    log::info!(
        "read_voting: {} total reads, {} passing filters, {} with variant observations",
        n_reads_total, n_reads_passing, profiles.len(),
    );

    profiles
        .into_iter()
        .map(|(name, observations)| ReadProfile { name, observations })
        .filter(|p| !p.observations.is_empty())
        .collect()
}

/// Log10-likelihood of a read under a given allele hypothesis.
fn score_read_allele(
    profile: &ReadProfile,
    allele_variants: &HashSet<usize>,
    sites: &[VariantSite],
) -> f64 {
    let mut score = 0.0;

    for &(site_idx, observed_base, qual) in &profile.observations {
        let site = &sites[site_idx];
        let allele_carries_variant = allele_variants.contains(&site_idx);

        if site.is_deletion {
            // For deletion sites:
            // observed_base == 0 means the base is deleted in the read
            // observed_base != 0 means the base is present
            let read_has_deletion = observed_base == 0;
            let expect_deletion = allele_carries_variant;

            let p_err = 10.0_f64.powf(-(qual as f64) / 10.0).clamp(1e-10, 0.5);
            if read_has_deletion == expect_deletion {
                score += (1.0 - p_err).log10();
            } else {
                score += (p_err / 3.0).log10();
            }
        } else {
            // SNP site
            let expected_base = if allele_carries_variant {
                site.alt_base
            } else {
                site.ref_base
            };

            let p_err = 10.0_f64.powf(-(qual as f64) / 10.0).clamp(1e-10, 0.5);

            if observed_base.to_ascii_uppercase() == expected_base {
                score += (1.0 - p_err).log10();
            } else {
                score += (p_err / 3.0).log10();
            }
        }
    }

    score
}

/// Find the best diplotypes by scoring all allele pairs.
/// Returns top_n diplotypes sorted by score (descending).
pub fn find_best_diplotypes(
    profiles: &[ReadProfile],
    allele_defs: &[(String, HashSet<usize>)],
    sites: &[VariantSite],
    top_n: usize,
) -> Vec<(String, String, f64)> {
    let n_alleles = allele_defs.len();
    let n_reads = profiles.len();

    if n_reads == 0 || n_alleles == 0 {
        return Vec::new();
    }

    // Precompute per-read per-allele log10-likelihoods
    let read_scores: Vec<Vec<f64>> = profiles
        .iter()
        .map(|profile| {
            allele_defs
                .iter()
                .map(|(_, variants)| score_read_allele(profile, variants, sites))
                .collect()
        })
        .collect();

    // Score each diplotype pair
    // P(read | diplotype i,j) = 0.5 * P(read|i) + 0.5 * P(read|j)
    // log10 P = log10(0.5) + log10(10^s_i + 10^s_j)
    //         = -0.301 + max(s_i,s_j) + log10(1 + 10^(-|s_i-s_j|))
    let log10_half = 0.5_f64.log10();

    let mut best: Vec<(String, String, f64)> = Vec::new();
    let mut min_score = f64::NEG_INFINITY;

    for i in 0..n_alleles {
        for j in i..n_alleles {
            let mut total = 0.0;
            for r in 0..n_reads {
                let s_i = read_scores[r][i];
                let s_j = read_scores[r][j];
                let max_s = s_i.max(s_j);
                let lse = max_s
                    + (10.0_f64.powf(s_i - max_s) + 10.0_f64.powf(s_j - max_s)).log10();
                total += log10_half + lse;
            }

            if best.len() < top_n || total > min_score {
                best.push((
                    allele_defs[i].0.clone(),
                    allele_defs[j].0.clone(),
                    total,
                ));
                best.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap());
                best.truncate(top_n);
                min_score = best.last().map_or(f64::NEG_INFINITY, |x| x.2);
            }
        }
    }

    best
}

/// Extract the base allele number from a star allele name.
/// "*10.002" → "*10", "*4" → "*4", "*36" → "*36"
fn base_allele(name: &str) -> String {
    if let Some(dot_pos) = name.find('.') {
        // Check if everything after the dot is digits (sub-allele number)
        let suffix = &name[dot_pos + 1..];
        if !suffix.is_empty() && suffix.chars().all(|c| c.is_ascii_digit()) {
            return name[..dot_pos].to_string();
        }
    }
    name.to_string()
}

/// Parse a haplotype string into its component alleles.
/// "*36+*10" → ["*36", "*10"], "*4x2" → ["*4"], "*1" → ["*1"]
fn parse_haplotype_components(hap: &str) -> Vec<String> {
    let mut components = Vec::new();
    for tandem_part in hap.split('+') {
        let trimmed = tandem_part.trim();
        // Strip duplication markers like "x2"
        let base = if let Some(pos) = trimmed.find('x') {
            let suffix = &trimmed[pos + 1..];
            if !suffix.is_empty() && suffix.chars().all(|c| c.is_ascii_digit()) {
                &trimmed[..pos]
            } else {
                trimmed
            }
        } else {
            trimmed
        };
        if !base.is_empty() {
            components.push(base.to_string());
        }
    }
    components
}

/// Check if an allele is a deletion (*5).
fn is_deletion_allele(allele: &str) -> bool {
    allele == "*5"
}

/// Check if voting allele `vote` is compatible with called allele `called`,
/// considering variant-set superset relationships.
/// Two alleles are compatible if:
/// 1. Their base allele names match, OR
/// 2. The voting allele's variant set is a superset of the called allele's variant set
///    (meaning voting can't distinguish them — the superset always scores ≥ the subset)
fn allele_compatible(
    vote: &str,
    called: &str,
    allele_defs: &[(String, HashSet<usize>)],
) -> bool {
    if base_allele(vote) == base_allele(called) {
        return true;
    }
    // If the called allele is not in allele_defs (e.g., hybrid like *68),
    // it has no variant definition for the D6 portion — treat as reference-like.
    // Any voting allele that's a superset of the empty set would match,
    // so just accept any vote for unknown called alleles.
    let called_vars = allele_defs.iter().find(|(n, _)| n == called).map(|(_, v)| v);
    if called_vars.is_none() {
        return true;
    }
    // Check superset: if vote's variants ⊇ called's variants, they're indistinguishable
    let vote_vars = allele_defs.iter().find(|(n, _)| n == vote).map(|(_, v)| v);
    if let (Some(vv), Some(cv)) = (vote_vars, called_vars) {
        if cv.is_subset(vv) {
            return true;
        }
    }
    false
}

/// Check if a voting diplotype (va, vb) is compatible with the called genotype.
/// Compatibility means: the voting alleles are base-allele matches for
/// some valid interpretation of the called haplotypes (accounting for
/// tandems, deletions, sub-alleles, and variant-set supersets).
fn is_compatible_diplotype(
    vote_a: &str,
    vote_b: &str,
    called_genotype: &str,
    allele_defs: &[(String, HashSet<usize>)],
) -> bool {
    let haps: Vec<&str> = called_genotype.split('/').collect();
    if haps.len() != 2 {
        return false;
    }

    // Get all allele components from each called haplotype
    let hap0_components = parse_haplotype_components(haps[0]);
    let hap1_components = parse_haplotype_components(haps[1]);

    // Handle deletions: if a haplotype is *5, its reads don't exist,
    // so voting sees only the other haplotype (as homozygous)
    let hap0_is_del = hap0_components.iter().any(|c| is_deletion_allele(c));
    let hap1_is_del = hap1_components.iter().any(|c| is_deletion_allele(c));

    // Build the set of expected alleles that voting might see
    let mut expected_alleles: Vec<String> = Vec::new();

    if hap0_is_del {
        // Only hap1 produces reads; voting sees it doubled
        for c in &hap1_components {
            if !is_deletion_allele(c) {
                expected_alleles.push(c.clone());
            }
        }
        // Voting sees homozygous for the non-deletion hap
        expected_alleles.extend(expected_alleles.clone());
    } else if hap1_is_del {
        for c in &hap0_components {
            if !is_deletion_allele(c) {
                expected_alleles.push(c.clone());
            }
        }
        expected_alleles.extend(expected_alleles.clone());
    } else {
        // Both haplotypes produce reads
        for c in &hap0_components {
            expected_alleles.push(c.clone());
        }
        for c in &hap1_components {
            expected_alleles.push(c.clone());
        }
    }

    // Check if (vote_a, vote_b) can be explained by the expected alleles.
    // A vote allele matches an expected allele if base alleles match OR
    // if the vote allele is a variant-set superset of the expected allele.
    for (i, ea) in expected_alleles.iter().enumerate() {
        if allele_compatible(vote_a, ea, allele_defs) {
            for (j, eb) in expected_alleles.iter().enumerate() {
                if j != i && allele_compatible(vote_b, eb, allele_defs) {
                    return true;
                }
            }
        }
    }

    false
}

/// Run read-level voting and compare with the pileup-based call.
pub fn validate_call(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    region_start: i64,
    region_end: i64,
    var_list: &[String],
    star_combinations: &StarCombinations,
    called_genotype: Option<&str>,
) -> VotingResult {
    // Skip voting for homozygous deletions — no CYP2D6 reads exist
    if let Some(geno) = called_genotype {
        let haps: Vec<&str> = geno.split('/').collect();
        if haps.len() == 2 {
            let h0 = parse_haplotype_components(haps[0]);
            let h1 = parse_haplotype_components(haps[1]);
            let both_del = h0.iter().all(|c| is_deletion_allele(c))
                && h1.iter().all(|c| is_deletion_allele(c));
            if both_del {
                log::info!("read_voting: skipping homozygous deletion '{}'", geno);
                return VotingResult {
                    top_diplotypes: Vec::new(),
                    n_informative_reads: 0,
                    call_rank: Some(0),
                    call_score_gap: Some(0.0),
                    flag: None,
                };
            }
        }
    }

    let sites = parse_variant_sites(var_list);
    if sites.is_empty() {
        return VotingResult {
            top_diplotypes: Vec::new(),
            n_informative_reads: 0,
            call_rank: None,
            call_score_gap: None,
            flag: None,
        };
    }

    let allele_defs = build_allele_defs(&star_combinations.dstar, &sites);

    log::info!(
        "read_voting: {} variant sites, {} allele definitions",
        sites.len(),
        allele_defs.len(),
    );

    let profiles = collect_read_profiles(reader, nchr, region_start, region_end, &sites);
    let n_informative = profiles.len();

    log::info!("read_voting: {} informative reads collected", n_informative);

    if n_informative == 0 {
        return VotingResult {
            top_diplotypes: Vec::new(),
            n_informative_reads: 0,
            call_rank: None,
            call_score_gap: None,
            flag: None,
        };
    }

    let top = find_best_diplotypes(&profiles, &allele_defs, &sites, 20);

    // Compare with the called genotype using compatibility matching.
    // Handles deletions (*5), tandems (*36+*10), sub-alleles (*10.002),
    // and duplications (*4x2).
    let (call_rank, call_score_gap) = if let Some(geno) = called_genotype {
        let rank = top.iter().position(|(ta, tb, _)| {
            is_compatible_diplotype(ta, tb, geno, &allele_defs)
        });
        let gap = rank.map(|r| {
            if r == 0 {
                0.0
            } else {
                top[0].2 - top[r].2
            }
        });
        (rank, gap)
    } else {
        (None, None)
    };

    // Flag based on score gap, not rank. Many alleles tie in score
    // (sub-alleles like *10.002, *10.009 differ by variants not covered by any read).
    // A score gap < 5.0 log10 units means the called diplotype is essentially tied
    // with the best. Only flag when there's a significant gap.
    const SCORE_GAP_THRESHOLD: f64 = 5.0;

    let flag = match call_score_gap {
        Some(gap) if gap < SCORE_GAP_THRESHOLD => {
            if gap > 0.0 {
                log::info!(
                    "read_voting: called '{}' compatible at rank {} (gap={:.1}, within threshold)",
                    called_genotype.unwrap_or("?"),
                    call_rank.map_or("?".to_string(), |r| (r + 1).to_string()),
                    gap,
                );
            }
            None // score gap is small, call is consistent
        }
        Some(gap) => {
            let best_name = if !top.is_empty() {
                format!("{}/{}", top[0].0, top[0].1)
            } else {
                "?".to_string()
            };
            log::info!(
                "read_voting: called '{}' disagrees with voting (rank={:?}, gap={:.1}); best is '{}'",
                called_genotype.unwrap_or("?"),
                call_rank,
                gap,
                best_name,
            );
            Some(format!(
                "Read_voting_disagrees_best_{}",
                best_name.replace('/', "_"),
            ))
        }
        None => {
            // Called genotype not found in top 20 at all
            let best_name = if !top.is_empty() {
                format!("{}/{}", top[0].0, top[0].1)
            } else {
                "?".to_string()
            };
            log::info!(
                "read_voting: called '{}' not found in top 20; best is '{}'",
                called_genotype.unwrap_or("?"),
                best_name,
            );
            Some(format!(
                "Read_voting_disagrees_best_{}",
                best_name.replace('/', "_"),
            ))
        }
    };

    // Log top 5
    for (i, (a, b, s)) in top.iter().take(5).enumerate() {
        log::info!("read_voting: #{} {}/{} score={:.1}", i + 1, a, b, s);
    }

    VotingResult {
        top_diplotypes: top,
        n_informative_reads: n_informative,
        call_rank,
        call_score_gap,
        flag,
    }
}

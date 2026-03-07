//! Read-level phasing constraints for diplotype validation.
//!
//! Inspired by Aldy's approach: mutations observed on the same read must belong
//! to the same allele. This module extracts per-read variant patterns from the BAM,
//! aggregates them into phasing modes, and scores candidate diplotypes for
//! consistency with the observed phasing evidence.
//!
//! Usage: after the pileup-based caller produces candidate diplotypes, run
//! `score_phasing()` to get a phasing consistency score that can break ties
//! or flag inconsistent calls.

use std::collections::{HashMap, HashSet};

use rust_htslib::bam::{self, Read as BamRead};

use crate::types::StarCombinations;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// A phasing mode: a set of (variant_name, allele) pairs observed on the same reads.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct PhasingMode {
    /// Sorted (variant_name, is_alt) pairs for variants covered by the reads.
    pub pattern: Vec<(String, bool)>,
}

/// Result of phasing analysis for a candidate diplotype.
#[derive(Debug, Clone)]
pub struct PhasingResult {
    /// Candidate diplotype (e.g., "*1/*4")
    pub diplotype: String,
    /// Number of phasing modes consistent with this diplotype
    pub consistent_modes: usize,
    /// Number of phasing modes inconsistent with this diplotype
    pub inconsistent_modes: usize,
    /// Total reads supporting consistent modes
    pub consistent_reads: usize,
    /// Total reads supporting inconsistent modes
    pub inconsistent_reads: usize,
    /// Phasing score: log-ratio of consistent vs inconsistent evidence
    pub phasing_score: f64,
}

/// Full output from the read-phasing module.
#[derive(Debug, Clone)]
pub struct ReadPhasingOutput {
    /// Number of distinct phasing modes (multi-variant read patterns)
    pub n_modes: usize,
    /// Total reads contributing to phasing modes
    pub n_phasing_reads: usize,
    /// Per-candidate phasing results
    pub candidate_scores: Vec<PhasingResult>,
    /// Best candidate by phasing score
    pub best_candidate: Option<String>,
    /// Flag if phasing contradicts the pileup call
    pub flag: Option<String>,
}

// ---------------------------------------------------------------------------
// Variant site for phasing
// ---------------------------------------------------------------------------

/// A variant position with REF and ALT bases for phasing extraction.
#[derive(Debug, Clone)]
struct PhaseSite {
    /// 0-based genomic position
    pos_0based: i64,
    /// Variant name (e.g., "g.42128945C>T")
    name: String,
    /// Expected REF base (upper-case)
    ref_base: u8,
    /// Expected ALT base (upper-case)
    alt_base: u8,
    /// True if deletion variant
    is_deletion: bool,
}

// ---------------------------------------------------------------------------
// Extract phase sites from target variant file
// ---------------------------------------------------------------------------

fn build_phase_sites(var_list: &[String]) -> Vec<PhaseSite> {
    let mut sites = Vec::new();
    for name in var_list {
        if let Some(site) = parse_variant_for_phasing(name) {
            sites.push(site);
        }
    }
    sites.sort_by_key(|s| s.pos_0based);
    sites
}

fn parse_variant_for_phasing(name: &str) -> Option<PhaseSite> {
    let s = name.strip_prefix("g.")?;
    let digit_end = s.find(|c: char| !c.is_ascii_digit())?;
    let pos: i64 = s[..digit_end].parse().ok()?;
    let rest = &s[digit_end..];

    // Single-base SNP: "C>T"
    if rest.len() == 3 && rest.as_bytes()[1] == b'>' {
        let ref_base = rest.as_bytes()[0].to_ascii_uppercase();
        let alt_base = rest.as_bytes()[2].to_ascii_uppercase();
        if ref_base.is_ascii_alphabetic() && alt_base.is_ascii_alphabetic() {
            return Some(PhaseSite {
                pos_0based: pos - 1,
                name: name.to_string(),
                ref_base,
                alt_base,
                is_deletion: false,
            });
        }
    }

    // Single-base deletion: "delT"
    if rest.starts_with("del") && rest.len() == 4 {
        let ref_base = rest.as_bytes()[3].to_ascii_uppercase();
        return Some(PhaseSite {
            pos_0based: pos - 1,
            name: name.to_string(),
            ref_base,
            alt_base: 0,
            is_deletion: true,
        });
    }

    None
}

// ---------------------------------------------------------------------------
// Extract per-read variant observations from BAM
// ---------------------------------------------------------------------------

/// Extract per-fragment phasing patterns from the BAM.
/// Returns aggregated modes: pattern → read count.
fn extract_phasing_modes(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    sites: &[PhaseSite],
) -> HashMap<Vec<(String, bool)>, usize> {
    if sites.is_empty() {
        return HashMap::new();
    }

    let tid = match reader.header().tid(nchr.as_bytes()) {
        Some(t) => t,
        None => return HashMap::new(),
    };

    // Build position lookup: 0-based pos → site index
    let mut pos_to_site: HashMap<i64, usize> = HashMap::new();
    for (i, site) in sites.iter().enumerate() {
        pos_to_site.insert(site.pos_0based, i);
    }

    let region_start = sites.first().unwrap().pos_0based;
    let region_end = sites.last().unwrap().pos_0based + 1;

    // Collect per-fragment observations
    let mut fragment_obs: HashMap<String, Vec<(usize, bool)>> = HashMap::new();

    if reader
        .fetch(bam::FetchDefinition::Region(
            tid as i32,
            region_start,
            region_end,
        ))
        .is_err()
    {
        return HashMap::new();
    }

    for result in reader.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };

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

        let ref_start = record.pos();
        let seq = record.seq().as_bytes();
        let cigar = record.cigar();
        let ref_end = cigar.end_pos();

        let fragment_name = String::from_utf8_lossy(record.qname()).to_string();

        for (&pos, &site_idx) in &pos_to_site {
            if pos < ref_start || pos >= ref_end {
                continue;
            }

            let site = &sites[site_idx];

            if site.is_deletion {
                let read_pos = cigar.read_pos(pos as u32, false, false).ok().flatten();
                let is_alt = read_pos.is_none(); // deletion = no base at position
                fragment_obs
                    .entry(fragment_name.clone())
                    .or_default()
                    .push((site_idx, is_alt));
            } else {
                let read_pos = match cigar.read_pos(pos as u32, false, false) {
                    Ok(Some(rp)) => rp as usize,
                    _ => continue,
                };
                if read_pos >= seq.len() {
                    continue;
                }
                let observed = seq[read_pos].to_ascii_uppercase();
                let is_alt = observed == site.alt_base;
                let is_ref = observed == site.ref_base;
                if !is_alt && !is_ref {
                    continue; // sequencing error or third allele
                }
                fragment_obs
                    .entry(fragment_name.clone())
                    .or_default()
                    .push((site_idx, is_alt));
            }
        }
    }

    // Aggregate into modes (patterns with 2+ sites)
    let mut modes: HashMap<Vec<(String, bool)>, usize> = HashMap::new();
    for (_fragment, obs) in &fragment_obs {
        if obs.len() < 2 {
            continue;
        }

        // Deduplicate (same site may be seen from mate)
        let mut deduped: HashMap<usize, bool> = HashMap::new();
        for &(site_idx, is_alt) in obs {
            deduped.entry(site_idx).or_insert(is_alt);
        }
        if deduped.len() < 2 {
            continue;
        }

        let mut pattern: Vec<(String, bool)> = deduped
            .iter()
            .map(|(&idx, &is_alt)| (sites[idx].name.clone(), is_alt))
            .collect();
        pattern.sort_by(|a, b| a.0.cmp(&b.0));

        *modes.entry(pattern).or_insert(0) += 1;
    }

    modes
}

// ---------------------------------------------------------------------------
// Score candidate diplotypes against phasing modes
// ---------------------------------------------------------------------------

/// Get the set of defining variants for a star allele (or tandem like *36+*10).
fn get_allele_variants(star: &str, dstar: &HashMap<String, String>) -> HashSet<String> {
    let mut all_variants = HashSet::new();

    // Handle tandem alleles: "*36+*10" → look up "*36" and "*10" separately
    // Handle duplication markers: "*4x2" → "*4", "*36x2+*10x2" → "*36" + "*10"
    let components: Vec<&str> = star.split('+').collect();
    for component in &components {
        // Strip xN duplication suffix (e.g., "*4x2" → "*4")
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
        let lookup_name = if base.starts_with('*') {
            base.to_string()
        } else {
            format!("*{}", base)
        };
        if let Some(var_key) = dstar.get(&lookup_name) {
            if var_key != "NA" {
                for v in var_key.split('_') {
                    all_variants.insert(v.to_string());
                }
            }
        }
    }

    all_variants
}

/// Check if a phasing mode is consistent with a diplotype.
///
/// A mode is consistent if ALL its observations can be explained by
/// a single allele in the diplotype. That is, every ALT observation in
/// the mode must come from the same haplotype, and every REF observation
/// must come from positions NOT in that haplotype.
fn mode_consistent_with_diplotype(
    mode: &[(String, bool)],
    hap1_variants: &HashSet<String>,
    hap2_variants: &HashSet<String>,
) -> bool {
    // Check if mode matches haplotype 1
    let matches_h1 = mode.iter().all(|(var, is_alt)| {
        let on_h1 = hap1_variants.contains(var);
        if *is_alt { on_h1 } else { !on_h1 }
    });

    // Check if mode matches haplotype 2
    let matches_h2 = mode.iter().all(|(var, is_alt)| {
        let on_h2 = hap2_variants.contains(var);
        if *is_alt { on_h2 } else { !on_h2 }
    });

    // A mode where all observations are REF is always consistent (uninformative)
    let all_ref = mode.iter().all(|(_, is_alt)| !is_alt);
    if all_ref {
        return true;
    }

    // Check if mode matches reference (no ALT matches any haplotype)
    let matches_ref = mode.iter().all(|(var, is_alt)| {
        if *is_alt {
            // ALT must not be on either haplotype for a "reference" match
            false
        } else {
            true
        }
    });

    matches_h1 || matches_h2 || matches_ref
}

/// Score a single candidate diplotype against all phasing modes.
fn score_candidate(
    diplotype: &str,
    modes: &HashMap<Vec<(String, bool)>, usize>,
    dstar: &HashMap<String, String>,
) -> PhasingResult {
    let parts: Vec<&str> = diplotype.split('/').collect();
    if parts.len() != 2 {
        return PhasingResult {
            diplotype: diplotype.to_string(),
            consistent_modes: 0,
            inconsistent_modes: 0,
            consistent_reads: 0,
            inconsistent_reads: 0,
            phasing_score: 0.0,
        };
    }

    let hap1_variants = get_allele_variants(parts[0], dstar);
    let hap2_variants = get_allele_variants(parts[1], dstar);

    let mut consistent_modes = 0usize;
    let mut inconsistent_modes = 0usize;
    let mut consistent_reads = 0usize;
    let mut inconsistent_reads = 0usize;

    for (pattern, &count) in modes {
        // Only evaluate modes that contain at least one ALT observation
        // (all-REF modes are uninformative)
        let has_alt = pattern.iter().any(|(_, is_alt)| *is_alt);
        if !has_alt {
            continue;
        }

        if mode_consistent_with_diplotype(pattern, &hap1_variants, &hap2_variants) {
            consistent_modes += 1;
            consistent_reads += count;
        } else {
            inconsistent_modes += 1;
            inconsistent_reads += count;
        }
    }

    let total_reads = (consistent_reads + inconsistent_reads) as f64;
    let phasing_score = if total_reads > 0.0 {
        let frac_consistent = consistent_reads as f64 / total_reads;
        // Log-odds ratio: positive = more consistent, negative = more inconsistent
        if frac_consistent > 0.0 && frac_consistent < 1.0 {
            (frac_consistent / (1.0 - frac_consistent)).ln()
        } else if frac_consistent >= 1.0 {
            10.0
        } else {
            -10.0
        }
    } else {
        0.0
    };

    PhasingResult {
        diplotype: diplotype.to_string(),
        consistent_modes,
        inconsistent_modes,
        consistent_reads,
        inconsistent_reads,
        phasing_score,
    }
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Analyze read-level phasing for the given pileup call and its alternatives.
///
/// `pileup_call`: the main call from match_star (e.g., "*1/*4" or "*1/*4;*2/*21")
/// `raw_candidates`: raw star allele combinations from match_star
pub fn score_phasing(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    var_list: &[String],
    star_combinations: &StarCombinations,
    pileup_call: Option<&str>,
) -> ReadPhasingOutput {
    let sites = build_phase_sites(var_list);
    if sites.is_empty() {
        return ReadPhasingOutput {
            n_modes: 0,
            n_phasing_reads: 0,
            candidate_scores: Vec::new(),
            best_candidate: None,
            flag: None,
        };
    }

    let modes = extract_phasing_modes(reader, nchr, &sites);
    let n_modes = modes.values().filter(|&&c| c > 0).count();
    let n_phasing_reads: usize = modes.values().sum();

    log::info!(
        "read_phasing: {} modes from {} reads ({} variant sites)",
        n_modes,
        n_phasing_reads,
        sites.len(),
    );

    if n_modes == 0 {
        return ReadPhasingOutput {
            n_modes: 0,
            n_phasing_reads: 0,
            candidate_scores: Vec::new(),
            best_candidate: None,
            flag: None,
        };
    }

    // Build candidate list from pileup call
    let candidates: Vec<String> = match pileup_call {
        Some(call) => call.split(';').map(|s| s.trim().to_string()).collect(),
        None => Vec::new(),
    };

    let dstar = &star_combinations.dstar;

    let mut candidate_scores: Vec<PhasingResult> = candidates
        .iter()
        .map(|c| score_candidate(c, &modes, dstar))
        .collect();

    // Sort by phasing score (descending)
    candidate_scores.sort_by(|a, b| {
        b.phasing_score
            .partial_cmp(&a.phasing_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Log results
    for r in &candidate_scores {
        log::info!(
            "read_phasing: {} — consistent={}/{} inconsistent={}/{} score={:.2}",
            r.diplotype,
            r.consistent_modes,
            r.consistent_reads,
            r.inconsistent_modes,
            r.inconsistent_reads,
            r.phasing_score,
        );
    }

    let best_candidate = candidate_scores.first().map(|r| r.diplotype.clone());

    // Flag if the best phasing candidate differs from the first pileup call
    let flag = if candidate_scores.len() >= 2 {
        let best = &candidate_scores[0];
        let second = &candidate_scores[1];
        if best.phasing_score > second.phasing_score + 1.0
            && best.consistent_reads >= 3
            && candidates.len() > 1
        {
            // Check if best differs from the first (default) candidate
            if candidates.first().map(|s| s.as_str()) != Some(best.diplotype.as_str()) {
                Some(format!(
                    "Phasing_prefers_{}",
                    best.diplotype.replace('/', "_")
                ))
            } else {
                None
            }
        } else {
            None
        }
    } else {
        None
    };

    ReadPhasingOutput {
        n_modes,
        n_phasing_reads,
        candidate_scores,
        best_candidate,
        flag,
    }
}

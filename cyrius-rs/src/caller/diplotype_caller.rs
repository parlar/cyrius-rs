//! Likelihood scoring of candidate diplotypes, for the resolver.
//!
//! `score_candidates` scores an explicit set of candidate diplotypes (the
//! `;`-ambiguous calls the matcher left) against per-read evidence. Reads are
//! collected from both the CYP2D6 and CYP2D7 regions; CYP2D7 is modelled as an
//! additional pseudo-allele so that mismapped reads are properly accounted for.
//!
//! The model for a single read under a diplotype given as D6 alleles with copy
//! counts summing to CN, with one copy of D7 always present:
//!
//!   P(read | dip) = w_d7 * P(read | D7) + sum_k w_k * P(read | allele_k)
//!
//! where w_d7 = 1/(CN+1) and w_k = copies_k/(CN+1) (for CN=2 with two
//! distinct alleles: w_d7 = w_i = w_j = 1/3).

use std::collections::{HashMap, HashSet};

use rust_htslib::bam::{self, Read as BamRead};

use crate::depth_calling::utilities::fetch_region;
use crate::types::StarCombinations;

// ---------------------------------------------------------------------------
// Precomputed quality-score lookup table
// ---------------------------------------------------------------------------

/// Precomputed log-likelihood values for each Phred quality score (0..=93).
/// Avoids repeated calls to 10.0_f64.powf() in the hot scoring loop.
struct QualLookup {
    /// log10(1 - p_err) for match
    log10_match: [f64; 94],
    /// log10(p_err / 3) for mismatch
    log10_mismatch: [f64; 94],
}

impl QualLookup {
    fn new() -> Self {
        let mut log10_match = [0.0f64; 94];
        let mut log10_mismatch = [0.0f64; 94];
        for q in 0..94u8 {
            let p_err = 10.0_f64.powf(-(q as f64) / 10.0).clamp(1e-10, 0.5);
            log10_match[q as usize] = (1.0 - p_err).log10();
            log10_mismatch[q as usize] = (p_err / 3.0).log10();
        }
        QualLookup {
            log10_match,
            log10_mismatch,
        }
    }

    #[inline(always)]
    fn score(&self, is_match: bool, qual: u8) -> f64 {
        let q = (qual as usize).min(93);
        if is_match {
            self.log10_match[q]
        } else {
            self.log10_mismatch[q]
        }
    }
}

// ---------------------------------------------------------------------------
// Public result types
// ---------------------------------------------------------------------------

/// One scored candidate diplotype (for the resolver).
#[derive(Debug, Clone)]
pub struct CandidateScore {
    pub diplotype: String,
    pub log10_likelihood: f64,
    /// Posterior normalised over the candidate set only (0 for unscorable).
    pub posterior: f64,
    /// Reads whose likelihood under this candidate beats every other
    /// candidate by more than DECISIVE_LOG10_LR.
    pub decisive_reads: usize,
}

// ---------------------------------------------------------------------------
// Variant site representation (shared between D6 target variants and D6/D7
// diagnostic SNPs)
// ---------------------------------------------------------------------------

/// A genomic position where we observe evidence.
/// Covers both target-variant sites (star-allele defining) and D6/D7
/// diagnostic SNPs.
#[derive(Debug, Clone)]
struct UnifiedSite {
    /// 1-based genomic position in the D6 coordinate space.
    d6_pos: i64,
    /// 1-based genomic position in D7 (for fetching reads mapped there).
    d7_pos: Option<i64>,
    /// The D6-reference base (upper-case, single byte for SNPs).
    d6_ref_base: u8,
    /// The D7 base at this position (if known from SNP file).
    d7_base: Option<u8>,
    /// For star-allele-defining SNPs: the ALT base on D6.
    /// None for pure D6/D7 diagnostic positions that aren't star-allele variants.
    d6_alt_base: Option<u8>,
    /// Human-readable name (e.g. "g.42128945C>T" or "diag_42127634").
    name: String,
    /// True if this is a single-base deletion variant.
    is_deletion: bool,
}

/// Per-read observations at unified sites.
#[derive(Debug)]
struct ReadObs {
    #[allow(dead_code)]
    name: String,
    /// (site_index, observed_base, base_quality).
    /// observed_base == 0 means the base was deleted.
    observations: Vec<(usize, u8, u8)>,
}

/// An allele definition: which unified-site indices carry the ALT allele.
#[derive(Debug, Clone)]
struct AlleleDef {
    name: String,
    /// Indices into the unified sites vector where this allele carries ALT.
    alt_sites: HashSet<usize>,
}

// ---------------------------------------------------------------------------
// Site construction
// ---------------------------------------------------------------------------

/// Build the unified site list from the target-variant names and the SNP file.
///
/// * `var_list`  — target variant names (e.g. "g.42128945C>T").
/// * `snp_file`  — raw content of CYP2D6_SNP_38.txt.
///
/// Returns sites covering:
/// 1. Every parseable single-base SNP from var_list (star-allele defining).
/// 2. Every single-base D6/D7 diagnostic position from the SNP file that is
///    NOT already covered by (1).
fn build_unified_sites(var_list: &[String], snp_file: &str) -> Vec<UnifiedSite> {
    let mut sites: Vec<UnifiedSite> = Vec::new();
    let mut seen_d6_pos: HashSet<i64> = HashSet::new();

    // 1. Target variants (star-allele defining)
    for name in var_list {
        if let Some(site) = parse_target_variant(name) {
            if seen_d6_pos.insert(site.d6_pos) {
                sites.push(site);
            }
        }
    }

    // 2. D6/D7 diagnostic SNPs — only single-base ones
    for line in snp_file.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 6 {
            continue;
        }
        let d6_pos: i64 = match cols[1].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let d6_base_str = cols[2];
        let d7_pos: i64 = match cols[3].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let d7_base_str = cols[4];

        // Only single-base SNPs
        if d6_base_str.len() != 1 || d7_base_str.len() != 1 {
            continue;
        }
        let d6_base = d6_base_str.as_bytes()[0].to_ascii_uppercase();
        let d7_base = d7_base_str.as_bytes()[0].to_ascii_uppercase();
        if d6_base == d7_base {
            continue;
        }

        if seen_d6_pos.insert(d6_pos) {
            sites.push(UnifiedSite {
                d6_pos,
                d7_pos: Some(d7_pos),
                d6_ref_base: d6_base,
                d7_base: Some(d7_base),
                d6_alt_base: None,
                name: format!("diag_{}", d6_pos),
                is_deletion: false,
            });
        } else {
            // Already have this position from target variants — enrich with D7 info
            if let Some(existing) = sites.iter_mut().find(|s| s.d6_pos == d6_pos) {
                if existing.d7_pos.is_none() {
                    existing.d7_pos = Some(d7_pos);
                }
                if existing.d7_base.is_none() {
                    existing.d7_base = Some(d7_base);
                }
            }
        }
    }

    sites
}

/// Parse a single-base SNP from a variant name like "g.42128945C>T".
fn parse_target_variant(name: &str) -> Option<UnifiedSite> {
    let s = name.strip_prefix("g.")?;
    let digit_end = s.find(|c: char| !c.is_ascii_digit())?;
    let pos: i64 = s[..digit_end].parse().ok()?;
    let rest = &s[digit_end..];

    // Single-base SNP: "C>T"
    if rest.len() == 3 && rest.as_bytes()[1] == b'>' {
        let ref_base = rest.as_bytes()[0].to_ascii_uppercase();
        let alt_base = rest.as_bytes()[2].to_ascii_uppercase();
        if ref_base.is_ascii_alphabetic() && alt_base.is_ascii_alphabetic() {
            return Some(UnifiedSite {
                d6_pos: pos,
                d7_pos: None,
                d6_ref_base: ref_base,
                d7_base: None,
                d6_alt_base: Some(alt_base),
                name: name.to_string(),
                is_deletion: false,
            });
        }
    }

    // Single-base deletion: "delT"
    if let Some(del_seq) = rest.strip_prefix("del") {
        if del_seq.len() == 1 && del_seq.as_bytes()[0].is_ascii_alphabetic() {
            let ref_base = del_seq.as_bytes()[0].to_ascii_uppercase();
            return Some(UnifiedSite {
                d6_pos: pos,
                d7_pos: None,
                d6_ref_base: ref_base,
                d7_base: None,
                d6_alt_base: Some(0), // deletion
                name: name.to_string(),
                is_deletion: true,
            });
        }
    }

    None
}

// ---------------------------------------------------------------------------
// D6 allele definitions (CYP2D7 is scored directly by `score_read_d7`)
// ---------------------------------------------------------------------------

/// Build D6 star-allele definitions from the dstar map.
fn build_d6_allele_defs(
    dstar: &HashMap<String, String>,
    sites: &[UnifiedSite],
) -> Vec<AlleleDef> {
    let name_to_idx: HashMap<&str, usize> = sites
        .iter()
        .enumerate()
        .filter(|(_, s)| s.d6_alt_base.is_some()) // only star-allele sites
        .map(|(i, s)| (s.name.as_str(), i))
        .collect();

    let mut defs: Vec<AlleleDef> = dstar
        .iter()
        .map(|(allele_name, variant_key)| {
            let mut alt_sites = HashSet::new();
            if variant_key != "NA" {
                for var_name in variant_key.split('_') {
                    if let Some(&idx) = name_to_idx.get(var_name) {
                        alt_sites.insert(idx);
                    }
                }
            }
            AlleleDef {
                name: allele_name.clone(),
                alt_sites,
            }
        })
        .collect();

    defs.sort_by(|a, b| a.name.cmp(&b.name));
    defs
}

// ---------------------------------------------------------------------------
// Read collection from both D6 and D7 regions
// ---------------------------------------------------------------------------

/// Collect read observations from both D6 and D7 genomic regions.
fn collect_reads_dual_region(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    sites: &[UnifiedSite],
    d6_start: i64,
    d6_end: i64,
    d7_start: i64,
    d7_end: i64,
) -> Vec<ReadObs> {
    let tid = match reader.header().tid(nchr.as_bytes()) {
        Some(t) => t,
        None => return Vec::new(),
    };

    // Build lookup: genomic_pos (0-based) → (site_index, is_d7_coord)
    let mut pos_lookup: HashMap<i64, Vec<(usize, bool)>> = HashMap::new();
    for (i, site) in sites.iter().enumerate() {
        // D6 coordinate (always present)
        pos_lookup
            .entry(site.d6_pos - 1) // convert to 0-based
            .or_default()
            .push((i, false));
        // D7 coordinate (if known)
        if let Some(d7p) = site.d7_pos {
            pos_lookup
                .entry(d7p - 1)
                .or_default()
                .push((i, true));
        }
    }

    let mut profiles: HashMap<String, Vec<(usize, u8, u8)>> = HashMap::new();

    // Fetch from both regions
    let regions = [(d6_start, d6_end), (d7_start, d7_end)];
    for &(start, end) in &regions {
        // A failed fetch means part of the evidence was never looked at;
        // scoring the rest could narrow a call it should not (E1). No reads
        // -> not scorable -> the call stays ambiguous.
        if let Err(e) = fetch_region(reader, tid, start, end, &format!("{nchr}:{start}-{end}")) {
            log::warn!("resolver: {e}; leaving the call unresolved");
            return Vec::new();
        }

        for result in reader.records() {
            let mut record = match result {
                Ok(r) => r,
                Err(e) => {
                    log::warn!("resolver: BAM read failed ({e}); leaving the call unresolved");
                    return Vec::new();
                }
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

            let ref_start = record.pos(); // 0-based
            record.cache_cigar();
            let cigar = match record.cigar_cached() {
                Some(c) => c,
                None => continue,
            };
            let ref_end = cigar.end_pos(); // 0-based exclusive

            // Only allocate read name once per read, and only if it overlaps any site
            let seq = record.seq().as_bytes();
            let qual = record.qual().to_vec();
            let mut read_name: Option<String> = None;

            // Iterate over variant sites (not over every ref position)
            // and check if each site falls within the read's span.
            for (&ref_pos_0, entries) in &pos_lookup {
                if ref_pos_0 < ref_start || ref_pos_0 >= ref_end {
                    continue;
                }

                for &(site_idx, _is_d7_coord) in entries {
                    let site = &sites[site_idx];

                    let read_pos_result =
                        cigar.read_pos(ref_pos_0 as u32, false, false).ok().flatten();

                    if site.is_deletion {
                        let (base, q) = match read_pos_result {
                            None => (0u8, 30u8),
                            Some(rp) => {
                                let rp = rp as usize;
                                if rp < seq.len() {
                                    (seq[rp], if rp < qual.len() { qual[rp] } else { 30 })
                                } else {
                                    continue;
                                }
                            }
                        };
                        let name = read_name.get_or_insert_with(|| {
                            String::from_utf8_lossy(record.qname()).to_string()
                        });
                        profiles
                            .entry(name.clone())
                            .or_default()
                            .push((site_idx, base, q));
                    } else if let Some(read_pos) = read_pos_result {
                        let rp = read_pos as usize;
                        if rp < seq.len() {
                            let base = seq[rp];
                            let q = if rp < qual.len() { qual[rp] } else { 30 };
                            let name = read_name.get_or_insert_with(|| {
                                String::from_utf8_lossy(record.qname()).to_string()
                            });
                            profiles
                                .entry(name.clone())
                                .or_default()
                                .push((site_idx, base, q));
                        }
                    }
                }
            }
        }
    }

    profiles
        .into_iter()
        .map(|(name, observations)| ReadObs { name, observations })
        .filter(|p| !p.observations.is_empty())
        .collect()
}

// ---------------------------------------------------------------------------
// Scoring
// ---------------------------------------------------------------------------

/// Log10-likelihood of a read under a D6 allele hypothesis.
/// For star-allele-defining sites: checks if observed matches expected (ref or alt).
/// For diagnostic sites: D6 alleles always carry the D6-ref base.
fn score_read_d6(
    obs: &ReadObs,
    allele: &AlleleDef,
    sites: &[UnifiedSite],
    lut: &QualLookup,
) -> f64 {
    let mut score = 0.0;
    for &(site_idx, observed, qual) in &obs.observations {
        let site = &sites[site_idx];
        let carries_alt = allele.alt_sites.contains(&site_idx);

        let expected = if site.is_deletion {
            if carries_alt { 0u8 } else { site.d6_ref_base }
        } else if let Some(alt) = site.d6_alt_base {
            if carries_alt { alt } else { site.d6_ref_base }
        } else {
            site.d6_ref_base
        };

        score += base_log10_lik(observed, expected, qual, site.is_deletion, lut);
    }
    score
}

/// Log10-likelihood of a read under the CYP2D7 hypothesis.
/// At diagnostic sites: expects D7 base.
/// At star-allele sites without D7 info: expects D6 ref (D7 is reference-like
/// for D6-specific variants).
fn score_read_d7(obs: &ReadObs, sites: &[UnifiedSite], lut: &QualLookup) -> f64 {
    let mut score = 0.0;
    for &(site_idx, observed, qual) in &obs.observations {
        let site = &sites[site_idx];

        let expected = if let Some(d7b) = site.d7_base {
            d7b
        } else {
            site.d6_ref_base
        };

        score += base_log10_lik(observed, expected, qual, site.is_deletion, lut);
    }
    score
}

/// Log10-likelihood of observing `observed` when `expected` is the true base.
/// Uses precomputed lookup table for speed.
#[inline(always)]
fn base_log10_lik(observed: u8, expected: u8, qual: u8, is_deletion: bool, lut: &QualLookup) -> f64 {
    let is_match = if is_deletion {
        (observed == 0) == (expected == 0)
    } else {
        observed.to_ascii_uppercase() == expected
    };
    lut.score(is_match, qual)
}

// ---------------------------------------------------------------------------
// Mismapped-read filter
// ---------------------------------------------------------------------------

/// Filter out reads that score poorly under every hypothesis (D6 or D7).
/// These are likely mismapped from elsewhere in the genome.
/// Returns (kept_reads, n_filtered).
fn filter_mismapped(
    reads: Vec<ReadObs>,
    d6_alleles: &[AlleleDef],
    sites: &[UnifiedSite],
    min_avg_log10_per_site: f64,
    lut: &QualLookup,
) -> (Vec<ReadObs>, usize) {
    let mut kept = Vec::with_capacity(reads.len());
    let mut n_filtered = 0usize;

    for obs in reads {
        let n_sites = obs.observations.len() as f64;
        if n_sites == 0.0 {
            continue;
        }

        let threshold = min_avg_log10_per_site * n_sites;

        // Check D7 first (single allele, fast)
        let d7_score = score_read_d7(&obs, sites, lut);
        if d7_score >= threshold {
            kept.push(obs);
            continue;
        }

        // Check D6 alleles with early exit once any passes threshold
        let mut passed = false;
        for a in d6_alleles {
            let score = score_read_d6(&obs, a, sites, lut);
            if score >= threshold {
                passed = true;
                break;
            }
        }

        if passed {
            kept.push(obs);
        } else {
            n_filtered += 1;
        }
    }

    (kept, n_filtered)
}

// ---------------------------------------------------------------------------
// Allele-label helpers
// ---------------------------------------------------------------------------

/// Split "*36x2" -> ("*36", 2); "*36" -> ("*36", 1). Only strips a trailing
/// 'x' followed by all ASCII digits (so "*4.013" and "*1x" are left intact).
fn split_copy_suffix(s: &str) -> (String, usize) {
    if let Some(pos) = s.rfind('x') {
        let num = &s[pos + 1..];
        if !num.is_empty() && num.chars().all(|c| c.is_ascii_digit()) {
            if let Ok(n) = num.parse::<usize>() {
                return (s[..pos].to_string(), n.max(1));
            }
        }
    }
    (s.to_string(), 1)
}

// ---------------------------------------------------------------------------
// Posterior normalization helper
// ---------------------------------------------------------------------------

/// Convert a slice of log10-likelihoods into a `Vec<f64>` of posteriors,
/// normalised using the numerically stable max-shift formula:
///   posterior_i = 10^(ll_i - (max_ll + log10(Σ 10^(ll_j - max_ll))))
///
/// Returns `None` if `max_ll` is not finite (e.g. all entries are `-inf` or
/// `NaN`); the caller is responsible for deciding what to do in that case.
fn posteriors_from_log10(lls: &[f64]) -> Option<Vec<f64>> {
    let max_ll = lls.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !max_ll.is_finite() {
        return None;
    }
    let sum: f64 = lls.iter().map(|&ll| 10.0_f64.powf(ll - max_ll)).sum();
    let log10_sum = max_ll + sum.log10();
    Some(lls.iter().map(|&ll| 10.0_f64.powf(ll - log10_sum)).collect())
}

// ---------------------------------------------------------------------------
// Explicit-multiset scoring (for the diplotype resolver)
// ---------------------------------------------------------------------------

/// Per-read log10-likelihoods under one explicit diplotype given as
/// (allele_index, copies) pairs summing to total_cn. General k-allele mixture:
/// one D7 copy + the specified D6 alleles, each weighted by copies/(cn+1).
/// Their sum is the diplotype's log10-likelihood. A read with no finite
/// score contributes 0.
fn read_log10_likelihoods(
    reads: &[ReadObs],
    d6_alleles: &[AlleleDef],
    sites: &[UnifiedSite],
    total_cn: u32,
    assignment: &[(usize, usize)],
    lut: &QualLookup,
) -> Vec<f64> {
    let ln10 = std::f64::consts::LN_10;
    let denom = (total_cn as usize + 1) as f64;
    let ln_w_d7 = (1.0 / denom).ln();
    let mut out = Vec::with_capacity(reads.len());
    for obs in reads {
        let s_d7 = ln_w_d7 + score_read_d7(obs, sites, lut) * ln10;
        // Per-read Vec allocation is fine here: the resolver scores only a tiny
        // explicit candidate set, not an O(n^2) enumeration of allele pairs.
        let mut terms: Vec<f64> = Vec::with_capacity(assignment.len() + 1);
        terms.push(s_d7);
        for &(idx, copies) in assignment {
            if copies == 0 {
                // map_candidates_to_raw never yields a zero count; purely defensive.
                continue;
            }
            let ln_w = ((copies as f64) / denom).ln();
            terms.push(ln_w + score_read_d6(obs, &d6_alleles[idx], sites, lut) * ln10);
        }
        let max_s = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        if !max_s.is_finite() {
            out.push(0.0);
            continue;
        }
        let lse = max_s + terms.iter().map(|t| (t - max_s).exp()).sum::<f64>().ln();
        out.push(lse / ln10);
    }
    out
}

/// Per-read log10 likelihood ratio a read needs over every other candidate
/// to count as clear evidence for one candidate (10:1).
const DECISIVE_LOG10_LR: f64 = 1.0;

/// For each candidate, the number of reads that favour it over every other
/// candidate by more than DECISIVE_LOG10_LR. `per_candidate[c][r]` is read
/// r's log10-likelihood under candidate c.
fn count_decisive_reads(per_candidate: &[Vec<f64>]) -> Vec<usize> {
    let mut counts = vec![0usize; per_candidate.len()];
    let n_reads = per_candidate.first().map_or(0, |v| v.len());
    for r in 0..n_reads {
        let mut best = 0usize;
        for c in 1..per_candidate.len() {
            if per_candidate[c][r] > per_candidate[best][r] {
                best = c;
            }
        }
        let runner_up = (0..per_candidate.len())
            .filter(|&c| c != best)
            .map(|c| per_candidate[c][r])
            .fold(f64::NEG_INFINITY, f64::max);
        if per_candidate[best][r] - runner_up > DECISIVE_LOG10_LR {
            counts[best] += 1;
        }
    }
    counts
}

/// Map a parsed multiset to (allele_index, copies) pairs, or None if ANY
/// allele is unknown to the scorer or the copy total != cn. Pure + testable.
fn build_assignment(
    multiset: &[(String, usize)],
    name_to_idx: &std::collections::HashMap<&str, usize>,
    cn: usize,
) -> Option<Vec<(usize, usize)>> {
    let mut assignment: Vec<(usize, usize)> = Vec::new();
    let mut total_copies = 0usize;
    for (name, copies) in multiset {
        let idx = *name_to_idx.get(name.as_str())?;
        assignment.push((idx, *copies));
        total_copies += copies;
    }
    if total_copies != cn {
        return None;
    }
    Some(assignment)
}

/// Allele name -> copy count, sorted by name.
type Multiset = Vec<(String, usize)>;

/// Main-allele name: strip any sub-allele suffix ("*4.013" -> "*4").
fn main_allele(name: &str) -> &str {
    name.split('.').next().unwrap_or(name)
}

/// Main-allele multiset of a candidate label. Accepts "/"-diplotypes with
/// "+"-tandems and "xN" copies, and "_"-joined allele lists.
fn label_main_multiset(label: &str) -> Multiset {
    use std::collections::BTreeMap;
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for comp in label.split(['/', '+', '_']) {
        let comp = comp.trim();
        if comp.is_empty() {
            continue;
        }
        let (name, mult) = split_copy_suffix(comp);
        *counts.entry(main_allele(&name).to_string()).or_insert(0) += mult;
    }
    counts.into_iter().collect()
}

/// Sub-allele multiset of a raw "_"-joined matcher candidate, plus its
/// main-allele form.
fn raw_multisets(raw: &str) -> (Multiset, Multiset) {
    use std::collections::BTreeMap;
    let mut sub: BTreeMap<String, usize> = BTreeMap::new();
    let mut main: BTreeMap<String, usize> = BTreeMap::new();
    for name in raw.split('_').filter(|n| !n.is_empty()) {
        *sub.entry(name.to_string()).or_insert(0) += 1;
        *main.entry(main_allele(name).to_string()).or_insert(0) += 1;
    }
    (sub.into_iter().collect(), main.into_iter().collect())
}

/// For each candidate label, the raw sub-allele multisets (from the matcher's
/// raw candidates) whose main alleles equal the label's, plus the D6 copy
/// count they share. None (fail closed) if any label has no matching raw
/// candidate or the matched raw candidates disagree on copy count.
fn map_candidates_to_raw(
    candidates: &[String],
    raw_candidates: &[String],
) -> Option<(Vec<Vec<Multiset>>, usize)> {
    let raws: Vec<_> = raw_candidates.iter().map(|r| raw_multisets(r)).collect();
    let mut cn: Option<usize> = None;
    let mut out = Vec::with_capacity(candidates.len());
    for cand in candidates {
        let want = label_main_multiset(cand);
        let matched: Vec<Multiset> = raws
            .iter()
            .filter(|(_, main)| *main == want)
            .map(|(sub, _)| sub.clone())
            .collect();
        if matched.is_empty() {
            return None;
        }
        for m in &matched {
            let copies: usize = m.iter().map(|(_, c)| c).sum();
            match cn {
                None => cn = Some(copies),
                Some(c) if c != copies => return None,
                _ => {}
            }
        }
        out.push(matched);
    }
    Some((out, cn?))
}

/// Score a set of candidate diplotype strings against per-read evidence.
/// Each candidate is scored through the exact sub-alleles the matcher found
/// for it (`raw_candidates`), not the main-allele definitions, so candidates
/// that explain the same variants differ only by how reads phase them. A
/// candidate with several raw realizations takes the best-scoring one.
/// Returns one CandidateScore per input candidate, with posteriors normalised
/// across the set. Returns an EMPTY vec to signal "cannot score" — no reads,
/// no sites, OR (FAIL-CLOSED) ANY candidate that has no matching raw
/// candidate or names an unknown allele — so the caller falls back to the
/// cascade's own choice. Scoring only the representable subset is forbidden:
/// it could make a confident wrong pick by dropping one side of the ambiguity.
pub fn score_candidates(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    var_list: &[String],
    snp_file: &str,
    star_combinations: &StarCombinations,
    candidates: &[String],
    raw_candidates: &[String],
) -> Vec<CandidateScore> {
    if candidates.is_empty() {
        return Vec::new();
    }
    let (per_candidate, cn) = match map_candidates_to_raw(candidates, raw_candidates) {
        Some(m) => m,
        None => return Vec::new(),
    };
    let lut = QualLookup::new();
    let sites = build_unified_sites(var_list, snp_file);
    if sites.is_empty() {
        return Vec::new();
    }
    let d6_alleles = build_d6_allele_defs(&star_combinations.dstar, &sites);
    let mut name_to_idx: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
    for (i, a) in d6_alleles.iter().enumerate() {
        name_to_idx.insert(a.name.as_str(), i);
    }
    let ((d6_start, d6_end), (d7_start, d7_end)) = fetch_windows(&sites);
    let all_reads = collect_reads_dual_region(
        reader, nchr, &sites, d6_start, d6_end, d7_start, d7_end,
    );
    let (reads, _nf) = filter_mismapped(all_reads, &d6_alleles, &sites, -1.5, &lut);
    if reads.is_empty() {
        return Vec::new();
    }

    let mut assignments: Vec<Vec<Vec<(usize, usize)>>> = Vec::with_capacity(candidates.len());
    for raws in &per_candidate {
        let mut cand_assignments = Vec::with_capacity(raws.len());
        for multiset in raws {
            match build_assignment(multiset, &name_to_idx, cn) {
                Some(a) => cand_assignments.push(a),
                None => return Vec::new(),
            }
        }
        assignments.push(cand_assignments);
    }

    // Each candidate takes its best-scoring raw realization; keep that
    // realization's per-read scores for the decisive-read count.
    let mut out: Vec<CandidateScore> = Vec::with_capacity(candidates.len());
    let mut per_read: Vec<Vec<f64>> = Vec::with_capacity(candidates.len());
    for (cand, cand_assignments) in candidates.iter().zip(assignments.iter()) {
        let mut best: Option<(f64, Vec<f64>)> = None;
        for a in cand_assignments {
            let reads_ll = read_log10_likelihoods(&reads, &d6_alleles, &sites, cn as u32, a, &lut);
            let ll: f64 = reads_ll.iter().sum();
            if best.as_ref().is_none_or(|(b, _)| ll > *b) {
                best = Some((ll, reads_ll));
            }
        }
        let (ll, reads_ll) = best.unwrap_or((f64::NEG_INFINITY, vec![0.0; reads.len()]));
        out.push(CandidateScore {
            diplotype: cand.clone(),
            log10_likelihood: ll,
            posterior: 0.0,
            decisive_reads: 0,
        });
        per_read.push(reads_ll);
    }
    for (c, n) in out.iter_mut().zip(count_decisive_reads(&per_read)) {
        c.decisive_reads = n;
    }

    let lls: Vec<f64> = out.iter().map(|c| c.log10_likelihood).collect();
    match posteriors_from_log10(&lls) {
        None => return Vec::new(), // no finite score — fail closed
        Some(posts) => {
            for (c, p) in out.iter_mut().zip(posts) {
                c.posterior = p;
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Fetch windows
// ---------------------------------------------------------------------------

/// CYP2D6 region coordinates (GRCh38).
const D6_START: i64 = 42123192;
const D6_END: i64 = 42132032;
const D7_START: i64 = 42139676;
const D7_END: i64 = 42145745;

/// Fetch windows (0-based, half-open) for the D6 and D7 regions, widened
/// from the gene windows above to cover every site's coordinate (some
/// downstream-exon9 D7 positions lie outside D7_START..D7_END).
fn fetch_windows(sites: &[UnifiedSite]) -> ((i64, i64), (i64, i64)) {
    let mut d6 = (D6_START, D6_END);
    let mut d7 = (D7_START, D7_END);
    for s in sites {
        d6 = (d6.0.min(s.d6_pos - 1), d6.1.max(s.d6_pos));
        if let Some(p) = s.d7_pos {
            d7 = (d7.0.min(p - 1), d7.1.max(p));
        }
    }
    (d6, d7)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_log10_lik_match() {
        let lut = QualLookup::new();
        // High-quality matching base should give score close to 0
        let score = base_log10_lik(b'A', b'A', 30, false, &lut);
        assert!(score > -0.01);
        assert!(score <= 0.0);
    }

    #[test]
    fn test_base_log10_lik_mismatch() {
        let lut = QualLookup::new();
        // High-quality mismatching base should give very negative score
        let score = base_log10_lik(b'A', b'T', 30, false, &lut);
        assert!(score < -2.0);
    }

    #[test]
    fn test_base_log10_lik_deletion() {
        let lut = QualLookup::new();
        // Deleted base matches deletion expectation
        let score_match = base_log10_lik(0, 0, 30, true, &lut);
        assert!(score_match > -0.01);

        // Present base does not match deletion
        let score_mismatch = base_log10_lik(b'A', 0, 30, true, &lut);
        assert!(score_mismatch < -2.0);
    }

    #[test]
    fn test_parse_target_variant_snp() {
        let site = parse_target_variant("g.42128945C>T").unwrap();
        assert_eq!(site.d6_pos, 42128945);
        assert_eq!(site.d6_ref_base, b'C');
        assert_eq!(site.d6_alt_base, Some(b'T'));
        assert!(!site.is_deletion);
    }

    #[test]
    fn test_parse_target_variant_del() {
        let site = parse_target_variant("g.42128945delC");
        // Single-base deletion
        assert!(site.is_some());
        let site = site.unwrap();
        assert_eq!(site.d6_ref_base, b'C');
        assert_eq!(site.d6_alt_base, Some(0));
        assert!(site.is_deletion);
    }

    #[test]
    fn test_parse_target_variant_multi_base_del() {
        // Multi-base deletion should be skipped (only single-base supported)
        let site = parse_target_variant("g.42128945delCT");
        assert!(site.is_none());
    }

    fn test_site(name: &str, pos: i64, d6: u8, d7: u8) -> UnifiedSite {
        UnifiedSite {
            d6_pos: pos, d7_pos: Some(pos + 1000), d6_ref_base: d6, d7_base: Some(d7),
            d6_alt_base: Some(d7), name: name.to_string(), is_deletion: false,
        }
    }

    #[test]
    fn test_read_log10_likelihoods_prefer_matching_reads() {
        let sites = vec![
            test_site("g.1A>T", 1, b'A', b'T'),
            test_site("g.2C>G", 2, b'C', b'G'),
        ];
        let allele_a = AlleleDef { name: "*A".into(), alt_sites: [0usize].into_iter().collect() };
        let allele_b = AlleleDef { name: "*B".into(), alt_sites: [1usize].into_iter().collect() };
        let alleles = vec![allele_a, allele_b];
        let reads: Vec<ReadObs> = (0..10)
            .map(|i| ReadObs { name: format!("r{i}"), observations: vec![(0, b'T', 40), (1, b'C', 40)] })
            .collect();
        let lut = QualLookup::new();
        let s_a: f64 = read_log10_likelihoods(&reads, &alleles, &sites, 2, &[(0usize, 2usize)], &lut).iter().sum();
        let s_b: f64 = read_log10_likelihoods(&reads, &alleles, &sites, 2, &[(1usize, 2usize)], &lut).iter().sum();
        assert!(s_a > s_b, "reads favor *A: {s_a} !> {s_b}");
    }

    #[test]
    fn test_build_assignment_fail_closed() {
        let mut idx: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
        idx.insert("*10", 0);
        idx.insert("*36", 1);
        let known = ms(&[("*10", 1), ("*36", 2)]);
        assert!(build_assignment(&known, &idx, 3).is_some());
        let mixed = ms(&[("*10", 1), ("*36", 1), ("*4.013", 1)]);
        assert!(build_assignment(&mixed, &idx, 3).is_none());
        assert!(build_assignment(&known, &idx, 2).is_none());
    }

    fn ms(v: &[(&str, usize)]) -> Vec<(String, usize)> {
        v.iter().map(|(n, c)| (n.to_string(), *c)).collect()
    }

    #[test]
    fn test_label_main_multiset() {
        let expected = ms(&[("*10", 1), ("*36", 2)]);
        assert_eq!(label_main_multiset("*36/*36+*10"), expected);
        assert_eq!(label_main_multiset("*10/*36x2"), expected);
        assert_eq!(label_main_multiset("*10.002_*36_*36"), expected);
    }

    fn strs(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn fetch_windows_cover_every_site() {
        use crate::depth_calling::utilities::get_var_names;
        let mut var_list = get_var_names(crate::data::TARGET_VARIANT);
        var_list.extend(get_var_names(crate::data::TARGET_VARIANT_HOMO));
        let sites = build_unified_sites(&var_list, crate::data::SNP_FILE);
        let ((d6s, d6e), (d7s, d7e)) = fetch_windows(&sites);
        for s in &sites {
            assert!(d6s < s.d6_pos && s.d6_pos <= d6e, "{} D6 {}", s.name, s.d6_pos);
            if let Some(p) = s.d7_pos {
                assert!(d7s < p && p <= d7e, "{} D7 {}", s.name, p);
            }
        }
    }

    #[test]
    fn decisive_reads_counts_reads_clearly_favoring_each_candidate() {
        // Per-candidate, per-read log10-likelihoods. A read counts for a
        // candidate only if it beats every other candidate by > 1 (10:1).
        let per_read = vec![
            vec![0.0, 0.0, 0.0, -0.5, -3.0],
            vec![-2.0, -1.5, 0.0, 0.0, 0.0],
        ];
        assert_eq!(count_decisive_reads(&per_read), vec![2, 1]);
    }

    #[test]
    fn candidates_map_to_their_matched_suballeles_with_d6_copy_count() {
        let (per_cand, cn) = map_candidates_to_raw(
            &strs(&["*1/*21", "*2/*21"]),
            &strs(&["*1_*21.002", "*2.001_*21"]),
        )
        .unwrap();
        assert_eq!(cn, 2, "D6 copies, not the D6+D7 total");
        assert_eq!(per_cand[0], vec![ms(&[("*1", 1), ("*21.002", 1)])]);
        assert_eq!(per_cand[1], vec![ms(&[("*2.001", 1), ("*21", 1)])]);
    }

    #[test]
    fn arrangements_of_one_multiset_share_the_same_raw_candidate() {
        let (per_cand, cn) = map_candidates_to_raw(
            &strs(&["*36/*36+*10", "*10/*36x2"]),
            &strs(&["*10.002_*36_*36"]),
        )
        .unwrap();
        assert_eq!(cn, 3);
        let expected = vec![ms(&[("*10.002", 1), ("*36", 2)])];
        assert_eq!(per_cand[0], expected);
        assert_eq!(per_cand[1], expected);
    }

    #[test]
    fn underscore_joined_candidates_are_mapped() {
        let (per_cand, cn) = map_candidates_to_raw(
            &strs(&["*10_*10_*4", "*10_*36_*4"]),
            &strs(&["*10_*10.002_*4.013", "*10_*36_*4"]),
        )
        .unwrap();
        assert_eq!(cn, 3);
        assert_eq!(per_cand[0], vec![ms(&[("*10", 1), ("*10.002", 1), ("*4.013", 1)])]);
        assert_eq!(per_cand[1], vec![ms(&[("*10", 1), ("*36", 1), ("*4", 1)])]);
    }

    #[test]
    fn candidate_without_a_matching_raw_fails_closed() {
        // *68 is part of the label but not of any raw candidate.
        assert!(map_candidates_to_raw(
            &strs(&["*68+*4/*2", "*68+*4/*1"]),
            &strs(&["*2_*4", "*1_*4"]),
        )
        .is_none());
    }

    #[test]
    fn test_split_copy_suffix() {
        assert_eq!(split_copy_suffix("*36"), ("*36".to_string(), 1));
        assert_eq!(split_copy_suffix("*36x2"), ("*36".to_string(), 2));
        assert_eq!(split_copy_suffix("*4.013"), ("*4.013".to_string(), 1));
        assert_eq!(split_copy_suffix("*1x"), ("*1x".to_string(), 1));
    }

    // E1 (decided by the human 2026-10-07): a failed BAM fetch must not be
    // read as evidence. The resolver only narrows `;`-ambiguous calls, so its
    // safe answer to a failed fetch is "no reads" -- not scorable, the call
    // stays ambiguous -- rather than scoring the half of the reads it did get.
    #[test]
    fn collect_reads_dual_region_returns_nothing_when_one_fetch_fails() {
        use crate::depth_calling::utilities::get_var_names;
        let bam = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../benchmark/HG00096_pgx_regions.bam");
        if !bam.exists() {
            eprintln!("skipping: benchmark BAM not present");
            return;
        }
        let mut var_list = get_var_names(crate::data::TARGET_VARIANT);
        var_list.extend(get_var_names(crate::data::TARGET_VARIANT_HOMO));
        let sites = build_unified_sites(&var_list, crate::data::SNP_FILE);
        let ((d6s, d6e), (d7s, d7e)) = fetch_windows(&sites);
        let mut reader = bam::IndexedReader::from_path(&bam).unwrap();

        let both = collect_reads_dual_region(&mut reader, "chr22", &sites, d6s, d6e, d7s, d7e);
        assert!(!both.is_empty(), "control: the real windows must yield reads");

        // D7 window inverted: its fetch fails, the D6 fetch succeeds.
        let d7_failed = collect_reads_dual_region(&mut reader, "chr22", &sites, d6s, d6e, d7e, d7s);
        assert!(
            d7_failed.is_empty(),
            "a failed D7 fetch must yield no reads, not the D6 half: got {}",
            d7_failed.len()
        );
    }
}

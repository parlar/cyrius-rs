//! Likelihood-based diplotype caller.
//!
//! Enumerates all valid diplotype combinations (constrained by copy number)
//! and scores each against per-read evidence. Reads are collected from both
//! the CYP2D6 and CYP2D7 regions; CYP2D7 is modelled as an additional
//! pseudo-allele so that mismapped reads are properly accounted for.
//!
//! The model for a single read under a diplotype (allele_i, allele_j) with
//! one copy of D7 always present:
//!
//!   P(read | i, j) = w_d7 * P(read | D7) + w_i * P(read | i) + w_j * P(read | j)
//!
//! where w_d7 = 1/(CN+1), w_d6 = CN/(CN+1) split equally between the two
//! D6 alleles (for CN=2: w_d7 = 1/3, w_i = w_j = 1/3).

use std::collections::{HashMap, HashSet};

use rayon::prelude::*;
use rust_htslib::bam::{self, Read as BamRead};

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
}

/// A scored diplotype candidate.
#[derive(Debug, Clone)]
pub struct ScoredDiplotype {
    pub allele_a: String,
    pub allele_b: String,
    /// Joint log10-likelihood across all reads.
    pub log10_likelihood: f64,
    /// Posterior probability, normalised over the retained top-N candidates
    /// (an approximation of the full posterior when mass lies outside top-N).
    pub posterior: f64,
}

/// Order scored diplotypes best-first: descending log10-likelihood, NaN-safe,
/// with a deterministic tiebreak on allele names so equal-likelihood
/// candidates (which are common) resolve identically across runs/threads.
fn sort_scored(a: &ScoredDiplotype, b: &ScoredDiplotype) -> std::cmp::Ordering {
    b.log10_likelihood
        .total_cmp(&a.log10_likelihood)
        .then_with(|| a.allele_a.cmp(&b.allele_a))
        .then_with(|| a.allele_b.cmp(&b.allele_b))
}

/// Full result from the diplotype caller.
#[derive(Debug, Clone)]
pub struct DiplotypeLikelihoodResult {
    /// Top diplotypes sorted by likelihood (descending).
    pub top_diplotypes: Vec<ScoredDiplotype>,
    /// Number of reads used (informative, covering ≥1 variant site).
    pub n_informative_reads: usize,
    /// Number of reads filtered as likely mismapped (low lik under all hypotheses).
    pub n_filtered_reads: usize,
    /// Total number of reads fetched before filtering.
    pub n_total_reads: usize,
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
    if rest.starts_with("del") {
        let del_seq = &rest[3..];
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
// D7 pseudo-allele + D6 allele definitions
// ---------------------------------------------------------------------------

/// Build the CYP2D7 pseudo-allele definition.
/// At D6/D7 diagnostic sites, D7 carries the D7 base (not the D6 ref).
/// At star-allele-defining sites, D7 carries the D6 reference (since those
/// variants are specific to D6 sub-alleles, not D6-vs-D7 differences).
fn build_d7_allele(_sites: &[UnifiedSite]) -> AlleleDef {
    // D7 doesn't carry any star-allele ALT variants — those are D6-internal.
    // At diagnostic sites where d7_base differs from d6_ref, we handle
    // scoring specially in score_read_for_allele().
    AlleleDef {
        name: "CYP2D7".to_string(),
        alt_sites: HashSet::new(), // D7 carries none of the D6 ALT alleles
    }
}

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
        if reader
            .fetch(bam::FetchDefinition::Region(tid as i32, start, end))
            .is_err()
        {
            continue;
        }

        for result in reader.records() {
            let mut record = match result {
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
// Diplotype enumeration + scoring
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


/// Score all diplotype pairs and return the top_n.
///
/// The mixture model includes one copy of D7 plus the D6 alleles:
///   P(read | i, j, copies_i, copies_j) =
///       w_d7 * P(read|D7) + w_i * P(read|i) + w_j * P(read|j)
///
/// For CN=2: w_d7 = 1/3, w_i = 1/3, w_j = 1/3.
/// For CN=3 with assignment (2 copies of i, 1 of j):
///   w_d7 = 1/4, w_i = 2/4 = 0.5, w_j = 1/4 = 0.25.
///
/// For CN >= 3, we try all valid copy-count splits and take the best.
fn score_diplotypes(
    reads: &[ReadObs],
    d6_alleles: &[AlleleDef],
    sites: &[UnifiedSite],
    total_cn: u32,
    top_n: usize,
    lut: &QualLookup,
) -> Vec<ScoredDiplotype> {
    let n_alleles = d6_alleles.len();
    let n_reads = reads.len();
    if n_reads == 0 || n_alleles == 0 {
        return Vec::new();
    }

    let cn = total_cn as usize;
    let denom = (cn + 1) as f64; // +1 for D7

    // Precompute per-read per-allele scores as flat array [read * n_alleles + allele]
    // Work in natural-log domain: exp() is ~3x faster than 10.0_f64.powf()
    let ln10 = std::f64::consts::LN_10;
    let mut d6_scores_ln = vec![0.0f64; n_reads * n_alleles];
    for (r, obs) in reads.iter().enumerate() {
        let base = r * n_alleles;
        for (a, allele) in d6_alleles.iter().enumerate() {
            d6_scores_ln[base + a] = score_read_d6(obs, allele, sites, lut) * ln10;
        }
    }

    // Precompute per-read D7 scores (ln domain)
    let d7_scores_ln: Vec<f64> = reads
        .iter()
        .map(|obs| score_read_d7(obs, sites, lut) * ln10)
        .collect();

    let ln_w_d7 = (1.0 / denom).ln();

    // Build all (i, j) pairs for parallel iteration
    let pairs: Vec<(usize, usize)> = (0..n_alleles)
        .flat_map(|i| (i..n_alleles).map(move |j| (i, j)))
        .collect();

    // Score diplotype pairs in parallel using thread-local top-N lists
    let thread_results: Vec<Vec<ScoredDiplotype>> = pairs
        .par_chunks(256)
        .map(|chunk| {
            let mut local_best: Vec<ScoredDiplotype> = Vec::new();
            let mut local_min = f64::NEG_INFINITY;

            for &(i, j) in chunk {
                // A heterozygous pair needs at least one copy of each allele.
                // When cn < 2 no such split exists, so skip it rather than
                // leaving best_total at -inf (which would poison the posterior).
                if i != j && cn < 2 {
                    continue;
                }
                let splits: Vec<(usize, usize)> = if i == j {
                    vec![(cn, 0)]
                } else {
                    (1..cn).map(|ci| (ci, cn - ci)).collect()
                };

                let mut best_total = f64::NEG_INFINITY;
                for &(copies_i, copies_j) in &splits {
                    let ln_w_i = (copies_i as f64 / denom).ln();
                    let ln_w_j = if copies_j > 0 {
                        (copies_j as f64 / denom).ln()
                    } else {
                        f64::NEG_INFINITY
                    };

                    let mut total = 0.0;
                    for r in 0..n_reads {
                        let base = r * n_alleles;
                        let s_d7 = ln_w_d7 + d7_scores_ln[r];
                        let s_i = ln_w_i + d6_scores_ln[base + i];
                        let s_j = if copies_j > 0 {
                            ln_w_j + d6_scores_ln[base + j]
                        } else {
                            f64::NEG_INFINITY
                        };

                        // log-sum-exp in natural log (exp() is much faster than powf())
                        let max_s = s_d7.max(s_i).max(s_j);
                        let lse = max_s
                            + ((s_d7 - max_s).exp()
                                + (s_i - max_s).exp()
                                + (s_j - max_s).exp())
                            .ln();
                        total += lse;
                    }

                    if total > best_total {
                        best_total = total;
                    }
                }

                // Convert back to log10 for storage
                let best_total_log10 = best_total / ln10;

                if local_best.len() < top_n || best_total_log10 > local_min {
                    local_best.push(ScoredDiplotype {
                        allele_a: d6_alleles[i].name.clone(),
                        allele_b: d6_alleles[j].name.clone(),
                        log10_likelihood: best_total_log10,
                        posterior: 0.0,
                    });
                    local_best.sort_by(|a, b| sort_scored(a, b));
                    local_best.truncate(top_n);
                    local_min = local_best
                        .last()
                        .map_or(f64::NEG_INFINITY, |x| x.log10_likelihood);
                }
            }
            local_best
        })
        .collect();

    // Merge thread-local results
    let mut best: Vec<ScoredDiplotype> = thread_results.into_iter().flatten().collect();
    best.sort_by(|a, b| sort_scored(a, b));
    best.truncate(top_n);

    // Normalise posteriors over the top-N (approximate)
    if !best.is_empty() {
        let lls: Vec<f64> = best.iter().map(|d| d.log10_likelihood).collect();
        if let Some(posts) = posteriors_from_log10(&lls) {
            for (d, p) in best.iter_mut().zip(posts) {
                d.posterior = p;
            }
        }
    }

    best
}

// ---------------------------------------------------------------------------
// Shared posterior normalization helper
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

/// Score one explicit diplotype given as (allele_index, copies) pairs summing
/// to total_cn. General k-allele mixture: one D7 copy + the specified D6
/// alleles, each weighted by copies/(cn+1). Returns log10-likelihood.
fn score_explicit_multiset(
    reads: &[ReadObs],
    d6_alleles: &[AlleleDef],
    sites: &[UnifiedSite],
    total_cn: u32,
    assignment: &[(usize, usize)],
    lut: &QualLookup,
) -> f64 {
    let ln10 = std::f64::consts::LN_10;
    let denom = (total_cn as usize + 1) as f64;
    let ln_w_d7 = (1.0 / denom).ln();
    let mut total = 0.0f64;
    for obs in reads {
        let s_d7 = ln_w_d7 + score_read_d7(obs, sites, lut) * ln10;
        // Per-read Vec allocation is fine here: the resolver scores only a tiny
        // explicit candidate set, not the full O(n^2) enumeration of score_diplotypes.
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
            continue;
        }
        let lse = max_s + terms.iter().map(|t| (t - max_s).exp()).sum::<f64>().ln();
        total += lse;
    }
    total / ln10
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
    let all_reads = collect_reads_dual_region(
        reader, nchr, &sites, D6_START, D6_END, D7_START, D7_END,
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

    let mut out: Vec<CandidateScore> = Vec::with_capacity(candidates.len());
    for (cand, cand_assignments) in candidates.iter().zip(assignments.iter()) {
        let ll = cand_assignments
            .iter()
            .map(|a| score_explicit_multiset(&reads, &d6_alleles, &sites, cn as u32, a, &lut))
            .fold(f64::NEG_INFINITY, f64::max);
        out.push(CandidateScore { diplotype: cand.clone(), log10_likelihood: ll, posterior: 0.0 });
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
// Public entry point
// ---------------------------------------------------------------------------

/// CYP2D6 region coordinates (GRCh38).
const D6_START: i64 = 42123192;
const D6_END: i64 = 42132032;
const D7_START: i64 = 42139676;
const D7_END: i64 = 42145745;

/// Run the likelihood-based diplotype caller.
///
/// This is an independent caller that uses per-read evidence from both D6
/// and D7 regions. It does not replace the existing pileup-based caller
/// but can be used alongside it for validation or as a primary caller.
pub fn call_diplotype(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    var_list: &[String],
    snp_file: &str,
    star_combinations: &StarCombinations,
    total_cn: u32,
    top_n: usize,
) -> DiplotypeLikelihoodResult {
    // Build qual lookup table once
    let lut = QualLookup::new();

    // 1. Build unified site list
    let sites = build_unified_sites(var_list, snp_file);
    log::info!(
        "diplotype_caller: {} unified sites ({} with D7 info)",
        sites.len(),
        sites.iter().filter(|s| s.d7_base.is_some()).count(),
    );

    if sites.is_empty() {
        return DiplotypeLikelihoodResult {
            top_diplotypes: Vec::new(),
            n_informative_reads: 0,
            n_filtered_reads: 0,
            n_total_reads: 0,
        };
    }

    // 2. Build allele definitions
    let d6_alleles = build_d6_allele_defs(&star_combinations.dstar, &sites);
    let _d7_allele = build_d7_allele(&sites);

    log::info!(
        "diplotype_caller: {} D6 alleles + 1 D7 pseudo-allele",
        d6_alleles.len(),
    );

    // 3. Collect reads from both regions
    let all_reads = collect_reads_dual_region(
        reader, nchr, &sites, D6_START, D6_END, D7_START, D7_END,
    );
    let n_total = all_reads.len();

    log::info!("diplotype_caller: {} reads collected from D6+D7 regions", n_total);

    // 4. Filter mismapped reads
    // Threshold: average log10-lik per site must be > -1.5
    // (i.e., average error rate per site < ~3%)
    let (reads, n_filtered) = filter_mismapped(all_reads, &d6_alleles, &sites, -1.5, &lut);
    let n_informative = reads.len();

    log::info!(
        "diplotype_caller: {} informative reads, {} filtered as mismapped",
        n_informative,
        n_filtered,
    );

    if n_informative == 0 {
        return DiplotypeLikelihoodResult {
            top_diplotypes: Vec::new(),
            n_informative_reads: 0,
            n_filtered_reads: n_filtered,
            n_total_reads: n_total,
        };
    }

    // 5. Score diplotypes
    let top = score_diplotypes(&reads, &d6_alleles, &sites, total_cn, top_n, &lut);

    // Log top 5
    for (i, d) in top.iter().take(5).enumerate() {
        log::info!(
            "diplotype_caller: #{} {}/{} ll={:.1} posterior={:.4}",
            i + 1,
            d.allele_a,
            d.allele_b,
            d.log10_likelihood,
            d.posterior,
        );
    }

    DiplotypeLikelihoodResult {
        top_diplotypes: top,
        n_informative_reads: n_informative,
        n_filtered_reads: n_filtered,
        n_total_reads: n_total,
    }
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

    #[test]
    fn test_d7_allele_has_no_alt_sites() {
        let sites = vec![
            UnifiedSite {
                d6_pos: 100,
                d7_pos: Some(200),
                d6_ref_base: b'A',
                d7_base: Some(b'G'),
                d6_alt_base: Some(b'T'),
                name: "test".to_string(),
                is_deletion: false,
            },
        ];
        let d7 = build_d7_allele(&sites);
        assert!(d7.alt_sites.is_empty());
        assert_eq!(d7.name, "CYP2D7");
    }

    fn test_site(name: &str, pos: i64, d6: u8, d7: u8) -> UnifiedSite {
        UnifiedSite {
            d6_pos: pos, d7_pos: Some(pos + 1000), d6_ref_base: d6, d7_base: Some(d7),
            d6_alt_base: Some(d7), name: name.to_string(), is_deletion: false,
        }
    }

    #[test]
    fn test_score_explicit_multiset_prefers_matching_reads() {
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
        let s_a = score_explicit_multiset(&reads, &alleles, &sites, 2, &[(0usize, 2usize)], &lut);
        let s_b = score_explicit_multiset(&reads, &alleles, &sites, 2, &[(1usize, 2usize)], &lut);
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
}

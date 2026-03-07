//! Physical phasing disambiguator for ambiguous diplotype calls.
//!
//! When match_star returns "more_than_one_match" for CN=2, this module uses
//! read-level phasing to pick the correct diplotype. It examines which defining
//! variants co-occur on the same reads (cis) vs different reads (trans).

use crate::depth_calling::snp_count::passing_read;
use crate::types::StarCombinations;
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::{HashMap, HashSet};

/// Variant info from the target variant file, mapping variant names to pileup positions
/// and ALT/REF base patterns.
#[derive(Debug, Clone)]
pub struct VariantLookup {
    /// variant_name -> (pileup_pos_0based, alt_bases, ref_bases)
    entries: HashMap<String, (i64, String, String)>,
}

impl VariantLookup {
    /// Parse the target variant file content into a lookup table.
    pub fn from_target_file(content: &str) -> Self {
        let mut entries = HashMap::new();
        for line in content.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 7 {
                continue;
            }
            let pos: i64 = match fields[1].parse::<i64>() {
                Ok(p) => p - 1, // convert 1-based to 0-based
                Err(_) => continue,
            };
            let alt_bases = fields[2].to_string();
            // REF may have comma-separated alternatives; take the first
            let ref_bases = fields[4].split(',').next().unwrap_or("").to_string();
            let variant_name = fields[6].to_string();
            entries.insert(variant_name, (pos, alt_bases, ref_bases));
        }
        VariantLookup { entries }
    }

    /// Look up a variant by name. Returns (pileup_pos_0based, alt_bases, ref_bases).
    fn get(&self, name: &str) -> Option<&(i64, String, String)> {
        self.entries.get(name)
    }
}

/// A diplotype option with raw (suballele) star names and clean (main allele) names.
#[derive(Debug, Clone)]
struct Diplotype {
    raw_hap1: String,
    raw_hap2: String,
    clean_hap1: String,
    clean_hap2: String,
}

/// Parse raw candidates (e.g., ["*1_*21.002", "*2.001_*21"]) and clean call
/// (e.g., "*1/*21;*2/*21") into diplotype options.
fn parse_diplotypes(raw_candidates: &[String], clean_call: &str) -> Vec<Diplotype> {
    let clean_parts: Vec<&str> = clean_call.split(';').collect();

    // Build diplotypes from raw candidates, matching with clean counterparts
    let mut diplotypes = Vec::new();
    for (i, raw) in raw_candidates.iter().enumerate() {
        let raw_parts: Vec<&str> = raw.split('_').collect();
        if raw_parts.len() != 2 {
            continue;
        }
        let (clean_h1, clean_h2) = if i < clean_parts.len() {
            let cp: Vec<&str> = clean_parts[i].split('/').collect();
            if cp.len() == 2 {
                (cp[0].to_string(), cp[1].to_string())
            } else {
                continue;
            }
        } else {
            // Fallback: strip suballele suffix
            (
                raw_parts[0].split('.').next().unwrap().to_string(),
                raw_parts[1].split('.').next().unwrap().to_string(),
            )
        };
        diplotypes.push(Diplotype {
            raw_hap1: raw_parts[0].to_string(),
            raw_hap2: raw_parts[1].to_string(),
            clean_hap1: clean_h1,
            clean_hap2: clean_h2,
        });
    }
    diplotypes
}

/// Get the set of defining variant names for a star allele (using raw suballele name).
/// Returns variant names like ["g.42128945C>T", "g.42129033G>A"].
fn get_star_variants(star: &str, dstar: &HashMap<String, String>) -> Vec<String> {
    if let Some(var_key) = dstar.get(star) {
        if var_key == "NA" {
            return vec![];
        }
        return var_key.split('_').map(|s| s.to_string()).collect();
    }
    // *1 has no defining variants (it's the reference)
    if star == "*1" {
        return vec![];
    }
    vec![]
}


/// Resolved variant info for phasing: pileup position (0-based) + ALT/REF patterns.
#[derive(Debug, Clone)]
struct ResolvedVariant {
    pos_0based: i64,
    alt_bases: String,
    ref_bases: String,
}

/// Check if two variants at nearby positions co-occur on the same reads (cis)
/// or on different reads (trans).
///
/// When `use_read_pairs` is true and positions are >120bp apart, performs
/// two separate fetches so that mate pairs spanning the distance can both
/// contribute observations. Read pairs share the same qname, so a mate
/// covering pos1 and its partner covering pos2 will be linked.
///
/// Returns (cis_count, trans_count).
fn check_phasing(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    rv1: &ResolvedVariant,
    rv2: &ResolvedVariant,
    use_read_pairs: bool,
) -> (usize, usize) {
    let pos1 = rv1.pos_0based;
    let pos2 = rv2.pos_0based;
    let tid = match reader.header().tid(nchr.as_bytes()) {
        Some(t) => t,
        None => return (0, 0),
    };

    let mut read_alleles: HashMap<String, (Option<bool>, Option<bool>)> = HashMap::new();

    let distance = (pos1 - pos2).unsigned_abs();
    let use_two_fetches = use_read_pairs && distance > 120;

    if use_two_fetches {
        // Fetch reads overlapping pos1. Use a wide window upstream because
        // reads starting up to ~150bp before pos1 can still cover it.
        let read_len = 150i64;
        reader
            .fetch(bam::FetchDefinition::Region(
                tid as i32,
                (pos1 - read_len).max(0),
                pos1 + 1,
            ))
            .unwrap();
        for record_result in reader.records() {
            let record = match record_result {
                Ok(r) => r,
                Err(_) => continue,
            };
            if record.is_secondary() || record.is_supplementary() || record.is_duplicate() {
                continue;
            }
            if !passing_read(&record, false, false) {
                continue;
            }
            let read_name = String::from_utf8_lossy(record.qname()).to_string();
            let seq = record.seq().as_bytes();
            let allele1 = check_allele_pileup(&record, &seq, rv1);
            let allele2 = check_allele_pileup(&record, &seq, rv2);
            let entry = read_alleles.entry(read_name).or_insert((None, None));
            if allele1.is_some() {
                entry.0 = allele1;
            }
            if allele2.is_some() {
                entry.1 = allele2;
            }
        }

        // Fetch reads overlapping pos2
        reader
            .fetch(bam::FetchDefinition::Region(
                tid as i32,
                (pos2 - read_len).max(0),
                pos2 + 1,
            ))
            .unwrap();
        for record_result in reader.records() {
            let record = match record_result {
                Ok(r) => r,
                Err(_) => continue,
            };
            if record.is_secondary() || record.is_supplementary() || record.is_duplicate() {
                continue;
            }
            if !passing_read(&record, false, false) {
                continue;
            }
            let read_name = String::from_utf8_lossy(record.qname()).to_string();
            let seq = record.seq().as_bytes();
            let allele1 = check_allele_pileup(&record, &seq, rv1);
            let allele2 = check_allele_pileup(&record, &seq, rv2);
            let entry = read_alleles.entry(read_name).or_insert((None, None));
            if allele1.is_some() {
                entry.0 = allele1;
            }
            if allele2.is_some() {
                entry.1 = allele2;
            }
        }
    } else {
        // Single fetch spanning both positions
        let start = pos1.min(pos2) - 1;
        let end = pos1.max(pos2);
        reader
            .fetch(bam::FetchDefinition::Region(tid as i32, start, end))
            .unwrap();

        for record_result in reader.records() {
            let record = match record_result {
                Ok(r) => r,
                Err(_) => continue,
            };
            if record.is_secondary() || record.is_supplementary() || record.is_duplicate() {
                continue;
            }
            if !passing_read(&record, false, false) {
                continue;
            }
            let read_name = String::from_utf8_lossy(record.qname()).to_string();
            let seq = record.seq().as_bytes();
            let allele1 = check_allele_pileup(&record, &seq, rv1);
            let allele2 = check_allele_pileup(&record, &seq, rv2);
            let entry = read_alleles.entry(read_name).or_insert((None, None));
            if allele1.is_some() {
                entry.0 = allele1;
            }
            if allele2.is_some() {
                entry.1 = allele2;
            }
        }
    }

    let mut cis = 0;
    let mut trans = 0;
    for (has_alt1, has_alt2) in read_alleles.values() {
        match (has_alt1, has_alt2) {
            (Some(true), Some(true)) => cis += 1,   // both ALT → same haplotype
            (Some(true), Some(false)) => trans += 1, // one ALT, one REF → different haplotypes
            (Some(false), Some(true)) => trans += 1,
            // both REF: uninformative for phasing — skip
            _ => {}
        }
    }

    (cis, trans)
}

/// Check whether a read carries the ALT allele at a pileup position.
/// Reads multiple bases from the read starting at rv.pos_0based and compares
/// against ALT and REF patterns from the target variant file.
/// Returns Some(true) for ALT, Some(false) for REF, None if position not covered.
fn check_allele_pileup(
    record: &rust_htslib::bam::Record,
    seq: &[u8],
    rv: &ResolvedVariant,
) -> Option<bool> {
    let check_len = rv.alt_bases.len().max(rv.ref_bases.len());
    let bases = read_bases_at_pos(record, seq, rv.pos_0based, check_len)?;
    let bases_upper = bases.to_uppercase();
    let alt_upper = rv.alt_bases.to_uppercase();
    let ref_upper = rv.ref_bases.to_uppercase();

    if bases_upper.starts_with(&alt_upper) || alt_upper.starts_with(&bases_upper) {
        Some(true)
    } else if bases_upper.starts_with(&ref_upper) || ref_upper.starts_with(&bases_upper) {
        Some(false)
    } else {
        None // doesn't match either pattern
    }
}

/// Read `n` bases from a read starting at a 0-based genomic position.
/// Walks the CIGAR to find the query offset, then extracts bases.
fn read_bases_at_pos(
    record: &rust_htslib::bam::Record,
    seq: &[u8],
    genomic_pos_0based: i64,
    n: usize,
) -> Option<String> {
    let record_start = record.pos(); // 0-based
    let ref_offset = genomic_pos_0based - record_start;
    if ref_offset < 0 {
        return None;
    }

    let mut ref_pos = 0i64;
    let mut query_pos = 0usize;

    for &cigar_elem in record.cigar().iter() {
        use rust_htslib::bam::record::Cigar::*;
        match cigar_elem {
            Match(len) | Equal(len) | Diff(len) => {
                let len = len as i64;
                if ref_pos + len > ref_offset {
                    let qp = query_pos + (ref_offset - ref_pos) as usize;
                    if qp + n <= seq.len() {
                        let bases: String = seq[qp..qp + n]
                            .iter()
                            .map(|&b| b as char)
                            .collect();
                        return Some(bases);
                    } else if qp < seq.len() {
                        let bases: String = seq[qp..seq.len()]
                            .iter()
                            .map(|&b| b as char)
                            .collect();
                        return Some(bases);
                    }
                    return None;
                }
                ref_pos += len;
                query_pos += len as usize;
            }
            Ins(len) => { query_pos += len as usize; }
            Del(len) => {
                let len = len as i64;
                if ref_pos + len > ref_offset { return None; }
                ref_pos += len;
            }
            SoftClip(len) => { query_pos += len as usize; }
            HardClip(_) | Pad(_) => {}
            RefSkip(len) => { ref_pos += len as i64; }
        }
    }
    None
}

/// Attempt to disambiguate an ambiguous CN=2 diplotype call using physical phasing.
///
/// `raw_candidates` are the raw star allele combinations (with suballeles), e.g.,
/// ["*1_*21.002", "*2.001_*21"]. `clean_call` is the main-allele form like "*1/*21;*2/*21".
///
/// Returns `Some(resolved_call)` if disambiguation succeeds, `None` otherwise.
pub fn disambiguate(
    clean_call: &str,
    raw_candidates: &[String],
    star_combinations: &StarCombinations,
    variant_lookup: &VariantLookup,
    bam_path: &str,
    nchr: &str,
    reference_fasta: Option<&str>,
    index_name: Option<&str>,
    use_read_pairs: bool,
) -> Option<String> {
    let diplotypes = parse_diplotypes(raw_candidates, clean_call);
    if diplotypes.len() != 2 {
        log::debug!(
            "phase_disambiguate: expected 2 diplotype options, got {} (raw: {:?})",
            diplotypes.len(),
            raw_candidates
        );
        return None;
    }

    let dt1 = &diplotypes[0];
    let dt2 = &diplotypes[1];

    let dstar = &star_combinations.dstar;

    // Get variants for each haplotype using raw (suballele) names
    let vars_dt1_h1 = get_star_variants(&dt1.raw_hap1, dstar);
    let vars_dt1_h2 = get_star_variants(&dt1.raw_hap2, dstar);
    let vars_dt2_h1 = get_star_variants(&dt2.raw_hap1, dstar);
    let vars_dt2_h2 = get_star_variants(&dt2.raw_hap2, dstar);

    log::debug!(
        "phase_disambiguate: dt1={}({})/{}({}) dt2={}({})/{}({})",
        dt1.raw_hap1,
        vars_dt1_h1.join(","),
        dt1.raw_hap2,
        vars_dt1_h2.join(","),
        dt2.raw_hap1,
        vars_dt2_h1.join(","),
        dt2.raw_hap2,
        vars_dt2_h2.join(","),
    );

    // Collect all unique variants
    let all_vars: HashSet<&str> = vars_dt1_h1
        .iter()
        .chain(vars_dt1_h2.iter())
        .chain(vars_dt2_h1.iter())
        .chain(vars_dt2_h2.iter())
        .map(|s| s.as_str())
        .collect();

    if all_vars.len() < 2 {
        log::debug!("phase_disambiguate: too few variants ({}) to phase", all_vars.len());
        return None;
    }

    let dt1_h1_set: HashSet<&str> = vars_dt1_h1.iter().map(|s| s.as_str()).collect();
    let dt1_h2_set: HashSet<&str> = vars_dt1_h2.iter().map(|s| s.as_str()).collect();
    let dt2_h1_set: HashSet<&str> = vars_dt2_h1.iter().map(|s| s.as_str()).collect();
    let dt2_h2_set: HashSet<&str> = vars_dt2_h2.iter().map(|s| s.as_str()).collect();

    // Find pairs of variants where the two options predict different phasing
    let all_vars_vec: Vec<&str> = all_vars.into_iter().collect();
    let mut phasing_pairs: Vec<(&str, &str, bool)> = Vec::new(); // (var_a, var_b, cis_in_dt1)

    for i in 0..all_vars_vec.len() {
        for j in (i + 1)..all_vars_vec.len() {
            let va = all_vars_vec[i];
            let vb = all_vars_vec[j];

            let cis_dt1 = are_cis_in_diplotype(va, vb, &dt1_h1_set, &dt1_h2_set);
            let cis_dt2 = are_cis_in_diplotype(va, vb, &dt2_h1_set, &dt2_h2_set);

            match (cis_dt1, cis_dt2) {
                (Some(true), Some(false)) => phasing_pairs.push((va, vb, true)),
                (Some(false), Some(true)) => phasing_pairs.push((va, vb, false)),
                _ => {} // Same phasing in both options, or indeterminate
            }
        }
    }

    if phasing_pairs.is_empty() {
        log::debug!("phase_disambiguate: no differentiating variant pairs found");
        return None;
    }

    log::info!(
        "phase_disambiguate: found {} differentiating variant pairs for {}",
        phasing_pairs.len(),
        clean_call
    );

    // Open BAM and check phasing for each pair
    let mut reader = match crate::depth_calling::utilities::open_alignment_file_with_index(
        bam_path,
        reference_fasta,
        index_name,
    ) {
        Ok(r) => r,
        Err(e) => {
            log::warn!("phase_disambiguate: failed to open BAM: {}", e);
            return None;
        }
    };

    let mut dt1_score = 0i32;
    let mut dt2_score = 0i32;

    for &(va, vb, cis_in_dt1) in &phasing_pairs {
        let rv_a = match variant_lookup.get(va) {
            Some((pos, alt, ref_b)) => ResolvedVariant {
                pos_0based: *pos,
                alt_bases: alt.clone(),
                ref_bases: ref_b.clone(),
            },
            None => {
                log::debug!("phase_disambiguate: no lookup for variant {}", va);
                continue;
            }
        };
        let rv_b = match variant_lookup.get(vb) {
            Some((pos, alt, ref_b)) => ResolvedVariant {
                pos_0based: *pos,
                alt_bases: alt.clone(),
                ref_bases: ref_b.clone(),
            },
            None => {
                log::debug!("phase_disambiguate: no lookup for variant {}", vb);
                continue;
            }
        };

        let (cis_count, trans_count) =
            check_phasing(&mut reader, nchr, &rv_a, &rv_b, use_read_pairs);

        if cis_count + trans_count == 0 {
            continue;
        }

        log::debug!(
            "phase_disambiguate: {}(pos={}) vs {}(pos={}): cis={}, trans={}, expected_cis_dt1={}",
            va, rv_a.pos_0based, vb, rv_b.pos_0based, cis_count, trans_count, cis_in_dt1
        );

        // Only count a pair's evidence when the imbalance is clear (>2:1 ratio)
        let total = cis_count + trans_count;
        let majority = cis_count.max(trans_count);
        if total < 3 || (majority as f64) / (total as f64) < 0.7 {
            log::debug!(
                "phase_disambiguate: skipping weak pair (cis={}, trans={})",
                cis_count, trans_count
            );
            continue;
        }

        if cis_count > trans_count {
            if cis_in_dt1 {
                dt1_score += (cis_count - trans_count) as i32;
            } else {
                dt2_score += (cis_count - trans_count) as i32;
            }
        } else {
            if cis_in_dt1 {
                dt2_score += (trans_count - cis_count) as i32;
            } else {
                dt1_score += (trans_count - cis_count) as i32;
            }
        }
    }

    log::info!(
        "phase_disambiguate: dt1_score={}, dt2_score={} for {}",
        dt1_score, dt2_score, clean_call
    );

    // Require clear winner: score >= 3 and at least 3x the other
    let min_score = 3;
    if dt1_score >= min_score && (dt2_score == 0 || dt1_score >= 3 * dt2_score) {
        let resolved = format!("{}/{}", dt1.clean_hap1, dt1.clean_hap2);
        log::info!("phase_disambiguate: resolved {} -> {}", clean_call, resolved);
        Some(resolved)
    } else if dt2_score >= min_score && (dt1_score == 0 || dt2_score >= 3 * dt1_score) {
        let resolved = format!("{}/{}", dt2.clean_hap1, dt2.clean_hap2);
        log::info!("phase_disambiguate: resolved {} -> {}", clean_call, resolved);
        Some(resolved)
    } else {
        log::info!("phase_disambiguate: inconclusive for {}", clean_call);
        None
    }
}

/// Check if two variants are on the same haplotype (cis) in a given diplotype.
/// Returns Some(true) for cis, Some(false) for trans, None if indeterminate.
fn are_cis_in_diplotype(
    var_a: &str,
    var_b: &str,
    hap1_vars: &HashSet<&str>,
    hap2_vars: &HashSet<&str>,
) -> Option<bool> {
    let a_on_h1 = hap1_vars.contains(var_a);
    let a_on_h2 = hap2_vars.contains(var_a);
    let b_on_h1 = hap1_vars.contains(var_b);
    let b_on_h2 = hap2_vars.contains(var_b);

    // Homozygous variant (on both haplotypes) is not informative for phasing
    if (a_on_h1 && a_on_h2) || (b_on_h1 && b_on_h2) {
        return None;
    }

    match (a_on_h1, b_on_h1, a_on_h2, b_on_h2) {
        (true, true, false, false) => Some(true),  // Both on hap1
        (false, false, true, true) => Some(true),  // Both on hap2
        (true, false, false, true) => Some(false), // a on hap1, b on hap2
        (false, true, true, false) => Some(false), // a on hap2, b on hap1
        _ => None,
    }
}

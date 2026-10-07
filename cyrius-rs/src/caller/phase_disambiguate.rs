//! Physical phasing disambiguator for ambiguous diplotype calls.
//!
//! When match_star returns "more_than_one_match" for CN=2, this module uses
//! read-level phasing to pick the correct diplotype. It examines which defining
//! variants co-occur on the same reads (cis) vs different reads (trans).

use crate::caller::match_star_allele::main_allele_name;
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
    /// NOTE: this is a module-local variant-name parser with rules specific to
    /// physical phasing; it is not a general-purpose parser.
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

    // Build diplotypes from raw candidates, matching each with the clean
    // counterpart that has the same main alleles. The two lists are sorted
    // differently, so pairing by position would mislabel them.
    let mut diplotypes = Vec::new();
    for raw in raw_candidates {
        let raw_parts: Vec<&str> = raw.split('_').collect();
        if raw_parts.len() != 2 {
            continue;
        }
        let main1 = main_allele_name(raw_parts[0]);
        let main2 = main_allele_name(raw_parts[1]);
        let mut raw_mains = [main1, main2];
        raw_mains.sort();
        let matched = clean_parts.iter().find_map(|part| {
            let cp: Vec<&str> = part.split('/').collect();
            if cp.len() != 2 {
                return None;
            }
            let mut clean_mains = [cp[0], cp[1]];
            clean_mains.sort();
            (clean_mains == raw_mains).then(|| (cp[0].to_string(), cp[1].to_string()))
        });
        // Fallback: the stripped suballele names
        let (clean_h1, clean_h2) =
            matched.unwrap_or_else(|| (main1.to_string(), main2.to_string()));
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
    // Anything not in the table (including *1, the reference) has no defining
    // variants for phasing purposes.
    vec![]
}


/// Resolved variant info for phasing: pileup position (0-based) + ALT/REF patterns.
#[derive(Debug, Clone)]
struct ResolvedVariant {
    pos_0based: i64,
    alt_bases: String,
    ref_bases: String,
}

/// Half-open BAM fetch window `[start, end)` covering both variant positions
/// for the single-fetch phasing path.
///
/// `end` is one past the later position: the region is half-open, so stopping
/// at the position itself drops any read whose leftmost mapped base is exactly
/// the later variant.
///
/// `start` is `min(pos) - 1`, one base wider than it needs to be. An htslib
/// region query returns every read *overlapping* `[beg, end)`, not only reads
/// starting inside it — a query at a single position comes back with reads
/// that start ~150 bp earlier — so `start = min(pos)` would already fetch
/// everything covering the earlier variant. The `-1` is a harmless leftover of
/// the very misconception that produced the `end` bug. It is kept deliberately:
/// narrowing a live fetch window is a behaviour change, not a comment fix.
///
/// `start` is clamped at 0 so a variant at position 0 cannot produce a negative
/// coordinate. That clamp is unreachable on real CYP2D6 data, where positions
/// are ~42 Mbp; it guards a bad coordinate from a malformed variant file, not
/// any live case.
///
/// Returns `(start, end)`.
fn single_fetch_window(pos1: i64, pos2: i64) -> (i64, i64) {
    let start = (pos1.min(pos2) - 1).max(0);
    let end = pos1.max(pos2) + 1;
    (start, end)
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
        if reader
            .fetch(bam::FetchDefinition::Region(
                tid as i32,
                (pos1 - read_len).max(0),
                pos1 + 1,
            ))
            .is_err()
        {
            return (0, 0);
        }
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
        if reader
            .fetch(bam::FetchDefinition::Region(
                tid as i32,
                (pos2 - read_len).max(0),
                pos2 + 1,
            ))
            .is_err()
        {
            return (0, 0);
        }
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
        let (start, end) = single_fetch_window(pos1, pos2);
        if reader
            .fetch(bam::FetchDefinition::Region(tid as i32, start, end))
            .is_err()
        {
            return (0, 0);
        }

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

    // Find pairs of variants where the two options predict different phasing.
    //
    // `all_vars` is a HashSet, and std's RandomState seeds its iteration order
    // per process, so draining it unsorted makes both the order of
    // `phasing_pairs` and the (va, vb) assignment within each pair vary from
    // run to run. Every consumer downstream happens to be order-symmetric
    // today — `are_cis_in_diplotype` is symmetric in a/b, `single_fetch_window`
    // is symmetric, the two-fetch branch merges both fetches so swapping
    // rv1/rv2 only swaps `entry.0`/`entry.1` and cis/trans survive that swap,
    // and the scores are order-independent sums — so the scores and the
    // resolved call are stable. The debug log is not, and any future
    // order-sensitive tiebreak would flip calls between runs. Sort for a
    // reproducible order.
    let mut all_vars_vec: Vec<&str> = all_vars.into_iter().collect();
    all_vars_vec.sort_unstable();
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_diplotypes_labels_each_raw_candidate_with_its_own_clean_form() {
        // Raw order (*1 before *10) differs from clean order ("*10_" sorts
        // before "*1_"), so pairing by position would swap the labels.
        let raw = vec!["*1_*4.004".to_string(), "*10_*4.012".to_string()];
        let dts = parse_diplotypes(&raw, "*10/*4;*1/*4");
        assert_eq!(dts.len(), 2);
        let by_raw: HashMap<&str, (&str, &str)> = dts
            .iter()
            .map(|d| (d.raw_hap1.as_str(), (d.clean_hap1.as_str(), d.clean_hap2.as_str())))
            .collect();
        assert_eq!(by_raw["*1"], ("*1", "*4"));
        assert_eq!(by_raw["*10"], ("*10", "*4"));
    }

    #[test]
    fn single_fetch_window_reaches_past_the_later_variant() {
        // The fetch region is half-open, so `end` must be one past the later
        // position; at `end == pos` a read whose leftmost mapped base IS the
        // later variant falls outside the window and is never fetched.
        let (start, end) = single_fetch_window(42_126_611, 42_126_666);
        // Exact equality: `start <= 42_126_610` also passes for -1, so it does
        // not guard the clamp, and `end > 42_126_666` also passes for
        // 99_999_999.
        assert_eq!(
            start, 42_126_610,
            "start must reach exactly one base before the earlier variant, got {start}"
        );
        assert_eq!(
            end, 42_126_667,
            "end must be exactly one past the later variant 42126666 so a read \
             starting exactly there is still fetched, got {end}"
        );
        // Argument order must not matter.
        assert_eq!(single_fetch_window(42_126_666, 42_126_611), (start, end));
    }

    #[test]
    fn single_fetch_window_clamps_a_negative_start() {
        // A variant at position 0 must not produce a negative fetch start.
        let (start, end) = single_fetch_window(0, 7);
        assert_eq!(start, 0, "start must be clamped to 0, got {start}");
        assert_eq!(end, 8, "end must be one past the later variant, got {end}");
    }

    /// 0-based position of the earlier of the fixture's two phasing variants.
    const FIXTURE_POS_EARLY: i64 = 1000;
    /// 0-based position of the later one: 100 bp on, so `distance <= 120` and
    /// `check_phasing` takes the single-fetch path under either flag setting.
    const FIXTURE_POS_LATE: i64 = 1100;
    /// 0-based start of the mate. With CIGAR 50M its last mapped base is
    /// exactly `FIXTURE_POS_EARLY`.
    const FIXTURE_MATE_POS: i64 = 951;

    /// Write an indexed BAM holding one two-record pair for the single-fetch
    /// phasing window.
    ///
    /// Record 1 (the mate) is 50M from `FIXTURE_MATE_POS`, so its last mapped
    /// base is exactly `FIXTURE_POS_EARLY`, where it carries the ALT. Record 2
    /// starts at exactly `FIXTURE_POS_LATE` — the later variant — and carries
    /// the ALT there. They share one qname, so the two observations land in the
    /// same `read_alleles` entry and tally as a single cis read, but only if
    /// the fetch window actually reaches record 2.
    ///
    /// NOTE: `bam::Record::new()` leaves `flags = 4` (unmapped), and htslib's
    /// `bam_endpos` ignores the CIGAR for an unmapped record, so an unflagged
    /// record spans only `[pos, pos + 1)` and silently drops out of region
    /// queries. The `set_flags` calls below are load-bearing.
    fn write_phasing_pair_bam(dir: &std::path::Path) -> std::path::PathBuf {
        let mut header = bam::header::Header::new();
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr22");
        sq.push_tag(b"LN", 50_818_468);
        header.push_record(&sq);
        let path = dir.join("phasing_pair.bam");
        {
            let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
            let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(50)]);
            let qual = vec![30u8; 50];
            let span = FIXTURE_POS_LATE + 50 - FIXTURE_MATE_POS;

            // The mate: ALT 'C' against REF 'G' at FIXTURE_POS_EARLY.
            let mut mate_seq = vec![b'G'; 50];
            mate_seq[(FIXTURE_POS_EARLY - FIXTURE_MATE_POS) as usize] = b'C';
            let mut mate = bam::Record::new();
            mate.set(b"pair1", Some(&cigar), &mate_seq, &qual);
            mate.set_tid(0);
            mate.set_pos(FIXTURE_MATE_POS);
            mate.set_mapq(60);
            mate.set_mtid(0);
            mate.set_mpos(FIXTURE_POS_LATE);
            mate.set_insert_size(span);
            // PAIRED | PROPER_PAIR | MATE_REVERSE — mapped, primary, not a dup.
            mate.set_flags(0x1 | 0x2 | 0x20);
            writer.write(&mate).unwrap();

            // The read starting exactly at the later variant: ALT 'T' against
            // REF 'A' at FIXTURE_POS_LATE, query offset 0.
            let mut read_seq = vec![b'A'; 50];
            read_seq[0] = b'T';
            let mut read = bam::Record::new();
            read.set(b"pair1", Some(&cigar), &read_seq, &qual);
            read.set_tid(0);
            read.set_pos(FIXTURE_POS_LATE);
            read.set_mapq(60);
            read.set_mtid(0);
            read.set_mpos(FIXTURE_MATE_POS);
            read.set_insert_size(-span);
            // PAIRED | PROPER_PAIR | REVERSE
            read.set_flags(0x1 | 0x2 | 0x10);
            writer.write(&read).unwrap();
        }
        bam::index::build(
            &path,
            None::<&std::path::PathBuf>,
            bam::index::Type::Bai,
            1,
        )
        .unwrap();
        path
    }

    #[test]
    fn check_phasing_fetches_a_read_starting_exactly_at_the_later_variant() {
        // Finding #52, end to end. The single-fetch region used to be
        // [min(pos) - 1, max(pos)); it is half-open, so the read whose leftmost
        // mapped base is exactly the later variant was never fetched. Its
        // mate's ALT at the earlier variant then sat alone in the entry as
        // (Some(true), None) and was dropped as uninformative, so the pair
        // tallied (0, 0) instead of (1, 0).
        //
        // The two window-arithmetic tests above cannot catch that: reverting
        // the call site to `let end = pos1.max(pos2);` leaves both of them
        // green and only this test red. That is what makes the fix
        // demonstrably non-vacuous.
        let dir = std::env::temp_dir()
            .join(format!("cyrius_phase_window_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bam_path = write_phasing_pair_bam(&dir);
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();

        let rv_early = ResolvedVariant {
            pos_0based: FIXTURE_POS_EARLY,
            alt_bases: "C".to_string(),
            ref_bases: "G".to_string(),
        };
        let rv_late = ResolvedVariant {
            pos_0based: FIXTURE_POS_LATE,
            alt_bases: "T".to_string(),
            ref_bases: "A".to_string(),
        };

        let single = check_phasing(&mut reader, "chr22", &rv_early, &rv_late, false);
        // 100 bp apart, so --phase-readpair takes the same single-fetch path.
        let readpair = check_phasing(&mut reader, "chr22", &rv_early, &rv_late, true);

        // Clean up before the assertions so a failure still leaves nothing behind.
        drop(reader);
        let _ = std::fs::remove_dir_all(&dir);

        assert_eq!(
            single,
            (1, 0),
            "the ALT-carrying read starting exactly at {FIXTURE_POS_LATE} plus \
             its mate's ALT at {FIXTURE_POS_EARLY} are one cis observation; got \
             (cis, trans) = {single:?}, which is what a window ending at the \
             later position returns"
        );
        assert_eq!(
            readpair,
            (1, 0),
            "the two variants are only 100 bp apart, so --phase-readpair must \
             take the same single-fetch path; got {readpair:?}"
        );
    }
}

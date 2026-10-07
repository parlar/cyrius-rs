use crate::depth_calling::copy_number_call::{
    call_reg1_cn, process_raw_call_denovo, process_raw_call_gc, CnProbResult,
};
use crate::depth_calling::haplotype::{
    extract_hap, get_haplotypes_from_bam, get_haplotypes_from_bam_single_region,
};
use crate::depth_calling::snp_count::{
    get_supporting_reads, get_supporting_reads_single_region,
};
use crate::fisher;
use crate::stats;
use crate::types::{CnRegions, SnpLookup};
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::HashMap;

const INTRON1_BP_APPROX: i64 = 42130500;
const EXON9_BP_APPROX: i64 = 42126611;
const P_CUTOFF: f64 = 0.05;

/// The CNV tag vocabulary: (tag, total_cn, exon9_and_downstream, exon9_to_intron1,
/// intron1_upstream). This is the single source of truth for both the accepted-tag
/// vocabulary (`CNV_ACCEPTED`) and the CN-pattern lookup table
/// (`cnvtag_lookup_table`), so the two can never drift apart (finding #26).
pub const CNVTAG_ENTRIES: [(&str, u32, u32, u32, u32); 28] = [
    ("star5_star5", 2, 0, 0, 0),
    ("star13_star13", 2, 2, 0, 0),
    ("star13intron1_star13intron1", 2, 2, 2, 0),
    ("star5", 3, 1, 1, 1),
    ("star13", 3, 2, 1, 1),
    ("star13intron1", 3, 2, 2, 1),
    ("star5_star5_star68", 3, 0, 0, 1),
    ("star5_star68", 4, 1, 1, 2),
    ("cn2", 4, 2, 2, 2),
    ("exon9hyb_star5", 4, 2, 2, 2),
    ("dup_star13", 4, 3, 2, 2),
    ("dup_star13intron1", 4, 3, 3, 2),
    ("star13_star68", 4, 2, 1, 2),
    ("cn3", 5, 3, 3, 3),
    ("exon9hyb", 5, 2, 3, 3),
    ("star68", 5, 2, 2, 3),
    ("cn4", 6, 4, 4, 4),
    ("exon9hyb_exon9hyb", 6, 2, 4, 4),
    ("star68_star68", 6, 2, 2, 4),
    ("dup_exon9hyb", 6, 3, 4, 4),
    ("dup_star68", 6, 3, 3, 4),
    ("exon9hyb_star68", 6, 2, 3, 4),
    ("cn5", 7, 5, 5, 5),
    ("exon9hyb_exon9hyb_exon9hyb", 7, 2, 5, 5),
    ("star68_star68_star68", 7, 2, 2, 5),
    ("cn6", 8, 6, 6, 6),
    ("exon9hyb_exon9hyb_exon9hyb_exon9hyb", 8, 2, 6, 6),
    ("star68_star68_star68_star68", 8, 2, 2, 6),
];

/// Pull just the tag names out of a `CNVTAG_ENTRIES`-shaped array, preserving order.
const fn tags_of<'a, const N: usize>(entries: &[(&'a str, u32, u32, u32, u32); N]) -> [&'a str; N] {
    let mut tags = [""; N];
    let mut i = 0;
    while i < N {
        tags[i] = entries[i].0;
        i += 1;
    }
    tags
}

/// The accepted CNV tag vocabulary, derived from `CNVTAG_ENTRIES` so it can never
/// list a tag the lookup table doesn't also know about, or vice versa (finding #26).
pub const CNV_ACCEPTED: [&str; 28] = tags_of(&CNVTAG_ENTRIES);

/// CNVTAG lookup table mapping cnv tags to CN regions.
pub fn cnvtag_lookup_table() -> HashMap<String, CnRegions> {
    let mut m = HashMap::new();
    for (tag, total, exon9, exon9_intron1, intron1) in CNVTAG_ENTRIES {
        m.insert(
            tag.to_string(),
            CnRegions {
                total_cn: total,
                exon9_and_downstream: exon9,
                exon9_to_intron1: exon9_intron1,
                intron1_upstream: intron1,
            },
        );
    }
    m
}

/// Clean variant sites
pub const CLEAN_VAR: &[&str] = &[
    "g.42129809T>C",
    "g.42129819G>T",
    "g.42128945C>T",
    "g.42126611C>G",
    "g.42130692G>A",
    "g.42127941G>A",
];

/// Noisy variant sites (may have strand bias)
pub const NOISY_VAR: &[&str] = &[
    "g.42127473C>T",
    "g.42129042T>C",
    "g.42129174C>A",
    "g.42129180A>T",
    "g.42127526C>T",
    "g.42128325A>G",
    "g.42126877G>A",
    "g.42127973T>C",
    "g.42127556T>C",
];

/// Get total expected CN at each variant site based on CNV configuration.
pub fn get_total_cn_per_site(
    cnvtag: &str,
    var_db: &SnpLookup,
    var_list: &[String],
) -> Option<Vec<u32>> {
    let lookup = cnvtag_lookup_table();
    let num_var_sites = var_db.dsnp1.len();
    let variant_names = &var_list[..num_var_sites];

    let mut num_before_intron1 = 0usize;
    let mut num_after_exon9 = 0usize;
    for var_name in variant_names {
        if var_name.len() >= 10 {
            let var_pos: i64 = var_name[2..10].parse().unwrap_or(0);
            if var_pos >= INTRON1_BP_APPROX {
                num_before_intron1 += 1;
            }
            if var_pos <= EXON9_BP_APPROX {
                num_after_exon9 += 1;
            }
        }
    }

    let cn_pattern = lookup.get(cnvtag)?;
    let mut cn_list = Vec::new();
    for _ in 0..num_after_exon9 {
        cn_list.push(cn_pattern.exon9_and_downstream);
    }
    for _ in 0..(num_var_sites - num_before_intron1 - num_after_exon9) {
        cn_list.push(cn_pattern.exon9_to_intron1);
    }
    for _ in 0..num_before_intron1 {
        cn_list.push(cn_pattern.intron1_upstream);
    }
    Some(cn_list)
}

/// Call CN for SNP sites between CYP2D6 and CYP2D7.
pub fn call_cn_snp(
    total_cn: u32,
    lsnp1: &[usize],
    lsnp2: &[usize],
    threshold: f64,
) -> Vec<Option<u32>> {
    let cn_prob: Vec<CnProbResult> = lsnp1
        .iter()
        .zip(lsnp2.iter())
        .map(|(&c1, &c2)| {
            call_reg1_cn(Some(total_cn), c1 as f64, c2 as f64, 0.0)[0].clone().into()
        })
        .collect();
    process_raw_call_gc(&cn_prob, threshold, true)
}

/// Call CN for variant sites in homology regions.
pub fn call_cn_var_homo(
    total_cn: u32,
    lsnp1: &[usize],
    lsnp2: &[usize],
) -> Vec<Option<u32>> {
    let cn_prob: Vec<CnProbResult> = lsnp1
        .iter()
        .zip(lsnp2.iter())
        .map(|(&c1, &c2)| call_reg1_cn(Some(total_cn), c1 as f64, c2 as f64, 4.0)[0].clone().into())
        .collect();
    let calls = process_raw_call_denovo(&cn_prob, 0.8, 0.65, None, true);
    calls
        .into_iter()
        .map(|c| c.map(|v| v.min(total_cn.saturating_sub(2))))
        .collect()
}

/// Call CN for variant sites in non-homology regions.
pub fn call_cn_var(
    cnvtag: &str,
    var_alt: &mut Vec<usize>,
    var_ref: &[usize],
    alt_forward: &[usize],
    alt_reverse: &[usize],
    var_list: &[String],
    var_db: &SnpLookup,
) -> Vec<Option<u32>> {
    let total_cn = get_total_cn_per_site(cnvtag, var_db, var_list).unwrap();
    let mut cn_prob: Vec<CnProbResult> = Vec::new();

    for (i, &forward) in alt_forward.iter().enumerate() {
        let reverse = alt_reverse[i];
        let total_ref = var_ref[i];
        let mut total_var = var_alt[i];

        if total_var > 0 && NOISY_VAR.contains(&var_list[i].as_str()) {
            let ntotal = (forward + reverse) as u64;
            let half = ntotal as f64 / 2.0;
            let (_, pvalue) =
                fisher::fisher_exact([[forward as u64, reverse as u64], [half as u64, half as u64]]);
            if pvalue < P_CUTOFF || forward <= 1 || reverse <= 1 {
                total_var = 0;
                var_alt[i] = 0;
            }
        }

        let min_read = if CLEAN_VAR.contains(&var_list[i].as_str()) {
            2
        } else if NOISY_VAR.contains(&var_list[i].as_str()) {
            7
        } else {
            4
        };

        cn_prob.push(
            call_reg1_cn(Some(total_cn[i]), total_var as f64, total_ref as f64, min_read as f64)[0]
                .clone()
                .into(),
        );
    }

    process_raw_call_denovo(&cn_prob, 0.8, 0.65, Some(&total_cn), true)
}

/// Define read filters for sequence search.
fn good_read(record: &bam::Record) -> bool {
    !record.is_secondary() && !record.is_supplementary()
}

/// Search for insertions at 42128936 defining *30/*40/*58.
pub fn get_allele_counts_var42128936(
    reader: &mut bam::IndexedReader,
    genome: &str,
) -> (usize, usize, usize) {
    let (chr, start, end) = if genome == "chr38" {
        ("chr22", 42128848i64, 42128978i64)
    } else {
        ("22", 42128848i64, 42128978i64)
    };

    let mut long_ins_read = 0usize;
    let mut short_ins_read = 0usize;
    let mut ref_read = 0usize;

    let tid = match reader.header().tid(chr.as_bytes()) {
        Some(t) => t,
        None => return (0, 0, 0),
    };
    if reader
        .fetch(bam::FetchDefinition::Region(tid as i32, start, end))
        .is_err()
    {
        log::warn!(
            "get_allele_counts_var42128936: failed to fetch region {}:{}-{}",
            chr,
            start,
            end
        );
        return (0, 0, 0);
    }

    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        if result.is_err() {
            continue;
        }
        if good_read(&record) {
            let seq = String::from_utf8_lossy(&record.seq().as_bytes()).to_string();
            if seq.contains("TGGGGCGAAAGGGGCGAAAGGGGCGAAAGGGGCGT") {
                long_ins_read += 1;
            } else if seq.contains("TTGGGGCGAAAGGGGCGAAAGGGGCGTC") {
                short_ins_read += 1;
            } else if seq.contains("TTGGGGCGAAAGGGGCGTC") {
                ref_read += 1;
            }
        }
    }
    (ref_read, long_ins_read, short_ins_read)
}

/// Update variant read counts for g42128936.
pub fn update_var42128936(
    var_list: &[String],
    var_alt: &mut Vec<usize>,
    var_ref: &mut Vec<usize>,
    ref_read: usize,
    long_ins_read: usize,
    short_ins_read: usize,
) {
    if let Some(idx) = var_list
        .iter()
        .position(|v| v == "g.42128936-42128937insGGGGCGAAAGGGGCGAAA")
    {
        var_alt[idx] = long_ins_read;
        var_ref[idx] = short_ins_read + ref_read;
    }
    if let Some(idx) = var_list
        .iter()
        .position(|v| v == "g.42128936-42128937insGGGGCGAAA")
    {
        var_alt[idx] = short_ins_read;
        var_ref[idx] = long_ins_read + ref_read;
    }
}

/// Call exon 9 conversion.
pub fn call_exon9gc(
    d6_count: &[usize],
    d7_count: &[usize],
    full_length_cn: Option<u32>,
) -> Option<u32> {
    let cn = full_length_cn? as u64;

    let mut d6_values = Vec::new();
    let mut cn_prob: Vec<CnProbResult> = Vec::new();
    for (i, &c1) in d6_count.iter().enumerate() {
        let c2 = d7_count[i];
        let sum = c1 + c2;
        // A zero-depth site carries no evidence about the D6 fraction at this
        // position; skip it rather than recording a phantom D6 value of 0.0
        // (which would otherwise always satisfy the min_d6 < 1.2 guard below).
        if sum > 0 {
            d6_values.push(cn as f64 * c1 as f64 / sum as f64);
        }
        cn_prob.push(call_reg1_cn(Some(cn as u32), c1 as f64, c2 as f64, 3.0)[0].clone().into());
    }

    let cn_prob_stringent = process_raw_call_gc(&cn_prob, 0.88, true);
    let cn_calls: Vec<u32> = cn_prob_stringent
        .iter()
        .filter_map(|&a| a)
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();

    if cn_calls.len() == 1 {
        // Round (rather than truncate) the mean counts before handing them to
        // call_reg1_cn, which casts its count arguments with `as u64`. A raw
        // mean like 44.5 would otherwise be silently truncated to 44.
        let mean_c1 = stats::mean(
            &d6_count.iter().map(|&x| x as f64).collect::<Vec<_>>(),
        )
        .round();
        let mean_c2 = stats::mean(
            &d7_count.iter().map(|&x| x as f64).collect::<Vec<_>>(),
        )
        .round();
        // Relaxed from 0.75 to 0.70 to handle *10D cases where the duplication
        // breakpoint downstream of exon9 causes borderline D6/D7 ratios at exon9
        // (e.g., HG00421: average probability 0.729 was failing the 0.75 cutoff).
        let ave_call = process_raw_call_gc(
            &[call_reg1_cn(Some(cn as u32), mean_c1, mean_c2, 3.0)[0].clone().into()],
            0.70,
            true,
        );
        if let Some(ave_cn) = ave_call[0] {
            if ave_cn == cn_calls[0] {
                if ave_cn == 1 {
                    let min_d6 = d6_values.iter().cloned().fold(f64::INFINITY, f64::min);
                    let max_d6 = d6_values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    if min_d6 < 1.2 && max_d6 < 1.3 {
                        return Some(ave_cn);
                    }
                } else {
                    return Some(ave_cn);
                }
            }
        }
    }
    None
}

/// Call variant g.42126938C>T based on read depth and phased haplotypes.
pub fn call_var42126938(
    reader: &mut bam::IndexedReader,
    full_length_cn: u32,
    base_db: &SnpLookup,
) -> (Vec<usize>, Vec<String>, bool) {
    let mut var_called = Vec::new();
    let mut g_haplotype = false;

    let (snp_d6, snp_d7, _) = get_supporting_reads(reader, base_db);
    if snp_d6.is_empty() || snp_d7.is_empty() {
        return (vec![], Vec::new(), false);
    }
    let d6_d7_base_count = vec![*snp_d6.last().unwrap(), *snp_d7.last().unwrap()];

    let d6_cn = call_cn_snp(full_length_cn, &[d6_d7_base_count[0]], &[d6_d7_base_count[1]], 0.8)[0];

    if let Some(cn) = d6_cn {
        if full_length_cn >= 2 && cn < full_length_cn - 2 {
            let n_positions = base_db.dsnp1.len();
            let positions: Vec<usize> = (0..n_positions).collect();
            let haplotype_per_read = get_haplotypes_from_bam(reader, base_db, &positions);
            let recombinant_read_count = extract_hap(&haplotype_per_read, &[0, 2]);
            if let Some(count_12) = recombinant_read_count.get("12") {
                if count_12.iter().sum::<u32>() > 1 {
                    let g_hap_count = extract_hap(&haplotype_per_read, &[1, 2]);
                    for _ in 0..(full_length_cn - 2 - cn) {
                        var_called.push("g.42126938C>T".to_string());
                    }
                    if let Some(count_g) = g_hap_count.get("12") {
                        if count_g.iter().sum::<u32>() > 1 {
                            g_haplotype = true;
                        }
                    }
                }
            }
        }
    }

    (d6_d7_base_count, var_called, g_haplotype)
}

/// Call variant g.42127526C>T based on read depth and phased haplotypes.
pub fn call_var42127526_var42127556(
    reader: &mut bam::IndexedReader,
    cnvtag: &str,
    base_db: &SnpLookup,
) -> (Vec<usize>, Vec<usize>, Vec<String>) {
    let mut var_called = Vec::new();
    let lookup = cnvtag_lookup_table();

    let (var_ref, var_alt, _var_ref_forward, _var_ref_reverse) =
        get_supporting_reads_single_region(reader, base_db, None);

    // This caller expects exactly the two diagnostic sites (42127526, 42127556).
    // If the panel returns fewer, bail out with no calls rather than panicking.
    if var_ref.len() < 2 || var_alt.len() < 2 {
        return (Vec::new(), Vec::new(), Vec::new());
    }

    let var7526_count = vec![var_ref[0], var_alt[0]];
    let var7556_count = vec![var_ref[1], var_alt[1]];

    if let Some(cn_pattern) = lookup.get(cnvtag) {
        let d6_cn = cn_pattern.exon9_to_intron1;
        let var7526_cn =
            call_cn_snp(d6_cn, &[var7526_count[1]], &[var7526_count[0]], 0.6)[0].unwrap_or(0);
        let var7556_cn =
            call_cn_snp(d6_cn, &[var7556_count[1]], &[var7556_count[0]], 0.6)[0].unwrap_or(0);

        let n_positions = base_db.dsnp1.len();
        let positions: Vec<usize> = (0..n_positions).collect();
        let haplotype_per_read =
            get_haplotypes_from_bam_single_region(reader, base_db, &positions);
        let recombinant_read_count = extract_hap(&haplotype_per_read, &[0, 1, 2]);

        let has_211 = recombinant_read_count
            .get("211")
            .map_or(false, |c| c.iter().sum::<u32>() > 1);
        let has_221 = recombinant_read_count
            .get("221")
            .map_or(false, |c| c.iter().sum::<u32>() > 1);

        if has_211 {
            for _ in 0..var7526_cn {
                var_called.push("g.42127526C>T".to_string());
            }
        } else if has_221 {
            let min_cn = std::cmp::min(var7526_cn, var7556_cn);
            for _ in 0..min_cn {
                var_called.push("g.42127526C>T".to_string());
                var_called.push("g.42127556T>C".to_string());
            }
        }
    }

    (var7526_count, var7556_count, var_called)
}

/// Call haplotype with regard to g.42127803C>T and g.42127941G>A.
pub fn call_var42127803hap(
    reader: &mut bam::IndexedReader,
    cnvtag: &str,
    base_db: &SnpLookup,
) -> bool {
    let mut diff_haplotype = false;
    if cnvtag == "cn2" {
        let n_positions = base_db.dsnp1.len();
        let positions: Vec<usize> = (0..n_positions).collect();
        let haplotype_per_read =
            get_haplotypes_from_bam_single_region(reader, base_db, &positions);
        let recombinant_read_count = extract_hap(&haplotype_per_read, &[0, 1]);
        if let (Some(c12), Some(c21)) = (
            recombinant_read_count.get("12"),
            recombinant_read_count.get("21"),
        ) {
            if c12.iter().sum::<u32>() > 1 && c21.iter().sum::<u32>() > 1 {
                diff_haplotype = true;
            }
        }
    }
    diff_haplotype
}

/// Call haplotype for g.42130655-42130656insA (for *15.003).
pub fn call_var42130655ins_a(
    reader: &mut bam::IndexedReader,
    _full_length_cn: u32,
    base_db: &SnpLookup,
) -> (Vec<usize>, Vec<String>) {
    let mut var_called = Vec::new();
    let length_between_d67_regions: i64 = 13000;
    let mut no_aa_on_d6: usize;
    let mut no_a_on_d6: usize = 0;

    let (var_aa_d6, var_aa_d7, _, _) = get_supporting_reads_single_region(
        reader,
        &SnpLookup {
            dsnp1: base_db.dsnp2.clone(),
            dsnp2: base_db.dsnp2.clone(),
            nchr: base_db.nchr.clone(),
            dindex: base_db.dindex.clone(),
        },
        Some(length_between_d67_regions),
    );
    no_aa_on_d6 = *var_aa_d7.last().unwrap_or(&0);

    let has_long_insert_size_reads = *var_aa_d7.last().unwrap_or(&0) >= 5
        && *var_aa_d6.last().unwrap_or(&0) <= 2;

    if has_long_insert_size_reads {
        let (var_aa_d6_2, _var_aa_d7_2, _, _) = get_supporting_reads_single_region(
            reader,
            &SnpLookup {
                dsnp1: base_db.dsnp2.clone(),
                dsnp2: base_db.dsnp2.clone(),
                nchr: base_db.nchr.clone(),
                dindex: base_db.dindex.clone(),
            },
            None,
        );
        no_aa_on_d6 += *var_aa_d6_2.last().unwrap_or(&0);

        let (var_a_d6, _var_a_d7, _, _) =
            get_supporting_reads_single_region(reader, base_db, None);
        no_a_on_d6 = *var_a_d6.last().unwrap_or(&0);

        if no_a_on_d6 >= no_aa_on_d6 {
            var_called.push("g.42130655-42130656insA".to_string());
        } else {
            var_called.push("g.42130655-42130656insA".to_string());
            var_called.push("g.42130655-42130656insA".to_string());
        }
    }

    (vec![no_a_on_d6, no_aa_on_d6], var_called)
}

/// Return called variants based on called copy number and list of variant names.
pub fn get_called_variants(
    var_list: &[String],
    cn_prob_processed: &[Option<u32>],
    starting_index: usize,
) -> Vec<String> {
    let mut total_callset = Vec::new();
    for (i, &cn_called) in cn_prob_processed.iter().enumerate() {
        if let Some(cn) = cn_called {
            if cn != 0 {
                for _ in 0..cn {
                    total_callset.push(var_list[i + starting_index].clone());
                }
            }
        }
    }
    total_callset
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Finding #24: a zero-depth exon9 site was pushed into d6_values as 0.0,
    // which trivially satisfies the `min_d6 < 1.2` guard for the CN-1 case.
    // Site 8 has no reads; site 7 calls CN 1 with a D6 value of 1.25. Two live
    // sites at 1.25 would correctly fail the min check, so this must be None.
    #[test]
    fn exon9gc_zero_depth_site_does_not_vote_as_d6_cn_zero() {
        let d6_count = vec![25usize, 0usize];
        let d7_count = vec![15usize, 0usize];
        assert_eq!(call_exon9gc(&d6_count, &d7_count, Some(2)), None);
    }

    // Characterization test (not a red-to-green regression test for this
    // task's change): when no exon9 site has any depth at all, the call must
    // be None. This already held before this task touched the function —
    // `cn_prob` is all `None`, so `cn_calls` has length 0 (not 1), and the
    // `cn_calls.len() == 1` gate a few lines up falls through to the final
    // `None` without ever reaching the per-site `d6_values` vector. Kept
    // because the "no depth anywhere -> None" behavior is worth pinning, even
    // though nothing in this task's diff is what makes it true.
    #[test]
    fn exon9gc_no_depth_anywhere_falls_through_to_none() {
        let d6_count = vec![0usize, 0usize];
        let d7_count = vec![0usize, 0usize];
        assert_eq!(call_exon9gc(&d6_count, &d7_count, Some(2)), None);
    }

    // Finding #24: the average D6/D7 counts were truncated (`as u64`) rather
    // than rounded before being scored. Here mean_c1 = (11+2)/2 = 6.5 and
    // mean_c2 = (0+2)/2 = 1.0; site 0 alone (11 vs 0 reads) confidently calls
    // CN 2. Truncating the mean yields no confident average call (None);
    // rounding it to the nearest integer yields Some(2), matching site 0.
    #[test]
    fn exon9gc_rounds_mean_counts_instead_of_truncating() {
        let d6_count = vec![11usize, 2usize];
        let d7_count = vec![0usize, 2usize];
        assert_eq!(call_exon9gc(&d6_count, &d7_count, Some(2)), Some(2));
    }

    // Finding #26: CNV_ACCEPTED (main.rs) and cnvtag_lookup_table()'s keys used to be
    // two independently maintained lists that happened to agree. Both now derive from
    // the single CNVTAG_ENTRIES const below, so the three tests that follow pin the
    // vocabulary's shape rather than re-asserting the (now tautological) "do the two
    // lists match" fact.

    // Pins the vocabulary size as a literal so an accidental addition or deletion in
    // CNVTAG_ENTRIES shows up as a test failure rather than silently changing which
    // cnvtags main.rs's no-call gate accepts.
    #[test]
    fn cnv_accepted_has_the_documented_28_tags() {
        assert_eq!(CNV_ACCEPTED.len(), 28);
        assert_eq!(CNVTAG_ENTRIES.len(), 28);
    }

    // Detects a duplicate tag in CNVTAG_ENTRIES: HashMap::insert silently collapses a
    // repeated key, so the map would end up shorter than the entries list even though
    // CNV_ACCEPTED (derived from the same list) still carries the duplicate twice.
    // The call_cn_var unwrap (finding #26) depends on every accepted tag having
    // exactly one row, so a collapse here is exactly the failure mode to catch.
    #[test]
    fn cnvtag_lookup_table_has_one_entry_per_const_entry() {
        assert_eq!(cnvtag_lookup_table().len(), CNVTAG_ENTRIES.len());
    }

    // The invariant finding #26 actually cares about: get_total_cn_per_site's
    // lookup.get(cnvtag)? must resolve for every tag main.rs's CNV_ACCEPTED gate lets
    // through, or call_cn_var's unwrap (call_variants.rs:173) panics.
    #[test]
    fn every_accepted_cnvtag_resolves_in_lookup_table() {
        let table = cnvtag_lookup_table();
        for tag in CNV_ACCEPTED {
            assert!(table.get(tag).is_some(), "accepted tag {tag:?} has no lookup table row");
        }
    }
}

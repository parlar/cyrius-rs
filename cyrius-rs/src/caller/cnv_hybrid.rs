use crate::types::CnConsensus;

const REP_END_POSITION: usize = 3;
const EXON9_END_POSITION: usize = 9;
const INTRON4_BP_POSITION: usize = 40;
const INTRON1_BP_POSITION: usize = 74;
const REP_SITES_CN: u32 = 2;
const EXON9_TO_INTRON4_SITES_MIN: usize = 18;
const INTRON4_TO_INTRON1_SITES_MIN: usize = 19;
const INTRON1_UPSTREAM_SITES_MIN: usize = 25;
const EXON9REGION_SITES_MIN: usize = 4;
const INTRON1_UPSTREAM_SITES_MIN_LOOSE: usize = 15;

/// Get the consensus (mode) CN from a list of calls, requiring a minimum count.
///
/// Matches Python's `sorted(Counter(values).items(), key=count, reverse=True)`:
/// `Counter` preserves first-insertion order and Python's `sorted` is stable,
/// so on a tie the value seen first in `values` wins. We mirror that here by
/// iterating `values` in order to build a first-seen-ordered (value, count)
/// list, then stably sorting by count descending.
fn get_consensus(values: &[u32], min_count: usize) -> Option<u32> {
    if values.is_empty() {
        return None;
    }
    let mut counts: Vec<(u32, usize)> = Vec::new();
    for &v in values {
        if let Some(entry) = counts.iter_mut().find(|(k, _)| *k == v) {
            entry.1 += 1;
        } else {
            counts.push((v, 1));
        }
    }
    counts.sort_by(|a, b| b.1.cmp(&a.1));
    if counts[0].1 >= min_count {
        Some(counts[0].0)
    } else {
        None
    }
}

/// Return a tag for the called CNV/hybrid group based on detected CN switching point.
pub fn get_cnvtag(
    total_cn: u32,
    _raw_d6_cn: &[f64],
    cn_call_per_site: &[Option<u32>],
    exon9gc_call_stringent: Option<u32>,
    spacer_cn: Option<u32>,
) -> (Option<String>, CnConsensus) {
    // exon9 to intron4 sites
    let exon9_intron4_sites: Vec<u32> = cn_call_per_site
        [EXON9_END_POSITION..INTRON4_BP_POSITION.min(cn_call_per_site.len())]
        .iter()
        .filter_map(|&a| a)
        .collect();
    let mut exon9_intron4_sites_consensus = get_consensus(&exon9_intron4_sites, EXON9_TO_INTRON4_SITES_MIN);

    // intron4 to intron1 sites
    let intron4_intron1_sites: Vec<u32> = cn_call_per_site
        [INTRON4_BP_POSITION.min(cn_call_per_site.len())
            ..INTRON1_BP_POSITION.min(cn_call_per_site.len())]
        .iter()
        .filter_map(|&a| a)
        .collect();
    let mut intron4_intron1_sites_consensus =
        get_consensus(&intron4_intron1_sites, INTRON4_TO_INTRON1_SITES_MIN);

    // Fill in missing consensus from neighbor
    if exon9_intron4_sites_consensus.is_none() && intron4_intron1_sites_consensus.is_some() {
        exon9_intron4_sites_consensus = intron4_intron1_sites_consensus;
    } else if intron4_intron1_sites_consensus.is_none()
        && exon9_intron4_sites_consensus.is_some()
    {
        intron4_intron1_sites_consensus = exon9_intron4_sites_consensus;
    }

    // intron1 upstream sites
    let intron1_upstream_sites: Vec<u32> = cn_call_per_site
        [INTRON1_BP_POSITION.min(cn_call_per_site.len())..]
        .iter()
        .filter_map(|&a| a)
        .collect();
    let mut intron1_upstream_sites_consensus =
        get_consensus(&intron1_upstream_sites, INTRON1_UPSTREAM_SITES_MIN);

    if intron1_upstream_sites_consensus.is_none() && total_cn >= 2 {
        let expected = total_cn - 2;
        let count_expected = intron1_upstream_sites
            .iter()
            .filter(|&&v| v == expected)
            .count();
        if count_expected >= INTRON1_UPSTREAM_SITES_MIN_LOOSE {
            intron1_upstream_sites_consensus = Some(expected);
        }
    }

    // Exon9 region consensus
    // Python computes `total_cn - spacer_cn` as a signed int; an impossible
    // spacer_cn > total_cn yields a negative which then fails the downstream
    // intron1_upstream consistency check at line ~127. We mirror that by
    // treating underflow as no-call here (which propagates the same way).
    let mut exon9region_sites_consensus = if let Some(sp_cn) = spacer_cn {
        let mut consensus_val = total_cn.checked_sub(sp_cn);
        if let (Some(exon9gc), Some(exon9_intron4), Some(current)) = (
            exon9gc_call_stringent,
            exon9_intron4_sites_consensus,
            consensus_val,
        ) {
            if current < exon9gc && exon9gc <= exon9_intron4 {
                consensus_val = Some(exon9gc);
            } else if current > exon9gc && exon9gc >= exon9_intron4 {
                consensus_val = Some(exon9gc);
            }
        }
        consensus_val
    } else {
        let exon9region_sites: Vec<u32> = cn_call_per_site
            [REP_END_POSITION..EXON9_END_POSITION.min(cn_call_per_site.len())]
            .iter()
            .filter_map(|&a| a)
            .collect();
        get_consensus(&exon9region_sites, EXON9REGION_SITES_MIN)
    };

    if exon9region_sites_consensus.is_none() && exon9gc_call_stringent.is_some() {
        exon9region_sites_consensus = exon9gc_call_stringent;
    }

    let consensus = CnConsensus {
        rep: REP_SITES_CN,
        exon9_and_downstream: exon9region_sites_consensus,
        exon9_to_intron4: exon9_intron4_sites_consensus,
        intron4_to_intron1: intron4_intron1_sites_consensus,
        intron1_upstream: intron1_upstream_sites_consensus,
    };

    let cn_increase = ["dup", "exon9hyb", "star68"];
    let cn_decrease = ["star13intron1", "star13", "star5"];

    // Check intron1_upstream must equal total_cn - 2
    if total_cn < 2 || consensus.intron1_upstream.is_none() || consensus.intron1_upstream.unwrap() != total_cn - 2 {
        return (None, consensus);
    }

    let mut change_point: Vec<String> = Vec::new();

    // Assign CNV events based on CNs of the different regions
    if let (Some(intron4_intron1), Some(intron1_upstream)) =
        (consensus.intron4_to_intron1, consensus.intron1_upstream)
    {
        for _ in 0..intron4_intron1.saturating_sub(intron1_upstream) {
            change_point.push("star13intron1".to_string());
        }
        for _ in 0..intron1_upstream.saturating_sub(intron4_intron1) {
            change_point.push("star68".to_string());
        }
    }
    if let (Some(exon9_intron4), Some(intron4_intron1)) =
        (consensus.exon9_to_intron4, consensus.intron4_to_intron1)
    {
        for _ in 0..exon9_intron4.saturating_sub(intron4_intron1) {
            change_point.push("star13intron1".to_string());
        }
    }
    if let (Some(exon9_downstream), Some(exon9_intron4)) =
        (consensus.exon9_and_downstream, consensus.exon9_to_intron4)
    {
        for _ in 0..exon9_downstream.saturating_sub(exon9_intron4) {
            change_point.push("star13".to_string());
        }
        for _ in 0..exon9_intron4.saturating_sub(exon9_downstream) {
            change_point.push("exon9hyb".to_string());
        }
    }
    if let Some(exon9_downstream) = consensus.exon9_and_downstream {
        let rep = consensus.rep;
        for _ in 0..rep.saturating_sub(exon9_downstream) {
            change_point.push("star5".to_string());
        }
        for _ in 0..exon9_downstream.saturating_sub(rep) {
            change_point.push("dup".to_string());
        }
    }

    if check_cn_match(
        &change_point,
        &cn_increase,
        &cn_decrease,
        consensus.intron1_upstream.unwrap(),
    ) {
        change_point.sort();
        let joined = change_point.join("_");
        let sv_call = transform_cnvtag(&joined);
        if sv_call.is_empty() {
            // The exon9hyb/star5 and dup/star13 cancellation rules in
            // transform_cnvtag can cancel every detected event, leaving no
            // tag. Route that into the existing no-call path instead of
            // reporting an empty-string CNV_group (finding #20).
            log::info!(
                "get_cnvtag: CNV events cancelled out (change points \"{}\") -- no CNV tag",
                joined
            );
            return (None, consensus);
        }
        return (Some(sv_call), consensus);
    }

    // No CNV — check if it's cn2
    if consensus.exon9_to_intron4 == Some(2)
        && consensus.intron4_to_intron1 == Some(2)
        && consensus.intron1_upstream == Some(2)
    {
        return (Some("cn2".to_string()), consensus);
    }

    (None, consensus)
}

/// Rename some CNV tags for downstream processing.
fn transform_cnvtag(cnvtag: &str) -> String {
    let mut split_call: Vec<String> = cnvtag.split('_').map(|s| s.to_string()).collect();

    if cnvtag != "exon9hyb_star5" {
        while split_call.contains(&"exon9hyb".to_string())
            && split_call.contains(&"star5".to_string())
        {
            if let Some(pos) = split_call.iter().position(|s| s == "exon9hyb") {
                split_call.remove(pos);
            }
            if let Some(pos) = split_call.iter().position(|s| s == "star5") {
                split_call.remove(pos);
            }
        }
    }
    if cnvtag != "dup_star13" {
        while split_call.contains(&"dup".to_string())
            && split_call.contains(&"star13".to_string())
        {
            if let Some(pos) = split_call.iter().position(|s| s == "dup") {
                split_call.remove(pos);
            }
            if let Some(pos) = split_call.iter().position(|s| s == "star13") {
                split_call.remove(pos);
            }
        }
    }

    if split_call.iter().all(|s| s == "dup") && !split_call.is_empty() {
        return format!("cn{}", split_call.len() + 2);
    }
    if cnvtag == "dup_dup_exon9hyb_star13intron1" {
        return "cn4".to_string();
    }
    split_call.join("_")
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- transform_cnvtag: direct coverage of the pure helper ----
    //
    // finding #20: total_cn=4, all sites CN 2, exon9gc=None, spacer=Some(4)
    // drives get_cnvtag to build change_point events exon9hyb x2 + star5 x2,
    // which sort/join to exactly the string used below. transform_cnvtag
    // correctly cancels every matched pair down to "" -- that part of the
    // pure function is not the bug. The bug is one level up, in get_cnvtag,
    // which wrapped that "" in `Some("")` instead of treating a fully
    // cancelled result as a no-call. That is covered by the test below.

    #[test]
    fn transform_cnvtag_cancels_matched_pairs_to_empty_string() {
        assert_eq!(transform_cnvtag("exon9hyb_exon9hyb_star5_star5"), "");
    }

    /// Normal path: an uncancelled "dup" still becomes a real tag.
    #[test]
    fn transform_cnvtag_still_produces_a_real_tag_for_uncancelled_input() {
        assert_eq!(transform_cnvtag("dup"), "cn3");
    }

    /// Line ~210 special-cases this exact pair so it is NOT cancelled --
    /// it denotes a specific *5-with-exon9-hybrid call, not a cancellation.
    #[test]
    fn transform_cnvtag_does_not_cancel_the_exon9hyb_star5_special_case() {
        assert_eq!(transform_cnvtag("exon9hyb_star5"), "exon9hyb_star5");
    }

    // ---- get_cnvtag: the actual bug -- Some("") must become a no-call ----

    #[test]
    fn get_cnvtag_reports_no_call_instead_of_an_empty_tag_when_events_cancel() {
        // finding #20's exact traced scenario: total_cn=4, all sites CN 2,
        // exon9gc_call_stringent=None, spacer_cn=Some(4). This produces
        // change_point events exon9hyb x2 + star5 x2, which cancel to ""
        // in transform_cnvtag. Before the fix, get_cnvtag returned
        // `Some("".to_string())`, which main.rs would print as CNV_group "".
        let cn_call_per_site = vec![Some(2u32); 120];
        let raw_d6_cn = vec![2.0f64; 120];
        let (tag, _consensus) = get_cnvtag(4, &raw_d6_cn, &cn_call_per_site, None, Some(4));
        assert_eq!(
            tag, None,
            "an all-events-cancel result must be a no-call, not Some(\"\")"
        );
    }
}

/// Check that the CNV combination produces the right final copy number.
fn check_cn_match(
    sv_list: &[String],
    cn_increase: &[&str],
    cn_decrease: &[&str],
    final_cn: u32,
) -> bool {
    if sv_list.is_empty() {
        return false;
    }
    let mut initial_cn: i32 = 2;
    for sv in sv_list {
        if cn_increase.contains(&sv.as_str()) {
            initial_cn += 1;
        }
        if cn_decrease.contains(&sv.as_str()) {
            initial_cn -= 1;
        }
    }
    initial_cn == final_cn as i32
}

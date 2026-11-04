use crate::caller::fuzzy_match;
use crate::types::{FeatureFlags, RawStarCall, StarCall, StarCombinations};
use std::collections::HashMap;

/// Direct genotype calls for certain CNV tags.
fn cnvtag_to_genotype() -> HashMap<&'static str, &'static str> {
    let mut m = HashMap::new();
    m.insert("star5_star5", "*5/*5");
    m.insert("star13_star13", "*13/*13");
    m.insert("star13intron1_star13intron1", "*13/*13");
    m.insert("star5_star5_star68", "*5/*68");
    m
}

/// Sub-alleles reported under their own name instead of their main allele.
/// *4.013 is *4 with the exon 9 gene conversion; the exon9hyb rules in
/// get_final_call_clean build its tandem notation (e.g. "*1/*4.013+*4").
const KEPT_SUBALLELES: &[&str] = &["*4.013"];
const RARE_ALLELES: &[&str] = &["*34", "*39", "*4.009", "*139"];

/// Main-allele name used in calls: strips the sub-allele suffix ("*2.001" ->
/// "*2") unless the sub-allele is in KEPT_SUBALLELES.
pub(crate) fn main_allele_name(star: &str) -> &str {
    if KEPT_SUBALLELES.contains(&star) {
        star
    } else {
        star.split('.').next().unwrap_or(star)
    }
}

fn get_var_list(var_observed: &[String]) -> String {
    if var_observed.is_empty() {
        "NA".to_string()
    } else {
        var_observed.join("_")
    }
}

pub(crate) fn convert_to_main_allele(list_of_star: &[String]) -> Vec<String> {
    let mut converted_set = std::collections::HashSet::new();
    for stars in list_of_star {
        let converted: Vec<String> = stars
            .split('_')
            .map(|star| main_allele_name(star).to_string())
            .collect();
        let mut sorted = converted;
        sorted.sort();
        converted_set.insert(sorted.join("_"));
    }
    // Sorted so downstream picks (e.g. call_star68's last-tag fallback) are
    // deterministic rather than following HashSet iteration order.
    let mut out: Vec<String> = converted_set.into_iter().collect();
    out.sort();
    out
}

/// Copy counts of every variant that is not an unknown site.
fn known_site_counts<'a>(
    vars: impl Iterator<Item = &'a str>,
    unknown: &[String],
) -> std::collections::BTreeMap<&'a str, usize> {
    let mut counts = std::collections::BTreeMap::new();
    for v in vars {
        if !unknown.iter().any(|u| u == v) {
            *counts.entry(v).or_insert(0) += 1;
        }
    }
    counts
}

/// Keys whose variants agree with `var_observed` at every site not in
/// `unknown`; an unknown site may carry any number of copies. A full scan:
/// measured 0.80 s on the largest table (1,551,523 keys), 0.14 s or less on
/// the others. Only called when `unknown` is non-empty.
fn compatible_keys<'k>(
    var_observed: &[String],
    unknown: &[String],
    keys: impl Iterator<Item = &'k String>,
) -> Vec<&'k String> {
    let target = known_site_counts(var_observed.iter().map(String::as_str), unknown);
    keys.filter(|k| {
        let vars = if k.as_str() == "NA" { "" } else { k.as_str() };
        known_site_counts(vars.split('_').filter(|v| !v.is_empty()), unknown) == target
    })
    .collect()
}

/// get_star when some sites are unknown: every compatible candidate is kept
/// and no frequency rule narrows the set.
fn get_star_masked(var_observed: &[String], unknown: &[String], dic: &HashMap<String, Vec<String>>) -> RawStarCall {
    let mut raw: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    for k in compatible_keys(var_observed, unknown, dic.keys()) {
        raw.extend(dic[k].iter().cloned());
    }
    raw_to_call(raw.into_iter().collect())
}

fn get_star_cn1_masked(var_observed: &[String], unknown: &[String], dic: &HashMap<String, String>) -> RawStarCall {
    let raw: std::collections::BTreeSet<String> =
        compatible_keys(var_observed, unknown, dic.keys()).into_iter().map(|k| dic[k].clone()).collect();
    raw_to_call(raw.into_iter().collect())
}

/// Labels a masked candidate list the way get_star labels an exact one,
/// minus the rare-allele rule.
fn raw_to_call(raw: Vec<String>) -> RawStarCall {
    if raw.is_empty() {
        return RawStarCall { call_info: Some("no_match".to_string()), candidate: vec![], star_call: vec![] };
    }
    let processed = convert_to_main_allele(&raw);
    let info = if raw.len() == 1 {
        "unique_match"
    } else if processed.len() == 1 {
        "unique_star"
    } else {
        "more_than_one_match"
    };
    RawStarCall { call_info: Some(info.to_string()), candidate: raw, star_call: processed }
}

/// Return the star allele call and assign tags.
fn get_star(var_observed: &[String], unknown: &[String], dic: &HashMap<String, Vec<String>>) -> RawStarCall {
    if !unknown.is_empty() {
        return get_star_masked(var_observed, unknown, dic);
    }
    let mut sorted = var_observed.to_vec();
    sorted.sort();
    let var_list = get_var_list(&sorted);

    if !dic.contains_key(&var_list) {
        return RawStarCall {
            call_info: Some("no_match".to_string()),
            candidate: vec![],
            star_call: vec![],
        };
    }

    let matches = &dic[&var_list];

    // Check if it's a single star allele (CN=1 dhap returns String, not Vec)
    // In our implementation, dhap is separate. For Vec<String>:
    if matches.len() == 1 {
        let raw_stars = matches.clone();
        let processed = convert_to_main_allele(&raw_stars);
        return RawStarCall {
            call_info: Some("unique_match".to_string()),
            candidate: raw_stars,
            star_call: processed,
        };
    }

    // More than one match
    let raw_stars = matches.clone();
    let processed = convert_to_main_allele(&raw_stars);
    if processed.len() == 1 {
        return RawStarCall {
            call_info: Some("unique_star".to_string()),
            candidate: raw_stars,
            star_call: processed,
        };
    }

    // Filter rare alleles
    let non_rare: Vec<String> = raw_stars
        .iter()
        .filter(|hap| !RARE_ALLELES.iter().any(|rare| hap.contains(rare)))
        .cloned()
        .collect();
    // If every candidate is rare, filtering would empty the call; keep them all.
    if non_rare.is_empty() {
        return RawStarCall {
            call_info: Some("more_than_one_match".to_string()),
            candidate: raw_stars,
            star_call: processed,
        };
    }
    let processed_filtered = convert_to_main_allele(&non_rare);
    if processed_filtered.len() == 1 {
        return RawStarCall {
            call_info: Some("pick_common_allele".to_string()),
            candidate: raw_stars,
            star_call: processed_filtered,
        };
    }

    RawStarCall {
        call_info: Some("more_than_one_match".to_string()),
        candidate: raw_stars,
        star_call: processed_filtered,
    }
}

/// Get star for CN=1 (dhap returns a single String, not Vec).
fn get_star_cn1(var_observed: &[String], unknown: &[String], dic: &HashMap<String, String>) -> RawStarCall {
    if !unknown.is_empty() {
        return get_star_cn1_masked(var_observed, unknown, dic);
    }
    let mut sorted = var_observed.to_vec();
    sorted.sort();
    let var_list = get_var_list(&sorted);

    if !dic.contains_key(&var_list) {
        return RawStarCall {
            call_info: Some("no_match".to_string()),
            candidate: vec![],
            star_call: vec![],
        };
    }

    let star = dic[&var_list].clone();
    let processed = convert_to_main_allele(&[star.clone()]);
    RawStarCall {
        call_info: Some("unique_match".to_string()),
        candidate: vec![star],
        star_call: processed,
    }
}

/// Star allele call for star68-related CN=1 cases (dhap dictionary).
fn call_star68_cn1(
    var_observed: &mut Vec<String>,
    unknown: &[String],
    cnvcall: &str,
    dic: &HashMap<String, String>,
) -> RawStarCall {
    let mut matchtag = vec![get_star_cn1(var_observed, unknown, dic)];
    let cn = cnvcall.split('_').filter(|&s| s == "star68").count();
    let count0692 = var_observed
        .iter()
        .filter(|v| v.as_str() == "g.42130692G>A")
        .count();

    for _ in 0..std::cmp::min(cn, count0692) {
        if let Some(pos) = var_observed.iter().position(|v| v == "g.42130692G>A") {
            var_observed.remove(pos);
        }
        matchtag.push(get_star_cn1(var_observed, unknown, dic));
    }

    if matchtag.len() == 1 {
        return matchtag.remove(0);
    }

    let num_match = matchtag
        .iter()
        .filter(|t| t.call_info.as_deref() != Some("no_match"))
        .count();
    if num_match == 0 {
        return RawStarCall {
            call_info: Some("no_match".to_string()),
            candidate: vec![],
            star_call: vec![],
        };
    }

    let mut hap_list: Vec<String> = Vec::new();
    for tag in &matchtag {
        hap_list.extend(tag.star_call.iter().cloned());
    }
    let unique: std::collections::HashSet<&String> = hap_list.iter().collect();
    if unique.len() == 1 {
        return RawStarCall {
            call_info: Some("unique_match".to_string()),
            candidate: hap_list.clone(),
            star_call: vec![hap_list[0].clone()],
        };
    }

    if !unknown.is_empty() {
        // Unknown sites: report every candidate rather than picking by rule.
        let mut all: Vec<String> = unique.into_iter().cloned().collect();
        all.sort();
        return RawStarCall {
            call_info: Some("more_than_one_match".to_string()),
            candidate: hap_list.clone(),
            star_call: all,
        };
    }

    let non_rare: Vec<String> = hap_list
        .iter()
        .filter(|hap| !RARE_ALLELES.iter().any(|rare| hap.contains(rare)))
        .cloned()
        .collect();
    let unique_nr: std::collections::HashSet<&String> = non_rare.iter().collect();
    if unique_nr.len() == 1 {
        return RawStarCall {
            call_info: Some("pick_common_allele".to_string()),
            candidate: hap_list,
            star_call: vec![non_rare[0].clone()],
        };
    }

    RawStarCall {
        call_info: Some("pick_common_allele".to_string()),
        candidate: hap_list.clone(),
        star_call: vec![hap_list.last().unwrap().clone()],
    }
}

/// Star allele call for star68-related cases.
fn call_star68(
    var_observed: &mut Vec<String>,
    unknown: &[String],
    cnvcall: &str,
    dic: &HashMap<String, Vec<String>>,
) -> RawStarCall {
    let mut matchtag = vec![get_star(var_observed, unknown, dic)];
    let cn = cnvcall.split('_').filter(|&s| s == "star68").count();
    let count0692 = var_observed
        .iter()
        .filter(|v| v.as_str() == "g.42130692G>A")
        .count();

    for _ in 0..std::cmp::min(cn, count0692) {
        if let Some(pos) = var_observed.iter().position(|v| v == "g.42130692G>A") {
            var_observed.remove(pos);
        }
        matchtag.push(get_star(var_observed, unknown, dic));
    }

    if matchtag.len() == 1 {
        return matchtag.remove(0);
    }

    let num_match = matchtag
        .iter()
        .filter(|t| t.call_info.as_deref() != Some("no_match"))
        .count();
    if num_match == 0 {
        return RawStarCall {
            call_info: Some("no_match".to_string()),
            candidate: vec![],
            star_call: vec![],
        };
    }

    let mut hap_list: Vec<String> = Vec::new();
    for tag in &matchtag {
        hap_list.extend(tag.star_call.iter().cloned());
    }
    if hap_list.is_empty() {
        return RawStarCall {
            call_info: Some("no_match".to_string()),
            candidate: vec![],
            star_call: vec![],
        };
    }
    let unique: std::collections::HashSet<&String> = hap_list.iter().collect();
    if unique.len() == 1 {
        return RawStarCall {
            call_info: Some("unique_match".to_string()),
            candidate: hap_list.clone(),
            star_call: vec![hap_list[0].clone()],
        };
    }

    if !unknown.is_empty() {
        // Unknown sites: report every candidate rather than picking by rule.
        let mut all: Vec<String> = unique.into_iter().cloned().collect();
        all.sort();
        return RawStarCall {
            call_info: Some("more_than_one_match".to_string()),
            candidate: hap_list.clone(),
            star_call: all,
        };
    }

    let non_rare: Vec<String> = hap_list
        .iter()
        .filter(|hap| !RARE_ALLELES.iter().any(|rare| hap.contains(rare)))
        .cloned()
        .collect();
    let unique_nr: std::collections::HashSet<&String> = non_rare.iter().collect();
    if unique_nr.len() == 1 {
        return RawStarCall {
            call_info: Some("pick_common_allele".to_string()),
            candidate: hap_list,
            star_call: vec![non_rare[0].clone()],
        };
    }

    // The one with g.42130692G>A removed is the most likely
    RawStarCall {
        call_info: Some("pick_common_allele".to_string()),
        candidate: hap_list.clone(),
        star_call: vec![hap_list.last().unwrap().clone()],
    }
}

/// Return the variant table to use for each CNV/hybrid group. The tag -> table map
/// is `construct_star_table::combination_table_for`, which the builder also reads so
/// it can build just this one table (finding #41).
fn get_dic<'a>(cnvcall: &str, star_combinations: &'a StarCombinations) -> Option<&'a HashMap<String, Vec<String>>> {
    super::construct_star_table::combination_table_for(cnvcall).map(|kind| kind.table(star_combinations))
}

/// Check if cnvcall uses CN=1 (dhap).
fn uses_cn1_dic(cnvcall: &str) -> bool {
    matches!(
        cnvcall,
        "star5" | "star13" | "star5_star68" | "star13_star68"
    )
}

/// Clean up final call to report diplotypes.
pub fn get_final_call_clean(
    final_call: &[String],
    cnvcall: &str,
    spacer_cn: Option<u32>,
) -> Option<String> {
    let mut final_call = final_call.to_vec();
    final_call.sort();

    if final_call.len() >= 2 && cnvcall == "cn2" {
        let diplotypes: Vec<String> = final_call.iter().map(|c| c.replace('_', "/")).collect();
        return Some(diplotypes.join(";"));
    }

    if final_call.is_empty() {
        return Some(String::new());
    }
    if final_call.len() > 1 {
        if final_call == ["*10_*10_*4.013", "*10_*36_*4"] {
            // Two structurally-consistent decompositions with DIFFERENT allele
            // multisets ({*10,*36,*4} vs {*10,*10,*4.013}); emit both as proper
            // diplotype arrangements and let the likelihood resolver pick (or
            // report honest ambiguity). Was previously hardcoded to *4/*36+*10.
            return Some("*4/*36+*10;*4.013/*10+*10".to_string());
        }
        // Format each candidate the way a single answer is formatted. Joining
        // the raw multisets ("*1_*4_*4") gave strings match_phenotype cannot
        // read. A candidate that cannot be formatted keeps its raw form, which
        // reads as unknown effect: the safe direction.
        let mut out: Vec<String> = Vec::new();
        for one in &final_call {
            let formatted = get_final_call_clean(std::slice::from_ref(one), cnvcall, spacer_cn)
                .unwrap_or_else(|| one.clone());
            for d in formatted.split(';') {
                if !out.iter().any(|x| x == d) {
                    out.push(d.to_string());
                }
            }
        }
        return Some(out.join(";"));
    }

    let called_stars = &final_call[0];
    let split_call: Vec<&str> = called_stars.split('_').collect();

    if cnvcall == "star5_star68" {
        return if called_stars == "*4" {
            Some("*5/*68+*4".to_string())
        } else {
            Some(format!("*68/{}", called_stars))
        };
    }

    if cnvcall == "star13_star68" {
        return if called_stars == "*4" {
            Some("*13/*68+*4".to_string())
        } else {
            Some(format!("*13_{}_*68", called_stars))
        };
    }

    if cnvcall == "cn2" {
        return Some(split_call.join("/"));
    }

    if cnvcall == "star13intron1" {
        if split_call[0] == "*2" {
            return Some(format!("*13/{}", split_call[1]));
        }
        if split_call[1] == "*2" {
            return Some(format!("*13/{}", split_call[0]));
        }
        return None;
    }

    if cnvcall == "dup_star13intron1" {
        let unique: std::collections::HashSet<&&str> = split_call.iter().collect();
        let unique_vec: Vec<&&str> = unique.into_iter().collect();
        if unique_vec.len() == 2 {
            let u0 = unique_vec[0];
            let u1 = unique_vec[1];
            let count0 = split_call.iter().filter(|&&s| s == *u0).count();
            let count1 = split_call.iter().filter(|&&s| s == *u1).count();
            if count0 == 2 {
                return Some(format!("*13+{}/{}", u0, u1));
            } else if count1 == 2 {
                return Some(format!("*13+{}/{}", u1, u0));
            }
        } else if unique_vec.len() == 1 {
            return Some(format!("*13+{}/{}", unique_vec[0], unique_vec[0]));
        }
        return None;
    }

    if cnvcall == "dup_star13" {
        if spacer_cn == Some(1) {
            let mut sc = split_call.iter().map(|s| s.to_string()).collect::<Vec<_>>();
            sc.push("*13".to_string());
            return Some(sc.join("_"));
        }
        return Some(split_call.join("/"));
    }

    if cnvcall == "star5" || cnvcall == "star13" {
        let star_num = &cnvcall[4..];
        return Some(format!("{}/*{}", called_stars, star_num));
    }

    if cnvcall.starts_with("cn") && cnvcall.len() <= 3 {
        let mut dup_allele: Vec<String> = Vec::new();
        let mut call_set: Vec<String> = Vec::new();
        for &star_allele in &split_call {
            let count = split_call.iter().filter(|&&s| s == star_allele).count();
            if count > 1 {
                if !dup_allele.contains(&star_allele.to_string()) {
                    call_set.push(format!("{}x{}", star_allele, count));
                    dup_allele.push(star_allele.to_string());
                }
            } else {
                call_set.push(star_allele.to_string());
            }
        }
        if call_set.len() == 2 {
            return Some(call_set.join("/"));
        }
        if call_set.len() == 1 {
            let var = split_call[0];
            return match cnvcall {
                "cn3" => Some(format!("{}/{}x2", var, var)),
                "cn4" => Some(format!("{}x2/{}x2", var, var)),
                "cn5" => Some(format!("{}x2/{}x3", var, var)),
                "cn6" => Some(format!("{}x3/{}x3", var, var)),
                _ => Some(call_set.join("_")),
            };
        }
        return Some(call_set.join("_"));
    }

    // exon9hyb_star5: the tandem arrangement is determined by spacer copy
    // number (a depth signal). The likelihood resolver is blind to depth, so
    // this structural rule is kept here rather than delegated to the resolver.
    if cnvcall == "exon9hyb_star5" {
        if split_call == ["*10", "*10"] {
            return if spacer_cn.is_some_and(|cn| cn > 2) {
                Some("*5/*36+*10".to_string())
            } else {
                None
            };
        }
        if split_call == ["*4", "*4"] && spacer_cn.is_some_and(|cn| cn > 2) {
            return Some("*5/*4.013+*4".to_string());
        }
        return Some(split_call.join("/"));
    }

    if cnvcall == "exon9hyb" {
        if split_call.contains(&"*4") && split_call.contains(&"*4.013") {
            let idx_4 = split_call.iter().position(|&s| s == "*4").unwrap();
            let idx_4013 = split_call.iter().position(|&s| s == "*4.013").unwrap();
            let remain: Vec<usize> =
                (0..split_call.len()).filter(|&n| n != idx_4 && n != idx_4013).collect();
            if remain.len() != 1 {
                return None;
            }
            return Some(format!("{}/*4.013+*4", split_call[remain[0]]));
        }
        if split_call.contains(&"*10") && split_call.contains(&"*36") {
            let idx_10 = split_call.iter().position(|&s| s == "*10").unwrap();
            let idx_36 = split_call.iter().position(|&s| s == "*36").unwrap();
            let remain: Vec<usize> =
                (0..split_call.len()).filter(|&n| n != idx_10 && n != idx_36).collect();
            if remain.len() != 1 {
                return None;
            }
            let remaining = split_call[remain[0]];
            // When remaining is also *36, both tandem arrangements are
            // indistinguishable from short reads (*36 is a superset of *10).
            // Output both alternatives.
            if remaining == "*36" {
                return Some("*36/*36+*10;*10/*36x2".to_string());
            }
            return Some(format!("{}/*36+*10", remaining));
        }
        if split_call.iter().filter(|&&s| s == "*36").count() == 2 {
            let remain: Vec<&&str> = split_call.iter().filter(|&&s| s != "*36").collect();
            return Some(format!("*36x2/{}", remain[0]));
        }
    }

    if cnvcall == "dup_exon9hyb" {
        let mut var = Vec::new();
        if split_call.contains(&"*36") {
            var = split_call
                .iter()
                .filter(|&&s| s != "*10" && s != "*36")
                .cloned()
                .collect();
        } else if split_call.contains(&"*4.013") {
            var = split_call
                .iter()
                .filter(|&&s| s != "*4" && s != "*4.013")
                .cloned()
                .collect();
        }
        if var.len() == 2 {
            let mut remaining: Vec<String> = split_call.iter().map(|s| s.to_string()).collect();
            remaining.retain(|s| s != var[0] && s != var[1]);
            remaining.sort_by(|a, b| b.cmp(a));
            if var[0] == var[1] {
                return Some(format!("{}x2/{}", var[0], remaining.join("+")));
            } else {
                return Some(format!("{}+{}/{}", var[0], var[1], remaining.join("+")));
            }
        }
        if var.is_empty() {
            if *called_stars == "*10_*10_*10_*36" {
                return Some("*10x2/*36+*10".to_string());
            }
            if *called_stars == "*4_*4_*4_*4.013" {
                return Some("*4x2/*4.013+*4".to_string());
            }
        }
    }

    if cnvcall.contains("exon9hyb") && cnvcall != "exon9hyb" && !cnvcall.contains("star") && !cnvcall.contains("dup") {
        // Multiple exon9hyb
        // Special case: *10_*10_*36_*36 with elevated spacer CN indicates
        // a *5 deletion on one haplotype with tandem *36+*10 on the other.
        // Expected spacer CN: *36+*10/*36+*10 = 4, *5/*36x2+*10x2 = 5.
        if cnvcall == "exon9hyb_exon9hyb" && called_stars == "*10_*10_*36_*36" {
            return if spacer_cn.is_some_and(|cn| cn > 4) {
                Some("*5/*36x2+*10x2".to_string())
            } else {
                Some("*36+*10/*36+*10".to_string())
            };
        }
        let specific_cases: Vec<(&str, &str)> = vec![
            ("*4_*4_*4.013_*4.013", "*4.013+*4/*4.013+*4"),
            ("*10_*36_*36_*36", "*36+*10/*36x2"),
            ("*10_*10_*36_*36_*36", "*36+*10/*36x2+*10"),
            ("*10_*10_*36_*36_*36_*36", "*36x2+*10/*36x2+*10"),
        ];
        for (pattern, result) in &specific_cases {
            if called_stars.as_str() == *pattern {
                return Some(result.to_string());
            }
        }
        if cnvcall == "exon9hyb_exon9hyb_exon9hyb" {
            let mut sc = split_call.to_vec();
            if sc.contains(&"*10")
                && sc.contains(&"*83")
                && sc.iter().filter(|&&s| s == "*36").count() == 2
            {
                // Remove exactly: one *10, two *36, one *83 (matching Python's .remove())
                if let Some(pos) = sc.iter().position(|&s| s == "*10") { sc.remove(pos); }
                if let Some(pos) = sc.iter().position(|&s| s == "*36") { sc.remove(pos); }
                if let Some(pos) = sc.iter().position(|&s| s == "*36") { sc.remove(pos); }
                if let Some(pos) = sc.iter().position(|&s| s == "*83") { sc.remove(pos); }
                if sc.len() == 1 {
                    return Some(format!("{}/*36x2+*83+*10", sc[0]));
                }
            }
        }
        let var: Vec<&&str> = split_call
            .iter()
            .filter(|&&s| s != "*10" && s != "*36" && s != "*83")
            .collect();
        if var.len() == 1 {
            // Python: split_call.remove(var[0]) removes only first occurrence
            let mut remaining: Vec<String> = split_call.iter().map(|s| s.to_string()).collect();
            if let Some(pos) = remaining.iter().position(|s| s == *var[0]) {
                remaining.remove(pos);
            }
            remaining.sort_by(|a, b| b.cmp(a));
            return Some(format!("{}/{}", var[0], remaining.join("+")));
        }
    }

    if cnvcall.contains("star68") {
        let cn = cnvcall.split('_').filter(|&s| s == "star68").count();
        let parts: std::collections::HashSet<&str> = cnvcall.split('_').collect();

        if parts.len() == 1 {
            // All star68 (single distinct component)
            if split_call.contains(&"*4") {
                let var: Vec<&&str> = split_call.iter().filter(|&&s| s != "*4").collect();
                if var.len() == 1 {
                    let mut genotype = format!("{}/", var[0]);
                    for _ in 0..cn {
                        genotype.push_str("*68+");
                    }
                    genotype.push_str("*4");
                    return Some(genotype);
                } else if split_call == ["*4", "*4"] {
                    // Deterministic *68/*4 tandem-notation construction for the
                    // given star68 copy number — formatting, not a tie-break.
                    // The resolver is blind to *68 hybrid arrangement, so this
                    // canonical notation stays here.
                    return match cn {
                        1 => Some("*4/*68+*4".to_string()),
                        2 => Some("*68+*4/*68+*4".to_string()),
                        3 => Some("*68+*4/*68+*68+*4".to_string()),
                        4 => Some("*68+*68+*4/*68+*68+*4".to_string()),
                        _ => None,
                    };
                }
            }
            if split_call.contains(&"*45") {
                let var: Vec<&&str> = split_call.iter().filter(|&&s| s != "*45").collect();
                if var.len() == 1 {
                    let mut genotype = format!("{}/", var[0]);
                    for _ in 0..cn {
                        genotype.push_str("*68+");
                    }
                    genotype.push_str("*45");
                    return Some(genotype);
                }
            }
        }

        if cnvcall == "dup_star68" {
            let var: Vec<&&str> = split_call
                .iter()
                .filter(|&&s| s != "*4" && s != "*68")
                .collect();
            if var.len() == 2 && var[0] == var[1] {
                return Some(format!("{}x2/*68+*4", var[0]));
            }
            if var.is_empty() && *called_stars == "*4_*4_*4" {
                return Some("*4x2/*68+*4".to_string());
            }
        }

        if cnvcall == "exon9hyb_star68" && *called_stars == "*4_*4_*4.013" {
            return Some("*4.013+*4/*68+*4".to_string());
        }

        // Default: add *68 to the call
        let mut sc: Vec<String> = split_call.iter().map(|s| s.to_string()).collect();
        for _ in 0..cn {
            sc.push("*68".to_string());
        }
        return Some(sc.join("_"));
    }

    Some(called_stars.clone())
}

/// Update variants based on called CNV.
pub fn update_variants(
    var_observed: &mut Vec<String>,
    cnvcall: &str,
    exon9: &crate::types::Exon9Values,
) {
    // g.42129809T>C and g.42129819G>T should have the same CN
    if var_observed.contains(&"g.42129809T>C".to_string()) {
        let count_809 = var_observed
            .iter()
            .filter(|v| v.as_str() == "g.42129809T>C")
            .count();
        let count_819 = var_observed
            .iter()
            .filter(|v| v.as_str() == "g.42129819G>T")
            .count();
        for _ in 0..(count_809.saturating_sub(count_819)) {
            var_observed.push("g.42129819G>T".to_string());
        }
    }

    // g.42127556T>C included in g.42127565T>C definition for *108
    if var_observed.contains(&"g.42127565T>C".to_string()) {
        let count_565 = var_observed
            .iter()
            .filter(|v| v.as_str() == "g.42127565T>C")
            .count();
        let count_556 = var_observed
            .iter()
            .filter(|v| v.as_str() == "g.42127556T>C")
            .count();
        for _ in 0..(count_565.saturating_sub(count_556)) {
            var_observed.push("g.42127556T>C".to_string());
        }
    }

    // g.42126611C>G in D6 part of hybrid gene
    if cnvcall.contains("star13") && !cnvcall.contains("intron1") {
        if let Some(pos) = var_observed.iter().position(|v| v == "g.42126611C>G") {
            var_observed.remove(pos);
        }
    }

    // Exon9hyb handling
    if cnvcall.contains("exon9hyb") && cnvcall != "exon9hyb_star5" {
        let cn = cnvcall.split('_').filter(|&s| s == "exon9hyb").count();
        if !var_observed.contains(&"exon9gc".to_string()) {
            for _ in 0..cn {
                var_observed.push("exon9gc".to_string());
                var_observed.push("g.42126611C>G".to_string());
            }
        }
        // Add variants if not called to sufficient CN
        {
            let var_to_add = &"g.42130692G>A";
            let count = var_observed.iter().filter(|v| v.as_str() == *var_to_add).count();
            if count > 0 && count <= cn {
                var_observed.push(var_to_add.to_string());
            }
        }
        {
            let var_to_add = &"g.42129754G>A";
            let count = var_observed.iter().filter(|v| v.as_str() == *var_to_add).count();
            if count > 0
                && count <= cn
                && !var_observed.contains(&"g.42128945C>T".to_string())
            {
                var_observed.push(var_to_add.to_string());
            }
        }
        for var_to_add in &["g.42128945C>T", "g.42129809T>C", "g.42129819G>T"] {
            let count = var_observed.iter().filter(|v| v.as_str() == *var_to_add).count();
            if count > 0
                && count <= cn
                && !var_observed.contains(&"g.42129754G>A".to_string())
            {
                var_observed.push(var_to_add.to_string());
            }
        }
    }

    // Exon 9 gene conversion by itself
    if cnvcall.contains("exon9hyb") || cnvcall == "cn2" {
        if let (Some(exon9cn_consensus), Some(exon9_cn)) =
            (exon9.exon9cn_in_consensus, exon9.exon9_cn)
        {
            if exon9cn_consensus > exon9_cn && exon9_cn <= 1
                && (!var_observed.contains(&"g.42126611C>G".to_string())
                    || (exon9.exon9_raw_site1.min(exon9.exon9_raw_site2) < 1.15
                        && exon9.exon9_raw_site1.max(exon9.exon9_raw_site2) < 1.2))
            {
                for _ in 0..(exon9cn_consensus - exon9_cn) {
                    var_observed.push("exon9gc".to_string());
                }
                if !var_observed.contains(&"g.42126611C>G".to_string()) {
                    for _ in 0..(exon9cn_consensus - exon9_cn) {
                        var_observed.push("g.42126611C>G".to_string());
                    }
                }
            }
        }
    }

    var_observed.retain(|v| v != "NA");
}

/// Targeted read-phase rules for two cn2 ambiguities. Each rule fires when
/// both of its diplotypes are among the candidates and drops the one its
/// haplotype signal excludes; any other candidates are kept.
///
/// Each signal is three-state, not a boolean. `Some(true)` and `Some(false)`
/// are positive evidence for one arrangement; `None` means the phase step
/// did not decide (too few linked reads, no evidence either way, or it never
/// ran). On `None` the rule is skipped and neither diplotype is dropped, so
/// the call stays ambiguous rather than being resolved by an absence of
/// evidence.
fn apply_cn2_phase_rules(
    clean_call: &str,
    var42126938_g_haplotype: Option<bool>,
    var42127803_diff_haplotype: Option<bool>,
) -> String {
    let mut diplotypes: Vec<&str> = clean_call.split(';').collect();
    let rules = [
        ("*1/*32", "*27/*41", var42126938_g_haplotype),
        ("*119/*2", "*1/*41", var42127803_diff_haplotype),
    ];
    for (if_true, if_false, signal) in rules {
        let Some(signal) = signal else { continue };
        if diplotypes.contains(&if_true) && diplotypes.contains(&if_false) {
            let excluded = if signal { if_false } else { if_true };
            diplotypes.retain(|d| *d != excluded);
        }
    }
    diplotypes.join(";")
}

/// For a cnN (N >= 3) no_match, try adding one more copy of each variant
/// currently seen N-2 times. Returns the result if exactly one distinct
/// main-allele call comes out of all retries that matched unambiguously
/// (unique_match, or unique_star: several sub-allele candidates with the
/// same main alleles).
fn retry_add_one_copy(
    var_observed: &[String],
    cnvcall: &str,
    dic: &HashMap<String, Vec<String>>,
) -> Option<RawStarCall> {
    let cn_num: usize = cnvcall[2..].parse().unwrap_or(0);
    // cn_num is always >= 2 for valid cnN tags; guard the subtraction.
    let target_count = cn_num.saturating_sub(2);
    let mut variant_tried: Vec<&String> = Vec::new();
    let mut matched_calls: Vec<RawStarCall> = Vec::new();
    for variant_to_try in var_observed {
        let count = var_observed.iter().filter(|v| *v == variant_to_try).count();
        if count == target_count && !variant_tried.contains(&variant_to_try) {
            variant_tried.push(variant_to_try);
            let mut var_new = var_observed.to_vec();
            var_new.push(variant_to_try.clone());
            let matchtag_new = get_star(&var_new, &[], dic);
            let unambiguous = matches!(
                matchtag_new.call_info.as_deref(),
                Some("unique_match") | Some("unique_star")
            );
            if unambiguous
                && !matched_calls.iter().any(|m| m.star_call == matchtag_new.star_call)
            {
                matched_calls.push(matchtag_new);
            }
        }
    }
    if matched_calls.len() == 1 {
        matched_calls.pop()
    } else {
        None
    }
}

/// Main star allele matching function.
pub fn match_star(
    var_observed: &mut Vec<String>,
    cnvcall: &str,
    spacer_cn: Option<u32>,
    star_combinations: &StarCombinations,
    exon9: &crate::types::Exon9Values,
    var42126938_g_haplotype: Option<bool>,
    var42127803_diff_haplotype: Option<bool>,
    features: &FeatureFlags,
) -> StarCall {
    let direct = cnvtag_to_genotype();
    if let Some(&genotype) = direct.get(cnvcall) {
        return StarCall {
            call_info: Some("unique_match".to_string()),
            variants_called: Some(String::new()),
            raw_call: Some(vec![genotype.to_string()]),
            clean_call: Some(genotype.to_string()),
        };
    }

    update_variants(var_observed, cnvcall, exon9);

    if uses_cn1_dic(cnvcall) {
        if cnvcall.contains("star68") {
            // star5_star68, star13_star68: use call_star68_cn1 for iterative matching
            let mut matchtag = call_star68_cn1(var_observed, &[], cnvcall, &star_combinations.dhap);
            // Feature: fuzzy_match — try fuzzy matching on no_match for star68 CN=1
            if features.fuzzy_match && matchtag.call_info.as_deref() == Some("no_match") {
                if let Some(fuzzy) = fuzzy_match::fuzzy_match_star_cn1(var_observed, &star_combinations.dhap) {
                    matchtag = fuzzy;
                }
            }
            let final_call = &matchtag.star_call;
            let final_call_clean = get_final_call_clean(final_call, cnvcall, spacer_cn);
            return StarCall {
                call_info: matchtag.call_info,
                variants_called: Some(var_observed.join(" ")),
                raw_call: Some(matchtag.candidate),
                clean_call: final_call_clean,
            };
        }
        let mut matchtag = get_star_cn1(var_observed, &[], &star_combinations.dhap);
        // Feature: fuzzy_match — try fuzzy matching on no_match
        if features.fuzzy_match && matchtag.call_info.as_deref() == Some("no_match") {
            if let Some(fuzzy) = fuzzy_match::fuzzy_match_star_cn1(var_observed, &star_combinations.dhap) {
                matchtag = fuzzy;
            }
        }
        let final_call = &matchtag.star_call;
        let final_call_clean = get_final_call_clean(final_call, cnvcall, spacer_cn);
        return StarCall {
            call_info: matchtag.call_info,
            variants_called: Some(var_observed.join(" ")),
            raw_call: Some(matchtag.candidate),
            clean_call: final_call_clean,
        };
    }

    let dic = match get_dic(cnvcall, star_combinations) {
        Some(d) => d,
        None => {
            return StarCall {
                call_info: None,
                variants_called: None,
                raw_call: None,
                clean_call: None,
            }
        }
    };

    // Every one of the ten combination tables is non-empty once built (sizes 1 to
    // 1,551,523), so an empty one here is never the data -- it means this
    // `StarCombinations` was built without the `TableKind` this tag selects. Since
    // finding #41 that is a live hazard: `d6_star_caller` holds a locally built
    // `StarCombinations` with one table in it while `params.star_combinations`, with
    // all ten empty, is still in scope, and nothing in the type distinguishes them.
    // Reaching for the wrong one would leave `get_star` with no matching key, so it
    // would return `no_match` -- no panic, no log, one silently changed call.
    // Checked only in debug builds: `cargo test` without `--release` and debug runs
    // catch the mistake; release builds, `cargo test --release` included, compile
    // it out, and production does not panic on it.
    debug_assert!(
        !dic.is_empty(),
        "combination table for CNV tag {} is empty: this StarCombinations was built \
         without TableKind {:?}",
        cnvcall,
        super::construct_star_table::combination_table_for(cnvcall)
    );

    if !cnvcall.contains("star68") {
        let mut matchtag = get_star(var_observed, &[], dic);

        // For cn3-cn6, try adding a variant on no_match
        if cnvcall.starts_with("cn")
            && cnvcall.len() <= 3
            && matchtag.call_info.as_deref() == Some("no_match")
        {
            if let Some(mt) = retry_add_one_copy(var_observed, cnvcall, dic) {
                let final_call_clean =
                    get_final_call_clean(&mt.star_call, cnvcall, spacer_cn);
                return StarCall {
                    call_info: mt.call_info,
                    variants_called: Some(var_observed.join(" ")),
                    raw_call: Some(mt.candidate),
                    clean_call: final_call_clean,
                };
            }
        }

        // Feature: fuzzy_match — try fuzzy matching on no_match (after cn3+ retry)
        if features.fuzzy_match && matchtag.call_info.as_deref() == Some("no_match") {
            if let Some(fuzzy) = fuzzy_match::fuzzy_match_star(var_observed, dic) {
                matchtag = fuzzy;
            }
        }

        let final_call = &matchtag.star_call;
        let mut final_call_clean = get_final_call_clean(final_call, cnvcall, spacer_cn);
        let call_info = matchtag.call_info.clone();

        // Resolve ambiguous cn2 cases with haplotype info
        // These two cn2 ambiguities are resolved by targeted read-phase
        // signals (var42126938 / var42127803 three-state haplotype signals,
        // where `None` means undecided and drops nothing), which are
        // more specific than the general likelihood resolver for these exact
        // allele sets. Kept here pending benchmark coverage that could
        // validate delegating them to the resolver.
        if call_info.as_deref() == Some("more_than_one_match") && cnvcall == "cn2" {
            if let Some(fcc) = final_call_clean.as_deref() {
                final_call_clean = Some(apply_cn2_phase_rules(
                    fcc,
                    var42126938_g_haplotype,
                    var42127803_diff_haplotype,
                ));
            }
        }

        StarCall {
            call_info,
            variants_called: Some(var_observed.join(" ")),
            raw_call: Some(matchtag.candidate),
            clean_call: final_call_clean,
        }
    } else {
        let mut var_observed_68 = var_observed.clone();
        let matchtag = call_star68(&mut var_observed_68, &[], cnvcall, dic);
        let final_call_clean =
            get_final_call_clean(&matchtag.star_call, cnvcall, spacer_cn);
        StarCall {
            call_info: matchtag.call_info,
            variants_called: Some(var_observed_68.join(" ")),
            raw_call: Some(matchtag.candidate),
            clean_call: final_call_clean,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unknown_phase_keeps_both_diplotypes() {
        assert_eq!(apply_cn2_phase_rules("*1/*41;*119/*2", None, None), "*1/*41;*119/*2");
        assert_eq!(apply_cn2_phase_rules("*1/*41;*119/*2", None, Some(true)), "*119/*2");
        assert_eq!(apply_cn2_phase_rules("*1/*41;*119/*2", None, Some(false)), "*1/*41");
        assert_eq!(apply_cn2_phase_rules("*1/*32;*27/*41", Some(true), None), "*1/*32");
        assert_eq!(apply_cn2_phase_rules("*1/*32;*27/*41", Some(false), None), "*27/*41");
        assert_eq!(apply_cn2_phase_rules("*1/*32;*27/*41", None, None), "*1/*32;*27/*41");
    }

    fn dic(entries: &[(&str, &[&str])]) -> HashMap<String, Vec<String>> {
        entries
            .iter()
            .map(|(k, v)| (k.to_string(), v.iter().map(|s| s.to_string()).collect()))
            .collect()
    }

    fn obs(v: &[&str]) -> Vec<String> { v.iter().map(|s| s.to_string()).collect() }

    #[test]
    fn an_unknown_site_never_acts_as_a_reference_call() {
        // Today "B unseen" means "no B": only *x. Masked, *y (which has B) stays.
        let d = dic(&[("g.A", &["*x"]), ("g.A_g.B", &["*y"])]);
        let call = get_star(&obs(&["g.A"]), &obs(&["g.B"]), &d);
        assert_eq!(call.candidate, obs(&["*x", "*y"]));
        assert_eq!(get_star(&obs(&["g.A"]), &[], &d).candidate, obs(&["*x"]));
    }

    #[test]
    fn unknown_sites_are_masked_together() {
        let d = dic(&[("g.A", &["*x"]), ("g.A_g.B_g.C", &["*y"]), ("g.A_g.C", &["*z"])]);
        let call = get_star(&obs(&["g.A"]), &obs(&["g.B", "g.C"]), &d);
        assert_eq!(call.candidate, obs(&["*x", "*y", "*z"]));
    }

    #[test]
    fn a_known_site_still_excludes_candidates() {
        // Masking B must not admit *w, which lacks the observed A.
        let d = dic(&[("g.A", &["*x"]), ("g.B", &["*w"])]);
        assert_eq!(get_star(&obs(&["g.A"]), &obs(&["g.B"]), &d).candidate, obs(&["*x"]));
    }

    #[test]
    fn a_frequency_rule_does_not_narrow_a_masked_set() {
        // *34 is in RARE_ALLELES. Unmasked, the rare rule picks the common one.
        let d = dic(&[("g.A", &["*1_*2"]), ("g.A_g.B", &["*1_*34"])]);
        let masked = get_star(&obs(&["g.A"]), &obs(&["g.B"]), &d);
        assert_eq!(masked.star_call, obs(&["*1_*2", "*1_*34"]));
        assert_eq!(masked.call_info.as_deref(), Some("more_than_one_match"));
    }

    #[test]
    fn an_irrelevant_unknown_site_gives_the_exact_answer_on_the_real_table() {
        // Scan path vs exact path, every 100th key of the real cn2 table.
        let c = crate::caller::construct_star_table::get_hap_table(
            crate::data::STAR_TABLE, &[crate::caller::construct_star_table::TableKind::Dhap2]);
        let t = &c.dhap2;
        let mut keys: Vec<&String> = t.keys().collect();
        keys.sort();
        let none = obs(&["g.not_a_site"]);
        let mut checked = 0;
        for k in keys.iter().step_by(100) {
            let vars: Vec<String> = if *k == "NA" { vec![] } else { k.split('_').map(String::from).collect() };
            let mut exact = t[*k].clone();
            exact.sort();
            exact.dedup(); // the masked path returns a set
            assert_eq!(get_star(&vars, &none, t).candidate, exact, "key {k}");
            checked += 1;
        }
        assert!(checked >= 200, "only {checked} keys checked");
    }

    #[test]
    fn damaged_na19920_keeps_its_true_genotype_among_the_candidates() {
        // Measured in gate B1: 5 unknown sites, 9 candidate multisets.
        let c = crate::caller::construct_star_table::get_hap_table(
            crate::data::STAR_TABLE, &[crate::caller::construct_star_table::TableKind::Dhap3Pair]);
        let observed = obs(&["g.42126611C>G", "g.42126611C>G", "g.42129130C>G", "g.42129130C>G",
                             "g.42129754G>A", "g.42129754G>A", "g.42130692G>A", "g.42130692G>A"]);
        let unknown = obs(&["g.42128936-42128937insGGGGCGAAAGGGGCGAAA", "g.42128936-42128937insGGGGCGAAA",
                            "g.42128932C>T", "g.42128933G>A", "g.42128945C>T"]);
        let call = get_star(&observed, &unknown, &c.dhap3pair);
        assert_eq!(call.star_call, obs(&["*10_*10_*145", "*10_*10_*4", "*145_*147_*147", "*145_*4_*4",
            "*147_*147_*4", "*1_*10_*10", "*1_*147_*147", "*1_*4_*4", "*4_*4_*4"]));
        // Unmasked, today's answer: the single wrong *1_*10_*10.
        assert_eq!(get_star(&observed, &[], &c.dhap3pair).star_call, obs(&["*1_*10_*10"]));
    }

    #[test]
    fn get_star_keeps_candidates_when_all_are_rare() {
        // Every candidate carries a rare allele and they differ after main-allele
        // conversion: the rare filter must not empty the call.
        let d = dic(&[("g.X", &["*34_*1", "*39_*1"])]);
        let out = get_star(&["g.X".to_string()], &[], &d);
        assert_eq!(out.call_info.as_deref(), Some("more_than_one_match"));
        assert_eq!(out.star_call, vec!["*1_*34".to_string(), "*1_*39".to_string()]);
    }

    #[test]
    fn star_4_013_keeps_its_name() {
        let d = dic(&[("g.X", &["*1_*4.013"])]);
        let out = get_star(&["g.X".to_string()], &[], &d);
        assert_eq!(out.star_call, vec!["*1_*4.013".to_string()]);
    }

    #[test]
    fn cn2_with_three_candidates_reports_diplotypes() {
        let calls: Vec<String> =
            ["*1_*32", "*119_*162", "*27_*41"].iter().map(|s| s.to_string()).collect();
        assert_eq!(
            get_final_call_clean(&calls, "cn2", None).as_deref(),
            Some("*119/*162;*1/*32;*27/*41")
        );
    }

    #[test]
    fn cn2_phase_rules_drop_only_the_excluded_diplotype() {
        assert_eq!(
            apply_cn2_phase_rules("*119/*162;*1/*32;*27/*41", Some(true), Some(false)),
            "*119/*162;*1/*32"
        );
        assert_eq!(apply_cn2_phase_rules("*1/*32;*27/*41", Some(false), Some(false)), "*27/*41");
        assert_eq!(apply_cn2_phase_rules("*119/*2;*1/*41", Some(false), Some(true)), "*119/*2");
        assert_eq!(apply_cn2_phase_rules("*1/*2;*1/*3", Some(true), Some(true)), "*1/*2;*1/*3");
    }

    #[test]
    fn cn3_retry_accepts_a_single_main_allele_result() {
        // Adding a second copy of g.A gives two raw candidates with the same
        // main alleles (unique_star); that is an unambiguous answer.
        let d = dic(&[("g.A_g.A_g.B", &["*1_*2.001_*2", "*1_*2_*2.002"])]);
        let out = retry_add_one_copy(&["g.A".to_string(), "g.B".to_string()], "cn3", &d).unwrap();
        assert_eq!(out.star_call, vec!["*1_*2_*2".to_string()]);
    }

    #[test]
    fn the_star68_paths_do_not_pick_by_habit_when_sites_are_unknown() {
        // Two different answers after removing one g.42130692G>A copy. Unmasked,
        // the "last one is most likely" rule picks one; masked, both are kept.
        let d = dic(&[("g.42130692G>A_g.X", &["*68_*4"]), ("g.X", &["*68_*2"])]);
        let mut v = obs(&["g.42130692G>A", "g.X"]);
        let picked = call_star68(&mut v.clone(), &[], "star68", &d);
        assert_eq!(picked.star_call.len(), 1);
        let kept = call_star68(&mut v, &obs(&["g.Y"]), "star68", &d);
        assert_eq!(kept.star_call.len(), 2, "{:?}", kept.star_call);
        assert_eq!(kept.call_info.as_deref(), Some("more_than_one_match"));
    }

    #[test]
    fn the_cn1_star68_path_does_not_pick_by_habit_when_sites_are_unknown() {
        // Single-copy table: values are one allele each.
        let d: HashMap<String, String> = [("g.42130692G>A_g.X", "*68"), ("g.X", "*2")]
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect();
        let v = obs(&["g.42130692G>A", "g.X"]);
        let picked = call_star68_cn1(&mut v.clone(), &[], "star5_star68", &d);
        assert_eq!(picked.star_call.len(), 1);
        let kept = call_star68_cn1(&mut v.clone(), &obs(&["g.Y"]), "star5_star68", &d);
        assert_eq!(kept.star_call, obs(&["*2", "*68"]));
        assert_eq!(kept.call_info.as_deref(), Some("more_than_one_match"));
    }

    #[test]
    fn call_star68_does_not_panic_when_matching_tag_is_all_rare() {
        // Tag 0 (with g.42130692G>A) is no_match; tag 1 (after removing it)
        // matches only rare-allele candidates.
        let d = dic(&[("g.X", &["*34_*1", "*39_*1"])]);
        let mut vars = vec!["g.42130692G>A".to_string(), "g.X".to_string()];
        let out = call_star68(&mut vars, &[], "star68", &d);
        assert!(!out.star_call.is_empty());
    }

    #[test]
    fn several_candidates_outside_cn2_are_each_formatted() {
        // HG02373 at 25 % depth (spacer CN 3). At full depth its single
        // candidate *10_*14_*36 is formatted *14/*36+*10 (measured).
        let call = vec![
            "*10_*14_*36".to_string(),
            "*10_*36_*65".to_string(),
            "*10_*114_*36".to_string(),
        ];
        // Each candidate formatted alone (measured): *14/*36+*10, *65/*36+*10,
        // *114/*36+*10. Today all three together give the raw
        // "*10_*114_*36;*10_*14_*36;*10_*36_*65" (measured).
        assert_eq!(
            get_final_call_clean(&call, "exon9hyb", Some(3)).as_deref(),
            Some("*114/*36+*10;*14/*36+*10;*65/*36+*10")
        );
    }

    #[test]
    fn the_hand_written_two_candidate_case_is_unchanged() {
        let call = vec!["*10_*10_*4.013".to_string(), "*10_*36_*4".to_string()];
        assert_eq!(
            get_final_call_clean(&call, "exon9hyb", Some(3)).as_deref(),
            Some("*4/*36+*10;*4.013/*10+*10")
        );
    }
}

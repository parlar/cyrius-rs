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

const KEPT_SUBALLELES: &[&str] = &[];
const RARE_ALLELES: &[&str] = &["*34", "*39", "*4.009", "*139"];

fn get_var_list(var_observed: &[String]) -> String {
    if var_observed.is_empty() {
        "NA".to_string()
    } else {
        var_observed.join("_")
    }
}

fn convert_to_main_allele(list_of_star: &[String]) -> Vec<String> {
    let mut converted_set = std::collections::HashSet::new();
    for stars in list_of_star {
        let star_split: Vec<&str> = stars.split('_').collect();
        let converted: Vec<String> = star_split
            .iter()
            .map(|star| {
                if !KEPT_SUBALLELES.contains(star) {
                    star.split('.').next().unwrap().to_string()
                } else {
                    star.to_string()
                }
            })
            .collect();
        let mut sorted = converted;
        sorted.sort();
        converted_set.insert(sorted.join("_"));
    }
    converted_set.into_iter().collect()
}

/// Return the star allele call and assign tags.
fn get_star(var_observed: &[String], dic: &HashMap<String, Vec<String>>) -> RawStarCall {
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
fn get_star_cn1(var_observed: &[String], dic: &HashMap<String, String>) -> RawStarCall {
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
    cnvcall: &str,
    dic: &HashMap<String, String>,
) -> RawStarCall {
    let mut matchtag = vec![get_star_cn1(var_observed, dic)];
    let cn = cnvcall.split('_').filter(|&s| s == "star68").count();
    let count0692 = var_observed
        .iter()
        .filter(|v| v.as_str() == "g.42130692G>A")
        .count();

    for _ in 0..std::cmp::min(cn, count0692) {
        if let Some(pos) = var_observed.iter().position(|v| v == "g.42130692G>A") {
            var_observed.remove(pos);
        }
        matchtag.push(get_star_cn1(var_observed, dic));
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
    cnvcall: &str,
    dic: &HashMap<String, Vec<String>>,
) -> RawStarCall {
    let mut matchtag = vec![get_star(var_observed, dic)];
    let cn = cnvcall.split('_').filter(|&s| s == "star68").count();
    let count0692 = var_observed
        .iter()
        .filter(|v| v.as_str() == "g.42130692G>A")
        .count();

    for _ in 0..std::cmp::min(cn, count0692) {
        if let Some(pos) = var_observed.iter().position(|v| v == "g.42130692G>A") {
            var_observed.remove(pos);
        }
        matchtag.push(get_star(var_observed, dic));
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

/// Return the variant table to use for each CNV/hybrid group.
fn get_dic<'a>(cnvcall: &str, star_combinations: &'a StarCombinations) -> Option<&'a HashMap<String, Vec<String>>> {
    match cnvcall {
        "star5" | "star13" | "star5_star68" | "star13_star68" => None, // Uses dhap (CN=1)
        "cn2" | "star13intron1" | "dup_star13" | "exon9hyb_star5" => {
            Some(&star_combinations.dhap2)
        }
        s if s.starts_with("star68") && !s.contains("exon9hyb") && !s.contains("dup") => {
            Some(&star_combinations.dhap2)
        }
        "exon9hyb" | "exon9hyb_star68" => Some(&star_combinations.dhap3),
        "cn3" | "dup_star13intron1" | "dup_star68" => Some(&star_combinations.dhap3pair),
        "cn4" => Some(&star_combinations.dhap4pair),
        "cn5" => Some(&star_combinations.dhap5pair),
        "cn6" => Some(&star_combinations.dhap6pair),
        "exon9hyb_exon9hyb" => Some(&star_combinations.dhap_exon9_x2),
        "exon9hyb_exon9hyb_exon9hyb" => Some(&star_combinations.dhap_exon9_x3),
        "exon9hyb_exon9hyb_exon9hyb_exon9hyb" => Some(&star_combinations.dhap_exon9_x4),
        "dup_exon9hyb" => Some(&star_combinations.dhap_dup_exon9),
        _ => None,
    }
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

    if final_call.len() == 2 && cnvcall == "cn2" {
        let d1: Vec<&str> = final_call[0].split('_').collect();
        let d2: Vec<&str> = final_call[1].split('_').collect();
        return Some(format!("{};{}", d1.join("/"), d2.join("/")));
    }

    if final_call.is_empty() || final_call.len() > 1 {
        if final_call == ["*10_*10_*4.013", "*10_*36_*4"] {
            return Some("*4/*36+*10".to_string());
        }
        return Some(final_call.join(";"));
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

    if cnvcall == "exon9hyb_star5" {
        if split_call == ["*10", "*10"] {
            return if spacer_cn.map_or(false, |cn| cn > 2) {
                Some("*5/*36+*10".to_string())
            } else {
                None
            };
        }
        if split_call == ["*4", "*4"] {
            if spacer_cn.map_or(false, |cn| cn > 2) {
                return Some("*5/*4.013+*4".to_string());
            }
        }
        return Some(split_call.join("/"));
    }

    if cnvcall == "exon9hyb" {
        if split_call.contains(&"*4") && split_call.contains(&"*4.013") {
            let idx_4 = split_call.iter().position(|&s| s == "*4").unwrap();
            let idx_4013 = split_call.iter().position(|&s| s == "*4.013").unwrap();
            let remain: Vec<usize> = (0..3).filter(|&n| n != idx_4 && n != idx_4013).collect();
            assert_eq!(remain.len(), 1);
            return Some(format!("{}/*4.013+*4", split_call[remain[0]]));
        }
        if split_call.contains(&"*10") && split_call.contains(&"*36") {
            let idx_10 = split_call.iter().position(|&s| s == "*10").unwrap();
            let idx_36 = split_call.iter().position(|&s| s == "*36").unwrap();
            let remain: Vec<usize> = (0..3).filter(|&n| n != idx_10 && n != idx_36).collect();
            assert_eq!(remain.len(), 1);
            let remaining = split_call[remain[0]];
            // When remaining is also *36, both tandem arrangements are
            // indistinguishable from short reads (*36 is a superset of *10).
            // Output both alternatives.
            if remaining == "*36" {
                return Some(format!("*36/*36+*10;*10/*36x2"));
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
            return if spacer_cn.map_or(false, |cn| cn > 4) {
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

        if parts.len() == 1 || (parts.contains("star68") && parts.len() == 1) {
            // All star68
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

        if cnvcall == "exon9hyb_star68" {
            if *called_stars == "*4_*4_*4.013" {
                return Some("*4.013+*4/*68+*4".to_string());
            }
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
        for var_to_add in &["g.42130692G>A"] {
            let count = var_observed.iter().filter(|v| v.as_str() == *var_to_add).count();
            if count > 0 && count <= cn {
                var_observed.push(var_to_add.to_string());
            }
        }
        for var_to_add in &["g.42129754G>A"] {
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
            if exon9cn_consensus > exon9_cn && exon9_cn <= 1 {
                if !var_observed.contains(&"g.42126611C>G".to_string())
                    || (exon9.exon9_raw_site1.min(exon9.exon9_raw_site2) < 1.15
                        && exon9.exon9_raw_site1.max(exon9.exon9_raw_site2) < 1.2)
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
    }

    var_observed.retain(|v| v != "NA");
}

/// Main star allele matching function.
pub fn match_star(
    var_observed: &mut Vec<String>,
    cnvcall: &str,
    spacer_cn: Option<u32>,
    star_combinations: &StarCombinations,
    exon9: &crate::types::Exon9Values,
    var42126938_g_haplotype: bool,
    var42127803_diff_haplotype: bool,
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
            let mut matchtag = call_star68_cn1(var_observed, cnvcall, &star_combinations.dhap);
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
        let mut matchtag = get_star_cn1(var_observed, &star_combinations.dhap);
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

    if !cnvcall.contains("star68") {
        let mut matchtag = get_star(var_observed, dic);

        // For cn3-cn6, try adding a variant on no_match
        if cnvcall.starts_with("cn")
            && cnvcall.len() <= 3
            && matchtag.call_info.as_deref() == Some("no_match")
        {
            let cn_num: usize = cnvcall[2..].parse().unwrap_or(0);
            let mut variant_tried = Vec::new();
            let mut matched_calls = Vec::new();
            for variant_to_try in var_observed.iter().cloned().collect::<Vec<_>>() {
                let count = var_observed
                    .iter()
                    .filter(|v| **v == variant_to_try)
                    .count();
                if count == cn_num - 2 && !variant_tried.contains(&variant_to_try) {
                    variant_tried.push(variant_to_try.clone());
                    let mut var_new = var_observed.clone();
                    var_new.push(variant_to_try);
                    let matchtag_new = get_star(&var_new, dic);
                    if matchtag_new.call_info.as_deref() == Some("unique_match") {
                        matched_calls.push(matchtag_new);
                    }
                }
            }
            if matched_calls.len() == 1 {
                let mt = &matched_calls[0];
                let final_call_clean =
                    get_final_call_clean(&mt.star_call, cnvcall, spacer_cn);
                return StarCall {
                    call_info: mt.call_info.clone(),
                    variants_called: Some(var_observed.join(" ")),
                    raw_call: Some(mt.candidate.clone()),
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
        if call_info.as_deref() == Some("more_than_one_match") && cnvcall == "cn2" {
            if let Some(fcc) = final_call_clean.clone() {
                let re = regex::Regex::new(r"[;/]+").unwrap();
                let mut parts: Vec<&str> = re.split(&fcc).collect();
                parts.sort();
                if parts == ["*1", "*27", "*32", "*41"] {
                    final_call_clean = if var42126938_g_haplotype {
                        Some("*1/*32".to_string())
                    } else {
                        Some("*27/*41".to_string())
                    };
                }
                if parts == ["*1", "*119", "*2", "*41"] {
                    final_call_clean = if var42127803_diff_haplotype {
                        Some("*119/*2".to_string())
                    } else {
                        Some("*1/*41".to_string())
                    };
                }
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
        let matchtag = call_star68(&mut var_observed_68, cnvcall, dic);
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

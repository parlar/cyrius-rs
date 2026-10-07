//! Fuzzy star allele matching.
//!
//! When exact dictionary lookup fails (no_match), this module scores candidate
//! star allele combinations by overlap with observed variants. Reports the best
//! match as "fuzzy_match" rather than "no_match" when the score is high enough.
//!
//! Uses sorted Vec (multiset) comparison to correctly handle duplicate variants
//! (e.g., a variant present at CN=2 appears twice in both observed and candidate).

use std::collections::HashMap;

use crate::caller::match_star_allele::convert_to_main_allele;
use crate::types::RawStarCall;

/// Minimum score threshold for a fuzzy match to be reported.
/// Score = matched_count - 0.5 * candidate_only_count - 0.3 * observed_only_count
const MIN_FUZZY_SCORE: f64 = 1.0;

/// Maximum number of missing or extra variants allowed for a fuzzy match.
const MAX_MISMATCH: usize = 2;

/// Scores within this distance of the best are treated as tied.
const TIE_EPS: f64 = 1e-9;

/// Score a candidate variant multiset against observed variant multiset.
/// Both inputs must be sorted.
fn score_match(candidate: &[&str], observed: &[&str]) -> f64 {
    // Count matches, candidate-only, observed-only using merge-style comparison
    let mut matched = 0usize;
    let mut ci = 0;
    let mut oi = 0;
    let mut candidate_only = 0usize;
    let mut observed_only = 0usize;

    while ci < candidate.len() && oi < observed.len() {
        match candidate[ci].cmp(observed[oi]) {
            std::cmp::Ordering::Equal => {
                matched += 1;
                ci += 1;
                oi += 1;
            }
            std::cmp::Ordering::Less => {
                candidate_only += 1;
                ci += 1;
            }
            std::cmp::Ordering::Greater => {
                observed_only += 1;
                oi += 1;
            }
        }
    }
    candidate_only += candidate.len() - ci;
    observed_only += observed.len() - oi;

    if candidate_only > MAX_MISMATCH || observed_only > MAX_MISMATCH {
        return f64::NEG_INFINITY;
    }

    matched as f64 - 0.5 * candidate_only as f64 - 0.3 * observed_only as f64
}

/// Attempt fuzzy matching against a CN>=2 dictionary (HashMap<String, Vec<String>>).
/// Returns a RawStarCall with call_info "fuzzy_match" if a good match is found.
pub fn fuzzy_match_star(
    var_observed: &[String],
    dic: &HashMap<String, Vec<String>>,
) -> Option<RawStarCall> {
    let mut sorted_observed: Vec<&str> = var_observed.iter().map(|s| s.as_str()).collect();
    sorted_observed.sort();

    if sorted_observed.is_empty() || sorted_observed.iter().all(|&v| v == "NA") {
        return None;
    }

    let mut best_score = f64::NEG_INFINITY;
    let mut best_keys: Vec<&String> = Vec::new();

    for key in dic.keys() {
        if key == "NA" {
            continue;
        }
        // Key is already sorted (alphabetical) from construct_star_table
        let candidate: Vec<&str> = key.split('_').collect();

        let score = score_match(&candidate, &sorted_observed);
        if score > best_score + TIE_EPS {
            best_score = score;
            best_keys = vec![key];
        } else if (score - best_score).abs() <= TIE_EPS {
            best_keys.push(key);
        }
    }

    // Report every tied best key rather than whichever HashMap iteration
    // happened to reach first.
    best_keys.sort();
    let mut best_stars: Vec<String> = Vec::new();
    for key in &best_keys {
        for star in &dic[*key] {
            if !best_stars.contains(star) {
                best_stars.push(star.clone());
            }
        }
    }

    if best_score >= MIN_FUZZY_SCORE && !best_stars.is_empty() {
        for key in &best_keys {
            let candidate: Vec<&str> = key.split('_').collect();
            let missing: Vec<&&str> = candidate.iter().filter(|c| !sorted_observed.contains(c)).collect();
            let extra: Vec<&&str> = sorted_observed.iter().filter(|o| !candidate.contains(o)).collect();
            log::info!(
                "fuzzy_match: score={:.1}, missing={:?}, extra={:?}, match={}",
                best_score, missing, extra, dic[*key].join(";")
            );
        }

        let processed = convert_to_main_allele(&best_stars);

        Some(RawStarCall {
            call_info: Some(format!("fuzzy_match(score={:.1})", best_score)),
            candidate: best_stars,
            star_call: processed,
        })
    } else {
        None
    }
}

/// Attempt fuzzy matching against a CN=1 dictionary (HashMap<String, String>).
pub fn fuzzy_match_star_cn1(
    var_observed: &[String],
    dic: &HashMap<String, String>,
) -> Option<RawStarCall> {
    let mut sorted_observed: Vec<&str> = var_observed.iter().map(|s| s.as_str()).collect();
    sorted_observed.sort();

    if sorted_observed.is_empty() || sorted_observed.iter().all(|&v| v == "NA") {
        return None;
    }

    let mut best_score = f64::NEG_INFINITY;
    let mut best_stars: Vec<String> = Vec::new();

    for (key, star) in dic {
        if key == "NA" {
            continue;
        }
        let candidate: Vec<&str> = key.split('_').collect();

        let score = score_match(&candidate, &sorted_observed);
        if score > best_score + TIE_EPS {
            best_score = score;
            best_stars = vec![star.clone()];
        } else if (score - best_score).abs() <= TIE_EPS && !best_stars.contains(star) {
            best_stars.push(star.clone());
        }
    }
    best_stars.sort();

    if best_score >= MIN_FUZZY_SCORE && !best_stars.is_empty() {
        log::info!(
            "fuzzy_match_cn1: score={:.1}, match={}",
            best_score, best_stars.join(";")
        );

        let processed = convert_to_main_allele(&best_stars);

        Some(RawStarCall {
            call_info: Some(format!("fuzzy_match(score={:.1})", best_score)),
            candidate: best_stars,
            star_call: processed,
        })
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn obs(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn tied_best_keys_are_all_reported() {
        // Both keys score 3 - 0.5 = 2.5 against the observed set.
        let dic: HashMap<String, Vec<String>> = [
            ("a_b_c_d".to_string(), vec!["*10_*20".to_string()]),
            ("a_b_c_e".to_string(), vec!["*20_*4".to_string()]),
        ]
        .into_iter()
        .collect();
        let out = fuzzy_match_star(&obs(&["a", "b", "c"]), &dic).unwrap();
        let mut calls = out.star_call.clone();
        calls.sort();
        assert_eq!(calls, vec!["*10_*20".to_string(), "*20_*4".to_string()]);
    }

    #[test]
    fn tied_best_keys_are_all_reported_cn1() {
        let dic: HashMap<String, String> = [
            ("a_b_d".to_string(), "*10".to_string()),
            ("a_b_e".to_string(), "*4".to_string()),
        ]
        .into_iter()
        .collect();
        let out = fuzzy_match_star_cn1(&obs(&["a", "b"]), &dic).unwrap();
        let mut calls = out.star_call.clone();
        calls.sort();
        assert_eq!(calls, vec!["*10".to_string(), "*4".to_string()]);
    }
}

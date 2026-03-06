//! Fuzzy star allele matching.
//!
//! When exact dictionary lookup fails (no_match), this module scores candidate
//! star allele combinations by overlap with observed variants. Reports the best
//! match as "fuzzy_match" rather than "no_match" when the score is high enough.
//!
//! Uses sorted Vec (multiset) comparison to correctly handle duplicate variants
//! (e.g., a variant present at CN=2 appears twice in both observed and candidate).

use std::collections::HashMap;

use crate::types::RawStarCall;

/// Minimum score threshold for a fuzzy match to be reported.
/// Score = matched_count - 0.5 * candidate_only_count - 0.3 * observed_only_count
const MIN_FUZZY_SCORE: f64 = 1.0;

/// Maximum number of missing or extra variants allowed for a fuzzy match.
const MAX_MISMATCH: usize = 2;

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
    let mut best_key = String::new();
    let mut best_stars: Vec<String> = Vec::new();

    for (key, star_list) in dic {
        if key == "NA" {
            continue;
        }
        // Key is already sorted (alphabetical) from construct_star_table
        let candidate: Vec<&str> = key.split('_').collect();

        let score = score_match(&candidate, &sorted_observed);
        if score > best_score {
            best_score = score;
            best_key = key.clone();
            best_stars = star_list.clone();
        }
    }

    if best_score >= MIN_FUZZY_SCORE && !best_stars.is_empty() {
        let candidate: Vec<&str> = best_key.split('_').collect();
        let missing: Vec<&&str> = candidate.iter().filter(|c| !sorted_observed.contains(c)).collect();
        let extra: Vec<&&str> = sorted_observed.iter().filter(|o| !candidate.contains(o)).collect();

        log::info!(
            "fuzzy_match: score={:.1}, missing={:?}, extra={:?}, match={}",
            best_score, missing, extra, best_stars.join(";")
        );

        let processed = convert_to_main_allele_simple(&best_stars);

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
    let mut best_star = String::new();

    for (key, star) in dic {
        if key == "NA" {
            continue;
        }
        let candidate: Vec<&str> = key.split('_').collect();

        let score = score_match(&candidate, &sorted_observed);
        if score > best_score {
            best_score = score;
            best_star = star.clone();
        }
    }

    if best_score >= MIN_FUZZY_SCORE && !best_star.is_empty() {
        log::info!(
            "fuzzy_match_cn1: score={:.1}, match={}",
            best_score, best_star
        );

        let processed = convert_to_main_allele_simple(&[best_star.clone()]);

        Some(RawStarCall {
            call_info: Some(format!("fuzzy_match(score={:.1})", best_score)),
            candidate: vec![best_star],
            star_call: processed,
        })
    } else {
        None
    }
}

/// Simple main-allele conversion (strip suballele suffixes like .013).
fn convert_to_main_allele_simple(stars: &[String]) -> Vec<String> {
    let mut result_set = std::collections::HashSet::new();
    for star_combo in stars {
        let parts: Vec<&str> = star_combo.split('_').collect();
        let converted: Vec<String> = parts
            .iter()
            .map(|s| s.split('.').next().unwrap().to_string())
            .collect();
        let mut sorted = converted;
        sorted.sort();
        result_set.insert(sorted.join("_"));
    }
    result_set.into_iter().collect()
}

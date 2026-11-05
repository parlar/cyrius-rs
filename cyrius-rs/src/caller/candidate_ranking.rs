//! Ranked candidate genotypes with a site-wise decomposition, and decision provenance.
//! Spec: docs/superpowers/specs/2026-11-05-candidate-ranking-design.md (Stage B).
//! Output only: nothing in this module changes the call.

#![allow(unused_imports, dead_code)]
use std::collections::{HashMap, HashSet};

use rayon::prelude::*;

use crate::caller::construct_star_table::combination_table_for;
use crate::caller::count_model::{site_ln_likelihood, SiteClass, MODEL_ID};
use crate::caller::match_star_allele::{convert_to_main_allele, get_final_call_clean, main_allele_name};
use crate::depth_calling::copy_number_call::{
    call_reg1_cn, cn_posterior, process_raw_call_denovo, CnProb, CnProbResult,
};
use crate::phenotype::sort_diplotypes;
use crate::types::StarCombinations;

/// Table tokens with no read evidence in this model (spec "Unmodelled tokens").
pub const UNMODELLED: [&str; 5] = [
    "exon9gc", "g.42126938C>T", "g.42127526C>T", "g.42127556T>C", "g.42130655-42130656insA",
];
/// Sites at or below this position are in the exon-9 region (`get_total_cn_per_site`).
pub const EXON9_BOUNDARY: i64 = 42126611;
/// Entries listed, plus the call's entry when it ranks lower.
pub const RANKING_LIST_SIZE: usize = 5;
/// Copy counts are tabulated for n < MAX_N (D6 copies, plus 2 D7 at homology sites).
const MAX_N: usize = 10;

/// One scored read site, with the counts the call used before any strand filter.
#[derive(Debug, Clone)]
pub struct ReadSite {
    pub name: String,
    pub pos: i64,
    pub class: SiteClass,
    pub homology: bool,
    pub alt: u32,
    pub ref_: u32,
    /// The filter that zeroed this site's alternate count in the caller. Such a site is not scored.
    pub mask: Option<&'static str>,
}

impl ReadSite {
    fn scored(&self) -> bool {
        self.mask.is_none() && self.alt + self.ref_ > 0
    }
}

/// (name, pos, class, homology) of every scored read site, in `var_list` order. The position
/// is parsed from the name, as `get_total_cn_per_site` does. Unmodelled tokens are skipped.
pub fn site_specs(target: &str, homology: &str) -> Vec<(String, i64, SiteClass, bool)> {
    let mut out = Vec::new();
    for (content, homo) in [(target, false), (homology, true)] {
        for line in content.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() < 7 || UNMODELLED.contains(&f[6]) {
                continue;
            }
            let pos: i64 = f[6].get(2..10).and_then(|p| p.parse().ok()).unwrap_or(0);
            if let Some(class) = SiteClass::from_variant_type(f[5], homo) {
                out.push((f[6].to_string(), pos, class, homo));
            }
        }
    }
    out
}

/// Per allele: the read sites it carries, whether it has the exon-9 conversion, and its
/// count of each unmodelled token (in `UNMODELLED` order).
#[derive(Debug, Clone, Default)]
struct AlleleDef {
    sites: Vec<usize>,
    gc: bool,
    unmodelled: [u32; 5],
}

fn allele_defs(dstar: &HashMap<String, String>, sites: &[ReadSite]) -> HashMap<String, AlleleDef> {
    let idx: HashMap<&str, usize> = sites.iter().enumerate().map(|(i, s)| (s.name.as_str(), i)).collect();
    dstar
        .iter()
        .map(|(allele, key)| {
            let mut def = AlleleDef::default();
            for tok in key.split('_').filter(|t| !t.is_empty() && *t != "NA") {
                if let Some(&i) = idx.get(tok) {
                    def.sites.push(i);
                }
                if let Some(u) = UNMODELLED.iter().position(|u| *u == tok) {
                    def.unmodelled[u] += 1;
                    if u == 0 {
                        def.gc = true;
                    }
                }
            }
            def.sites.sort_unstable();
            (allele.clone(), def)
        })
        .collect()
}

fn copies_of<'d>(raw: &str, defs: &'d HashMap<String, AlleleDef>) -> Option<Vec<&'d AlleleDef>> {
    raw.split('_').map(|a| defs.get(a)).collect()
}

/// Spec "Dosage translation": (k, n) at every read site for one profile. A copy with the
/// exon-9 conversion does not count at pos <= EXON9_BOUNDARY; homology sites add 2 D7 copies.
fn dosage_into(copies: &[&AlleleDef], sites: &[ReadSite], out: &mut [(u8, u8)]) {
    let m = copies.len() as u8;
    let g = copies.iter().filter(|c| c.gc).count() as u8;
    for (o, s) in out.iter_mut().zip(sites) {
        let n6 = if s.pos <= EXON9_BOUNDARY { m - g } else { m };
        *o = (0, n6 + if s.homology { 2 } else { 0 });
    }
    for c in copies {
        for &i in &c.sites {
            if !(c.gc && sites[i].pos <= EXON9_BOUNDARY) {
                out[i].0 += 1;
            }
        }
    }
}

fn dosage(copies: &[&AlleleDef], sites: &[ReadSite]) -> Vec<(u8, u8)> {
    let mut v = vec![(0, 0); sites.len()];
    dosage_into(copies, sites, &mut v);
    v
}

/// Per-site ln-likelihoods indexed `n * MAX_N + k`. All zero for a site that is not scored,
/// and for n = 0.
fn ll_table(sites: &[ReadSite]) -> Vec<Vec<f64>> {
    sites
        .iter()
        .map(|s| {
            let mut v = vec![0.0; MAX_N * MAX_N];
            if s.scored() {
                for n in 1..MAX_N {
                    for k in 0..=n {
                        v[n * MAX_N + k] = site_ln_likelihood(s.class, s.alt, s.ref_, k as u32, n as u32);
                    }
                }
            }
            v
        })
        .collect()
}

/// Composite ln-likelihood of a profile: the site terms summed in site order.
fn profile_score(d: &[(u8, u8)], ll: &[Vec<f64>]) -> f64 {
    d.iter().zip(ll).map(|(&(k, n), t)| t[n as usize * MAX_N + k as usize]).sum()
}

/// Canonical identity of a genotype label (spec "Genotype identity"): its gene copies as
/// sorted main alleles. `*A(*B)` becomes `*B`, `xN` is expanded, and `*5` is dropped as no copy.
pub fn copy_multiset(label: &str) -> Vec<String> {
    static PARENT: std::sync::OnceLock<regex::Regex> = std::sync::OnceLock::new();
    let parent = PARENT.get_or_init(|| regex::Regex::new(r"\*[0-9.]+\((\*[0-9.]+)\)").unwrap());
    let text = parent.replace_all(label, "$1");
    let mut out = Vec::new();
    for part in text.split(['/', '+', '_']) {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let (name, n) = match part.rsplit_once('x') {
            Some((a, b)) if !b.is_empty() && b.bytes().all(|c| c.is_ascii_digit()) => {
                (a, b.parse::<usize>().unwrap_or(1))
            }
            _ => (part, 1),
        };
        let name = main_allele_name(name);
        if name != "*5" {
            out.extend(std::iter::repeat_n(name.to_string(), n));
        }
    }
    out.sort();
    out
}

/// The entries one raw profile produces: the formatter's label(s), split on `;`, each marked
/// formatted. A profile the formatter cannot format keeps its raw main-allele form.
fn entry_labels(raw: &str, cnvcall: &str, spacer_cn: Option<u32>) -> Vec<(String, bool)> {
    let main = convert_to_main_allele(&[raw.to_string()]);
    match get_final_call_clean(&main, cnvcall, spacer_cn) {
        Some(s) if !s.is_empty() && !s.contains('_') => {
            s.split(';').map(|d| (sort_diplotypes(d.trim()), true)).collect()
        }
        _ => main.into_iter().map(|m| (m, false)).collect(),
    }
}

#[cfg(test)]
mod tests_sites {
    use super::*;
    use crate::caller::construct_star_table::get_hap_table;

    fn real_sites() -> Vec<ReadSite> {
        site_specs(crate::data::TARGET_VARIANT, crate::data::TARGET_VARIANT_HOMO)
            .into_iter()
            .map(|(name, pos, class, homology)| ReadSite { name, pos, class, homology, alt: 0, ref_: 30, mask: None })
            .collect()
    }

    fn at(sites: &[ReadSite], d: &[(u8, u8)], name: &str) -> (u8, u8) {
        d[sites.iter().position(|s| s.name == name).unwrap()]
    }

    #[test]
    fn there_are_158_scored_sites_19_of_them_homology() {
        let s = real_sites();
        assert_eq!(s.len(), 158);
        assert_eq!(s.iter().filter(|x| x.homology).count(), 19);
        assert!(s.iter().all(|x| x.name != "g.42130655-42130656insA"));
    }

    #[test]
    fn exon9_converted_copies_do_not_count_in_the_exon9_region() {
        let sc = get_hap_table(crate::data::STAR_TABLE, &[]);
        let sites = real_sites();
        let defs = allele_defs(&sc.dstar, &sites);
        let d = dosage(&copies_of("*10_*36_*36", &defs).unwrap(), &sites);
        assert_eq!(at(&sites, &d, "g.42126611C>G"), (1, 1));
        assert_eq!(at(&sites, &d, "g.42129754G>A"), (3, 3));
    }

    #[test]
    fn homology_sites_add_two_d7_copies() {
        let sc = get_hap_table(crate::data::STAR_TABLE, &[]);
        let sites = real_sites();
        let defs = allele_defs(&sc.dstar, &sites);
        let d = dosage(&copies_of("*1_*1", &defs).unwrap(), &sites);
        assert_eq!(at(&sites, &d, "g.42126578C>T"), (0, 4));
    }

    #[test]
    fn copy_multiset_matches_across_formats() {
        assert_eq!(copy_multiset("*10/*2"), copy_multiset("*2/*10"));
        assert_eq!(copy_multiset("*4.013/*10+*10"), copy_multiset("*10_*10_*4.013"));
        assert_eq!(copy_multiset("*36/*36+*10"), copy_multiset("*10/*36x2"));
        assert_eq!(copy_multiset("*5/*36x2+*10x2"), vec!["*10", "*10", "*36", "*36"]);
        assert_eq!(copy_multiset("*2(*45)/*10"), vec!["*10", "*45"]);
        assert_eq!(copy_multiset("*2.001/*1"), vec!["*1", "*2"]);
    }

    #[test]
    fn arrangements_are_separate_entries_with_one_multiset() {
        let e = entry_labels("*10_*36_*36", "exon9hyb", None);
        assert_eq!(e.len(), 2, "{e:?}");
        assert!(e.iter().all(|(_, formatted)| *formatted));
        assert_eq!(copy_multiset(&e[0].0), copy_multiset(&e[1].0));
    }

    #[test]
    fn a_cn2_profile_formats_to_its_main_alleles() {
        assert_eq!(entry_labels("*1_*2.001", "cn2", None), vec![("*1/*2".to_string(), true)]);
    }

    #[test]
    fn unformatted_profiles_are_kept_and_flagged() {
        let e = entry_labels("*10_*10_*4.013", "exon9hyb", None);
        assert_eq!(e, vec![("*10_*10_*4.013".to_string(), false)]);
    }
}

// ---------------------------------------------------------------------------
// Decision provenance (spec "Decision provenance"), by replaying the caller's own
// functions per site. A replay that disagrees with the value the call used reads
// "unexplained", so drift cannot pass silently.
// ---------------------------------------------------------------------------

/// One provenance record.
#[derive(Debug, Clone, PartialEq)]
pub struct ProvenanceRecord {
    pub site: String,
    pub mechanism: &'static str,
    pub original_alt: Option<u32>,
    pub original_ref: Option<u32>,
    /// The per-site copy count before token rules (None: the site was not called).
    pub per_site_copies: Option<u32>,
    /// The copy count in the call's final variant list (`Variants_called`).
    pub called_copies: u32,
}

fn argmax(v: &[f64]) -> Option<usize> {
    let mut best: Option<usize> = None;
    for (i, x) in v.iter().enumerate() {
        if best.is_none_or(|b| *x > v[b]) {
            best = Some(i);
        }
    }
    best
}

/// The rule that set one site's per-site copy count, or None when the plain posterior did.
/// `actual` is the value the call used. It is compared with the replay first.
pub fn site_mechanism(
    full_cn: u32,
    alt: f64,
    ref_: f64,
    min_read: f64,
    list_total_cn: Option<u32>,
    cap: Option<u32>,
    actual: Option<u32>,
) -> Option<&'static str> {
    let raw = call_reg1_cn(Some(full_cn), alt, ref_, min_read)[0];
    let total = list_total_cn.map(|t| vec![t]);
    let decided = process_raw_call_denovo(&[CnProbResult(raw)], 0.8, 0.65, total.as_deref(), true)[0];
    let replay = match cap {
        Some(c) => decided.map(|v| v.min(c)),
        None => decided,
    };
    if replay != actual {
        return Some("unexplained");
    }
    if cap.is_some() && replay != decided {
        return Some("homology_cap");
    }
    match raw {
        None => Some("no_call"),
        Some(CnProb::Single(v)) => {
            if full_cn == 1 && alt > min_read {
                return Some("one_copy_shortcut");
            }
            let top = cn_posterior(full_cn, alt, ref_).and_then(|p| argmax(&p));
            if v == 0 && top.is_some_and(|t| t != 0) && alt <= min_read && ref_ >= min_read {
                Some("min_read_rule")
            } else {
                None
            }
        }
        Some(CnProb::TwoOptions { cn1, prob1, .. }) => match decided {
            None => Some("no_call"),
            Some(d) if prob1 <= 0.65 && d != cn1 => Some("min_cn_rule"),
            _ => None,
        },
    }
}

/// Values `d6_star_caller` already has for the read sites and their per-site calls.
pub struct SiteCalls<'a> {
    /// Target sites, then homology sites (`params.var_list`).
    pub var_list: &'a [String],
    /// Alternate counts after `update_var42128936`, before any strand filter.
    pub alt_original: &'a [usize],
    /// After `--strand-bias-all` (equal to `alt_original` without the flag).
    pub alt_after_sba: &'a [usize],
    /// After `call_cn_var`'s NOISY_VAR rule: the counts the per-site call used.
    pub alt_final: &'a [usize],
    pub ref_: &'a [usize],
    pub homo_alt: &'a [usize],
    pub homo_ref: &'a [usize],
    /// `get_total_cn_per_site` for the target sites (empty if the tag has no pattern).
    pub total_cn_sites: &'a [u32],
    pub total_cn: u32,
    /// `cn_call_var` and `cn_call_var_homo`.
    pub decisions: &'a [Option<u32>],
    pub homo_decisions: &'a [Option<u32>],
    /// `call_variants::min_read_for`.
    pub min_read: fn(&str) -> f64,
}

fn count_of(tokens: &[String], t: &str) -> u32 {
    tokens.iter().filter(|x| *x == t).count() as u32
}

fn token_counts(tokens: &[String]) -> HashMap<String, u32> {
    let mut m = HashMap::new();
    for t in tokens.iter().filter(|t| t.as_str() != "NA") {
        *m.entry(t.clone()).or_insert(0) += 1;
    }
    m
}

/// Token counts of one raw profile (all its copies' definitions, `NA` dropped).
fn profile_token_counts(raw: &str, dstar: &HashMap<String, String>) -> Option<HashMap<String, u32>> {
    let mut m = HashMap::new();
    for allele in raw.split('_') {
        for t in dstar.get(allele)?.split('_').filter(|t| !t.is_empty() && *t != "NA") {
            *m.entry(t.to_string()).or_insert(0) += 1;
        }
    }
    Some(m)
}

/// Read sites with original counts and filter masks, plus per-site and token-level
/// provenance records. `before_rules` is the call set before `match_star`.
pub fn sites_and_provenance(
    c: &SiteCalls,
    before_rules: &[String],
    variants_called: &[String],
    cnvcall: &str,
) -> (Vec<ReadSite>, Vec<ProvenanceRecord>) {
    let specs: HashMap<String, (i64, SiteClass, bool)> =
        site_specs(crate::data::TARGET_VARIANT, crate::data::TARGET_VARIANT_HOMO)
            .into_iter()
            .map(|(n, p, cl, h)| (n, (p, cl, h)))
            .collect();
    let n_target = c.alt_original.len();
    let mut sites = Vec::new();
    let mut prov = Vec::new();
    for (i, name) in c.var_list.iter().enumerate() {
        let homo = i >= n_target;
        let (alt0, ref0, decision) = if homo {
            let j = i - n_target;
            (c.homo_alt[j], c.homo_ref[j], c.homo_decisions[j])
        } else {
            (c.alt_original[i], c.ref_[i], c.decisions[i])
        };
        let mask = if homo {
            None
        } else if c.alt_original[i] > 0 && c.alt_after_sba[i] == 0 {
            Some("strand_bias_all")
        } else if c.alt_after_sba[i] > 0 && c.alt_final[i] == 0 {
            Some("noisy_var_strand_filter")
        } else {
            None
        };
        let replay_of = |i: usize| {
            c.total_cn_sites.get(i).and_then(|&tcn| {
                site_mechanism(tcn, c.alt_final[i] as f64, ref0 as f64, (c.min_read)(name), Some(tcn), None, decision)
            })
        };
        let mechanism = if homo {
            site_mechanism(c.total_cn, alt0 as f64, ref0 as f64, 4.0, None, Some(c.total_cn.saturating_sub(2)), decision)
        } else if c.total_cn_sites.get(i).is_none() {
            Some("unexplained")
        } else if mask.is_some() {
            if replay_of(i) == Some("unexplained") { Some("unexplained") } else { mask }
        } else {
            replay_of(i)
        };
        if let Some(m) = mechanism {
            prov.push(ProvenanceRecord {
                site: name.clone(),
                mechanism: m,
                original_alt: Some(alt0 as u32),
                original_ref: Some(ref0 as u32),
                per_site_copies: decision,
                called_copies: count_of(variants_called, name),
            });
        }
        if let Some(&(pos, class, h)) = specs.get(name) {
            sites.push(ReadSite { name: name.clone(), pos, class, homology: h, alt: alt0 as u32, ref_: ref0 as u32, mask });
        }
    }
    let mut tokens: Vec<&String> = before_rules.iter().chain(variants_called).filter(|t| t.as_str() != "NA").collect();
    tokens.sort();
    tokens.dedup();
    for t in tokens {
        let (before, after) = (count_of(before_rules, t), count_of(variants_called, t));
        if before != after {
            let mechanism = if cnvcall.contains("star68") && t == "g.42130692G>A" { "star68_step" } else { "update_variants" };
            let site = sites.iter().find(|s| &s.name == t);
            prov.push(ProvenanceRecord {
                site: t.clone(),
                mechanism,
                original_alt: site.map(|s| s.alt),
                original_ref: site.map(|s| s.ref_),
                per_site_copies: Some(before),
                called_copies: after,
            });
        }
    }
    (sites, prov)
}

/// Spec `Matching_path`: the whole-call step that produced the call. `matched` holds the
/// token counts of the call's first raw candidate.
pub fn matching_path(
    cnvcall: &str,
    call_info: Option<&str>,
    variants_called: &[String],
    uncallable: &[String],
    matched: Option<&HashMap<String, u32>>,
) -> &'static str {
    if call_info.is_some_and(|c| c.starts_with("fuzzy_match")) {
        return "fuzzy_match";
    }
    if cnvcall.contains("star68") {
        return "star68_step";
    }
    let unknown: HashSet<&str> = uncallable
        .iter()
        .map(String::as_str)
        .filter(|u| !variants_called.iter().any(|v| v == u))
        .collect();
    if let Some(m) = matched {
        let observed = token_counts(variants_called);
        let diffs: Vec<&String> = m
            .keys()
            .chain(observed.keys())
            .filter(|t| !unknown.contains(t.as_str()))
            .filter(|t| m.get(*t).copied().unwrap_or(0) != observed.get(*t).copied().unwrap_or(0))
            .collect::<HashSet<_>>()
            .into_iter()
            .collect();
        if !diffs.is_empty() {
            // retry_add_one_copy adds one copy of one variant observed exactly N - 2 times, N >= 3.
            let n = if cnvcall.starts_with("cn") && cnvcall.len() <= 3 { cnvcall[2..].parse::<u32>().ok() } else { None };
            let retry = n.is_some_and(|n| {
                n >= 3
                    && diffs.len() == 1
                    && observed.get(diffs[0]).copied().unwrap_or(0) == n - 2
                    && m.get(diffs[0]).copied().unwrap_or(0) == observed.get(diffs[0]).copied().unwrap_or(0) + 1
            });
            return if retry { "retry_add_one_copy" } else { "unexplained" };
        }
    }
    if unknown.is_empty() { "exact" } else { "masked" }
}

#[cfg(test)]
mod tests_provenance {
    use super::*;

    fn s(v: &[&str]) -> Vec<String> {
        v.iter().map(|x| x.to_string()).collect()
    }

    // Read counts chosen from a Python replay of call_reg1_cn and process_raw_call_denovo
    // (PREDICTED until this test runs).
    #[test]
    fn each_rule_is_named_and_an_untouched_site_is_not() {
        assert_eq!(site_mechanism(2, 20.0, 20.0, 4.0, Some(2), None, Some(1)), None); // control
        assert_eq!(site_mechanism(1, 5.0, 95.0, 4.0, Some(1), None, Some(1)), Some("one_copy_shortcut"));
        assert_eq!(site_mechanism(2, 3.0, 4.0, 4.0, Some(2), None, Some(0)), Some("min_read_rule"));
        assert_eq!(site_mechanism(2, 5.0, 1.0, 4.0, Some(2), None, Some(1)), Some("min_cn_rule"));
        assert_eq!(site_mechanism(4, 30.0, 10.0, 4.0, None, Some(2), Some(2)), Some("homology_cap"));
        assert_eq!(site_mechanism(2, 0.0, 0.0, 4.0, Some(2), None, None), Some("no_call"));
    }

    #[test]
    fn a_replay_that_disagrees_with_the_call_reads_unexplained() {
        assert_eq!(site_mechanism(2, 20.0, 20.0, 4.0, Some(2), None, Some(2)), Some("unexplained"));
    }

    #[test]
    fn matching_paths() {
        let m = |v: &[&str]| token_counts(&s(v));
        let obs = s(&["g.A", "g.B"]);
        assert_eq!(matching_path("cn2", Some("unique_match"), &obs, &[], Some(&m(&["g.A", "g.B"]))), "exact");
        assert_eq!(matching_path("cn2", Some("more_than_one_match"), &obs, &s(&["g.C"]), Some(&m(&["g.A", "g.B", "g.C"]))), "masked");
        assert_eq!(matching_path("cn3", Some("unique_match"), &obs, &[], Some(&m(&["g.A", "g.B", "g.B"]))), "retry_add_one_copy");
        assert_eq!(matching_path("cn2", Some("fuzzy_match(score=1.0)"), &obs, &[], Some(&m(&["g.A"]))), "fuzzy_match");
        assert_eq!(matching_path("cn2", Some("unique_match"), &obs, &[], Some(&m(&["g.A"]))), "unexplained");
        // retry tightening
        let cn3 = |mm: &[&str]| matching_path("cn3", Some("unique_match"), &obs, &[], Some(&m(mm)));
        assert_eq!(cn3(&["g.A", "g.B", "g.B"]), "retry_add_one_copy");
        assert_eq!(cn3(&["g.A", "g.B", "g.B", "g.B"]), "unexplained");
        assert_eq!(cn3(&["g.A"]), "unexplained");
        assert_eq!(cn3(&["g.A", "g.A", "g.B", "g.B"]), "unexplained");
        // cn4 needs the variant observed twice (N - 2); B is observed once, so +1 is not a retry.
        assert_eq!(matching_path("cn4", Some("unique_match"), &obs, &[], Some(&m(&["g.A", "g.B", "g.B"]))), "unexplained");
    }

    #[test]
    fn matching_path_without_a_matched_candidate() {
        let obs = s(&["g.A"]);
        assert_eq!(matching_path("cn2", None, &obs, &[], None), "exact");
        assert_eq!(matching_path("cn2", None, &obs, &s(&["g.C"]), None), "masked");
    }

    fn one_site<'a>(
        var_list: &'a [String],
        alt_final: &'a [usize],
        alt_sba: &'a [usize],
        total_cn_sites: &'a [u32],
        decisions: &'a [Option<u32>],
    ) -> SiteCalls<'a> {
        SiteCalls {
            var_list,
            alt_original: &[5],
            alt_after_sba: alt_sba,
            alt_final,
            ref_: &[95],
            homo_alt: &[],
            homo_ref: &[],
            total_cn_sites,
            total_cn: 4,
            decisions,
            homo_decisions: &[],
            min_read: crate::caller::call_variants::min_read_for,
        }
    }

    #[test]
    fn unmasked_target_site_reaches_site_mechanism() {
        let vl = s(&["g.42129754G>A"]);
        let c = one_site(&vl, &[5], &[5], &[1], &[Some(1)]);
        let (_, prov) = sites_and_provenance(&c, &[], &[], "cn2");
        assert_eq!(prov.iter().find(|p| p.site == vl[0]).map(|p| p.mechanism), Some("one_copy_shortcut"));
    }

    #[test]
    fn homology_cap_through_sites_and_provenance() {
        let vl = s(&["g.42126578C>T"]);
        let c = SiteCalls {
            var_list: &vl,
            alt_original: &[],
            alt_after_sba: &[],
            alt_final: &[],
            ref_: &[],
            homo_alt: &[30],
            homo_ref: &[10],
            total_cn_sites: &[],
            total_cn: 4,
            decisions: &[],
            homo_decisions: &[Some(2)],
            min_read: crate::caller::call_variants::min_read_for,
        };
        let (_, prov) = sites_and_provenance(&c, &[], &[], "cn2");
        assert_eq!(prov.iter().find(|p| p.site == vl[0]).map(|p| p.mechanism), Some("homology_cap"));
    }

    #[test]
    fn masked_site_whose_decision_disagrees_with_the_replay_is_unexplained() {
        let vl = s(&["g.42129754G>A"]);
        // strand_bias_all mask (alt 5 -> 0); replay of alt 0, ref 95 at tcn 1 is not Some(1).
        let c = one_site(&vl, &[0], &[0], &[1], &[Some(1)]);
        let (_, prov) = sites_and_provenance(&c, &[], &[], "cn2");
        assert_eq!(prov.iter().find(|p| p.site == vl[0]).map(|p| p.mechanism), Some("unexplained"));
    }

    #[test]
    fn missing_copy_number_pattern_is_unexplained() {
        let vl = s(&["g.42129754G>A"]);
        let c = one_site(&vl, &[5], &[5], &[], &[Some(1)]);
        let (_, prov) = sites_and_provenance(&c, &[], &[], "cn2");
        assert_eq!(prov.iter().find(|p| p.site == vl[0]).map(|p| p.mechanism), Some("unexplained"));
    }

    #[test]
    fn star68_token_diff_is_star68_step() {
        let before = s(&["g.42130692G>A", "g.42130692G>A"]);
        let after = s(&["g.42130692G>A"]);
        let c = SiteCalls {
            var_list: &[],
            alt_original: &[],
            alt_after_sba: &[],
            alt_final: &[],
            ref_: &[],
            homo_alt: &[],
            homo_ref: &[],
            total_cn_sites: &[],
            total_cn: 2,
            decisions: &[],
            homo_decisions: &[],
            min_read: crate::caller::call_variants::min_read_for,
        };
        let (_, prov) = sites_and_provenance(&c, &before, &after, "star68");
        assert_eq!(prov.iter().find(|p| p.site == "g.42130692G>A").map(|p| p.mechanism), Some("star68_step"));
    }

    #[test]
    fn filters_mask_sites_and_token_rules_are_recorded() {
        let var_list = s(&["g.42129754G>A", "g.42127473C>T", "g.42126578C>T"]);
        let calls = SiteCalls {
            var_list: &var_list,
            alt_original: &[9, 10],
            alt_after_sba: &[0, 10],
            alt_final: &[0, 0],
            ref_: &[11, 12],
            homo_alt: &[0],
            homo_ref: &[30],
            total_cn_sites: &[2, 2],
            total_cn: 4,
            decisions: &[Some(0), Some(0)],
            homo_decisions: &[Some(0)],
            min_read: crate::caller::call_variants::min_read_for,
        };
        let before = s(&["g.42129809T>C"]);
        let after = s(&["g.42129809T>C", "g.42129819G>T"]);
        let (sites, prov) = sites_and_provenance(&calls, &before, &after, "cn2");
        let mech = |n: &str| prov.iter().find(|p| p.site == n).map(|p| p.mechanism);
        assert_eq!(mech("g.42129754G>A"), Some("strand_bias_all"));
        assert_eq!(mech("g.42127473C>T"), Some("noisy_var_strand_filter"));
        assert_eq!(mech("g.42126578C>T"), None);
        assert_eq!(mech("g.42129819G>T"), Some("update_variants"));
        let site = |n: &str| sites.iter().find(|x| x.name == n).unwrap().clone();
        assert_eq!((site("g.42129754G>A").alt, site("g.42129754G>A").mask), (9, Some("strand_bias_all")));
        assert!(site("g.42126578C>T").homology);
    }
}

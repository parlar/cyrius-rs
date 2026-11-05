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

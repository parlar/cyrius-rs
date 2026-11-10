# Candidate Ranking (Stage B) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Behind a new opt-in flag `--candidate-ranking`, add a `Candidate_ranking` JSON field.
It holds the ranked genotypes from the call's CNV group, scored by Stage A's count model, with
a per-site decomposition of each gap and the decision provenance of the call.

**Architecture:**
- A new `count_model` module holds Stage A's chosen model M4 (exp-H20 result `0a4b020`).
- A new `candidate_ranking` module is a pure function of values `d6_star_caller` already has:
  read sites, table, call and per-site decisions. It does no BAM I/O.
- Provenance is found by **replaying** the caller's own decision functions per site, and is
  checked at run time against the value the call actually used. Two behaviour-preserving
  refactors: `min_read_for`, and `call_reg1_cn` computing its posterior through the new
  `cn_posterior` (user decision 2026-11-05: one copy of the math, not two). B3 proves
  neither changes any output.
- `main.rs` keeps three snapshots of the alternate counts, calls the module once, and inserts
  the field only when the flag is set.

**Tech Stack:** Rust (MSRV 1.85), statrs 0.17 (`ln_gamma`), rayon, serde_json, regex. The
checks use Python 3.

**Spec:** `docs/superpowers/specs/2026-11-05-candidate-ranking-design.md` (Stage B). Stage A
plan and result: `~/cyrius-research/exp-H20` (`PLAN.md`, result `0a4b020`).

## Global Constraints

- Branch `feat/candidate-ranking`. Commit per task. Never merge or push.
- Commit trailer: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Output only: no change to any call, filter, `Confidence`, resolver result or TSV, with or
  without the flag.
- Without the flag, every JSON field except `Tool_commit` and `Flags` is identical to master
  `baf3459`, and there is no `Candidate_ranking` key.
- `cargo clippy --all-targets --locked`: the count of lines matching `^(warning|error)` stays
  at or below 22 (`~/cyrius-research/unknown-sites/lint_baseline.txt`).
- `cargo test --locked` passes in full after every task.
- Model and scope come from Stage A only. Ranked tags are `cn2`, `cn3`, `exon9hyb`,
  `exon9hyb_exon9hyb` and `dup_exon9hyb`, the tags Stage A's passing groups received
  (measured from the exp-H20 runs: cn2 169, cn3 6, exon9hyb 7, exon9hyb_exon9hyb 5,
  dup_exon9hyb 1).
  - `star5` → `group_failed_transfer_test`.
  - Tags containing `star13` or `star68` → `hybrid_dosage_undefined`.
  - `star5_star5`, `star13_star13`, `star13intron1_star13intron1`, `star5_star5_star68` →
    `direct_genotype_tag`.
  - Every other tag → `group_untested_in_stage_a`.
- `Evidence_scale` is always `"uncalibrated_composite"`. Fragment reuse in Stage A was a
  median of 96.1%, so the gaps overstate independent evidence.

## Locked checks (this commit, before any Stage B run)

These are run by Task 8 on 160 files: the 53 GeT-RM benchmark BAMs, the 53 at 25% and 53 at
10% depth (`~/cyrius-research/unknown-sites/thinned/{p25,p10}`), and the damaged NA19920. The
flags are `--strand-bias-all --fuzzy-match --changepoint-hybrid --het-check`. There are four
configurations:
- **base:** `~/cyrius-research/exp-H20/frozen/cyrius-rs`, i.e. `baf3459`, sha256
  `1f9bcead…e85e`;
- **noflag:** the branch, without the flag;
- **flag1, flag2:** the branch with `--candidate-ranking`, run twice.

- **B1.** For every listed candidate in every flag1 file, the sum of its
  `Site_log10_ratio` rows equals its `Delta_composite_log10_likelihood` within 1e-9.
- **B2.** `Candidate_ranking` is identical between flag1 and flag2 for all 160 files.
- **B3a.** noflag vs base: every field except `Tool_commit` and `Flags` is identical, and no
  noflag file has `Candidate_ranking`.
- **B3b.** flag1 vs base: every field except `Tool_commit`, `Flags` and `Candidate_ranking`
  is identical.
- **B4.**
  - No `Decision_provenance` record and no `Matching_path` reads `unexplained`. That would
    mean the replay disagreed with the value the call used.
  - Every `update_variants` or `star68_step` record has `Called_copies` equal to that
    token's count in `Variants_called`.
- **Vacuity guards** (a check that saw nothing proves nothing). Each of these fails the run:
  - fewer than 160 files in any configuration;
  - fewer than 80 `ranked` files in flag1. PREDICTED (not run): about 120;
  - fewer than 10 listed candidates with at least one discriminating row;
  - no provenance record at all;
  - `Genotype` missing anywhere in base.
- **Descriptive (no decision):**
  - `Status` and `Reason` counts;
  - `Matching_path` counts;
  - the `Call_rank` distribution;
  - wall time per configuration.

The result commit records the check output verbatim. No threshold changes after the run.

---

### Task 1: The count model

**Files:**
- Create: `cyrius-rs/src/caller/count_model.rs`
- Modify: `cyrius-rs/src/caller/mod.rs` (add `pub mod count_model;` after `pub mod consistency;`)

**Interfaces:**
- Produces:
  - `pub enum SiteClass { Snp, SnpMulti, Indel, IndelMulti, GeneConversion, Homology }`
  - `SiteClass::from_variant_type(&str, bool) -> Option<SiteClass>`
  - `SiteClass::params(self) -> (f64, f64, f64)`
  - `pub fn site_ln_likelihood(class, alt: u32, ref_: u32, k: u32, n: u32) -> f64`
  - `pub const MODEL_ID: &str`

- [ ] **Step 1: Write the failing test** — create `count_model.rs` containing only the test module and a stub, so it compiles and fails:

```rust
//! Per-site count model for the candidate ranking: Stage A's chosen model M4
//! (exp-H20 result 0a4b020). Beta-binomial with a per-class error rate, allele
//! bias and overdispersion, fitted on 169 cn2 samples (144 HPRC + 25 GeT-RM).
//! The parameters are the `all_cn2` rows of exp-H20 `results/params.tsv`.

use statrs::function::gamma::ln_gamma;

/// Reported as `Candidate_ranking.Model`.
pub const MODEL_ID: &str = "M4 (exp-H20 0a4b020)";

/// Site class: the `variant_type` column, with the homology file as its own class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiteClass {
    Snp,
    SnpMulti,
    Indel,
    IndelMulti,
    GeneConversion,
    Homology,
}

impl SiteClass {
    pub fn from_variant_type(_variant_type: &str, _homology: bool) -> Option<SiteClass> {
        None
    }
}

pub fn site_ln_likelihood(_class: SiteClass, _alt: u32, _ref: u32, _k: u32, _n: u32) -> f64 {
    0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference values from exp-H20 `scripts/models.py` (`ll_vec("M4", ...)`) with the
    /// `all_cn2` parameters, printed by Python before this plan was written.
    #[test]
    fn matches_the_stage_a_python_model() {
        let cases = [
            (SiteClass::Snp, 20, 20, 1, 2, -2.424344032803589),
            (SiteClass::Snp, 0, 30, 0, 2, -0.008584880431968145),
            (SiteClass::Indel, 7, 33, 1, 2, -3.5155621349832145),
            (SiteClass::Homology, 10, 50, 1, 4, -3.1241100424440447),
            (SiteClass::GeneConversion, 15, 15, 1, 2, -1.954377886442014),
            (SiteClass::Snp, 3, 40, 0, 3, -7.284606459391441),
            (SiteClass::IndelMulti, 0, 25, 0, 1, -0.000961107793626681),
            (SiteClass::SnpMulti, 12, 0, 2, 2, -0.006228342306030754),
        ];
        for (class, a, r, k, n, want) in cases {
            let got = site_ln_likelihood(class, a, r, k, n);
            assert!((got - want).abs() < 1e-9, "{class:?} {a}/{r} k={k} n={n}: {got} vs {want}");
        }
    }

    #[test]
    fn is_a_distribution_over_the_alt_count() {
        use SiteClass::*;
        for class in [Snp, SnpMulti, Indel, IndelMulti, GeneConversion, Homology] {
            let total: f64 = (0..=30).map(|a| site_ln_likelihood(class, a, 30 - a, 1, 2).exp()).sum();
            assert!((total - 1.0).abs() < 1e-9, "{class:?}: {total}");
        }
    }

    #[test]
    fn variant_types_map_to_classes() {
        assert_eq!(SiteClass::from_variant_type("snp", false), Some(SiteClass::Snp));
        assert_eq!(SiteClass::from_variant_type("snp_multi", false), Some(SiteClass::SnpMulti));
        assert_eq!(SiteClass::from_variant_type("indel", false), Some(SiteClass::Indel));
        assert_eq!(SiteClass::from_variant_type("indel_multi", false), Some(SiteClass::IndelMulti));
        assert_eq!(SiteClass::from_variant_type("gene_conversion", false), Some(SiteClass::GeneConversion));
        assert_eq!(SiteClass::from_variant_type("snp", true), Some(SiteClass::Homology));
        assert_eq!(SiteClass::from_variant_type("weird", false), None);
    }
}
```

Add `pub mod count_model;` to `cyrius-rs/src/caller/mod.rs` after `pub mod consistency;`.

- [ ] **Step 2: Run, expect failure**

Run: `cd /home/parlar_ai/dev/cyrius-rs/cyrius-rs && cargo test --locked --lib count_model 2>&1 | tail -5`
Expected: 3 failed (`matches_the_stage_a_python_model`, `is_a_distribution_over_the_alt_count`, `variant_types_map_to_classes`).

- [ ] **Step 3: Implement** — replace the stubs with:

```rust
impl SiteClass {
    pub fn from_variant_type(variant_type: &str, homology: bool) -> Option<SiteClass> {
        if homology {
            return Some(SiteClass::Homology);
        }
        match variant_type {
            "snp" => Some(SiteClass::Snp),
            "snp_multi" => Some(SiteClass::SnpMulti),
            "indel" => Some(SiteClass::Indel),
            "indel_multi" => Some(SiteClass::IndelMulti),
            "gene_conversion" => Some(SiteClass::GeneConversion),
            _ => None,
        }
    }

    /// (error rate, allele bias, overdispersion): exp-H20 `results/params.tsv`, M4, all_cn2.
    pub fn params(self) -> (f64, f64, f64) {
        match self {
            SiteClass::Snp => (0.0004015441770681196, 0.9339356268269258, 0.023072497450790153),
            SiteClass::SnpMulti => (0.0005162070342784222, 0.980262902642141, 0.0027068788553678198),
            SiteClass::Indel => (0.0016441425094134746, 0.754160272105761, 0.31571547734888167),
            SiteClass::IndelMulti => (0.0001, 0.8434766740495846, 0.11761644655153744),
            SiteClass::GeneConversion => (0.0002807889779143907, 0.9301523446404545, 0.0),
            SiteClass::Homology => (0.0006937702284259667, 1.0072875994763135, 0.0033292053680482486),
        }
    }
}

fn ln_choose(n: f64, a: f64) -> f64 {
    ln_gamma(n + 1.0) - ln_gamma(a + 1.0) - ln_gamma(n - a + 1.0)
}

/// Natural-log probability of `alt` alternate reads among `alt + ref_` at a site with `n`
/// copies, `k` of them carrying the variant. Same formula as exp-H20 `models.ll_vec("M4")`:
/// p = (k/n)(1 - eps) + (1 - k/n) eps, biased by beta, beta-binomial with rho (binomial when
/// rho < 1e-6).
pub fn site_ln_likelihood(class: SiteClass, alt: u32, ref_: u32, k: u32, n: u32) -> f64 {
    let (eps, beta, rho) = class.params();
    let (a, r) = (alt as f64, ref_ as f64);
    let f = k as f64 / n as f64;
    let p = f * (1.0 - eps) + (1.0 - f) * eps;
    let p = (beta * p / (beta * p + 1.0 - p)).clamp(1e-12, 1.0 - 1e-12);
    if rho < 1e-6 {
        ln_choose(a + r, a) + a * p.ln() + r * (1.0 - p).ln()
    } else {
        let s = (1.0 - rho) / rho;
        let (al, be) = (p * s, (1.0 - p) * s);
        ln_choose(a + r, a) + ln_gamma(a + al) + ln_gamma(r + be) - ln_gamma(a + r + al + be)
            - ln_gamma(al) - ln_gamma(be) + ln_gamma(al + be)
    }
}
```

- [ ] **Step 4: Run, expect pass**

Run: `cargo test --locked --lib count_model 2>&1 | tail -3`
Expected: `3 passed`.

- [ ] **Step 5: Mutation check.** Change the Indel β `0.754160272105761` to `1.0`, run Step 4,
  and expect `matches_the_stage_a_python_model` to FAIL. Restore it and expect a pass again.

- [ ] **Step 6: Commit**

```bash
cd /home/parlar_ai/dev/cyrius-rs && git add cyrius-rs/src/caller/count_model.rs cyrius-rs/src/caller/mod.rs
git commit -m "feat(ranking): Stage A's count model M4 as a per-site likelihood (exp-H20 0a4b020)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 2: Replay helpers in the caller's own modules

**Files:**
- Modify: `cyrius-rs/src/depth_calling/copy_number_call.rs` (add `cn_posterior` after `call_reg1_cn`, plus a test)
- Modify: `cyrius-rs/src/caller/call_variants.rs` (extract `min_read_for` and use it in `call_cn_var`, plus a test)

**Interfaces:**
- Produces:
  - `pub fn cn_posterior(cn: u32, count_reg1: f64, count_reg2: f64) -> Option<Vec<f64>>`
  - `pub fn min_read_for(var: &str) -> f64`

- [ ] **Step 1: Write the failing tests.** Append to the test module of `copy_number_call.rs`
  (create `#[cfg(test)] mod tests { use super::*; ... }` if none exists):

```rust
    /// cn_posterior is the posterior call_reg1_cn decides from. Wherever call_reg1_cn
    /// returns two options or a stringent single value, they must agree. min_read = -1
    /// keeps the min-read rule from firing; cn >= 2 keeps the one-copy shortcut out.
    #[test]
    fn cn_posterior_is_the_posterior_call_reg1_cn_decides_from() {
        let mut checked = 0;
        for cn in 2..=4u32 {
            for c1 in 0..30 {
                for c2 in 0..30 {
                    let (c1, c2) = (c1 as f64, c2 as f64);
                    let post = cn_posterior(cn, c1, c2);
                    match call_reg1_cn(Some(cn), c1, c2, -1.0)[0] {
                        None => assert!(post.is_none()),
                        Some(CnProb::Single(v)) => {
                            let p = post.unwrap();
                            let top = (0..p.len()).fold(0, |b, i| if p[i] > p[b] { i } else { b });
                            assert_eq!(top as u32, v, "cn={cn} {c1}/{c2}");
                            checked += 1;
                        }
                        Some(CnProb::TwoOptions { cn1, prob1, .. }) => {
                            let p = post.unwrap();
                            assert!((p[cn1 as usize] - prob1).abs() <= 0.0005 + 1e-12, "cn={cn} {c1}/{c2}");
                            checked += 1;
                        }
                    }
                }
            }
        }
        assert!(checked > 2000);
    }
```

Append to the test module of `call_variants.rs`:

```rust
    #[test]
    fn min_read_follows_the_site_lists() {
        assert_eq!(min_read_for("g.42126611C>G"), 2.0); // CLEAN_VAR
        assert_eq!(min_read_for("g.42127473C>T"), 7.0); // NOISY_VAR
        assert_eq!(min_read_for("g.42129754G>A"), 4.0); // neither
    }
```

- [ ] **Step 2: Run, expect failure**

Run: `cargo test --locked --lib 2>&1 | grep -E "cannot find function|error\[E0425\]" | head -3`
Expected: `cannot find function` for `cn_posterior` and `min_read_for`.

- [ ] **Step 3: Implement.** In `copy_number_call.rs`, after `call_reg1_cn`:

```rust
/// The normalised posterior over 0..=cn that `call_reg1_cn` decides from: the same expected
/// depths, the same smaller-count Poisson term and the same reversal. `None` where
/// `call_reg1_cn` computes no Poisson term. Used only to report decision provenance.
pub fn cn_posterior(cn: u32, count_reg1: f64, count_reg2: f64) -> Option<Vec<f64>> {
    let nsum = count_reg1 + count_reg2;
    if cn == 0 || nsum == 0.0 {
        return None;
    }
    let count = if count_reg1 <= count_reg2 { count_reg1 as u64 } else { count_reg2 as u64 };
    let mut prob = Vec::with_capacity(cn as usize + 1);
    for i in 0..=cn {
        let depth_expected = if i == 0 {
            (ERROR_RATE / 3.0) * nsum
        } else if i == cn {
            nsum - ERROR_RATE * nsum
        } else {
            nsum * i as f64 / cn as f64
        };
        prob.push(Poisson::new(depth_expected).unwrap().pmf(count));
    }
    let sum: f64 = prob.iter().sum();
    if sum == 0.0 {
        return None;
    }
    let mut post: Vec<f64> = prob.iter().map(|p| p / sum).collect();
    if count_reg2 < count_reg1 {
        post.reverse();
    }
    Some(post)
}
```

In `call_variants.rs`, add before `call_cn_var`:

```rust
/// The read minimum `call_cn_var` gives `call_reg1_cn` for a site (2 clean, 7 noisy, else 4).
pub fn min_read_for(var: &str) -> f64 {
    if CLEAN_VAR.contains(&var) {
        2.0
    } else if NOISY_VAR.contains(&var) {
        7.0
    } else {
        4.0
    }
}
```

In `call_cn_var`, delete the `let min_read = if CLEAN_VAR... else { 4 };` block and change the
push to:

```rust
        cn_prob.push(
            call_reg1_cn(Some(total_cn[i]), total_var as f64, total_ref as f64, min_read_for(&var_list[i]))[0]
                .into(),
        );
```

**Shared, not copied (user decision 2026-11-05):** after adding `cn_posterior`, make
`call_reg1_cn`'s `Some(cn)` arm use it, so the math exists once. Replace everything from
`let nsum = count_reg1 + count_reg2;` through the `post_prob.reverse();` block with:

```rust
            let Some(post_prob) = cn_posterior(cn, count_reg1, count_reg2) else {
                return vec![None];
            };
```

Leave the rest of the arm (from `let mut post_prob_sorted` on) unchanged. The `Some(0)` and
`Some(1) if count_reg1 > min_read` arms stay before it, as they are.

- [ ] **Step 4: Run, expect pass**

Run: `cargo test --locked 2>&1 | grep -E "^test result" | head -5`
Expected: every `test result: ok`.

- [ ] **Step 5: Mutation check.** Remove `post.reverse();` from `cn_posterior`, run, and
  expect `cn_posterior_is_the_posterior_call_reg1_cn_decides_from` to FAIL. Restore it.

- [ ] **Step 6: Commit**

```bash
git add cyrius-rs/src/depth_calling/copy_number_call.rs cyrius-rs/src/caller/call_variants.rs
git commit -m "feat(ranking): cn_posterior and min_read_for, for replaying per-site decisions

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 3: Read sites, dosage translation and genotype identity

**Files:**
- Create: `cyrius-rs/src/caller/candidate_ranking.rs`
- Modify: `cyrius-rs/src/caller/mod.rs` (add `pub mod candidate_ranking;` after `pub mod call_variants;`)

**Interfaces:**
- Consumes: `count_model::{SiteClass, site_ln_likelihood, MODEL_ID}`,
  `match_star_allele::{convert_to_main_allele, get_final_call_clean, main_allele_name}`,
  `phenotype::sort_diplotypes`
- Produces:
  - `pub const UNMODELLED: [&str; 5]`, `pub const EXON9_BOUNDARY: i64`,
    `pub const RANKING_LIST_SIZE: usize`
  - `pub struct ReadSite`
  - `pub fn site_specs(target, homology) -> Vec<(String, i64, SiteClass, bool)>`
  - `pub fn copy_multiset(label) -> Vec<String>`
  - private: `AlleleDef`, `allele_defs`, `dosage_into`, `dosage`, `ll_table`,
    `profile_score`, `copies_of`, `entry_labels`

- [ ] **Step 1: Write the module with its tests** — `cyrius-rs/src/caller/candidate_ranking.rs`:

```rust
//! Ranked candidate genotypes with a site-wise decomposition, and decision provenance.
//! Spec: docs/superpowers/specs/2026-11-05-candidate-ranking-design.md (Stage B).
//! Output only: nothing in this module changes the call.

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
```

Add `pub mod candidate_ranking;` to `cyrius-rs/src/caller/mod.rs` after `pub mod call_variants;`.

The `use` lines for `HashSet`, `rayon`, `combination_table_for`, `MODEL_ID`,
`call_reg1_cn`/`cn_posterior`/`process_raw_call_denovo`/`CnProb`/`CnProbResult` and
`StarCombinations` are used from Task 4 on. Until then add
`#![allow(unused_imports)]` as the first line under the doc comment, and remove it in Task 5,
Step 5.

- [ ] **Step 2: Run the tests**

Run: `cargo test --locked --lib candidate_ranking 2>&1 | tail -4`
Expected: `7 passed`. If one of the formatter-dependent tests
(`a_cn2_profile_formats...`, `unformatted_profiles...`, `arrangements...`) fails, print the
actual `entry_labels` output first and report it. Do not edit the expected value to match: a
formatter surprise changes the spec's identity rules.

- [ ] **Step 3: Mutation check.** In `dosage_into`, change `if !(c.gc && ...)` to
  `if !(false && ...)`, and expect `exon9_converted_copies...` to FAIL. Restore it.

- [ ] **Step 4: Commit**

```bash
git add cyrius-rs/src/caller/candidate_ranking.rs cyrius-rs/src/caller/mod.rs
git commit -m "feat(ranking): read sites, dosage translation and copy-multiset identity

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 4: Decision provenance by replay

**Files:**
- Modify: `cyrius-rs/src/caller/candidate_ranking.rs` (append the provenance section and its tests)

**Interfaces:**
- Produces:
  - `pub struct ProvenanceRecord { site, mechanism, original_alt, original_ref, per_site_copies, called_copies }`
  - `pub fn site_mechanism(full_cn, alt, ref_, min_read, list_total_cn, cap, actual) -> Option<&'static str>`
  - `pub struct SiteCalls<'a>`
  - `pub fn sites_and_provenance(&SiteCalls, before_rules, variants_called, cnvcall) -> (Vec<ReadSite>, Vec<ProvenanceRecord>)`
  - `pub fn matching_path(cnvcall, call_info, variants_called, uncallable, matched) -> &'static str`
  - `fn token_counts`, `fn profile_token_counts`

- [ ] **Step 1: Append the code and the tests**

```rust
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
        let mechanism = if mask.is_some() {
            mask
        } else if homo {
            site_mechanism(c.total_cn, alt0 as f64, ref0 as f64, 4.0, None, Some(c.total_cn.saturating_sub(2)), decision)
        } else {
            c.total_cn_sites.get(i).and_then(|&tcn| {
                site_mechanism(tcn, c.alt_final[i] as f64, ref0 as f64, (c.min_read)(name), Some(tcn), None, decision)
            })
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
        let differs = m.keys().chain(observed.keys()).filter(|t| !unknown.contains(t.as_str())).any(|t| {
            m.get(t).copied().unwrap_or(0) != observed.get(t).copied().unwrap_or(0)
        });
        if differs {
            // retry_add_one_copy needs target_count = N - 2 >= 1, so it never fires for cn2
            // (corrected during Task 4: the first draft also let cn2 through).
            let retry_possible = cnvcall.starts_with("cn") && cnvcall.len() <= 3
                && cnvcall[2..].parse::<u32>().is_ok_and(|n| n >= 3);
            return if retry_possible { "retry_add_one_copy" } else { "unexplained" };
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
```

- [ ] **Step 2: Run, expect pass**

Run: `cargo test --locked --lib candidate_ranking 2>&1 | tail -4`
Expected: `11 passed`. If `each_rule_is_named...` fails on a rule, print that replay's
`call_reg1_cn` output and report it. Do not change counts until the cause is understood: the
Python replay that chose them is a prediction.

- [ ] **Step 3: Mutation check.** In `site_mechanism`, change `if replay != actual` to
  `if false`, and expect `a_replay_that_disagrees...` to FAIL. Restore it.

- [ ] **Step 4: Commit**

```bash
git add cyrius-rs/src/caller/candidate_ranking.rs
git commit -m "feat(ranking): decision provenance by replay, checked against the value the call used

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 5: The ranking

**Files:**
- Modify: `cyrius-rs/src/caller/candidate_ranking.rs` (append the ranking section and its tests)

**Interfaces:**
- Produces:
  - `pub fn ranking_scope(cnvcall) -> Result<(), &'static str>`
  - `pub struct RankingInput<'a>`
  - `pub fn rank(input: RankingInput) -> serde_json::Value`
  - `fn rank_with_table(input, table) -> serde_json::Value`

- [ ] **Step 1: Append the code and the tests**

```rust
// ---------------------------------------------------------------------------
// Ranking (spec "Profiles and scores", "Ranking", "Decomposition", "Call rank", "Output").
// ---------------------------------------------------------------------------

/// Which CNV tags Stage A tested and passed (exp-H20 0a4b020). Err is the not_ranked reason.
pub fn ranking_scope(cnvcall: &str) -> Result<(), &'static str> {
    const DIRECT: [&str; 4] = ["star5_star5", "star13_star13", "star13intron1_star13intron1", "star5_star5_star68"];
    const RANKED: [&str; 5] = ["cn2", "cn3", "exon9hyb", "exon9hyb_exon9hyb", "dup_exon9hyb"];
    if DIRECT.contains(&cnvcall) {
        Err("direct_genotype_tag")
    } else if cnvcall.contains("star13") || cnvcall.contains("star68") {
        Err("hybrid_dosage_undefined")
    } else if cnvcall == "star5" {
        Err("group_failed_transfer_test")
    } else if RANKED.contains(&cnvcall) {
        Ok(())
    } else {
        Err("group_untested_in_stage_a")
    }
}

/// Everything the ranking reads. Built in `d6_star_caller`.
pub struct RankingInput<'a> {
    pub cnvcall: &'a str,
    pub spacer_cn: Option<u32>,
    pub star_combinations: &'a StarCombinations,
    pub sites: &'a [ReadSite],
    pub provenance: Vec<ProvenanceRecord>,
    pub variants_called: &'a [String],
    pub uncallable: &'a [String],
    /// The final `Genotype` (after the resolver and the phase steps).
    pub call: Option<&'a str>,
    pub call_info: Option<&'a str>,
    pub raw_call: Option<&'a [String]>,
}

/// One scored profile: (composite ln-likelihood, table key index, raw diplotype index).
type ScoredProfile = (f64, u32, u32);

#[derive(Debug, Clone)]
struct Entry {
    score: f64,
    key: u32,
    raw: u32,
    label: String,
    formatted: bool,
    multiset: Vec<String>,
    rank: usize,
}

fn not_ranked(reason: &str) -> serde_json::Value {
    serde_json::json!({"Status": "not_ranked", "Reason": reason})
}

/// The `Candidate_ranking` value for one sample.
pub fn rank(input: RankingInput) -> serde_json::Value {
    if let Err(reason) = ranking_scope(input.cnvcall) {
        return not_ranked(reason);
    }
    match combination_table_for(input.cnvcall).map(|k| k.table(input.star_combinations)) {
        Some(t) if !t.is_empty() => rank_with_table(input, t),
        _ => not_ranked("no_table"),
    }
}

fn rank_with_table(input: RankingInput, table: &HashMap<String, Vec<String>>) -> serde_json::Value {
    let Some(call) = input.call.filter(|c| !c.is_empty()) else {
        return not_ranked("no_call");
    };
    let sites = input.sites;
    let defs = allele_defs(&input.star_combinations.dstar, sites);
    let ll = ll_table(sites);
    let mut keys: Vec<&String> = table.keys().collect();
    keys.sort();

    // Score every profile (each raw diplotype of each key).
    let scored: Vec<(Vec<ScoredProfile>, usize)> = keys
        .par_iter()
        .enumerate()
        .map(|(ki, key)| {
            let mut buf = vec![(0u8, 0u8); sites.len()];
            let mut out = Vec::new();
            let mut bad = 0usize;
            for (ri, raw) in table[*key].iter().enumerate() {
                match copies_of(raw, &defs) {
                    Some(c) if c.len() + 2 < MAX_N => {
                        dosage_into(&c, sites, &mut buf);
                        out.push((profile_score(&buf, &ll), ki as u32, ri as u32));
                    }
                    _ => bad += 1,
                }
            }
            (out, bad)
        })
        .collect();
    let unscorable: usize = scored.iter().map(|x| x.1).sum();
    let mut profs: Vec<ScoredProfile> = scored.into_iter().flat_map(|x| x.0).collect();
    if profs.is_empty() {
        return not_ranked("no_scorable_profile");
    }
    profs.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    let raw_of = |key: u32, raw: u32| -> &String { &table[keys[key as usize]][raw as usize] };

    // Walk in score order. Each distinct (multiset, label) is one entry.
    let call_sets: Vec<Vec<String>> = call.split(';').map(|g| copy_multiset(g.trim())).collect();
    let mut seen: HashSet<(Vec<String>, String)> = HashSet::new();
    let (mut listed, mut call_entry): (Vec<Entry>, Option<Entry>) = (Vec::new(), None);
    let (mut n_entries, mut tied_at_cutoff, mut current_rank) = (0usize, 0usize, 1usize);
    let mut last_score: Option<f64> = None;
    let mut i = 0;
    while i < profs.len() {
        let (score, key, _) = profs[i];
        let mut j = i;
        while j < profs.len() && profs[j].0 == score && profs[j].1 == key {
            j += 1;
        }
        if last_score != Some(score) {
            current_rank = n_entries + 1;
            last_score = Some(score);
        }
        let fifth = listed.get(RANKING_LIST_SIZE - 1).map(|e| e.score);
        if fifth.is_some_and(|f| score < f) && call_entry.is_some() {
            break;
        }
        let mut labels: Vec<(String, bool, u32)> = Vec::new();
        for p in &profs[i..j] {
            for (l, f) in entry_labels(raw_of(p.1, p.2), input.cnvcall, input.spacer_cn) {
                labels.push((l, f, p.2));
            }
        }
        labels.sort();
        labels.dedup_by(|a, b| a.0 == b.0);
        for (label, formatted, raw) in labels {
            let multiset = copy_multiset(&label);
            if !seen.insert((multiset.clone(), label.clone())) {
                continue;
            }
            n_entries += 1;
            let e = Entry { score, key, raw, label, formatted, multiset, rank: current_rank };
            let is_call = call_sets.contains(&e.multiset);
            if is_call && call_entry.is_none() {
                call_entry = Some(e.clone());
            }
            if listed.len() < RANKING_LIST_SIZE {
                listed.push(e);
            } else if listed[RANKING_LIST_SIZE - 1].score == score {
                tied_at_cutoff += 1;
            }
        }
        i = j;
    }

    // Decomposition against the rank-1 entry's profile.
    let ln10 = std::f64::consts::LN_10;
    let best = listed[0].clone();
    let d_best = dosage(&copies_of(raw_of(best.key, best.raw), &defs).unwrap(), sites);
    let observed = token_counts(input.variants_called);
    let entry_json = |e: &Entry| -> serde_json::Value {
        let copies = copies_of(raw_of(e.key, e.raw), &defs).unwrap();
        let d = dosage(&copies, sites);
        let rows: Vec<serde_json::Value> = sites
            .iter()
            .enumerate()
            .filter(|(i, s)| s.scored() && d[*i] != d_best[*i])
            .map(|(i, s)| {
                let ((kb, nb), (kc, nc)) = (d_best[i], d[i]);
                let r = (ll[i][nb as usize * MAX_N + kb as usize] - ll[i][nc as usize * MAX_N + kc as usize]) / ln10;
                serde_json::json!({"Site": s.name, "Alt_reads": s.alt, "Ref_reads": s.ref_,
                    "Copies_best": [kb, nb], "Copies_candidate": [kc, nc], "Site_log10_ratio": r})
            })
            .collect();
        let mut unmod = [0u32; 5];
        for c in &copies {
            for (u, v) in unmod.iter_mut().zip(c.unmodelled) {
                *u += v;
            }
        }
        let um: Vec<serde_json::Value> = UNMODELLED
            .iter()
            .zip(unmod)
            .filter_map(|(t, cand)| {
                let in_call = observed.get(*t).copied().unwrap_or(0);
                (in_call != cand).then(|| serde_json::json!({"Token": t, "Copies_call": in_call, "Copies_candidate": cand}))
            })
            .collect();
        serde_json::json!({"Rank": e.rank, "Genotype": e.label, "Copy_multiset": e.multiset,
            "Label_formatted": e.formatted, "Best_profile_key": keys[e.key as usize],
            "Delta_composite_log10_likelihood": (best.score - e.score) / ln10,
            "Discriminating_sites": rows, "Unmodelled_differences": um})
    };
    let mut candidates: Vec<serde_json::Value> = listed.iter().map(&entry_json).collect();
    if let Some(ce) = &call_entry {
        if !listed.iter().any(|e| e.key == ce.key && e.label == ce.label) {
            candidates.push(entry_json(ce));
        }
    }

    // Matching path, and records for retry / fuzzy differences from the observed profile.
    let matched = input
        .raw_call
        .and_then(|r| r.first())
        .and_then(|raw| profile_token_counts(raw, &input.star_combinations.dstar));
    let path = matching_path(input.cnvcall, input.call_info, input.variants_called, input.uncallable, matched.as_ref());
    let mut provenance = input.provenance;
    if let (Some(m), "retry_add_one_copy" | "fuzzy_match") = (&matched, path) {
        let mut toks: Vec<&String> = m.keys().chain(observed.keys()).collect();
        toks.sort();
        toks.dedup();
        for t in toks {
            let (o, k) = (observed.get(t).copied().unwrap_or(0), m.get(t).copied().unwrap_or(0));
            if o != k {
                let site = sites.iter().find(|s| &s.name == t);
                provenance.push(ProvenanceRecord {
                    site: t.clone(),
                    mechanism: path,
                    original_alt: site.map(|s| s.alt),
                    original_ref: site.map(|s| s.ref_),
                    per_site_copies: Some(o),
                    called_copies: k,
                });
            }
        }
    }
    let prov_json: Vec<serde_json::Value> = provenance
        .iter()
        .map(|p| serde_json::json!({"Site": p.site, "Mechanism": p.mechanism, "Original_alt": p.original_alt,
            "Original_ref": p.original_ref, "Per_site_copies": p.per_site_copies, "Called_copies": p.called_copies}))
        .collect();
    let masked: Vec<&str> = sites.iter().filter(|s| s.mask.is_some()).map(|s| s.name.as_str()).collect();
    serde_json::json!({
        "Status": "ranked",
        "Model": MODEL_ID,
        "Evidence_scale": "uncalibrated_composite",
        "Hypothesis_space": input.cnvcall,
        "Profiles_evaluated": profs.len(),
        "Profiles_unscorable": unscorable,
        "Matching_path": path,
        "Masked_sites": masked,
        "Unmodelled_tokens": UNMODELLED,
        "Call_rank": call_entry.as_ref().map(|e| e.rank),
        "Call_delta_composite_log10_likelihood": call_entry.as_ref().map(|e| (best.score - e.score) / ln10),
        "Entries_tied_at_cutoff": tied_at_cutoff,
        "Candidates": candidates,
        "Decision_provenance": prov_json,
    })
}

#[cfg(test)]
mod tests_rank {
    use super::*;
    use crate::caller::construct_star_table::{get_hap_table, TableKind};

    /// Sites with reads drawn from one profile: 40 reads, alt = round(40 k / n).
    fn sites_from(raw: &str, sc: &StarCombinations) -> Vec<ReadSite> {
        let blank: Vec<ReadSite> = site_specs(crate::data::TARGET_VARIANT, crate::data::TARGET_VARIANT_HOMO)
            .into_iter()
            .map(|(name, pos, class, homology)| ReadSite { name, pos, class, homology, alt: 0, ref_: 40, mask: None })
            .collect();
        let defs = allele_defs(&sc.dstar, &blank);
        let d = dosage(&copies_of(raw, &defs).unwrap(), &blank);
        blank
            .into_iter()
            .zip(d)
            .map(|(mut s, (k, n))| {
                let alt = if n == 0 { 0 } else { (40.0 * k as f64 / n as f64).round() as u32 };
                s.alt = alt;
                s.ref_ = 40 - alt;
                s
            })
            .collect()
    }

    fn tokens(raw: &str, sc: &StarCombinations) -> Vec<String> {
        let mut v: Vec<String> = profile_token_counts(raw, &sc.dstar)
            .unwrap()
            .into_iter()
            .flat_map(|(t, n)| std::iter::repeat_n(t, n as usize))
            .collect();
        v.sort();
        v
    }

    fn input<'a>(sc: &'a StarCombinations, sites: &'a [ReadSite], vc: &'a [String], call: &'a str, raw: &'a [String]) -> RankingInput<'a> {
        RankingInput { cnvcall: "cn2", spacer_cn: Some(2), star_combinations: sc, sites, provenance: Vec::new(),
            variants_called: vc, uncallable: &[], call: Some(call), call_info: Some("unique_match"), raw_call: Some(raw) }
    }

    fn check_sums(v: &serde_json::Value) {
        for c in v["Candidates"].as_array().unwrap() {
            let s: f64 = c["Discriminating_sites"].as_array().unwrap().iter().map(|r| r["Site_log10_ratio"].as_f64().unwrap()).sum();
            assert!((s - c["Delta_composite_log10_likelihood"].as_f64().unwrap()).abs() < 1e-9, "{c}");
        }
    }

    #[test]
    fn an_exact_cn2_call_ranks_first_and_its_equivalent_ties_with_it() {
        let sc = get_hap_table(crate::data::STAR_TABLE, &[TableKind::Dhap2]);
        let sites = sites_from("*1_*2", &sc);
        let vc = tokens("*1_*2", &sc);
        let raw = vec!["*1_*2".to_string()];
        let v = rank(input(&sc, &sites, &vc, "*1/*2", &raw));
        assert_eq!(v["Status"], "ranked", "{v}");
        assert_eq!(v["Call_rank"], 1);
        assert_eq!(v["Matching_path"], "exact");
        let c = v["Candidates"].as_array().unwrap();
        let by_ms = |ms: &[&str]| c.iter().filter(|e| e["Copy_multiset"] == serde_json::json!(ms)).count();
        assert_eq!(by_ms(&["*1", "*2"]), 1, "each genotype appears once");
        assert_eq!(by_ms(&["*34", "*39"]), 1);
        assert!(c.iter().filter(|e| e["Copy_multiset"] == serde_json::json!(["*34", "*39"])).all(|e| e["Rank"] == 1));
        check_sums(&v);
        // Lookalikes that differ only at unmodelled tokens (*27 g.42126938C>T, *127/*139 g.42127526C>T)
        // or only by phase (*34/*39) tie with the call: the scored read sites cannot tell them apart.
        assert_eq!(by_ms(&["*2", "*27"]), 1, "{v}");
        assert!(c.iter().all(|e| e["Rank"] == 1 && e["Discriminating_sites"].as_array().unwrap().is_empty()), "{v}");
    }

    #[test]
    fn ranking_is_deterministic() {
        let sc = get_hap_table(crate::data::STAR_TABLE, &[TableKind::Dhap2]);
        let sites = sites_from("*2_*4", &sc);
        let vc = tokens("*2_*4", &sc);
        let raw = vec!["*2_*4".to_string()];
        assert_eq!(rank(input(&sc, &sites, &vc, "*2/*4", &raw)), rank(input(&sc, &sites, &vc, "*2/*4", &raw)));
    }

    #[test]
    fn ties_share_a_rank_and_are_counted_at_the_cutoff() {
        let sc = get_hap_table(crate::data::STAR_TABLE, &[]);
        let mut sites = sites_from("*1_*1", &sc);
        sites.iter_mut().for_each(|s| s.mask = Some("strand_bias_all")); // nothing scored: all tie
        let table: HashMap<String, Vec<String>> = ["*1_*1", "*1_*2", "*1_*3", "*1_*4", "*1_*6", "*1_*9", "*1_*10"]
            .iter()
            .enumerate()
            .map(|(i, r)| (format!("k{i}"), vec![r.to_string()]))
            .collect();
        let vc: Vec<String> = Vec::new();
        let raw = vec!["*1_*1".to_string()];
        let v = rank_with_table(input(&sc, &sites, &vc, "*1/*1", &raw), &table);
        let c = v["Candidates"].as_array().unwrap();
        assert_eq!(c.len(), 5);
        assert!(c.iter().all(|e| e["Rank"] == 1));
        assert_eq!(v["Entries_tied_at_cutoff"], 2);
        assert_eq!(v["Call_rank"], 1);
    }

    #[test]
    fn a_call_below_the_list_is_appended() {
        let sc = get_hap_table(crate::data::STAR_TABLE, &[]);
        let sites = sites_from("*1_*1", &sc);
        let table: HashMap<String, Vec<String>> = ["*1_*1", "*1_*2", "*1_*3", "*1_*4", "*1_*6", "*1_*9", "*4_*10"]
            .iter()
            .enumerate()
            .map(|(i, r)| (format!("k{i}"), vec![r.to_string()]))
            .collect();
        let vc = tokens("*4_*10", &sc);
        let raw = vec!["*4_*10".to_string()];
        let v = rank_with_table(input(&sc, &sites, &vc, "*4/*10", &raw), &table);
        let c = v["Candidates"].as_array().unwrap();
        assert_eq!(c.len(), 6, "{v}");
        assert!(v["Call_rank"].as_u64().unwrap() > 5);
        assert_eq!(c[5]["Rank"], v["Call_rank"]);
        assert!(!c[5]["Discriminating_sites"].as_array().unwrap().is_empty(), "rows to sum: {v}");
        check_sums(&v);
    }

    #[test]
    fn not_ranked_reasons() {
        for (tag, reason) in [("star5_star5", "direct_genotype_tag"), ("star13", "hybrid_dosage_undefined"),
            ("exon9hyb_star68", "hybrid_dosage_undefined"), ("star5", "group_failed_transfer_test"),
            ("cn4", "group_untested_in_stage_a"), ("exon9hyb_star5", "group_untested_in_stage_a")] {
            assert_eq!(ranking_scope(tag), Err(reason), "{tag}");
        }
        for tag in ["cn2", "cn3", "exon9hyb", "exon9hyb_exon9hyb", "dup_exon9hyb"] {
            assert_eq!(ranking_scope(tag), Ok(()), "{tag}");
        }
        let sc = get_hap_table(crate::data::STAR_TABLE, &[TableKind::Dhap2]);
        let sites = sites_from("*1_*1", &sc);
        let mut i = input(&sc, &sites, &[], "", &[]);
        i.call = None;
        assert_eq!(rank(i)["Reason"], "no_call");
    }
}
```

- [ ] **Step 2: Run, expect pass**

Run: `cargo test --locked --lib candidate_ranking 2>&1 | tail -4`
Expected: `22 passed` (17 from Tasks 3-4 after the Task 4 fix, plus 5 here; measured at 5a7a5bb).

- [ ] **Step 3: Mutation checks**, restoring the code after each:
  - Delete `tied_at_cutoff += 1;` → `ties_share_a_rank...` FAILS.
  - Change `current_rank = n_entries + 1;` to `current_rank = 1;` →
    `a_call_below_the_list_is_appended` FAILS.

- [ ] **Step 4: Full suite and lint count**

Run:
- `cargo test --locked 2>&1 | grep "^test result" | grep -vc " 0 failed"` → `0`.
- `cargo clippy --all-targets --locked 2>&1 | grep -cE "^(warning|error)"` → at most `22`.

- [ ] **Step 5:** Remove the `#![allow(unused_imports, dead_code)]` line Task 3 added. Add no new allow attributes. Re-run Step 4.

- [ ] **Step 6: Commit**

```bash
git add cyrius-rs/src/caller/candidate_ranking.rs
git commit -m "feat(ranking): ranked entries, decomposition, call rank, matching path and output

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 6: Wiring, the flag and the JSON field

**Files:**
- Modify: `cyrius-rs/src/types.rs` (`FeatureFlags`, `D6Call`, `D6Call::no_call`)
- Modify: `cyrius-rs/src/main.rs` (CLI, `feature_flags`, `d6_star_caller`, JSON)
- Modify: `cyrius-rs/src/caller/match_star_allele.rs` (the two test `FeatureFlags` literals)

- [ ] **Step 1: Types.**
  - In `FeatureFlags`, after `kmer_validation`, add:
    ```rust
    /// Add `Candidate_ranking` to the JSON (output only; uncalibrated)
    pub candidate_ranking: bool,
    ```
  - In `D6Call`, after `uncallable_target_sites`, add:
    ```rust
    /// `Candidate_ranking` JSON (only with --candidate-ranking).
    pub candidate_ranking: Option<serde_json::Value>,
    ```
  - In `no_call`, add `candidate_ranking: None,`.
  - Add `candidate_ranking: false,` to both `FeatureFlags { ... }` literals in
    `match_star_allele.rs` tests.

- [ ] **Step 2: CLI.**
  - In `struct Cli`, after `kmer_validation: bool,`, add:
    ```rust
    #[arg(long = "candidate-ranking", help = "Add Candidate_ranking to the JSON: ranked genotypes in the call's CNV group, a site-wise decomposition and decision provenance (output only; uncalibrated)")]
    candidate_ranking: bool,
    ```
  - In `feature_flags`, add `candidate_ranking: cli.candidate_ranking,`.

- [ ] **Step 3: Snapshots in `d6_star_caller`:**
  - Directly after the `update_var42128936(...);` statement:
    `let var_alt_original = var_alt.clone();`.
  - Directly after the `if features.strand_bias_all { ... }` block:
    `let var_alt_after_sba = var_alt.clone();`.
  - Directly after `total_callset.extend(var42130655ins_a);`:
    `let callset_before_rules = total_callset.clone();`.

- [ ] **Step 4: The call.** Directly before the final `D6Call {` literal of `d6_star_caller`
  (the one with `total_cn: Some(total_cn),`):

```rust
    // Feature: candidate_ranking. Output only (spec 2026-11-05-candidate-ranking-design.md).
    let candidate_ranking = if features.candidate_ranking {
        use cyrius_rs::caller::candidate_ranking as cr;
        let variants_called: Vec<String> = star_called
            .variants_called
            .as_deref()
            .unwrap_or("")
            .split_whitespace()
            .map(str::to_string)
            .collect();
        let total_cn_sites = cyrius_rs::caller::call_variants::get_total_cn_per_site(
            &cnvtag_str, &params.var_db, &params.var_list,
        )
        .unwrap_or_default();
        let calls = cr::SiteCalls {
            var_list: &params.var_list,
            alt_original: &var_alt_original,
            alt_after_sba: &var_alt_after_sba,
            alt_final: &var_alt,
            ref_: &var_ref,
            homo_alt: &var_homo_alt,
            homo_ref: &var_homo_ref,
            total_cn_sites: &total_cn_sites,
            total_cn,
            decisions: &cn_call_var,
            homo_decisions: &cn_call_var_homo,
            min_read: cyrius_rs::caller::call_variants::min_read_for,
        };
        let (sites, provenance) =
            cr::sites_and_provenance(&calls, &callset_before_rules, &variants_called, &cnvtag_str);
        Some(cr::rank(cr::RankingInput {
            cnvcall: &cnvtag_str,
            spacer_cn: raw_cn_call.spacer_cn,
            star_combinations: &star_combinations,
            sites: &sites,
            provenance,
            variants_called: &variants_called,
            uncallable: &uncallable_target_sites,
            call: final_star_allele_call.as_deref(),
            call_info: star_called.call_info.as_deref(),
            raw_call: star_called.raw_call.as_deref(),
        }))
    } else {
        None
    };
```

  Then:
  - In that final `D6Call` literal, add `candidate_ranking,`.
  - Add `candidate_ranking: None,` to every other full `D6Call { ... }` literal in `main.rs`.
    Run `cargo build --locked 2>&1 | grep -A3 "missing field"` to list them. The literals that
    use `..D6Call::no_call(...)` need nothing.

- [ ] **Step 5: JSON.**
  - Change `let call_json = serde_json::json!({` to `let mut call_json = serde_json::json!({`.
  - Directly before `json_map.insert(sample_id.clone(), call_json);`, add:

```rust
    if let Some(cr) = &cyp2d6_call.candidate_ranking {
        call_json
            .as_object_mut()
            .expect("call_json is a JSON object")
            .insert("Candidate_ranking".to_string(), cr.clone());
    }
```

- [ ] **Step 6: Build and smoke-run**

```bash
cd /home/parlar_ai/dev/cyrius-rs/cyrius-rs && cargo build --release --locked 2>&1 | tail -1
T=/tmp/claude-1066/-home-parlar-ai-dev-cyrius-rs/b2ebb398-d73c-4f9f-8aa6-379fd1929e31/scratchpad/stageb_smoke
rm -rf $T && mkdir -p $T/on $T/off
./target/release/cyrius-rs -i ../benchmark/NA19920_pgx_regions.bam -o $T/on --id S --strand-bias-all --fuzzy-match --changepoint-hybrid --het-check --candidate-ranking 2> $T/on.log
./target/release/cyrius-rs -i ../benchmark/NA19920_pgx_regions.bam -o $T/off --id S --strand-bias-all --fuzzy-match --changepoint-hybrid --het-check 2> $T/off.log
python3 -c "import json;d=json.load(open('$T/on/S.json'))['S'];c=d['Candidate_ranking'];print(d['Genotype'],c['Status'],c.get('Call_rank'),c.get('Matching_path'),[x['Genotype'] for x in c.get('Candidates',[])])"
python3 -c "import json;print('Candidate_ranking' in json.load(open('$T/off/S.json'))['S'])"
```
Expected: one line with the genotype and `ranked` (or a `not_ranked` reason), then `False`.
This is a smoke run only. The checks are Task 8.

- [ ] **Step 7: Full suite and lint count** — as Task 5, Step 4.

- [ ] **Step 8: Commit**

```bash
git add cyrius-rs/src/types.rs cyrius-rs/src/main.rs cyrius-rs/src/caller/match_star_allele.rs
git commit -m "feat(ranking): --candidate-ranking flag adds Candidate_ranking to the JSON (output only)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 6b: Hybrid count fixed to the call's (amendment 2026-11-10)

Why: spec "Profiles and scores", **Hybrid count**. In the Task 6 smoke run, profiles with
more `exon9gc` copies than the call had `n_s = 0` at exon-9 sites, where the site term is 0
(the maximum), and topped 7 of 9 hybrid-family or exon-9 samples. Design choice this patches:
"Profiles = every key of the called group's table". It is the first patch to that choice.

**Files:**
- Modify: `cyrius-rs/src/caller/candidate_ranking.rs` (`rank_with_table`, one struct, one test)

**Interfaces:**
- Consumes: `count_of(tokens: &[String], t: &str) -> u32`, `copies_of`, `AlleleDef.gc`.
- Produces: JSON field `Profiles_other_exon9gc_count` (int), and the `not_ranked` reason
  `no_profile_with_called_exon9gc_count`.

- [ ] **Step 1: Write the failing test** (in `mod tests_rank`):

```rust
    #[test]
    fn profiles_with_another_exon9gc_count_are_not_scored() {
        let sc = get_hap_table(crate::data::STAR_TABLE, &[]);
        let sites = sites_from("*10_*36_*36", &sc);
        let vc = tokens("*10_*36_*36", &sc);
        assert_eq!(count_of(&vc, "exon9gc"), 2);
        let raw = vec!["*10_*36_*36".to_string()];
        let table: HashMap<String, Vec<String>> = ["*10_*36_*36", "*36_*36_*36", "*36_*36_*83"]
            .iter()
            .enumerate()
            .map(|(i, r)| (format!("k{i}"), vec![r.to_string()]))
            .collect();
        let mut i = input(&sc, &sites, &vc, "*36/*36+*10", &raw);
        i.cnvcall = "exon9hyb";
        let v = rank_with_table(i, &table);
        assert_eq!(v["Profiles_evaluated"], 1, "{v}");
        assert_eq!(v["Profiles_other_exon9gc_count"], 2, "{v}");
        assert_eq!(v["Call_rank"], 1, "{v}");
        let only_other: HashMap<String, Vec<String>> = [("k0".to_string(), vec!["*36_*36_*83".to_string()])].into();
        let mut i = input(&sc, &sites, &vc, "*36/*36+*10", &raw);
        i.cnvcall = "exon9hyb";
        assert_eq!(rank_with_table(i, &only_other)["Reason"], "no_profile_with_called_exon9gc_count");
    }
```

- [ ] **Step 2: Run it, expect FAIL**

Run: `cargo test --locked --lib profiles_with_another_exon9gc_count_are_not_scored 2>&1 | tail -5`
Expected: FAIL on `Profiles_evaluated` (3, not 1).

- [ ] **Step 3: Implement.** Directly above `fn rank_with_table`, add:

```rust
/// One table key's scored profiles, and how many of its raws were left out.
struct KeyScores {
    profiles: Vec<ScoredProfile>,
    unscorable: usize,
    other_exon9gc: usize,
}
```

In `rank_with_table`, replace the block from `// Score every profile` through the
`if profs.is_empty() { ... }` check with:

```rust
    // Score every profile (each raw diplotype of each key) that has the called number of
    // exon9gc copies (spec "Profiles and scores", Hybrid count).
    let called_gc = count_of(input.variants_called, "exon9gc") as usize;
    let scored: Vec<KeyScores> = keys
        .par_iter()
        .enumerate()
        .map(|(ki, key)| {
            let mut buf = vec![(0u8, 0u8); sites.len()];
            let mut ks = KeyScores { profiles: Vec::new(), unscorable: 0, other_exon9gc: 0 };
            for (ri, raw) in table[*key].iter().enumerate() {
                match copies_of(raw, &defs) {
                    Some(c) if c.len() + 2 < MAX_N => {
                        if c.iter().filter(|d| d.gc).count() != called_gc {
                            ks.other_exon9gc += 1;
                        } else {
                            dosage_into(&c, sites, &mut buf);
                            ks.profiles.push((profile_score(&buf, &ll), ki as u32, ri as u32));
                        }
                    }
                    _ => ks.unscorable += 1,
                }
            }
            ks
        })
        .collect();
    let unscorable: usize = scored.iter().map(|x| x.unscorable).sum();
    let other_exon9gc: usize = scored.iter().map(|x| x.other_exon9gc).sum();
    let mut profs: Vec<ScoredProfile> = scored.into_iter().flat_map(|x| x.profiles).collect();
    if profs.is_empty() {
        return not_ranked(if other_exon9gc > 0 { "no_profile_with_called_exon9gc_count" } else { "no_scorable_profile" });
    }
```

In the final `serde_json::json!({ ... })`, directly after `"Profiles_unscorable": unscorable,`, add:

```rust
        "Profiles_other_exon9gc_count": other_exon9gc,
```

- [ ] **Step 4: Run, expect PASS.** `cargo test --locked --lib candidate_ranking 2>&1 | tail -4` → 24 passed.

- [ ] **Step 5: Mutation check.** Replace `!= called_gc` with `!= usize::MAX` (no profile excluded).
  Expect `profiles_with_another_exon9gc_count_are_not_scored` to FAIL. Restore.

- [ ] **Step 6: Full suite and lint count** (as Task 5, Step 4): 0 failing lines; clippy at most 22,
  with no new allow attributes.

- [ ] **Step 7: Re-run the 9 samples (descriptive).** Build release. For NA18526, HG00463,
  NA18563, NA18564, NA18565, NA18572, HG02373, NA18545 and HG00421 (benchmark BAMs, flags
  `--strand-bias-all --fuzzy-match --changepoint-hybrid --het-check --candidate-ranking`), report:
  Genotype, CNV_group, Status, Call_rank, Call_delta, the top entry's genotype and its
  `Best_profile_key` exon9gc count, the call's exon9gc count, `Profiles_evaluated`,
  `Profiles_other_exon9gc_count`, and wall time. Before (15c57b4): call ranks 2, 6, 38, 1, 2, 5,
  1, 6, 1. This step is descriptive: the fix is judged by Steps 2-5, not by these ranks.

- [ ] **Step 8: Commit**

```bash
git add cyrius-rs/src/caller/candidate_ranking.rs
git commit -m "fix(ranking): score only profiles with the call's exon9gc copy count

A profile with more exon9gc copies had n_s = 0 at exon-9 sites and a free
site term of 0 (spec amendment 2026-11-10, Hybrid count). Patches design
choice: hypothesis space = every key of the called group's table (patch 1).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 7: Spec and README

**Files:**
- Modify: `docs/superpowers/specs/2026-11-05-candidate-ranking-design.md`
- Modify: `README.md` (flag table and a short `Candidate_ranking` paragraph)

- [ ] **Step 1: Spec edits** (Stage B section):
  - **Decision provenance:** replace "The functions involved report the mechanism alongside
    each value…" with this: provenance is found by replaying `call_reg1_cn`,
    `process_raw_call_denovo` and the homology cap per site, using the caller's own functions
    plus `cn_posterior`. A replay that disagrees with the value the call used reads
    `unexplained`, and B4 requires none. Rule-decided values are recorded whether or not the
    plain posterior would have agreed.
  - **Not-ranked reasons:** add `group_untested_in_stage_a` and `no_scorable_profile`.
    Replace "Groups where it does not fit are also `not_ranked`" with the ranked-tag list and
    the star5 outcome from this plan's Global Constraints.
  - **Matching_path:** derived from the call's first raw candidate against `Variants_called`.
    `retry_add_one_copy` and `fuzzy_match` add one provenance record per differing token.
  - **Model:** M4 (exp-H20 `0a4b020`), constants in `count_model.rs`. Fragment reuse median
    96.1%.
- [ ] **Step 2: README.** Add `--candidate-ranking` to the flag list: "Adds
  `Candidate_ranking` (output only, uncalibrated). See the spec." Add one paragraph under the
  output fields: what `Candidate_ranking` holds, that its numbers are an uncalibrated
  composite score and not probabilities, and which groups are ranked.
- [ ] **Step 3: Commit**

```bash
git add docs/superpowers/specs/2026-11-05-candidate-ranking-design.md README.md
git commit -m "docs: candidate ranking Stage B as built (provenance by replay, ranked tags, reasons)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 8: Run the locked checks and record the result

**Files:**
- Create (research): `~/cyrius-research/candidate-ranking/run_one_b.sh`,
  `~/cyrius-research/candidate-ranking/run_b.sh`,
  `~/cyrius-research/candidate-ranking/check_b.py`
- Modify (repo): this plan, appending "Result".

- [ ] **Step 1: Scripts.** `run_one_b.sh`:

```bash
#!/usr/bin/env bash
# run_one_b.sh BIN OUTDIR ID BAM [EXTRA...]
set -u
B=$1; O=$2; S=$3; BAM=$4; shift 4
"$B" -i "$BAM" -o "$O" --id "$S" --strand-bias-all --fuzzy-match --changepoint-hybrid --het-check "$@" 2> "$O/$S.log"
echo "$S $?" >> "$O/exit.txt"
```

`run_b.sh`:

```bash
#!/usr/bin/env bash
# Stage B checks: base (baf3459), noflag, flag1, flag2 on 160 files.
set -euo pipefail
C=/home/parlar_ai/cyrius-research/candidate-ranking
REPO=/home/parlar_ai/dev/cyrius-rs
U=/home/parlar_ai/cyrius-research/unknown-sites
BASE=/home/parlar_ai/cyrius-research/exp-H20/frozen/cyrius-rs
BIN=$REPO/cyrius-rs/target/release/cyrius-rs
[ -z "$(git -C $REPO status --porcelain -- cyrius-rs)" ] || { echo "crate dir dirty: refusing"; exit 1; }
(cd /home/parlar_ai/cyrius-research/exp-H20 && sha256sum -c <(head -1 data/checksums.txt))
(cd $REPO/cyrius-rs && cargo build --release --locked 2>&1 | tail -1)
{ git -C $REPO rev-parse --short=7 HEAD; sha256sum $BIN; } > $C/commit.txt
: > $C/runtime.txt
for cfg in base noflag flag1 flag2; do
  case $cfg in base) B=$BASE; X="";; noflag) B=$BIN; X="";; *) B=$BIN; X="--candidate-ranking";; esac
  rm -rf $C/$cfg; mkdir -p $C/$cfg/bench $C/$cfg/p25 $C/$cfg/p10 $C/$cfg/dmg
  {
    grep -v '^#' $REPO/benchmark/getrm_cyp2d6_truth.tsv | grep -v '^$' | cut -f1 \
      | while read -r s; do echo "$B $C/$cfg/bench $s $REPO/benchmark/${s}_pgx_regions.bam $X"; done
    for p in p25 p10; do
      while read -r s; do echo "$B $C/$cfg/$p $s $U/thinned/$p/${s}_pgx_regions.bam $X"; done < $U/thinned/samples.txt
    done
    echo "$B $C/$cfg/dmg dmg /home/parlar_ai/cyrius-research/review-fixes/na19920_missing_splice/NA19920_missing_splice.bam $X"
  } > $C/$cfg/list.txt
  start=$(date +%s)
  xargs -P 8 -L 1 bash $C/run_one_b.sh < $C/$cfg/list.txt
  echo "$cfg $(( $(date +%s) - start ))s $(cat $C/$cfg/*/exit.txt | awk '$2!=0' | wc -l) nonzero-exit" | tee -a $C/runtime.txt
done
```

`check_b.py`:

```python
#!/usr/bin/env python3
"""Stage B checks B1-B4 and vacuity guards (plan 2026-11-05-candidate-ranking-stage-b.md).
Exit 1 on any failure."""
import collections, glob, json, os, sys
C = "/home/parlar_ai/cyrius-research/candidate-ranking"
SETS = ("bench", "p25", "p10", "dmg")

def load(cfg):
    out = {}
    for st in SETS:
        for f in sorted(glob.glob(f"{C}/{cfg}/{st}/*.json")):
            s = os.path.basename(f)[:-5]
            out[(st, s)] = json.load(open(f))[s]
    return out

base, nof, f1, f2 = (load(c) for c in ("base", "noflag", "flag1", "flag2"))
fails = []
counts = [len(x) for x in (base, nof, f1, f2)]
if counts != [160] * 4:
    fails.append(f"file counts {counts}, want 160 each")
if not base or any("Genotype" not in v for v in base.values()):
    fails.append("Genotype missing in base: refusing a vacuous B3")

def diff(a, b, ignore):
    return [(k, f) for k in sorted(a) for f in sorted(set(a[k]) | set(b.get(k, {})))
            if f not in ignore and a[k].get(f) != b.get(k, {}).get(f)]

d = diff(base, nof, {"Tool_commit", "Flags"})
if d: fails.append(f"B3a: {len(d)} differences without the flag, e.g. {d[:3]}")
if any("Candidate_ranking" in v for v in nof.values()): fails.append("B3a: Candidate_ranking present without the flag")
d = diff(base, f1, {"Tool_commit", "Flags", "Candidate_ranking"})
if d: fails.append(f"B3b: {len(d)} existing-field differences with the flag, e.g. {d[:3]}")
d = [k for k in sorted(f1) if f1[k].get("Candidate_ranking") != f2.get(k, {}).get("Candidate_ranking")]
if d: fails.append(f"B2: {len(d)} files differ between the two flagged runs, e.g. {d[:3]}")

status, reasons, paths, ranks = (collections.Counter() for _ in range(4))
rows_nonempty, n_rec, worst = 0, 0, 0.0
for k, v in sorted(f1.items()):
    cr = v.get("Candidate_ranking")
    if cr is None:
        fails.append(f"no Candidate_ranking in flagged file {k}")
        continue
    status[cr["Status"]] += 1
    if cr["Status"] != "ranked":
        reasons[cr["Reason"]] += 1
        continue
    paths[cr["Matching_path"]] += 1
    ranks[str(cr["Call_rank"])] += 1
    if cr["Matching_path"] == "unexplained":
        fails.append(f"B4: unexplained matching path at {k}")
    for c in cr["Candidates"]:
        s = sum(r["Site_log10_ratio"] for r in c["Discriminating_sites"])
        worst = max(worst, abs(s - c["Delta_composite_log10_likelihood"]))
        rows_nonempty += bool(c["Discriminating_sites"])
    called = v.get("Variants_called") or []
    for r in cr["Decision_provenance"]:
        n_rec += 1
        if r["Mechanism"] == "unexplained":
            fails.append(f"B4: unexplained provenance at {k} {r['Site']}")
        if r["Mechanism"] in ("update_variants", "star68_step") and r["Called_copies"] != called.count(r["Site"]):
            fails.append(f"B4: Called_copies mismatch at {k} {r['Site']}")
if worst > 1e-9: fails.append(f"B1: worst |sum(rows) - delta| = {worst}")
if status["ranked"] < 80: fails.append(f"vacuous: {status['ranked']} ranked files (< 80)")
if rows_nonempty < 10: fails.append(f"vacuous: {rows_nonempty} candidates with discriminating rows (< 10)")
if n_rec == 0: fails.append("vacuous: no provenance record")
report = {"fails": fails, "status": status, "not_ranked_reasons": reasons, "matching_paths": paths,
          "call_ranks": ranks, "b1_worst": worst, "candidates_with_rows": rows_nonempty,
          "provenance_records": n_rec, "files": counts}
json.dump(report, open(f"{C}/check_b.json", "w"), indent=1)
print(json.dumps(report, indent=1))
sys.exit(1 if fails else 0)
```

- [ ] **Step 2: Run**

```bash
mkdir -p /home/parlar_ai/cyrius-research/candidate-ranking && cd /home/parlar_ai/cyrius-research/candidate-ranking
bash run_b.sh && python3 check_b.py; echo "check exit=$?"
```
Expected: `check exit=0`, with `fails` empty.

- [ ] **Step 3: Result commit.** Append a "Result" section to this plan, with the commit
  measured, `runtime.txt`, and `check_b.json` verbatim. Commit it as:
  `result: candidate ranking Stage B checks — B1-B4 <pass|fail>`.
  On any failure, record it, stop and report. No threshold is changed.

# Diplotype Resolver + Scorer Consolidation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the probabilistic `diplotype_caller` the single tie-breaking scorer that resolves `;`-ambiguous diplotype calls, emit calibrated confidence + honest ambiguity, migrate the 6 hardcoded special-cases reproduce-or-promote, and delete four logging-only scorer modules.

**Architecture:** A thin resolver layer sits between the existing cnvtag→`match_star` cascade and the output. When the cascade yields a `;`-joined candidate set, the resolver scores each candidate's *allele multiset* with `diplotype_caller` and either picks confidently (top posterior clears a margin) or reports the near-tied set as ambiguous. Unambiguous samples bypass the resolver entirely, protecting the 53/53 benchmark.

**Tech Stack:** Rust, rust-htslib, rayon. Spec: `docs/superpowers/specs/2026-06-04-diplotype-resolver-design.md`.

**Working directory:** `/home/parlar_ai/dev/cyrius-rs/cyrius-rs` (the crate). Branch: `feat/diplotype-resolver`.

**Build/test/benchmark commands (used throughout):**
- Build: `cargo build --release`
- Unit tests for a module: `cargo test --release <module>::`
- Full benchmark (primary gate): from repo root `/home/parlar_ai/dev/cyrius-rs`, run
  `bash benchmark/run_benchmark.sh --strand-bias-all --fuzzy-match --changepoint-hybrid --het-check` — expect `Correct: 53 / 53`.

**Key facts (verified against code):**
- `diplotype_caller` scores a diplotype as a *multiset of base alleles* (pair of distinct alleles + copy split summing to total CN), blind to arrangement. Different multisets are separable by read likelihood; same-multiset-different-arrangement are not (→ correct ambiguity signal).
- The cascade entry is `match_star(...)` → `star_called.clean_call: Option<String>` (possibly `;`-joined). Consumed in `d6_star_caller` at `main.rs:777-827`.
- `phase_disambiguate` already hooks `;`-calls for `cn2 more_than_one_match` at `main.rs:801-827`; it stays untouched. The resolver runs *after* it.
- The existing `Confidence` TSV field is a computed score from `caller::confidence` — DO NOT clobber it. The resolver's posterior goes to JSON as a new field; the genotype-string decision flows into the existing `phasing_ambiguous`/`long_read_recommended` logic at `main.rs:1200-1226`.
- `StarCombinations.dstar: HashMap<String,String>` maps allele name (e.g. `"*36"`, `"*4.013"`, `"*1"`) → underscore-joined variant key. `build_d6_allele_defs(&dstar, &sites)` produces `AlleleDef{name, alt_sites}` whose names equal dstar keys.

---

## Task 1: Add `score_candidates` (explicit-multiset scoring) to `diplotype_caller`

Pure addition. `call_diplotype` is untouched; the scorer still only logs after this task.

**Files:**
- Modify: `src/caller/diplotype_caller.rs` (add `CandidateScore`, `parse_multiset`, `split_copy_suffix`, `score_candidates`; add `#[cfg(test)]` tests)

- [ ] **Step 1: Write failing tests for the multiset parser**

Add to the existing `#[cfg(test)] mod tests` block at the bottom of `src/caller/diplotype_caller.rs`:

```rust
#[test]
fn test_split_copy_suffix() {
    assert_eq!(split_copy_suffix("*36"), ("*36".to_string(), 1));
    assert_eq!(split_copy_suffix("*36x2"), ("*36".to_string(), 2));
    assert_eq!(split_copy_suffix("*4.013"), ("*4.013".to_string(), 1));
    // 'x' that is not a copy suffix must not be stripped
    assert_eq!(split_copy_suffix("*1x"), ("*1x".to_string(), 1));
}

#[test]
fn test_parse_multiset() {
    // *36/*36+*10  -> {*36:2, *10:1}
    let m = parse_multiset("*36/*36+*10");
    assert_eq!(m.iter().find(|(n, _)| n == "*36").map(|(_, c)| *c), Some(2));
    assert_eq!(m.iter().find(|(n, _)| n == "*10").map(|(_, c)| *c), Some(1));
    // *10/*36x2 -> {*10:1, *36:2}  (same multiset as above)
    let m2 = parse_multiset("*10/*36x2");
    assert_eq!(m2.iter().find(|(n, _)| n == "*36").map(|(_, c)| *c), Some(2));
    assert_eq!(m2.iter().find(|(n, _)| n == "*10").map(|(_, c)| *c), Some(1));
}
```

- [ ] **Step 2: Run tests, verify they fail to compile (functions undefined)**

Run: `cargo test --release diplotype_caller::tests::test_parse_multiset 2>&1 | tail -5`
Expected: compile error `cannot find function split_copy_suffix` / `parse_multiset`.

- [ ] **Step 3: Implement the parser + the public types/function**

Add near the top of `src/caller/diplotype_caller.rs` (after the existing `use` lines), the public result type:

```rust
/// One scored candidate diplotype (for the resolver).
#[derive(Debug, Clone)]
pub struct CandidateScore {
    pub diplotype: String,
    pub log10_likelihood: f64,
    /// Posterior normalised over the candidate set only (0 for unscorable).
    pub posterior: f64,
}
```

Add these functions in the "Diplotype enumeration + scoring" section (near `score_diplotypes`):

```rust
/// Split "*36x2" -> ("*36", 2); "*36" -> ("*36", 1). Only strips a trailing
/// 'x' followed by all ASCII digits (so "*4.013" and "*1x" are left intact).
fn split_copy_suffix(s: &str) -> (String, usize) {
    if let Some(pos) = s.rfind('x') {
        let num = &s[pos + 1..];
        if !num.is_empty() && num.chars().all(|c| c.is_ascii_digit()) {
            if let Ok(n) = num.parse::<usize>() {
                return (s[..pos].to_string(), n.max(1));
            }
        }
    }
    (s.to_string(), 1)
}

/// Parse a diplotype string into a base-allele -> copies multiset.
/// e.g. "*36/*36+*10" -> [("*10",1),("*36",2)] (order by allele name).
fn parse_multiset(diplotype: &str) -> Vec<(String, usize)> {
    use std::collections::BTreeMap;
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for hap in diplotype.split('/') {
        for comp in hap.split('+') {
            let comp = comp.trim();
            if comp.is_empty() {
                continue;
            }
            let (name, mult) = split_copy_suffix(comp);
            *counts.entry(name).or_insert(0) += mult;
        }
    }
    counts.into_iter().collect()
}
```

- [ ] **Step 4: Run parser tests, verify pass**

Run: `cargo test --release diplotype_caller::tests::test_parse_multiset diplotype_caller::tests::test_split_copy_suffix 2>&1 | tail -5`
Expected: 2 passed.

- [ ] **Step 5: Write a failing test for `score_candidates` ranking using synthetic data**

This test builds two candidate strings over alleles known to exist in a tiny dstar, with synthetic reads favoring one multiset. Because constructing real `bam::IndexedReader` reads is heavy, this test exercises the *scoring core* via a thin internal helper. Add the helper signature test first:

```rust
#[test]
fn test_score_explicit_multiset_prefers_matching_reads() {
    // Two sites; allele A has ALT at site0, allele B has ALT at site1.
    let sites = vec![
        test_site("g.1A>T", 1, b'A', b'T'),
        test_site("g.2C>G", 2, b'C', b'G'),
    ];
    let allele_a = AlleleDef { name: "*A".into(), alt_sites: [0usize].into_iter().collect() };
    let allele_b = AlleleDef { name: "*B".into(), alt_sites: [1usize].into_iter().collect() };
    let alleles = vec![allele_a, allele_b];
    // 10 reads all showing ALT at site0 (supports *A), REF at site1.
    let reads: Vec<ReadObs> = (0..10)
        .map(|i| ReadObs { name: format!("r{i}"), observations: vec![(0, b'T', 40), (1, b'C', 40)] })
        .collect();
    let lut = QualLookup::new();
    // Homozygous *A vs homozygous *B at cn=2 — *A must score higher.
    let s_a = score_explicit_multiset(&reads, &alleles, &sites, 2, &[(0usize, 2usize)], &lut);
    let s_b = score_explicit_multiset(&reads, &alleles, &sites, 2, &[(1usize, 2usize)], &lut);
    assert!(s_a > s_b, "reads favor *A: {s_a} !> {s_b}");
}
```

Add this test-only helper to the test module (constructs a single-base SNP `UnifiedSite`):

```rust
fn test_site(name: &str, pos: i64, d6: u8, d7: u8) -> UnifiedSite {
    UnifiedSite {
        d6_pos: pos, d7_pos: Some(pos + 1000), d6_ref_base: d6, d7_base: Some(d7),
        d6_alt_base: Some(d7), name: name.to_string(), is_deletion: false,
    }
}
```

- [ ] **Step 6: Run it, verify failure (helper undefined)**

Run: `cargo test --release diplotype_caller::tests::test_score_explicit_multiset 2>&1 | tail -5`
Expected: compile error `cannot find function score_explicit_multiset`.

- [ ] **Step 7: Implement `score_explicit_multiset` and `score_candidates`**

Add to `src/caller/diplotype_caller.rs`:

```rust
/// Score one explicit diplotype given as (allele_index, copies) pairs summing
/// to total_cn. General k-allele mixture: one D7 copy + the specified D6
/// alleles, each weighted by copies/(cn+1). Returns log10-likelihood.
fn score_explicit_multiset(
    reads: &[ReadObs],
    d6_alleles: &[AlleleDef],
    sites: &[UnifiedSite],
    total_cn: u32,
    assignment: &[(usize, usize)],
    lut: &QualLookup,
) -> f64 {
    let ln10 = std::f64::consts::LN_10;
    let denom = (total_cn as usize + 1) as f64;
    let ln_w_d7 = (1.0 / denom).ln();
    let mut total = 0.0f64;
    for obs in reads {
        let s_d7 = ln_w_d7 + score_read_d7(obs, sites, lut) * ln10;
        let mut terms: Vec<f64> = Vec::with_capacity(assignment.len() + 1);
        terms.push(s_d7);
        for &(idx, copies) in assignment {
            if copies == 0 {
                continue;
            }
            let ln_w = ((copies as f64) / denom).ln();
            terms.push(ln_w + score_read_d6(obs, &d6_alleles[idx], sites, lut) * ln10);
        }
        let max_s = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        if !max_s.is_finite() {
            continue;
        }
        let lse = max_s + terms.iter().map(|t| (t - max_s).exp()).sum::<f64>().ln();
        total += lse;
    }
    total / ln10
}

/// Map a parsed multiset to (allele_index, copies) pairs, or None if ANY
/// allele is unknown to the scorer or the copy total != cn. Pure + testable.
fn build_assignment(
    multiset: &[(String, usize)],
    name_to_idx: &std::collections::HashMap<&str, usize>,
    cn: usize,
) -> Option<Vec<(usize, usize)>> {
    let mut assignment: Vec<(usize, usize)> = Vec::new();
    let mut total_copies = 0usize;
    for (name, copies) in multiset {
        let idx = *name_to_idx.get(name.as_str())?;
        assignment.push((idx, *copies));
        total_copies += copies;
    }
    if total_copies != cn {
        return None;
    }
    Some(assignment)
}

/// Score a set of candidate diplotype strings against per-read evidence.
/// Returns one CandidateScore per input candidate, with posteriors normalised
/// across the set. Returns an EMPTY vec to signal "cannot score" — no reads,
/// no sites, OR (FAIL-CLOSED) ANY candidate that is unparseable, names an
/// unknown allele, or whose copy total != cn — so the caller falls back to the
/// cascade's own choice. Scoring only the representable subset is forbidden:
/// it could make a confident wrong pick by dropping one side of the ambiguity.
pub fn score_candidates(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    var_list: &[String],
    snp_file: &str,
    star_combinations: &StarCombinations,
    total_cn: u32,
    candidates: &[String],
) -> Vec<CandidateScore> {
    if candidates.is_empty() {
        return Vec::new();
    }
    let lut = QualLookup::new();
    let sites = build_unified_sites(var_list, snp_file);
    if sites.is_empty() {
        return Vec::new();
    }
    let d6_alleles = build_d6_allele_defs(&star_combinations.dstar, &sites);
    let mut name_to_idx: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
    for (i, a) in d6_alleles.iter().enumerate() {
        name_to_idx.insert(a.name.as_str(), i);
    }
    let all_reads = collect_reads_dual_region(
        reader, nchr, &sites, D6_START, D6_END, D7_START, D7_END,
    );
    let (reads, _nf) = filter_mismapped(all_reads, &d6_alleles, &sites, -1.5, &lut);
    if reads.is_empty() {
        return Vec::new();
    }

    let cn = total_cn as usize;
    // First pass: validate ALL candidates (fail-closed). Bail out entirely if
    // any candidate can't be represented.
    let mut assignments: Vec<Vec<(usize, usize)>> = Vec::with_capacity(candidates.len());
    for cand in candidates {
        let multiset = parse_multiset(cand);
        match build_assignment(&multiset, &name_to_idx, cn) {
            Some(a) => assignments.push(a),
            None => return Vec::new(), // fail-closed -> caller falls back
        }
    }

    // Second pass: score each validated candidate.
    let mut out: Vec<CandidateScore> = Vec::with_capacity(candidates.len());
    for (cand, assignment) in candidates.iter().zip(assignments.iter()) {
        let ll = score_explicit_multiset(&reads, &d6_alleles, &sites, total_cn, assignment, &lut);
        out.push(CandidateScore { diplotype: cand.clone(), log10_likelihood: ll, posterior: 0.0 });
    }

    // Normalise posteriors across the candidate set.
    let max_ll = out
        .iter()
        .map(|c| c.log10_likelihood)
        .fold(f64::NEG_INFINITY, f64::max);
    if !max_ll.is_finite() {
        return Vec::new();
    }
    let sum: f64 = out
        .iter()
        .map(|c| 10.0_f64.powf(c.log10_likelihood - max_ll))
        .sum();
    let log10_sum = max_ll + sum.log10();
    for c in &mut out {
        c.posterior = 10.0_f64.powf(c.log10_likelihood - log10_sum);
    }
    out
}
```

Note: `ReadObs.name` is currently `#[allow(dead_code)]`; the test constructs it directly, which is fine.

- [ ] **Step 7b: Add a fail-closed unit test for `build_assignment`**

Add to the test module (verifies the FAIL-CLOSED rule: a mixed set with one
unknown allele yields None, so `score_candidates` would return empty → fallback):

```rust
#[test]
fn test_build_assignment_fail_closed() {
    let mut idx: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
    idx.insert("*10", 0);
    idx.insert("*36", 1);
    // Known alleles, copies sum to cn=3 -> Some.
    let known = parse_multiset("*36/*36+*10");
    assert!(build_assignment(&known, &idx, 3).is_some());
    // One unknown allele (*4.013) -> None (fail-closed), even mixed with known.
    let mixed = parse_multiset("*4.013/*10+*36");
    assert!(build_assignment(&mixed, &idx, 3).is_none());
    // Copy total mismatch (cn=2 but multiset has 3 copies) -> None.
    assert!(build_assignment(&known, &idx, 2).is_none());
}
```

- [ ] **Step 8: Run the scoring test + full module tests, verify pass**

Run: `cargo test --release diplotype_caller:: 2>&1 | tail -8`
Expected: all pass (including the new `test_score_explicit_multiset_prefers_matching_reads`).

- [ ] **Step 9: Build the whole crate to confirm no breakage**

Run: `cargo build --release 2>&1 | tail -3`
Expected: `Finished`, no errors.

- [ ] **Step 10: Commit**

```bash
git add src/caller/diplotype_caller.rs
git commit -m "feat(diplotype_caller): add score_candidates explicit-multiset scoring"
```

---

## Task 2: Create the resolver decision rule (`resolver.rs`)

Pure function over `CandidateScore`; no BAM, fully unit-testable.

**Files:**
- Create: `src/caller/resolver.rs`
- Modify: `src/caller/mod.rs` (add `pub mod resolver;`)

- [ ] **Step 1: Write failing tests**

Create `src/caller/resolver.rs` with only the test module first:

```rust
use crate::caller::diplotype_caller::CandidateScore;

#[cfg(test)]
mod tests {
    use super::*;

    fn cs(name: &str, post: f64) -> CandidateScore {
        CandidateScore { diplotype: name.into(), log10_likelihood: post.ln(), posterior: post }
    }

    #[test]
    fn single_candidate_is_confident() {
        let out = resolve(&[cs("*1/*2", 1.0)]).unwrap();
        assert_eq!(out.genotype, "*1/*2");
        assert!(!out.ambiguous);
    }

    #[test]
    fn clear_winner_picks_top() {
        let out = resolve(&[cs("*1/*2", 0.98), cs("*1/*3", 0.02)]).unwrap();
        assert_eq!(out.genotype, "*1/*2");
        assert!(!out.ambiguous);
    }

    #[test]
    fn near_tie_reports_ambiguous_set() {
        let out = resolve(&[cs("*36/*36+*10", 0.5), cs("*10/*36x2", 0.5)]).unwrap();
        assert!(out.ambiguous);
        assert!(out.genotype.contains(';'));
        assert!(out.genotype.contains("*36/*36+*10") && out.genotype.contains("*10/*36x2"));
    }

    #[test]
    fn all_unscorable_returns_none() {
        let mut a = cs("*1/*2", 0.0);
        a.log10_likelihood = f64::NEG_INFINITY;
        assert!(resolve(&[a]).is_none());
    }
}
```

- [ ] **Step 2: Run tests, verify fail (resolve undefined)**

Run: `cargo test --release resolver:: 2>&1 | tail -5`
Expected: compile error `cannot find function resolve`.

- [ ] **Step 3: Implement the resolver**

Prepend to `src/caller/resolver.rs` (above the test module):

```rust
//! Likelihood-based tie-breaker for ;-ambiguous diplotype candidate sets.
//! Picks a single call when one candidate's posterior clearly dominates,
//! otherwise reports the near-tied set as honest ambiguity.

use crate::caller::diplotype_caller::CandidateScore;

/// Minimum posterior gap between the top two candidates required for a
/// confident single call. Calibrated on the GeT-RM benchmark ONLY and set
/// conservatively: prefer declaring ambiguity over a confident wrong pick.
/// (Starting value; re-confirm during Task 3's benchmark calibration.)
pub const RESOLVE_MARGIN: f64 = 0.90;

#[derive(Debug, Clone)]
pub struct ResolverOutcome {
    /// Single call, or ";"-joined near-tied set.
    pub genotype: String,
    pub ambiguous: bool,
    pub top_posterior: f64,
}

/// Decide a genotype from scored candidates. Returns None when nothing was
/// scorable (every score non-finite) — the caller should then keep the
/// cascade's own choice.
pub fn resolve(scores: &[CandidateScore]) -> Option<ResolverOutcome> {
    let mut finite: Vec<&CandidateScore> =
        scores.iter().filter(|c| c.log10_likelihood.is_finite()).collect();
    if finite.is_empty() {
        return None;
    }
    finite.sort_by(|a, b| {
        b.posterior
            .total_cmp(&a.posterior)
            .then_with(|| a.diplotype.cmp(&b.diplotype))
    });
    let top = finite[0];
    if finite.len() == 1 {
        return Some(ResolverOutcome {
            genotype: top.diplotype.clone(),
            ambiguous: false,
            top_posterior: top.posterior,
        });
    }
    let second = finite[1];
    if top.posterior - second.posterior >= RESOLVE_MARGIN {
        Some(ResolverOutcome {
            genotype: top.diplotype.clone(),
            ambiguous: false,
            top_posterior: top.posterior,
        })
    } else {
        let kept: Vec<String> = finite
            .iter()
            .filter(|c| top.posterior - c.posterior < RESOLVE_MARGIN)
            .map(|c| c.diplotype.clone())
            .collect();
        Some(ResolverOutcome {
            ambiguous: kept.len() > 1,
            genotype: kept.join(";"),
            top_posterior: top.posterior,
        })
    }
}
```

Then add to `src/caller/mod.rs`, keeping alphabetical-ish order (after `read_voting;`):

```rust
pub mod resolver;
```

- [ ] **Step 4: Run tests, verify pass**

Run: `cargo test --release resolver:: 2>&1 | tail -6`
Expected: 4 passed.

- [ ] **Step 5: Commit**

```bash
git add src/caller/resolver.rs src/caller/mod.rs
git commit -m "feat(resolver): add likelihood-based diplotype tie-breaker"
```

---

## Task 3: Wire the resolver into `d6_star_caller` for `;`-ambiguous calls

Acts ONLY on calls that already contain `;` (already reported ambiguous today), so it can only improve those — the ~47 unambiguous samples are untouched. This is the safe entry point.

**Files:**
- Modify: `src/main.rs` (add a `resolver_posterior` variable and a resolver block after the `phase_disambiguate` block at ~`main.rs:827`; thread `resolver_posterior` into the JSON output and the `D6Call`)
- Modify: `src/types.rs` (add `resolver_posterior: Option<f64>` to `D6Call`, defaulting to `None` in `D6Call::no_call`)

- [ ] **Step 1: Add the `resolver_posterior` field to `D6Call`**

In `src/types.rs`, add to `struct D6Call` (after `kmer_validation`):

```rust
    /// Top posterior from the likelihood resolver (None if it didn't run).
    pub resolver_posterior: Option<f64>,
```

And in `D6Call::no_call(...)`, add `resolver_posterior: None,` to the returned struct.

- [ ] **Step 2: Build to find every D6Call construction site that now needs the field**

Run: `cargo build --release 2>&1 | grep -E "missing field|error" | head`
There are **three** `D6Call { ... }` literals: two no-call returns at `main.rs:398` and `main.rs:598`, and the success literal at `main.rs:1027`. Add `resolver_posterior: None,` to **both** no-call literals (398 and 598). (The success literal at 1027 gets the populated value in Step 5.) Re-run the build until only the success literal remains flagged (it's fixed in Step 5). Confirm exact lines with `grep -n "D6Call {" src/main.rs`.

- [ ] **Step 3: Confirm the total-CN variable**

The total CN is `total_cn` (defined at `main.rs:410` as `let total_cn = raw_cn_call.d67_cn.unwrap();`, a `u32`), in scope at the insertion point. The existing `call_diplotype` invocation (main.rs:1003) uses `&params.snp_db.nchr`, `&params.var_list`, `data::SNP_FILE`, `&params.star_combinations`, `total_cn` — the resolver block below mirrors exactly these argument sources.

- [ ] **Step 4: Add the resolver block after the `phase_disambiguate` block**

Locate the closing brace of the `phase_disambiguate` block (currently ends at `main.rs:827`, the `}` before `if let Some(ref call) = final_star_allele_call {` at line 829). Insert a `let mut resolver_posterior: Option<f64> = None;` before the `if star_called.call_info.is_some()` block (~`main.rs:795`), and insert this resolver block immediately AFTER the phase_disambiguate block and BEFORE the existing `if let Some(ref call) = final_star_allele_call {` filter at line 829:

```rust
        // Feature: likelihood resolver — when the cascade leaves a ;-ambiguous
        // candidate set, score each candidate's allele multiset and either pick
        // confidently or report honest ambiguity. Falls back silently to the
        // cascade's own choice if the BAM can't be scored.
        if let Some(ref call) = final_star_allele_call {
            if call.contains(';') {
                let candidates: Vec<String> =
                    call.split(';').map(|s| s.trim().to_string()).collect();
                if let Ok(mut r) =
                    open_alignment_file_with_index(bam_path, reference_fasta, index_name)
                {
                    let scores = cyrius_rs::caller::diplotype_caller::score_candidates(
                        &mut r,
                        &params.snp_db.nchr,
                        &params.var_list,
                        data::SNP_FILE,
                        &params.star_combinations,
                        total_cn,
                        &candidates,
                    );
                    if let Some(outcome) = cyrius_rs::caller::resolver::resolve(&scores) {
                        log::info!(
                            "resolver: {} -> {} (ambiguous={}, top_posterior={:.4})",
                            call, outcome.genotype, outcome.ambiguous, outcome.top_posterior
                        );
                        resolver_posterior = Some(outcome.top_posterior);
                        final_star_allele_call = Some(outcome.genotype);
                    }
                }
            }
        }
```

- [ ] **Step 5: Thread `resolver_posterior` into the returned `D6Call`**

At the success-path `D6Call { ... }` literal (`main.rs:1027`, the one whose `genotype: final_star_allele_call,` field is set), add:

```rust
        resolver_posterior,
```

- [ ] **Step 6: Build, run unit tests**

Run: `cargo build --release 2>&1 | tail -3 && cargo test --release 2>&1 | grep "test result" | head -3`
Expected: build `Finished`; all tests pass.

- [ ] **Step 7: Benchmark gate (CRITICAL)**

From repo root `/home/parlar_ai/dev/cyrius-rs`:
Run: `bash benchmark/run_benchmark.sh --strand-bias-all --fuzzy-match --changepoint-hybrid --het-check 2>&1 | tail -5`
Expected: `Correct: 53 / 53`.

If it drops below 53: inspect which sample moved (the `Different` rows). If the resolver flipped a previously-correct `;`-call to a wrong single call, the `RESOLVE_MARGIN` is too low — raise it (e.g. 0.90 → 0.95) and re-run. If a previously-correct single call moved, that's a bug (the resolver should only touch `;`-calls) — debug before proceeding. Do NOT proceed past 53/53 here; this task must not regress.

- [ ] **Step 8: Commit**

```bash
git add src/main.rs src/types.rs
git commit -m "feat: run likelihood resolver on ;-ambiguous diplotype calls"
```

---

## Task 4: Surface the resolver posterior in JSON output

**Files:**
- Modify: `src/main.rs` (JSON assembly around `main.rs:1239-1276`)

- [ ] **Step 1: Add the field to the JSON `call_json` object**

In the `serde_json::json!({...})` call that builds `call_json` (the block containing `"Kmer_validation": kmer_validation_json,` near `main.rs:1275`), add a sibling field:

```rust
        "Resolver_posterior": cyp2d6_call.resolver_posterior,
```

- [ ] **Step 2: Build + benchmark gate**

Run: `cargo build --release 2>&1 | tail -2`
Then from repo root: `bash benchmark/run_benchmark.sh --strand-bias-all --fuzzy-match --changepoint-hybrid --het-check 2>&1 | tail -4`
Expected: build clean; `Correct: 53 / 53` (output-only change, no call movement).

- [ ] **Step 3: Spot-check the JSON for one ambiguous sample**

Run (NA18565 is a known ambiguous sample):
`./target/release/cyrius-rs -i /home/parlar_ai/dev/cyrius-rs/benchmark/<NA18565-bam> --strand-bias-all --fuzzy-match --changepoint-hybrid --het-check -o /tmp/res 2>/dev/null; grep -o 'Resolver_posterior[^,]*' /tmp/res/*.json`
Expected: a numeric `Resolver_posterior` value present (resolver ran on the `;`-call). Exact BAM filename: `ls /home/parlar_ai/dev/cyrius-rs/benchmark/ | grep -i NA18565`.

- [ ] **Step 4: Commit**

```bash
git add src/main.rs
git commit -m "feat: surface resolver posterior in JSON output"
```

---

## Task 5: Migrate the 6 hardcoded special-cases (reproduce-or-promote)

This is the careful, gated part. Process the 6 sites ONE AT A TIME. For each: expose its candidate set so the resolver (Task 3) sees it, then verify the benchmark. Each site is its own commit and its own benchmark run.

**Sites (all in `src/caller/match_star_allele.rs`):**
1. line 319 — `if final_call == ["*10_*10_*4.013", "*10_*36_*4"]` → `*4/*36+*10`
2. line 422 — `if split_call == ["*10", "*10"]` (exon9hyb_star5)
3. line 429 — `if split_call == ["*4", "*4"]` (exon9hyb_star5)
4. line 575 — `else if split_call == ["*4", "*4"]` (star68)
5. line 862 — `if parts == ["*1", "*27", "*32", "*41"]`
6. line 869 — `if parts == ["*1", "*119", "*2", "*41"]`

**Per-site procedure (repeat for each):**

- [ ] **Step A: Read the site and classify it**

Read ~20 lines around the site. Determine:
- (i) What single diplotype string it currently returns.
- (ii) The set of *structurally distinct candidate diplotypes* this case is choosing among. For multiset-different candidates the resolver can rank them; for same-multiset-different-arrangement candidates it will (correctly) tie → ambiguous.
- (iii) Whether the pick is driven by a *real structural signal* (e.g. `spacer_cn` at sites 422/429) rather than an arbitrary/frequency choice.

- [ ] **Step B: Decide reproduce vs promote vs keep**

- If the candidates are **multiset-different** (e.g. site 319: `{*10,*10,*4.013}` vs `{*10,*36,*4}`; sites 862/869: 4 distinct alleles): change the branch to RETURN the candidates `;`-joined instead of the hardcoded single string, so the resolver scores them. Go to Step C.
- If the pick is driven by a **real structural signal** (sites 422/429 spacer-CN, possibly 575): do NOT delete. Instead, refactor the literal equality into a clearly-named, documented helper/rule (e.g. `// structural rule: exon9hyb_star5 with spacer_cn>2 ⇒ *5/*36+*10`) so the intent is explicit. Leave the call unchanged. Commit as a documentation/clarity refactor; skip Steps C–E for this site.
- If unsure, default to the safe path: keep the hardcode, add a `// TODO(resolver): candidates not yet exposed` note, and move on (do not regress).

- [ ] **Step C: Write a regression test for the expected post-migration call**

Add a test to `match_star_allele.rs`'s test module (or `resolver.rs`) asserting the behavior you expect after migration: either the single resolved call (if you expect the scorer to reproduce it) or that the output is the `;`-ambiguous set (if you expect a promote). If you can't construct the scorer inputs in a unit test, rely on the benchmark gate in Step E and add an integration note instead.

- [ ] **Step D: Make the edit**

Replace the hardcoded `return Some("<single>".to_string());` with the candidate set. For site 319 the two `final_call` decompositions are the two candidate **multisets**: `["*10","*10","*4.013"]` = `{*10,*10,*4.013}` and `["*10","*36","*4"]` = `{*10,*36,*4}`. Each must become a diplotype string whose multiset matches exactly:

```rust
        if final_call == ["*10_*10_*4.013", "*10_*36_*4"] {
            // Expose both structurally-consistent arrangements; the resolver
            // (main.rs) scores them by read likelihood and picks or reports
            // honest ambiguity.
            //   {*10,*36,*4}     -> "*4/*36+*10"      (the old hardcoded pick)
            //   {*10,*10,*4.013} -> "*10/*4.013+*10"
            return Some("*4/*36+*10;*10/*4.013+*10".to_string());
        }
```

IMPORTANT: the exact arrangement notation must follow this codebase's conventions (e.g. `+` for tandem, `xN` for duplication). Verify in Step A that each string's multiset (via `parse_multiset`) equals its `final_call` decomposition; the example above is illustrative — do NOT introduce an allele (like `*4`) that isn't in the decomposition.

- [ ] **Step E: Build + benchmark gate for THIS site**

Run: `cargo build --release 2>&1 | tail -2`
From repo root: `bash benchmark/run_benchmark.sh --strand-bias-all --fuzzy-match --changepoint-hybrid --het-check 2>&1 | tail -5`

Decision:
- **53/53** → the resolver reproduced the call (or correctly resolved it). Keep the change. Commit.
- **52/53 (or fewer) and the moved sample is now reported as `;`-ambiguous** → acceptable principled regression (the spec authorizes this). Update `benchmark/getrm_cyp2d6_truth.tsv` expectations ONLY if the benchmark harness treats the `;`-set as a miss AND the user has confirmed; otherwise record the sample in the plan's notes and keep the change. Commit.
- **A sample moved to a *confidently wrong single* call** → REVERT only this site's edit, keeping any unrelated working-tree changes intact. Use `git restore -p src/caller/match_star_allele.rs` and unstage only this site's hunk (or manually re-apply the original hardcoded branch). Do NOT use `git checkout -- <file>` — it would discard unrelated edits. Keep the hardcode and log it as a scorer limitation.

- [ ] **Step F: Commit this site**

```bash
git add -A
git commit -m "refactor(match_star_allele): migrate hardcoded case at <LINE> to resolver (reproduce|promote|kept)"
```

Repeat A–F for all 6 sites. After the last site, run the full benchmark once more and record the final tally + any documented principled regressions in `docs/superpowers/specs/2026-06-04-diplotype-resolver-design.md` under a new "## Migration outcomes" section.

---

## Task 6: Delete the four logging-only scorer modules

Only after Task 5 confirms none of them was needed to reproduce a hardcoded case (i.e. no Step-B decision kept a module alive). If one was kept, exclude it from the deletions below.

**Files:**
- Delete: `src/caller/read_voting.rs`, `src/caller/read_phasing.rs`, `src/caller/af_phasing.rs`, `src/caller/cn_classifier.rs`
- Modify: `src/caller/mod.rs` (remove the 4 `pub mod` lines)
- Modify: `src/main.rs` (remove the 4 CLI flag fields, their `FeatureFlags` fields/assignments, and the logging-only invocation blocks)
- Modify: `src/types.rs` (remove the 4 `FeatureFlags` fields)
- Modify: `../README.md` (remove the 4 flag rows and the 4 module-tree entries — they document deleted features)

- [ ] **Step 1: Find every reference to the four modules**

Run:
```bash
grep -rn "read_voting\|read_phasing\|af_phasing\|cn_classifier" src/ | grep -v "src/caller/read_voting.rs\|src/caller/read_phasing.rs\|src/caller/af_phasing.rs\|src/caller/cn_classifier.rs"
```
This lists `mod.rs`, `main.rs`, and `types.rs` references to remove.

- [ ] **Step 2: Remove the module files and declarations**

```bash
git rm src/caller/read_voting.rs src/caller/read_phasing.rs src/caller/af_phasing.rs src/caller/cn_classifier.rs
```
Remove the 4 `pub mod ...;` lines from `src/caller/mod.rs`.

- [ ] **Step 3: Remove CLI flags, FeatureFlags fields, and invocation blocks in main.rs**

For each of the 4 features: delete the `#[arg(long = "...")]` field in `Cli`, the corresponding `FeatureFlags` field in `src/types.rs`, the `<feature>: cli.<feature>,` assignment, the experimental-features log-line references, and the `if features.<feature> { ... }` invocation block (each is logging-only; search `grep -n "features.read_voting\|features.read_phasing\|features.af_phasing\|features.cn_classifier" src/main.rs`).

- [ ] **Step 3b: Remove the dead flags + modules from README.md**

The four flags are documented in `README.md` (the option table rows at lines ~204 `--read-voting`, ~213 `--cn-classifier`, ~216 `--af-phasing`, ~217 `--read-phasing`) and the module tree (entries around lines ~291-295: `cn_classifier.rs`, `read_phasing.rs`, `read_voting.rs`, and the `af_phasing.rs` entry). Delete those rows/entries so users don't see dead flags. Confirm with `grep -n "read-voting\|read-phasing\|af-phasing\|cn-classifier\|read_voting\|read_phasing\|af_phasing\|cn_classifier" ../README.md` and remove every match that refers to a deleted module. (If Task 5 kept one module alive, leave that one's README entry.)

- [ ] **Step 4: Build until clean**

Run: `cargo build --release 2>&1 | grep -E "error|warning: unused" | head`
Fix any dangling references until no errors/unused warnings. Re-run until clean.

- [ ] **Step 5: Full test + benchmark gate (both flag configs)**

Run: `cargo test --release 2>&1 | grep "test result" | head -3`
From repo root, standard config:
`bash benchmark/run_benchmark.sh --strand-bias-all --fuzzy-match --changepoint-hybrid --het-check 2>&1 | tail -4`
Expected: matches the Task 5 final tally (53/53, or the documented principled-regression tally). The deleted flags must no longer be passed; if any benchmark script references them, remove those references.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "refactor: delete four logging-only scorer modules consolidated into diplotype_caller"
```

---

## Task 7: Final verification + memory update

- [ ] **Step 1: Clean build, all tests, both benchmark configs**

```bash
cargo build --release 2>&1 | tail -2
cargo test --release 2>&1 | grep "test result"
```
From repo root:
```bash
bash benchmark/run_benchmark.sh --strand-bias-all --fuzzy-match --changepoint-hybrid --het-check 2>&1 | tail -4
```
Expected: clean build, all tests pass, final tally as documented.

- [ ] **Step 2: Record durable outcomes IN THE REPO**

Add a `## Migration outcomes` section to `docs/superpowers/specs/2026-06-04-diplotype-resolver-design.md` capturing: the final benchmark tally, each hardcoded-case decision (reproduced / promoted-to-ambiguous / kept), the calibrated `RESOLVE_MARGIN` value, and any module kept alive. This is the committed, durable record (the `.claude` memory file below is outside the repo and is NOT committed).

- [ ] **Step 3: Final commit (repo files only)**

```bash
git add -A
git commit -m "docs: record resolver consolidation outcomes and calibration"
```

- [ ] **Step 4: Update project memory (non-commit follow-up)**

Separately (NOT part of any commit — this path is outside the git repo), update `/home/parlar_ai/.claude/projects/-home-parlar-ai-dev-cyrius-rs/memory/MEMORY.md`: note that `diplotype_caller` is now the active tie-breaking resolver (no longer logging-only), that `read_voting`/`read_phasing`/`af_phasing`/`cn_classifier` were deleted, and record `RESOLVE_MARGIN`'s value.

---

## Notes / deviations from spec

- **Confidence field:** the spec said "populate the existing Confidence field as tiers." Implementation discovered `Confidence` is already a computed score from `caller::confidence`. To avoid clobbering it, the resolver instead (a) drives the genotype string — which flows into the existing `phasing_ambiguous`/`long_read_recommended` logic — and (b) adds a new `Resolver_posterior` JSON field. Net effect matches the spec's intent (calibrated confidence + first-class ambiguity) without destroying existing signal.
- **HG00421 threshold** (`call_exon9gc` 0.75→0.70) is out of scope (CN threshold, not a diplotype pick) — untouched.
- **`phase_disambiguate`** runs before the resolver and is untouched.

# Design: Likelihood-based diplotype resolver, calibrated confidence, and scorer consolidation

**Date:** 2026-06-04
**Status:** Implemented (see Migration outcomes)
**Scope:** short-read only

## Problem

The CYP2D6 caller's decision core is a heuristic cascade (`depth → GMM → diagnostic-SNP
ratio → cnvtag → match_star_allele`). At its hard cases it relies on **6 hardcoded
special-case branches** in `match_star_allele.rs` plus a sample-tuned threshold, and it
ships **five experimental "scorer" modules** (`diplotype_caller`, `read_voting`,
`read_phasing`, `af_phasing`, `cn_classifier`) that compute signal and then **only log**
— none influence the call. The probabilistic `diplotype_caller` is the most principled of
them but is dead weight.

This produces three problems addressed here:
1. Structural ambiguity is resolved by hardcoded guesses instead of evidence or honest
   "ambiguous" reporting.
2. The same read-vs-allele scoring is reimplemented several times, all unused.
3. The hardcoded branches generalize poorly and are the maintainability ceiling.

## Goals

- Turn `diplotype_caller`'s (now-fixed) posterior into a **calibrated per-call
  confidence** and a **first-class ambiguity output**.
- Make `diplotype_caller` the **single** scorer that **breaks ties** among
  structurally-consistent candidate diplotypes.
- Replace the hardcoded special-cases with general scorer-driven reasoning, **reproduce-
  or-promote** per case.

## Non-goals (explicitly out of scope)

- Long-read support, multi-evidence weighting (Approach C), LD/population phasing.
- Any change to depth / GMM / diagnostic-SNP / cnvtag logic.
- `phase_disambiguate` (separate, already-active CN=2 disambiguator — untouched).
- The HG00421 `call_exon9gc` CN threshold (0.75→0.70) — a CN-calling threshold, not a
  diplotype pick; out of scope.
- Validation-set expansion; PharmCAT-format output.

## Success criterion

53/53 on the GeT-RM benchmark is **strongly preferred**. Small principled regressions are
acceptable **only** when the miss is a genuinely ambiguous case the tool now reports
honestly as ambiguous (rather than hardcoding a guess). Each such regression must be
documented.

## Key feasibility finding

`diplotype_caller` scores a diplotype as a **multiset of base alleles** (a pair of
distinct alleles + a copy split summing to total CN). It is **blind to arrangement**.
This maps exactly onto the real biology:

- **Different multisets** (e.g. `{*10,*10,*4.013}` vs `{*10,*36,*4}`) → the scorer **can**
  separate them by read likelihood.
- **Same multiset, different arrangement** (e.g. `*36/*36+*10` vs `*10/*36x2`, both
  `{*10,*36,*36}`) → the scorer gives **identical** scores → which is the **correct**
  signal to report ambiguity, because short reads genuinely cannot resolve arrangement.

## Architecture

A thin **resolver layer** sits between the existing cascade and the output. Nothing
upstream of the cnvtag changes.

```
BAM ─► [depth+GMM+SNP cascade] ─► cnvtag + candidate diplotype set ─┐
                                                                     │
                          ┌──────────────────────────────────────────┘
                          ▼
              [resolver]  if |candidates| == 1 ─► emit, high confidence
                          if |candidates| > 1  ─► diplotype_caller scores
                                                  each candidate's MULTISET
                                                     │
                                   ┌─────────────────┴───────────────────┐
                          top posterior clears margin           within margin
                                   │                                     │
                          single call + confidence        ambiguous ;-set + low
                                                          confidence + Phasing_ambiguous=true
```

- The cascade becomes a candidate **generator** (not decider) at the ambiguity points.
- The scorer runs **once per ambiguous sample**. For unambiguous samples the resolver is
  a no-op, so their calls cannot move — this is what protects the 53/53.

## Components

### C1. Candidate exposure — `caller/match_star_allele.rs`
The 6 hardcoded `if final_call == [...]` / `if split_call == [...]` sites stop *deciding*
and instead emit their candidate set (`Vec<String>`). The public entry can still return a
string for the unambiguous path. The existing `;`-joining at `main.rs:1205`
(`phasing_ambiguous = count_of_diplotypes > 1`) already consumes a multi-candidate set.

### C2. Scorer interface — `caller/diplotype_caller.rs`
Add one public function:
```
score_candidates(reader, nchr, var_list, snp_file, star_combinations,
                 total_cn, candidates: &[String]) -> Vec<CandidateScore>
// CandidateScore { diplotype: String, log10_likelihood: f64, posterior: f64 }
```
Reuses the existing read-collection / per-read-per-allele machinery, but scores the
**explicit multiset** parsed from each candidate string (strip `xN`, split on `/` and
`+` → base-allele→copies map), then normalizes posteriors **across the candidate set
only**. This refactor makes per-candidate scores exact (vs. the current "best split for a
pair").

### C3. Decision rule — new `caller/resolver.rs`
Given scored candidates:
- 1 candidate → emit it.
- `top1.posterior − top2.posterior ≥ MARGIN` → emit top-1 (confident).
- else → emit the near-tied set `;`-joined; `Phasing_ambiguous=true`,
  `Long_read_recommended=true`.

`MARGIN` is a single named constant, calibrated **once** on the benchmark (not per-sample),
set conservatively to favor declaring ambiguity over a confident wrong pick. Documented as
GeT-RM-calibrated.

### C4. Confidence mapping
Populate the **existing** `Confidence` field as tiers: `HIGH` (single candidate or margin
cleared) / `AMBIGUOUS` (within margin). Surface the numeric top posterior in JSON output
(not TSV). No new output columns.

### C5. Consolidation
Keep `diplotype_caller`. Delete `read_voting`, `read_phasing`, `af_phasing`,
`cn_classifier` and their CLI flags/wiring (all confirmed logging-only) — **unless** the
per-case validation finds one reproduces a hardcoded result the scorer cannot.
`phase_disambiguate` is untouched.

## Hardcoded-case migration (the 53/53 safety mechanism)

Rule: **no hardcoded branch is deleted until the scorer is proven to reproduce its result,
or the case is shown genuinely ambiguous and promoted to an honest ambiguous call.**
Per-case, not big-bang.

| Site (match_star_allele.rs) | Candidates | Multiset test → expected |
|---|---|---|
| 319 `*4/*36+*10` | `{*10,*10,*4.013}` vs `{*10,*36,*4}` | Different multisets → scorer reproduces; remove hardcode once confirmed. |
| 422 / 429 exon9hyb_star5 | spacer-CN-gated arrangement | Likely same multiset → scorer can't separate; spacer-CN is the real discriminator → **keep as a named, documented structural rule**. |
| 575 star68 `*4`/`*4` | classify | Reproduce-or-promote. |
| 862 / 869 4-allele sets | `{*1,*27,*32,*41}` / `{*1,*119,*2,*41}` | Different-multiset → scorer ranks; reproduce-or-promote. |

**Per-case workflow (TDD):**
1. Test: run scorer on the candidate set, assert expected ranking.
2. Scorer reproduces → delete branch, keep test as regression guard.
3. True tie (same multiset, no other signal) → delete branch, flow to ambiguous-output
   path, update expected benchmark line to the `;`-set (authorized principled regression).
4. Scorer confidently wrong → **do not delete**; keep hardcode, log as scorer limitation.

After all six: full benchmark. Bar = 53/53, or fewer with every miss documented as a now-
correctly-ambiguous case.

**Honest caveat:** the exon9hyb_star5 spacer-CN rule (and possibly star68) encode real
structural signal, not arbitrary picks. Those are **refactored into named documented
rules**, not deleted. Realistic outcome: remove ~3–4 sample-style hardcodes, convert ~2
into documented structural rules.

## Error handling

- `score_candidates` treats **any** failure (BAM open/fetch error, zero informative reads,
  unparseable candidate, or a candidate whose base allele isn't in `dstar`) as **"cannot
  score"** → resolver falls back to the cascade's current behavior for that sample. Never
  panics, never silently drops a sample. Worst case = today's behavior.
- Scorer BAM access reuses the hardened `open_alignment_file_with_index` → clean no-call on
  bad BAM.

## Testing (TDD)

1. Unit: candidate-string → multiset parser (`*36+*10`, `*10x2`, `*4.013`, `*1/*2`, plus
   malformed / unknown-allele cases).
2. Unit: `score_candidates` on synthetic reads — one clearly-favored multiset (assert
   ranking + posterior normalization); one same-multiset pair (assert equal scores →
   ambiguity).
3. Unit: resolver decision rule at/above/below `MARGIN`; 1-candidate and scorer-failure
   fallbacks.
4. Per-case regression tests (one per migrated hardcode).
5. Integration: full 53-sample benchmark via `run_benchmark.sh` as the final gate.

## Rollout order

Each step ends green (builds, tests pass, benchmark ≥ target):
1. Refactor `diplotype_caller` to expose `score_candidates` — pure addition, no behavior
   change (scorer still only logs).
2. Add `resolver.rs` + multiset parser + decision rule, gated to act **only** on already-
   `;`-ambiguous samples (safest entry: those already report ambiguity, so we can only
   improve them).
3. Migrate the 6 hardcoded cases one at a time (reproduce-or-promote), each behind its own
   test + benchmark check.
4. Populate `Confidence` + JSON posterior from the resolver.
5. Delete the four logging-only modules and their flags.
6. Final full benchmark + both flag configurations
   (`--strand-bias-all --fuzzy-match --changepoint-hybrid --het-check`, and the full
   experimental set).

## Risks

- **Margin overfitting:** the one new tunable, fit on 53 samples. Mitigated by conservative
  setting (prefer ambiguity) and explicit documentation that it's GeT-RM-calibrated.
- **Scorer can't reproduce a hardcode:** mitigated by reproduce-or-promote — the hardcode
  stays if the scorer would regress the call.

## Migration outcomes (implemented 2026-06-04)

**Status:** Implemented on branch `feat/diplotype-resolver`. Final benchmark: **53/53**
GeT-RM (both the standard config and the full surviving-experimental config). 44 unit
tests pass; clean build, zero new warnings.

**Resolver:** `diplotype_caller::score_candidates` (explicit allele-multiset scoring,
fail-closed) + `resolver::resolve` are live. The resolver runs on any `;`-ambiguous
genotype the cascade produces (after `phase_disambiguate`), scoring each candidate's
multiset and either picking confidently or keeping an honest ambiguous set. Its top
posterior is surfaced as `Resolver_posterior` in JSON.

> **Correction (2026-10-07 review).** As implemented here the resolver never scored
> anything: `main.rs` passed the D6+D7 total CN (4 for a normal genome) where
> `build_assignment` needs the D6 copy count (2), so every candidate failed closed
> and `Resolver_posterior` was always null. It also scored main-allele definitions
> (`*2`, `*21`), which, once the CN was fixed, called HG00589 and NA18973 as `*2/*21`
> at posterior 1.0 (truth `*1/*21`). Fixed: `score_candidates` now takes the
> matcher's raw candidates, scores each answer through the exact sub-alleles that
> matched it (e.g. `*1`+`*21.002` vs `*2.001`+`*21`), and takes the D6 copy count from
> them. On GeT-RM all four `;`-ambiguous samples now score as honest 0.5/0.5 ties
> (the reads carry no phase signal between the options); no call changed. The
> site-319 "migration" below is unreachable: `KEPT_SUBALLELES` is empty, so `.013`
> is stripped before that comparison. `RESOLVE_MARGIN` has not been calibrated on
> any sample where the resolver separates candidates.

**`RESOLVE_MARGIN` = 0.10** (near-tie posterior-gap epsilon, not the original "confidence
gap"). During code review the gap-based rule was found to over-include non-competitive
candidates for 3+ candidates; it was replaced with a single coherent near-tie rule: a
candidate joins the reported set iff its posterior is within `RESOLVE_MARGIN` of the top,
and the call is confident iff the top stands alone. Calibrated on GeT-RM only; left at 0.10
because no benchmark sample required tuning it.

**Hardcoded-case dispositions (Task 5).** The plan anticipated migrating several; the
honest result is **1 migrated, 5 kept as documented rules** — most "hardcodes" encode
signals the arrangement-blind multiset resolver genuinely cannot reproduce:

| Site (match_star_allele.rs) | Signal it encodes | Disposition |
|---|---|---|
| `get_final_call_clean` ~319 (`*10_*10_*4.013` vs `*10_*36_*4`) | pick between two **different-multiset** decompositions | **Migrated** → emits `*4/*36+*10;*4.013/*10+*10`; resolver decides. Benchmark-neutral (not in GeT-RM), held 53/53. |
| `exon9hyb_star5` ~422 / ~429 | **spacer copy number** (a depth signal) | **Kept**, documented — resolver is blind to depth. |
| `star68` ~575 | deterministic `*68+*4` tandem **notation** per CN | **Kept**, documented — notation construction, not a tie-break. |
| cn2 `more_than_one_match` ~862 / ~869 | targeted **read-phase booleans** (`var42126938`/`var42127803`) | **Kept**, documented — purpose-built phase signals; no benchmark coverage to validate delegating to the resolver. |

**Consolidation (Task 6).** Deleted four logging-only scorer modules now subsumed by the
single resolver path: `read_voting`, `read_phasing`, `af_phasing`, `cn_classifier` (~1800
lines + their CLI flags + README docs). `diplotype_caller`, `resolver`, `phase_disambiguate`,
`kmer_validation`, `clip_evidence`, `d7_depth` and the rest remain.

**Honest limitations.** (1) The site-319 migration is benchmark-neutral, not
benchmark-verified — no GeT-RM sample exercises it, so the resolver's pick there is
principled but untested on real data. (2) `RESOLVE_MARGIN` is fit to 53 samples. (3) The
five kept rules remain hardcoded; revisiting them needs targeted test data (or long reads).

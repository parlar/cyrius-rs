# Design: Ranked candidate genotypes — read evidence and decision provenance

**Date:** 2026-11-05 (revised the same day after review)
**Status:** Design revised; not implemented
**Scope:** short-read only; output only (no call, filter, confidence, or TSV changes)

## Problem

A call reports one genotype (or a `;` set) and a filter. It does not report which
other genotypes were close, or what separates the call from them. A reader
cannot tell a call that beats every rival by many reads from one that wins at a
single site with a handful of reads, or one that wins only because a rule or a
filter overrode the reads.

## Goal

For each called sample, report two separate things:

1. **Read evidence.** A ranked list of individual genotypes from the call's
   copy-number group. Each is scored by a composite log-likelihood from a count
   model, and each gap to the best genotype is decomposed into one term per site.
2. **Decision provenance.** For every site where a rule or a filter set the
   call's count, the rule, the original read counts, and the count it produced.

`Call_rank`, where the call sits in the read-evidence ranking, is descriptive.
A call can rank below 1. That happens when a rule, a filter, or evidence this
model does not see decided it, and the provenance shows which.

## Review of the first draft (2026-11-05)

The first draft scored candidates with the caller's own Poisson term and left
rule-touched sites out of the score, so that the call would rank 1 by
construction. A review raised five points and two gaps. Each was checked before
changing the design:

| # | Point | Checked how | Result | Change |
|---|---|---|---|---|
| 1 | The caller's Poisson ("score the smaller count, then reverse") is asymmetric | Recomputed `call_reg1_cn`'s terms in Python | Reproduced. CN=2 at 20/20: het vs hom-ref 15.415, het vs hom-alt 1.843. At 19/21 vs 21/19 the large penalty moves from one homozygote to the other (14.239 ↔ 2.099). | The count model is chosen by measurement (Stage A); the Poisson term is the baseline to beat |
| 2 | Leaving rule sites out of the score makes it self-confirming | Recomputed | Reproduced. CN=1 at 5/95 takes the one-copy shortcut (min_read 4); the Poisson term favours 0 copies by 30.48 log10. | Read evidence is kept at every site; rules go to provenance; "call ranks 1" is no longer an invariant |
| 3 | Strand-filtered sites become reference evidence | Recomputed | Reproduced. 0/20 at CN=2 favours 0 over 1 copy by 4.053 log10, and the draft would not have flagged it. | Original counts kept; filtered sites masked and listed in provenance |
| 4 | Summing sites is a composite likelihood, not a joint one | Prior measurement | A median 48.7% of fragments are counted at 2 or more of D6/D7 sites 3–6 (`~/cyrius-research/REPORT.md`, H6 input assumptions). Not yet measured on target sites. | Named a composite log-likelihood; no evidential scale until Stage C; reuse on target sites measured in Stage A |
| 5 | Labels overlap, so rank 2 can repeat rank 1's genotype | Rebuilt the cn2 table in Python (22,192 keys) | Reproduced. `*1/*2` is the label of 13 keys and is inside `*1/*2;*34/*39` for 2 more. | Each individual genotype is ranked once, by its best sub-allele profile |
| — | Tie order can pick different representative keys | Code reading | Holds | The canonical key breaks ties |
| — | `get_final_call_clean` can return `None` | Code reading: four `return None` paths | Holds | Such keys are excluded and counted |

Earlier work bears on the model as well. A beta-binomial with one global error
rate was rejected in H1 stage 3, because one sample drove the rate and it would
be overconfident at indels. So error parameters are set per site class and
fitted across many samples.

## Non-goals

- Copy number and structure. Only genotypes in the table of the called CNV group
  are ranked.
- Priors. No allele or diplotype frequencies (refuted here: H15/H19).
- Probabilities, and any evidential reading of the gap, until Stage C passes.
- A support or resampling measure (bootstrap, site jackknife). That comes later,
  built on this.
- Any change to the call, the filter, `Confidence`, the resolver, or the TSV.

## Stage A — choose the count model by measurement

This stage comes before any feature code, with a locked plan (plan → code →
result). Its result fixes the model that Stage B uses. If it fails, this design
stops.

**Candidate models.** Each gives `P(a | n, k)`: `a` alternate reads out of `n`
reads at a site with `n_s` total copies, `k` of them carrying the variant.

- **Baseline:** the caller's Poisson term, normalised over `k`.
- **Binomial:** `a ~ Binomial(n, p_k)` with
  `p_k = (k/n_s)(1 − ε_c) + (1 − k/n_s) ε_c`, then adjusted for an allele bias
  `β_c` (alternate reads sampled at relative rate `β_c`). This is the standard
  genotype-likelihood form (Li 2011).
- **Beta-binomial:** the binomial plus an overdispersion `ρ_c`.

`c` is a site class from the `variant_type` column. The target file has snp 90,
snp_multi 13, indel 16, indel_multi 5, gene_conversion 16, and there are 19
homology sites. Stage A decides which classes share parameters.

**Truth.** Use benchmark samples whose truth has two alleles without copy-number
change. At **core sites** the expected copies come from the truth. A site is
core for an allele if every sub-allele of that main allele carries it, or none
does. Sites that differ among sub-alleles are left out. The call is not used to
define truth, which avoids circularity.

**Held out.** Two folds over the 53, split by sorted sample ID, alternating.
Fit on one fold and score on the other. Scores:

- mean held-out log score of the true `k` under a flat prior over `k` (a proper
  scoring rule);
- per-site accuracy of the argmax.

**Decision rule** (exact margins locked in the Stage A plan before any run):

- The binomial family is adopted only if it beats the Poisson baseline on
  held-out log score in both folds. Otherwise stop.
- `β_c` and `ρ_c` are kept only if each improves held-out log score by more than
  the locked margin. Otherwise drop the parameter (the simplest model that fits).

**Also measured (descriptive):** the fraction of fragments counted at 2 or more
target sites. This sets how strongly the composite sum overstates independent
evidence.

## Stage B — the feature (opt-in)

It sits behind a new flag, `--candidate-ranking`. Without the flag the JSON field
is absent, and the output is byte-identical to today's.

### Read sites and counts

The 159 target sites are 140 in `data/CYP2D6_target_variant_38.txt` and 19 in
`data/CYP2D6_target_variant_homology_region_38.txt`.

- Counts are the **original** counts, before the NOISY_VAR strand filter in
  `call_cn_var`
  ([call_variants.rs:229-236](../../../cyrius-rs/src/caller/call_variants.rs#L229-L236))
  and before `--strand-bias-all`.
- A site where either filter zeroed the alternate count is **masked**. It
  contributes 0 to every candidate and is listed in provenance. Masking is
  preferred over modelling the filter: the filter exists because those reads are
  suspected artifacts, and scoring them either way would assert something
  unmeasured.
- `n_s` is the same `full_cn` the caller uses at that site.

### Unmodelled tokens

`exon9gc`, `g.42126938C>T`, `g.42127526C>T` and `g.42127556T>C` have no read
site. `g.42130655-42130656insA` takes its copies from its own caller. They
contribute 0. Each candidate lists where it differs from the call at these
tokens (`Unmodelled_differences`). So two genotypes that differ only there tie,
which is the honest result for this model.

### Hypotheses and scores

- **Profiles.** Every key of the table that `match_star` uses for the called
  CNV group (`get_dic(cnvcall)`, or `dhap` for one-copy groups). The largest
  table has 1,551,523 keys
  ([match_star_allele.rs:79](../../../cyrius-rs/src/caller/match_star_allele.rs#L79)).
- **Profile score** `S(K)`: the sum of `log10 P(a_s | n, k_s(K))` over unmasked
  read sites, in `var_list` order, under the Stage A model. This is a
  **composite** log-likelihood.
- **Genotypes.** Each raw diplotype in a key is converted to main alleles and
  labelled by `get_final_call_clean` on that single diplotype, with the call's
  `cnvcall` and `spacer_cn`. A diplotype whose label is `None` is excluded and
  counted (`Hypotheses_unlabelled`).
- **Genotype score.** The maximum `S(K)` over the keys that contain the genotype.
  This is a profile score: the best-fitting sub-allele realisation. Its key is
  the genotype's `Best_profile_key`.

### Ranking

- Genotypes are sorted by score descending (`total_cmp`), then
  `Best_profile_key`, then label.
- Each genotype appears once.
- Ranks use competition ranking ("1, 2, 2, 4"). Genotypes share a rank exactly
  when their scores are equal. Genotypes that share a `Best_profile_key` are
  indistinguishable by these reads.
- `Delta_composite_log10_likelihood` = rank-1 score − genotype score ≥ 0.

### Decomposition

For each listed genotype, against the rank-1 genotype's best key:

- `Discriminating_sites`: one row per unmasked read site where the two profiles'
  expected copies differ. Each row has the reads, both expected copies, and the
  `Site_log10_ratio`. The rows sum to the delta (invariant; tolerance 1e-9).
- `Unmodelled_differences`: as above.

### Decision provenance

Every read site whose call came from something other than the plain argmax gets
a record:

- `Site`
- `Mechanism`, one of:
  - `min_read_rule` ([copy_number_call.rs:54](../../../cyrius-rs/src/depth_calling/copy_number_call.rs#L54));
  - `one_copy_shortcut` ([:14](../../../cyrius-rs/src/depth_calling/copy_number_call.rs#L14));
  - `min_cn_rule` ([:164](../../../cyrius-rs/src/depth_calling/copy_number_call.rs#L164));
  - `homology_cap` ([call_variants.rs:207](../../../cyrius-rs/src/caller/call_variants.rs#L207));
  - `noisy_var_strand_filter`, `strand_bias_all`;
  - `update_variants` ([match_star_allele.rs:761](../../../cyrius-rs/src/caller/match_star_allele.rs#L761));
  - `star68_step` ([:314-326](../../../cyrius-rs/src/caller/match_star_allele.rs#L314-L326));
  - `no_call`.
- `Original_alt`, `Original_ref`
- `Called_copies`

The decision functions report which mechanism set each value. The caller keeps
using only the value, and check B3 proves that no call changes.

### Call rank

The call is the final `Genotype` after the resolver and the phase steps, split
on `;`.

- `Call_rank` is the best rank of any of its genotypes, or `null` if none is in
  the hypothesis space.
- `Call_delta_composite_log10_likelihood` is that genotype's delta.

### Output

```
"Candidate_ranking": {
  "Status": "ranked",
  "Model": "<Stage A model id>",
  "Evidence_scale": "uncalibrated_composite",
  "Hypothesis_space": "<cnv group>",
  "Profiles_evaluated": <int>, "Hypotheses_unlabelled": <int>,
  "Masked_sites": [...], "Unmodelled_tokens": [...],
  "Call_rank": <int or null>, "Call_delta_composite_log10_likelihood": <float or null>,
  "Entries_tied_at_cutoff": <int>,
  "Candidates": [
    {"Rank": 1, "Genotype": "...", "Best_profile_key": "...",
     "Delta_composite_log10_likelihood": 0.0,
     "Discriminating_sites": [], "Unmodelled_differences": []},
    {"Rank": 2, "Genotype": "...", "Best_profile_key": "...",
     "Delta_composite_log10_likelihood": <float>,
     "Discriminating_sites": [
       {"Site": "g.…", "Alt_reads": <int>, "Ref_reads": <int>,
        "Expected_copies_best": <int>, "Expected_copies_candidate": <int>,
        "Site_log10_ratio": <float>}],
     "Unmodelled_differences": [
       {"Token": "…", "Copies_call": <int>, "Copies_candidate": <int>}]}
  ],
  "Decision_provenance": [
    {"Site": "g.…", "Mechanism": "…", "Original_alt": <int>,
     "Original_ref": <int>, "Called_copies": <int or null>}]
}
```

- `Candidates` holds the first `RANKING_LIST_SIZE = 5` genotypes, plus the
  call's best genotype if it ranks lower.
- `Entries_tied_at_cutoff` counts unlisted genotypes that share the fifth
  entry's rank.
- If no ranking is possible, the field is
  `{"Status": "not_ranked", "Reason": "direct_genotype_tag" | "no_table" | "no_call"}`.

### Code layout

- New module `cyrius-rs/src/caller/candidate_ranking.rs`. It is a pure function
  with no BAM I/O.
- The count model goes in its own small module, with the Stage A parameters as
  constants and the measurement that set them cited.
- The decision functions in `copy_number_call.rs` / `call_variants.rs` gain a
  mechanism alongside each value. The values are unchanged.
- `main.rs` keeps the original counts, calls the module once after the final
  call is set (only with the flag), and writes the field.

### Tests (unit, by name)

- `binomial_model_is_symmetric_in_ref_and_alt`: swapping `a` and `n − a`
  mirrors `k` and `n_s − k` when `β = 1`. This pins the fix for point 1.
- `read_evidence_is_kept_where_a_rule_decided`: the 5/95 CN=1 case scores
  against 1 copy, and provenance records `one_copy_shortcut`.
- `strand_filtered_site_is_masked_not_scored`: the 20/20 → 0/20 case contributes
  0 and is in `Masked_sites` and provenance.
- `each_genotype_appears_once`: on the real cn2 table, `*1/*2` is one entry, and
  `*34/*39` is a separate entry.
- `canonical_key_breaks_ties`, `ranking_is_deterministic`.
- `unlabelled_diplotypes_are_counted`.
- `decomposition_sums_to_delta`.
- `ties_share_a_rank_and_are_counted_at_the_cutoff`.
- `resolver_narrowed_call_finds_its_rank`.
- `mechanism_is_reported_for_each_rule`: one case per mechanism, each with a
  control site that is not touched.
- `not_ranked_reasons`: one case per reason.

### Implementation checks (locked before the Stage B run)

These establish consistency, arithmetic and backward compatibility only.

- **B1.** In every output file (53 benchmark, 25% and 10% thinned sets, damaged
  NA19920), every listed candidate's rows sum to its delta within 1e-9.
- **B2.** Two runs of the same binary give identical `Candidate_ranking`.
- **B3.** Without the flag, every field except the build provenance
  (`Tool_commit`, `Flags`) is identical to master `baf3459`. With the flag, every
  existing field is identical as well (`compare_all_fields.py`).
- **B4.** Provenance agrees with the call: for every read site, `Called_copies`
  equals the copy count the call used.

## Stage C — validation before any evidential reading

These checks are required before the gap is described as evidence, and before
the flag can become a default. Each gets its own locked plan.

1. **Held-out truth ranking.** The rank and gap of the true genotype on the 7
   GeT-RM samples not yet used anywhere (HG00423, NA10855, NA10860, NA19109,
   NA19174, NA19226, NA19790; 30x NYGC), and in 2-fold on the 53.
2. **Depth sensitivity.** The same on the 25% and 10% thinned sets. The gap
   should shrink with depth, and the truth must not sit behind a large gap.
3. **Perturbation.** Vary the model parameters within their Stage A uncertainty,
   and vary the caller thresholds (min_read, the 0.8/0.65 cutoffs). Rank-1
   genotypes should be stable where the gap is large.
4. **Known mapping failures.** The damaged NA19920, the multimapper-depleted
   files, and the contamination mixtures. Where the call is wrong, the truth must
   not trail by a large gap. Such a gap would be false confidence.

Until all four pass, `Evidence_scale` stays `uncalibrated_composite` and the
flag stays opt-in.

## Limits

- The ranking can only be as good as the table. A novel allele cannot appear.
- The composite sum counts a fragment once at each site it covers, so it
  overstates independent evidence. Stage A measures by how much on target sites.
- Rules and filters are reported, not judged. Provenance shows where the call
  rests on a rule; it does not measure whether the rule is right.

## References

- Li H. (2011). A statistical framework for SNP calling, mutation discovery,
  association mapping and population genetical parameter estimation from
  sequencing data. *Bioinformatics* 27(21):2987–2993.
- Varin C., Reid N., Firth D. (2011). An overview of composite likelihood
  methods. *Statistica Sinica* 21(1):5–42.

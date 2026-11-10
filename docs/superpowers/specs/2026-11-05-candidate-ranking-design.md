# Design: Ranked candidate genotypes — read evidence and decision provenance

**Date:** 2026-11-05 (revised twice the same day after review)
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

1. **Read evidence.** A ranked list of genotypes from the call's copy-number
   group. Each is scored by a composite log-likelihood from a count model, and
   each gap to the best genotype is decomposed into one term per site.
2. **Decision provenance.** For every site where a rule or a filter set the
   call's count, the rule, the original read counts, and the count it produced.

`Call_rank`, where the call sits in the read-evidence ranking, is descriptive.
A call can rank below 1. That happens when a rule, a filter, or evidence this
model does not see decided it, and the provenance shows which.

## Review record

### Round 1 (first draft)

The first draft scored candidates with the caller's own Poisson term and left
rule-touched sites out of the score, so that the call would rank 1 by
construction.

| # | Point | Checked how | Result | Change |
|---|---|---|---|---|
| 1 | The caller's Poisson ("score the smaller count, then reverse") is asymmetric | Recomputed `call_reg1_cn` in Python | Reproduced. CN=2 at 20/20: het vs hom-ref 15.415, vs hom-alt 1.843. 19/21 vs 21/19 swaps 14.239 ↔ 2.099. | Count model chosen by measurement (Stage A); Poisson is the baseline to beat |
| 2 | Leaving rule sites out makes the score self-confirming | Recomputed | Reproduced. CN=1 at 5/95 takes the one-copy shortcut; Poisson favours 0 copies by 30.48 log10. | Read evidence kept at every site; rules go to provenance; "call ranks 1" dropped |
| 3 | Strand-filtered sites become reference evidence | Recomputed | Reproduced. 0/20 at CN=2 favours 0 over 1 copy by 4.053 log10. | Original counts kept; filtered sites masked and listed |
| 4 | A site sum is a composite likelihood | Prior measurement | A median 48.7% of fragments are counted at 2 or more of D6/D7 sites 3–6 (`~/cyrius-research/REPORT.md`, H6). | Named composite; no evidential scale until Stage C |
| 5 | Labels overlap | Rebuilt the cn2 table in Python (22,192 keys) | Reproduced. `*1/*2` labels 13 keys and is inside `*1/*2;*34/*39` for 2. | Each genotype ranked once |
| — | Tie order; `get_final_call_clean` returns `None` | Code reading | Both hold | Canonical key breaks ties; see Genotype identity |

### Round 2 (second draft)

| # | Point | Checked how | Result | Change |
|---|---|---|---|---|
| 1 | Hybrid profiles can exceed the site's total copies | `star_table.txt`; Dhap3 builder; `CNVTAG_ENTRIES` | Holds. `*10` and `*36` both carry `g.42126611C>G`; Dhap3 (used for exon9hyb) has `*10_*36_*36` with 3 copies of it; that site's `full_cn` for exon9hyb is 2. | Dosage translation section; copies come from the profile, not from `full_cn` |
| 2 | One raw diplotype can format to several genotypes | `get_final_call_clean`, exon9hyb branch | Holds. `*10_*36_*36` → `*36/*36+*10;*10/*36x2`. | Formatter output is split on `;` before dedup; arrangements share a score |
| 3 | Candidate labels need not match final-call labels | Special two-candidate case; `phenotype::sort_diplotypes` | Holds. `*10_*10_*4.013` alone is unformatted, but in the special case it becomes `*4.013/*10+*10`. Output sorting turns `*10/*2` into `*2/*10`. | Matching uses a canonical copy multiset, not the label |
| 4 | Provenance omits `retry_add_one_copy` | Code reading | Holds. It promotes a site from N−2 to N−1 copies and discards which one. | The retry returns the promoted site as metadata, recorded in provenance |
| 5 | Four "untouched" validation samples were development data | `exp-H11/PLAN.md:14`; grep of research manifests | Holds, and worse. HG00423, NA10860, NA19109 and NA19226 are H11 development samples. NA10855, NA19174 and NA19790 appear in exp-H1, H4, C1, R1 and guard manifests. | All 7 become regression cases; an untouched set is a Stage C prerequisite |

### Prior work this design must build on

- **H11** (`~/cyrius-research/exp-H11`; plan `176c0be`, code `e7ea2f5`, result
  `d3caa22`) already scored genotypes by read counts. It used a beta-binomial
  with ρ = 0.01, a per-site leave-one-out background, top-5 labels and a log10
  margin, in the frozen `count_lik.py`. It was tested on 161 held-out HPRC
  samples × 5 depth sets, in cn2, star5, cn3 and cn4 only.
  - As a caller it was "inconclusive". At native depth it gave 118 exact, 40
    ambiguous and 3 wrong, against the cascade's 157/1/3.
  - 180 of its 198 `;` runs sat at margin exactly 0, between options with
    identical variant multisets.
  - It found that "allele fraction = copies / D6 copies" is **false for
    conversion alleles** (NA19777 `*82`, allele fraction 0.19–0.27).
- **H1 stage 3** rejected a beta-binomial with one global error rate. One sample
  drove the rate, and it would be overconfident at indels.

## Non-goals

- Copy number and structure. Only genotypes in the table of the called CNV group
  are ranked.
- Priors (refuted here: H15/H19).
- Probabilities, and any evidential reading of the gap, until Stage C passes.
- Telling tandem arrangements apart. Read counts cannot, so arrangements of one
  profile always tie.
- A support or resampling measure. That comes later, built on this.
- Any change to the call, the filter, `Confidence`, the resolver, or the TSV.

## Dosage translation

A table profile is a multiset of allele copies. At each read site `s`, the model
needs the number of D6-mapped copies `n_s(K)` and how many of them carry the
variant, `k_s(K)`. Both come from the **profile**, not from the caller's
`full_cn`:

- **Exon-9 conversion.** A copy whose allele carries `exon9gc` has D7-like
  sequence in the exon-9 region, and its reads there map to CYP2D7. At read
  sites in that region it contributes to neither `n_s` nor `k_s`. The region is
  `pos ≤ 42126611`, the same boundary `get_total_cn_per_site` uses. So the
  `g.42126611C>G` tokens on `*36`-type copies are not counted there. For
  `*10_*36_*36` at `g.42126611C>G`: `n_s = 1`, `k_s = 1`.
- **Homology sites.** `n_s` = contributing D6 copies + 2 D7 copies. Every one of
  the 28 CNV tags has `total_cn − intron1_upstream = 2` (measured from
  `CNVTAG_ENTRIES`).
- **Everything else.** Every copy contributes. `k_s` is the number of copies
  whose allele carries `s`.
- **Out of scope for now.** Tags with `star13`, `star13intron1` or `star68`.
  These hybrids are not table alleles (no `*13` or `*68` line in
  `star_table.txt`). Their copies are implied by the tag, so the profile cannot
  say which sites they cover. These samples are `not_ranked` with reason
  `hybrid_dosage_undefined`.

Stage A tests whether this translation fits real counts in each tag group.
Groups where it does not fit are also `not_ranked`.

## Stage A — choose the count model by measurement

This stage comes before any feature code, with a locked plan (plan → code →
result). Its result fixes the model for Stage B. If it fails, this design stops.

**Locked plan:** `~/cyrius-research/exp-H20/PLAN.md` (H20, plan commit `5dc4e0f`).
The exact models, folds, thresholds and counts live there, not here.

**Start from H11.** The beta-binomial candidate is H11's frozen `count_lik.py`
model. It is not re-derived.

**Candidate models.** Each gives `P(a | n, k_s, n_s)` for `a` alternate reads of
`n`:

- **Poisson baseline:** the caller's term, normalised over `k`.
- **Binomial:** `p = (k_s/n_s)(1 − ε_c) + (1 − k_s/n_s) ε_c`, then adjusted for
  an allele bias `β_c`. This is the standard genotype-likelihood form (Li 2011).
- **Beta-binomial:** H11's model (ρ, per-site background).

`c` is a site class from `variant_type`: snp 90, snp_multi 13, indel 16,
indel_multi 5, gene_conversion 16, and 19 homology sites. gene_conversion gets
its own `β_c` because of the H11 `*82` finding.

**Data and truth.**
- The 53 GeT-RM samples (main-allele truth).
- The 161 HPRC samples with assembly typing from H11 (`exp-H8/data/truth_T1.tsv`).

Both sets have been seen before. That is acceptable for fitting, but rules both
out for Stage C. Expected copies come from the truth at core sites: sites that
every sub-allele of a main allele carries, or none does. The call is never used
as truth.

**Held out.** Two folds, by sorted sample ID. Fit on one fold and score on the
other. Scores:
- the mean held-out log score of the true `k_s` (a proper scoring rule);
- argmax accuracy.

**Decision rule.** Exact margins are locked in the Stage A plan.
- A candidate model is adopted only if it beats the Poisson baseline on held-out
  log score in both folds. Otherwise stop.
- Each extra parameter (`β_c`, ρ, the background) stays only if it improves the
  held-out log score by more than the locked margin.
- **Transfer test.** Parameters are fitted on cn2 only. Then each other in-scope
  tag group (star5, cn3–cn6, the exon9hyb family) is scored with the dosage
  translation. A group is ranked in Stage B only if its held-out log score per
  site is within the locked margin of cn2's.

**Also measured (descriptive):**
- the fraction of fragments counted at 2 or more target sites;
- each group's sample count, so groups with too few samples are reported as
  untested rather than passed.

## Stage B — the feature (opt-in)

It sits behind a new flag, `--candidate-ranking`. Without the flag the field is
absent, and everything except build provenance is identical to today.

### Read sites and counts

- The 159 target sites: 140 in `data/CYP2D6_target_variant_38.txt` and 19 in
  `data/CYP2D6_target_variant_homology_region_38.txt`.
- Counts are the **original** counts, before the NOISY_VAR strand filter
  ([call_variants.rs:229-236](../../../cyrius-rs/src/caller/call_variants.rs#L229-L236))
  and before `--strand-bias-all`.
- A site where a filter zeroed the alternate count is **masked**. It contributes
  0 and is listed in provenance.

### Unmodelled tokens

`exon9gc`, `g.42126938C>T`, `g.42127526C>T` and `g.42127556T>C` have no read
site. `g.42130655-42130656insA` takes its copies from its own caller. They
contribute 0. Each candidate lists where it differs from the call at them
(`Unmodelled_differences`).

### Profiles and scores

- **Profiles.** Every key of the table `match_star` uses for the called group
  (`get_dic(cnvcall)`, or `dhap` for one-copy groups). There are up to 1,551,523
  ([match_star_allele.rs:79](../../../cyrius-rs/src/caller/match_star_allele.rs#L79)).
- **Hybrid count (amendment 2026-11-10).** Only profiles with as many `exon9gc`
  copies as `Variants_called` has `exon9gc` tokens are scored. The others are
  counted in `Profiles_other_exon9gc_count`. Why:
  - The hybrid count comes from depth, as the CNV group does. The read model
    cannot score it.
  - A profile with more `exon9gc` copies has fewer counted copies at exon-9 sites
    (dosage translation). At `n_s = 0` the site term is 0, the largest possible,
    so that profile ignores the site's reads.
  - Measured at `15c57b4` on 9 hybrid-family or exon-9 benchmark samples: in 7,
    the top entry had more `exon9gc` copies than the call. NA18526 (truth
    `*1/*36+*36+*10`) was topped by `*36_*36_*36_*83`, with 4 copies against 2, by
    +1.82 log10, all at `g.42126611C>G` (17 alt, 33 ref).
  - Every profile in a ranked table has the same number of copies (measured: cn2
    2, cn3 3, exon9hyb 3, exon9hyb_exon9hyb 4, dup_exon9hyb 4). So with the
    hybrid count fixed, `n_s` is the same for every profile at every site.
- **Profile score** `S(K)`: the sum over unmasked read sites, in `var_list`
  order, of `log10 P(a_s | n_s, k_s(K), n_s(K))` under the Stage A model, after
  dosage translation. This is a composite log-likelihood.

### Genotype identity

- **Entries.** Each raw diplotype in a key is converted to main alleles and
  formatted by `get_final_call_clean` on that single diplotype, with the call's
  `cnvcall` and `spacer_cn`. The result is split on `;`. Each piece is one
  **entry** (one arrangement). An unformattable diplotype (`None`, or a raw
  underscore string) is shown in raw main-allele form with
  `Label_formatted: false`. It is not dropped.
- **Display.** Entry labels go through `phenotype::sort_diplotypes`, as the
  `Genotype` field does.
- **Identity.** An entry's identity is its **copy multiset**:
  - split the label on `/`, `+` and `_`;
  - expand `xN` into N copies;
  - map `*A(*B)` to `*B`;
  - apply `main_allele_name` (which keeps `*4.013`);
  - sort.

  So `*10/*2` ≡ `*2/*10`, and `*4.013/*10+*10` ≡ `*10_*10_*4.013`.
  `*36/*36+*10` and `*10/*36x2` share one multiset but are separate entries.
- **Entry score.** The maximum `S(K)` over the keys that produce the entry (a
  profile score: the best sub-allele realisation). Its key is the
  `Best_profile_key`.

### Ranking

- Entries are sorted by score descending (`total_cmp`), then
  `Best_profile_key`, then label.
- Each (multiset, label) pair appears once.
- Ranks use competition ranking. Entries share a rank exactly when their scores
  are equal, so arrangements of one profile always share a rank.
- `Delta_composite_log10_likelihood` = rank-1 score − entry score ≥ 0.

### Decomposition

For each listed entry, against the rank-1 entry's best key:

- `Discriminating_sites`: one row per unmasked read site where the two profiles'
  translated `(k_s, n_s)` differ. Each row has the reads, both `(k_s, n_s)`, and
  `Site_log10_ratio`. The rows sum to the delta (tolerance 1e-9).
- `Unmodelled_differences`: as above.

### Decision provenance

One record per read site whose call came from anything other than the plain
argmax:

- `Site`
- `Mechanism`, one of:
  - `min_read_rule` ([copy_number_call.rs:54](../../../cyrius-rs/src/depth_calling/copy_number_call.rs#L54));
  - `one_copy_shortcut` ([:14](../../../cyrius-rs/src/depth_calling/copy_number_call.rs#L14));
  - `min_cn_rule` ([:164](../../../cyrius-rs/src/depth_calling/copy_number_call.rs#L164));
  - `homology_cap` ([call_variants.rs:207](../../../cyrius-rs/src/caller/call_variants.rs#L207));
  - `noisy_var_strand_filter`, `strand_bias_all`;
  - `update_variants` ([match_star_allele.rs:761](../../../cyrius-rs/src/caller/match_star_allele.rs#L761));
  - `star68_step` ([:314-326](../../../cyrius-rs/src/caller/match_star_allele.rs#L314-L326));
  - `retry_add_one_copy` ([:901](../../../cyrius-rs/src/caller/match_star_allele.rs#L901)), recorded at the promoted site;
  - `no_call`.
- `Original_alt`, `Original_ref`
- `Called_copies`

Whole-call steps are recorded once, in `Matching_path`: `exact`, `masked`,
`retry_add_one_copy`, `fuzzy_match` or `star68_step`.

The functions involved report the mechanism alongside each value (for
`retry_add_one_copy`, the promoted site). The caller keeps using only the value,
and check B3 proves that no call changes.

### Call rank

- The call is the final `Genotype` after the resolver and the phase steps, split
  on `;`.
- `Call_rank` is the best rank of any entry whose copy multiset equals the
  multiset of one of the call's genotypes. It is `null` if none does.
- `Call_delta_composite_log10_likelihood` is that entry's delta.

### Output

```
"Candidate_ranking": {
  "Status": "ranked",
  "Model": "<Stage A model id>",
  "Evidence_scale": "uncalibrated_composite",
  "Hypothesis_space": "<cnv group>",
  "Profiles_evaluated": <int>, "Profiles_other_exon9gc_count": <int>,
  "Matching_path": "<exact | masked | …>",
  "Masked_sites": [...], "Unmodelled_tokens": [...],
  "Call_rank": <int or null>, "Call_delta_composite_log10_likelihood": <float or null>,
  "Entries_tied_at_cutoff": <int>,
  "Candidates": [
    {"Rank": 1, "Genotype": "...", "Copy_multiset": ["*1", "*2"],
     "Label_formatted": true, "Best_profile_key": "...",
     "Delta_composite_log10_likelihood": 0.0,
     "Discriminating_sites": [], "Unmodelled_differences": []},
    {"Rank": 2, ...,
     "Discriminating_sites": [
       {"Site": "g.…", "Alt_reads": <int>, "Ref_reads": <int>,
        "Copies_best": [<k>, <n>], "Copies_candidate": [<k>, <n>],
        "Site_log10_ratio": <float>}],
     "Unmodelled_differences": [
       {"Token": "…", "Copies_call": <int>, "Copies_candidate": <int>}]}
  ],
  "Decision_provenance": [
    {"Site": "g.…", "Mechanism": "…", "Original_alt": <int>,
     "Original_ref": <int>, "Called_copies": <int or null>}]
}
```

- `Candidates` holds the first `RANKING_LIST_SIZE = 5` entries, plus the call's
  best entry if it ranks lower.
- `Entries_tied_at_cutoff` counts unlisted entries that share the fifth entry's
  rank.
- If no ranking is possible, the field is `{"Status": "not_ranked", "Reason": …}`,
  where the reason is one of:
  - `direct_genotype_tag`;
  - `no_table`;
  - `no_call`;
  - `hybrid_dosage_undefined`;
  - `group_failed_transfer_test`;
  - `no_profile_with_called_exon9gc_count` (amendment 2026-11-10: every
    scorable profile has another `exon9gc` count).

### Code layout

- New module `cyrius-rs/src/caller/candidate_ranking.rs`. It is a pure function
  with no BAM I/O.
- The count model and the dosage translation each get their own small module.
  The Stage A parameters are constants, with the measurement that set them cited.
- The decision functions gain a mechanism alongside each value. The values are
  unchanged.
- `main.rs` keeps the original counts, calls the module once after the final
  call is set (only with the flag), and writes the field.

### Tests (unit, by name)

- `binomial_model_is_symmetric_in_ref_and_alt` (with `β = 1`).
- `exon9_converted_copies_do_not_count_in_the_exon9_region`: `*10_*36_*36` at
  `g.42126611C>G` gives `(k, n) = (1, 1)`.
- `homology_sites_add_two_d7_copies`.
- `star13_and_star68_tags_are_not_ranked`.
- `read_evidence_is_kept_where_a_rule_decided`: the 5/95 CN=1 case, with
  provenance `one_copy_shortcut`.
- `strand_filtered_site_is_masked_not_scored`: the 20/20 → 0/20 case.
- `arrangements_are_separate_entries_with_one_rank`: `*36/*36+*10` and
  `*10/*36x2`.
- `copy_multiset_matches_across_formats`: `*10/*2` ≡ `*2/*10`;
  `*4.013/*10+*10` ≡ `*10_*10_*4.013`.
- `unformatted_diplotypes_are_kept_and_flagged`.
- `each_genotype_appears_once`: `*1/*2` and `*34/*39` are separate entries.
- `retry_records_the_promoted_site`.
- `canonical_key_breaks_ties`, `ranking_is_deterministic`.
- `decomposition_sums_to_delta`.
- `ties_share_a_rank_and_are_counted_at_the_cutoff`.
- `resolver_narrowed_call_finds_its_rank`.
- `mechanism_is_reported_for_each_rule`: one case per mechanism, each with an
  untouched control site.
- `not_ranked_reasons`: one case per reason.

### Implementation checks (locked before the Stage B run)

These establish consistency, arithmetic and backward compatibility only.

- **B1.** In every output file (53 benchmark, 25% and 10% thinned, damaged
  NA19920, and the 7 regression samples once available), every listed entry's
  rows sum to its delta within 1e-9.
- **B2.** Two runs of the same binary give identical `Candidate_ranking`.
- **B3.** Without the flag, every field except build provenance (`Tool_commit`,
  `Flags`) is identical to master `baf3459`. With the flag, every existing field
  is identical as well.
- **B4.** For every read site, provenance `Called_copies` equals the count the
  call used.

## Stage C — validation before any evidential reading

These checks are required before the gap is described as evidence, and before
the flag can become a default.

**Prerequisite: an untouched validation set.** No known set qualifies:
- the 53 GeT-RM samples and all 7 extra GeT-RM samples have been used in
  development or earlier experiments;
- the 161 HPRC samples evaluated a closely related model (H11);
- all of these are fitting data in Stage A.

The Stage C plan must name a new set with star-allele truth, and freeze its IDs
and checksums before any result on it is read. Sources to investigate
(unverified): assembly-typed samples outside HPRC R2, and long-read 1000G
samples. Until such a set exists, Stage C cannot pass.

**Checks** (each in its own locked plan):

1. **Held-out truth ranking.** The rank and gap of the true genotype on the
   untouched set. The 7 GeT-RM extras are reported as regression cases only.
2. **Depth sensitivity.** The same on thinned versions of the untouched set,
   plus the existing 25% and 10% sets as regression. The truth must not trail a
   wrong rank 1 by a large gap.
3. **Perturbation.** Vary the model parameters within their Stage A uncertainty,
   and vary the caller thresholds. Rank 1 should be stable where the gap is
   large.
4. **Known mapping failures.** The damaged NA19920, the multimapper-depleted
   files, the contamination mixtures, and conversion alleles such as NA19777's
   `*82`. Where the call is wrong, the truth must not trail by a large gap.

Until all four pass, `Evidence_scale` stays `uncalibrated_composite` and the
flag stays opt-in.

## Limits

- The ranking can only be as good as the table and the dosage translation. A
  novel allele, or a hybrid the translation does not cover, cannot be scored
  correctly.
- The composite sum counts a fragment at each site it covers, so it overstates
  independent evidence.
- Rules and filters are reported, not judged.
- The per-site parameters come from NYGC 30x data. Other pipelines are
  unverified (a Stage C item).

## References

- Li H. (2011). A statistical framework for SNP calling, mutation discovery,
  association mapping and population genetical parameter estimation from
  sequencing data. *Bioinformatics* 27(21):2987–2993.
- Varin C., Reid N., Firth D. (2011). An overview of composite likelihood
  methods. *Statistica Sinica* 21(1):5–42.

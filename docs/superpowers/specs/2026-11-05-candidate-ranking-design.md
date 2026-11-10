# Design: Ranked candidate genotypes — read evidence and decision provenance

**Date:** 2026-11-05 (revised twice the same day after review)
**Status:** Stage A done (M4). Stage B implemented as described below, with the hybrid-count amendment of 2026-11-10. Stage C open.
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
  `hybrid_dosage_undefined`, except the direct-genotype tags listed below.

Stage A tests whether this translation fits real counts in each tag group.
Only the CNV groups that passed are ranked: `cn2`, `cn3`, `exon9hyb`,
`exon9hyb_exon9hyb` and `dup_exon9hyb` (`ranking_scope` in
`candidate_ranking.rs`; exp-H20 `0a4b020`). Every other group is `not_ranked`:
- the direct-genotype tags `star5_star5`, `star13_star13`,
  `star13intron1_star13intron1` and `star5_star5_star68`: reason
  `direct_genotype_tag`;
- other tags with `star13` or `star68`: `hybrid_dosage_undefined`;
- `star5`: `group_failed_transfer_test`. This is the locked Stage A verdict.
  A post hoc analysis, written after the result was read, attributes the
  failure to a depth artefact in the transfer statistic rather than a misfit
  (−0.0035 nats per reading; exp-H20 `results/POSTHOC_star5_transfer.md`). The
  verdict was not changed;
- any other tag: `group_untested_in_stage_a`.

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
  `data/CYP2D6_target_variant_homology_region_38.txt`. One of them,
  `g.42130655-42130656insA`, is an unmodelled token (below), so 158 are scored.
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
  order, of `ln P(a_s | n_s, k_s(K), n_s(K))` under the Stage A model (M4,
  constants in `count_model.rs`), after dosage translation. This is a composite
  log-likelihood. A site with no reads adds 0. The output reports differences
  in log10.
- **Unscorable profiles.** A profile with an allele that has no definition, or
  with too many copies for the likelihood table, is not scored. It is counted in
  `Profiles_unscorable`.

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
- `Unmodelled_differences`: one row per unmodelled token whose count in the
  entry's best profile (`Copies_candidate`) differs from its count in
  `Variants_called` (`Copies_call`). These rows compare with the call's variant
  list, not with the rank-1 entry.

### Decision provenance

Provenance is found by replaying `call_reg1_cn`, `process_raw_call_denovo` and
the homology cap per site. The replay uses the caller's own functions plus
`cn_posterior`. A replay that disagrees with the value the call used reads
`unexplained`, and check B4 requires none. When a rule is recorded (code
behaviour, `site_mechanism`):
- `one_copy_shortcut` and `no_call` are recorded whenever they apply, whether or
  not the plain posterior would have agreed.
- `min_read_rule` applies only when the posterior's top value is not 0, so it
  always changes the value.
- `min_cn_rule` and `homology_cap` are recorded only when they change the
  value. `min_cn_rule` takes the smaller of the two options, so it is recorded
  only when the smaller option is not the top option.

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
  - `star68_step` ([:314-326](../../../cyrius-rs/src/caller/match_star_allele.rs#L314-L326)),
    which cannot appear in a ranked output (below);
  - `retry_add_one_copy` ([:901](../../../cyrius-rs/src/caller/match_star_allele.rs#L901)),
    recorded at the promoted site, per genotype option;
  - `fuzzy_match`, recorded at each token that differs, per genotype option;
  - `no_call`;
  - `unexplained`.
- `Original_alt`, `Original_ref`
- `Per_site_copies`: the per-site copy count before token rules, or `null` when
  the site was not called.
- `Called_copies`: the copy count in the call's final variant list
  (`Variants_called`).
- `Genotype_option`: the genotype options that a whole-call record holds for
  (below), or `null`.

For `retry_add_one_copy` and `fuzzy_match` records the two counts change
meaning. `Per_site_copies` is the count in `Variants_called`, and
`Called_copies` is the count in the genotype option's profile (Whole-call
records, below).

Details:
- **Masked target sites** are also replayed. If the replay disagrees with the
  value the call used, the site reads `unexplained`. Otherwise it carries the
  filter mechanism (`strand_bias_all` or `noisy_var_strand_filter`).
- A target site with no copy-number pattern entry reads `unexplained`.
- **`retry_add_one_copy`** is labelled only when all of these hold:
  - the tag is `cnN` with N >= 3;
  - exactly one token differs between `Variants_called` and the call's first
    raw candidate, not counting uncallable tokens that are absent from
    `Variants_called`;
  - that token's count in the raw candidate is its observed count + 1;
  - its observed count is N - 2.

  Any other difference gives `unexplained`.

Whole-call steps are named in `Matching_path`: `exact`, `masked`,
`retry_add_one_copy`, `fuzzy_match`, `star68_step` or `unexplained`. The path is
derived from the call's first raw candidate against `Variants_called`.

**star68.** Tags that contain `star68` are not ranked (`direct_genotype_tag` or
`hybrid_dosage_undefined`, Output). So `star68_step` cannot appear in a ranked
output, as a record or as a `Matching_path`, and the B4 clause on `star68_step`
records cannot fail. The value stays in both lists for completeness.

**Whole-call records, per genotype option (amendment 2026-11-10b).** For
`retry_add_one_copy` and `fuzzy_match` calls:
- The options are the raw candidates in `Raw_star_allele` whose copy multiset
  is one of the final `Genotype`'s options. If none is, all raw candidates are
  used.
- Raw candidates with the same token profile share one set of records.
- For each distinct profile, there is one record per token whose count differs
  from `Variants_called`:
  - `Per_site_copies` is the count in `Variants_called`;
  - `Called_copies` is the count in that profile;
  - `Genotype_option` is the sorted list of the option labels with that profile.
- Every other record has `Genotype_option: null`.
- Basis: the final review found p25 HG02373 (fuzzy; options `*14/*36+*10`,
  `*36+*10/*65` and `*36+*10/*114`) with a single record, built from the first
  raw candidate, that holds for `*14` only. This replaces "one record per
  differing token of the first raw candidate". It was the third patch on that
  design choice, which was reopened with the user; the user chose per-option
  records.

The replay runs in `candidate_ranking.rs` and calls the caller's functions. It
does not change them. Check B3 proves that no call changes.

### Call rank

- The call is the final `Genotype` after the resolver and the phase steps, split
  on `;`.
- `Call_rank` is the best rank of any entry whose copy multiset equals the
  multiset of one of the call's genotypes. It is `null` if none does.
- An entry is scored by its best profile (`Best_profile_key`, Genotype
  identity). So `Call_rank` ranks the genotype's best entry, which may be
  another sub-allele version than the called raw candidate. Example (run 2,
  benchmark HG00421, `*2/*10x2`): the raw call `*10.001_*10.001_*2.002` has
  `g.42129754G>A` once, as `Variants_called` does, while the entry's best
  profile has it twice.
- `Call_delta_composite_log10_likelihood` is that entry's delta. It is `null`
  when `Call_rank` is.
- `Call_rank` is `null` when no scored entry has the call's copy multiset. One
  case: the call's own profile has another `exon9gc` count than
  `Variants_called`, so it was not scored (Hybrid count).

### Output

```
"Candidate_ranking": {
  "Status": "ranked",
  "Model": "M4 (exp-H20 0a4b020)",
  "Evidence_scale": "uncalibrated_composite",
  "Hypothesis_space": "<cnv group>",
  "Profiles_evaluated": <int>, "Profiles_unscorable": <int>,
  "Profiles_other_exon9gc_count": <int>,
  "Matching_path": "<exact | masked | retry_add_one_copy | fuzzy_match | star68_step | unexplained>",
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
    {"Site": "g.…", "Mechanism": "…", "Original_alt": <int or null>,
     "Original_ref": <int or null>, "Per_site_copies": <int or null>,
     "Called_copies": <int>, "Genotype_option": [<label>, …] or null}]
}
```

- `Profiles_evaluated` counts scored profiles. `Profiles_unscorable` counts
  profiles that could not be scored. `Profiles_other_exon9gc_count` counts
  profiles left out by the hybrid-count rule.
- `Unmodelled_tokens` lists `exon9gc`, `g.42126938C>T`, `g.42127526C>T`,
  `g.42127556T>C` and `g.42130655-42130656insA`.
- `Original_alt` and `Original_ref` are `null` for a token-level record whose
  token is not a read site.
- `Candidates` holds the first `RANKING_LIST_SIZE = 5` entries, plus the call's
  best entry if it ranks lower.
- `Entries_tied_at_cutoff` counts unlisted entries that share the fifth entry's
  rank. An unlisted call entry tied at the cutoff is listed and not counted.
- With the flag, every output carries the field. A sample that stops before a
  CNV group or total copy number gets `not_ranked` with reason `no_call`
  (amendment 2026-11-10).
- If no ranking is possible, the field is `{"Status": "not_ranked", "Reason": …}`,
  where the reason is one of:
  - `direct_genotype_tag`;
  - `hybrid_dosage_undefined`;
  - `group_failed_transfer_test`;
  - `group_untested_in_stage_a`;
  - `no_table`;
  - `no_call`;
  - `no_scorable_profile` (no profile could be scored);
  - `no_profile_with_called_exon9gc_count` (amendment 2026-11-10: every
    scorable profile has another `exon9gc` count).

### Code layout

- New module `cyrius-rs/src/caller/candidate_ranking.rs`. It is a pure function
  with no BAM I/O.
- The count model and the dosage translation each get their own small module.
  The Stage A parameters are constants, with the measurement that set them cited.
- The decision functions are unchanged. The module replays them per site to
  find the mechanism (Decision provenance).
- `main.rs` keeps the original counts, calls the module once after the final
  call is set (only with the flag), and writes the field.

### Tests (unit, by name)

The tests as built, in `candidate_ranking.rs` unless noted. Where an earlier
version of this list named a test that was folded into another, the earlier
name follows in brackets.

Count model (`count_model.rs`):
- `matches_the_stage_a_python_model`, `is_a_distribution_over_the_alt_count`,
  `variant_types_map_to_classes`.
- `binomial_model_is_symmetric_in_ref_and_alt`: with `β = 1`, for ρ = 0 and
  ρ > 0. The fitted SNP bias breaks the symmetry.

Sites and genotype identity:
- `there_are_158_scored_sites_19_of_them_homology`.
- `exon9_converted_copies_do_not_count_in_the_exon9_region`: `*10_*36_*36` at
  `g.42126611C>G` gives `(k, n) = (1, 1)`.
- `homology_sites_add_two_d7_copies`.
- `copy_multiset_matches_across_formats`: `*10/*2` ≡ `*2/*10`;
  `*4.013/*10+*10` ≡ `*10_*10_*4.013`.
- `arrangements_are_separate_entries_with_one_multiset`
  [`arrangements_are_separate_entries_with_one_rank`]: `*10_*36_*36` gives two
  formatted entries with one multiset. They share a rank because they share a
  profile.
- `a_cn2_profile_formats_to_its_main_alleles`.
- `unformatted_profiles_are_kept_and_flagged`
  [`unformatted_diplotypes_are_kept_and_flagged`].

Decision provenance:
- `each_rule_is_named_and_an_untouched_site_is_not`
  [`mechanism_is_reported_for_each_rule`, `read_evidence_is_kept_where_a_rule_decided`]:
  `one_copy_shortcut` (5/95 at CN 1), `min_read_rule`, `min_cn_rule`,
  `homology_cap` and `no_call`, with an untouched control site. It tests the
  mechanism, not that the site keeps its counts.
- `a_replay_that_disagrees_with_the_call_reads_unexplained`,
  `masked_site_whose_decision_disagrees_with_the_replay_is_unexplained`,
  `missing_copy_number_pattern_is_unexplained`.
- `unmasked_target_site_reaches_site_mechanism`,
  `homology_cap_through_sites_and_provenance`.
- `filters_mask_sites_and_token_rules_are_recorded`
  [`strand_filtered_site_is_masked_not_scored`]: `strand_bias_all`,
  `noisy_var_strand_filter` and `update_variants`. A masked site keeps its
  original alternate count and its mask.
- `star68_token_diff_is_star68_step`.
- `matching_paths`, `matching_path_without_a_matched_candidate`.
- `retry_record_names_the_promoted_site_and_its_option`
  [`retry_records_the_promoted_site`]: a cn3 retry gives exactly one record.
- `fuzzy_records_are_made_per_genotype_option`: the p25 HG02373 shape, with one
  record per option, and no record for a raw candidate the call does not name.
- `fuzzy_raws_with_one_profile_share_records_and_an_unnamed_call_uses_all_raws`.

Ranking:
- `an_exact_cn2_call_ranks_first_and_its_equivalent_ties_with_it`
  [`each_genotype_appears_once`]: also checks that the rows sum to the delta.
- `ranking_is_deterministic`.
- `canonical_key_breaks_ties`: in an all-tied table, the listed order follows
  the sorted key names.
- `ties_share_a_rank_and_are_counted_at_the_cutoff`,
  `a_call_tied_beyond_the_cutoff_is_listed_not_counted`.
- `a_call_below_the_list_is_appended` [`decomposition_sums_to_delta`, with the
  test above]: an appended entry with rows that sum to its delta.
- `resolver_narrowed_call_finds_its_rank`: `Call_rank` is the best rank over
  the call's `;` options, in either order.
- `profiles_with_another_exon9gc_count_are_not_scored`: also gives
  `no_profile_with_called_exon9gc_count`.
- `not_ranked_reasons` [`star13_and_star68_tags_are_not_ranked`]: every reason
  that `ranking_scope` gives, and `no_call`. `no_table` and
  `no_scorable_profile` have no test.

Integration (`tests/output_provenance.rs`):
- `candidate_ranking_is_present_on_an_early_no_call_only_with_the_flag`.

### Implementation checks (locked before the Stage B run)

These establish consistency, arithmetic and backward compatibility only. The
locked wording, with the vacuity guards, is in the Stage B plan
(`docs/superpowers/plans/2026-11-05-candidate-ranking-stage-b.md`, Locked
checks). Where the two differ, the plan governs.

- **B1.** In every output file (53 benchmark, 25% and 10% thinned, damaged
  NA19920, and the 7 regression samples once available), every listed entry's
  rows sum to its delta within 1e-9.
- **B2.** Two runs of the same binary give identical `Candidate_ranking`.
- **B3.** Without the flag, every field except build provenance (`Tool_commit`,
  `Flags`) is identical to master `baf3459`. With the flag, every existing field
  is identical as well.
- **B4.** No `Decision_provenance` record and no `Matching_path` reads
  `unexplained`. Every `update_variants` or `star68_step` record has
  `Called_copies` equal to that token's count in `Variants_called`.

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

## Known limits (Stage B measurements)

Source: the outputs of run 2 of the locked checks, with the binary built at
`7d68260`. That is 160 files: 53 benchmark, 53 at 25% depth, 53 at 10% depth,
and the damaged NA19920. Run 2 predates the per-option records (`63ee1ab`),
which change only `Decision_provenance`. Re-run at `63ee1ab`, p25 HG02373
differs from run 2 only in `Decision_provenance` and `Tool_commit`. No model
change followed any of these measurements.

### Calls below rank 1

54 of the 160 outputs are ranked. The call ranks 1 in 42 of them and below 1 in
12. No ranked output has `Call_rank: null`. In 11 of the 12, the call names the
truth among its options. The exception is p25 NA18973.

For each of the 12, the table gives the site that carries most of the gap, read
from the call entry's `Discriminating_sites` against the rank-1 entry. Gap and
site ratio are in log10. Where the entry has a second row, that row is given
after a semicolon.

| Group | Set | Sample | Truth | Call | Rank | Gap | Site | Alt/ref | (k, n) rank 1 / call | Site ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | benchmark | HG00463 | `*36+*10/*36+*10` | `*36+*10/*36+*10` | 4 | 1.68 | `g.42127634C>A` | 8/90 | (1, 4) / (0, 4) | 1.68 |
| 1 | benchmark | NA18545 | `*5/*36x2+*10x2` | `*5/*36x2+*10x2` | 4 | 1.86 | `g.42127634C>A` | 9/95 | (1, 4) / (0, 4) | 1.86 |
| 1 | benchmark | NA18563 | `*1/*36+*10` | `*1/*36+*10` | 14 | 1.87 | `g.42127634C>A` | 9/70 | (1, 3) / (0, 3) | 1.87 |
| 1 | benchmark | NA18572 | `*36+*10/*41` | `*36+*10/*41` | 3 | 0.31 | `g.42127634C>A` | 9/102 | (1, 3) / (0, 3) | 0.31 |
| 1 | benchmark | NA18617 | `*36+*10/*36+*10` | `*36+*10/*36+*10` | 4 | 0.78 | `g.42127634C>A` | 5/79 | (1, 4) / (0, 4) | 0.78 |
| 1 | benchmark | NA18632 | `*36+*36+*10/*52` | `*36+*36+*10/*52` | 2 | 1.19 | `g.42127634C>A` | 7/93 | (1, 4) / (0, 4) | 1.19 |
| 1 | 25% | HG02373 | `*14/*36+*10` | `*14/*36+*10;*36+*10/*65;*36+*10/*114` | 2 | 1.11 | `g.42127634C>A` | 2/23 | (1, 3) / (0, 3) | 1.11 |
| 1 | 25% | NA18565 | `*10/*36x2` | `*10/*36x2;*36/*36+*10` | 2 | 1.11 | `g.42127634C>A` | 2/23 | (1, 3) / (0, 3) | 1.11 |
| 2 | 25% | NA18959 | `*2/*36+*10` | `*2/*36+*10` | 13 | 2.72 | `g.42126634A>C` | 2/29 | (1, 5) / (0, 5) | 1.48; `g.42127634C>A` 3/31: 1.24 |
| 2 | 25% | NA18973 | `*1/*21` | `*2/*21` | 2 | 1.58 | `g.42126634A>C` | 2/23 | (1, 4) / (0, 4) | 1.58 |
| 3 | benchmark | HG01094 | `*1/*31` | `*1/*31` | 2 | 1.10 | `g.42128218-42128219insG` | 1/41 | (1, 2) / (0, 2) | 1.10 |
| 4 | benchmark | NA18959 | `*2/*36+*10` | `*2/*36+*10` | 2 | 0.13 | `g.42130692G>A` | 39/43 | (1, 3) / (2, 3) | 0.13 |

Groups by mechanism:

1. **Alternate-read background at g.42127634C>A (8 calls).** The site is the only
   row, so it carries the whole gap. Each rank-1 entry carries `*177`, the only
   allele in `star_table.txt` with this variant. Across run 2's 53 benchmark
   outputs (`Variant_raw_count`), no truth carries `*177`. Even so, all 13
   samples whose truth has `*36` show alt reads at this site, at 3.4-12.6%. Of
   the other 39 with reads, 13 show alt reads, up to 26.3% (HG01190,
   `*68+*4/*5`). The M4 SNP error rate (ε ≈ 0.0004) does not model this
   background. It is left to Stage C.
2. **Alternate-read background at the homology site g.42126634A>C, at 25% depth
   (2 calls).** The rank-1 entries carry `*85`, the only allele with this
   variant: `*21/*85` for NA18973, and `*177_*36_*85` for NA18959, which also
   carries `*177` (the second row).
   - At full depth the same samples have 3 alt reads of 98 and 3 of 124 there.
     Across the 53 benchmark outputs the alt fraction at this site is at most
     3.2%, and no truth carries `*85`.
   - At 25% depth, 2 alt reads remain of 25 and of 31. With the two D7 copies
     (n = 4 and 5) and the homology class (ε ≈ 0.0007), this favours one
     carrier copy over none by 1.58 and 1.48.
   - Here the ranking does not support the truth. The call `*2/*21` is wrong,
     and the truth `*1/*21` is not among the listed entries.
3. **One stray read at an indel site (1 call).** HG01094 is the only benchmark
   sample without `*21` that has an alt read at g.42128218-42128219insG.
   - Under M4's indel class (ρ = 0.316), its 1 read in 42 favours one copy over
     none by 1.10. 0 alt reads in 42 would favour none by only 1.44.
   - With the same ε and β but ρ = 0, the same 1 of 42 favours none by 7.44.
     The overdispersion is what lets one read reverse the order.
4. **A rule decides against the reads (1 call).** For benchmark NA18959, the
   per-site call at g.42130692G>A was 1 copy, and `update_variants` set it to
   2. The provenance record shows both counts (`Per_site_copies` 1,
   `Called_copies` 2).
   - The reads, 39 alt and 43 ref at n = 3, favour 1 copy by 0.13.
   - The rank-1 entry `*2_*36_*39` keeps the per-site count. The call agrees
     with the truth. This works as designed.

The site ratios of groups 2-4, and NA18563's, are reproduced by a Python
transcription of `site_ln_likelihood` that matches the reference values in
`matches_the_stage_a_python_model`. The values for 0 of 42, and for ρ = 0, come
from the same transcription.

### Other limits

- **Lookalike ties.** Entries that differ from the call only at unmodelled
  tokens, or only by phase, tie with it at delta 0. On the Dhap2 table, `*1/*2`
  ties with `*2/*27`, `*2/*127`, `*2/*139` and `*34/*39`.
  - Basis: `*27` = g.42126938C>T, `*139` = g.42127526C>T, and `*127` =
    g.42127556T>C + g.42127526C>T. All of these tokens are unmodelled. `*34/*39`
    carries the same variants as `*1/*2`, in a different phase.
  - In run 2, 29 of the 42 calls at rank 1 share that rank with at least one
    other listed entry, and none of those entries has a discriminating row.
    Examples: benchmark NA19920 (`*1/*4x2`, 5 entries at rank 1) and HG00337
    (`*2x2/*22`, 3).
- **Hybrid count, before and after the amendment.** Nine samples, in the order
  NA18526, HG00463, NA18563, NA18564, NA18565, NA18572, HG02373, NA18545,
  HG00421.
  - Call ranks before (commit `15c57b4`): 2, 6, 38, 1, 2, 5, 1, 6, 1.
  - After (`98cb61d`): 1, 4, 14, 1, 1, 3, 1, 4, 1.
  - Run 2 gives the same ranks as after.
- **Runtime.** Measured on the release build of the code committed as
  `63ee1ab`, built from the working tree just before that commit.
  - Flags: `--strand-bias-all --fuzzy-match --changepoint-hybrid --het-check`.
  - Two runs each, on a shared 72-core machine at a load average of about 5.
  - Wall time and maximum resident memory, flag off to on:
    - NA18563 (exon9hyb): 4.27-4.69 s to 7.16-7.19 s (+2.5 to 2.9 s), and
      558 MiB to 644-646 MiB.
    - HG01094 (cn2): 1.15-1.16 s to 1.22 s (+0.06 to 0.07 s), and 18-20 MiB
      in both.
- **Scores are not probabilities.** They are an uncalibrated composite (fragment
  reuse median 96.1%, exp-H20). Candidate ranks are descriptive and do not
  change the call.

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

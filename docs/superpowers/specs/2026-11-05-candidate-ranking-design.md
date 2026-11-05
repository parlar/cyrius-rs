# Design: Ranked candidate genotypes with a site-wise likelihood decomposition

**Date:** 2026-11-05
**Status:** Design agreed; not implemented
**Scope:** short-read only; output only (no call, filter, confidence, or TSV changes)

## Problem

A call reports one genotype (or a `;` set) and a filter. It does not report which
other genotypes were close, or what evidence separates the call from them. A
reader cannot tell a call that beats every rival by many reads from one that
wins at a single site with a handful of reads, or one that wins only because of
a rule.

## Goal

For each called sample, report:

1. a ranked list of candidate genotypes from the call's own copy-number group;
2. for each listed candidate, the log10 likelihood ratio to the best candidate,
   decomposed exactly into one term per discriminating site;
3. where the reported call sits in that ranking (`Call_rank`).

The ranking uses the evidence and decision rules the call itself used, so that
it explains the call rather than offering a second opinion.

## Non-goals

- Copy number and structure. Only genotypes in the table of the called CNV group
  are ranked. Whether the copy number is right is not tested here.
- Priors. No allele or diplotype frequencies (refuted here before: H15/H19 in
  `~/cyrius-research/REPORT.md`).
- Probabilities. The output reports likelihood ratios, not posteriors. Nothing in
  this model is calibrated.
- A support or resampling measure (bootstrap, site jackknife). That is the next
  step, and it will be validated against known-wrong calls before it is reported.
- Any change to the call, the filter, `Confidence`, the resolver, or the TSV.

## Definitions

### Read sites and their likelihoods

The 159 target sites in `data/CYP2D6_target_variant_38.txt` (140) and
`data/CYP2D6_target_variant_homology_region_38.txt` (19) are the **read sites**.
Each has an alternate read count and a reference read count. The counts are the
ones the call used:

- `var_alt` / `var_ref` after the `--strand-bias-all` filter and after the
  NOISY_VAR strand filter in `call_cn_var`
  ([call_variants.rs:229-236](../../../cyrius-rs/src/caller/call_variants.rs#L229-L236));
- `var_homo_alt` / `var_homo_ref` for the 19 homology sites.

Each site has a total copy number `n_s`. It is `get_total_cn_per_site(cnvtag)` for
the 140 sites and the total D6+D7 copy number for the 19 homology sites, the
same `full_cn` the call passes to `call_reg1_cn`.

For a site with `a` alternate and `r` reference reads, the **site
log-likelihood** of `k` variant copies (`0 ≤ k ≤ n_s`) is the natural-log
Poisson term that `call_reg1_cn`
([copy_number_call.rs:10](../../../cyrius-rs/src/depth_calling/copy_number_call.rs#L10))
computes before it normalises. That means the same expected depth for each `k`
(with `ERROR_RATE` at `k = 0` and `k = n_s`), the same "score the smaller of the
two counts" step, and the same reversal. A new helper
`site_ln_likelihoods(n_s, a, r) -> Option<Vec<f64>>` computes it with `ln_pmf`, so
high depth does not underflow. It returns `None` when `a + r = 0` and when
`n_s = 0`, which are the cases where `call_reg1_cn` computes no Poisson term.
`call_reg1_cn` itself is not changed.

### Rule-determined tokens

Some copy counts in the call are set by a rule, not by the site likelihood. Such
a token is **rule-determined**. It adds nothing to a candidate's likelihood.
Instead, a candidate's count at that token is compared with the counts the call
allowed.

A token is rule-determined if any of these holds:

1. It has no read site: `exon9gc`, `g.42126938C>T`, `g.42127526C>T`,
   `g.42127556T>C`. These are the four table tokens in `data/star_table.txt`
   with no line in either target file, apart from `NA`.
2. It is `g.42130655-42130656insA`, whose copies come from its own caller
   (`call_var42130655ins_a`).
3. It is a read site whose per-site decision differs from the argmax of its site
   log-likelihoods. This covers:
   - the min-read rule ([copy_number_call.rs:54](../../../cyrius-rs/src/depth_calling/copy_number_call.rs#L54));
   - the one-copy shortcut ([copy_number_call.rs:14](../../../cyrius-rs/src/depth_calling/copy_number_call.rs#L14));
   - the `min(cn1, cn2)` rule ([copy_number_call.rs:164](../../../cyrius-rs/src/depth_calling/copy_number_call.rs#L164));
   - the homology cap ([call_variants.rs:207](../../../cyrius-rs/src/caller/call_variants.rs#L207));
   - `n_s = 0`.
4. It is a read site whose count in the final variant list (`Variants_called`)
   differs from its per-site decision. This covers `update_variants`
   ([match_star_allele.rs:761](../../../cyrius-rs/src/caller/match_star_allele.rs#L761))
   and the `*68` step.

The **allowed counts** at a rule-determined token are the single count it has in
`Variants_called`, with one exception. In `*68` groups, the allowed counts at
`g.42130692G>A` are every count the `*68` step tried: the observed count down to
the observed count minus `min(cn, count)`
([match_star_allele.rs:314-326](../../../cyrius-rs/src/caller/match_star_allele.rs#L314-L326)).

A read site whose per-site decision is "unknown" (no call) is **not**
rule-determined. Its site log-likelihoods enter the score as soft evidence. A
site with no reads contributes zero to every candidate.

### Hypothesis space

The hypothesis space is every key of the table that `match_star` uses for the
called CNV group: `get_dic(cnvcall)`, or `dhap` for the one-copy groups. A key's
**expected copies** at a token is the number of times the token occurs in the
key (`NA` = no tokens).

Its size ranges up to the largest table, which has 1,551,523 keys
([match_star_allele.rs:79](../../../cyrius-rs/src/caller/match_star_allele.rs#L79)).

### Score, rule differences, and label

For each key `K`:

- **Score** `S(K)` is the sum over read sites that are not rule-determined of
  `site_ln_likelihoods(...)[c_s(K)]`, converted to log10. Sites are summed in
  `var_list` index order, so equal inputs give bit-identical sums.
- **Rule differences** `D(K)` is the number of rule-determined tokens at which
  `c(K)` is not an allowed count.
- **Label** is the clean genotype string `get_final_call_clean` produces for
  `convert_to_main_allele(table[K])`, called with the same `cnvcall` and
  `spacer_cn` as the call. Sub-alleles therefore merge (`*2` and `*2.001` are one
  genotype). A label that contains `;` is a set of genotypes that this evidence
  cannot separate (observationally equivalent).

The fragment-phase step, the cn2 phase rules, and the resolver are not applied
to labels. They are separate evidence, and they can narrow the call after
matching.

### Ranking

Keys are sorted by `S` descending (`total_cmp`), then `D` ascending, then
label. Walking that order, each distinct label becomes one **entry**, taking its
first (best) key. The entry's

    Delta_log10_likelihood = S(best key of rank-1 entry) − S(entry's key) ≥ 0

Ranks use standard competition ranking ("1, 2, 2, 4"). Two entries share a rank
only when their `S` and `D` are exactly equal.

For an exact match (the final variant list is itself a table key and there are
no unknown sites), the call's key has the highest possible `S` and `D = 0`. Any
key with equal `S` differs from the call only at rule-determined tokens, so it
has `D ≥ 1`. **So for an exact match the call ranks 1 by construction.** If
check C1 below fails, it means there is a bug, not a finding.

### Decomposition

For each listed entry, against the rank-1 entry's key:

- `Discriminating_sites`: one row per read site (not rule-determined) where the
  two keys' expected copies differ. Each row has the reads, both expected
  copies, and the site's `Log10_likelihood_ratio` (rank-1 term minus candidate
  term). The rows sum to `Delta_log10_likelihood` (invariant; tolerance 1e-9).
- `Rule_determined_differences`: one row per rule-determined token where the
  entry's count is not allowed. Each row has the allowed counts and the entry's
  count.

### Call rank

The call is the final `Genotype` after the resolver and the phase steps, split
on `;`. `Call_rank` is the smallest rank of an entry whose label shares at least
one genotype with the call. It is `null` if no entry does.

## Output

There is one new JSON field, `Candidate_ranking`. No existing field changes.

```
"Candidate_ranking": {
  "Status": "ranked",
  "Model": "site_poisson",
  "Hypothesis_space": "<cnv group>",
  "Hypotheses_evaluated": <number of table keys>,
  "Observed_profile_in_table": <bool: final variant list is a table key>,
  "Rule_determined_tokens": ["exon9gc", ...],
  "Call_rank": <int or null>,
  "Entries_tied_at_cutoff": <int>,
  "Candidates": [
    {"Rank": 1, "Genotype": "...", "Delta_log10_likelihood": 0.0,
     "Discriminating_sites": [], "Rule_determined_differences": []},
    {"Rank": 2, "Genotype": "...", "Delta_log10_likelihood": <float>,
     "Discriminating_sites": [
       {"Site": "g.…", "Alt_reads": <int>, "Ref_reads": <int>,
        "Expected_copies_best": <int>, "Expected_copies_candidate": <int>,
        "Log10_likelihood_ratio": <float>}],
     "Rule_determined_differences": [
       {"Token": "…", "Allowed_copies": [<int>, …], "Expected_copies_candidate": <int>}]}
  ]
}
```

- `Candidates` holds the first `RANKING_LIST_SIZE = 5` entries in rank order,
  plus the call's entry if it ranks lower.
- `Entries_tied_at_cutoff` counts unlisted entries that share the fifth entry's
  rank, so nothing is cut silently.
- When no ranking is possible, the field is
  `{"Status": "not_ranked", "Reason": "<reason>"}`, where the reason is one of:
  - `direct_genotype_tag` (`star5_star5`, `star13_star13`,
    `star13intron1_star13intron1`, `star5_star5_star68`);
  - `no_table`;
  - `no_call`.

## Code layout

- New module `cyrius-rs/src/caller/candidate_ranking.rs`. It is a pure function
  of the per-site evidence, the table, `cnvcall`, `spacer_cn`, the final variant
  list, the `*68` allowed counts, and the final call. It does no BAM I/O.
- `site_ln_likelihoods` goes in `copy_number_call.rs`, next to `call_reg1_cn`.
- `main.rs` builds the per-site evidence from values it already holds, calls the
  module once after the final call is set, and writes the JSON field.
- `match_star` exposes, without changing its result, whether the final variant
  list is a table key, and the `*68` allowed counts.

## Tests (unit, by name)

- `site_likelihoods_match_call_reg1_cn`: over a grid of `(n_s, a, r)`, wherever
  `call_reg1_cn` returns two options, the argmax and the normalised top two
  probabilities of `site_ln_likelihoods` match `cn1`, `cn2`, `prob1` and `prob2`
  to three decimals. Wherever it returns a single value from the 0.9 cutoff, the
  argmax matches that value.
- One test per rule-determined clause (1–4) and one for the `*68` allowed set.
  Each builds a site the rule touches and checks it is classified
  rule-determined. Each also has a control: an untouched site that is not.
- `exact_match_ranks_first`: on the real table, an exact-match profile gives
  `Call_rank == 1`.
- `decomposition_sums_to_delta`.
- `sub_alleles_merge_into_one_entry`.
- `ranking_is_deterministic`: two runs give identical output, including tie
  order.
- `ties_share_a_rank_and_are_counted_at_the_cutoff`.
- `resolver_narrowed_call_finds_its_entry`.
- `not_ranked_reasons`: one case per reason.

## Validation

The analysis plan, including these checks, is committed before any run
(plan → code → result).

**Pass/fail checks**

- **C1.** In the 53 benchmark samples, every sample with
  `Observed_profile_in_table = true` and empty `Uncallable_target_sites` has
  `Call_rank == 1`. Any failure is a bug, to be fixed before the result commit.
- **C2.** For every listed candidate in every output file (53 benchmark, 25% and
  10% thinned sets, damaged NA19920), the `Discriminating_sites` rows sum to
  `Delta_log10_likelihood` within 1e-9.
- **C3.** Every existing JSON field is identical to master `baf3459` on the same
  files (`compare_all_fields.py`, new field ignored).

**Descriptive (reported, not judged)**

- `Call_rank` for samples that fail C1's precondition (fuzzy, retry, masked).
- The rank of the true genotype in the 25% set and the damaged NA19920, wherever
  the call is wrong or a `;` set.
- The distribution of the rank-2 `Delta_log10_likelihood` across the 53.
- Wall time per sample against master. PREDICTED (not run): the largest table
  dominates.

## Limits

- The ranking can only be as good as the table. A genotype that is not in the
  table (a novel allele) cannot appear, and the best table entry still ranks 1.
- A large `Delta_log10_likelihood` means the reads at the scored sites favour
  rank 1 under this model. It does not mean the call is right. Errors that are
  not noise (missing reads at a site, contamination, a wrong copy number) can
  give a large gap to a wrong call.
- Rule-determined tokens show where the call rests on a rule. They do not
  measure how reliable the rule is.

## Follow-up

The support measure (site jackknife or read bootstrap) builds on this
decomposition. Before any support number is reported, it must separate known
right from known wrong calls on the benchmark, the thinned sets, and the damaged
files.

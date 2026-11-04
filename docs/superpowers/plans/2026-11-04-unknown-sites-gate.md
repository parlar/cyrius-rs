# Unknown Target Sites in Allele Matching — Gate Record

> This is a **gate record**, not an implementation plan. It locks the cheap
> experiments that can kill the idea and their decision rules **before** they
> run. The implementation plan is written only if the idea survives them.
> Commit order: this file (`plan:`) → throwaway experiment code → `result:`.

**The defect (review #1, measured):** an uncallable target site is treated as
"no variant here". On NA19920 with the 68 alignments over the `*4` splice site
removed, the call moves from `*1/*4x2` (activity 1.0, **Intermediate
Metabolizer**) to `*1/*10x2` at 0.95 HIGH, `PASS`, activity 1.5, **Normal
Metabolizer**. Measured output: `~/cyrius-research/review-fixes/na19920_missing_splice/out/dmg.json`.
Since `e035eb4` the five unknown sites are listed in `Uncallable_target_sites`,
but the call is unchanged.

**The direction (agreed during the P1 review):** keep the uncertainty inside
allele matching. An unknown site matches any copy count; every compatible
candidate is kept; a clinical result is reported only when all candidates
agree; no-call is the fallback.

---

## Measured facts this record rests on

| Fact | How measured |
|---|---|
| The matcher is an exact key lookup: observed variants are sorted, joined with `_`, and looked up in the CNV group's table. A missing site is simply absent from the key. | Read `get_star`, `match_star_allele.rs:58-124`. |
| Reporting already exists for multi-candidate calls. `clinical_equivalence` returns `same_effect` only when every candidate has the same activity score **and** phenotype, and `unknown_effect` if any is `n/a`. | Read `phenotype.rs:212-225`. |
| The resolver runs on any `;`-joined candidate set and may pick one. | `docs/superpowers/plans/2026-06-04-diplotype-resolver.md`, line 7. |
| Routing of the 53 benchmark samples: **45** use the main `get_star` path (cn2 25, exon9hyb 7, cn3 6, exon9hyb_exon9hyb 5, dup_star13intron1 1, dup_exon9hyb 1); **4** the single-copy path (star5 3, star5_star68 1); **3** the `*68` path (star68); **1** a fixed genotype (star5_star5). | `CNV_group` counted over `~/cyrius-research/review-fixes/accepted/*.json`, against `cnvtag_to_genotype`, `uses_cn1_dic` and the `star68` branch in `match_star`. |
| **0 of 53** benchmark samples have any uncallable site. | `Uncallable_target_sites` non-empty on 0/53 (e035eb4 run; all 33 fields unchanged through f21e8fe). |
| Damaged NA19920 is in group `cn3`, main path, `unique_match` `*1_*10.002_*10.002`, with 5 unknown sites: the two `g.42128936-42128937` insertions, `g.42128932C>T`, `g.42128933G>A`, `g.42128945C>T`. Three have zero reads; `g.42128932C>T` has 4 and `g.42128933G>A` 3, too few at total CN 5. | `dmg.json`; read counts from the P1 measurement. |
| Everything an offline experiment needs is public: `get_hap_table`, `combination_table_for`, `get_final_call_clean`, `match_phenotype`, `clinical_equivalence`. | `lib.rs`, `caller/mod.rs`, `construct_star_table.rs:111,144`. |
| The project's case file exists as a draft outside the repo, `~/cyrius-research/judgment-gate-cases.draft.md`. Relevant entries: **B7/B8** (a gate with no power, a filter that cannot fire — the shape of the vacuous kill test below), **C3/C4** (results measured only on clean input; rows no sample reaches). | Read. |
| Case **C4**'s off-by-one target `g.42128445delCTGCC` ("0,0 in all 62") is fixed: it now has reference reads in 52 of 53 samples. The one `0,0` sample is NA19317 (`*5/*5`, 117 zero-read sites), where copy number 0 is the correct reading, and it takes the fixed-genotype path. So on the benchmark, zero reads occur only where there is no gene. | `Variant_raw_count` over `accepted/*.json`. |
| The largest combination table has 1,551,523 keys. | **Code comment** at `match_star_allele.rs` (~L886), not re-measured. B1 will measure it. |

---

## Gate A — revised

**1. Governing principle.** Evidence that was never observed must not narrow
a call. Unknown is a third state, and it reaches the clinical answer only when
it cannot change it.

**2. What would kill it.** The kill test written in the P1 plan was "on fully
observed inputs the candidate set equals today's call on all 53". That test is
**vacuous**: 0 of 53 samples have an unknown site, so masking hides nothing
and the test passes by construction. It is kept as a regression invariant,
not as a kill test. The real kill tests are K1, K3 and K4 below.

**3. Already refuted here?** No. A search of the git log and the design docs
for masking, unknown dosage and candidate sets finds no earlier attempt. The
resolver work uses the word "candidate set" for `;`-ambiguous calls, which is
a different thing and is relevant only through K4.

**4. Simplest thing that could work.** The trivial baseline is: **no-call
whenever any target site is unknown.** It fixes the measured failure with a
few lines. The candidate-set mechanism is compared with it **case by case,
not by a rate threshold** (user decision, 2026-11-04). In every case the
candidate set either establishes a phenotype or reports the alternatives; the
baseline always reports nothing. So the candidate set is never worse than the
baseline **unless it establishes a wrong phenotype**. That one outcome is the
kill test (K3). How often it establishes an answer is measured and reported,
but decides nothing.

**Never-wrong by design, not by hope.** Masking only removes constraints, so
every key the evidence allows stays compatible, the true one included. An
established phenotype can then be wrong only if a step **after** the lookup
drops a candidate without evidence. One such step exists: `get_star`'s
rare-allele rule, which keeps only common alleles when there is a choice
(`*34`, `*39`, `*4.009`, `*139`). That is a frequency prior, and priors were
measured to entrench wrong calls here (H19, `~/cyrius-research/REPORT.md`).
**When any site is unknown, the rare-allele rule is not applied.** The
simulation then tests the implementation, not the principle.

The guarantee is conditional, and the condition is stated: the true genotype
must be in the table, and every site the matcher **does** trust must be
right. Masking protects only against sites flagged unknown. A site that
damage makes confidently wrong, instead of unknown, is not caught. K1 tests
that condition on the one real damaged sample.

**5. What must be true of the inputs.**

| Assumption | Status |
|---|---|
| A site in `Uncallable_target_sites` was not observed, rather than observed and rejected. | **Partly false, known.** `call_cn_var` reaches `None` by more than one route: zero reads, zero probability mass, posterior below cutoff (`process_raw_call_denovo`), and the `NOISY_VAR` strand filter. The posterior route leaves two copy-number values open, so a full mask is wider than the evidence. B1 measures both a full mask and a mask limited to the two open values on NA19920. |
| `Variants_called` in the output is exactly what `get_star` saw. | True for the main path: it is `var_observed` after `update_variants`, inside `match_star` (`main.rs:1046`, `match_star_allele.rs`). Not true for the `*68` path, which edits its own copy. |
| A target site listed as uncallable has no other caller that resolved it. | **False for at least one site.** `g.42130655-42130656insA` has a dedicated caller. Such sites must be excluded from the mask when that caller decided. Not relevant to NA19920. |

---

## Gate B — the experiments, with rules locked before running

### B1. Offline: lookup and dropout simulation (no pipeline change)

A throwaway example program reads the tables and the existing benchmark
outputs. It changes nothing in `src/`.

**Masked lookup.** For an observed multiset `O` and an unknown-site set `U`,
a table key `K` is compatible when, for every variant `v` not in `U`, the
count of `v` in `K` equals its count in `O`. The candidates are the union of
the table entries of every compatible key, followed by the same
post-processing `get_star` applies **except the rare-allele rule** (see
Gate A #4), then `get_final_call_clean`, `match_phenotype` and
`clinical_equivalence`. The run also reports how many cases the rare-allele
rule *would* have changed, so its effect is visible rather than assumed.

**On NA19920 (damaged).** Report the candidate set, each candidate's activity
and phenotype, and the equivalence verdict, with a full mask and with the
two-value mask.

**Dropout simulation.** For each of the 45 main-path samples, and for each
target-list site, one case: remove every copy of that site's variant from
`Variants_called` (if present) and mark it unknown. Classify each case as
**established-correct**, **established-wrong**, **unresolved**, or **empty**
(no compatible key). "Correct" means the established phenotype equals the
phenotype of the truth genotype, computed with the same functions. Report
cases where the dropped variant was present (A) and absent (B) separately,
and in total. Also report what today's exact lookup does in the A cases:
right, wrong, or no match. That is the size of the problem.

The 8 samples on other paths are excluded and listed by name in the output.
No silent cap.

**What trivially passes these, and why it does not here.** The truth's own
key is always compatible, because masking only removes constraints. So
"truth is in the set" cannot fail on a sample today's matcher already gets
right. It **can** fail on damaged NA19920, because more than the masked sites
might have changed. "Established-wrong = 0" can also fail: the rare-allele
rule and `get_final_call_clean` can drop a candidate after the lookup.

### B2. End-to-end on damaged NA19920 (throwaway pipeline hack)

A scratch branch adds the masked lookup to the main path behind an
environment variable, and runs the damaged BAM through the whole pipeline,
resolver included. The branch is deleted afterwards and never merged.

---

## Locked decision rules

| Rule | Kills the idea when | Then |
|---|---|---|
| **K1** | The truth `*1/*4x2` is **not** among the NA19920 candidates, under the full mask. | Masking alone cannot fix the measured failure. Stop; ship the no-call baseline. |
| **K3** | **Any** dropout case is established-wrong. Zero tolerance. | Masking can produce a confident wrong answer, so the never-wrong design does not hold. Stop and diagnose before anything else. |
| **K4** | End to end, damaged NA19920 is reported with a **single phenotype other than Intermediate** — e.g. the resolver collapses the set back to `*1/*10x2`. | The resolver undoes the masking. Scope grows to the resolver; re-gate before building. |

**There is no rate threshold.** An earlier draft had one (K2: established-
correct under 25 % kills the idea). It was dropped before anything ran: the
case-by-case comparison above makes the candidate set never worse than the
no-call baseline whenever K3 holds, so a rate measures how much it helps, not
whether it is safe. The rate is reported as **established-correct /
unresolved / empty**, for present and absent dropped variants separately.

**Single-site dropout is the mildest kind of coverage loss.** The real
NA19920 damage hid five sites at once. The implementation plan must also
measure multi-site loss before any claim about how often it answers.

**If K1, K3 and K4 all pass:** write the implementation plan. It must include the
regression invariant (no unknown sites → exactly today's call on all 53, and
the masked lookup with an empty mask returns exactly the table entry for
every key of every table), the reviewer's five required tests, and decisions
on the single-copy and `*68` paths, the empty-set fallback, and the two-value
mask.

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

---

## Addendum (locked before running): does the defect occur in thinner real data?

**Why.** The user asked whether validating the BAM up front would beat
handling the consequences. For whole-file problems it would, and several
such checks exist. But the damaged NA19920 is a valid file with a coverage
hole, and every hole tested in B1 was made by hand. On the 53 real files,
none has an uncallable site. So whether this defect occurs in real data at
all is **unmeasured**. This addendum measures it before B2.

**Input.** Each of the 53 truth samples' BAMs, thinned with
`samtools view -s 42.25` and `samtools view -s 42.10` (seed 42; 25 % and
10 % of read names, so mates stay together). Both levels are fixed here,
before any run.

**Run.** The `a1b919b` binary (its `src/` is identical to this branch), with
the benchmark's flags: `--strand-bias-all --fuzzy-match --changepoint-hybrid --het-check`.

**Per run, measured:** genotype, filter, `Uncallable_target_sites`, and the
phenotype of the call against the phenotype of the truth genotype (same
functions as B1). A call is **confidently wrong** when it establishes one
phenotype (all its alternatives agree) and that phenotype is not the
truth's.

**Decision rule — no tuned numbers:**

| Rule | Outcome | Then |
|---|---|---|
| **S1** | **No** run, at either level, is confidently wrong while listing at least one uncallable site. | The defect does not arise even in thin data. Stop the masking plan; keep the reporting field; rely on up-front checks. |
| **S2** | At least one such run exists. | The defect is real in realistic data. Continue to B2. |
| **K3 on real holes** | On every S2 run on the main path, masking its own uncallable sites establishes a wrong phenotype. Zero tolerance. | Stop and diagnose before B2. |

**Reported, deciding nothing:** wrong calls that list **no** uncallable
site (other causes, such as copy number at low depth); how many of the S2
runs carry `PASS`; and how often the truth falls out of the masked
candidate set, meaning a site the matcher still trusted was itself wrong.
That last number is the limit of masking.

---

## Results so far (B1 and the thinned-data addendum; B2 not yet run)

Raw outputs: `~/cyrius-research/unknown-sites/b1_output.txt`,
`real_full.txt`, `real_p25.txt`, `real_p10.txt`. Thinned BAMs and their
caller outputs: `~/cyrius-research/unknown-sites/thinned/`.

**Controls (all passed before any result was read).** Exact lookup finds the
same raw candidates as the pipeline on 45/45 main-path samples; the formatted
call matches the reported genotype on 43/43 (2 changed by fragment phasing,
skipped). Masking by enumeration and by full scan agree on NA19920 (9 keys
each). The real-holes mode reports 0 holes and 0 confidently wrong calls on
the 53 full-depth outputs (36 calls established, so the comparison ran), and
it **does** flag the damaged NA19920 as confidently wrong at `PASS`.

**K1 — PASS.** Damaged NA19920, 5 unknown sites masked: 9 candidates
(`*1/*4x2` among them), activity 0.0 to 1.5 plus several of unknown
function. Verdict: unresolved. All 9 differ only at the two zero-read sites
(the 9-base insertion and the splice site); the two low-read sites are 0 in
every candidate, so the two-value mask would change nothing here.

**K3 — PASS, and not vacuously.** 7,155 single-site dropouts on the 45
main-path samples: 0 established-wrong. In 167 of them masking **added**
candidates and still established a phenotype, all correct. Today's exact
lookup on the 197 dropouts of a present variant: 28 right, **17 confidently
wrong**, 16 unresolved, 124 no exact match. Restricted to the 30 samples
whose unmasked call establishes a phenotype: present-variant dropouts 98
correct / 37 unresolved / 0 wrong; absent-variant dropouts 3,514 / 1,121 / 0.
15 of 45 samples never establish a phenotype even unmasked: each carries an
allele with no function row.

**Skipping the rare-allele rule was not shown to be necessary.** It would
have changed 66 cases and made 0 of them wrong. It stays skipped on
principle (H19); this is caution, not a measured need.

**Addendum — S2, so continue.** At 25 % of reads, 12 of 53 samples get a
genotype; 41 are refused at the copy-number step (37 with no CNV tag, 4 with
no total CN). 7 runs list uncallable sites. **One** is confidently wrong with
holes: NA21105, called `*2/*3` (activity 1.0, Intermediate) against truth
`*3/*111` (Indeterminate: `*111` has no function row), 5 holes, filter
`Fuzzy_match`, **not PASS**. Masking its 5 holes gives
`*111/*3;*2/*3;*3/*30`, unresolved, truth included. K3 on the 6 real-hole
runs masking can handle: 0 established-wrong. Truth fell out of the masked
set twice: NA18973 (`*2/*21` for truth `*1/*21`; a trusted site was wrong,
phenotype still right) and NA19908 (no compatible key, so the fallback). At
10 %, 2 of 53 get a genotype and neither is wrong. Confidently wrong with
**no** holes: 0 at either level.

**Honest reading.** The defect is real in thinned data but rare, and the one
real instance was already flagged `Fuzzy_match`. The only `PASS` instance is
still the hand-made NA19920. Thin WGS is mostly refused before matching.
Data where copy number succeeds but single sites are missing (for example
capture gaps) has not been measured.

**Found along the way, must be fixed in any implementation:**
`get_final_call_clean` prints raw multisets for any non-cn2 group with more
than one candidate. It occurs on real data at 25 %: HG02373 is printed
`*10_*14_*36;*10_*36_*65;*10_*114_*36`, so its phenotype reads as
Indeterminate although the truth `*14/*36+*10` is the first candidate.

---

## B2 result — K4 PASS. All kill tests passed.

Prototype: `~/cyrius-research/unknown-sites/b2_prototype.patch` (applies cleanly
to `da0b66a`; the scratch branch was deleted unmerged, as locked). It masks the
uncallable sites in the main matching path when `CYRIUS_MASK_PROTO` is set,
skips the rare-allele rule, and formats each candidate on its own. Outputs:
`~/cyrius-research/unknown-sites/b2/`.

**K4 — PASS.** Damaged NA19920, prototype on: all 9 candidates are reported;
the resolver scores them exactly equal (log10 likelihood -727.49 each, 0
decisive reads) and leaves the set ambiguous. Filter `PASS` ->
`More_than_one_possible_genotype`, confidence 0.95 HIGH -> 0.80 MEDIUM,
`Clinical_equivalence` `unknown_effect`, `Long_read_recommended` true. The
phenotype column lists every option. Control, same binary with the prototype
off: `*1/*10x2`, Normal, `PASS`, 0.95, unchanged from before.

**Invariant held.** All 53 full-depth benchmark files with the prototype on: 0
differing field values against `accepted/` (every field except `Tool_commit`
and `Flags`). Prototype off on the 25 % set: 0 differences against `a1b919b`.

**On real thinned holes, prototype on (25 %): 2 of 53 change, both better.**
- NA21105: `*2/*3` [`Fuzzy_match`] -> `*3/*111` [`PASS`], picked by the
  resolver from 3 masked candidates with 5 decisive reads (log10 likelihood
  -185.00 vs -203.46). The 5 reads carry `g.42129759C>T`, `*111`'s defining
  variant, which had 5 alt / 3 ref reads: too few for the dosage step, so today
  it is treated as absent. The result equals the truth.
- NA20289: `*6/*11` unchanged; filter `Fuzzy_match` -> `LowQ_coverage_MAD`. The
  call no longer needs fuzzy matching, which exposes a coverage warning the
  fuzzy filter had been hiding. Still not `PASS`; genotype equals the truth.

At 10 %: no change on any sample.

**Caveat carried forward.** NA21105's `PASS` rests on the resolver, whose
model and posterior are not calibrated (review #8). Here the evidence is five
reads carrying the allele's own defining variant, which is direct; but the
implementation plan should not let a resolver pick on a masked set reach
`PASS` without a stated rule for it.

**Next:** the implementation plan.

## Result: implementation (Task 5 of the implementation plan) — measured at 176cb35

Runs: `~/cyrius-research/unknown-sites/final/run.sh` (same flags as the benchmark),
outputs in `~/cyrius-research/unknown-sites/final/{dmg,p25,p10}`, audit with
`mask_audit --real` in `final/real_p25.txt`, `final/real_p10.txt`, field diffs
against the `a1b919b` thinned outputs in `final/diff_p25.txt`, `final/diff_p10.txt`.

| Input | Locked expectation | Measured | |
|---|---|---|---|
| 53 full-depth samples | 0 differing field-values | gates 1, 2, 3, 4 and 4b: 0 of 53 each | met |
| Damaged NA19920 | 9 diplotypes as unordered pairs; `More_than_one_possible_genotype`; `unknown_effect`; long-read true; 0.80 MEDIUM | `*1/*10x2;*1/*147x2;*1/*4x2;*145/*147x2;*4x2/*145;*4/*4x2;*10x2/*145;*4/*10x2;*4/*147x2`, same 9 pairs; `More_than_one_possible_genotype`; `unknown_effect`; true; 0.80 MEDIUM; resolver posterior 0.111 (a 9-way tie) | met |
| NA21105, 25 % | `*3/*111`, `PASS`, resolver | `*3/*111`, `PASS`, `Resolved_by` resolver (was `*2/*3` `Fuzzy_match`) | met |
| NA20289, 25 % | `*6/*11`, `LowQ_coverage_MAD` | `*6/*11`, `LowQ_coverage_MAD` (was `Fuzzy_match`) | met |
| HG02373, 25 % | no `_`; `*14/*36+*10` among options | `*14/*36+*10;*36+*10/*65;*36+*10/*114` | met |
| NA19908, 25 % (D1 case) | unchanged | unchanged: `*1/*46;*43/*45;*12/*43`, `More_than_one_possible_genotype` | met |
| HG03259, 25 % (single-copy path) | PREDICTED (not run) | unchanged: `*5/*106`, `PASS` (1 unknown site; masked lookup on the single-copy table still finds only `*106`). At 10 %: unchanged, `LowQ_coverage_MAD`, 5 unknown sites | recorded |

**Changed samples.** 25 %: 3 of 53 (NA21105, NA20289, HG02373), 11 field-values,
every one in the table above. 10 %: 0 of 53.

**Audit (`mask_audit --real`).** 25 %: confidently wrong with holes 1 -> **0**
(NA21105 fixed); confidently wrong without holes 0 -> 0; ESTABLISHED-WRONG 0.
10 %: 0 -> 0 throughout.

**Stop rules:** none triggered (no full-depth difference; no run became
confidently wrong; the NA19920 set equals the locked 9).

**Not exercised by any run.** The new filter value `Unknown_sites_unmatched`:
NA19908 is the only D1 case and its ambiguity filter comes first. Pinned instead
by the unit test `an_unmatched_sample_keeps_its_fallback_answer_and_is_flagged`
(NA19908's measured inputs).

**Open, carried to the user.** The gate record's B2 caveat (a resolver pick on a
masked set reaching `PASS` needs a stated rule) was not addressed by the
implementation plan; NA21105 reaches `PASS` that way. See the final report.

**Found while measuring (not in the plan).**
- *Tie order is not reproducible.* The same binary (176cb35) run 6 times on the
  damaged NA19920 gave the same order in runs 1-5 and a different order in run
  6; the set is the same 9 diplotypes, `Confidence` 0.80 and `Resolver_posterior`
  0.11111111111112557 every time (raw runs: `~/cyrius-research/unknown-sites/tie_runs/r1..r6/`).
  The earlier comparison was across two builds (0cf2b80 and 176cb35, a
  behaviour-neutral refactor), not one binary. The resolver sorts by posterior,
  then by name (`resolver.rs` ~L38), but the 9 posteriors are equal only to
  about 1e-13 and differ in the last bits between runs, so the name tie-break
  never fires. Only the order changes: the set, the filter and the
  ambiguous/resolved decision (margin 0.10) do not. `confidence.rs` (~L203)
  scores completeness from the first `;` alternative only, so `Confidence`
  could change with the order (not observed). The resolver behaved this way
  before this branch; masking makes many-way ties reachable. 25 % and 10 %
  sets: identical between the two builds.
- *Size of a masked set with many unknown sites* (probe on the real tables,
  empty variant list): 2-copy table 14 / 35 / 70 / 140 unknown sites ->
  15 / 136 / 378 / 10,153 diplotypes; 3-copy table -> 25 / 256 / 729 / 20,164
  (a 219 KB genotype string). Scan <= 0.11 s, formatting <= 0.25 s. Nothing
  downstream caps the set; resolver time on such a set is not measured. Not
  reached in any run: real thinned files had at most 5 unknown sites, and
  NA19920 with reads removed at every 10th or every 4th target site (14 or 35
  sites, `~/cyrius-research/unknown-sites/many_holes/`) is refused at the
  copy-number step (no call, ~1 s).

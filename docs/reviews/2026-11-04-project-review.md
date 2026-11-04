# Project review, 4 November 2026

Reviewed code revision: `3807326`. The active Rust caller is the main scope; the Python implementation, benchmark harness, embedded resources, design notes, methods inventory, and README were used for comparison. Production source files were not edited. Local README edits and the existing untracked presentation PDF were left intact.

The Cyrius approach is a reasonable foundation for short-read CYP2D6 calling. Summing D6/D7 depth, using paralog-differentiating sites, and preserving structural ambiguity are useful choices for this locus. The immediate problem is the handling of uncertainty: missing observations, conflicting observations, and technical failures can become confident genotype evidence. These problems should be fixed before adding more inference modules or presenting confidence as a measure of accuracy.

## Verification performed

| Check | Result |
|---|---|
| `cargo test --all-targets --locked` | 179 tests passed: 161 library, 12 binary, 6 integration |
| `cargo build --release --locked` | Passed; benchmark JSON records identify commit `3807326` |
| Fresh default benchmark | 51 exact, 2 ambiguous, 0 wrong, 0 no-call, across 53 truth samples |
| Fresh benchmark with `--strand-bias-all --fuzzy-match --changepoint-hybrid --het-check` | Same counts; genotype, filter, CNV group, confidence, label, and resolution method identical for all 53 samples |
| `cargo clippy --all-targets --locked -- -D warnings` | Failed on existing lint findings; normal build and tests pass |
| Targeted coverage removal in a real benchmark sample | Produced an incorrect PASS/HIGH call and changed the predicted phenotype |
| Synthetic fragment, activity-score, and corrupt-BAM probes | Reproduced the defects described below |

The two unresolved benchmark calls are NA18565 and NA19908. Fragment phasing resolves two other samples. Neither remaining ambiguity is resolved by the likelihood resolver; both have a top conditional posterior approximately 0.5. This benchmark therefore does not validate a successful real-data decision by the likelihood resolver.

Evidence from this review is retained locally in `/tmp/cyrius_review_29pdv5po/`. It contains `probe.rs`, `probe_output.txt`, synthetic BAMs, the modified NA19920 BAM and output, and the corrupt-BAM probe. Fresh benchmark results are in `/tmp/bench_9DHkn2/` and `/tmp/bench_4JDUgZ/`. These temporary directories are not permanent test fixtures.

## Confirmed correctness defects

### 1. [P1] Uncallable target sites become reference evidence

**Status:** fixed on `feat/unknown-sites` (commits `ad4f123` formatting, `cec2334` masked lookup, `4b64bbb` `*68` paths, `0cf2b80` + `176cb35` wiring and the D1-a test, `841303e` results). Uncallable sites are now matched as unknown. Damaged NA19920 now gives 9 possible diplotypes including the truth `*1/*4x2`, `More_than_one_possible_genotype`, confidence 0.80 MEDIUM, instead of `*1/*10x2` `PASS` 0.95 HIGH; no field changes on the 53 full-depth samples. The original finding follows.

Locations: [call_variants.rs](../../cyrius-rs/src/caller/call_variants.rs), `get_called_variants`, around lines 575–593; [main.rs](../../cyrius-rs/src/main.rs), variant-set construction around line 817; [confidence.rs](../../cyrius-rs/src/caller/confidence.rs), around lines 237–250.

`call_cn_var` can return `None` for a site without sufficient evidence. `get_called_variants` discards this state. The matcher receives the same variant multiset for a missing site and a confidently reference site. Consequently, a unique dictionary match can be reported even when a variant needed to distinguish that genotype was never observed.

**Real-data reproduction:** From `NA19920_pgx_regions.bam`, remove only alignments on chr22 overlapping zero-based position 42128944, the main `g.42128945C>T` splice variant defining *4. This removes 68 alignments. Run the freshly built binary with default flags.

| Input | Genotype | Filter | Confidence | Activity | Phenotype |
|---|---|---|---|---|---|
| Original | `*1/*4x2` | PASS | 0.87 HIGH | 1.0 | Intermediate |
| Splice-site coverage removed | `*1/*10x2` | PASS | 0.95 HIGH | 1.5 | Normal |

The modified output explicitly contains `0,0` at the splice site. It reports `Phasing_ambiguous=false` and `Long_read_recommended=false`. The wrong answer receives a higher confidence score than the original answer.

A separate function-level probe makes all 140 non-homology target sites uncallable. The matcher returns `*1/*1`, `unique_match`. Given otherwise ideal depth and diagnostic-SNP inputs, the confidence code returns 1.0 HIGH. This is a component probe, not a claim that a completely empty BAM passes the full pipeline.

**Required change:** Retain per-site states for reference, alternate dosage, uncertain dosage, and unavailable evidence. Evaluate candidates against these states or site likelihoods. An absent defining variant may exclude an allele only when there is adequate reference evidence. Missing differentiating sites should preserve alternatives or produce an explicit insufficient-evidence result. Wild-type calls need reference callability across the supported target panel.

### 2. [P1] Multiplication activity lookup changes the metabolizer phenotype

**Status (fix/review-p1):** fixed in 146ad3d; test `an_exact_copy_count_multiplies_the_single_copy_activity`.

Location: [phenotype.rs](../../cyrius-rs/src/phenotype.rs), lines 15–16 and 67–70.

The loader removes `≥` from names and activity values. The fallback then uses any `*Ax3` row for all copy numbers of three or more. The embedded `*41x3` row is an exact three-copy row, with activity 0.75; it is not a `*41x≥3` bucket.

Reproduced through `match_phenotype` using the embedded table:

| Diplotype | Current result | Expected using the embedded single-copy activities |
|---|---|---|
| `*41x4/*10` | 1.0, Intermediate | 1.25, Normal |
| `*41x5/*5` | 0.75, Intermediate | 1.25, Normal |
| `*2x4/*5` | 3.0, Ultrarapid | 4.0, Ultrarapid |

The *41 activity value of 0.25 is also documented in the [ClinPGx function update](https://blog.clinpgx.org/cyp2d6-allele-function-update/). CPIC classifies activity 1.25 as Normal and activity 1.0 as Intermediate; see the [CPIC beta-blocker guideline](https://pmc.ncbi.nlm.nih.gov/articles/PMC11502236/).

**Required change:** Preserve exact counts and lower bounds as different data types. For an exact multiplication, calculate the activity from the exact copy count and single-copy value unless a specific validated exception applies. Preserve bounds when the copy count itself is unknown. Add boundary tests that check phenotype as well as numeric activity.

### 3. [P1] Disagreeing overlapping mates can reverse the phasing result when reordered

**Status (fix/review-p1):** fixed in 9593eb0; test `disagreeing_mates_give_the_same_result_in_either_order`.

Location: [phase_disambiguate.rs](../../cyrius-rs/src/caller/phase_disambiguate.rs), `observe_read`, lines 241–251.

Observations are grouped by fragment name, but a later nonempty observation overwrites the earlier allele at each site. Conflicting observations are not reconciled or marked uncertain.

**Reproduction:** Three fragments, each represented by two overlapping, high-quality primary mate records. One mate observes ALT/ALT at two sites; the other observes ALT/REF. With the exact same records in one order, `disambiguate` returns the trans candidate `*2/*3`. Reversing their order within each equal-coordinate pair returns the cis candidate `*1/*4`.

This is a synthetic counterexample, not an observed error in the 53 baseline samples. It demonstrates that valid changes to record order can change a biological conclusion.

**Required change:** Merge mate observations using base qualities and an explicit conflict policy. Equal-quality disagreement should become uninformative at that site. A low-quality observation must not overwrite a stronger one. Deduplicate overlapping observations and test invariance to mate order and equal-coordinate BAM record order.

### 4. [P1] Low-quality bases can supply decisive phasing evidence

**Status (fix/review-p1):** fixed in 323a8c7; test `base_quality_zero_reads_do_not_decide_phase`. The main allele-counting path (`snp_count.rs`) is unchanged; see Out of scope. Omitted qualities (`0xFF`) leave phase unresolved. The short-quality-array guards at both call sites are unreachable rather than live rejection conditions: `Record::qual()` returns exactly `seq_len()` bytes (rust-htslib 0.47.1, `record.rs:582-585`) and both sites already bound their slice by `seq.len()`. They are kept so the slice cannot panic, and are labelled as such in the code.

Locations: [phase_disambiguate.rs](../../cyrius-rs/src/caller/phase_disambiguate.rs), `phase_read_passes` and `check_allele_pileup`, around lines 213 and 456; [snp_count.rs](../../cyrius-rs/src/depth_calling/snp_count.rs), around lines 117–146.

The phasing path tests sequence patterns and some record flags, but does not inspect the base qualities of the supporting observations. The main allele-counting path similarly checks mapping quality and sequence matches without using base quality. Mapping quality and base quality represent different uncertainties.

**Reproduction:** Three distinct synthetic fragments with MAPQ 60, ALT at both differentiating sites, and base quality zero throughout are sufficient for `disambiguate` to select `*1/*4` over `*2/*3`.

**Required change:** Use base-quality-aware likelihoods or a validated quality threshold at all allele-discriminating bases. Multi-base allele patterns need a policy for the quality of the entire informative pattern. Keep the deliberate treatment of low-MAPQ paralog reads separate: discarding every low-MAPQ read would damage the core Cyrius dosage method.

### 5. [P1] The targeted phase rules interpret “not demonstrated” as the opposite phase

**Status (fix/review-p1):** partly fixed in 79cc8ff (plus 3705503, a clippy type alias): the unsupported negative branch was removed, not replaced; test `phase_signal_42126938_only_ever_says_yes` and `unknown_phase_keeps_both_diplotypes`. `phase_signal_42126938` can no longer return `Some(false)`; no affirmative non-G pattern is established. Open question, no plan.

Locations: [match_star_allele.rs](../../cyrius-rs/src/caller/match_star_allele.rs), `apply_cn2_phase_rules`, lines 755–769; [call_variants.rs](../../cyrius-rs/src/caller/call_variants.rs), `call_var42126938` and `call_var42127803hap`.

The specialized haplotype functions initialize their signals to `false`. Too few linked reads leaves the signal false, as does an absence of evidence for the target arrangement. The matcher nevertheless excludes one candidate on either boolean value.

For example, observations matching the *41 variant set produce raw candidates `*1_*41` and `*119_*2.001`. Passing the default false signal makes `match_star` report only `*1/*41`; passing true reports only `*119/*2`. There is no third state for unavailable or inconclusive phase. These rules execute before the general fragment stage and resolver, including when the general phasing feature is disabled.

**Required change:** Return a typed result such as cis, trans, and unknown, with fragment counts and a support criterion. Eliminate a candidate only with affirmative evidence. Keep both candidates when phase is unknown. Report these specialized resolutions in provenance as well.

### 6. [P1] BAM decoding failures are accepted as successful partial counts

**Status (fix/review-p1):** fixed in 158dbdd; test `a_corrupt_block_is_an_error_not_a_short_count`. Related, not one of the six: `Phase_pairs_truncated` (2968087, extended to an unopenable BAM in 832acd8) reports when a read error cut the pair list short.

Locations: [bin_count.rs](../../cyrius-rs/src/depth_calling/bin_count.rs), lines 23–25; [snp_count.rs](../../cyrius-rs/src/depth_calling/snp_count.rs), pileup error handling around lines 105–109; analogous record loops in the haplotype, phasing, insertion, and likelihood paths.

Fetch failures are now handled, but failures encountered while reading fetched records are often skipped with `continue`. An error after reading a prefix of a region therefore produces a seemingly valid count for that prefix.

**Reproduction:** Write an indexed synthetic BAM containing 2,000 mapped reads, then invalidate the CRC of one later BGZF data block without changing offsets or the index. `get_read_count` returns `Ok(242)` while htslib reports a CRC32 mismatch and failed block read.

**Required change:** Propagate decoding and pileup errors. Treat mandatory evidence as unavailable and return a structured technical failure. Optional scoring should stop and preserve the prior ambiguity. Do not score a candidate or normalize depth using an incomplete region. A successful indexed fetch does not establish successful reading of its contents.

## Methodological and validation findings

### 7. [P2] Confidence is a heuristic quality index, not calibrated call accuracy

Location: [confidence.rs](../../cyrius-rs/src/caller/confidence.rs), especially lines 66–88 and 188–286.

The weighted geometric mean, component weights, and HIGH/MEDIUM cutoffs have no demonstrated out-of-sample calibration. Several components summarize the same reads or the consequences of the same matching decision. They are not independent confirmation channels.

Specific weaknesses compound the missing-data defect:

- Expected-variant completeness requires only one ALT observation, regardless of copy number or uncertainty.
- `*1/*1` returns perfect completeness and specificity immediately, including when target counts are all missing.
- Expected variants come from main-allele definitions, although exact matched suballeles can have different variants.
- Only the first reported diplotype is evaluated, so an ambiguous set is not assessed symmetrically.
- The score does not estimate the correctness of structural assignment or phenotype separately.

Retain the component metrics, make callability and mandatory QC prerequisites explicit, and label the aggregate as an uncalibrated quality score until it is evaluated on independent samples. Calibrate genotype, structural, and phenotype error separately, using reliability plots and error versus call-rate curves. A complex calibration model fitted to 53 mostly correct training calls would give little trustworthy information about rare failures.

### 8. [P2] The resolver’s generative assumptions need correction and independent validation

Locations: [diplotype_caller.rs](../../cyrius-rs/src/caller/diplotype_caller.rs), around lines 432–471, 587–616, and 802–833.

The resolver is a useful architectural step, but its output should remain a conditional model score:

- It adds exactly one D7 copy to a mixture whose other entries are diploid-sample D6 copy counts. Ordinary diploid D7 dosage is two, and hybrid structures require position-dependent sequence contributions. Filtering on MAPQ further changes effective sampling proportions, so replacing one constant with two alone would not complete the model.
- At diagnostic-only sites, every D6 allele expects D6-reference sequence. This does not represent the D7-derived segments of hybrid alleles.
- Overlapping mate observations are appended and multiplied without a consensus policy, which can exaggerate evidence from one fragment.
- The implementation takes the maximum likelihood over a label’s matched raw suballeles, then normalizes those maxima as posteriors. This is a profile likelihood, not marginalization with explicit suballele priors.
- The candidate set has already been restricted by the upstream CNV choice, exact/fuzzy matching, rare-allele rules, and targeted phase rules. Normalization cannot measure uncertainty about candidates that were removed or never generated.

The decisive-fragment safeguards are useful, but do not validate these assumptions. In the fresh benchmark no sample is resolved by this stage. Validate it on synthetic informative examples and independent real samples where different arrangements genuinely produce distinguishable fragments before giving it broader authority.

### 9. [P2] The benchmark harness cannot serve as a regression gate

Location: [bench/run_benchmark.sh](../../bench/run_benchmark.sh), lines 39–42 and 117–129.

An existing release binary is reused without rebuilding. The script can therefore test stale code. It also returns failure only for execution errors, allowing wrong calls, no-calls, and missing BAMs to pass. Pointing it at a directory without BAMs reproduced `0 / 0` exact, 53 missing BAMs, and exit status zero.

Add a strict regression mode that rebuilds the intended revision, records input/resource provenance, checks the expected sample count, and fails on unapproved per-sample result changes. Explicitly permitted ambiguous calls should be versioned, rather than treating every ambiguity or every truth-containing candidate set as equivalent. Exploratory scoring can remain available separately.

### 10. [P2] The proposed long-read reflex description overstates what PCR read counts establish

Locations: [README.md](../../README.md), long-read proposal around lines 377–385; [main.rs](../../cyrius-rs/src/main.rs), reflex decision around lines 1257–1265.

Long-read phasing is a sensible direction. However, PCR amplicon counts are not automatically proportional to genomic copy number, and an absent product cannot by itself distinguish a deletion from amplification failure or allele dropout. Resolving an individual gene sequence also does not establish which duplicated genes are linked on the same chromosome without suitable junction or flanking evidence. These assay limitations and structural testing approaches are discussed in the [PharmVar structural-variation tutorial](https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/cpt.3044).

Keep dosage confirmation and haplotype sequencing as separate assay objectives unless the combined assay is quantitatively validated. Design explicit deletion/duplication/hybrid coverage, with appropriate controls. Reflex criteria should consider insufficient callability, technical failure, unknown function, and structural inconsistency as well as reported phase ambiguity. The observed 1/53 clinically uncertain ambiguity rate is not an established future reflex rate.

## Assessment of the methods and alternatives

### Depth normalization and the GMM: keep as a baseline, validate transferability

The depth-plus-paralog-ratio design is well motivated for homologous genes and is the basis of the [original Cyrius method](https://pmc.ncbi.nlm.nih.gov/articles/PMC7997805/). Rejecting unusable controls and retaining low-MAPQ dosage reads are sensible. A generic small-variant VCF workflow would lose information the caller needs.

The fixed GMM parameters and narrow control set should nevertheless be tested across library preparation, aligner, read length, duplicate burden, and coverage. A familiar integer-like CN is not sufficient evidence that a new pipeline has the same noise or bias distribution. Preserve plausible CN alternatives near boundaries rather than committing immediately to one integer and one downstream structure.

Consider cohort- or pipeline-calibrated depth parameters with robust outlier handling. Do not refit against the truth labels later used to report accuracy. Report effective control coverage and target callability alongside MAD.

### Site dosage: compare a beta-binomial model with the current Poisson rule

[copy_number_call.rs](../../cyrius-rs/src/depth_calling/copy_number_call.rs) scores the smaller allele count with a Poisson PMF, using the observed count sum to set expected depth. Given total informative depth, a binomial distribution is a more direct starting model for an allele split; a beta-binomial can model excess variation from alignment and sampling bias. This is a proposed improvement, not an assertion that a new model will automatically improve benchmark accuracy.

For site i, score D6-supporting fragments conditional on total informative fragments and the D6 dosage fraction, with learned bias and dispersion where the validation data support them. Distinguish D6-specific target sites from pooled D6/D7 diagnostic sites. Use log probabilities and retain the full dosage likelihood vector. Account for correlations when the same fragments cover adjacent sites.

### Structural inference: prefer a count-based HMM experiment over more fallback rules

The current consensus method is interpretable and should remain the comparator. A beta-binomial HMM over diagnostic counts is a reasonable next experiment because it can use depth information directly, retain uncertainty, and constrain transitions to biologically plausible configurations.

The current optional HMM uses Gaussian emissions on rounded CN ratios, a common site-depth threshold, and a constant transition probability per marker. Its emission variance does not otherwise distinguish a site with ten observations from one with hundreds, and marker spacing is uneven. Consider count-dependent emissions and distance/structure-aware transitions, with held-out tuning. Summarize posterior support as well as a Viterbi path.

The changepoint module is recursive binary segmentation using a maximum Welch statistic; it is described as CBS but does not implement a full circular/permutation-calibrated CBS procedure. Its fixed statistic thresholds should not be read as controlled error probabilities. Compare consensus, the count-based HMM, and segmentation on the same locked challenge set. A fallback that changes only previously uncallable cases still requires validation of those rescued calls.

### Candidate generation: stop discarding uncertainty before scoring it

Exact multiset lookup is useful for well-observed, supported configurations. Its limitations should be explicit: many high-CN tables assume at most two distinct base alleles, some arrangements are assigned canonical notation, and rare candidates are removed by a fixed list in `get_star`.

Use the table as a constrained candidate generator, retaining uncertain dosage alternatives and data-consistent rare alleles. If population frequencies are used, expose them as optional, documented priors; retain likelihood-based ambiguity. A rarity rule is not direct read evidence.

Fuzzy matching currently permits up to two missing and two extra variant copies, scoring every matched variant equally. It does not distinguish a low-coverage missing variant from a well-covered contradictory reference genotype, or a benign marker from a function-changing splice variant. The README already documents incorrect fuzzy calls on some other pipelines. Keep fuzzy suggestions distinct from supported genotype conclusions, or replace their scoring with explicit read evidence and contradiction rules. Increasing heuristic score weights would not repair the lost information.

If the candidate tables become a scaling or expressiveness limit, evaluate constrained integer optimization, as used by [Aldy 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC9977157/), against the same independent data. That is a useful comparator and possible architecture, not a reason to replace this project before fixing the demonstrated defects.

### Fragment evidence: consolidate around one quality-aware representation

Use one fragment-observation implementation for special-site phasing, generic phasing, and the likelihood scorer. It should preserve base quality, mapping uncertainty, allele-pattern support, overlap conflicts, and technical errors. A single molecule should contribute once to an evidence term. Keep unknown phase explicit, and model hybrid sequence segments where needed.

Longer inserts can extend the range of short-read phasing. More elaborate scoring cannot establish physical phase between disconnected variant blocks without adding assumptions or new data. Population phasing, RNA evidence, or a graph aligner should not silently convert these cases into certain individual haplotypes.

### K-mers: useful corroboration, limited independence

The optional k-mer module avoids base-to-reference alignment when matching the extracted read sequences, which is useful. It still uses the same alignment-selected BAM reads and the depth caller’s total CN. Its estimated D6 dosage is explicitly `total_cn * overall_ratio`, and agreement is checked against an expectation using that same total CN. It cannot independently verify the absolute CN or recover reads discarded upstream.

Keep it as a paralog-ratio consistency check. For a more independent measurement, count suitable sequence markers from raw reads or a broadly recruited read set, normalize against independent control markers, and validate marker uniqueness and sequence-context bias. Report unavailable QC separately from measured disagreement, including partial-fetch failures.

### Long reads and graph methods: use when they add missing information

Long reads spanning the relevant variant blocks or structural junctions are the appropriate way to resolve some remaining ambiguities. Evaluate a validated assay and an appropriate platform-specific caller; [StarPhase](https://pmc.ncbi.nlm.nih.gov/articles/PMC11661245/) is one published HiFi comparator. Long reads still require careful locus reconstruction and copy-number interpretation.

Graph or local multi-reference realignment can be useful for paralog and indel alignment bias. The existing results on graph-aligned data show that it is not enough to substitute an aligner and assume the caller’s previous calibration transfers. Benchmark candidate reference sequences, alignment handling, and inference together.

## Validation and implementation priorities

1. Add the reproduced missing-site, exact-copy activity, fragment-order, base-quality, unknown-phase, and corrupt-record cases as permanent regression tests; fix those behaviors first.
2. Define explicit evidence and result types: unknown site dosage, uncertain phase, technical failure, supported genotype, and unresolved candidate set. Replace string-formatted counts and positional tuples at inference boundaries with structured fields.
3. Correct confidence semantics and reference-callability checks. Keep genotype confidence, structural certainty, and clinical-effect equivalence distinct.
4. Make the benchmark a reproducible gate, and add CI that can run synthetic tests without local BAM symlinks. Retain a separately provisioned real-data suite. The existing tests pass but missed all six reproduced correctness issues.
5. Freeze rules and resources before evaluating new donors. Split by donor and family; multiple alignments or platforms for the same individual are portability tests, not additional independent genotype truths. Use orthogonal truth for difficult structural cases, with adjudication recorded independently of this caller’s answer.
6. Challenge callability and robustness deliberately: downsample by fragment, vary insert length and read length, remove coverage at defining sites, introduce conflicting mates and low-quality observations, test high CN and hybrid configurations, and retain unsupported or novel patterns as such. Assess no-calls and false confident phenotype assignments as well as exact diplotype concordance.
7. Then compare count-based dosage/HMM models and a corrected fragment likelihood model in shadow mode. Predefine the error and call-rate criteria before choosing a method. Do not add more independent-looking scores derived from the same evidence without accounting for dependence.

Version the embedded PharmVar/CPIC snapshots and resource-generation process explicitly, including reference assembly identity and checksums. The tool build hash is useful provenance, but it does not by itself document database source releases or validate that an input uses the expected GRCh38 coordinates. The functionality table lacks rows for 24 main alleles present in the star table; missing function is appropriately indeterminate, but the supported interpretation scope should be visible.

Historical documentation also needs reconciliation. `PDX_CALLER_METHODS.md` says GC correction is absent even though the live caller performs LOWESS correction. `ANALYSIS_3_WRONG_SAMPLES.md` contains sample truth labels that disagree with the versioned current truth table and speculative explanations of partial duplications. Treat it as an archived hypothesis record, not as current validation evidence. Full-copy star notation should not be justified as a representation of a truncated copy solely from a depth gradient.

The recommended direction is to retain the working Cyrius baseline, make evidence loss and uncertainty explicit, establish independent validation, and then evaluate richer inference. The demonstrated failures require changes to evidence handling even if every statistical method remains otherwise unchanged.

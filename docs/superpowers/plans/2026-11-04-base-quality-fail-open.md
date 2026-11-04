# Base-Quality Fail-Open Paths Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Settle by measurement whether the base-quality gate added in
`fix/review-p1` can silently turn withheld evidence into a reference call,
and if it can, make that condition reported instead of silent — then give the
resolver the same view of an unrecorded base quality.

**Architecture:** `record.qual()` is read in exactly three places in the crate
and each place has a different rule. Two of them (`haplotype.rs:251`,
`phase_disambiguate.rs:512`) now reject a base below 13 or equal to `0xFF`.
The third (`diplotype_caller.rs:53`) clamps the value to 93, so `0xFF` becomes
the *highest* confidence. This plan does not unify them by refactor; it makes
the fail-open consequence visible and fixes the clamp.

**Tech Stack:** Rust, rust-htslib 0.47.1, `cargo test`, the 53-sample GeT-RM
benchmark, and the gate script at `~/cyrius-research/review-fixes/gate.sh`.

## Global Constraints

- No threshold is chosen after seeing a distribution. Any threshold is written
  into this plan, and committed, before the run that measures it.
- Every expected output change goes into
  `~/cyrius-research/review-fixes/expected_task_<N>.tsv` and
  `expected_cumulative_after_<N>.tsv` **before** the gate is run. Empty file
  means "no change permitted".
- `cargo clippy --all-targets --locked` must report at most **65** messages
  (the measured baseline in `lint_baseline.txt`).
- `cargo test` must emit at least one `test result` line and every such line
  must show `0 failed`.
- Commit order per claim: `plan:` → `code:` → `result:`. Commit 1 and commit 3
  are never the same commit.
- Beyond the gate's 10 fields, compare all 33 output fields across all 53
  samples directly, and assert field presence before trusting equality. A
  comparison that finds every field absent is vacuous, not passing.

---

## Gate A — answered before any code

**1. Governing principle.** Evidence that the caller decided not to trust must
be reported, not silently converted into a confident reference call. Abstaining
is allowed; pretending is not.

**2. What measurement would kill this.** Construct an input on which a variant
*is* called with good base qualities, then set those qualities to `0xFF` and
nothing else. If the variant is still called, or if no such input can be
constructed at all, there is nothing here to fix and the plan stops at Task 1.

**3. Has this already been refuted here?** The specific claim has not been
proposed before, but the review's stated evidence for it **has been measured
and does not hold**:

> Review I1: "Before this branch the same BAM called them."

MEASURED on HG01190 — the only benchmark sample with non-zero alt counts at
both `g.42126938C>T` (16,64) and `g.42127526C>T` (21,29) — with every `QUAL`
field set to `*` (1 279 889 reads in, 1 279 889 out, verified):

| binary | input | Genotype | Filter | Confidence |
|---|---|---|---|---|
| 3807326 (before the branch) | original | `*5/*68+*4` | Fuzzy_match | 0.61 |
| 3807326 | qualities stripped | `*5/*68+*4` | Fuzzy_match | 0.61 |
| e035eb4 (branch HEAD) | original | `*5/*68+*4` | Fuzzy_match | 0.61 |
| e035eb4 | qualities stripped | `*5/*68+*4` | Fuzzy_match | 0.61 |

No variant is lost by stripping, under either binary. 0 of 31 comparable
fields differ. Control: the HEAD original run reproduces the accepted
`Genotype` exactly. Reproduction kept in
`~/cyrius-research/review-fixes/qual_stripped/`.

The reason, MEASURED: **no sample of the 53 calls any of the three variants.**
The grep hits for those positions in the output JSON are in
`Variant_raw_count` (the raw depth table), not `Variants_called`. Both
producers are gated on a haplotype pattern that is already closed on this
cohort:

- `call_var42126938` (call_variants.rs:433) needs
  `full_length_cn >= 2 && cn < full_length_cn - 2` **and**
  `recombinant_read_count["12"].sum() > 1`.
- `call_var42127526_var42127556` (call_variants.rs:475) computes
  `var7526_cn` and `var7556_cn` from **depth**, then discards both unless
  `has_211` or `has_221` — each derived from the quality-gated haplotype path.

So the fail-open **structure** is verified by reading the code, and the
branch's own test proves every read is dropped at those sites when qualities
are absent. What is **not** shown is that this changes any call on any input
yet tried. The benchmark gives zero coverage of it. Task 1 exists solely to
settle that, and the plan is expected to die there if it cannot.

**4. Simplest thing that could work.** One check when the BAM is opened: if the
base qualities are absent, say so in the output. That covers the whole class in
one place, needs no per-variant plumbing and no new threshold. Per-site
reporting (Task 3) is justified only if Task 1 produces a call change that this
check would not flag.

**5. What must be true of the inputs.**

| Assumption | How it was verified (not why it is plausible) |
|---|---|
| `record.qual().len() == record.seq().as_bytes().len()` | Read in rust-htslib 0.47.1: `Record::seq()` sets `len: self.seq_len()` (record.rs:572-577), `Seq::as_bytes()` is `(0..self.len())` (1636-1638), `Record::qual()` returns `[..self.seq_len()]` (582-585). |
| htslib materialises SAM `QUAL = *` as `0xFF` | Built such a BAM and read it back; `samtools view` shows `*`, and the branch's `QUAL_UNAVAILABLE` path fires on it. |
| The benchmark BAMs carry real base qualities | Control run on the unmodified HG01190 reproduces the accepted call; `samtools view` on it shows a real quality string. |
| `0xFF` reaching the resolver is scored as high quality | Read `QualLookup::score` (diplotype_caller.rs:52-59): `let q = (qual as usize).min(93)`, and index 93 is the largest match / smallest mismatch likelihood. |

---

## Gate B — the kill experiment comes first

Task 1 is the cheapest experiment that can end this plan, and it runs before
anything is changed in `src/`. If `has_211` cannot be made true by any
fixture, I1 is unreachable: record that, delete Tasks 2 and 3, and keep only
Task 4 (which stands on its own).

**Decision point:** stop after Task 1 and report the result to the user before
starting Task 2.

---

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `cyrius-rs/src/caller/call_variants.rs` | the two variant producers and their tests | Task 1 adds tests only; Task 3 may add reporting |
| `cyrius-rs/src/depth_calling/haplotype.rs` | `get_bases_per_read` and the fixture helpers | Task 1 extends `haplotype_fixture_with_quality` to 3 sites |
| `cyrius-rs/src/main.rs` | pipeline, output assembly | Task 2 adds the input check and its output field |
| `cyrius-rs/src/caller/diplotype_caller.rs` | resolver scoring | Task 4 fixes the clamp |
| `README.md` | the user-facing claim | Tasks 2 and 4 update it |

---

## Task 1: the reproduction (kill gate — no `src/` behaviour changes)

**Files:**
- Modify: `cyrius-rs/src/depth_calling/haplotype.rs` (test module only)
- Modify: `cyrius-rs/src/caller/call_variants.rs` (test module only)

**Interfaces:**
- Consumes: `write_hap_bam_with_quality(path, read_len, &[(name, pos, &[(offset, base)])], q)`,
  `hap_lookup(&[(d6_pos_1based, d7_pos_1based, "REF_ALT")])`,
  `Scratch::new`, `extract_hap` — all already in `haplotype.rs`'s test module.
- Produces: `hap3_fixture_with_quality(q) -> (Scratch, bam::IndexedReader, SnpLookup)`,
  a three-site version of the existing two-site helper.

- [ ] **Step 1: write the failing test — three sites, pattern `211`, good quality**

Add to `haplotype.rs`'s test module. Three 1 bp sites mirroring the real
`call_var42127526_var42127556` panel shape; four reads, two of which spell
`211` (site 0 → D7 base, sites 1 and 2 → D6 base).

```rust
    /// Three-site version of `haplotype_fixture_with_quality`, shaped like the
    /// panel `call_var42127526_var42127556` reads. Two reads spell `211`,
    /// which is the pattern that caller gates its push on.
    fn hap3_fixture_with_quality(q: u8) -> (Scratch, bam::IndexedReader, SnpLookup) {
        let dir = Scratch::new(&format!("hap3_{q}"));
        let bam_path = dir.join("hap3.bam");
        write_hap_bam_with_quality(
            &bam_path,
            150,
            &[
                ("a211", 1000, &[(0, "T"), (30, "G"), (60, "C")]),
                ("b211", 1000, &[(0, "T"), (30, "G"), (60, "C")]),
                ("c111", 1000, &[(0, "G"), (30, "G"), (60, "C")]),
                ("d111", 1000, &[(0, "G"), (30, "G"), (60, "C")]),
            ],
            q,
        );
        let base_db = hap_lookup(&[
            (1001, 2001, "G_T"),
            (1031, 2031, "G_A"),
            (1061, 2061, "C_T"),
        ]);
        let reader = bam::IndexedReader::from_path(&bam_path).unwrap();
        (dir, reader, base_db)
    }

    #[test]
    fn three_site_fixture_spells_211_at_good_quality() {
        // Control for the test below: the fixture must produce the pattern the
        // caller gates on, or the 0xFF test proves nothing.
        let (_dir, mut reader, base_db) = hap3_fixture_with_quality(30);
        let haps = get_haplotypes_from_bam_single_region(&mut reader, &base_db, &[0, 1, 2]).unwrap();
        let counts = extract_hap(&haps, &[0, 1, 2]);
        assert!(!counts.is_empty(), "fixture produced no haplotypes at all");
        assert!(
            counts.get("211").map_or(0, |v| v.iter().sum::<u32>()) > 1,
            "need >1 read spelling 211 for has_211 to fire, got {counts:?}"
        );
    }

    #[test]
    fn three_site_fixture_spells_nothing_without_qualities() {
        let (_dir, mut reader, base_db) = hap3_fixture_with_quality(QUAL_UNAVAILABLE);
        let haps = get_haplotypes_from_bam_single_region(&mut reader, &base_db, &[0, 1, 2]).unwrap();
        let counts = extract_hap(&haps, &[0, 1, 2]);
        assert_eq!(
            counts.get("211").map_or(0, |v| v.iter().sum::<u32>()),
            0,
            "unrecorded qualities must drop the reads, got {counts:?}"
        );
    }
```

- [ ] **Step 2: run them**

Run: `cd cyrius-rs && cargo test --release --locked three_site_fixture -- --nocapture`

Expected: **PREDICTED (not run)** — the base letters and offsets above are a
first attempt at the `211` encoding and the `hap_lookup` REF/ALT order. If
`three_site_fixture_spells_211_at_good_quality` fails, read the printed
`counts` map and correct the fixture's bases to match the encoding
`extract_hap` actually uses. **Do not** weaken the assertion to match the
output; fix the fixture. The second test must pass only once the first does.

- [ ] **Step 3: lift it to the caller — the decisive test**

Add to `call_variants.rs`'s test module. This is the measurement that decides
the plan, because it calls the real producer.

```rust
    #[test]
    fn unrecorded_base_quality_silently_drops_a_called_variant() {
        // The fail-open claim, end to end. The depth counts come from
        // `snp_count`, which never looks at base quality, so they are
        // unchanged; the pushed variant comes from the quality-gated
        // haplotype path. If the variant disappears while the depth counts
        // stay put, evidence was withheld and nothing reported it.
        let cnvtag = /* a tag present in cnvtag_lookup_table() whose
                        exon9_to_intron1 is >= 1 */;
        let (_d1, mut r_good, db) = hap3_fixture_with_quality(30);
        let (good_7526, good_7556, good_vars) =
            call_var42127526_var42127556(&mut r_good, cnvtag, &db).unwrap();
        assert!(
            good_vars.contains(&"g.42127526C>T".to_string()),
            "control failed: the good-quality fixture must call the variant, got {good_vars:?}"
        );

        let (_d2, mut r_bad, db2) = hap3_fixture_with_quality(QUAL_UNAVAILABLE);
        let (bad_7526, bad_7556, bad_vars) =
            call_var42127526_var42127556(&mut r_bad, cnvtag, &db2).unwrap();
        assert_eq!(bad_7526, good_7526, "depth counts must not depend on quality");
        assert_eq!(bad_7556, good_7556, "depth counts must not depend on quality");
        assert!(
            bad_vars.is_empty(),
            "expected the variant to vanish, got {bad_vars:?}"
        );
    }
```

The `hap3_fixture_with_quality` helper and `QUAL_UNAVAILABLE` must be reachable
from this module; re-export from the test module or move the helper to a shared
`#[cfg(test)]` location rather than copying it.

- [ ] **Step 4: run it and record the verdict**

Run: `cd cyrius-rs && cargo test --release --locked unrecorded_base_quality -- --nocapture`

Three outcomes, and each one is a result worth committing:

| Outcome | Meaning | Next |
|---|---|---|
| Both the control and the assertion pass | I1 is real and reachable | continue to Task 2 |
| The control cannot be made to pass | I1 is unreachable through this caller | stop; keep only Task 4 |
| The variant is still called without qualities | the gate does not reach the push | stop; I1 is refuted |

- [ ] **Step 5: commit the result, separately from the plan**

```bash
git add -A
git commit -m "result: base-quality fail-open — <supported|refuted|unreachable> at the caller

<the actual test output, pasted>"
```

- [ ] **Step 6: STOP and report to the user before Task 2.**

---

## Task 2: report a BAM whose base qualities are absent

**Only if Task 1 came out "supported".**

**Files:**
- Modify: `cyrius-rs/src/main.rs`
- Modify: `README.md`
- Test: `cyrius-rs/src/main.rs` test module

**Interfaces:**
- Produces: a `Base_qualities_absent: bool` output field, and a `Filter` value.

The threshold, fixed here before measurement: inspect the **first 1000 primary,
non-secondary, non-supplementary alignments** in the CYP2D6 region. If **every
one** of them has `record.qual()` entirely `0xFF`, the input is reported as
having no base qualities. "Every one" rather than a fraction, because a
partially-stripped BAM is a different condition and is Task 3's business, not
this check's.

- [ ] **Step 1: write the failing test**

```rust
    #[test]
    fn a_bam_without_base_qualities_is_reported() {
        let bam = write_bam_with_all_quals(QUAL_UNAVAILABLE, 1000);
        assert!(base_qualities_absent(&bam).unwrap());
    }

    #[test]
    fn a_bam_with_base_qualities_is_not_reported() {
        // Control: the check must fail for the qualities, not for the fixture.
        let bam = write_bam_with_all_quals(30, 1000);
        assert!(!base_qualities_absent(&bam).unwrap());
    }

    #[test]
    fn one_good_read_is_enough_to_not_report() {
        // Pins "every one", not "most". A single recorded quality means the
        // file is not quality-stripped, whatever the rest looks like.
        let bam = write_bam_mixed_quals(QUAL_UNAVAILABLE, 999, 30, 1);
        assert!(!base_qualities_absent(&bam).unwrap());
    }
```

- [ ] **Step 2: run them; expect failure with "cannot find function `base_qualities_absent`"**

- [ ] **Step 3: implement `base_qualities_absent`** in `main.rs`, next to the
existing header check, and wire its result into the output record and `Filter`
exactly as `Header_drops_multimappers` is wired. Follow that field: it is the
established pattern for "an input property that changes `Filter`".

- [ ] **Step 4: run the three tests; all pass**

- [ ] **Step 5: write the expected output change BEFORE running the gate**

The new field appears on every sample, so it IS an output change. Write it to
`~/cyrius-research/review-fixes/expected_task_2.tsv` as
`<sample>\tBase_qualities_absent\tfalse` for all 53, then run
`bash ~/cyrius-research/review-fixes/gate.sh 2`. A gate run before the
expectations are written does not count.

- [ ] **Step 6: update the README** — the bullet added in e035eb4 says this
release does not detect a quality-stripped input. It now does; say so, and say
what the check does not cover (partial or binned low quality).

- [ ] **Step 7: commit**

---

## Task 3: report partially gated-out evidence

**Only if Task 1 produced a call change that Task 2's whole-file check would
not flag** — that is, an input with real qualities where some are below 13.
If Task 2's check covers every case Task 1 found, skip this task and say so;
do not build it speculatively.

Scope to decide at that point, with the user: `Uncallable_target_sites` cannot
carry these sites, because they are not in `CYP2D6_target_variant_38.txt` and
the field's name would then be false. A separate field is the honest option.

---

## Task 4: the resolver must not score an unrecorded quality as the best one

**Stands on its own; independent of Task 1's verdict.**

**Files:**
- Modify: `cyrius-rs/src/caller/diplotype_caller.rs`
- Test: same file's test module

**Interfaces:**
- Consumes: `QUAL_UNAVAILABLE` and `PHASE_MIN_BASE_QUALITY` from
  `crate::depth_calling` (already re-exported there).
- Produces: no signature change. `QualLookup::score` keeps
  `fn score(&self, is_match: bool, qual: u8) -> f64`.

- [ ] **Step 1: write the failing test**

```rust
    #[test]
    fn an_unrecorded_quality_is_not_the_most_confident_one() {
        let lut = QualLookup::new();
        // 0xFF means "no quality recorded". Clamping it to 93 made it the most
        // confident base in the table, which is backwards.
        assert!(
            lut.score(true, QUAL_UNAVAILABLE) < lut.score(true, 93),
            "an unrecorded quality must not outscore a measured q93 match"
        );
        assert!(
            lut.score(true, QUAL_UNAVAILABLE) <= lut.score(true, PHASE_MIN_BASE_QUALITY),
            "an unrecorded quality must not outscore the phasing floor"
        );
    }

    #[test]
    fn a_measured_quality_still_scores_by_its_value() {
        // Control: the change must only touch 0xFF, not the ordinary range.
        let lut = QualLookup::new();
        assert!(lut.score(true, 30) > lut.score(true, 10));
        assert!(lut.score(false, 30) < lut.score(false, 10));
    }
```

- [ ] **Step 2: run them; the first must fail, the second must pass**

If the second fails, the table itself is wrong and that is a separate finding —
stop and report it rather than folding it into this task.

- [ ] **Step 3: implement**

Map `QUAL_UNAVAILABLE` to the lowest-information slot rather than clamping it
upward. Do not silently extend the array; `0xFF` is a sentinel, not a Phred
value, and the code should say so at the mapping site.

- [ ] **Step 4: run the whole suite and the gate**

Write `expected_task_4.tsv` first. **PREDICTED (not run):** empty, because no
benchmark BAM has unrecorded qualities (verified: the control run on the
unmodified HG01190 reproduces the accepted call). If the gate shows any
difference, that prediction was wrong and the difference is the finding —
investigate it, do not widen the expectations file to absorb it.

- [ ] **Step 5: measure whether it changes a stripped-BAM call**

Re-run the kept reproduction and compare against the recorded table in Gate A:

```bash
~/cyrius-research/review-fixes/qual_stripped/HG01190_noqual.bam
```

- [ ] **Step 6: commit**

---

## Self-review notes

- Task 1's fixture bases are labelled `PREDICTED (not run)` because the
  `extract_hap` pattern encoding and the `hap_lookup` REF/ALT order were read
  but not exercised with three sites. Step 2 corrects the fixture, not the
  assertion.
- The `cnvtag` in Task 1 Step 3 is left as a comment on purpose: the right tag
  must be read out of `cnvtag_lookup_table()` at implementation time, and
  guessing it here would be a placeholder dressed as code.
- Every control test above pairs a rejecting input with an accepting one, so
  no zero-count assertion can pass vacuously.

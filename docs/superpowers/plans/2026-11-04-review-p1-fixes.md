# Review P1 Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the six P1 defects from `docs/reviews/2026-11-04-project-review.md`, so that missing, conflicting, low-quality or unreadable evidence can no longer turn into a confident genotype.

**Architecture:** There are six independent fixes. Each one keeps the current Cyrius-derived pipeline and changes only how evidence is handled at one boundary: phenotype lookup, BAM read loops, target-site callability, fragment merging, phasing base quality, and the two targeted cn2 phase rules. Every task is test-first. Every task ends with a per-sample benchmark comparison against a frozen baseline.

**Tech Stack:** Rust (cyrius-rs crate, rust-htslib), bash benchmark harness, Python 3 for the comparison script.

## Global Constraints

- **G1 Branch:** work in a worktree at `/home/parlar_ai/dev/cyrius-rs-p1` on a new branch `fix/review-p1`, created from `feat/diplotype-resolver` at `3807326`. Never merge or push without the user.
- **G2 Commits:** commit steps run only after the user has approved executing this plan. End every commit message with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **G3 No tuning on the benchmark:** no threshold in this plan may be changed after seeing benchmark results. If a gate fails, stop and report; do not adjust the threshold.
- **G4 Regression gate.** Run `~/cyrius-research/review-fixes/gate.sh <task-number>` after every task, after its commit (the binary records the commit it was built from). Task 0 writes it.

  The gate is a script with `set -euo pipefail`, not a chain of piped commands, because piped commands hide failures. Measured 2026-11-04:
  - `( exit 101 ) 2>&1 | tail -1` exits **0** — a failing `cargo build` piped to `tail` looks like a pass.
  - `grep -E "^test result"` over a compiler error prints nothing, so "every `test result` line shows 0 failed" is **trivially true** when nothing compiled.

  ```bash
  #!/usr/bin/env bash
  # ~/cyrius-research/review-fixes/gate.sh <task-number>
  set -euo pipefail
  TASK="${1:?task number}"
  R="$HOME/cyrius-research/review-fixes"
  L="$R/logs/task_$TASK"; mkdir -p "$L"
  fail() { echo "GATE FAIL (task $TASK): $*" >&2; exit 1; }

  cd /home/parlar_ai/dev/cyrius-rs-p1/cyrius-rs
  # Exit status is checked directly; full logs are kept for inspection.
  cargo build --release --locked > "$L/build.log" 2>&1 || fail "build failed, see $L/build.log"
  cargo test  --release --locked > "$L/test.log"  2>&1 || fail "tests failed, see $L/test.log"
  # The summary lines must EXIST (guards the vacuous pass) and all show 0 failed.
  n_result=$(grep -cE "^test result" "$L/test.log" || true)
  [ "$n_result" -ge 1 ] || fail "no 'test result' line in $L/test.log"
  if grep -E "^test result" "$L/test.log" | grep -qv "0 failed"; then
    fail "a test failed, see $L/test.log"
  fi

  # G5 lint. Clippy must have COMPILED before its message count means anything.
  cargo clippy --all-targets --locked > "$L/clippy.log" 2>&1 || true
  if grep -qE "^error(\[|:)" "$L/clippy.log"; then
    fail "clippy did not compile, see $L/clippy.log"
  fi
  n_lint=$(grep -cE "^(warning|error)" "$L/clippy.log" || true)
  [ "$n_lint" -le "$(cat "$R/lint_baseline.txt")" ] || fail "clippy $n_lint > baseline"

  cd /home/parlar_ai/dev/cyrius-rs-p1
  BENCH_KEEP=1 BENCH_BAM_DIR=/home/parlar_ai/dev/cyrius-rs/benchmark \
    bash bench/run_benchmark.sh --strand-bias-all --fuzzy-match --changepoint-hybrid --het-check \
    > "$L/bench.log" 2>&1 || fail "benchmark run failed, see $L/bench.log"
  # run_benchmark.sh prints "Error:      0" (padded), so match on whitespace, anchored.
  grep -qE "^Error:[[:space:]]+0$" "$L/bench.log" || fail "benchmark reported errors"
  if grep -q "No BAM" "$L/bench.log"; then fail "a BAM was missing"; fi
  NEW=$(sed -n 's/^Results kept in //p' "$L/bench.log")
  [ -n "$NEW" ] && [ -d "$NEW" ] || fail "no results directory in $L/bench.log"

  # Two comparisons, per G8.
  # build.rs stamps `git rev-parse --short=7 HEAD` (+ "-dirty" if the CRATE dir has changes); match it exactly.
  HEAD_SHORT="$(git rev-parse --short=7 HEAD)"
  python3 "$R/compare_runs.py" "$R/accepted" "$NEW" \
    --expect-commit "$HEAD_SHORT" --expected "$R/expected_task_$TASK.tsv" \
    || fail "step comparison failed"
  python3 "$R/compare_runs.py" "$R/baseline" "$NEW" \
    --expect-commit "$HEAD_SHORT" --expected "$R/expected_cumulative_after_$TASK.tsv" \
    || fail "cumulative comparison failed"

  echo "$NEW" > "$L/results_dir.txt"
  echo "GATE PASS (task $TASK): results in $NEW"
  ```

  The gate passes only when the script prints `GATE PASS` and exits 0. **On a pass, promote the run** so the next task compares against it:
  ```bash
  R=~/cyrius-research/review-fixes
  rm -rf "$R/accepted" && cp -r "$(cat "$R/logs/task_N/results_dir.txt")" "$R/accepted"
  ```
- **G5 Lint:** folded into `gate.sh` above; the baseline lives in `~/cyrius-research/review-fixes/lint_baseline.txt`.

  **Measured 2026-11-04, and it invalidated the old figure.** At `3807326`, `cargo clippy --all-targets` *could not compile* the `output_provenance` test: `clippy::option_env_unwrap` is deny-by-default and fires on `option_env!("CYRIUS_GIT_COMMIT").expect(..)` at `cyrius-rs/tests/output_provenance.rs:19`. The old G5 command piped clippy into `grep -c` and ignored its exit status, so it counted the output of an **aborted** run and reported `68`. A clippy run that actually compiles reports **65**. The plan's recorded baseline was an artifact of the failure it failed to notice — the reviewer's point about piped exit statuses, demonstrated on committed code.

  Resolved by prerequisite commit `401e324` (user chose this on 2026-11-04 over rewriting the helper to use `env!`): a narrow `#[allow(clippy::option_env_unwrap)]` on that one function, with the reasoning recorded in a doc comment. The lint's premise — don't turn a compile-time-knowable absence into a run-time panic — is exactly the trade the helper wants, because the test must be able to *fail* rather than refuse to build, and `build.rs` always emits the variable anyway.

  **Baseline is therefore 65**, and `gate.sh` checks that clippy compiled before trusting any count. A count below the baseline never passes when clippy failed to compile.
- **G8 Expected changes are cumulative and valued.** Two failures the first draft of this plan had:
  - *Chaining.* Task 3 accepts changes to `Filter`, `Confidence` and `Confidence_label`. Comparing Task 4 against the frozen `3807326` baseline with "no changes allowed" would then fail on Task 3's accepted output. So every task compares twice: against `accepted/` (the previous passing run) with **only that task's** expected changes, and against `baseline/` (frozen at `3807326`) with the **cumulative** expected changes.
  - *Values.* An expected-changes file naming only `sample<TAB>field` lets any new value through; the reviewer's probe got `"Confidence": "not-a-number"` past it. Every line must name the expected new value: `sample<TAB>field<TAB>value`. A file must always exist, even when empty: an empty file means no change is allowed. A listed change that did not happen is also a failure.
- **G9 Honest status in Task 7.** Review #1 is recorded as *diagnosed and reported; fix deferred to its own plan* — this plan adds the `Uncallable_target_sites` field and the measurement, and changes no call. Review #5's negative-phase branch is recorded as *removed, not replaced*. Neither gets a tick.
- **G6 Divergence from bcyrius:** fixes 1, 4 and 5 change behaviour inherited from Cyrius/bcyrius. This is deliberate: the user decided on 2026-06-03 that the Rust code is reviewed on its own merits. Task 7 lists each change in the README's "Changes from bcyrius" section.
- **G7 Decision D1 — resolved, and both of its options dropped.** D1 asked whether an uncallable site should force a no-call or a warned-but-provisional genotype. Neither: the answer is to keep the uncertainty inside allele matching and report a candidate set (see Task 3's SCOPE box). That change lives in its own plan; this plan keeps only the reporting field. No D1 confirmation is needed any more.

## Measured facts this plan relies on

- At `3807326`, 179 tests pass and 1 is ignored (`cargo test --release --locked`).
- The NA19920 reproduction was run on 2026-11-04 with `3807326`:

  | Input | Genotype | Filter | Confidence | Phenotype | `g.42128945C>T` reads |
  |---|---|---|---|---|---|
  | Original | `*1/*4x2` | PASS | 0.87 HIGH | Intermediate | 35,27 |
  | 68 reads at the site removed | `*1/*10x2` | PASS | 0.95 HIGH | Normal | 0,0 |

- **The damage is not confined to the splice site** (measured 2026-11-04, `3807326` release binary, on the reviewer's `NA19920_missing_splice.bam`). Removing those 68 alignments takes **three** target sites to zero reads, not one:

  | Site | Original `alt,ref` | Damaged |
  |---|---|---|
  | `g.42128945C>T` | `35,27` | `0,0` |
  | `g.42128936-42128937insGGGGCGAAA` | `0,63` | `0,0` |
  | `g.42128936-42128937insGGGGCGAAAGGGGCGAAA` | `0,63` | `0,0` |

  and **31** target sites lose depth in all, several nearly so: `g.42128933G>A` `0,62 -> 0,3`, `g.42128932C>T` `0,62 -> 0,4`, `g.42128926T>C` `0,62 -> 0,6`, `g.42128903-42128904insTA` `0,58 -> 0,11`. At total CN 5 a handful of reads may not survive the posterior cutoff either.

  **Now measured, by 3a's own code path (commit `ed5be20`): the count is 5.**
  Three sites from zero coverage, and two more that keep a few reads and still
  fail the cutoff:

  | Site | reads left (alt,ref) | why |
  |---|---|---|
  | `g.42128945C>T` | `0,0` | no coverage |
  | `g.42128936-42128937insGGGGCGAAA` | `0,0` | no coverage |
  | `g.42128936-42128937insGGGGCGAAAGGGGCGAAA` | `0,0` | no coverage |
  | `g.42128932C>T` | `0,4` | too few reads at total CN 5 |
  | `g.42128933G>A` | `0,3` | too few reads at total CN 5 |

  This reproduces the external reviewer's prototype figure of five uncallable
  dosages, which this plan had recorded as the reviewer's output rather than a
  fact. Both the number and the mechanism now come from the shipped code.

- **Confidence rose when evidence was removed:** `0.87 HIGH` → `0.95 HIGH` (same measurement). The damaged input scores *higher* than the intact one, because the missing site stops contradicting `*10`. This is the defect in one line.
- `call_reg1_cn` returns `None` when a site has no reads, and `Some(0)` when the expected copy number there is 0 ([copy_number_call.rs:13-18](../../../cyrius-rs/src/depth_calling/copy_number_call.rs)). So a homozygous deletion (NA19317, `*5/*5`) gives 0, not `None`.
- Benchmark read depth at the 140 target sites (a75638a outputs): only NA19317 has sites with zero reads, which is expected for `*5/*5`. Every other sample has at least 7 reads at every site.
- None of the 53 benchmark calls is one of the four diplotypes that the targeted cn2 rules choose between (`*1/*32`, `*27/*41`, `*119/*2`, `*1/*41`). So Task 6 cannot change a benchmark call.
- Haplotype codes in `haplotype.rs`: `'1'` = CYP2D6 base, `'2'` = CYP2D7 base, `'x'` = not seen ([haplotype.rs:286-293](../../../cyrius-rs/src/depth_calling/haplotype.rs)).

---

### Task 0: Worktree, frozen baseline and comparison script

**Files:**
- Create: `~/cyrius-research/review-fixes/compare_runs.py`
- Create: `~/cyrius-research/review-fixes/baseline/` (benchmark outputs at `3807326`)

- [ ] **Step 1: Create the worktree**

```bash
cd /home/parlar_ai/dev/cyrius-rs
git worktree add ../cyrius-rs-p1 -b fix/review-p1 feat/diplotype-resolver
ln -s /home/parlar_ai/dev/cyrius-rs/benchmark /home/parlar_ai/dev/cyrius-rs-p1/benchmark
cd /home/parlar_ai/dev/cyrius-rs-p1 && git log --oneline -1
```
Expected: `3807326 feat: update README and .gitignore for graft's local graph cache`

- [ ] **Step 1b: Copy the untracked documents into the worktree.** `docs/reviews/2026-11-04-project-review.md` and this plan are untracked at `3807326`, so a fresh worktree does not contain them and Task 7 would have nothing to edit.

```bash
mkdir -p /home/parlar_ai/dev/cyrius-rs-p1/docs/reviews /home/parlar_ai/dev/cyrius-rs-p1/docs/superpowers/plans
cp /home/parlar_ai/dev/cyrius-rs/docs/reviews/2026-11-04-project-review.md \
   /home/parlar_ai/dev/cyrius-rs-p1/docs/reviews/
cp /home/parlar_ai/dev/cyrius-rs/docs/superpowers/plans/2026-11-04-review-p1-fixes.md \
   /home/parlar_ai/dev/cyrius-rs-p1/docs/superpowers/plans/
ls /home/parlar_ai/dev/cyrius-rs-p1/docs/reviews/ /home/parlar_ai/dev/cyrius-rs-p1/docs/superpowers/plans/
```
Expected: both files listed. They are committed in Task 7.

- [ ] **Step 2: Write the comparison script**

```python
#!/usr/bin/env python3
"""Compare two benchmark result directories sample by sample.

Usage: compare_runs.py BASE NEW [--expect-commit HASH] --expected FILE

FILE lists the changes this step is allowed to make, one per line:

    sample<TAB>field<TAB>expected_new_value

The value is required. A file naming only the field lets ANY new value
through: in the reviewer's probe a bare field list accepted
`"Confidence": "not-a-number"`. Use the literal `*` as the value only for a
field whose exact value genuinely cannot be known in advance, and say in the
task why.

The file must exist even when empty; an empty file means no change is allowed.

Exit 0 only if: both runs hold the same 53 samples; NEW was built from HASH;
every difference is listed with a matching new value; and every listed change
actually happened (a stale expectation is a failure, not a pass).
"""
import argparse, csv, glob, json, os, sys

JSON_FIELDS = ["Genotype", "Filter", "CNV_group", "Total_CN", "Spacer_CN",
               "Confidence", "Confidence_label", "Resolved_by"]
TSV_FIELDS = ["Activity score", "Predicted phenotype"]


def load(d):
    out = {}
    for f in sorted(glob.glob(os.path.join(d, "*.json"))):
        sample = os.path.basename(f)[:-5]
        record = next(iter(json.load(open(f)).values()))
        row = {k: record.get(k) for k in JSON_FIELDS}
        row["Tool_commit"] = record.get("Tool_commit")
        tsv = f[:-5] + ".tsv"
        rows = list(csv.DictReader(open(tsv), delimiter="\t")) if os.path.exists(tsv) else []
        for k in TSV_FIELDS:
            row[k] = rows[0].get(k) if rows else None
        out[sample] = row
    return out


def read_expected(path):
    if not os.path.exists(path):
        print(f"MISSING EXPECTED FILE: {path} (create it empty if no change is allowed)")
        sys.exit(1)
    want = {}
    for line in open(path):
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 3:
            print(f"BAD EXPECTED LINE (need sample<TAB>field<TAB>value): {line!r}")
            sys.exit(1)
        want[(parts[0], parts[1])] = parts[2]
    return want


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base"); ap.add_argument("new")
    ap.add_argument("--expect-commit")
    ap.add_argument("--expected", required=True)
    a = ap.parse_args()
    base, new = load(a.base), load(a.new)
    expected = read_expected(a.expected)
    ok = True
    if len(base) != 53 or set(base) != set(new):
        print(f"SAMPLE SET: base {len(base)}, new {len(new)}, "
              f"missing {sorted(set(base) - set(new))}, extra {sorted(set(new) - set(base))}")
        ok = False
    if a.expect_commit:
        bad = sorted(s for s, r in new.items() if r["Tool_commit"] != a.expect_commit)
        if bad:
            print(f"STALE BINARY: {len(bad)} samples not built from {a.expect_commit}: {bad[:5]}")
            ok = False
    seen, n_diff = set(), 0
    for s in sorted(set(base) & set(new)):
        for k in JSON_FIELDS + TSV_FIELDS:
            if base[s][k] == new[s][k]:
                continue
            n_diff += 1
            want = expected.get((s, k))
            if want is None:
                tag = "NOT EXPECTED"
            elif want == "*" or want == str(new[s][k]):
                tag = "expected"
                seen.add((s, k))
            else:
                tag = f"WRONG VALUE (wanted {want!r}, got {new[s][k]!r})"
            print(f"{s}\t{k}\t{base[s][k]!r} -> {new[s][k]!r}\t{tag}")
            ok = ok and tag == "expected"
    stale = sorted(set(expected) - seen)
    if stale:
        print(f"EXPECTED BUT UNCHANGED: {stale}")
        ok = False
    print(f"{n_diff} differences; {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Build the frozen baseline and the promotion copy**

```bash
mkdir -p ~/cyrius-research/review-fixes
cd /home/parlar_ai/dev/cyrius-rs-p1/cyrius-rs
cargo build --release --locked > /tmp/base_build.log 2>&1 || { echo "BUILD FAILED"; tail -20 /tmp/base_build.log; exit 1; }
cd /home/parlar_ai/dev/cyrius-rs-p1
BENCH_KEEP=1 BENCH_BAM_DIR=/home/parlar_ai/dev/cyrius-rs/benchmark \
  bash bench/run_benchmark.sh --strand-bias-all --fuzzy-match --changepoint-hybrid --het-check \
  > /tmp/bench_base.log 2>&1 || { echo "BENCH FAILED"; tail -20 /tmp/bench_base.log; exit 1; }
tail -8 /tmp/bench_base.log
R=~/cyrius-research/review-fixes
cp -r "$(sed -n 's/^Results kept in //p' /tmp/bench_base.log)" $R/baseline
cp -r $R/baseline $R/accepted     # the first task compares against this
```
Expected: the benchmark summary shows `Exact: 51 / 53`, `Ambiguous: 2`.

- [ ] **Step 4: Prove the comparator can fail.** A comparator that has only ever seen identical runs is not a verified gate. Each probe below must be rejected for the stated reason.

```bash
R=~/cyrius-research/review-fixes; B=$R/baseline; C="python3 $R/compare_runs.py"
: > $R/empty.tsv
# 1. identical runs, right commit, nothing expected -> PASS
$C $B $B --expect-commit 3807326 --expected $R/empty.tsv; echo "exit $?"
# 2. stale binary -> reject
$C $B $B --expect-commit 0000000 --expected $R/empty.tsv; echo "exit $?"
# 3. a changed call that nothing expected -> reject
T=$(mktemp -d); cp $B/* $T/; sed -i 's/"\*3\/\*3"/"*3\/*4"/' $T/HG00111.json
$C $B $T --expected $R/empty.tsv; echo "exit $?"
# 4. the field is listed but the VALUE is wrong -> reject (the old field-only
#    form accepted any value, including "not-a-number")
printf 'HG00111\tGenotype\t*3/*5\n' > $R/probe.tsv
$C $B $T --expected $R/probe.tsv; echo "exit $?"
# 5. the field is listed with the right value -> PASS
printf 'HG00111\tGenotype\t*3/*4\n' > $R/probe.tsv
$C $B $T --expected $R/probe.tsv; echo "exit $?"
# 6. an expectation that did not happen -> reject (stale expectation)
$C $B $B --expected $R/probe.tsv; echo "exit $?"
# 7. a missing expected file -> reject, never silently "no change allowed"
$C $B $B --expected $R/does_not_exist.tsv; echo "exit $?"
rm -rf $T $R/probe.tsv
```
Expected:
- 1 → `0 differences; PASS`, `exit 0`
- 2 → `STALE BINARY: 53 samples ...`, `exit 1`
- 3 → `HG00111	Genotype	'*3/*3' -> '*3/*4'	NOT EXPECTED`, `exit 1`
- 4 → `... WRONG VALUE (wanted '*3/*5', got '*3/*4')` **and** `EXPECTED BUT UNCHANGED: [('HG00111', 'Genotype')]`, `exit 1` (both lines appear: the listed change did not happen as written)
- 5 → `1 differences; PASS`, `exit 0`
- 6 → `EXPECTED BUT UNCHANGED: [('HG00111', 'Genotype')]`, `exit 1`
- 7 → `MISSING EXPECTED FILE: ...`, `exit 1`

If any probe passes when it should fail, stop: the gate is not measuring anything.

- [ ] **Step 5: Write `gate.sh`** exactly as given in G4, `chmod +x` it, and record the lint baseline **with the gate's own command** so the number and the check cannot drift apart:

```bash
R=~/cyrius-research/review-fixes
cd /home/parlar_ai/dev/cyrius-rs-p1/cyrius-rs
cargo clippy --all-targets --locked > /tmp/clippy_base.log 2>&1 || true
grep -qE "^error(\[|:)" /tmp/clippy_base.log && { echo "clippy did not compile"; exit 1; }
grep -cE "^(warning|error)" /tmp/clippy_base.log | tee $R/lint_baseline.txt
```
Expected: a number is written to `$R/lint_baseline.txt`. **Measured 2026-11-04: `65`** (with `401e324` applied; see G5). The `68` in earlier drafts came from a clippy run that aborted on a compile error, so it was never a valid baseline.

**Prerequisite:** this step cannot pass at `3807326` — clippy fails to compile the `output_provenance` test there. Commit `401e324` fixes that and must be in place first. If `grep -qE "^error(\[|:)"` finds anything in `clippy_base.log`, stop: the count that follows would be meaningless.

- [ ] **Step 6: Prove `gate.sh` can fail.** Break the build on purpose, run the gate, restore.

```bash
cd /home/parlar_ai/dev/cyrius-rs-p1/cyrius-rs
cp src/lib.rs /tmp/lib.rs.bak
echo 'fn deliberate_syntax_error( {' >> src/lib.rs
~/cyrius-research/review-fixes/gate.sh 0; echo "exit $?"
cp /tmp/lib.rs.bak src/lib.rs && rm /tmp/lib.rs.bak
```
Expected: `GATE FAIL (task 0): build failed, see .../build.log` and `exit 1`. A `GATE PASS` here means the gate is dead and the plan stops.

- [ ] **Step 7: Prove `gate.sh` can PASS.** A gate that always fails is as dead as one that always passes. On the clean tree (build restored), with empty expected files:

```bash
R=~/cyrius-research/review-fixes
: > $R/expected_task_0.tsv; : > $R/expected_cumulative_after_0.tsv
cd /home/parlar_ai/dev/cyrius-rs-p1/cyrius-rs && git status --short -- .   # crate dir only, as build.rs checks; must print nothing
$R/gate.sh 0; echo "exit $?"
```
Expected: `GATE PASS (task 0): results in ...` and `exit 0`, with both comparisons printing `0 differences; PASS`. Do not promote this run; `accepted/` stays the baseline copy.

---


### Task 1: Activity for exact copy counts (review #2)

**Files:**
- Modify: `cyrius-rs/src/phenotype.rs:5-31` (loader), `:42-87` (`lookup_functionality`), tests at `:388-396`

**Interfaces:**
- Produces: `lookup_functionality(name, table)`. Same signature. For `*AxN` with no exact row, it now returns N × the single-copy activity. Bucket rows keep their `≥` (key `*Ax≥3`) and only supply the function label for N ≥ 3.

- [ ] **Step 1: Write the failing tests**

In `cyrius-rs/src/phenotype.rs`, replace the test `three_or_more_copies_use_the_x3_row`:

```rust
    #[test]
    fn three_or_more_copies_use_the_x3_row() {
        assert_eq!(activity_of("*1/*2x4", "*2x4"), "3.0");
    }
```

with:

```rust
    #[test]
    fn an_exact_copy_count_multiplies_the_single_copy_activity() {
        // *41 is 0.25 per copy; the *41x3 row is an exact 3-copy row, not a
        // ">= 3" bucket, so 4 and 5 copies must not reuse its 0.75.
        assert_eq!(activity_of("*10/*41x4", "*41x4"), "1.0");
        assert_eq!(activity_of("*5/*41x5", "*41x5"), "1.25");
        // *2x≥3 is a lower bound ("≥3.0"); an exact count of 4 is 4.0.
        assert_eq!(activity_of("*5/*2x4", "*2x4"), "4.0");
    }

    #[test]
    fn an_exact_table_row_still_wins() {
        assert_eq!(activity_of("*1/*41x3", "*41x3"), "0.75");
        assert_eq!(activity_of("*1/*2x2", "*2x2"), "2.0");
    }

    #[test]
    fn exact_multiplication_sets_the_cpic_phenotype() {
        let t = table();
        let p = match_phenotype(Some("*41x4/*10"), &t);
        assert_eq!(p[0].total_activity, "1.25");
        assert_eq!(p[0].predicted_phenotype, "Normal Metabolizer");
        let p = match_phenotype(Some("*41x5/*5"), &t);
        assert_eq!(p[0].total_activity, "1.25");
        assert_eq!(p[0].predicted_phenotype, "Normal Metabolizer");
    }

    #[test]
    fn three_or_more_normal_copies_keep_the_increased_function_label() {
        let (details, _) = compare_genotypes("*1x4/*5", &table());
        let d = details.into_iter().find(|d| d.haplotype == "*1x4").unwrap();
        assert_eq!(d.activity, "4.0");
        assert_eq!(d.function, "Increased function");
    }
```

- [ ] **Step 2: Run them and check that they fail for the stated reason**

Run: `cd /home/parlar_ai/dev/cyrius-rs-p1/cyrius-rs && cargo test --release --locked phenotype:: 2>&1 | grep -E "^test |panicked|left|right"`
Expected:
- FAIL `an_exact_copy_count_multiplies_the_single_copy_activity`, with left `"0.75"` against right `"1.0"`;
- FAIL `exact_multiplication_sets_the_cpic_phenotype`;
- PASS `an_exact_table_row_still_wins`.

- FAIL `three_or_more_normal_copies_keep_the_increased_function_label`, with left `"3.0"` against right `"4.0"` (today `*1x4` borrows the `*1x≥3` row's "3.0").

- [ ] **Step 3: Implement**

In `load_haplotype_functionality`, keep the table text as it is:

```rust
        let haplotype = fields[0].to_string();
        let activity = fields[1].to_string();
```

Replace the doc comment and body of `lookup_functionality` from `let n = copies?;` to the end of the function with:

```rust
    let n = copies?;
    let single = table.get(base)?;
    let activity: f64 = single.activity.parse().ok()?;
    // The copy count is exact, so the activity is n times one copy. A
    // "*Ax≥3" row is a lower bound for an unknown count of 3 or more: it
    // supplies only the function label, never the activity.
    let function = match table.get(&format!("{}x≥3", base)) {
        Some(bucket) if n >= 3 => bucket.function.clone(),
        _ => single.function.clone(),
    };
    Some(HaplotypeFunctionality {
        activity: format_activity(activity * n as f64),
        function,
        evidence_strength: "n/a".to_string(),
        evidence_summary: format!(
            "No table entry for {}; activity is {} x the {} value.",
            name, n, base
        ),
    })
}
```

Update the doc comment above `lookup_functionality` to:

```rust
/// Look up a haplotype's functionality. When there is no exact row:
/// a sub-allele uses its main allele ("*4.013" -> "*4"); otherwise "*AxN"
/// gets N times the single-copy activity of "*A". The table's "*Ax≥3" rows
/// are lower bounds for an unknown count, so they only give the function
/// label for N >= 3.
```

- [ ] **Step 4: Run the tests**

Run: `cargo test --release --locked phenotype:: 2>&1 | grep -E "^test result|FAILED"`
Expected: `test result: ok.` with no `FAILED`.

- [ ] **Step 5: Commit**

```bash
cd /home/parlar_ai/dev/cyrius-rs-p1
git add cyrius-rs/src/phenotype.rs
git commit -m "fix(phenotype): exact copy counts multiply the single-copy activity

The loader stripped the ≥ from bucket rows such as *2x≥3, so they looked
like exact rows such as *41x3, and *41x4 borrowed the *41x3 activity
(0.75). *41x4/*10 was called Intermediate (1.0) instead of Normal (1.25).
Review #2.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

- [ ] **Step 6: Regression gate (G4).** No benchmark call has three or more copies of one allele, so no difference is allowed: both expected files are empty.

```bash
R=~/cyrius-research/review-fixes
: > $R/expected_task_1.tsv
: > $R/expected_cumulative_after_1.tsv
$R/gate.sh 1
# only after GATE PASS:
rm -rf $R/accepted && cp -r "$(cat $R/logs/task_1/results_dir.txt)" $R/accepted
```


---

### Task 2: BAM read errors are errors (review #6)

**Files:**
- Modify (required evidence, return `Err`):
  - `cyrius-rs/src/depth_calling/bin_count.rs:23-25, 48, 62, 65-83, 249-251`
  - `cyrius-rs/src/depth_calling/snp_count.rs:108-111`
  - `cyrius-rs/src/depth_calling/haplotype.rs:76-79`
  - `cyrius-rs/src/caller/call_variants.rs:286-288`
- Modify (optional evidence, stop and report unavailable):
  - `cyrius-rs/src/caller/phase_disambiguate.rs:265-345, 352-398, 404-440, 678-683`
  - `cyrius-rs/src/caller/diplotype_caller.rs:338-341`
  - `cyrius-rs/src/caller/d7_depth.rs:83-131, 164-167, 182-187`
  - `cyrius-rs/src/caller/kmer_validation.rs:167-222, 465-480`
  - `cyrius-rs/src/caller/clip_evidence.rs:101-104`
- Modify (call sites of the changed public signatures):
  - `cyrius-rs/examples/phase_audit.rs:67` — the only out-of-crate caller of `count_pair`. `cargo build` does not compile `examples/` by default, so this breaks only under `--all-targets`; G4 runs clippy with `--all-targets`, which is what catches it.
- Test: `cyrius-rs/src/depth_calling/bin_count.rs` (tests module), `cyrius-rs/src/caller/phase_disambiguate.rs` (tests module)

**Interfaces:**
- Produces:
  - `check_phasing(...) -> Result<PhaseCounts, String>`
  - `pair_evidence(...) -> Result<Option<PairEvidence>, String>`
  - `count_pair(...) -> Result<Option<PairEvidence>, String>`
  - `count_site(...) -> Option<(usize, usize, usize)>` (same type, but now `None` on a fetch or read error)
- Not touched: the three `Err(_) => continue` lines that parse text files (`diplotype_caller.rs:161`, `:166`, `phase_disambiguate.rs:38`).

- [ ] **Step 1: Write the failing tests**

Add to the tests module of `cyrius-rs/src/depth_calling/bin_count.rs`:

```rust
    /// An indexed BAM of 2,000 mapped 100 bp reads on chr22, starting at 0, 2, 4, ...
    fn write_test_bam(path: &std::path::Path) {
        let mut header = bam::header::Header::new();
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr22");
        sq.push_tag(b"LN", 50_818_468);
        header.push_record(&sq);
        {
            let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
            let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(100)]);
            for i in 0..2000i64 {
                let seq: Vec<u8> = (0..100).map(|j| b"ACGT"[((i + j) % 4) as usize]).collect();
                let mut r = bam::Record::new();
                r.set(format!("r{i}").as_bytes(), Some(&cigar), &seq, &[30u8; 100]);
                r.set_tid(0);
                r.set_pos(i * 2);
                r.set_mapq(60);
                r.set_flags(0);
                writer.write(&r).unwrap();
            }
        }
        bam::index::build(path, None::<&std::path::PathBuf>, bam::index::Type::Bai, 1).unwrap();
    }

    /// Flip one byte of the CRC32 of BGZF block `block` (block 0 holds the
    /// header). Offsets and the index stay valid, so a fetch succeeds and the
    /// failure only shows up while the records are read.
    fn corrupt_bgzf_block(path: &std::path::Path, block: usize) {
        let mut bytes = std::fs::read(path).unwrap();
        let mut off = 0usize;
        for _ in 0..block {
            off += u16::from_le_bytes([bytes[off + 16], bytes[off + 17]]) as usize + 1;
        }
        let size = u16::from_le_bytes([bytes[off + 16], bytes[off + 17]]) as usize + 1;
        assert!(size > 28, "block {block} is the EOF marker; write more reads");
        bytes[off + size - 8] ^= 0xFF;
        std::fs::write(path, bytes).unwrap();
    }

    fn test_bam_dir(tag: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("cyrius_bamio_{}_{}", tag, std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn a_corrupt_block_is_an_error_not_a_short_count() {
        let dir = test_bam_dir("corrupt");
        let path = dir.join("corrupt.bam");
        write_test_bam(&path);
        corrupt_bgzf_block(&path, 2);
        let region: Region = ("chr22".to_string(), 0, 5000, "test".to_string());
        let count = get_read_count(&mut bam::IndexedReader::from_path(&path).unwrap(), &region, 0);
        let mapq0 = get_mapq0_counts(&mut bam::IndexedReader::from_path(&path).unwrap(), &region);
        let length = get_read_length(&mut bam::IndexedReader::from_path(&path).unwrap());
        let _ = std::fs::remove_dir_all(&dir);
        assert!(count.is_err(), "get_read_count gave {count:?}");
        assert!(mapq0.is_err(), "get_mapq0_counts gave {mapq0:?}");
        assert!(length.is_err(), "get_read_length gave {length:?}");
    }

    #[test]
    fn the_same_bam_uncorrupted_counts_every_read() {
        // Control: the test above must fail because of the corrupt block.
        let dir = test_bam_dir("intact");
        let path = dir.join("ok.bam");
        write_test_bam(&path);
        let region: Region = ("chr22".to_string(), 0, 5000, "test".to_string());
        let count = get_read_count(&mut bam::IndexedReader::from_path(&path).unwrap(), &region, 0);
        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(count, Ok(2000));
    }
```

Add to the tests module of `cyrius-rs/src/caller/phase_disambiguate.rs`:

```rust
    #[test]
    fn check_phasing_reports_a_read_error_instead_of_partial_counts() {
        let dir = fixture_dir("readerr");
        // Enough fragments that the BAM spans several BGZF blocks.
        let mut records = Vec::new();
        for i in 0..1500 {
            records.extend(fixture_fragment(&format!("f{i}"), true, true, 60, 60));
        }
        let path = write_fixture_bam(&dir, "readerr.bam", records);
        let mut bytes = std::fs::read(&path).unwrap();
        let size0 = u16::from_le_bytes([bytes[16], bytes[17]]) as usize + 1;
        let size1 = u16::from_le_bytes([bytes[size0 + 16], bytes[size0 + 17]]) as usize + 1;
        let off2 = size0 + size1;
        let size2 = u16::from_le_bytes([bytes[off2 + 16], bytes[off2 + 17]]) as usize + 1;
        assert!(size2 > 28, "too few records for a third data block");
        bytes[off2 + size2 - 8] ^= 0xFF;
        std::fs::write(&path, bytes).unwrap();
        let lookup = fixture_lookup();
        let early = lookup.resolve("g.early").unwrap();
        let late = lookup.resolve("g.late").unwrap();
        let mut reader = bam::IndexedReader::from_path(&path).unwrap();
        let counts = check_phasing(&mut reader, "chr22", &early, &late, false);
        let _ = std::fs::remove_dir_all(&dir);
        assert!(counts.is_err(), "check_phasing gave {counts:?}");
    }
```

- [ ] **Step 2: Run the tests and check that they fail as expected**

Run: `cd /home/parlar_ai/dev/cyrius-rs-p1/cyrius-rs && cargo test --release --locked corrupt 2>&1 | grep -E "^test |error\[|gave"`
Expected:
- `a_corrupt_block_is_an_error_not_a_short_count` FAILS with `get_read_count gave Ok(<n>)` and n < 2000;
- `the_same_bam_uncorrupted_counts_every_read` PASSES.

Run: `cargo test --release --locked check_phasing_reports 2>&1 | grep -E "error\[E"`
Expected: a compile error, because `is_err` doesn't exist on `PhaseCounts`. That is the failing state for the phasing test.

- [ ] **Step 3: Required evidence returns `Err`**

`bin_count.rs`, `get_read_count` (line 23):
```rust
    while let Some(result) = reader.read(&mut record) {
        if let Err(e) = result {
            return Err(format!("BAM read failed in {}: {}", what, e));
        }
```
`bin_count.rs`, `get_read_length` (line 48) and `get_read_length_from_region` (line 62): replace `Ok(get_read_length_from_fetched(reader))` with `get_read_length_from_fetched(reader)`.
`bin_count.rs`, `get_read_length_from_fetched`:
```rust
fn get_read_length_from_fetched(reader: &mut bam::IndexedReader) -> Result<f64, String> {
    let mut read_lengths = Vec::new();
    let mut record = bam::Record::new();
    let mut counter = 0;
    while let Some(result) = reader.read(&mut record) {
        if let Err(e) = result {
            return Err(format!("BAM read failed while estimating read length: {}", e));
        }
        read_lengths.push(record.seq_len() as f64);
        counter += 1;
        if counter > 2000 {
            break;
        }
    }
    if read_lengths.is_empty() {
        return Ok(150.0);
    }
    Ok(stats::median(&read_lengths))
}
```
`bin_count.rs`, `get_mapq0_counts` (line 249):
```rust
    while let Some(result) = reader.read(&mut record) {
        if let Err(e) = result {
            return Err(format!("BAM read failed in {}: {}", what, e));
        }
```
`snp_count.rs`, `get_reads_by_region` (line 108):
```rust
            let pileup = match pileup_result {
                Ok(p) => p,
                Err(e) => return Err(format!("BAM pileup failed at {}:{}: {}", nchr, snp_position, e)),
            };
```
`haplotype.rs`, `get_bases_per_read` (line 76):
```rust
                let pileup = match pileup_result {
                    Ok(p) => p,
                    Err(e) => return Err(format!("BAM pileup failed while reading haplotypes: {}", e)),
                };
```
`call_variants.rs`, `get_allele_counts_var42128936` (line 286):
```rust
    while let Some(result) = reader.read(&mut record) {
        if let Err(e) = result {
            return Err(format!("BAM read failed at g.42128936: {}", e));
        }
```

- [ ] **Step 4: Optional evidence stops, and its step makes no decision**

`phase_disambiguate.rs`, `check_phasing`: change the return type to `Result<PhaseCounts, String>`.
- `None => return Ok(PhaseCounts::default()),` for the missing contig.
- Each `.is_err() { return PhaseCounts::default(); }` after a fetch becomes `.is_err() { return Err(format!("BAM fetch failed at {}:{}", nchr, pos1.min(pos2))); }`.
- Each record loop:
```rust
                let record = match record_result {
                    Ok(r) => r,
                    Err(e) => return Err(format!("BAM read failed at {}:{}: {}", nchr, pos1.min(pos2), e)),
                };
```
- End the function with `Ok(counts)`.

`pair_evidence` and `count_pair`: change the return type to `Result<Option<PairEvidence>, String>`. In `pair_evidence`:
- the two lookup misses return `Ok(None)`;
- `let counts = check_phasing(reader, nchr, &rv_a, &rv_b, use_read_pairs)?;`;
- the final value is `Ok(Some(PairEvidence { ... }))`.

In `disambiguate` (line 678):
```rust
        let evidence =
            match pair_evidence(&mut reader, nchr, va, vb, variant_lookup, use_read_pairs) {
                Ok(Some(e)) => e,
                Ok(None) => continue,
                Err(e) => {
                    log::warn!("phase_disambiguate: {e}; phase left unresolved");
                    return PhaseOutcome { resolved: None, pairs };
                }
            };
```
`count_site`:
- the fetch failure returns `None`;
- the record loop uses `Err(_) => return None,`.

Existing tests that call `check_phasing` or `count_pair` (lines 930, 933, 1066, 1086, 1111, 1250, 1253, 1254, 1255) get `.unwrap()` on the `Result`. The compiler lists each one.

`diplotype_caller.rs`, `collect_reads_dual_region` (line 338):
```rust
            let mut record = match result {
                Ok(r) => r,
                Err(e) => {
                    log::warn!("resolver: BAM read failed ({e}); leaving the call unresolved");
                    return Vec::new();
                }
            };
```
`d7_depth.rs`:
- `count_reads_in_region` returns `Option<u64>`. A fetch failure gives `return None;`, the record loop uses `if result.is_err() { return None; }`, and the function ends with `Some(count)`.
- `measure_region` returns `Option<RegionDepth>`, with `let read_count = count_reads_in_region(reader, tid, start, end, min_mapq)?;` and `Some(RegionDepth { ... })`.
- In `compute_d7_depth`: the read-length loop uses `if result.is_err() { return None; }`, and each `measure_region(...)` call at lines 182–187 ends in `?`.

`kmer_validation.rs`, `count_kmers_in_region`:
- returns `Option<KmerCounts>`;
- the missing-chromosome and fetch-failure branches log as now, then `return None;`;
- the record loop uses `Err(_) => return None,`;
- the function ends with `Some(KmerCounts { ... })`.

At the two callers (lines 465 and 474):
```rust
    let counts_d6 = match count_kmers_in_region(&mut reader, chrom, 42123000, 42132000, &lookup) {
        Some(c) => c,
        None => return KmerValidationResult::unavailable("BAM read failed in the CYP2D6 region"),
    };
    let counts_d7 = match count_kmers_in_region(&mut reader, chrom, 42135000, 42146000, &lookup) {
        Some(c) => c,
        None => return KmerValidationResult::unavailable("BAM read failed in the CYP2D7 region"),
    };
```
`clip_evidence.rs` (line 101):
```rust
        let record = match result {
            Ok(r) => r,
            Err(_) => return None,
        };
```

- [ ] **Step 5: Run all tests**

Run: `cargo test --release --locked 2>&1 | grep -E "^test result|FAILED|panicked"`
Expected: every `test result` is `ok`, with 0 failed, and both corrupt-BAM tests pass.

- [ ] **Step 6: Check that no BAM read error is still skipped**

Run: `grep -rnE "is_err\(\) *\{ *$|Err\(_\) *=> *continue" src --include=*.rs`
At `3807326` this command gives 19 hits (measured). **Measured afterwards: 7 hits**, all accounted for:
- Three are text parsing, not BAM reads: `diplotype_caller.rs:161,166` (`cols[...].parse()`) and `phase_disambiguate.rs:38` (`fields[1].parse`).
- Two are test function *names* in `bin_count.rs` ending in `_is_err() {`.
- Two are BAM reads in `d7_depth.rs:100,165` of the form `if result.is_err() { return None; }`.

**The expected count was first written as 5, and that was wrong** — a prediction about the *form* the d7-depth fix would take, not a missed site. It is corrected here rather than argued with.

The pattern over-matches, so the count alone is not the gate. `return None` discards the partial `count`/`lengths` and reports the optional evidence as unavailable, which is exactly what this task asks of an optional path; `continue` would keep counting past a bad record, which is the defect. The grep cannot tell them apart, so run the sharper check as well:

```bash
grep -rnA 2 -E "is_err\(\) *\{ *$|Err\(_\) *=> *" src --include=*.rs \
  | grep -E "continue;" | grep -v "_is_err"
```
**Measured: no output** — no read error is skipped-and-continued anywhere. That is the property the task is actually buying, and an empty result is the pass condition. If this prints anything, a read error is still being stepped over.

Read each hit of the broad grep anyway, and classify it into one of the three groups above. A hit that fits none of them is a finding.

- [ ] **Step 7: Commit**

```bash
git add cyrius-rs/src
git commit -m "fix(io): a BAM read error stops the step instead of giving a partial count

Read and pileup errors after a successful fetch were skipped with
continue, so a corrupt block gave a short count that was then used for
normalization and calling. Required steps now return an error; optional
steps (phasing, resolver, d7-depth, k-mer, clip evidence) stop and decide
nothing. Review #6.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

- [ ] **Step 8: Regression gate (G4).** No change is allowed.

```bash
R=~/cyrius-research/review-fixes
: > $R/expected_task_2.tsv
: > $R/expected_cumulative_after_2.tsv
$R/gate.sh 2
# only after GATE PASS:
rm -rf $R/accepted && cp -r "$(cat $R/logs/task_2/results_dir.txt)" $R/accepted
```


---

### Task 3: A site with no call is not a reference call (review #1)

> **SCOPE — read before implementing. Only 3a (reporting) belongs in this
> plan. The inference fix has outgrown it and needs its own plan and its own
> gate.**
>
> **The direction is settled, and it is neither of the two options this plan
> offered before.** Both of those picked a single genotype and then argued
> about how loudly to warn. The reviewer's third approach is better and is
> what we build: keep the uncertainty inside allele matching.
>
> 1. An uncallable dosage is **unknown**, not zero. A confident zero excludes
>    variants; `None` excludes nothing.
> 2. Mask all unknown sites **together** and keep **every** compatible
>    candidate, applying the existing structural constraints and whatever
>    fragment evidence passes the Task 5 / Task 6 quality gates. Match on raw
>    sub-alleles, simplify names only for reporting.
> 3. Report a clinical result **only when the candidate set establishes it**.
>    If every compatible candidate carries the same activity score and
>    phenotype, print it — that is the common case and must not be suppressed.
>    If they differ, or any is unknown, the interpretation is unresolved and
>    the alternatives are listed.
> 4. No-call is the **fallback**, for when nothing fits or enumeration cannot
>    be completed.
>
> Why this is right: it removes the need to choose between "Normal" and
> "Intermediate" on evidence that establishes neither, *without* blanking the
> clinical answer whenever the surviving candidates happen to agree.
>
> **Why it is not a step in this plan.** It changes `match_star` — the core
> matcher every call goes through — rather than adding a filter after it.
> Under the judgment gate that is a new mechanism, so it needs its own Gate A
> answers, its own trivial baseline, and its own cheapest-kill experiment
> before any code. Bolting it onto a plan whose other six tasks are
> self-contained bug fixes would also break this plan's own property that each
> task is independently gateable.
>
> **Gate A, answered here so the next plan can start from it:**
>
> - *Governing principle:* evidence that was never observed must not narrow a
>   call. Unknown is a third state, and it propagates to the clinical answer
>   only when it cannot change it.
> - *What would kill it:* on fully observed inputs the candidate set must
>   collapse to exactly today's single call on all 53 benchmark samples. One
>   sample that gains a spurious alternative kills it — the masking is then
>   admitting candidates the data already excludes.
> - *Already refuted here?* No. This is not panel tie-breaking (H15) or
>   priors (H19); it adds no new evidence source and no prior. It only stops
>   treating absence as a reference observation.
> - *Trivial baseline to beat:* today's behaviour plus 3a's report. The new
>   mechanism has to earn its place by changing an answer that the report
>   alone would leave wrong.
> - *Input assumptions to verify first:* that `None` really means "not
>   observed" at every site 3a flags, and not "observed and rejected" — the
>   two must not be conflated, and `call_cn_var` reaches `None` by more than
>   one route (zero reads, zero probability mass, posterior below cutoff).
>
> **Numbers that are NOT established yet.** The reviewer's prototype reports
> five uncallable dosages and nine compatible main-allele combinations for the
> damaged NA19920 (including `*1/*4x2` Intermediate, `*1/*10x2` Normal,
> `*4/*4x2` Poor). I could not reproduce either figure: my own measurement
> (above) finds three sites at zero reads and 31 degraded, and the true
> `None` count needs 3a's code path. Treat 5 and 9 as the reviewer's
> prototype output, not as facts, until 3a and the new plan's own Gate B
> measure them. The *direction* they establish — that more than one genotype
> is compatible, with conflicting phenotypes — is what matters here, and my
> measurement supports it.
>
> **Tests the new plan must require** (from the reviewer, kept verbatim in
> intent):
> - an unknown dosage never acts as a reference call;
> - all unknown sites are masked jointly, not one at a time;
> - the true NA19920 genotype stays in the candidate set, established without
>   consulting its truth label;
> - conflicting or unknown clinical effects block a definitive interpretation;
> - fully observed calls keep their current behaviour, on all 53.
>
> **So in this plan:** 3a proceeds (it is a pure addition and the new plan
> depends on its `None` reporting). 3b and 3c are **withdrawn** — the filter
> and the confidence cap were containment for a wrong answer that the new
> approach does not produce. Task 7 records review #1 as *diagnosed and
> reported, fix deferred to its own plan*, with this measurement attached.

**Files:**
- Modify: `cyrius-rs/src/caller/call_variants.rs` (new helpers next to `get_called_variants`, line ~577)
- Modify: `cyrius-rs/src/types.rs:134-195` (new `D6Call` field)
- Modify: `cyrius-rs/src/main.rs:453-462, 657-668, 815-818, 1025-1046, ~1300` — the `None` dosages are collected and reported. The `1196-1240` block (filter downgrade and confidence cap) is **not** touched: that was 3c.
- Test: `cyrius-rs/src/caller/call_variants.rs` (tests module)
- Create: `~/cyrius-research/review-fixes/make_na19920_missing_splice.sh` — the reproduction script. Still built here (3b Step 1) even though the rest of 3b is withdrawn: the reviewer's copy lives in `/tmp` and will vanish, and the next plan's Gate B needs this input.

**Interfaces:**
- Produces:
  - `pub fn uncallable_variants(var_list: &[String], cn_calls: &[Option<u32>], starting_index: usize) -> Vec<String>`
  - `D6Call.uncallable_target_sites: Option<Vec<String>>`
  - JSON field `Uncallable_target_sites`
  - *Withdrawn with 3c, do not build:* `UNCALLABLE_FILTER`, `downgrade_filter_for_uncallable_sites`, and the confidence cap (`cap_below_high`, `label_for_score`, `Confidence_uncapped`). The next plan decides what, if anything, replaces them.

#### 3a — Report the uncallable sites (no change to the call)

- [ ] **Step 1: Write the failing test** in the tests module of `call_variants.rs`:

```rust
    #[test]
    fn uncallable_variants_lists_sites_without_a_call_and_not_zero_calls() {
        let names: Vec<String> = ["g.a", "g.b", "g.c", "g.d"].iter().map(|s| s.to_string()).collect();
        let calls = [Some(0), None, Some(1)];
        assert_eq!(uncallable_variants(&names, &calls, 0), vec!["g.b".to_string()]);
        assert_eq!(uncallable_variants(&names, &[None], 3), vec!["g.d".to_string()]);
        assert!(uncallable_variants(&names, &[Some(0), Some(2)], 0).is_empty());
    }
```

- [ ] **Step 2: Run it**

Run: `cargo test --release --locked uncallable_variants 2>&1 | grep -E "error\[E0425\]|cannot find function"`
Expected: `cannot find function 'uncallable_variants'`.

- [ ] **Step 3: Implement the helper** directly below `get_called_variants` in `call_variants.rs`:

```rust
/// Target variants whose dosage could not be called (`None`): no usable
/// reads, or no copy number reached its posterior cutoff. A missing call is
/// not a reference call, so these are reported instead of being dropped the
/// way `get_called_variants` drops them.
pub fn uncallable_variants(
    var_list: &[String],
    cn_calls: &[Option<u32>],
    starting_index: usize,
) -> Vec<String> {
    cn_calls
        .iter()
        .enumerate()
        .filter(|(_, c)| c.is_none())
        .map(|(i, _)| var_list[i + starting_index].clone())
        .collect()
}
```

- [ ] **Step 4: Plumb the field.**

In `types.rs`, add to `D6Call` after `rep7_mapq0_fraction`:
```rust
    /// Target variants with no dosage call (see `call_variants::uncallable_variants`).
    pub uncallable_target_sites: Option<Vec<String>>,
```
Add `uncallable_target_sites: None,` to the `D6Call` literal in `D6Call::no_call` and to the two full literals at `main.rs:453` and `main.rs:657`.

In `main.rs`, after line 817 (`let called_var_homo = ...`):
```rust
    let mut uncallable_target_sites = uncallable_variants(&params.var_list, &cn_call_var, 0);
    uncallable_target_sites.extend(uncallable_variants(&params.var_list, &cn_call_var_homo, cn_call_var.len()));
```
Then add `uncallable_variants` to the `call_variants` import list at `main.rs:8`, and add `uncallable_target_sites: Some(uncallable_target_sites),` to the final `D6Call` literal (`main.rs:1025`). In the JSON object, after `"Header_drops_multimappers": header_drops_multimappers,`:
```rust
        "Uncallable_target_sites": cyp2d6_call.uncallable_target_sites,
```

- [ ] **Step 5: Run the tests and build:** `cargo test --release --locked 2>&1 | grep -E "^test result|FAILED"`. Expected: all ok.

- [ ] **Step 6: Commit** (`feat(output): report target sites without a dosage call`, with the G2 trailer).

- [ ] **Step 7: Regression gate (G4).** The compared fields don't include the new JSON field, so no change is allowed.

```bash
R=~/cyrius-research/review-fixes
: > $R/expected_task_3a.tsv
: > $R/expected_cumulative_after_3a.tsv
$R/gate.sh 3a
# only after GATE PASS:
rm -rf $R/accepted && cp -r "$(cat $R/logs/task_3a/results_dir.txt)" $R/accepted
```


#### 3b — one step kept, the rest withdrawn

> The reproduction below is still built: the reviewer's BAM lives in `/tmp`
> and will vanish, and the next plan's Gate B needs exactly this input. The
> necessity analysis that followed it is withdrawn — the candidate-set
> matcher answers that question structurally instead.

- [ ] **Step 1 (still executed): Rebuild the NA19920 reproduction** (the reviewer's copy in `/tmp` is temporary). Write `~/cyrius-research/review-fixes/make_na19920_missing_splice.sh`:

```bash
#!/usr/bin/env bash
# Remove only the alignments overlapping g.42128945 (the *4 splice variant)
# from NA19920. Measured 2026-11-04: 68 alignments removed, 1,298,717 left.
set -euo pipefail
O=/home/parlar_ai/dev/cyrius-rs/benchmark/NA19920_pgx_regions.bam
D=~/cyrius-research/review-fixes/na19920_missing_splice
mkdir -p "$D"
samtools view -h "$O" \
  | awk 'BEGIN{FS=OFS="\t"} NR==FNR{k[$1 FS $2 FS $4]=1; next} /^@/ || !(($1 FS $2 FS $4) in k)' \
      <(samtools view "$O" chr22:42128945-42128945) - \
  | samtools view -b -o "$D/NA19920_missing_splice.bam"
samtools index "$D/NA19920_missing_splice.bam"
echo "total $(samtools view -c "$D/NA19920_missing_splice.bam") at-site $(samtools view -c "$D/NA19920_missing_splice.bam" chr22:42128945-42128945)"
```
Run: `bash ~/cyrius-research/review-fixes/make_na19920_missing_splice.sh`
Expected: `total 1298717 at-site 0`


<details><summary>Withdrawn content</summary>

#### 3b — Measure before deciding (locked rule)

- [ ] **Step 2: Count the flagged benchmark samples**

```bash
NEW=$(cat ~/cyrius-research/review-fixes/logs/task_3a/results_dir.txt)
python3 - "$NEW" <<'EOF'
import glob, json, os, sys
flagged = {}
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*.json"))):
    r = next(iter(json.load(open(f)).values()))
    if r.get("Uncallable_target_sites"):
        flagged[os.path.basename(f)[:-5]] = r["Uncallable_target_sites"]
print(len(flagged), "of 53 samples have uncallable target sites")
for s, v in flagged.items():
    print(s, len(v), v[:6])
EOF
```

- [ ] **Step 3: Decide whether each flagged site could have changed the call.**

  The first draft of this plan said: "Every benchmark call is concordant, so
  each flagged sample is a false alarm." That inference is wrong (review #8). A
  genotype can match truth while the evidence was never sufficient to separate
  it from a rival — that is exactly what NA19920 with its splice-site reads
  removed demonstrates: the call becomes `*1/*10x2`, and had truth happened to
  be `*1/*10x2` the same missing evidence would have looked like a success.
  **Concordance and evidence sufficiency are separate questions and are
  recorded separately.**

  The criterion is necessity, not concordance: does the uncallable site
  actually distinguish the called genotype from another candidate the matcher
  would have considered?

```bash
NEW=$(cat ~/cyrius-research/review-fixes/logs/task_3a/results_dir.txt)
python3 ~/cyrius-research/review-fixes/site_necessity.py "$NEW" \
  /home/parlar_ai/dev/cyrius-rs/cyrius-rs/data/star_table.txt
```

  `site_necessity.py` (written in this step) reports, per flagged sample:
  `sample`, `n_uncallable`, and `n_decisive` — the number of uncallable sites
  at which at least one other star allele in the same CNV group differs from
  the called genotype's own dosage, i.e. sites whose true value could have
  selected a different candidate. It must also print, for each decisive site,
  the rival allele that differs there, so the number can be checked by hand.

  Prove the script can fail before trusting it: run it on the NA19920
  reproduction built in Step 1, where `g.42128945C>T` is known to be decisive (it
  is the difference between `*4` and `*10`). Expected: `n_decisive >= 1`. If it
  reports 0 there, the script is broken and the gate means nothing.

- [ ] **Step 4: Apply the rule fixed now, before the result is seen.**
  - **No flagged sample has `n_decisive >= 1`:** continue to 3c, and record in
    Task 7 that on this benchmark the filter is informational only.
  - **Any flagged sample has `n_decisive >= 1`:** STOP and report. A benchmark
    call then rests on a site that was never observed, which is a finding in its
    own right and changes what D1 has to do. Do not continue on the assumption
    that a concordant call is a safe one.
  - Either way, report the flagged list and the concordance of each flagged
    sample **as two separate columns**. Neither is evidence about the other.

</details>

#### 3c — WITHDRAWN

> Superseded. The filter and the confidence cap were containment for a wrong
> genotype that the candidate-set approach does not produce. `UNCALLABLE_FILTER`,
> `downgrade_filter_for_uncallable_sites` and the confidence cap are therefore
> **not built** — 3a's reporting field is all this plan adds. **Not executed.**

<details><summary>Withdrawn content</summary>

#### 3c — The filter (withdrawn; was: only after 3b passes and the user confirms D1)

- [ ] **Step 1: Write the failing test** in the tests module of `call_variants.rs`:

```rust
    #[test]
    fn an_uncallable_site_removes_pass_and_fuzzy_match_only() {
        let sites = vec!["g.42128945C>T".to_string()];
        let pass = Some("PASS".to_string());
        assert_eq!(downgrade_filter_for_uncallable_sites(pass.clone(), &sites).as_deref(), Some(UNCALLABLE_FILTER));
        assert_eq!(downgrade_filter_for_uncallable_sites(Some("Fuzzy_match".into()), &sites).as_deref(), Some(UNCALLABLE_FILTER));
        assert_eq!(downgrade_filter_for_uncallable_sites(Some("More_than_one_possible_genotype".into()), &sites).as_deref(), Some("More_than_one_possible_genotype"));
        assert_eq!(downgrade_filter_for_uncallable_sites(pass.clone(), &[]), pass);
        assert_eq!(downgrade_filter_for_uncallable_sites(None, &sites), None);
    }
```

- [ ] **Step 2: Run it.** Expected: `cannot find function 'downgrade_filter_for_uncallable_sites'`.

- [ ] **Step 3: Implement** below `uncallable_variants`:

```rust
/// Filter for a call made while some target sites had no dosage call.
pub const UNCALLABLE_FILTER: &str = "Uncallable_target_sites";

/// A PASS or Fuzzy_match call made with uncallable target sites is not a
/// PASS: the missing sites were treated as reference by the matcher. Other
/// filters already say the call is not clean and are kept.
pub fn downgrade_filter_for_uncallable_sites(
    filter: Option<String>,
    uncallable: &[String],
) -> Option<String> {
    match filter.as_deref() {
        Some("PASS" | "Fuzzy_match") if !uncallable.is_empty() => Some(UNCALLABLE_FILTER.to_string()),
        _ => filter,
    }
}
```

In `main.rs`, directly before the multi-mapper guard (`let header_drops_multimappers =`, line ~1199):
```rust
    let uncallable = cyp2d6_call.uncallable_target_sites.clone().unwrap_or_default();
    let filtered = cyrius_rs::caller::call_variants::downgrade_filter_for_uncallable_sites(
        cyp2d6_call.filter.clone(),
        &uncallable,
    );
    if filtered != cyp2d6_call.filter {
        log::warn!(
            "Sample {}: {} target site(s) have no dosage call ({}); a missing site is not a \
             reference call, so the call is not a PASS.",
            sample_id, uncallable.len(), uncallable.join(", "),
        );
        cyp2d6_call.filter = filtered;
    }
```
First, in `confidence.rs`, give the label mapping one home and add the cap.
The thresholds are moved, not changed (`0.85` / `0.5` as measured at
`3807326`), so G3 is not touched:

```rust
/// Label thresholds, named so the cap cannot drift away from them.
pub const CONFIDENCE_HIGH: f64 = 0.85;
pub const CONFIDENCE_MEDIUM: f64 = 0.5;

/// The score the cap reports. Any value below `CONFIDENCE_HIGH` would serve;
/// this is a reporting ceiling, not a decision threshold, so it is not a
/// tuned parameter under G3.
pub const CAPPED_SCORE: f64 = 0.84;

/// The label for a score. `compute_confidence` must use this too, instead of
/// its own inline comparison, so there is exactly one mapping.
pub fn label_for_score(score: f64) -> &'static str {
    if score >= CONFIDENCE_HIGH {
        "HIGH"
    } else if score >= CONFIDENCE_MEDIUM {
        "MEDIUM"
    } else {
        "LOW"
    }
}

/// Stop a score being HIGH, moving the number and the label together.
///
/// Capping the label alone would let the output read `0.95 (MEDIUM)` (review
/// #4). The uncapped value is not lost: the caller reports it as
/// `Confidence_uncapped`.
pub fn cap_below_high(c: &mut ConfidenceScore) -> Option<f64> {
    if c.score >= CONFIDENCE_HIGH {
        let was = c.score;
        c.score = CAPPED_SCORE;
        c.label = label_for_score(c.score).to_string();
        Some(was)
    } else {
        None
    }
}
```

Then in `main.rs`, change `let confidence = ...compute_confidence(` to
`let mut confidence = ...compute_confidence(`. After that call, add:

```rust
    // Missing evidence caps confidence. The trigger is the uncallable-site
    // list itself, NOT the filter string: `downgrade_filter_for_uncallable_sites`
    // deliberately leaves other filters in place, so a call that is both
    // ambiguous and missing target evidence keeps a filter like
    // `More_than_one_possible_genotype` — and keying the cap on the filter
    // would then leave it at HIGH (review #4).
    let confidence_uncapped = if uncallable.is_empty() {
        None
    } else {
        cyrius_rs::caller::confidence::cap_below_high(&mut confidence)
    };
```

Add two JSON fields next to `Confidence`: `"Confidence_uncapped": confidence_uncapped`
(null when nothing was capped) and `"Confidence_capped_by": if confidence_uncapped.is_some() { "Uncallable_target_sites" } else { "n/a" }`.

Note the ordering constraint: `uncallable` is bound in 3c's earlier block, which
must therefore stay above this one.

- [ ] **Step 4: Run the tests.** Expected: all ok.

- [ ] **Step 5: End-to-end check on the reproduction (must turn red → not PASS)**

```bash
cd /home/parlar_ai/dev/cyrius-rs-p1/cyrius-rs && cargo build --release --locked 2>&1 | tail -1
W=$(mktemp -d); B=target/release/cyrius-rs
$B -i /home/parlar_ai/dev/cyrius-rs/benchmark/NA19920_pgx_regions.bam -o $W/orig --id NA19920 2>/dev/null
$B -i ~/cyrius-research/review-fixes/na19920_missing_splice/NA19920_missing_splice.bam -o $W/mod --id NA19920 2>/dev/null
for d in orig mod; do python3 -c "
import json; r=next(iter(json.load(open('$W/$d/NA19920.json')).values()))
print('$d', r['Genotype'], r['Filter'], r['Confidence'], r['Confidence_label'], r['Uncallable_target_sites'])"; done
```
Expected:
- `orig *1/*4x2 PASS 0.87 HIGH []`
- `mod *1/*10x2 Uncallable_target_sites <score> MEDIUM-or-LOW [..., 'g.42128945C>T', ...]`

If `mod` still shows PASS or HIGH, the fix failed. Stop.

- [ ] **Step 6: Commit** (`fix(calling): a target site without a dosage call removes PASS and HIGH`, with the review #1 reference and the G2 trailer).

- [ ] **Step 7: Regression gate (G4).** This is the first task that changes
  benchmark output, so its expected file names the new values — not just the
  fields (G8). Every value is known **before** the run, from the baseline:
  `Filter` becomes `Uncallable_target_sites` for a sample whose baseline filter
  was `PASS` or `Fuzzy_match`; `Confidence` becomes `0.84` exactly when the
  baseline score was `>= 0.85`; `Confidence_label` becomes `MEDIUM` in that
  case. Write the file from the baseline, then run — do not write it from the
  new output, which would make the gate self-confirming.

```bash
R=~/cyrius-research/review-fixes
python3 - > $R/expected_task_3c.tsv <<'EOF'
import glob, json, os
# Flagged samples come from 3b; paste that list here before running.
FLAGGED = []   # e.g. ["NA19317"]
for s in FLAGGED:
    r = next(iter(json.load(open(f"{os.path.expanduser('~')}/cyrius-research/review-fixes/baseline/{s}.json")).values()))
    if r["Filter"] in ("PASS", "Fuzzy_match"):
        print(f"{s}\tFilter\tUncallable_target_sites")
    if float(r["Confidence"]) >= 0.85:
        print(f"{s}\tConfidence\t0.84")
        print(f"{s}\tConfidence_label\tMEDIUM")
EOF
cat $R/expected_task_3c.tsv
cp $R/expected_task_3c.tsv $R/expected_cumulative_after_3c.tsv
$R/gate.sh 3c
# only after GATE PASS:
rm -rf $R/accepted && cp -r "$(cat $R/logs/task_3c/results_dir.txt)" $R/accepted
```
  If 3b found no flagged sample, both files are empty and nothing may change.

---

</details>

---

### Task 4: Disagreeing mates make a site uninformative (review #3)

**Files:**
- Modify: `cyrius-rs/src/caller/phase_disambiguate.rs:224-252` (`Fragment`, `observe_read`), `:334-345` (tally)
- Test: `cyrius-rs/src/caller/phase_disambiguate.rs` (tests module)

**Interfaces:**
- Produces: `enum SiteObs { Unseen, Allele(bool), Conflict }`, private to the module. `Fragment.site1`/`Fragment.site2` replace `allele1`/`allele2`.

- [ ] **Step 1: Write the failing test**

```rust
    /// A 150 bp read starting at 990 covers both fixture sites (1000, 1100).
    fn spanning_record(qname: &str, early_alt: bool, late_alt: bool, qual: u8, flags: u16) -> bam::Record {
        spanning_record_quals(qname, early_alt, late_alt, vec![qual; 150], flags)
    }

    /// `spanning_record` with a per-base quality array, so a single base of a
    /// pattern can be made untrustworthy (review #5).
    fn spanning_record_quals(
        qname: &str,
        early_alt: bool,
        late_alt: bool,
        quals: Vec<u8>,
        flags: u16,
    ) -> bam::Record {
        let mut seq = vec![b'A'; 150];
        seq[(FIXTURE_POS_EARLY - 990) as usize] = if early_alt { b'C' } else { b'G' };
        seq[(FIXTURE_POS_LATE - 990) as usize] = if late_alt { b'T' } else { b'A' };
        let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(150)]);
        let mut record = bam::Record::new();
        record.set(qname.as_bytes(), Some(&cigar), &seq, &quals);
        record.set_tid(0);
        record.set_pos(990);
        record.set_mapq(60);
        record.set_mtid(0);
        record.set_flags(flags);
        record
    }

    #[test]
    fn disagreeing_mates_give_the_same_result_in_either_order() {
        // Each fragment: one mate ALT/ALT (cis), the other ALT/REF (trans),
        // both at the same coordinate. The late site conflicts.
        let order = |cis_first: bool| -> Vec<bam::Record> {
            let mut v = Vec::new();
            for i in 0..3 {
                let q = format!("f{i}");
                let cis = spanning_record(&q, true, true, 30, 0x1 | 0x40);
                let trans = spanning_record(&q, true, false, 30, 0x1 | 0x80);
                if cis_first { v.push(cis); v.push(trans); } else { v.push(trans); v.push(cis); }
            }
            v
        };
        let a = disambiguate_fixture("conflict_a", order(true));
        let b = disambiguate_fixture("conflict_b", order(false));
        assert_eq!(a.resolved, b.resolved, "record order changed the phase decision");
        assert_eq!(a.resolved, None, "a conflicting site is not evidence");
    }
```

- [ ] **Step 2: Run it.** Run: `cargo test --release --locked disagreeing_mates 2>&1 | grep -E "^test |panicked|left|right"`
  Expected: FAIL with `record order changed the phase decision`. The two orders give `Some("*B/*C")` and `Some("*1/*A")`.

- [ ] **Step 3: Implement.** Replace `Fragment` and `observe_read`:

```rust
/// What the reads of one fragment showed at one site.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
enum SiteObs {
    #[default]
    Unseen,
    /// ALT (`true`) or REF (`false`), with every read agreeing.
    Allele(bool),
    /// Reads of the same fragment disagree: the site is not evidence.
    Conflict,
}

impl SiteObs {
    /// Add one read's observation. Order-independent: any disagreement is a
    /// conflict for good, so no read can overwrite another.
    fn merge(self, seen: Option<bool>) -> SiteObs {
        match (self, seen) {
            (s, None) => s,
            (SiteObs::Unseen, Some(a)) => SiteObs::Allele(a),
            (SiteObs::Allele(old), Some(a)) if old == a => SiteObs::Allele(a),
            (SiteObs::Allele(_), Some(_)) | (SiteObs::Conflict, Some(_)) => SiteObs::Conflict,
        }
    }

    fn allele(self) -> Option<bool> {
        match self {
            SiteObs::Allele(a) => Some(a),
            _ => None,
        }
    }
}

/// What one fragment (all reads sharing a qname) showed at a pair's two sites.
#[derive(Debug, Default)]
struct Fragment {
    site1: SiteObs,
    site2: SiteObs,
    /// Highest MAPQ over every read of this qname that passed the filters.
    max_mapq: u8,
}

/// Merge one fetched read's observations at both sites into its fragment.
fn observe_read(
    fragments: &mut HashMap<String, Fragment>,
    record: &bam::Record,
    rv1: &ResolvedVariant,
    rv2: &ResolvedVariant,
) {
    let read_name = String::from_utf8_lossy(record.qname()).to_string();
    let seq = record.seq().as_bytes();
    let allele1 = check_allele_pileup(record, &seq, rv1);
    let allele2 = check_allele_pileup(record, &seq, rv2);
    let entry = fragments.entry(read_name).or_default();
    entry.site1 = entry.site1.merge(allele1);
    entry.site2 = entry.site2.merge(allele2);
    entry.max_mapq = entry.max_mapq.max(record.mapq());
}
```
In `check_phasing`'s tally, change `match (f.allele1, f.allele2) {` to `match (f.site1.allele(), f.site2.allele()) {`.

- [ ] **Step 4: Run all tests.** Expected: all ok, including the new test.

- [ ] **Step 5: Commit** (`fix(phasing): disagreeing reads of one fragment make the site uninformative`, with the review #3 reference and the G2 trailer).

- [ ] **Step 6: Regression gate (G4).** No change is allowed by this task. The cumulative file carries Task 3c's changes forward (G8) — comparing against the frozen baseline with an empty file here would fail on Task 3c's accepted output, which was the first draft's mistake. If HG00589 or NA18973 changes `Resolved_by` or `Genotype`, STOP and report; that would mean their decisive fragments contained conflicts.

```bash
R=~/cyrius-research/review-fixes
: > $R/expected_task_4.tsv
: > $R/expected_cumulative_after_4.tsv   # 3c withdrawn: this plan changes no benchmark output
$R/gate.sh 4
# only after GATE PASS:
rm -rf $R/accepted && cp -r "$(cat $R/logs/task_4/results_dir.txt)" $R/accepted
```


---

### Task 5: Base quality gates phasing evidence (review #4, phasing path)

**Files:**
- Modify: `cyrius-rs/src/caller/phase_disambiguate.rs:456-532` (`check_allele_pileup`, `read_bases_at_pos`)
- Test: `cyrius-rs/src/caller/phase_disambiguate.rs` (tests module; uses `spanning_record` from Task 4)

**Interfaces:**
- Produces, in `depth_calling/mod.rs`: `pub const PHASE_MIN_BASE_QUALITY: u8 = 13;` and `pub const QUAL_UNAVAILABLE: u8 = 0xFF;`, re-exported from `caller::phase_disambiguate` so that path still resolves. In `phase_disambiguate.rs`: `fn query_span_at_pos(record, seq_len, genomic_pos_0based, n) -> Option<(usize, usize)>` and the test helper `spanning_record_quals(qname, early_alt, late_alt, quals, flags)`.
- Layering note for Task 6: import these from `crate::depth_calling::{...}` inside `depth_calling/haplotype.rs`. Do **not** reach up into `crate::caller::` from `depth_calling` — nothing in that module does today.

**Threshold:** 13 is the default `--min-BQ` of `bcftools mpileup`. It is chosen before any benchmark run and must not be changed after one (G3).

- [ ] **Step 1: Write the failing tests**

```rust
    #[test]
    fn base_quality_zero_reads_do_not_decide_phase() {
        let records = (0..3).map(|i| spanning_record(&format!("q{i}"), true, true, 0, 0)).collect();
        assert_eq!(disambiguate_fixture("bq0", records).resolved, None);
    }

    #[test]
    fn the_same_reads_at_base_quality_30_do_decide_phase() {
        // Control: the test above must fail for quality, not for layout.
        let records = (0..3).map(|i| spanning_record(&format!("q{i}"), true, true, 30, 0)).collect();
        assert_eq!(disambiguate_fixture("bq30", records).resolved.as_deref(), Some("*1/*A"));
    }

    #[test]
    fn quality_12_is_below_the_gate_and_quality_13_is_not() {
        // The boundary, both sides. 13 is the threshold, so 13 must decide.
        let at = |q: u8| {
            let records = (0..3).map(|i| spanning_record(&format!("q{i}"), true, true, q, 0)).collect();
            disambiguate_fixture(&format!("bq{q}"), records).resolved
        };
        assert_eq!(at(12), None);
        assert_eq!(at(13).as_deref(), Some("*1/*A"));
    }

    #[test]
    fn reads_without_base_qualities_do_not_decide_phase() {
        // A BAM with qualities omitted stores 0xFF. That is "unknown", not
        // "excellent", so it must not decide (review #5).
        let records = (0..3)
            .map(|i| spanning_record(&format!("q{i}"), true, true, QUAL_UNAVAILABLE, 0))
            .collect();
        assert_eq!(disambiguate_fixture("bqmissing", records).resolved, None);
    }

    #[test]
    fn only_the_compared_bases_quality_matters() {
        // Low quality ON a discriminating base rejects; low quality on a
        // neighbouring base does not. This is what proves the gate reads the
        // compared span and not the whole read.
        let early = (FIXTURE_POS_EARLY - 990) as usize;
        let with_low_at = |offset: usize| {
            let records = (0..3)
                .map(|i| {
                    let mut quals = vec![30u8; 150];
                    quals[offset] = 2;
                    spanning_record_quals(&format!("q{i}"), true, true, quals, 0)
                })
                .collect();
            disambiguate_fixture(&format!("low{offset}"), records).resolved
        };
        assert_eq!(with_low_at(early), None, "low quality on the compared base must reject");
        assert_eq!(
            with_low_at(early + 20).as_deref(),
            Some("*1/*A"),
            "low quality on an unrelated base must not reject"
        );
    }
```

- [ ] **Step 1b: Write the multi-base-pattern test.** The fixture above has
single-base variants, so it cannot show the case the reviewer asked for: one
untrustworthy base *inside* a pattern matched over several bases (an indel
allele). That belongs to `check_allele_pileup`, which matches whole patterns.

```rust
    #[test]
    fn one_bad_base_inside_a_multi_base_pattern_rejects_the_whole_pattern() {
        // A 4-base ALT pattern where only the third base is poor. The pattern
        // still matches the sequence, so without the gate it would be
        // accepted; every compared base must reach the threshold.
        let rv = ResolvedVariant {
            pos_0based: 1000,
            alt_bases: "ACGT".to_string(),
            ref_bases: "A".to_string(),
        };
        let seq = b"ACGT".to_vec();
        let good = pileup_record_with_quals(1000, &seq, &[30, 30, 30, 30]);
        let bad = pileup_record_with_quals(1000, &seq, &[30, 30, 2, 30]);
        assert_eq!(check_allele_pileup(&good, &seq, &rv), Some(true));
        assert_eq!(check_allele_pileup(&bad, &seq, &rv), None);
    }

    #[test]
    fn a_pattern_whose_qualities_are_omitted_is_rejected() {
        let rv = ResolvedVariant {
            pos_0based: 1000,
            alt_bases: "ACGT".to_string(),
            ref_bases: "A".to_string(),
        };
        let seq = b"ACGT".to_vec();
        let omitted = pileup_record_with_quals(1000, &seq, &[QUAL_UNAVAILABLE; 4]);
        assert_eq!(check_allele_pileup(&omitted, &seq, &rv), None);
    }
```

`pileup_record_with_quals(pos, seq, quals)` builds a single aligned record at
`pos` with an all-`Match` CIGAR over `seq` and the given per-base qualities —
the same construction `spanning_record_quals` uses, without the two-site
layout. The `ResolvedVariant` field names above are verified against
`phase_disambiguate.rs:183-187` (`pos_0based: i64`, `alt_bases: String`,
`ref_bases: String`), and `check_allele_pileup(record, seq, rv) -> Option<bool>`
against `:456-460`.

- [ ] **Step 1c: Run the new tests.** Expected before Step 3:
`one_bad_base_inside_a_multi_base_pattern_rejects_the_whole_pattern` FAILS on
the `bad` assertion (left `Some(true)`, right `None`),
`a_pattern_whose_qualities_are_omitted_is_rejected` FAILS the same way, and
both `good` assertions pass. After Step 3, all pass.

- [ ] **Step 2: Run them.** Run: `cargo test --release --locked decide_phase 2>&1 | grep -E "^test |left|right"`
  Expected:
  - `base_quality_zero_reads_do_not_decide_phase` FAILS (left `Some("*1/*A")`, right `None`);
  - `the_same_reads_at_base_quality_30_do_decide_phase` PASSES.

- [ ] **Step 3: Implement.** Replace `read_bases_at_pos` with a span finder plus a thin wrapper. The CIGAR walk is unchanged.

In `depth_calling/mod.rs` — **not** in `phase_disambiguate.rs`. Task 6 gates
the same bases inside `depth_calling/haplotype.rs`, and `depth_calling`
currently imports nothing from `caller` (verified: zero `crate::caller`
references under `src/depth_calling/`). Defining these in `caller` would make
Task 6 introduce the first upward dependency in the crate. They belong in the
lower layer, where the bases are read:

```rust
/// Minimum base quality for an allele-discriminating base: the
/// `bcftools mpileup` default (`--min-BQ 13`). A pattern counts only if every
/// base it was matched on reaches it. Used by both the fragment-phasing path
/// (`caller::phase_disambiguate`) and the targeted-rule path
/// (`depth_calling::haplotype`), so it lives here, in the lower layer.
pub const PHASE_MIN_BASE_QUALITY: u8 = 13;

/// The value SAM/BAM stores when base qualities are **omitted**. It is not a
/// quality of 255: it means there is no quality information, so it can never
/// satisfy `PHASE_MIN_BASE_QUALITY` (review #5).
pub const QUAL_UNAVAILABLE: u8 = 0xFF;
```

Then in `phase_disambiguate.rs`, re-export so the documented public path and
every existing reference keep working:

```rust
pub use crate::depth_calling::{PHASE_MIN_BASE_QUALITY, QUAL_UNAVAILABLE};

/// Query offset and number of available bases (at most `n`) for a 0-based
/// genomic position, walking the CIGAR. `None` if the position is not on an
/// aligned base of the read.
fn query_span_at_pos(
    record: &rust_htslib::bam::Record,
    seq_len: usize,
    genomic_pos_0based: i64,
    n: usize,
) -> Option<(usize, usize)> {
    let ref_offset = genomic_pos_0based - record.pos();
    if ref_offset < 0 {
        return None;
    }
    let mut ref_pos = 0i64;
    let mut query_pos = 0usize;
    for &cigar_elem in record.cigar().iter() {
        use rust_htslib::bam::record::Cigar::*;
        match cigar_elem {
            Match(len) | Equal(len) | Diff(len) => {
                let len = len as i64;
                if ref_pos + len > ref_offset {
                    let qp = query_pos + (ref_offset - ref_pos) as usize;
                    return if qp < seq_len { Some((qp, n.min(seq_len - qp))) } else { None };
                }
                ref_pos += len;
                query_pos += len as usize;
            }
            Ins(len) => { query_pos += len as usize; }
            Del(len) => {
                let len = len as i64;
                if ref_pos + len > ref_offset { return None; }
                ref_pos += len;
            }
            SoftClip(len) => { query_pos += len as usize; }
            HardClip(_) | Pad(_) => {}
            RefSkip(len) => { ref_pos += len as i64; }
        }
    }
    None
}

/// Read up to `n` bases from a read starting at a 0-based genomic position.
fn read_bases_at_pos(
    record: &rust_htslib::bam::Record,
    seq: &[u8],
    genomic_pos_0based: i64,
    n: usize,
) -> Option<String> {
    let (qp, len) = query_span_at_pos(record, seq.len(), genomic_pos_0based, n)?;
    Some(seq[qp..qp + len].iter().map(|&b| b as char).collect())
}
```

Replace the body of `check_allele_pileup` (keep its doc comment, and add one line saying that bases below `PHASE_MIN_BASE_QUALITY` give `None`):

```rust
    let check_len = rv.alt_bases.len().max(rv.ref_bases.len());
    let (qp, len) = query_span_at_pos(record, seq.len(), rv.pos_0based, check_len)?;
    let bases_upper: String = seq[qp..qp + len].iter().map(|&b| b as char).collect::<String>().to_uppercase();
    let alt_upper = rv.alt_bases.to_uppercase();
    let ref_upper = rv.ref_bases.to_uppercase();
    let is_alt = match (
        bases_upper.starts_with(&alt_upper),
        bases_upper.starts_with(&ref_upper),
    ) {
        (true, true) => alt_upper.len() > ref_upper.len(),
        (true, false) => true,
        (false, true) => false,
        (false, false) => return None, // holds neither whole pattern
    };
    // Every base the winning pattern was matched on must be trustworthy.
    //
    // Unavailable quality is NOT high quality (review #5). SAM/BAM stores
    // "qualities omitted" as 0xFF, which would sail past a `q < 13` test, and
    // a quality array too short to cover the compared bases means those bases
    // have no quality at all. Both leave the phase unresolved rather than
    // deciding it. The earlier draft of this plan said 0xFF "passes"; that was
    // the same defect this task exists to remove.
    let compared = if is_alt { alt_upper.len() } else { ref_upper.len() };
    let quals = record.qual();
    if quals.len() < qp + compared {
        return None; // no quality for the bases compared
    }
    if quals[qp..qp + compared]
        .iter()
        .any(|&q| q == QUAL_UNAVAILABLE || q < PHASE_MIN_BASE_QUALITY)
    {
        return None;
    }
    Some(is_alt)
```

- [ ] **Step 4: Run all tests.** Expected: all ok, including the existing `check_allele_pileup_*` tests (their fixtures use quality 30) and the Step 1b pair.

- [ ] **Step 5: Commit** (`fix(phasing): bases below quality 13 are not phasing evidence`, with the review #4 reference and the G2 trailer).

- [ ] **Step 6: Regression gate (G4).** No change is allowed by this task. If HG00589 or NA18973 loses its fragment-phase resolution, STOP and report. Do not lower the threshold.

```bash
R=~/cyrius-research/review-fixes
: > $R/expected_task_5.tsv
: > $R/expected_cumulative_after_5.tsv   # 3c withdrawn: this plan changes no benchmark output
$R/gate.sh 5
# only after GATE PASS:
rm -rf $R/accepted && cp -r "$(cat $R/logs/task_5/results_dir.txt)" $R/accepted
```


---

### Task 6: The targeted cn2 phase rules get a third state (review #5)

**Files:**
- Modify: `cyrius-rs/src/caller/call_variants.rs:397-438` (`call_var42126938`), `:499-520` (`call_var42127803hap`)
- Modify: `cyrius-rs/src/depth_calling/haplotype.rs:31-270` (`get_bases_per_read`) — the quality gate for the targeted path, see Step 0 below
- Modify: `cyrius-rs/src/caller/match_star_allele.rs:755-769` (`apply_cn2_phase_rules`), `:821-822` (`match_star` parameters)
- Modify: `cyrius-rs/src/main.rs:763-790, 836-837` (types flow through; no logic change)
- Test: `cyrius-rs/src/depth_calling/haplotype.rs` (tests module, through a real BAM fixture), plus both caller modules' tests

**Interfaces:**
- Produces:
  - `call_var42126938(...) -> Result<(Vec<usize>, Vec<String>, Option<bool>), String>`
  - `call_var42127803hap(...) -> Result<Option<bool>, String>`
  - `match_star(..., var42126938_g_haplotype: Option<bool>, var42127803_diff_haplotype: Option<bool>, ...)`
  - `Some(true)` = evidence for the first diplotype of the rule, `Some(false)` = evidence for the second, `None` = unknown (both are kept).

**Rules** (fixed now):
- **g.42127803** (positions 42127803, 42127941; codes from `haplotype.rs`): trans (`*119/*2`) when reads show `"12"` ≥ 2 and `"21"` ≥ 2. Cis (`*1/*41`) when reads show `"11"` ≥ 2 and `"22"` ≥ 2. If both or neither hold, unknown.
- **g.42126938**: the positive rule is kept (`"12"` at positions 1, 2 seen in more than one read means `*1/*32`). Everything else is **unknown**. There is no `Some(false)`: no pattern is established that only `*27/*41` can produce, and reads that merely fail to show `"12"` may be covering sequence both candidates share, so counting them is not negative evidence (review #5, second pass). An earlier draft of this plan admitted `Some(false)` once 3 reads were "informative", where "informative" counted every pattern seen at those positions — the original defect with a threshold added. Task 7 records this as an open question with no plan.
- Consequence to accept deliberately: this rule can now only ever *narrow towards* `*1/*32`, never away from it. Ties it used to break the other way stay ties. No benchmark call is affected (measured: none of the 53 is one of the four diplotypes these rules choose between).

- [ ] **Step 0: Gate base quality on the targeted path too (review #3, second pass).**

Task 5 gates `phase_disambiguate.rs`. The targeted cn2 rules do not go through
it: they read bases via `get_haplotypes_from_bam` / `get_haplotypes_from_bam_single_region`
in `depth_calling/haplotype.rs`, which applies `passing_read` and a MAPQ floor
but **no base-quality check at all** (verified 2026-11-04: no `qual` reference
in either function or in their shared helper). The reviewer demonstrated the
consequence with four synthetic base-quality-zero reads: the extraction yields
`11 = 2`, `22 = 2`, and `phase_signal_42127803` then returns `Some(false)` —
a definitive cis call built entirely on untrustworthy bases.

Both public entry points funnel through one helper, `get_bases_per_read`, and
its only non-test callers are the three targeted cn2 rules
(`call_variants.rs:417`, `:472`, `:509`). So one gate in one place covers the
whole targeted path and nothing else. `call_var42127526_var42127556` gains the
gate as well as the two rules this task changes; that is intended.

Reuse Task 5's constants — do not define a second threshold. They live in
`depth_calling/mod.rs`, so add `use crate::depth_calling::{PHASE_MIN_BASE_QUALITY, QUAL_UNAVAILABLE};`
at the top of `haplotype.rs`. Do not import from `crate::caller::`: nothing in
`depth_calling` does, and this task must not be the first.

In `get_bases_per_read`, directly after the existing `if end_pos > seq.len() { continue; }`:

```rust
                    // Every extracted base must be trustworthy, by the same
                    // rule the fragment-phasing path uses (review #3/#5).
                    // A short or omitted quality array is "unknown", never
                    // "good": both drop the read at this site.
                    let quals = record.qual();
                    if quals.len() < end_pos
                        || quals[qpos..end_pos]
                            .iter()
                            .any(|&q| q == QUAL_UNAVAILABLE || q < PHASE_MIN_BASE_QUALITY)
                    {
                        continue;
                    }
```

- [ ] **Step 0b: Test it through a BAM, not through the count helper.** Testing
`phase_signal_42127803(c11, c12, c21, c22)` directly cannot catch this: the
defect is that bad bases reach it with respectable counts. The test must build
a BAM, run the extraction, and assert the counts.

```rust
    #[test]
    fn base_quality_zero_reads_are_not_haplotype_evidence() {
        // Four reads that would otherwise give 11=2, 22=2 — enough for
        // phase_signal_42127803 to answer Some(false) — are dropped because
        // their bases carry quality 0 (review #3).
        let (mut reader, base_db) = haplotype_fixture_with_quality(0);
        let haps = get_haplotypes_from_bam_single_region(&mut reader, &base_db, &[0, 1]).unwrap();
        let counts = extract_hap(&haps, &[0, 1]);
        assert!(
            counts.values().all(|v| v.iter().sum::<u32>() == 0),
            "quality-0 bases must not be counted, got {counts:?}"
        );
    }

    #[test]
    fn the_same_reads_at_quality_30_are_haplotype_evidence() {
        // Control: the test above must fail for quality, not for the fixture.
        let (mut reader, base_db) = haplotype_fixture_with_quality(30);
        let haps = get_haplotypes_from_bam_single_region(&mut reader, &base_db, &[0, 1]).unwrap();
        let counts = extract_hap(&haps, &[0, 1]);
        assert_eq!(counts.get("11").map(|v| v.iter().sum::<u32>()), Some(2));
        assert_eq!(counts.get("22").map(|v| v.iter().sum::<u32>()), Some(2));
    }

    #[test]
    fn reads_without_base_qualities_are_not_haplotype_evidence() {
        let (mut reader, base_db) = haplotype_fixture_with_quality(0xFF);
        let haps = get_haplotypes_from_bam_single_region(&mut reader, &base_db, &[0, 1]).unwrap();
        let counts = extract_hap(&haps, &[0, 1]);
        assert!(counts.values().all(|v| v.iter().sum::<u32>() == 0));
    }
```

`haplotype_fixture_with_quality(q)` writes a two-site BAM with four reads —
two showing the region-1 base at both sites (`"11"`) and two showing the
region-2 base at both (`"22"`) — every base at quality `q`, and returns an
indexed reader plus a matching `SnpLookup`. Build it from the existing
fixture helpers in this module's tests (`sed -n '470,520p' cyrius-rs/src/depth_calling/haplotype.rs`
shows the pattern the current tests use); the only new parameter is the
quality byte.

- [ ] **Step 0c: Run them.** Expected before the Step 0 change:
`base_quality_zero_reads_are_not_haplotype_evidence` and
`reads_without_base_qualities_are_not_haplotype_evidence` FAIL (counts 2 and 2,
not 0), `the_same_reads_at_quality_30_are_haplotype_evidence` PASSES. After the
change, all three pass.

- [ ] **Step 1: Write the failing tests**

In the `match_star_allele.rs` tests:
```rust
    #[test]
    fn unknown_phase_keeps_both_diplotypes() {
        assert_eq!(apply_cn2_phase_rules("*1/*41;*119/*2", None, None), "*1/*41;*119/*2");
        assert_eq!(apply_cn2_phase_rules("*1/*41;*119/*2", None, Some(true)), "*119/*2");
        assert_eq!(apply_cn2_phase_rules("*1/*41;*119/*2", None, Some(false)), "*1/*41");
        assert_eq!(apply_cn2_phase_rules("*1/*32;*27/*41", Some(true), None), "*1/*32");
        assert_eq!(apply_cn2_phase_rules("*1/*32;*27/*41", Some(false), None), "*27/*41");
        assert_eq!(apply_cn2_phase_rules("*1/*32;*27/*41", None, None), "*1/*32;*27/*41");
    }
```
In the `call_variants.rs` tests:
```rust
    #[test]
    fn phase_signal_42127803_needs_affirmative_reads() {
        assert_eq!(phase_signal_42127803(0, 5, 4, 0), Some(true));   // trans
        assert_eq!(phase_signal_42127803(6, 0, 0, 3), Some(false));  // cis
        assert_eq!(phase_signal_42127803(0, 5, 0, 0), None);         // one pattern only
        assert_eq!(phase_signal_42127803(0, 0, 0, 0), None);         // no reads
        assert_eq!(phase_signal_42127803(3, 3, 3, 3), None);         // both
    }

    #[test]
    fn phase_signal_42126938_never_says_no() {
        // Corrected 2026-11-04: an earlier draft of this plan gave this test
        // two arguments and expected `Some(false)` once three reads were
        // "informative". That contradicts the rule above and reinstates the
        // very defect review #5 is about -- absence of the G pattern is not
        // evidence of the opposite arrangement. The implementer caught the
        // contradiction and followed the rule, which was right.
        assert_eq!(phase_signal_42126938(2), Some(true));
        assert_eq!(phase_signal_42126938(5), Some(true));
        assert_eq!(phase_signal_42126938(1), None);
        assert_eq!(phase_signal_42126938(0), None);
    }
```

- [ ] **Step 2: Run them.** Expected: compile errors (`expected bool, found Option`, and `cannot find function 'phase_signal_42127803'`).

- [ ] **Step 3: Implement.**

In `call_variants.rs`, add above `call_var42126938`:
```rust
/// g.42126938 rule: `Some(true)` (G haplotype, *1/*32) when more than one
/// read shows it; otherwise **unknown**.
///
/// There is deliberately no `Some(false)`. Reads that fail to show the G
/// pattern do not demonstrate the opposite arrangement — they may cover
/// sequence that both candidates share — so a read count cannot stand in for
/// negative evidence, however large it is (review #5, second pass). An
/// earlier draft of this plan returned `Some(false)` once three reads were
/// "informative", where "informative" counted every pattern at those
/// positions; that is the original defect with a threshold bolted on.
///
/// Reinstating a `false` branch requires an affirmative pattern that only the
/// non-G arrangement can produce, measured against phase truth. That is not
/// in scope here and has no plan yet.
fn phase_signal_42126938(g_reads: u32) -> Option<bool> {
    if g_reads > 1 { Some(true) } else { None }
}

/// g.42127803 rule over read patterns at (42127803, 42127941), '1' = CYP2D6
/// base, '2' = CYP2D7 base: trans (*119/*2, `Some(true)`) needs "12" and "21"
/// in at least two reads each; cis (*1/*41, `Some(false)`) needs "11" and "22"
/// in at least two reads each; both or neither is unknown.
fn phase_signal_42127803(c11: u32, c12: u32, c21: u32, c22: u32) -> Option<bool> {
    let trans = c12 >= 2 && c21 >= 2;
    let cis = c11 >= 2 && c22 >= 2;
    match (trans, cis) {
        (true, false) => Some(true),
        (false, true) => Some(false),
        _ => None,
    }
}
```
In `call_var42126938`:
- set the return type's last element to `Option<bool>`;
- start with `let mut g_haplotype: Option<bool> = None;`;
- the early return becomes `return Ok((vec![], Vec::new(), None));`;
- replace the inner `if let Some(count_g) = g_hap_count.get("12") { ... }` block with:
```rust
                    let g_reads = g_hap_count.get("12").map_or(0, |v| v.iter().sum::<u32>());
                    g_haplotype = phase_signal_42126938(g_reads);
```
Replace the body of `call_var42127803hap`:
```rust
pub fn call_var42127803hap(
    reader: &mut bam::IndexedReader,
    cnvtag: &str,
    base_db: &SnpLookup,
) -> Result<Option<bool>, String> {
    if cnvtag != "cn2" {
        return Ok(None);
    }
    let n_positions = base_db.dsnp1.len();
    let positions: Vec<usize> = (0..n_positions).collect();
    let haplotype_per_read = get_haplotypes_from_bam_single_region(reader, base_db, &positions)?;
    let counts = extract_hap(&haplotype_per_read, &[0, 1]);
    let n = |k: &str| counts.get(k).map_or(0, |v| v.iter().sum::<u32>());
    Ok(phase_signal_42127803(n("11"), n("12"), n("21"), n("22")))
}
```
In `match_star_allele.rs`:
- change both `bool` parameters of `apply_cn2_phase_rules` and `match_star` to `Option<bool>`;
- in the rule loop, insert `let Some(signal) = signal else { continue };` as the first line of the `for` body.

In `main.rs`, the values flow through unchanged; fix any type annotation the compiler reports.

- [ ] **Step 4: Run all tests.** Expected: all ok. Fix existing `match_star(...)` test calls that pass `false`/`true` by passing `Some(false)`/`Some(true)`; the compiler lists them.

- [ ] **Step 5: Commit** (`fix(calling): targeted cn2 phase rules keep both diplotypes when phase is unknown`, with the review #5 reference and the G2 trailer).

- [ ] **Step 6: Regression gate (G4).** No change is allowed by this task. No benchmark call is one of these four diplotypes, and Step 0's quality gate touches only the three targeted cn2 rules, so a difference here is a bug.

```bash
R=~/cyrius-research/review-fixes
: > $R/expected_task_6.tsv
: > $R/expected_cumulative_after_6.tsv   # 3c withdrawn: this plan changes no benchmark output
$R/gate.sh 6
# only after GATE PASS:
rm -rf $R/accepted && cp -r "$(cat $R/logs/task_6/results_dir.txt)" $R/accepted
```


---

### Task 7: Documentation, final gate and reproduction

**Files:**
- Modify: `README.md` (Filter values table, JSON fields, "Changes from bcyrius")
- Modify: `docs/reviews/2026-11-04-project-review.md` (add a "Status" line under each P1 heading)

- [ ] **Step 1: README.**
  - **No Filter table row.** `Uncallable_target_sites` is a JSON field, not a filter value — 3c was withdrawn, so no filter string and no confidence cap exist. Adding that row would document behaviour the code does not have.
  - JSON fields: add `Uncallable_target_sites` — the list of target variants that got no dosage call. Say precisely what that means, because there is more than one route to it: no usable reads at the site, too few reads for any copy number to reach its posterior cutoff, or the strand-bias filter (`NOISY_VAR`) having zeroed the alt count. Add the caveat that a handful of sites also have a dedicated caller that can resolve them from elsewhere (`g.42130655-42130656insA` via the haplotype database), so a listed site means "the main dosage path could not call it", not "nothing is known about it" — the field errs toward over-reporting missing evidence, which is the safe direction for a warning. State plainly that it is **reported only**: it does not change the genotype, the filter or the confidence in this release, and that an unobserved site is still treated as reference by the matcher. That is the deferred fix.
  - JSON fields: add `Phase_pairs_truncated` (added by Task 8). Document all three states and be explicit that `null` is not `false`: `null` = the phase step did not run; `false` = it ran to completion, **including the legitimate case where it ran and had nothing to phase**, so an empty `Phase_pairs` with `false` is a real result and not a failure; `true` = the list was cut short, either by a BAM read error partway or by the file not opening at all. Say that every pair listed is a complete measurement even when the list is truncated — it is the list that is short, not the numbers.
  - "Changes from bcyrius": one line each for the five fixes that change behaviour (reviews #2, #3, #4, #5, #6), naming what changed. Review #1 changes no behaviour, so it belongs under the new field's description, not here. `Phase_pairs_truncated` is a new output field rather than a behaviour change, so it belongs with the JSON fields too.

- [ ] **Step 2: Review file.** Under each P1 heading add `**Status (fix/review-p1):** <status>; test <test name>.` Use these statuses exactly (per G9):
  - **#1:** `diagnosed and reported, fix deferred to its own plan` — this branch adds the `Uncallable_target_sites` field and the measurement; the matcher still treats an unobserved site as reference. Attach the measurement: the damaged NA19920 takes three target sites to zero reads and 31 sites lose depth, the genotype moves `*1/*4x2` → `*1/*10x2`, and confidence *rises* 0.87 → 0.95.
  - **#2, #3, #6:** `fixed in <short hash>`.
  - **#4:** `fixed in <short hash>` plus: "The main allele-counting path (`snp_count.rs`) is unchanged; see Out of scope." Also note that omitted qualities (`0xFF`) and a quality array too short to cover the compared bases both leave phase unresolved.
  - **#5:** `partly fixed in <short hash>: the unsupported negative branch was removed, not replaced` plus: "`phase_signal_42126938` can no longer return `Some(false)`; no affirmative non-G pattern is established. Open question, no plan."
  - Do not write "fixed" for #1 or a bare "fixed" for #5.

- [ ] **Step 3: Final gate.**
  - Run `gate.sh 7` with both `expected_task_7.tsv` and `expected_cumulative_after_7.tsv` empty. With 3c withdrawn, **no task in this plan may change any of the 53 outputs** — which makes the gate strictly stronger than before. The lint check is inside `gate.sh` and compares against `lint_baseline.txt` from Task 0; there is no separate clippy step and no hardcoded 68.
  - Re-run the NA19920 reproduction built in 3b Step 1, with the final binary, and record the output. Expected **unchanged from the baseline measurement**: genotype `*1/*10x2`, and `Uncallable_target_sites` non-empty listing at least `g.42128945C>T`, `g.42128936-42128937insGGGGCGAAA` and `g.42128936-42128937insGGGGCGAAAGGGGCGAAA`. The filter stays `PASS` and the confidence stays `0.95 HIGH`: that is the deferred fix, not a regression, and the review file must say so.

- [ ] **Step 4: Commit** (`docs: README and review status for the P1 fixes`, with the G2 trailer). Report to the user:
  - the per-sample comparison output;
  - the reproduction output;
  - the test totals;
  - the lint count.

---

### Task 8: A truncated phase-pair list says so (review follow-up)

**Added 2026-11-04 by user decision**, after the Task 2 review found the one
place where partial evidence still reaches the output. Task 2 made an aborted
phasing run stop deciding, which was right, but it still reports the variant
pairs it had already finished and the JSON array reads as complete. Every
number in it is a correct, finished measurement — the list is just shorter
than it would have been, and the only hint is that `Resolved_by` is absent.

The user chose to mark it rather than leave it to the docs. Run this **after
Task 6 and before Task 7**, so Task 7 documents the field.

**Files:**
- Modify: `cyrius-rs/src/caller/phase_disambiguate.rs` (`PhaseOutcome`, the
  abort return at `:687`, the normal return at `:750`)
- Modify: `cyrius-rs/src/types.rs` (`D6Call`, new field)
- Modify: `cyrius-rs/src/main.rs` (`phase_pairs_json` at `:184`, the
  `disambiguate` call site at `~:882`, the JSON record)
- Test: both modules' tests

**Interfaces:**
- Produces: `PhaseOutcome.pairs_truncated: bool`,
  `D6Call.phase_pairs_truncated: Option<bool>`, JSON key
  `Phase_pairs_truncated`.

**The three states, fixed now:**

| Situation | `Phase_pairs` | `Phase_pairs_truncated` |
|---|---|---|
| Phase step did not run | `null` | `null` |
| Ran to completion | array | `false` |
| Aborted on a read error | array of the pairs finished | `true` |

`Phase_pairs_truncated` is not one of the fields `compare_runs.py` compares,
so this task must still show 0 benchmark differences.

- [ ] **Step 1: Write the failing tests.**

In `main.rs` tests, extend the existing `phase_pairs_json` coverage:
```rust
    #[test]
    fn phase_pairs_truncated_is_null_when_the_phase_step_did_not_run() {
        assert_eq!(phase_pairs_truncated_json(None), serde_json::Value::Null);
    }

    #[test]
    fn phase_pairs_truncated_distinguishes_a_complete_list_from_an_aborted_one() {
        assert_eq!(phase_pairs_truncated_json(Some(false)), serde_json::json!(false));
        assert_eq!(phase_pairs_truncated_json(Some(true)), serde_json::json!(true));
    }
```

In `phase_disambiguate.rs` tests:
```rust
    #[test]
    fn a_completed_run_is_not_marked_truncated() {
        let records = (0..3).map(|i| spanning_record(&format!("q{i}"), true, true, 30, 0)).collect();
        let outcome = disambiguate_fixture("complete", records);
        assert_eq!(outcome.resolved.as_deref(), Some("*1/*A"));
        assert!(!outcome.pairs_truncated, "a run that finished must not be marked truncated");
    }
```

And the one that matters — reuse Task 2's corrupt-BAM construction so the
abort path is reached for real, not simulated:
```rust
    #[test]
    fn an_aborted_run_marks_its_pair_list_truncated() {
        // Same corrupt-BGZF fixture Task 2 used: the read fails partway, so
        // `disambiguate` stops. It still reports the pairs it finished, and
        // must say that the list is incomplete.
        let outcome = disambiguate_on_corrupt_fixture();
        assert_eq!(outcome.resolved, None, "an aborted run must not decide");
        assert!(outcome.pairs_truncated, "an aborted run must mark its list truncated");
    }
```

`disambiguate_on_corrupt_fixture()` builds the fixture the way
`check_phasing_reports_a_read_error_instead_of_partial_counts` already does in
this module (Task 2 added it) and calls `disambiguate` on it. Reuse that
helper rather than writing a second corruption routine; if it is not factored
out, factor it out and have both tests use it.

- [ ] **Step 2: Run them.**

Run: `cargo test --release --locked truncated 2>&1 | grep -E "^test |left|right|error\["`
Expected: the three new tests fail to compile (`no field 'pairs_truncated'`,
`cannot find function 'phase_pairs_truncated_json'`). That is the RED.

- [ ] **Step 3: Implement.**

In `phase_disambiguate.rs`, add to `PhaseOutcome`:
```rust
    /// True when the step stopped early (a read error) and `pairs` therefore
    /// holds only the pairs finished before it stopped. Each listed pair is
    /// still a complete measurement; the list is not complete. `false` on a
    /// run that evaluated every differentiating pair.
    pub pairs_truncated: bool,
```
`Default` derives `false`, so the existing `PhaseOutcome::default()` returns
stay correct: a step that never started has an empty list, not a truncated one.

At the abort return (`:687`) use `pairs_truncated: true`; at the normal return
(`:750`) use `pairs_truncated: false`.

In `types.rs`, add to `D6Call`:
```rust
    /// See `PhaseOutcome.pairs_truncated`. `None` when the phase step did not
    /// run.
    pub phase_pairs_truncated: Option<bool>,
```
and `None` in `D6Call::no_call`.

In `main.rs`, beside `phase_pairs_json`:
```rust
/// `Phase_pairs_truncated`: `null` when the phase step did not run, otherwise
/// whether `Phase_pairs` is cut short by a read error.
fn phase_pairs_truncated_json(truncated: Option<bool>) -> serde_json::Value {
    match truncated {
        None => serde_json::Value::Null,
        Some(t) => serde_json::json!(t),
    }
}
```
Set it from the outcome where `phase_pairs = Some(outcome.pairs);` is assigned
(capture `outcome.pairs_truncated` before moving `outcome.pairs`), and add
`"Phase_pairs_truncated": phase_pairs_truncated_json(cyp2d6_call.phase_pairs_truncated),`
immediately after the `"Phase_pairs"` entry.

- [ ] **Step 4: Run all tests.** Expected: all ok, output pristine.

- [ ] **Step 5: Commit** (`feat(output): Phase_pairs says when a read error cut it short`, with the G2 trailer).

- [ ] **Step 6: Regression gate (G4).** No change is allowed: the new key is not a compared field.

```bash
R=~/cyrius-research/review-fixes
: > $R/expected_task_8.tsv
: > $R/expected_cumulative_after_8.tsv
$R/gate.sh 8
# only after GATE PASS:
rm -rf $R/accepted && cp -r "$(cat $R/logs/task_8/results_dir.txt)" $R/accepted
```

---


## Out of scope (named so they aren't lost)

- Base quality in the main allele-counting path (`snp_count.rs`). It changes the core dosage counts inherited from Cyrius, so it needs its own measured experiment.
- Review P2 items 7–10: confidence calibration, resolver model, a strict benchmark mode in the repo, and the long-read reflex wording.
- Provenance (`Resolved_by`) for the targeted cn2 rules.
- Beta-binomial dosage, count-based HMM, and held-out donor validation. Each needs its own plan through the judgment gate.
- A pattern-based "not G" criterion for g.42126938.

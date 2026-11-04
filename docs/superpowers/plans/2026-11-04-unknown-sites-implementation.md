# Unknown Target Sites in Allele Matching — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A target site with no dosage call stops counting as "no variant". It matches any copy count; every compatible candidate is kept; a phenotype is reported only when all candidates agree.

**Architecture:** The matcher today looks up one exact key. When a sample has unknown sites, it instead scans the table for every key that agrees on all *known* sites, and skips the two frequency rules that narrow a candidate set by habit rather than evidence. The existing reporting (`clinical_equivalence`, the resolver, `More_than_one_possible_genotype`) then does the rest unchanged. A separate formatting bug that this makes reachable is fixed first.

**Tech Stack:** Rust (floor 1.85), rust-htslib 0.47.1, the 53-sample GeT-RM benchmark, the gate experiments in `docs/superpowers/plans/2026-11-04-unknown-sites-gate.md`.

## Global Constraints

- **No sample without an unknown site may change.** On all 53 full-depth benchmark samples, every output field except `Tool_commit` and `Flags` stays identical. (0 of 53 have an unknown site — measured.)
- **Unknown sites never narrow a call.** No frequency rule (the rare-allele rule; the `*68` "last candidate" pick) is applied to a masked candidate set.
- **Exact path untouched.** When a sample has no unknown site, every matcher function runs today's code, unchanged, by construction (an early branch, not an edited body).
- **Rust floor 1.85.** `cargo +1.85 check --all-targets --locked` must pass. No std API newer than 1.85.
- **Lint.** `cargo clippy --all-targets --locked` message count ≤ **22** (measured at `f215094`).
- **Tests.** `cargo test --release --locked` prints at least one `test result` line and every one shows `0 failed`.
- **Expectations before runs.** Expected output changes are written to the gate's expected files *before* the gate runs.
- **One commit per task**, message prefix `fix:` / `feat:` / `test:` / `docs:`, ending with the `Co-Authored-By` line.

## Decision needed before Task 4

**D1 — what happens when a sample has unknown sites and *no* table key fits even allowing for them.** The P1 plan agreed "no-call is the fallback". Measured since: in the 25 % thinned data, NA19908 is exactly this case, and today's fallback (fuzzy matching) returns `*1/*46;*43/*45;*12/*43` — three options **including the truth** — already labelled `More_than_one_possible_genotype`. A no-call would discard that.

- **D1-a (recommended):** keep today's fallback (retry, then fuzzy), but the result never reaches `PASS`: it gets the new filter value `Unknown_sites_unmatched`.
- **D1-b:** no-call, as first agreed.

Tasks 1–3 do not depend on D1. Task 4 is written for D1-a; the D1-b variant is given in Task 4 Step 6.

---

## File Structure

| File | Change |
|---|---|
| `cyrius-rs/src/caller/match_star_allele.rs` | Task 1 formatting fix; Task 2 masked lookup; Task 3 `*68` paths; Task 4 `match_star` takes the unknown sites |
| `cyrius-rs/src/types.rs` | Task 4: `StarCall` gains `unknown_sites_unmatched: bool` |
| `cyrius-rs/src/main.rs` | Task 4: pass the uncallable list; one filter value; two PASS-only checks |
| `README.md`, `docs/reviews/2026-11-04-project-review.md` | Task 6 |
| `~/cyrius-research/unknown-sites/gate.sh` | Task 0 (outside the repo) |

---

## Task 0: The gate for this branch

**Files:** Create `~/cyrius-research/unknown-sites/gate.sh`, `~/cyrius-research/unknown-sites/lint_baseline.txt`, `~/cyrius-research/unknown-sites/compare_all_fields.py`.

- [ ] **Step 1: copy the lint-branch gate and point it here**

```bash
R=~/cyrius-research/unknown-sites
sed -e 's#/home/parlar_ai/dev/cyrius-rs-lint#/home/parlar_ai/dev/cyrius-rs-unknown#g' \
    -e 's#R="\$HOME/cyrius-research/lint-cleanup"#R="$HOME/cyrius-research/unknown-sites"#' \
    ~/cyrius-research/lint-cleanup/gate.sh > $R/gate.sh && chmod +x $R/gate.sh
echo 22 > $R/lint_baseline.txt
grep -n "cyrius-rs-lint\|lint-cleanup" $R/gate.sh || echo "no stale paths"
```

Expected: `no stale paths`.

- [ ] **Step 2: add the all-fields comparison the gate lacks**

The gate's `compare_runs.py` checks 10 fields. This checks all of them, and refuses to pass when the fields it compares are absent (the vacuous pass that bit once in this branch).

```python
#!/usr/bin/env python3
"""compare_all_fields.py OLD_DIR NEW_DIR -> exit 1 on any difference.
Compares every field except Tool_commit and Flags, per sample. Fields are
nested under the sample name; refuses to pass if Genotype is absent anywhere."""
import json, glob, os, sys
old_dir, new_dir = sys.argv[1], sys.argv[2]
def load(d):
    out = {}
    for f in glob.glob(os.path.join(d, "*.json")):
        s = os.path.basename(f)[:-5]
        out[s] = json.load(open(f))[s]
    return out
a, b = load(old_dir), load(new_dir)
if set(a) != set(b) or not a:
    sys.exit(f"sample sets differ or empty: {len(a)} vs {len(b)}")
if any("Genotype" not in a[s] or "Genotype" not in b[s] for s in a):
    sys.exit("Genotype missing: refusing a vacuous comparison")
bad = [(s, k) for s in sorted(a) for k in sorted(set(a[s]) | set(b[s]))
       if k not in ("Tool_commit", "Flags") and a[s].get(k) != b[s].get(k)]
for s, k in bad:
    print(f"{s}\t{k}\t{json.dumps(a[s].get(k))[:80]}\t->\t{json.dumps(b[s].get(k))[:80]}")
print(f"{len(a)} samples, {len(bad)} differing field-values")
sys.exit(1 if bad else 0)
```

Append to `gate.sh`, before its final `echo "GATE PASS ..."`:

```bash
python3 "$R/compare_all_fields.py" "$REF/accepted" "$NEW" || fail "a field other than Tool_commit/Flags changed"
```

- [ ] **Step 3: prove the gate rejects, then passes**

```bash
R=~/cyrius-research/unknown-sites
: > $R/expected_task_0.tsv; : > $R/expected_cumulative_after_0.tsv
echo 21 > $R/lint_baseline.txt; bash $R/gate.sh 0 2>&1 | tail -1   # must FAIL: clippy 22 > baseline
echo 22 > $R/lint_baseline.txt; bash $R/gate.sh 0 2>&1 | tail -1   # must PASS
python3 $R/compare_all_fields.py ~/cyrius-research/review-fixes/accepted ~/cyrius-research/unknown-sites/thinned/out_p25; echo "exit=$?"   # must exit 1
```

Expected: `GATE FAIL (task 0): clippy 22 > baseline`, then `GATE PASS (task 0)`, then a list of differences and `exit=1`.

---

## Task 1: Format every candidate of a multi-candidate call

**The bug.** `get_final_call_clean` turns a single candidate multiset into a diplotype (`*10_*14_*36` → `*14/*36+*10`), and does so for several candidates only in group `cn2`. For any other group with more than one candidate it joins the raw multisets (`match_star_allele.rs:326-334`), which `match_phenotype` cannot read, so the phenotype becomes Indeterminate. No full-depth benchmark sample reaches it (0 of 53 genotypes contain `_`). It is reached in real data: HG02373 at 25 % depth prints `*10_*14_*36;*10_*36_*65;*10_*114_*36`, although the first candidate is its truth. Masking would make it common.

**Files:** Modify `cyrius-rs/src/caller/match_star_allele.rs` (`get_final_call_clean`, ~L326-334; tests module).

- [ ] **Step 1: write the failing test**

Add to the `tests` module:

```rust
    #[test]
    fn several_candidates_outside_cn2_are_each_formatted() {
        // HG02373 at 25 % depth (spacer CN 3). At full depth its single
        // candidate *10_*14_*36 is formatted *14/*36+*10 (measured).
        let call = vec![
            "*10_*14_*36".to_string(),
            "*10_*36_*65".to_string(),
            "*10_*114_*36".to_string(),
        ];
        // Each candidate formatted alone (measured): *14/*36+*10, *65/*36+*10,
        // *114/*36+*10. Today all three together give the raw
        // "*10_*114_*36;*10_*14_*36;*10_*36_*65" (measured).
        assert_eq!(
            get_final_call_clean(&call, "exon9hyb", Some(3)).as_deref(),
            Some("*114/*36+*10;*14/*36+*10;*65/*36+*10")
        );
    }

    #[test]
    fn the_hand_written_two_candidate_case_is_unchanged() {
        let call = vec!["*10_*10_*4.013".to_string(), "*10_*36_*4".to_string()];
        assert_eq!(
            get_final_call_clean(&call, "exon9hyb", Some(3)).as_deref(),
            Some("*4/*36+*10;*4.013/*10+*10")
        );
    }
```

- [ ] **Step 2: run; the first fails, the second passes**

Run: `cd cyrius-rs && cargo test --release --locked several_candidates_outside_cn2 the_hand_written_two -- --nocapture`
Expected: `several_candidates_outside_cn2_are_each_formatted ... FAILED`, left `"*10_*114_*36;*10_*14_*36;*10_*36_*65"`; the other `ok`.

- [ ] **Step 3: implement**

Replace

```rust
    if final_call.is_empty() || final_call.len() > 1 {
        if final_call == ["*10_*10_*4.013", "*10_*36_*4"] {
```

with

```rust
    if final_call.is_empty() {
        return Some(String::new());
    }
    if final_call.len() > 1 {
        if final_call == ["*10_*10_*4.013", "*10_*36_*4"] {
```

and replace that block's closing `return Some(final_call.join(";"));` with:

```rust
        // Format each candidate the way a single answer is formatted. Joining
        // the raw multisets ("*1_*4_*4") gave strings match_phenotype cannot
        // read. A candidate that cannot be formatted keeps its raw form, which
        // reads as unknown effect: the safe direction.
        let mut out: Vec<String> = Vec::new();
        for one in &final_call {
            let formatted = get_final_call_clean(std::slice::from_ref(one), cnvcall, spacer_cn)
                .unwrap_or_else(|| one.clone());
            for d in formatted.split(';') {
                if !out.iter().any(|x| x == d) {
                    out.push(d.to_string());
                }
            }
        }
        return Some(out.join(";"));
```

- [ ] **Step 4: run both tests, then the suite** — both pass; suite 0 failed.

- [ ] **Step 5: gate** — `: > ~/cyrius-research/unknown-sites/expected_task_1.tsv; : > ~/cyrius-research/unknown-sites/expected_cumulative_after_1.tsv; bash ~/cyrius-research/unknown-sites/gate.sh 1`. Expected `GATE PASS`: the full-depth benchmark does not reach this branch.

- [ ] **Step 6: commit** — `fix(matching): format each candidate of a multi-candidate call outside cn2`.

---

## Task 2: The masked lookup

**Files:** Modify `cyrius-rs/src/caller/match_star_allele.rs` (new helpers above `get_star`; `get_star`; `get_star_cn1`; `retry_add_one_copy`; existing test call sites).

**Interfaces:**
- Produces: `fn get_star(var_observed: &[String], unknown: &[String], dic: &HashMap<String, Vec<String>>) -> RawStarCall` and `fn get_star_cn1(var_observed: &[String], unknown: &[String], dic: &HashMap<String, String>) -> RawStarCall`. An empty `unknown` runs today's body unchanged.

- [ ] **Step 1: write the failing tests**

```rust
    fn obs(v: &[&str]) -> Vec<String> { v.iter().map(|s| s.to_string()).collect() }

    #[test]
    fn an_unknown_site_never_acts_as_a_reference_call() {
        // Today "B unseen" means "no B": only *x. Masked, *y (which has B) stays.
        let d = dic(&[("g.A", &["*x"]), ("g.A_g.B", &["*y"])]);
        let call = get_star(&obs(&["g.A"]), &obs(&["g.B"]), &d);
        assert_eq!(call.candidate, obs(&["*x", "*y"]));
        assert_eq!(get_star(&obs(&["g.A"]), &[], &d).candidate, obs(&["*x"]));
    }

    #[test]
    fn unknown_sites_are_masked_together() {
        let d = dic(&[("g.A", &["*x"]), ("g.A_g.B_g.C", &["*y"]), ("g.A_g.C", &["*z"])]);
        let call = get_star(&obs(&["g.A"]), &obs(&["g.B", "g.C"]), &d);
        assert_eq!(call.candidate, obs(&["*x", "*y", "*z"]));
    }

    #[test]
    fn a_known_site_still_excludes_candidates() {
        // Masking B must not admit *w, which lacks the observed A.
        let d = dic(&[("g.A", &["*x"]), ("g.B", &["*w"])]);
        assert_eq!(get_star(&obs(&["g.A"]), &obs(&["g.B"]), &d).candidate, obs(&["*x"]));
    }

    #[test]
    fn a_frequency_rule_does_not_narrow_a_masked_set() {
        // *34 is in RARE_ALLELES. Unmasked, the rare rule picks the common one.
        let d = dic(&[("g.A", &["*1_*2"]), ("g.A_g.B", &["*1_*34"])]);
        let masked = get_star(&obs(&["g.A"]), &obs(&["g.B"]), &d);
        assert_eq!(masked.star_call, obs(&["*1_*2", "*1_*34"]));
        assert_eq!(masked.call_info.as_deref(), Some("more_than_one_match"));
    }

    #[test]
    fn an_irrelevant_unknown_site_gives_the_exact_answer_on_the_real_table() {
        // Scan path vs exact path, every 100th key of the real cn2 table.
        let c = crate::caller::construct_star_table::get_hap_table(
            crate::data::STAR_TABLE, &[crate::caller::construct_star_table::TableKind::Dhap2]);
        let t = &c.dhap2;
        let mut keys: Vec<&String> = t.keys().collect();
        keys.sort();
        let none = obs(&["g.not_a_site"]);
        let mut checked = 0;
        for k in keys.iter().step_by(100) {
            let vars: Vec<String> = if *k == "NA" { vec![] } else { k.split('_').map(String::from).collect() };
            let mut exact = t[*k].clone();
            exact.sort();
            exact.dedup(); // the masked path returns a set
            assert_eq!(get_star(&vars, &none, t).candidate, exact, "key {k}");
            checked += 1;
        }
        assert!(checked >= 200, "only {checked} keys checked");
    }

    #[test]
    fn damaged_na19920_keeps_its_true_genotype_among_the_candidates() {
        // Measured in gate B1: 5 unknown sites, 9 candidate multisets.
        let c = crate::caller::construct_star_table::get_hap_table(
            crate::data::STAR_TABLE, &[crate::caller::construct_star_table::TableKind::Dhap3Pair]);
        let observed = obs(&["g.42126611C>G", "g.42126611C>G", "g.42129130C>G", "g.42129130C>G",
                             "g.42129754G>A", "g.42129754G>A", "g.42130692G>A", "g.42130692G>A"]);
        let unknown = obs(&["g.42128936-42128937insGGGGCGAAAGGGGCGAAA", "g.42128936-42128937insGGGGCGAAA",
                            "g.42128932C>T", "g.42128933G>A", "g.42128945C>T"]);
        let call = get_star(&observed, &unknown, &c.dhap3pair);
        assert_eq!(call.star_call, obs(&["*10_*10_*145", "*10_*10_*4", "*145_*147_*147", "*145_*4_*4",
            "*147_*147_*4", "*1_*10_*10", "*1_*147_*147", "*1_*4_*4", "*4_*4_*4"]));
        // Unmasked, today's answer: the single wrong *1_*10_*10.
        assert_eq!(get_star(&observed, &[], &c.dhap3pair).star_call, obs(&["*1_*10_*10"]));
    }
```

Update the existing tests that call `get_star(&x, &d)` to `get_star(&x, &[], &d)`.

- [ ] **Step 2: run; expect a compile error** (`get_star` takes 2 arguments) — that is the red for a signature change.

- [ ] **Step 3: implement**

Above `get_star`:

```rust
/// Copy counts of every variant that is not an unknown site.
fn known_site_counts<'a>(
    vars: impl Iterator<Item = &'a str>,
    unknown: &[String],
) -> std::collections::BTreeMap<&'a str, usize> {
    let mut counts = std::collections::BTreeMap::new();
    for v in vars {
        if !unknown.iter().any(|u| u == v) {
            *counts.entry(v).or_insert(0) += 1;
        }
    }
    counts
}

/// Keys whose variants agree with `var_observed` at every site not in
/// `unknown`; an unknown site may carry any number of copies. A full scan:
/// measured 0.80 s on the largest table (1,551,523 keys), 0.14 s or less on
/// the others. Only called when `unknown` is non-empty.
fn compatible_keys<'k>(
    var_observed: &[String],
    unknown: &[String],
    keys: impl Iterator<Item = &'k String>,
) -> Vec<&'k String> {
    let target = known_site_counts(var_observed.iter().map(String::as_str), unknown);
    keys.filter(|k| {
        let vars = if k.as_str() == "NA" { "" } else { k.as_str() };
        known_site_counts(vars.split('_').filter(|v| !v.is_empty()), unknown) == target
    })
    .collect()
}

/// get_star when some sites are unknown: every compatible candidate is kept
/// and no frequency rule narrows the set.
fn get_star_masked(var_observed: &[String], unknown: &[String], dic: &HashMap<String, Vec<String>>) -> RawStarCall {
    let mut raw: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    for k in compatible_keys(var_observed, unknown, dic.keys()) {
        raw.extend(dic[k].iter().cloned());
    }
    raw_to_call(raw.into_iter().collect())
}

fn get_star_cn1_masked(var_observed: &[String], unknown: &[String], dic: &HashMap<String, String>) -> RawStarCall {
    let raw: std::collections::BTreeSet<String> =
        compatible_keys(var_observed, unknown, dic.keys()).into_iter().map(|k| dic[k].clone()).collect();
    raw_to_call(raw.into_iter().collect())
}

/// Labels a masked candidate list the way get_star labels an exact one,
/// minus the rare-allele rule.
fn raw_to_call(raw: Vec<String>) -> RawStarCall {
    if raw.is_empty() {
        return RawStarCall { call_info: Some("no_match".to_string()), candidate: vec![], star_call: vec![] };
    }
    let processed = convert_to_main_allele(&raw);
    let info = if raw.len() == 1 {
        "unique_match"
    } else if processed.len() == 1 {
        "unique_star"
    } else {
        "more_than_one_match"
    };
    RawStarCall { call_info: Some(info.to_string()), candidate: raw, star_call: processed }
}
```

Change the first lines of `get_star` to

```rust
fn get_star(var_observed: &[String], unknown: &[String], dic: &HashMap<String, Vec<String>>) -> RawStarCall {
    if !unknown.is_empty() {
        return get_star_masked(var_observed, unknown, dic);
    }
```

leaving the rest of its body unchanged; likewise `get_star_cn1` with `get_star_cn1_masked`. In `retry_add_one_copy` change `get_star(&var_new, dic)` to `get_star(&var_new, &[], dic)` (it runs only in the fallback, which keeps today's behaviour).

- [ ] **Step 4: run the six tests and the suite.** All pass. If `damaged_na19920_keeps_its_true_genotype_among_the_candidates` lists different multisets, **stop**: the measured B1 set is the expectation and must not be edited to match.

- [ ] **Step 5: `cargo +1.85 check --all-targets --locked`** — passes.

- [ ] **Step 6: gate 2** (expected files empty). `match_star` still passes no unknown sites, so nothing can change. **Commit** `feat(matching): masked lookup — an unknown site matches any copy count`.

---

## Task 3: The `*68` paths

**Files:** `call_star68`, `call_star68_cn1` in `match_star_allele.rs`.

**Interfaces:** `fn call_star68(var_observed: &mut Vec<String>, unknown: &[String], cnvcall: &str, dic: &HashMap<String, Vec<String>>) -> RawStarCall`; `call_star68_cn1` likewise with `HashMap<String, String>`.

- [ ] **Step 1: failing test**

```rust
    #[test]
    fn the_star68_paths_do_not_pick_by_habit_when_sites_are_unknown() {
        // Two different answers after removing one g.42130692G>A copy. Unmasked,
        // the "last one is most likely" rule picks one; masked, both are kept.
        let d = dic(&[("g.42130692G>A_g.X", &["*68_*4"]), ("g.X", &["*68_*2"])]);
        let mut v = obs(&["g.42130692G>A", "g.X"]);
        let picked = call_star68(&mut v.clone(), &[], "star68", &d);
        assert_eq!(picked.star_call.len(), 1);
        let kept = call_star68(&mut v, &obs(&["g.Y"]), "star68", &d);
        assert_eq!(kept.star_call.len(), 2, "{:?}", kept.star_call);
        assert_eq!(kept.call_info.as_deref(), Some("more_than_one_match"));
    }
```

- [ ] **Step 2: run; compile error** (signature).

- [ ] **Step 3: implement.** Add the `unknown` parameter; pass it to every `get_star` / `get_star_cn1` call inside; and immediately after the existing `if unique.len() == 1 { ... }` block in **both** functions insert:

```rust
    if !unknown.is_empty() {
        // Unknown sites: report every candidate rather than picking by rule.
        let mut all: Vec<String> = unique.into_iter().cloned().collect();
        all.sort();
        return RawStarCall {
            call_info: Some("more_than_one_match".to_string()),
            candidate: hap_list.clone(),
            star_call: all,
        };
    }
```

Update both callers in `match_star` to pass `&[]` for now.

- [ ] **Step 4: test + suite + `cargo +1.85 check`. Step 5: gate 3 (empty expectations). Step 6: commit** `feat(matching): the *68 paths keep every candidate when sites are unknown`.

---

## Task 4: Wire the unknown sites through `match_star`

**Files:** `match_star_allele.rs` (`match_star`), `types.rs` (`StarCall`), `main.rs`.

**Interfaces:** `pub fn match_star(var_observed: &mut Vec<String>, unknown_sites: &[String], cnvcall: &str, ... )` (new second parameter); `StarCall.unknown_sites_unmatched: bool`.

- [ ] **Step 1: add the field** to `StarCall` in `types.rs`:

```rust
    /// True when the sample had unknown target sites and no table key fit even
    /// with them masked, so the call came from the retry or fuzzy fallback,
    /// which treats an unknown site as absent. Such a call never reaches PASS.
    pub unknown_sites_unmatched: bool,
```

`cargo build` then lists every `StarCall { .. }` literal (8 in `match_star_allele.rs`). Add `unknown_sites_unmatched: false` to each; Step 3 sets the fallback ones.

- [ ] **Step 2: failing test**

```rust
    #[test]
    fn a_listed_site_already_in_the_call_is_not_masked() {
        // g.B is listed uncallable but a dedicated caller (or update_variants)
        // put it in the call, so it is known; g.C stays unknown.
        let still = still_unknown(&obs(&["g.B", "g.C"]), &obs(&["g.A", "g.B"]));
        assert_eq!(still, obs(&["g.C"]));
        assert!(still_unknown(&[], &obs(&["g.A"])).is_empty());
    }
```

(End-to-end behaviour of `match_star` is tested by Task 5's runs.)

- [ ] **Step 3: implement in `match_star`.** After `update_variants(var_observed, cnvcall, exon9);`:

```rust
    let unknown = still_unknown(unknown_sites, var_observed);
```

with, next to `get_var_list`:

```rust
/// Sites the dosage step could not call, minus any that a dedicated caller or
/// `update_variants` has since placed in the call.
fn still_unknown(unknown_sites: &[String], var_observed: &[String]) -> Vec<String> {
    unknown_sites.iter().filter(|u| !var_observed.contains(u)).cloned().collect()
}
```

Pass `&unknown` to `get_star`, `get_star_cn1`, `call_star68`, `call_star68_cn1`. On the main path, after `let mut matchtag = get_star(var_observed, &unknown, dic);`:

```rust
        let unknown_sites_unmatched =
            !unknown.is_empty() && matchtag.call_info.as_deref() == Some("no_match");
```

and set `unknown_sites_unmatched` (instead of `false`) in the `StarCall` returned by the retry branch and in the final `StarCall` of that path. Do the same on the single-copy and `*68` paths, using their own `no_match` result.

- [ ] **Step 4: `main.rs`.** Pass `&uncallable_target_sites` as the new second argument of `match_star`. In the filter chain insert, before `} else if is_fuzzy {`:

```rust
            } else if star_called.unknown_sites_unmatched {
                genotype_filter = Some("Unknown_sites_unmatched".to_string());
```

In the two PASS-only checks (`spacer_cn_check` ~L976, `consistency_check` ~L1015) change `filter == "PASS" || filter == "Fuzzy_match"` to also accept `"Unknown_sites_unmatched"`, so those checks still see the calls they saw before.

- [ ] **Step 5: test + suite + `cargo +1.85 check` + gate 4 (empty expectations).** Commit `feat(matching): an uncallable target site matches any copy count`.

- [ ] **Step 6 (only if D1-b was chosen instead):** in Step 3, when `unknown_sites_unmatched` is true, return `StarCall { call_info: Some("no_match".into()), raw_call: Some(vec![]), clean_call: None, variants_called: Some(var_observed.join(" ")), unknown_sites_unmatched: true }` instead of running retry/fuzzy, and drop the Step 4 filter value.

---

## Task 5: Prove it on the inputs that have holes

No code change. Rules locked here, before running:

| Input | Expected | Source of the expectation |
|---|---|---|
| 53 full-depth samples | 0 differing field-values (gate, `compare_all_fields.py`) | 0 of 53 have an unknown site (measured) |
| Damaged NA19920 | The 9 diplotypes below as a **set**; filter `More_than_one_possible_genotype`; `Clinical_equivalence` `unknown_effect`; `Long_read_recommended` true; confidence 0.80 MEDIUM | Measured with the B2 prototype |
| NA21105, 25 % | `*3/*111`, `PASS`, `Resolved_by` `resolver` | Measured with the B2 prototype |
| NA20289, 25 % | `*6/*11`, filter `LowQ_coverage_MAD` | Measured with the B2 prototype |
| HG02373, 25 % | No `_` in the genotype; `*14/*36+*10` among its options | Task 1 test |
| HG03259, 25 % (single-copy path, not in the prototype) | **PREDICTED (not run)** — record what happens | — |

The 9 NA19920 diplotypes: `*1/*10x2`, `*1/*147x2`, `*1/*4x2`, `*145/*147x2`, `*4x2/*145`, `*4/*4x2`, `*10x2/*145`, `*4/*10x2`, `*4/*147x2`.

**Stop rules:** any full-depth difference; any thinned or damaged run that is confidently wrong with this code and was not before; the NA19920 set differing from the 9 above.

- [ ] **Step 1:** run the damaged NA19920 and the 25 % and 10 % sets with the new binary into `~/cyrius-research/unknown-sites/final/`.
- [ ] **Step 2:** `mask_audit --real` on each output directory: report confidently-wrong-with-holes and confidently-wrong-without-holes, against the counts recorded in the gate record.
- [ ] **Step 3:** diff the 25 % and 10 % outputs against `~/cyrius-research/unknown-sites/thinned/out_p25` and `out_p10`; every changed sample must be in the table above or explained in the result commit.
- [ ] **Step 4:** append the results to the gate record; commit `result:`.

---

## Task 6: Documentation

- [ ] **README.md:** `Uncallable_target_sites` is no longer "reported only": its sites are matched as unknown. Describe the new filter value `Unknown_sites_unmatched`, and that an ambiguous result caused by unknown sites appears as `More_than_one_possible_genotype` with the alternatives listed. Remove the sentence that names the NA19920 phenotype flip as a deferred defect and replace it with the measured after-state from Task 5.
- [ ] **docs/reviews/2026-11-04-project-review.md:** review #1 status → fixed, with the commit hashes and the Task 5 measurements.
- [ ] Commit `docs:`.

---

## Self-review notes

- Every expected value above is either measured (gate record, B1/B2 outputs, the scan timing, HG02373's full-depth formatting) or labelled PREDICTED.
- The invariant "no unknown site, no change" holds by construction (early return into the masked functions only when `unknown` is non-empty) and is checked by the gate on all fields.
- Not covered by this plan: the per-site two-value mask (no measured case needs it: NA19920's low-read sites are 0 in every candidate), and making the resolver itself aware of unknown sites (K4 measured that it does not collapse NA19920).

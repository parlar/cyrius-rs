//! Gate B1 for `docs/superpowers/plans/2026-11-04-unknown-sites-gate.md`.
//!
//! Throwaway experiment. Changes nothing in `src/`. It answers two questions
//! with the real tables and the real reporting functions:
//!
//! * K1: with the damaged NA19920's five unknown sites masked, is the truth
//!   `*1/*4x2` among the candidates?
//! * K3: in a single-site dropout simulation over the main-path benchmark
//!   samples, does masking ever ESTABLISH a wrong phenotype?
//!
//! Usage:
//!   cargo run --release --example mask_audit -- <accepted_dir> <truth.tsv> <damaged.json>
//!   cargo run --release --example mask_audit -- --real <out_dir> <truth.tsv>
//!     (addendum: real coverage holes in thinned BAMs; rules S1/S2 and K3)

use cyrius_rs::caller::construct_star_table::{combination_table_for, get_hap_table, TableKind};
use cyrius_rs::caller::match_star_allele::get_final_call_clean;
use cyrius_rs::data::{HAPLOTYPE_FUNC, STAR_TABLE, TARGET_VARIANT, TARGET_VARIANT_HOMO};
use cyrius_rs::phenotype::{clinical_equivalence, load_haplotype_functionality, match_phenotype, sort_diplotypes};
use cyrius_rs::types::{HaplotypeFunctionality, StarCombinations};
use std::collections::{BTreeMap, BTreeSet, HashMap};

type Table = HashMap<String, Vec<String>>;

// ---- copied from match_star_allele.rs (private there); kept byte-for-byte in meaning ----
const RARE_ALLELES: &[&str] = &["*34", "*39", "*4.009", "*139"];
fn main_allele_name(star: &str) -> &str {
    if star == "*4.013" { star } else { star.split('.').next().unwrap_or(star) }
}
fn convert_to_main_allele(list: &[String]) -> Vec<String> {
    let mut set = BTreeSet::new();
    for stars in list {
        let mut v: Vec<&str> = stars.split('_').map(main_allele_name).collect();
        v.sort();
        set.insert(v.join("_"));
    }
    set.into_iter().collect()
}
fn key_of(vars: &[String]) -> String {
    let mut v = vars.to_vec();
    v.sort();
    if v.is_empty() { "NA".to_string() } else { v.join("_") }
}
// ------------------------------------------------------------------------------------------

/// get_star's post-processing of a raw candidate list. `rare_rule` = today's
/// behaviour; the masked path runs it with `false` (gate record, Gate A #4).
fn post_process(raw: &[String], rare_rule: bool) -> Vec<String> {
    let processed = convert_to_main_allele(raw);
    if raw.len() <= 1 || processed.len() == 1 || !rare_rule {
        return processed;
    }
    let non_rare: Vec<String> = raw
        .iter()
        .filter(|h| !RARE_ALLELES.iter().any(|r| h.contains(r)))
        .cloned()
        .collect();
    if non_rare.is_empty() { processed } else { convert_to_main_allele(&non_rare) }
}

fn table_of(c: &StarCombinations, k: TableKind) -> &Table {
    match k {
        TableKind::Dhap2 => &c.dhap2,
        TableKind::Dhap3 => &c.dhap3,
        TableKind::Dhap3Pair => &c.dhap3pair,
        TableKind::Dhap4Pair => &c.dhap4pair,
        TableKind::Dhap5Pair => &c.dhap5pair,
        TableKind::Dhap6Pair => &c.dhap6pair,
        TableKind::Exon9X2 => &c.dhap_exon9_x2,
        TableKind::Exon9X3 => &c.dhap_exon9_x3,
        TableKind::Exon9X4 => &c.dhap_exon9_x4,
        TableKind::DupExon9 => &c.dhap_dup_exon9,
    }
}

/// Largest count of any one variant in any key of the table. MEASURED per
/// table, so the enumeration below cannot miss a key by assuming a bound.
fn max_variant_count(t: &Table) -> usize {
    let mut best = 0;
    for k in t.keys() {
        if k == "NA" { continue; }
        let mut c: HashMap<&str, usize> = HashMap::new();
        for v in k.split('_') { *c.entry(v).or_default() += 1; }
        best = best.max(c.values().copied().max().unwrap_or(0));
    }
    best
}

/// Masked lookup by enumeration: every assignment of 0..=max copies to the
/// unknown sites, each looked up exactly. Returns (counts, key) per hit.
fn masked_enum(observed: &[String], unknown: &[String], max: usize, t: &Table) -> Vec<(Vec<usize>, String)> {
    let mut hits = Vec::new();
    let mut counts = vec![0usize; unknown.len()];
    loop {
        let mut vars = observed.to_vec();
        for (u, &n) in unknown.iter().zip(&counts) {
            for _ in 0..n { vars.push(u.clone()); }
        }
        let k = key_of(&vars);
        if t.contains_key(&k) { hits.push((counts.clone(), k)); }
        // odometer
        let mut i = 0;
        loop {
            if i == counts.len() { return hits; }
            counts[i] += 1;
            if counts[i] <= max { break; }
            counts[i] = 0;
            i += 1;
        }
    }
}

/// Masked lookup by full scan, as an independent cross-check of `masked_enum`.
fn masked_scan(observed: &[String], unknown: &[String], t: &Table) -> BTreeSet<String> {
    let mut obs: BTreeMap<&str, usize> = BTreeMap::new();
    for v in observed { *obs.entry(v.as_str()).or_default() += 1; }
    let unk: BTreeSet<&str> = unknown.iter().map(|s| s.as_str()).collect();
    let mut out = BTreeSet::new();
    for k in t.keys() {
        let mut kc: BTreeMap<&str, usize> = BTreeMap::new();
        if k != "NA" { for v in k.split('_') { *kc.entry(v).or_default() += 1; } }
        let a: BTreeMap<_, _> = kc.iter().filter(|(v, _)| !unk.contains(*v)).collect();
        let b: BTreeMap<_, _> = obs.iter().filter(|(v, _)| !unk.contains(*v)).collect();
        if a == b { out.insert(k.clone()); }
    }
    out
}

#[derive(Debug, Clone, PartialEq)]
enum Verdict { Established(String, String), Unresolved, Empty, NoCall }

/// Format each candidate multiset ON ITS OWN, the way the pipeline formats a
/// single answer. `get_final_call_clean` given several multisets in a non-cn2
/// group just joins the raw lists (match_star_allele.rs:326-334), which
/// `match_phenotype` cannot read -- a pipeline gap no benchmark sample reaches.
fn format_each(processed: &[String], cnv: &str, spacer: Option<u32>) -> Option<String> {
    let mut out: Vec<String> = Vec::new();
    for m in processed {
        let one = get_final_call_clean(std::slice::from_ref(m), cnv, spacer)?;
        for d in one.split(';') {
            if !out.iter().any(|x| x == d) { out.push(d.to_string()); }
        }
    }
    Some(out.join(";"))
}

fn verdict(processed: &[String], cnv: &str, spacer: Option<u32>, func: &HashMap<String, HaplotypeFunctionality>) -> (Verdict, Option<String>) {
    if processed.is_empty() { return (Verdict::Empty, None); }
    let call = match format_each(processed, cnv, spacer) { Some(c) => c, None => return (Verdict::NoCall, None) };
    let preds = match_phenotype(Some(&call), func);
    let (eq, _) = clinical_equivalence(&preds);
    if eq == "same_effect" {
        (Verdict::Established(preds[0].total_activity.clone(), preds[0].predicted_phenotype.clone()), Some(call))
    } else {
        (Verdict::Unresolved, Some(call))
    }
}

fn truth_effect(truth: &str, func: &HashMap<String, HaplotypeFunctionality>) -> Option<(String, String)> {
    // The benchmark's own normalisation of `*parent(*child)` -> `*child`.
    let re_paren = |s: &str| -> String {
        let mut out = String::new();
        let mut rest = s;
        while let Some(i) = rest.find('(') {
            let head = &rest[..i];
            let j = rest[i..].find(')').map(|j| i + j).unwrap();
            let child = &rest[i + 1..j];
            let start = head.rfind('*').unwrap();
            out.push_str(&head[..start]);
            out.push_str(child);
            rest = &rest[j + 1..];
        }
        out.push_str(rest);
        out
    };
    let t = re_paren(truth);
    let p = match_phenotype(Some(&t), func);
    if p.len() == 1 && p[0].total_activity != "n/a" {
        Some((p[0].total_activity.clone(), p[0].predicted_phenotype.clone()))
    } else {
        None
    }
}

fn target_names() -> Vec<String> {
    let mut v = Vec::new();
    for content in [TARGET_VARIANT, TARGET_VARIANT_HOMO] {
        for line in content.lines() {
            if line.starts_with('#') || line.trim().is_empty() { continue; }
            v.push(line.split('\t').nth(6).unwrap().to_string());
        }
    }
    v
}

fn is_main_path(cnv: &str) -> bool {
    !matches!(cnv, "star5_star5" | "star13_star13" | "star13intron1_star13intron1" | "star5_star5_star68")
        && !matches!(cnv, "star5" | "star13" | "star5_star68" | "star13_star68")
        && !cnv.contains("star68")
}

/// Truth diplotype -> sorted main-allele multiset, e.g. "*1/*4x2" -> "*1_*4_*4".
fn truth_multiset(t: &str) -> String {
    let mut v: Vec<String> = Vec::new();
    for hap in t.split('/') {
        for part in hap.split('+') {
            let (name, n) = match part.split_once('x') {
                Some((a, b)) if b.chars().all(|c| c.is_ascii_digit()) && !b.is_empty() => (a, b.parse::<usize>().unwrap()),
                _ => (part, 1),
            };
            for _ in 0..n { v.push(main_allele_name(name).to_string()); }
        }
    }
    v.sort();
    v.join("_")
}

fn norm_truth(truth: &str) -> String {
    let mut out = String::new();
    let mut rest = truth;
    while let Some(i) = rest.find('(') {
        let head = &rest[..i];
        let j = rest[i..].find(')').map(|j| i + j).unwrap();
        let start = head.rfind('*').unwrap();
        out.push_str(&head[..start]);
        out.push_str(&rest[i + 1..j]);
        rest = &rest[j + 1..];
    }
    out.push_str(rest);
    out
}

/// (activity, phenotype) of a ';'-joined call if all alternatives agree and
/// the activity is known; None otherwise ("not established").
fn established(call: &str, func: &HashMap<String, HaplotypeFunctionality>) -> Option<(String, String)> {
    let preds = match_phenotype(Some(call), func);
    let (eq, _) = clinical_equivalence(&preds);
    (eq == "same_effect" && preds[0].total_activity != "n/a")
        .then(|| (preds[0].total_activity.clone(), preds[0].predicted_phenotype.clone()))
}

fn real_mode(out_dir: &str, truth_path: &str) {
    let func = load_haplotype_functionality(HAPLOTYPE_FUNC);
    let mut truth: HashMap<String, String> = HashMap::new();
    for line in std::fs::read_to_string(truth_path).unwrap().lines() {
        if line.starts_with('#') || line.trim().is_empty() { continue; }
        let f: Vec<&str> = line.split('\t').collect();
        truth.insert(f[0].to_string(), f[1].to_string());
    }
    let mut paths: Vec<_> = std::fs::read_dir(out_dir).unwrap().filter_map(|e| { let p = e.ok()?.path(); (p.extension()? == "json").then_some(p) }).collect();
    paths.sort();
    let (mut runs, mut called, mut with_holes, mut s2, mut s2_pass, mut wrong_no_holes) = (0, 0, 0, 0, 0, 0);
    let (mut masked_runs, mut k3_wrong, mut truth_out) = (0, 0, 0);
    let mut combos_cache: HashMap<String, StarCombinations> = HashMap::new();
    for p in paths {
        let s = p.file_stem().unwrap().to_string_lossy().to_string();
        let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&p).unwrap()).unwrap();
        let d = &v[&s];
        runs += 1;
        let Some(call) = d["Genotype"].as_str() else { continue };
        called += 1;
        let tr = norm_truth(&truth[&s]);
        let tp = match_phenotype(Some(&tr), &func);
        let t_eff = (tp[0].total_activity.clone(), tp[0].predicted_phenotype.clone());
        let c_eff = established(call, &func);
        let holes: Vec<String> = d["Uncallable_target_sites"].as_array().map(|a| a.iter().map(|x| x.as_str().unwrap().to_string()).collect()).unwrap_or_default();
        let conf_wrong = matches!(&c_eff, Some(e) if *e != t_eff);
        let filter = d["Filter"].as_str().unwrap_or("None");
        if !holes.is_empty() { with_holes += 1; }
        if conf_wrong && holes.is_empty() { wrong_no_holes += 1; }
        let tag = if conf_wrong && !holes.is_empty() { s2 += 1; if filter == "PASS" { s2_pass += 1; } "S2 CONFIDENTLY WRONG WITH HOLES" } else if conf_wrong { "confidently wrong, no holes" } else { "" };
        println!("{s:8} call={call:<40} filter={filter:<32} est={c_eff:?} truth={} {t_eff:?} holes={} {tag}", truth[&s], holes.len());
        if holes.is_empty() { continue; }
        // Masking on the real holes (main path only).
        let cnv = d["CNV_group"].as_str().unwrap_or("");
        if !(is_main_path(cnv) && combination_table_for(cnv).is_some()) { println!("    masking: not on the main path ({cnv}), skipped"); continue; }
        let kind = combination_table_for(cnv).unwrap();
        let combos = combos_cache.entry(format!("{kind:?}")).or_insert_with(|| get_hap_table(STAR_TABLE, &[kind]));
        let t = table_of(combos, kind);
        let observed: Vec<String> = d["Variants_called"].as_array().map(|a| a.iter().map(|x| x.as_str().unwrap().to_string()).collect()).unwrap_or_default();
        let unknown: Vec<String> = holes.iter().filter(|h| !observed.contains(h)).cloned().collect();
        if unknown.len() != holes.len() { println!("    {} listed site(s) were resolved by a dedicated caller and are not masked", holes.len() - unknown.len()); }
        let hits = masked_enum(&observed, &unknown, max_variant_count(t), t);
        let mut raw: BTreeSet<String> = BTreeSet::new();
        for (_, k) in &hits { raw.extend(t[k].iter().cloned()); }
        let raw: Vec<String> = raw.into_iter().collect();
        let processed = post_process(&raw, false);
        let (mv, mcall) = verdict(&processed, cnv, d["Spacer_CN"].as_u64().map(|x| x as u32), &func);
        masked_runs += 1;
        let tin = processed.iter().any(|m| *m == truth_multiset(&tr));
        if !tin { truth_out += 1; }
        let mwrong = matches!(&mv, Verdict::Established(a, ph) if (a.clone(), ph.clone()) != t_eff);
        if mwrong { k3_wrong += 1; }
        println!("    masked {} site(s): {} keys -> {} candidates -> {mv:?}{} | truth in set: {tin}{}",
            unknown.len(), hits.len(), processed.len(), if mwrong { "  <-- K3 VIOLATION" } else { "" },
            mcall.map(|c| format!(" | call {}", if c.len() > 90 { format!("{}...", &c[..90]) } else { c })).unwrap_or_default());
    }
    println!("\nruns {runs} | gave a genotype {called} | listed >=1 uncallable site {with_holes}");
    println!("confidently wrong WITH holes (S2): {s2} (of which PASS {s2_pass}) | confidently wrong with NO holes: {wrong_no_holes}");
    println!("masking on real holes: {masked_runs} runs | ESTABLISHED-WRONG {k3_wrong} | truth fell out of the masked set {truth_out}");
    println!("rule: {}", if s2 == 0 { "S1 (no confidently wrong run with holes)" } else { "S2 (defect occurs in this data)" });
    println!("K3 on real holes: {}", if k3_wrong == 0 { "PASS" } else { "FAIL" });
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args[1] == "--real" { return real_mode(&args[2], &args[3]); }
    let (acc_dir, truth_path, dmg_path) = (&args[1], &args[2], &args[3]);

    let func = load_haplotype_functionality(HAPLOTYPE_FUNC);
    let targets = target_names();
    println!("target sites: {}", targets.len());

    let mut truth: HashMap<String, String> = HashMap::new();
    for line in std::fs::read_to_string(truth_path).unwrap().lines() {
        if line.starts_with('#') || line.trim().is_empty() { continue; }
        let f: Vec<&str> = line.split('\t').collect();
        truth.insert(f[0].to_string(), f[1].to_string());
    }

    // Load each sample's outputs.
    let mut samples: Vec<(String, serde_json::Value)> = Vec::new();
    let mut names: Vec<_> = std::fs::read_dir(acc_dir).unwrap().filter_map(|e| {
        let p = e.ok()?.path();
        (p.extension()? == "json").then_some(p)
    }).collect();
    names.sort();
    for p in names {
        let s = p.file_stem().unwrap().to_string_lossy().to_string();
        let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&p).unwrap()).unwrap();
        samples.push((s.clone(), v[&s].clone()));
    }
    let dmg_v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(dmg_path).unwrap()).unwrap();
    let dmg = dmg_v.as_object().unwrap().values().next().unwrap().clone();

    // Build only the tables needed.
    let mut kinds: BTreeSet<String> = BTreeSet::new();
    let mut wanted: Vec<TableKind> = Vec::new();
    for (_, d) in samples.iter().map(|(a, b)| (a, b)).chain(std::iter::once((&"dmg".to_string(), &dmg))) {
        let cnv = d["CNV_group"].as_str().unwrap_or("");
        if !is_main_path(cnv) { continue; }
        if let Some(k) = combination_table_for(cnv) {
            if kinds.insert(format!("{k:?}")) { wanted.push(k); }
        }
    }
    let combos = get_hap_table(STAR_TABLE, &wanted);
    let mut maxc: HashMap<String, usize> = HashMap::new();
    for k in &wanted {
        let t = table_of(&combos, *k);
        let m = max_variant_count(t);
        println!("table {k:?}: {} keys, max copies of one variant in a key = {m}", t.len());
        maxc.insert(format!("{k:?}"), m);
    }

    let vars_of = |d: &serde_json::Value| -> Vec<String> {
        d["Variants_called"].as_array().map(|a| a.iter().map(|x| x.as_str().unwrap().to_string()).collect()).unwrap_or_default()
    };
    let spacer_of = |d: &serde_json::Value| d["Spacer_CN"].as_u64().map(|x| x as u32);

    // ================= K1: the damaged NA19920 =================
    println!("\n================ K1: damaged NA19920 ================");
    let cnv = dmg["CNV_group"].as_str().unwrap();
    let kind = combination_table_for(cnv).unwrap();
    let t = table_of(&combos, kind);
    let observed = vars_of(&dmg);
    let unknown: Vec<String> = dmg["Uncallable_target_sites"].as_array().unwrap().iter().map(|x| x.as_str().unwrap().to_string()).collect();
    assert!(unknown.iter().all(|u| !observed.contains(u)), "an unknown site is also in the observed list");
    println!("cnv {cnv}, table {kind:?}, observed {observed:?}");
    println!("unknown sites ({}): {unknown:?}", unknown.len());
    let hits = masked_enum(&observed, &unknown, maxc[&format!("{kind:?}")], t);
    let scan = masked_scan(&observed, &unknown, t);
    let enum_keys: BTreeSet<String> = hits.iter().map(|(_, k)| k.clone()).collect();
    println!("compatible keys: {} by enumeration, {} by full scan, identical: {}", enum_keys.len(), scan.len(), enum_keys == scan);
    let raw_counts = dmg["Variant_raw_count"].clone();
    println!("raw reads (alt,ref) at unknown sites: {}", unknown.iter().map(|u| format!("{u}={}", raw_counts[u].as_str().unwrap_or("?"))).collect::<Vec<_>>().join("  "));
    let mut raw: BTreeSet<String> = BTreeSet::new();
    for (counts, k) in &hits {
        let entries = &t[k];
        let mains = convert_to_main_allele(entries);
        println!("  counts at unknown sites {counts:?} -> {} raw entries -> main {mains:?}", entries.len());
        raw.extend(entries.iter().cloned());
    }
    let raw: Vec<String> = raw.into_iter().collect();
    let processed = post_process(&raw, false);
    let (v, call) = verdict(&processed, cnv, spacer_of(&dmg), &func);
    println!("masked call: {}", call.clone().unwrap_or_default());
    for dip in call.clone().unwrap_or_default().split(';') {
        let p = match_phenotype(Some(dip), &func);
        println!("    {dip:<22} activity {:<5} {}", p[0].total_activity, p[0].predicted_phenotype);
    }
    println!("masked verdict: {v:?}");
    let k1 = processed.iter().any(|m| m == "*1_*4_*4");
    println!("K1 (truth *1/*4x2 i.e. main multiset *1_*4_*4 among candidates): {}", if k1 { "PASS" } else { "FAIL" });
    let with_rare = post_process(&raw, true);
    println!("rare-allele rule would change the candidate list here: {}", with_rare != processed);

    // ================= control + K3: dropout simulation =================
    println!("\n================ control: today's lookup reproduces today's call ================");
    let mut main_path = Vec::new();
    let mut excluded = Vec::new();
    for (s, d) in &samples {
        let cnv = d["CNV_group"].as_str().unwrap_or("");
        if is_main_path(cnv) && combination_table_for(cnv).is_some() { main_path.push((s, d)); } else { excluded.push(format!("{s}({cnv})")); }
    }
    println!("main-path samples: {}; excluded {}: {}", main_path.len(), excluded.len(), excluded.join(" "));
    // Control (a): the exact lookup finds the same RAW candidates the pipeline
    // reported. Independent of formatting and of every later step.
    // Control (b): the formatted call matches the reported Genotype, ignoring
    // order, on samples where nothing after the matcher changed the call.
    let alt_set = |g: &str| -> BTreeSet<String> { sort_diplotypes(g).split(';').map(|x| x.to_string()).collect() };
    let (mut raw_bad, mut fmt_bad, mut fmt_n, mut fuzzy_n) = (0, 0, 0, 0);
    for (s, d) in &main_path {
        let cnv = d["CNV_group"].as_str().unwrap();
        let t = table_of(&combos, combination_table_for(cnv).unwrap());
        let obs = vars_of(d);
        let mine: BTreeSet<String> = t.get(&key_of(&obs)).cloned().unwrap_or_default().into_iter().collect();
        let theirs: BTreeSet<String> = d["Raw_star_allele"].as_array().map(|a| a.iter().map(|x| x.as_str().unwrap().to_string()).collect()).unwrap_or_default();
        let fuzzy = d["Call_info"].as_str().map_or(false, |c| c.starts_with("fuzzy"));
        if fuzzy { fuzzy_n += 1; }
        if mine != theirs && !fuzzy {
            raw_bad += 1;
            println!("  RAW differs: {s} cnv={cnv} mine={mine:?} pipeline={theirs:?} call_info={}", d["Call_info"]);
        }
        if fuzzy || !d["Resolved_by"].is_null() { continue; }
        fmt_n += 1;
        let call = format_each(&post_process(&mine.iter().cloned().collect::<Vec<_>>(), true), cnv, spacer_of(d));
        let json_call = d["Genotype"].as_str().unwrap_or("");
        if call.as_deref().map(alt_set) != Some(alt_set(json_call)) {
            fmt_bad += 1;
            println!("  FORMAT differs: {s} cnv={cnv} json={json_call:?} mine={call:?}");
        }
    }
    println!("control (a) raw candidates: {raw_bad} differ of {} ({fuzzy_n} fuzzy-matched samples skipped: their key is not in the table by definition)", main_path.len());
    println!("control (b) formatted call: {fmt_bad} differ of {fmt_n} (samples with Resolved_by set, or fuzzy, skipped)");

    println!("\n================ K3: single-site dropout ================");
    #[derive(Default)]
    struct Tally { cases: usize, est_ok: usize, est_wrong: usize, unresolved: usize, empty: usize, nocall: usize, truth_na: usize, grew: usize, rare_would_change: usize, grew_ok: usize, grew_wrong: usize, grew_unres: usize, rare_wrong: usize }
    let (mut a, mut b) = (Tally::default(), Tally::default());
    // Same tallies restricted to samples whose UNMASKED call is established, so
    // that alleles with no function row do not inflate "unresolved".
    let (mut ea, mut eb) = (Tally::default(), Tally::default());
    let mut est_samples = 0;
    let (mut today_ok, mut today_wrong, mut today_nomatch, mut today_unres) = (0, 0, 0, 0);
    let mut wrong_cases: Vec<String> = Vec::new();
    for (s, d) in &main_path {
        let cnv = d["CNV_group"].as_str().unwrap();
        let kind = combination_table_for(cnv).unwrap();
        let t = table_of(&combos, kind);
        let m = maxc[&format!("{kind:?}")];
        let obs = vars_of(d);
        let spacer = spacer_of(d);
        let te = truth.get(*s).and_then(|tr| truth_effect(tr, &func));
        let base_keys: BTreeSet<String> = t.get(&key_of(&obs)).map(|_| key_of(&obs)).into_iter().collect();
        let base_raw = t.get(&key_of(&obs)).cloned().unwrap_or_default();
        let (bv, bcall) = verdict(&post_process(&base_raw, false), cnv, spacer, &func);
        let base_est = matches!(bv, Verdict::Established(..));
        if base_est { est_samples += 1; } else { println!("  baseline (no site masked) not established: {s} cnv={cnv} call={bcall:?} -> {bv:?}"); }
        for site in &targets {
            let present = obs.iter().filter(|v| *v == site).count();
            let reduced: Vec<String> = obs.iter().filter(|v| *v != site).cloned().collect();
            let hits = masked_enum(&reduced, std::slice::from_ref(site), m, t);
            let keys: BTreeSet<String> = hits.iter().map(|(_, k)| k.clone()).collect();
            let mut raw: BTreeSet<String> = BTreeSet::new();
            for k in &keys { raw.extend(t[k].iter().cloned()); }
            let raw: Vec<String> = raw.into_iter().collect();
            let processed = post_process(&raw, false);
            let (v, call) = verdict(&processed, cnv, spacer, &func);
            for tl in [if present > 0 { &mut a } else { &mut b }] { tl.cases += 1; }
            if base_est { let e = if present > 0 { &mut ea } else { &mut eb }; e.cases += 1;
                match (&v, &te) {
                    (Verdict::Established(act, ph), Some((ta, tp))) if act == ta && ph == tp => e.est_ok += 1,
                    (Verdict::Established(..), Some(_)) => e.est_wrong += 1,
                    _ => e.unresolved += 1,
                }
            }
            let tl = if present > 0 { &mut a } else { &mut b };
            if keys != base_keys { tl.grew += 1; }
            let grew = keys != base_keys;
            if !raw.is_empty() && post_process(&raw, true) != processed {
                tl.rare_would_change += 1;
                // Counterfactual: WITH the rare-allele rule, would this be a confident wrong answer?
                if let (Verdict::Established(ra, rp), _) = verdict(&post_process(&raw, true), cnv, spacer, &func) {
                    if let Some((ta, tp)) = &te { if ra != *ta || rp != *tp {
                        tl.rare_wrong += 1;
                        println!("  RARE RULE WOULD BE WRONG: {s} dropped={site} est=({ra},{rp}) truth=({ta},{tp})");
                    } }
                }
            }
            if grew {
                match (&v, &te) {
                    (Verdict::Established(act, ph), Some((ta, tp))) if act == ta && ph == tp => tl.grew_ok += 1,
                    (Verdict::Established(..), Some(_)) => tl.grew_wrong += 1,
                    _ => tl.grew_unres += 1,
                }
            }
            match (&v, &te) {
                (Verdict::Empty, _) => tl.empty += 1,
                (Verdict::NoCall, _) => tl.nocall += 1,
                (Verdict::Unresolved, _) => tl.unresolved += 1,
                (Verdict::Established(..), None) => tl.truth_na += 1,
                (Verdict::Established(act, ph), Some((ta, tp))) => {
                    if act == ta && ph == tp { tl.est_ok += 1; } else {
                        tl.est_wrong += 1;
                        wrong_cases.push(format!("{s} cnv={cnv} dropped={site} (present {present}) call={call:?} est=({act},{ph}) truth={:?} ({ta},{tp})", truth.get(*s)));
                    }
                }
            }
            // What today's exact lookup does when the dropped variant was present.
            if present > 0 {
                match t.get(&key_of(&reduced)) {
                    None => today_nomatch += 1,
                    Some(r) => match verdict(&post_process(r, true), cnv, spacer, &func) {
                        (Verdict::Established(act, ph), _) => match &te {
                            Some((ta, tp)) if act == *ta && ph == *tp => today_ok += 1,
                            Some(_) => today_wrong += 1,
                            None => {}
                        },
                        _ => today_unres += 1,
                    },
                }
            }
        }
    }
    let show = |name: &str, t: &Tally| {
        println!("{name}: cases {} | established-correct {} | ESTABLISHED-WRONG {} | unresolved {} | empty {} | get_final_call_clean None {} | truth effect n/a {} | rare rule would have changed {} (and been WRONG {})",
            t.cases, t.est_ok, t.est_wrong, t.unresolved, t.empty, t.nocall, t.truth_na, t.rare_would_change, t.rare_wrong);
        println!("    of which the candidate set GREW: {} -> established-correct {} | ESTABLISHED-WRONG {} | unresolved {}", t.grew, t.grew_ok, t.grew_wrong, t.grew_unres);
    };
    show("A (dropped variant was present)", &a);
    show("B (dropped variant was absent) ", &b);
    let tot = Tally { cases: a.cases + b.cases, est_ok: a.est_ok + b.est_ok, est_wrong: a.est_wrong + b.est_wrong, unresolved: a.unresolved + b.unresolved, empty: a.empty + b.empty, nocall: a.nocall + b.nocall, truth_na: a.truth_na + b.truth_na, grew: a.grew + b.grew, rare_would_change: a.rare_would_change + b.rare_would_change,
        grew_ok: a.grew_ok + b.grew_ok, grew_wrong: a.grew_wrong + b.grew_wrong, grew_unres: a.grew_unres + b.grew_unres, rare_wrong: a.rare_wrong + b.rare_wrong };
    show("TOTAL                          ", &tot);
    println!("restricted to the {est_samples} samples whose unmasked call is established:");
    println!("    A: cases {} | established-correct {} | ESTABLISHED-WRONG {} | unresolved {}", ea.cases, ea.est_ok, ea.est_wrong, ea.unresolved);
    println!("    B: cases {} | established-correct {} | ESTABLISHED-WRONG {} | unresolved {}", eb.cases, eb.est_ok, eb.est_wrong, eb.unresolved);
    println!("today's exact lookup on the A cases: right {today_ok} | WRONG {today_wrong} | unresolved {today_unres} | no exact match (would go to retry/fuzzy) {today_nomatch}");
    println!("K3 (zero established-wrong): {}", if tot.est_wrong == 0 { "PASS" } else { "FAIL" });
    for w in &wrong_cases { println!("  WRONG: {w}"); }
}

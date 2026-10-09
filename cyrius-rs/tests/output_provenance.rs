//! Provenance in the output (condition C3): `--version` and every per-sample
//! JSON record name the tool version, the git commit it was built from, and
//! the effective value of every feature flag, so a result file says which
//! build and which configuration produced it.
//!
//! The JSON test runs the real binary on a synthetic, read-free BAM: the run
//! ends in a depth no-call, but the JSON record is still written in full, so
//! no benchmark BAM is needed.

use rust_htslib::bam;
use std::path::{Path, PathBuf};
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_cyrius-rs");

/// The commit `build.rs` stamped into this build. `option_env!` rather than
/// `env!` so a missing stamp is a test failure, not a compile error.
fn build_commit() -> &'static str {
    option_env!("CYRIUS_GIT_COMMIT").expect("build.rs must set CYRIUS_GIT_COMMIT")
}

/// A short hash, optionally `-dirty`, or `unknown` when git was unavailable.
fn looks_like_a_commit(s: &str) -> bool {
    let hash = s.strip_suffix("-dirty").unwrap_or(s);
    s == "unknown" || (hash.len() == 7 && hash.bytes().all(|b| b.is_ascii_hexdigit()))
}

/// Write an indexed BAM with a chr22 header and no reads.
fn write_empty_bam(dir: &Path) -> PathBuf {
    let mut header = bam::header::Header::new();
    let mut sq = bam::header::HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr22");
    sq.push_tag(b"LN", 50_818_468);
    header.push_record(&sq);
    let path = dir.join("empty.bam");
    {
        bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
    }
    bam::index::build(&path, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
    path
}

/// Run the binary on a read-free BAM with `flags`; return the sample's JSON
/// record and the TSV text.
fn run_on_empty_bam(tag: &str, flags: &[&str]) -> (serde_json::Value, String) {
    let dir = std::env::temp_dir().join(format!("cyrius_provenance_{}_{}", tag, std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let bam = write_empty_bam(&dir);
    let out = dir.join("out");
    let status = Command::new(BIN)
        .args(["-i", bam.to_str().unwrap(), "-o", out.to_str().unwrap(), "--id", "s1"])
        .args(flags)
        .stderr(std::process::Stdio::null())
        .status()
        .expect("spawn cyrius-rs");
    let json = std::fs::read_to_string(out.join("s1.json"));
    let tsv = std::fs::read_to_string(out.join("s1.tsv"));
    let _ = std::fs::remove_dir_all(&dir);

    assert!(status.success(), "cyrius-rs exited with {status}");
    let mut json: serde_json::Value = serde_json::from_str(&json.unwrap()).unwrap();
    (json["s1"].take(), tsv.unwrap())
}

/// Every `FeatureFlags` field, in declaration order.
const FLAG_NAMES: [&str; 12] = [
    "strand_bias_all",
    "fuzzy_match",
    "phase_disambiguate",
    "phase_readpair",
    "changepoint_hybrid",
    "het_check",
    "spacer_cn_check",
    "consistency_check",
    "hmm_cnv",
    "clip_evidence",
    "d7_depth",
    "kmer_validation",
];

fn flags_of(record: &serde_json::Value) -> Vec<(String, bool)> {
    record["Flags"]
        .as_object()
        .unwrap_or_else(|| panic!("no Flags object in {record}"))
        .iter()
        .map(|(k, v)| (k.clone(), v.as_bool().unwrap_or_else(|| panic!("{k} is not a bool: {v}"))))
        .collect()
}

#[test]
fn version_names_the_tool_version_and_the_commit() {
    let out = Command::new(BIN).arg("--version").output().expect("spawn cyrius-rs");
    assert!(out.status.success());
    let stdout = String::from_utf8(out.stdout).unwrap();
    assert_eq!(stdout, format!("BCyrius 1.0.5 ({})\n", build_commit()));
    assert!(looks_like_a_commit(build_commit()), "{:?}", build_commit());
}

#[test]
fn json_records_tool_version_commit_and_every_effective_flag() {
    let (record, tsv) = run_on_empty_bam("noflags", &[]);

    assert_eq!(record["Tool_version"], "1.0.5");
    assert_eq!(record["Tool_commit"], build_commit());
    // Every flag, with its effective value: fuzzy matching is the default.
    let expected: Vec<(String, bool)> = FLAG_NAMES
        .iter()
        .map(|&n| (n.to_string(), n == "fuzzy_match"))
        .collect();
    assert_eq!(flags_of(&record), expected);

    // The TSV stays as it was.
    assert!(
        !tsv.contains("Tool_") && !tsv.contains("Flags"),
        "provenance must not be added to the TSV: {tsv}"
    );
}

#[test]
fn json_flags_follow_the_command_line() {
    let (record, _) = run_on_empty_bam("someflags", &["--no-fuzzy-match", "--het-check"]);
    let expected: Vec<(String, bool)> = FLAG_NAMES
        .iter()
        .map(|&n| (n.to_string(), n == "het_check"))
        .collect();
    assert_eq!(flags_of(&record), expected);
}

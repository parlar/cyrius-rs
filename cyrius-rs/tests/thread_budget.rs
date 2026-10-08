//! End-to-end thread budget for `--threads 1`.
//!
//! Review finding #10 was a resource bug, so the honest red for it is a
//! measurement, not a compile error: run the real binary at `-t 1` and watch
//! how many OS threads it actually has. Before the global pool was sized from
//! the flag, `--diplotype-caller -t 1` ran rayon on a pool of one worker per
//! core (72 on the development node); after, it runs on one.
//!
//! Linux-only (it reads `/proc/<pid>/status`) and skipped elsewhere, skipped
//! when the benchmark BAM is absent, and skipped on machines too small for the
//! bound to discriminate.

use std::path::Path;
use std::process::{Command, Stdio};
use std::time::Duration;

/// Generous ceiling on concurrent threads at `-t 1`.
///
/// The run's real peak is set by `construct_star_table::get_hap_table`, which
/// spawns 10 scoped threads unconditionally (not gated on `--threads`), plus
/// the main thread: measured at exactly 11 across 10 repeats on a 72-CPU node.
/// With the pool unsized, the same run peaks at 75. An unbounded rayon pool
/// adds one worker per core, so on any machine with more than `BUDGET` cores
/// this bound separates the two states by a wide margin in both directions.
const BUDGET: usize = 32;

#[test]
fn diplotype_caller_at_one_thread_stays_within_the_thread_budget() {
    if !Path::new("/proc/self/status").exists() {
        eprintln!("skipping: needs /proc/<pid>/status (Linux)");
        return;
    }
    let available = std::thread::available_parallelism().map_or(1, |n| n.get());
    if available <= BUDGET {
        eprintln!(
            "skipping: {} CPU(s) available, so an all-cores rayon pool would \
             also fit inside the budget of {}",
            available, BUDGET
        );
        return;
    }

    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let bam = root.join("benchmark/NA18545_pgx_regions.bam");
    if !bam.exists() {
        eprintln!("skipping: benchmark BAM not present");
        return;
    }

    let out_dir = std::env::temp_dir().join(format!("cyrius_thread_budget_{}", std::process::id()));
    std::fs::create_dir_all(&out_dir).expect("create temp out dir");

    let mut child = Command::new(env!("CARGO_BIN_EXE_cyrius-rs"))
        .args([
            "-i",
            bam.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "-t",
            "1",
            "--diplotype-caller",
        ])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .expect("spawn cyrius-rs");

    let status_path = format!("/proc/{}/status", child.id());
    let mut peak = 0usize;
    let status = loop {
        // Sample before reaping, so the last sample is of a live process.
        if let Ok(status) = std::fs::read_to_string(&status_path) {
            if let Some(n) = status
                .lines()
                .find_map(|line| line.strip_prefix("Threads:"))
                .and_then(|value| value.trim().parse::<usize>().ok())
            {
                peak = peak.max(n);
            }
        }
        match child.try_wait().expect("wait for cyrius-rs") {
            Some(status) => break status,
            None => std::thread::sleep(Duration::from_millis(5)),
        }
    };

    let _ = std::fs::remove_dir_all(&out_dir);

    assert!(status.success(), "cyrius-rs exited with {}", status);
    assert!(peak > 0, "never sampled a thread count from {}", status_path);
    assert!(
        peak <= BUDGET,
        "--diplotype-caller -t 1 peaked at {} threads on a {}-CPU machine, \
         over the budget of {}: --threads is not bounding rayon's global pool",
        peak,
        available,
        BUDGET
    );
    eprintln!("peak threads at -t 1 = {} (budget {})", peak, BUDGET);
}

//! Stamp the git commit into the build as `CYRIUS_GIT_COMMIT`: the short hash
//! of HEAD, with `-dirty` when the crate directory has uncommitted changes, or
//! `unknown` when git is unavailable (e.g. a source tarball).

use std::path::{Path, PathBuf};
use std::process::Command;

/// Run git in the crate directory; `None` on any failure.
fn git(args: &[&str]) -> Option<String> {
    let out = Command::new("git")
        .args(args)
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .ok()?;
    if !out.status.success() {
        return None;
    }
    Some(String::from_utf8_lossy(&out.stdout).trim().to_string())
}

fn rerun_if_exists(path: &Path) {
    if path.exists() {
        println!("cargo:rerun-if-changed={}", path.display());
    }
}

fn main() {
    let commit = match git(&["rev-parse", "--short=7", "HEAD"]) {
        Some(hash) if !hash.is_empty() => {
            // `-- .` limits the check to this crate's directory.
            let dirty = git(&["status", "--porcelain", "--", "."]).is_some_and(|s| !s.is_empty());
            if dirty {
                format!("{hash}-dirty")
            } else {
                hash
            }
        }
        _ => "unknown".to_string(),
    };
    println!("cargo:rustc-env=CYRIUS_GIT_COMMIT={commit}");

    // Rerun after a commit or checkout. In a worktree HEAD lives in the
    // worktree's git dir, while branch refs and packed-refs live in the
    // common dir. The paths git prints are relative to the crate directory.
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    if let Some(git_dir) = git(&["rev-parse", "--git-dir"]) {
        let git_dir = manifest_dir.join(git_dir);
        let head = git_dir.join("HEAD");
        rerun_if_exists(&head);
        if let Some(common_dir) = git(&["rev-parse", "--git-common-dir"]) {
            let common_dir = manifest_dir.join(common_dir);
            rerun_if_exists(&common_dir.join("packed-refs"));
            if let Ok(head) = std::fs::read_to_string(&head) {
                if let Some(reference) = head.trim().strip_prefix("ref: ") {
                    rerun_if_exists(&common_dir.join(reference));
                }
            }
        }
    }
    // ...and when the crate's own files change, which is what flips `-dirty`.
    for path in ["build.rs", "Cargo.toml", "Cargo.lock", "src", "data", "tests"] {
        rerun_if_exists(&manifest_dir.join(path));
    }
}

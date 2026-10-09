//! Process-wide thread-pool setup.
//!
//! Review finding #10: `--threads` claimed to bound the tool's parallelism but
//! nothing in `src` ever called `rayon::ThreadPoolBuilder`, so every rayon user
//! -- the control-region loop in `depth_calling::bin_count` and, until the
//! `--diplotype-caller` path was deleted, the `par_chunks` in
//! `caller::diplotype_caller` -- ran on rayon's *global* pool,
//! which sizes itself to every core on the machine. On a 72-core node `-t 2`
//! spawned ~72 rayon workers. `init_rayon_global_pool` sizes that global pool
//! from the flag, so `--threads` bounds **rayon's** parallelism.
//!
//! It does not bound the process's total thread count, and the flag's help text
//! ("Number of threads") should not be read as promising that it does. Three
//! other thread sources are deliberately untouched by this module:
//!
//! * `caller::construct_star_table::get_hap_table` spawns one scoped thread per
//!   combination table asked for -- at most one per run since finding #41, and
//!   none from `prepare_resource` -- not gated on `--threads` at all.
//! * `main.rs`'s two `std::thread::scope` blocks spawn a fixed 4 threads each
//!   whenever `threads > 1`, regardless of the value.
//! * every alignment-file open calls `set_threads(1)`, adding an htslib worker.

use std::sync::OnceLock;

/// Size of the global pool *this module* successfully built, if it built one.
///
/// Recorded so the `Err` arm can report a real number without having to ask
/// rayon -- see [`init_rayon_global_pool`].
static POOL_THREADS: OnceLock<usize> = OnceLock::new();

/// Parallelism actually available to this process.
///
/// On Linux `std::thread::available_parallelism` honours both
/// `sched_getaffinity` and the cgroup v1/v2 CPU quota, so inside a container
/// limited to 2 CPUs this returns 2, not the host's core count. That is the
/// property that makes it the right ceiling here: the cap tracks the CPU budget
/// the process was actually given rather than the hardware it happens to sit
/// on. It can only fail on platforms that cannot report it, where 1 is the
/// safe answer.
fn available_parallelism() -> usize {
    std::thread::available_parallelism().map_or(1, |n| n.get())
}

/// Clamp a `--threads` request into `1..=cap`.
///
/// Split out from [`rayon_pool_threads`] so both ends can be pinned by tests
/// that do not depend on the machine they run on. `cap` is itself floored at 1
/// so a `cap` of 0 can never produce a request of 0.
fn clamp_pool_threads(requested: usize, cap: usize) -> usize {
    requested.clamp(1, cap.max(1))
}

/// Pool size to request from rayon for a given `--threads` value.
///
/// Both ends are clamped, because rayon forwards this number to the OS as-is:
///
/// * **Floor.** `rayon`'s `num_threads(0)` means "use the default", i.e. one
///   worker per core -- exactly the unbounded behaviour this module exists to
///   remove -- so a `--threads 0` (which clap accepts; `usize` parses it fine)
///   must never be forwarded verbatim. It becomes 1, matching the flag's own
///   `default_value = "1"`.
/// * **Ceiling.** `num_threads(n)` really does try to spawn `n` OS threads.
///   Forwarded verbatim, `--threads 100000` spawned >10000 threads and reserved
///   59-65 GB of address space in stacks (two measurements) before making
///   progress; the run had not finished 20 s in, where the same run without any
///   pool sizing at all finished in 5 s on 73 threads and 10 GB. Capping at the
///   available parallelism keeps an absurd request from being strictly worse
///   than the unbounded pool this module replaced.
///
/// The ceiling is `available_parallelism` rather than a large constant on
/// purpose: rayon's work here (diplotype scoring, control-region counting) is
/// CPU-bound, so oversubscribing buys nothing, and the control-region loop
/// chunks to `n_cores` anyway. The cost is that a cgroup-limited process cannot
/// oversubscribe its quota deliberately -- an acceptable trade, since a quota
/// of 1 genuinely cannot run 8 workers in parallel.
///
/// Clamping down is logged at `warn`: a user who asked for 100000 threads
/// should learn what they got instead of silently getting something else.
pub fn rayon_pool_threads(requested: usize) -> usize {
    let cap = available_parallelism();
    let n = clamp_pool_threads(requested, cap);
    if requested > n {
        log::warn!(
            "--threads {} exceeds the {} CPU(s) available to this process; \
             using {} rayon thread(s) instead. Rayon spawns one OS thread per \
             requested thread, so the original request would have cost thread \
             stacks and scheduler time without buying any parallelism.",
            requested,
            cap,
            n
        );
    } else if requested == 0 {
        log::warn!(
            "--threads 0 means \"one worker per core\" to rayon, which is the \
             unbounded behaviour --threads exists to avoid; using 1 rayon \
             thread instead."
        );
    }
    n
}

/// Size rayon's global thread pool from `--threads`.
///
/// Returns `Some(n)` when the global pool is known to have `n` workers, and
/// `None` when it could not be sized and its size is unknown.
///
/// `build_global` may only succeed once per process and returns `Err` if the
/// pool is already initialised (for instance because some earlier rayon call
/// built the default pool). That is a performance problem, not a correctness
/// one -- rayon's work is deterministic here either way -- so it is logged and
/// the run continues. Panicking at start-up would turn a resource-usage bug
/// into a failed run.
pub fn init_rayon_global_pool(requested: usize) -> Option<usize> {
    let n = rayon_pool_threads(requested);
    match rayon::ThreadPoolBuilder::new().num_threads(n).build_global() {
        Ok(()) => {
            // `build_global` succeeded, so the pool has exactly `n` workers.
            let _ = POOL_THREADS.set(n);
            log::debug!("Rayon global thread pool sized to {} thread(s).", n);
            Some(n)
        }
        // Deliberately *not* `rayon::current_num_threads()` here. That call
        // lazily initialises rayon's default, all-cores global registry when
        // none exists, so if `build_global` ever failed for a thread-spawn
        // reason rather than "already initialised", asking for the live size
        // would make this error handler build the very unbounded pool it is
        // warning about. Report only what we already know.
        Err(e) => match POOL_THREADS.get().copied() {
            Some(live) => {
                log::debug!(
                    "Rayon's global thread pool was already sized to {} thread(s) \
                     by this process, so the request for {} was a no-op ({}).",
                    live,
                    n,
                    e
                );
                Some(live)
            }
            None => {
                log::warn!(
                    "Could not size rayon's global thread pool to {} thread(s): {}. \
                     Continuing on whatever pool is live; results are unaffected, \
                     but --threads will not bound rayon in this run.",
                    n,
                    e
                );
                None
            }
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clamp_pins_both_ends() {
        // Machine-independent: the floor, the identity range and the ceiling.
        assert_eq!(clamp_pool_threads(0, 8), 1);
        assert_eq!(clamp_pool_threads(1, 8), 1);
        assert_eq!(clamp_pool_threads(7, 8), 7);
        assert_eq!(clamp_pool_threads(8, 8), 8);
        assert_eq!(clamp_pool_threads(9, 8), 8);
        assert_eq!(clamp_pool_threads(100_000, 8), 8);
        assert_eq!(clamp_pool_threads(usize::MAX, 8), 8);
        // A cap of 0 must still leave a usable request, not 0.
        assert_eq!(clamp_pool_threads(4, 0), 1);
        assert_eq!(clamp_pool_threads(0, 0), 1);
    }

    #[test]
    fn requests_the_machine_can_honour_are_passed_through() {
        let available = available_parallelism();
        for requested in [1usize, 2, 4, 8] {
            if requested <= available {
                assert_eq!(rayon_pool_threads(requested), requested);
            }
        }
        assert_eq!(rayon_pool_threads(available), available);
    }

    #[test]
    fn zero_threads_never_requests_rayons_all_cores_default() {
        // num_threads(0) means "every core" to rayon, which is the bug this
        // module fixes, so -t 0 must become 1 rather than be forwarded.
        assert_eq!(rayon_pool_threads(0), 1);
    }

    #[test]
    fn absurd_requests_are_capped_to_the_available_parallelism() {
        // Forwarded verbatim, -t 100000 made rayon spawn >10000 OS threads and
        // reserve ~59 GB of address space -- strictly worse than the unbounded
        // pool this module replaced. The ceiling is the whole point of the clamp.
        let available = available_parallelism();
        assert_eq!(rayon_pool_threads(100_000), available);
        assert_eq!(rayon_pool_threads(usize::MAX), available);
        assert!(rayon_pool_threads(usize::MAX) >= 1);
    }
}

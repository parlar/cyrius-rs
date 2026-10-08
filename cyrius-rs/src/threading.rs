//! Process-wide thread-pool setup.
//!
//! Review finding #10: `--threads` claimed to bound the tool's parallelism but
//! nothing in `src` ever called `rayon::ThreadPoolBuilder`, so every rayon user
//! -- the control-region loop in `depth_calling::bin_count` and the
//! `par_chunks` in `caller::diplotype_caller` -- ran on rayon's *global* pool,
//! which sizes itself to every core on the machine. On a 72-core node `-t 2`
//! spawned ~72 rayon workers. `init_rayon_global_pool` sizes that global pool
//! from the flag so `--threads` actually bounds the work it claims to bound.

/// Pool size to request from rayon for a given `--threads` value.
///
/// `rayon`'s `num_threads(0)` means "use the default", i.e. one worker per
/// core -- exactly the unbounded behaviour this fix exists to remove -- so a
/// `--threads 0` (which clap accepts, `usize` parses it fine) must never be
/// forwarded verbatim. It is clamped to a single thread, matching the flag's
/// own `default_value = "1"`: asking for no threads is treated as asking for
/// the minimum, not for everything.
pub fn rayon_pool_threads(requested: usize) -> usize {
    requested.max(1)
}

/// Size rayon's global thread pool from `--threads`.
///
/// Returns the number of threads the global pool actually has afterwards.
///
/// `build_global` may only succeed once per process and returns `Err` if the
/// pool is already initialised (for instance because some earlier rayon call
/// built the default pool). That is a performance problem, not a correctness
/// one -- rayon's work is deterministic here either way -- so it is logged and
/// the run continues. Panicking at start-up would turn a resource-usage bug
/// into a failed run.
pub fn init_rayon_global_pool(requested: usize) -> usize {
    let n = rayon_pool_threads(requested);
    match rayon::ThreadPoolBuilder::new().num_threads(n).build_global() {
        Ok(()) => log::debug!("Rayon global thread pool sized to {} thread(s).", n),
        Err(e) => log::warn!(
            "Could not size rayon's global thread pool to {} thread(s) ({}); \
             continuing with the existing pool of {}. Results are unaffected; \
             --threads will not bound rayon in this run.",
            n,
            e,
            rayon::current_num_threads()
        ),
    }
    rayon::current_num_threads()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn requested_thread_count_is_passed_through() {
        assert_eq!(rayon_pool_threads(1), 1);
        assert_eq!(rayon_pool_threads(2), 2);
        assert_eq!(rayon_pool_threads(4), 4);
        assert_eq!(rayon_pool_threads(72), 72);
    }

    #[test]
    fn zero_threads_never_requests_rayons_all_cores_default() {
        // num_threads(0) means "every core" to rayon, which is the bug this
        // module fixes, so -t 0 must become 1 rather than be forwarded.
        assert_eq!(rayon_pool_threads(0), 1);
    }
}

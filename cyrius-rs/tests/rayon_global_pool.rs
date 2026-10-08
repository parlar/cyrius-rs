//! The global rayon pool is process-wide and can only be initialised once, so
//! it cannot be exercised from the lib unit-test binary (other tests there
//! touch rayon and would initialise it first, non-deterministically). This
//! integration test owns its own process, so it can assert the real effect of
//! `init_rayon_global_pool`: that `--threads` actually bounds rayon.

/// `--threads N` must size the global pool to exactly N, not to every core on
/// the machine. Review finding #10: before this fix nothing ever called
/// `ThreadPoolBuilder`, so rayon sized itself to `nproc` whatever `-t` said.
///
/// The assertion only *means* that on a machine with more than 2 CPUs: an
/// unbounded pool on a 1- or 2-core runner is also of size 2, so the test would
/// pass with the bug present. Skip rather than assert something vacuous.
#[test]
fn init_rayon_global_pool_bounds_the_pool_to_the_request() {
    let available = std::thread::available_parallelism().map_or(1, |n| n.get());
    if available <= 2 {
        eprintln!(
            "skipping: {} CPU(s) available, so a pool of 2 cannot be told apart \
             from rayon's all-cores default",
            available
        );
        return;
    }

    // Sized to 2 (not `available`, which is > 2 here, as checked above).
    assert_eq!(cyrius_rs::threading::init_rayon_global_pool(2), Some(2));
    assert_eq!(rayon::current_num_threads(), 2);

    // Second call cannot re-size an already-built global pool; it must report
    // the live size rather than panic, and must do so without asking rayon
    // (which would build the default pool if none existed).
    assert_eq!(cyrius_rs::threading::init_rayon_global_pool(8), Some(2));
    assert_eq!(rayon::current_num_threads(), 2);
}

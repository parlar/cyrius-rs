//! The global rayon pool is process-wide and can only be initialised once, so
//! it cannot be exercised from the lib unit-test binary (other tests there
//! touch rayon and would initialise it first, non-deterministically). This
//! integration test owns its own process, so it can assert the real effect of
//! `init_rayon_global_pool`: that `--threads` actually bounds rayon.

/// `--threads N` must size the global pool to exactly N, not to every core on
/// the machine. Review finding #10: before this fix nothing ever called
/// `ThreadPoolBuilder`, so rayon sized itself to `nproc` whatever `-t` said.
#[test]
fn init_rayon_global_pool_bounds_the_pool_to_the_request() {
    // Sized to 2 (not `nproc`) even though this machine has many more cores.
    assert_eq!(cyrius_rs::threading::init_rayon_global_pool(2), 2);

    // Second call cannot re-size an already-built global pool; it must report
    // the live size rather than panic.
    assert_eq!(cyrius_rs::threading::init_rayon_global_pool(8), 2);
}

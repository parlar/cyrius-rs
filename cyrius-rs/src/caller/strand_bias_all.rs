//! Extended strand bias filtering for ALL variant sites.
//!
//! The original code only applies Fisher's exact strand bias test to NOISY_VAR sites.
//! This module applies it to non-NOISY sites with different thresholds:
//! - CLEAN_VAR: lenient (p < 0.01 AND one strand has zero reads)
//! - Other sites: moderate (p < 0.02 AND (forward <= 1 or reverse <= 1))
//!
//! NOISY_VAR sites are skipped here — they are already handled by call_cn_var().

use crate::caller::call_variants::CLEAN_VAR;
use crate::fisher;

const P_CUTOFF_MODERATE: f64 = 0.02;
const P_CUTOFF_LENIENT: f64 = 0.01;

/// Minimum total alt reads before strand bias check applies.
/// Sites with very few reads shouldn't be filtered by strand bias alone.
const MIN_ALT_FOR_BIAS_CHECK: usize = 3;

/// Check strand bias for a non-NOISY variant site.
/// Returns true if the variant should be zeroed out.
fn has_strand_bias(var_name: &str, forward: usize, reverse: usize) -> bool {
    let total = forward + reverse;
    if total < MIN_ALT_FOR_BIAS_CHECK {
        return false;
    }

    let half = total as f64 / 2.0;
    let (_, pvalue) = fisher::fisher_exact([
        [forward as u64, reverse as u64],
        [half as u64, half as u64],
    ]);

    if CLEAN_VAR.contains(&var_name) {
        // Lenient: only flag extreme bias (significant p AND completely one-sided)
        pvalue < P_CUTOFF_LENIENT && (forward == 0 || reverse == 0)
    } else {
        // Moderate: significant p AND nearly one-sided
        pvalue < P_CUTOFF_MODERATE && (forward <= 1 || reverse <= 1)
    }
}

/// Apply strand bias filtering to all non-NOISY variant sites.
/// Modifies var_alt in-place, zeroing out sites that fail the strand bias check.
/// Returns the number of sites that were filtered.
pub fn apply_strand_bias_all(
    var_alt: &mut [usize],
    alt_forward: &[usize],
    alt_reverse: &[usize],
    var_list: &[String],
) -> usize {
    use crate::caller::call_variants::NOISY_VAR;

    let mut filtered_count = 0;
    for i in 0..var_alt.len().min(alt_forward.len()) {
        if var_alt[i] == 0 {
            continue;
        }
        let var_name = &var_list[i];
        // Skip NOISY_VAR — they are already handled by call_cn_var()
        if NOISY_VAR.contains(&var_name.as_str()) {
            continue;
        }
        if has_strand_bias(var_name, alt_forward[i], alt_reverse[i]) {
            log::debug!(
                "strand_bias_all: filtering {} (fwd={}, rev={})",
                var_name, alt_forward[i], alt_reverse[i]
            );
            var_alt[i] = 0;
            filtered_count += 1;
        }
    }
    filtered_count
}

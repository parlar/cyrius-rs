/// Fisher exact test for 2x2 contingency tables.
/// Matches scipy.stats.fisher_exact output (odds ratio, two-sided p-value).

use statrs::function::factorial::ln_factorial;

/// Compute the hypergeometric probability for a 2x2 table.
fn hypergeom_pmf(a: u64, b: u64, c: u64, d: u64) -> f64 {
    let n = a + b + c + d;
    (ln_factorial(a + b) + ln_factorial(c + d) + ln_factorial(a + c) + ln_factorial(b + d)
        - ln_factorial(a)
        - ln_factorial(b)
        - ln_factorial(c)
        - ln_factorial(d)
        - ln_factorial(n))
    .exp()
}

/// Fisher exact test for a 2x2 contingency table [[a, b], [c, d]].
/// Returns (odds_ratio, two_sided_p_value) matching scipy.stats.fisher_exact.
pub fn fisher_exact(table: [[u64; 2]; 2]) -> (f64, f64) {
    let a = table[0][0];
    let b = table[0][1];
    let c = table[1][0];
    let d = table[1][1];

    // Odds ratio (matching scipy.stats.fisher_exact behavior)
    let numerator = a as f64 * d as f64;
    let denominator = b as f64 * c as f64;
    let odds_ratio = if denominator == 0.0 {
        if numerator == 0.0 {
            f64::NAN
        } else {
            f64::INFINITY
        }
    } else {
        numerator / denominator
    };

    // Two-sided p-value: sum probabilities of all tables as or less probable
    // than the observed table
    let p_observed = hypergeom_pmf(a, b, c, d);

    let row1 = a + b;
    let col1 = a + c;
    let n = a + b + c + d;

    let a_min = if col1 + row1 > n {
        (col1 + row1 - n) as u64
    } else {
        0
    };
    let a_max = std::cmp::min(row1, col1);

    let mut p_value = 0.0;
    for ai in a_min..=a_max {
        let bi = row1 - ai;
        let ci = col1 - ai;
        let di = n - ai - bi - ci;
        let p = hypergeom_pmf(ai, bi, ci, di);
        if p <= p_observed + 1e-12 {
            // small tolerance for floating point
            p_value += p;
        }
    }

    // Clamp to [0, 1]
    let p_value = p_value.min(1.0);

    (odds_ratio, p_value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fisher_exact_basic() {
        // scipy.stats.fisher_exact([[1, 9], [11, 3]]) -> (0.0303..., 0.00137...)
        let (or, pv) = fisher_exact([[1, 9], [11, 3]]);
        assert!((or - 1.0 * 3.0 / (9.0 * 11.0)).abs() < 1e-10);
        assert!(pv < 0.01);
    }

    #[test]
    fn test_fisher_exact_symmetric() {
        let (_, pv) = fisher_exact([[5, 5], [5, 5]]);
        assert!((pv - 1.0).abs() < 1e-6);
    }
}

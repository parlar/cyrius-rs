/// Statistical helper functions replacing numpy.median, numpy.mean, etc.

pub fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n % 2 == 0 {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    }
}

pub fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.iter().sum::<f64>() / values.len() as f64
}

/// Median Absolute Deviation (MAD) with the standard constant factor 1.4826.
pub const MAD_CONSTANT: f64 = 1.4826;

pub fn mad(values: &[f64]) -> f64 {
    let med = median(values);
    let abs_devs: Vec<f64> = values.iter().map(|&a| (a - med).abs()).collect();
    MAD_CONSTANT * median(&abs_devs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_median_odd() {
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
    }

    #[test]
    fn test_median_even() {
        assert_eq!(median(&[4.0, 1.0, 3.0, 2.0]), 2.5);
    }

    #[test]
    fn test_mean() {
        assert!((mean(&[1.0, 2.0, 3.0, 4.0]) - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_mad() {
        let vals = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let result = mad(&vals);
        // MAD = 1.4826 * median(|x - 3|) = 1.4826 * 1.0 = 1.4826
        assert!((result - 1.4826).abs() < 1e-10);
    }
}

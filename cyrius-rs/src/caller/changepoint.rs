//! Paralog ratio changepoint detection for hybrid identification.
//!
//! Runs circular binary segmentation (CBS) on the D6/(D6+D7) ratio profile
//! to detect step changes that indicate hybrid alleles (gene conversion events).
//!
//! When the standard `get_cnvtag` approach can't determine the CNV group,
//! changepoint detection may reveal hybrid structure from the spatial pattern
//! of D6 vs D7 signal across the gene.

/// A segment of roughly constant paralog ratio.
#[derive(Debug, Clone)]
pub struct Segment {
    pub start_idx: usize,
    pub end_idx: usize,
    pub mean_ratio: f64,
    pub n_valid: usize,
}

/// Result of changepoint analysis.
#[derive(Debug, Clone)]
pub struct ChangePointResult {
    pub segments: Vec<Segment>,
    pub n_changepoints: usize,
    pub is_hybrid: bool,
    pub overall_ratio: f64,
}

/// Compute paralog ratios at each diagnostic position.
fn compute_ratios(snp_d6: &[usize], snp_d7: &[usize], min_depth: usize) -> (Vec<f64>, Vec<bool>) {
    let n = snp_d6.len().min(snp_d7.len());
    let mut ratios = vec![f64::NAN; n];
    let mut mask = vec![false; n];

    for i in 0..n {
        let total = snp_d6[i] + snp_d7[i];
        if total >= min_depth {
            ratios[i] = snp_d6[i] as f64 / total as f64;
            mask[i] = true;
        }
    }

    (ratios, mask)
}

/// CBS recursive binary segmentation using Welch's t-test.
fn cbs_recurse(
    values: &[f64],
    start: usize,
    end: usize,
    t_thresh: f64,
    min_bins: usize,
) -> Vec<(usize, usize, f64)> {
    let segment = &values[start..end];
    let valid: Vec<f64> = segment.iter().filter(|v| !v.is_nan()).copied().collect();

    if valid.len() < 2 * min_bins {
        let mean = if valid.is_empty() {
            f64::NAN
        } else {
            valid.iter().sum::<f64>() / valid.len() as f64
        };
        return vec![(start, end, mean)];
    }

    let mut best_t = 0.0_f64;
    let mut best_split = None;

    for split in min_bins..=(segment.len() - min_bins) {
        let left: Vec<f64> = segment[..split].iter().filter(|v| !v.is_nan()).copied().collect();
        let right: Vec<f64> = segment[split..].iter().filter(|v| !v.is_nan()).copied().collect();

        if left.len() < min_bins || right.len() < min_bins {
            continue;
        }

        let t_abs = welch_t_stat(&left, &right).abs();
        if t_abs > best_t {
            best_t = t_abs;
            best_split = Some(split);
        }
    }

    if best_t < t_thresh || best_split.is_none() {
        let mean = valid.iter().sum::<f64>() / valid.len() as f64;
        return vec![(start, end, mean)];
    }

    let split = best_split.unwrap();
    let mut left_segs = cbs_recurse(values, start, start + split, t_thresh, min_bins);
    let right_segs = cbs_recurse(values, start + split, end, t_thresh, min_bins);
    left_segs.extend(right_segs);
    left_segs
}

/// Welch's t-statistic for two samples.
fn welch_t_stat(a: &[f64], b: &[f64]) -> f64 {
    let n1 = a.len() as f64;
    let n2 = b.len() as f64;
    if n1 < 2.0 || n2 < 2.0 {
        return 0.0;
    }

    let mean1 = a.iter().sum::<f64>() / n1;
    let mean2 = b.iter().sum::<f64>() / n2;

    let var1 = a.iter().map(|x| (x - mean1).powi(2)).sum::<f64>() / (n1 - 1.0);
    let var2 = b.iter().map(|x| (x - mean2).powi(2)).sum::<f64>() / (n2 - 1.0);

    let se = (var1 / n1 + var2 / n2).sqrt();
    if se < 1e-10 {
        return 0.0;
    }

    (mean1 - mean2) / se
}

/// Run changepoint detection on the paralog ratio profile.
///
/// Uses CBS with a two-pass approach: primary pass at t_thresh=3.0,
/// then a secondary pass at t_thresh=2.0 if the overall ratio suggests
/// imbalance but no changepoint was found.
pub fn detect_changepoints(
    snp_d6: &[usize],
    snp_d7: &[usize],
    min_depth: usize,
    t_thresh: f64,
    min_bins: usize,
) -> ChangePointResult {
    let (ratios, mask) = compute_ratios(snp_d6, snp_d7, min_depth);

    let valid_ratios: Vec<f64> = ratios.iter().zip(mask.iter())
        .filter(|(_, &m)| m)
        .map(|(&r, _)| r)
        .collect();

    let overall_ratio = if valid_ratios.is_empty() {
        0.0
    } else {
        valid_ratios.iter().sum::<f64>() / valid_ratios.len() as f64
    };

    let mut raw_segments = cbs_recurse(&ratios, 0, ratios.len(), t_thresh, min_bins);

    // Two-pass: if primary found only 1 segment and ratio suggests imbalance,
    // try again with lower threshold
    if raw_segments.len() == 1 && !(0.40..=0.60).contains(&overall_ratio) {
        let secondary = cbs_recurse(&ratios, 0, ratios.len(), (t_thresh - 1.0).max(2.0), min_bins);
        if secondary.len() > 1 {
            raw_segments = secondary;
        }
    }

    let segments: Vec<Segment> = raw_segments
        .iter()
        .map(|&(start, end, mean)| {
            let n_valid = mask[start..end].iter().filter(|&&m| m).count();
            Segment {
                start_idx: start,
                end_idx: end,
                mean_ratio: if mean.is_nan() { 0.0 } else { mean },
                n_valid,
            }
        })
        .collect();

    let n_changepoints = segments.len().saturating_sub(1);

    // Detect hybrid from segment ratio differences
    let min_ratio_diff = 0.15;
    let min_segment_valid = 5;
    let is_hybrid = segments.windows(2).any(|pair| {
        let diff = (pair[0].mean_ratio - pair[1].mean_ratio).abs();
        diff >= min_ratio_diff
            && pair[0].n_valid >= min_segment_valid
            && pair[1].n_valid >= min_segment_valid
    });

    ChangePointResult {
        segments,
        n_changepoints,
        is_hybrid,
        overall_ratio,
    }
}

/// Classify changepoint result into a suggested CNV modification.
///
/// Returns Some(suggested_cnv_tag) if the changepoint pattern suggests
/// a hybrid that the standard pipeline might have missed.
pub fn suggest_cnv_from_changepoints(
    result: &ChangePointResult,
    current_cnvtag: Option<&str>,
    total_cn: u32,
) -> Option<String> {
    if !result.is_hybrid {
        return None;
    }

    // Only suggest when the standard pipeline couldn't determine CNV group
    if current_cnvtag.is_some() {
        return None;
    }

    // A clean hybrid should produce exactly 2 segments (or 3 with a small
    // transition zone). More segments indicates noisy data, not a hybrid.
    if result.segments.len() < 2 || result.segments.len() > 3 {
        return None;
    }

    // Both the first and last segment must have substantial data
    let first = &result.segments[0];
    let last = &result.segments[result.segments.len() - 1];

    if first.n_valid < 10 || last.n_valid < 10 {
        return None;
    }

    // SNP array goes from 3' (downstream, exon9) to 5' (upstream, exon1).
    // *68 (gene-pseudo): D6 at 5' (last), D7 at 3' (first) -> first LOW, last HIGH -> diff < 0
    // *13 (pseudo-gene): D7 at 5' (last), D6 at 3' (first) -> first HIGH, last LOW -> diff > 0
    // Require a substantial step (>0.20) for confidence
    let diff = first.mean_ratio - last.mean_ratio;

    if diff > 0.20 {
        // 3' D6-rich, 5' D7-rich -> *13-like hybrid (pseudo-gene)
        match total_cn {
            2 => Some("star13".to_string()),
            3 => Some("dup_star13".to_string()),
            _ => None,
        }
    } else if diff < -0.20 {
        // 3' D7-rich, 5' D6-rich -> *68-like hybrid (gene-pseudo)
        match total_cn {
            2 => Some("star5_star68".to_string()),
            3 => Some("star68".to_string()),
            4 => Some("star68_star68".to_string()),
            _ => None,
        }
    } else {
        None
    }
}

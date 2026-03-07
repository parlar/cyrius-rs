//! HMM-based CNV segmentation and classification for CYP2D6.
//!
//! Uses a Hidden Markov Model on the per-SNP D6 copy number profile to:
//! 1. Segment the gene into regions of constant CYP2D6 copy number (Viterbi)
//! 2. Classify the segment pattern into a CNV group tag
//!
//! This is an alternative to the consensus-based approach in `cnv_hybrid.rs`.
//! The HMM naturally handles noisy observations and detects breakpoints from
//! the spatial pattern of D6/D7 signal across the gene.
//!
//! # Usage
//! ```ignore
//! let params = HmmParams::from_gmm_prior(total_cn, sd);
//! let result = hmm_segment(&raw_d6_cn, &snp_d6, &snp_d7, total_cn, &params);
//! let cnvtag = classify_segments(&result.segments, total_cn, spacer_cn);
//! ```

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// CYP2D6 gene structure: SNP site index boundaries
// These match the positions used in cnv_hybrid.rs
// ---------------------------------------------------------------------------
const REP_END: usize = 3;
const EXON9_END: usize = 9;
const INTRON4_BP: usize = 40;
const INTRON1_BP: usize = 74;

// ---------------------------------------------------------------------------
// HMM parameters
// ---------------------------------------------------------------------------

/// Parameters for the HMM.
#[derive(Debug, Clone)]
pub struct HmmParams {
    /// Number of CN states (0..max_cn inclusive)
    pub n_states: usize,
    /// Emission standard deviation per state (indexed by CN)
    pub emission_sd: Vec<f64>,
    /// Log prior probability per state
    pub log_prior: Vec<f64>,
    /// Log probability of staying in the same state
    pub log_stay: f64,
    /// Log probability of transitioning to a different state
    pub log_switch: f64,
    /// Minimum total read depth at a site to use the observation
    pub min_depth: usize,
}

impl HmmParams {
    /// Construct default parameters from the total CN call and an estimated SD.
    ///
    /// `sd_per_copy` is the standard deviation of the D6-CN observation per
    /// unit of copy number — typically ~0.3–0.5 for WGS data. The GMM's
    /// `sd_cn2` value (divided by 2 to get per-copy) is a reasonable source.
    pub fn from_gmm_prior(total_cn: u32, sd_per_copy: f64) -> Self {
        let max_cn = (total_cn + 2) as usize; // allow some room above total
        let n_states = max_cn + 1;

        // Emission SD scales with sqrt(CN) (matching GMM convention), floored
        let emission_sd: Vec<f64> = (0..n_states)
            .map(|cn| {
                if cn == 0 {
                    0.15 // tight distribution around 0
                } else {
                    sd_per_copy * (cn as f64).sqrt()
                }
            })
            .collect();

        // Flat prior — let the data speak
        let log_prior = vec![-(n_states as f64).ln(); n_states];

        // Transition: high self-transition probability (CN is locally constant)
        // p_stay = 0.995, p_switch = 0.005 / (n_states - 1)
        let p_stay: f64 = 0.995;
        let p_switch_each = (1.0 - p_stay) / (n_states - 1).max(1) as f64;

        HmmParams {
            n_states,
            emission_sd,
            log_prior,
            log_stay: p_stay.ln(),
            log_switch: p_switch_each.ln(),
            min_depth: 10,
        }
    }

    /// Build with custom stay probability (useful for tuning sensitivity).
    pub fn with_stay_prob(mut self, p_stay: f64) -> Self {
        let p_switch_each = (1.0 - p_stay) / (self.n_states - 1).max(1) as f64;
        self.log_stay = p_stay.ln();
        self.log_switch = p_switch_each.ln();
        self
    }
}

// ---------------------------------------------------------------------------
// HMM core: Viterbi decoding
// ---------------------------------------------------------------------------

/// Log-probability of observing `x` given CN state `cn`.
/// Emission model: Gaussian(mean=cn, sd=emission_sd[cn])
fn log_emission(x: f64, cn: usize, params: &HmmParams) -> f64 {
    let mean = cn as f64;
    let sd = params.emission_sd[cn];
    let diff = x - mean;
    // log N(x | mean, sd) = -0.5*ln(2*pi) - ln(sd) - 0.5*((x-mean)/sd)^2
    -0.5 * std::f64::consts::TAU.ln() - sd.ln() - 0.5 * (diff / sd).powi(2)
}

/// Log transition probability from state `from` to state `to`.
fn log_transition(from: usize, to: usize, params: &HmmParams) -> f64 {
    if from == to {
        params.log_stay
    } else {
        params.log_switch
    }
}

/// Run Viterbi decoding on a sequence of D6-CN observations.
///
/// `observations` is the per-SNP D6 copy number (e.g. `total_cn * d6_fraction`).
/// `valid` marks which sites have sufficient depth to trust.
///
/// Returns the most-likely state (CN) sequence, same length as `observations`.
pub fn viterbi(observations: &[f64], valid: &[bool], params: &HmmParams) -> Vec<usize> {
    let t = observations.len();
    if t == 0 {
        return vec![];
    }
    let n = params.n_states;

    // dp[t][s] = log prob of best path ending in state s at time t
    let mut dp = vec![vec![f64::NEG_INFINITY; n]; t];
    // backpointer[t][s] = previous state on best path
    let mut bp = vec![vec![0usize; n]; t];

    // Initialize
    for s in 0..n {
        dp[0][s] = params.log_prior[s]
            + if valid[0] {
                log_emission(observations[0], s, params)
            } else {
                0.0 // uninformative observation
            };
    }

    // Forward pass
    for i in 1..t {
        for s in 0..n {
            let e = if valid[i] {
                log_emission(observations[i], s, params)
            } else {
                0.0
            };
            let mut best_prev = 0;
            let mut best_score = f64::NEG_INFINITY;
            for p in 0..n {
                let score = dp[i - 1][p] + log_transition(p, s, params);
                if score > best_score {
                    best_score = score;
                    best_prev = p;
                }
            }
            dp[i][s] = best_score + e;
            bp[i][s] = best_prev;
        }
    }

    // Backtrace
    let mut path = vec![0usize; t];
    let mut best_final = 0;
    let mut best_score = f64::NEG_INFINITY;
    for s in 0..n {
        if dp[t - 1][s] > best_score {
            best_score = dp[t - 1][s];
            best_final = s;
        }
    }
    path[t - 1] = best_final;
    for i in (0..t - 1).rev() {
        path[i] = bp[i + 1][path[i + 1]];
    }

    path
}

// ---------------------------------------------------------------------------
// Segment extraction
// ---------------------------------------------------------------------------

/// A contiguous region of constant HMM-decoded CN.
#[derive(Debug, Clone, PartialEq)]
pub struct CnSegment {
    /// SNP site index where this segment starts (inclusive)
    pub start_idx: usize,
    /// SNP site index where this segment ends (inclusive)
    pub end_idx: usize,
    /// Decoded CN state for this segment
    pub cn: u32,
    /// Number of valid (high-depth) observations in this segment
    pub n_valid: usize,
    /// Mean observed D6-CN value in this segment (valid sites only)
    pub mean_obs: f64,
}

/// Extract contiguous segments from a Viterbi state path.
pub fn extract_segments(
    path: &[usize],
    observations: &[f64],
    valid: &[bool],
) -> Vec<CnSegment> {
    if path.is_empty() {
        return vec![];
    }

    let mut segments = Vec::new();
    let mut seg_start = 0;

    for i in 1..=path.len() {
        if i == path.len() || path[i] != path[seg_start] {
            // Close current segment
            let mut n_valid = 0;
            let mut sum_obs = 0.0;
            for j in seg_start..i {
                if valid[j] {
                    n_valid += 1;
                    sum_obs += observations[j];
                }
            }
            segments.push(CnSegment {
                start_idx: seg_start,
                end_idx: i - 1,
                cn: path[seg_start] as u32,
                n_valid,
                mean_obs: if n_valid > 0 {
                    sum_obs / n_valid as f64
                } else {
                    path[seg_start] as f64
                },
            });
            if i < path.len() {
                seg_start = i;
            }
        }
    }

    segments
}

// ---------------------------------------------------------------------------
// Full HMM segmentation pipeline
// ---------------------------------------------------------------------------

/// Result of HMM segmentation.
#[derive(Debug, Clone)]
pub struct HmmResult {
    /// Viterbi state path (one CN per SNP site)
    pub path: Vec<usize>,
    /// Extracted segments
    pub segments: Vec<CnSegment>,
    /// Per-site D6 CN observations used
    pub observations: Vec<f64>,
    /// Per-site validity mask
    pub valid: Vec<bool>,
}

/// Run full HMM segmentation on D6/D7 SNP data.
///
/// * `raw_d6_cn` — per-SNP D6 copy number (`total_cn * d6_fraction`), already computed
/// * `snp_d6` / `snp_d7` — raw read counts used only for depth filtering
/// * `total_cn` — GMM-called total D6+D7 copy number
/// * `params` — HMM parameters
pub fn hmm_segment(
    raw_d6_cn: &[f64],
    snp_d6: &[usize],
    snp_d7: &[usize],
    _total_cn: u32,
    params: &HmmParams,
) -> HmmResult {
    let n = raw_d6_cn.len().min(snp_d6.len()).min(snp_d7.len());

    let observations: Vec<f64> = raw_d6_cn[..n].to_vec();
    let valid: Vec<bool> = (0..n)
        .map(|i| snp_d6[i] + snp_d7[i] >= params.min_depth)
        .collect();

    let path = viterbi(&observations, &valid, params);
    let segments = extract_segments(&path, &observations, &valid);

    HmmResult {
        path,
        segments,
        observations,
        valid,
    }
}

// ---------------------------------------------------------------------------
// CNV classification from segments
// ---------------------------------------------------------------------------

/// Map the HMM segment pattern to a CNV group tag (same vocabulary as cnv_hybrid).
///
/// The classification logic mirrors the region-based approach:
/// - Look at the CN in each gene region (REP, exon9, intron4, intron1, upstream)
/// - Identify structural events from CN changes between adjacent regions
/// - Combine events into a CNV tag
///
/// Returns `None` if the pattern doesn't match any known configuration.
pub fn classify_segments(
    segments: &[CnSegment],
    total_cn: u32,
    _spacer_cn: Option<u32>,
) -> Option<String> {
    if segments.is_empty() {
        return None;
    }

    // Compute the dominant CN in each gene region by weighted vote
    let region_cn = get_region_cn(segments);

    let rep_cn = region_cn.get("rep").copied()?;
    let exon9_cn = region_cn.get("exon9").copied()?;
    let body_cn = region_cn.get("body").copied()?; // exon9-to-intron4
    let intron4_cn = region_cn.get("intron4").copied()?; // intron4-to-intron1
    let upstream_cn = region_cn.get("upstream").copied()?; // intron1-upstream

    // Basic sanity: upstream must be total_cn - 2 (REP region always = 2)
    if total_cn < 2 || upstream_cn != total_cn - 2 {
        return None;
    }

    // Collect structural events from CN transitions (same logic as cnv_hybrid)
    let mut events: Vec<String> = Vec::new();

    // intron4→intron1 vs upstream
    for _ in 0..intron4_cn.saturating_sub(upstream_cn) {
        events.push("star13intron1".to_string());
    }
    for _ in 0..upstream_cn.saturating_sub(intron4_cn) {
        events.push("star68".to_string());
    }

    // body (exon9-intron4) vs intron4-intron1
    for _ in 0..body_cn.saturating_sub(intron4_cn) {
        events.push("star13intron1".to_string());
    }

    // exon9 vs body
    for _ in 0..exon9_cn.saturating_sub(body_cn) {
        events.push("star13".to_string());
    }
    for _ in 0..body_cn.saturating_sub(exon9_cn) {
        events.push("exon9hyb".to_string());
    }

    // rep vs exon9
    for _ in 0..rep_cn.saturating_sub(exon9_cn) {
        events.push("star5".to_string());
    }
    for _ in 0..exon9_cn.saturating_sub(rep_cn) {
        events.push("dup".to_string());
    }

    if events.is_empty() {
        // Uniform CN across all regions
        if body_cn == 2 && intron4_cn == 2 && upstream_cn == 2 {
            return Some("cn2".to_string());
        }
        return None;
    }

    // Validate: events must produce the correct final CN
    let cn_increase = ["dup", "exon9hyb", "star68"];
    let cn_decrease = ["star13intron1", "star13", "star5"];
    let mut computed_cn: i32 = 2;
    for ev in &events {
        if cn_increase.contains(&ev.as_str()) {
            computed_cn += 1;
        }
        if cn_decrease.contains(&ev.as_str()) {
            computed_cn -= 1;
        }
    }
    if computed_cn != upstream_cn as i32 {
        return None;
    }

    events.sort();
    let tag = transform_tag(&events.join("_"));
    Some(tag)
}

/// Determine the dominant CN in each named gene region from HMM segments.
fn get_region_cn(segments: &[CnSegment]) -> HashMap<String, u32> {
    // Map each SNP index range to a region name
    let regions: &[(&str, usize, usize)] = &[
        ("rep", 0, REP_END),
        ("exon9", REP_END, EXON9_END),
        ("body", EXON9_END, INTRON4_BP),
        ("intron4", INTRON4_BP, INTRON1_BP),
        ("upstream", INTRON1_BP, usize::MAX), // extends to end
    ];

    let mut result = HashMap::new();

    for &(name, region_start, region_end) in regions {
        // Weighted vote: each segment contributes its overlap length
        let mut votes: HashMap<u32, usize> = HashMap::new();
        for seg in segments {
            let overlap_start = seg.start_idx.max(region_start);
            let overlap_end = (seg.end_idx + 1).min(region_end);
            if overlap_start < overlap_end {
                *votes.entry(seg.cn).or_insert(0) += overlap_end - overlap_start;
            }
        }
        if let Some((&cn, _)) = votes.iter().max_by_key(|&(_, &count)| count) {
            result.insert(name.to_string(), cn);
        }
    }

    result
}

/// Rename CNV tags to match the existing vocabulary (mirrors cnv_hybrid::transform_cnvtag).
fn transform_tag(tag: &str) -> String {
    let mut parts: Vec<String> = tag.split('_').map(|s| s.to_string()).collect();

    // Cancel out exon9hyb + star5 pairs (except the single exon9hyb_star5)
    if tag != "exon9hyb_star5" {
        while parts.contains(&"exon9hyb".to_string())
            && parts.contains(&"star5".to_string())
        {
            if let Some(pos) = parts.iter().position(|s| s == "exon9hyb") {
                parts.remove(pos);
            }
            if let Some(pos) = parts.iter().position(|s| s == "star5") {
                parts.remove(pos);
            }
        }
    }
    // Cancel out dup + star13 pairs (except the single dup_star13)
    if tag != "dup_star13" {
        while parts.contains(&"dup".to_string())
            && parts.contains(&"star13".to_string())
        {
            if let Some(pos) = parts.iter().position(|s| s == "dup") {
                parts.remove(pos);
            }
            if let Some(pos) = parts.iter().position(|s| s == "star13") {
                parts.remove(pos);
            }
        }
    }

    if parts.iter().all(|s| s == "dup") && !parts.is_empty() {
        return format!("cn{}", parts.len() + 2);
    }
    if tag == "dup_dup_exon9hyb_star13intron1" {
        return "cn4".to_string();
    }
    parts.join("_")
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_params(n_states: usize) -> HmmParams {
        HmmParams {
            n_states,
            emission_sd: (0..n_states)
                .map(|cn| if cn == 0 { 0.15 } else { 0.3 * (cn as f64).sqrt() })
                .collect(),
            log_prior: vec![-(n_states as f64).ln(); n_states],
            log_stay: 0.995_f64.ln(),
            log_switch: (0.005 / (n_states - 1) as f64).ln(),
            min_depth: 5,
        }
    }

    #[test]
    fn test_viterbi_uniform_cn2() {
        // All observations around CN=2 → should decode to all 2s
        let params = make_params(6);
        let obs: Vec<f64> = vec![2.1, 1.9, 2.0, 2.05, 1.95, 2.0, 2.1, 1.9];
        let valid = vec![true; obs.len()];
        let path = viterbi(&obs, &valid, &params);
        assert!(path.iter().all(|&s| s == 2), "Expected all CN=2, got {:?}", path);
    }

    #[test]
    fn test_viterbi_step_change() {
        // First half CN=2, second half CN=1 → should detect the step
        let params = make_params(6);
        let mut obs = Vec::new();
        // 20 sites at CN~2
        for _ in 0..20 {
            obs.push(2.0);
        }
        // 20 sites at CN~1
        for _ in 0..20 {
            obs.push(1.0);
        }
        let valid = vec![true; obs.len()];
        let path = viterbi(&obs, &valid, &params);

        // First half should be CN=2
        assert!(
            path[0..15].iter().all(|&s| s == 2),
            "First half should be CN=2, got {:?}",
            &path[0..15]
        );
        // Second half should be CN=1
        assert!(
            path[25..40].iter().all(|&s| s == 1),
            "Second half should be CN=1, got {:?}",
            &path[25..40]
        );
    }

    #[test]
    fn test_viterbi_with_noise() {
        // CN=2 with significant noise — should still decode as 2
        let params = make_params(6);
        let obs: Vec<f64> = vec![1.7, 2.3, 1.8, 2.2, 2.0, 1.6, 2.4, 2.1, 1.9, 2.0];
        let valid = vec![true; obs.len()];
        let path = viterbi(&obs, &valid, &params);
        assert!(path.iter().all(|&s| s == 2), "Expected all CN=2, got {:?}", path);
    }

    #[test]
    fn test_viterbi_invalid_sites() {
        // Some sites invalid — should not break decoding
        let params = make_params(6);
        let obs: Vec<f64> = vec![2.0, 999.0, 2.0, 2.0, 999.0, 2.0];
        let valid = vec![true, false, true, true, false, true];
        let path = viterbi(&obs, &valid, &params);
        assert!(path.iter().all(|&s| s == 2), "Expected all CN=2, got {:?}", path);
    }

    #[test]
    fn test_extract_segments_uniform() {
        let path = vec![2, 2, 2, 2, 2];
        let obs = vec![2.0, 2.1, 1.9, 2.0, 2.0];
        let valid = vec![true; 5];
        let segs = extract_segments(&path, &obs, &valid);
        assert_eq!(segs.len(), 1);
        assert_eq!(segs[0].cn, 2);
        assert_eq!(segs[0].start_idx, 0);
        assert_eq!(segs[0].end_idx, 4);
        assert_eq!(segs[0].n_valid, 5);
    }

    #[test]
    fn test_extract_segments_two_regions() {
        let path = vec![2, 2, 2, 1, 1, 1];
        let obs = vec![2.0, 2.0, 2.0, 1.0, 1.0, 1.0];
        let valid = vec![true; 6];
        let segs = extract_segments(&path, &obs, &valid);
        assert_eq!(segs.len(), 2);
        assert_eq!(segs[0].cn, 2);
        assert_eq!(segs[0].end_idx, 2);
        assert_eq!(segs[1].cn, 1);
        assert_eq!(segs[1].start_idx, 3);
    }

    #[test]
    fn test_classify_cn2() {
        // All regions at CN=2 → cn2
        let segments = vec![CnSegment {
            start_idx: 0,
            end_idx: 100,
            cn: 2,
            n_valid: 80,
            mean_obs: 2.0,
        }];
        let tag = classify_segments(&segments, 4, Some(2));
        assert_eq!(tag, Some("cn2".to_string()));
    }

    #[test]
    fn test_classify_star5() {
        // star5: rep=2 but exon9 onwards drops to CN=1
        // total_cn=3, upstream=1 → total_cn-2=1 ✓
        // rep(2) > exon9(1) → star5
        let segments = vec![
            CnSegment {
                start_idx: 0,
                end_idx: REP_END - 1,
                cn: 2,
                n_valid: 3,
                mean_obs: 2.0,
            },
            CnSegment {
                start_idx: REP_END,
                end_idx: 100,
                cn: 1,
                n_valid: 80,
                mean_obs: 1.0,
            },
        ];
        let tag = classify_segments(&segments, 3, Some(2));
        assert_eq!(tag, Some("star5".to_string()));
    }

    #[test]
    fn test_classify_dup() {
        // dup: rep=2, exon9 onwards = CN=3
        // total_cn=5, upstream=3 → total_cn-2=3 ✓
        // exon9(3) > rep(2) → dup
        let segments = vec![
            CnSegment {
                start_idx: 0,
                end_idx: REP_END - 1,
                cn: 2,
                n_valid: 3,
                mean_obs: 2.0,
            },
            CnSegment {
                start_idx: REP_END,
                end_idx: 100,
                cn: 3,
                n_valid: 80,
                mean_obs: 3.0,
            },
        ];
        let tag = classify_segments(&segments, 5, Some(2));
        assert_eq!(tag, Some("cn3".to_string()));
    }

    #[test]
    fn test_classify_exon9hyb() {
        // exon9hyb: exon9 region CN=2, body/intron4/upstream CN=3
        // total_cn=5, upstream=3 → total_cn-2=3 ✓
        // body(3) > exon9(2) → exon9hyb event
        // exon9(2) = rep(2) → no star5/dup
        let segments = vec![
            CnSegment {
                start_idx: 0,
                end_idx: EXON9_END - 1,
                cn: 2,
                n_valid: 9,
                mean_obs: 2.0,
            },
            CnSegment {
                start_idx: EXON9_END,
                end_idx: 100,
                cn: 3,
                n_valid: 80,
                mean_obs: 3.0,
            },
        ];
        let tag = classify_segments(&segments, 5, Some(2));
        assert_eq!(tag, Some("exon9hyb".to_string()));
    }

    #[test]
    fn test_classify_star68() {
        // star68: intron4-intron1=2, upstream=3
        // total_cn=5, upstream=3 → total_cn-2=3 ✓
        // upstream(3) > intron4(2) → star68
        let segments = vec![
            CnSegment {
                start_idx: 0,
                end_idx: INTRON1_BP - 1,
                cn: 2,
                n_valid: 60,
                mean_obs: 2.0,
            },
            CnSegment {
                start_idx: INTRON1_BP,
                end_idx: 100,
                cn: 3,
                n_valid: 25,
                mean_obs: 3.0,
            },
        ];
        let tag = classify_segments(&segments, 5, Some(2));
        assert_eq!(tag, Some("star68".to_string()));
    }

    #[test]
    fn test_classify_none_on_mismatch() {
        // upstream CN doesn't match total_cn - 2 → None
        let segments = vec![CnSegment {
            start_idx: 0,
            end_idx: 100,
            cn: 3,
            n_valid: 80,
            mean_obs: 3.0,
        }];
        let tag = classify_segments(&segments, 4, Some(2));
        assert_eq!(tag, None);
    }

    #[test]
    fn test_hmm_segment_full_pipeline() {
        // Simulate a star5 pattern: first few sites CN~2, rest CN~1
        let total_cn = 3u32;
        let params = HmmParams::from_gmm_prior(total_cn, 0.3);

        let n = 90;
        let mut raw_d6_cn = Vec::with_capacity(n);
        let mut snp_d6 = Vec::with_capacity(n);
        let mut snp_d7 = Vec::with_capacity(n);

        for i in 0..n {
            if i < REP_END {
                raw_d6_cn.push(2.0);
                snp_d6.push(30);
                snp_d7.push(15);
            } else {
                raw_d6_cn.push(1.0);
                snp_d6.push(15);
                snp_d7.push(30);
            }
        }

        let result = hmm_segment(&raw_d6_cn, &snp_d6, &snp_d7, total_cn, &params);

        assert_eq!(result.segments.len(), 2);
        assert_eq!(result.segments[0].cn, 2);
        assert_eq!(result.segments[1].cn, 1);

        let tag = classify_segments(&result.segments, total_cn, Some(2));
        assert_eq!(tag, Some("star5".to_string()));
    }

    #[test]
    fn test_transform_tag_dup_only() {
        assert_eq!(transform_tag("dup"), "cn3");
        assert_eq!(transform_tag("dup_dup"), "cn4");
        assert_eq!(transform_tag("dup_dup_dup"), "cn5");
    }

    #[test]
    fn test_transform_tag_cancellation() {
        // exon9hyb + star5 cancel (unless it's exactly the pair)
        assert_eq!(transform_tag("exon9hyb_star5"), "exon9hyb_star5");
        assert_eq!(transform_tag("dup_exon9hyb_star5"), "cn3");
    }

    #[test]
    fn test_empty_input() {
        let params = make_params(6);
        let path = viterbi(&[], &[], &params);
        assert!(path.is_empty());

        let segs = extract_segments(&[], &[], &[]);
        assert!(segs.is_empty());

        let tag = classify_segments(&[], 4, Some(2));
        assert_eq!(tag, None);
    }
}

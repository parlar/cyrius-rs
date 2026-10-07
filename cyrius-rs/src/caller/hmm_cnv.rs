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
//! let sd_per_copy = estimate_sd_per_copy(&raw_d6_cn, &snp_d6, &snp_d7);
//! let params = HmmParams::from_sd_per_copy(total_cn, sd_per_copy);
//! let result = hmm_segment(&raw_d6_cn, &snp_d6, &snp_d7, total_cn, &params);
//! let cnvtag = classify_segments(&result.segments, total_cn, spacer_cn);
//! ```

use std::collections::{BTreeMap, HashMap};

// ---------------------------------------------------------------------------
// CYP2D6 gene structure: SNP site index boundaries
// These match the positions used in cnv_hybrid.rs
// ---------------------------------------------------------------------------
const REP_END: usize = 3;
const EXON9_END: usize = 9;
const INTRON4_BP: usize = 40;
const INTRON1_BP: usize = 74;

// ---------------------------------------------------------------------------
// Emission-model constants
// ---------------------------------------------------------------------------

/// Minimum total (D6 + D7) read depth for a site's observation to be used.
const DEFAULT_MIN_DEPTH: usize = 10;

/// Lower bound on the per-copy emission SD.
///
/// A robust dispersion legitimately collapses to 0 on a flat or heavily tied
/// profile (a whole-gene deletion, for instance), and a zero SD makes the
/// Gaussian emission degenerate. The bound comes from the model's own
/// transition cost rather than from taste: the round-trip penalty for leaving
/// a state and returning is `2 * ln(p_stay / p_switch_each)`, which is 13.4
/// nats at `total_cn = 2` and 14.5 at `total_cn = 5`. For Viterbi not to split
/// a single site off for a deviation below half a CN unit — the midpoint
/// between adjacent CN states, where the observation genuinely becomes
/// ambiguous — the emission SD must satisfy `0.5 * (0.5 / sd)^2 < 13.4`, i.e.
/// `sd > 0.097`. Rounded up, 0.10.
///
/// Measured, not asserted: the floor binds on 22 of the 61 benchmark profiles
/// (the quietest ones, where the SD matters least), and segmentation behaviour
/// over the cohort is identical for any floor from 0.05 to 0.15 — 54/60
/// agreement with `get_cnvtag`, median 2 segments. Only a floor of exactly 0
/// degrades it, and then only via the one all-flat `star5_star5` profile. So
/// this value guards the degenerate case without doing the discriminating
/// work itself.
const MIN_SD_PER_COPY: f64 = 0.10;

/// Gaussian consistency factor converting a median absolute deviation to an SD.
const MAD_TO_SD: f64 = 1.4826;

/// Fewest usable sites a structural region must have for its dispersion to
/// enter the emission-SD estimate.
const MIN_REGION_SITES: usize = 10;

/// The structural regions of the gene, as half-open index ranges into the SNP
/// site array. CYP2D6 copy-number transitions occur at these boundaries, so
/// within one region the true CN is constant. Shared by [`estimate_sd_per_copy`]
/// and [`get_region_cn`] so the two cannot drift apart.
const REGION_BOUNDS: [(&str, usize, usize); 5] = [
    ("rep", 0, REP_END),
    ("exon9", REP_END, EXON9_END),
    ("body", EXON9_END, INTRON4_BP),
    ("intron4", INTRON4_BP, INTRON1_BP),
    ("upstream", INTRON1_BP, usize::MAX), // extends to the end
];

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
    /// Construct default parameters from the total CN call and a per-copy SD.
    ///
    /// `sd_per_copy` is the standard deviation of the per-site D6-CN
    /// observation at copy number 1; the emission SD for state `cn` is
    /// `sd_per_copy * sqrt(cn)`. Take it from [`estimate_sd_per_copy`], which
    /// measures it on the profile being segmented.
    ///
    /// Do **not** pass the depth GMM's `sd` here (review findings #15/#81):
    /// that value describes the normalised depth of the whole D6+D7 region,
    /// not the per-site D6/D7 read split this HMM runs over, and is 4-10x too
    /// small — tight enough that ordinary per-site noise outweighs the
    /// ~14-nat round-trip switch penalty and Viterbi splits off single sites.
    pub fn from_sd_per_copy(total_cn: u32, sd_per_copy: f64) -> Self {
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
            min_depth: DEFAULT_MIN_DEPTH,
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

/// Robust (MAD-based) dispersion of a set of observations, scaled to a
/// Gaussian SD. Returns 0.0 for fewer than two values.
fn robust_sd(values: &[f64]) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let level = crate::stats::median(values);
    let deviations: Vec<f64> = values.iter().map(|x| (x - level).abs()).collect();
    MAD_TO_SD * crate::stats::median(&deviations)
}

/// Estimate the per-copy emission SD from the per-site profile itself.
///
/// The emission model compares each site's observation to an integer CN, so
/// the SD it needs is the *marginal* spread of the observations about their
/// segment level. The depth GMM cannot supply that: its `sd` describes the
/// normalised depth of the whole D6+D7 region, and using it makes the
/// emissions 4-10x too tight (review findings #15/#81).
///
/// Method: take the robust dispersion of the usable observations **within
/// each structural region** of [`REGION_BOUNDS`], and use the median across
/// regions. CYP2D6 copy-number transitions occur at those boundaries, so a
/// within-region spread is observation noise rather than CN signal, and the
/// median across regions tolerates one region holding an off-boundary
/// breakpoint. The result is divided by `sqrt(median observation)` to undo the
/// `sqrt(cn)` scaling that [`HmmParams::from_sd_per_copy`] re-applies, and
/// floored at [`MIN_SD_PER_COPY`].
///
/// A whole-profile dispersion is the fallback when no region has
/// [`MIN_REGION_SITES`] usable sites.
///
/// Measured over the 61 benchmark profiles that have a d67 profile, before
/// flooring (per copy): q1 0.092, median 0.112, q3 0.136, max 0.238 — against
/// the single value 0.040 the GMM route produced. Replacing the GMM SD cuts the
/// median segment count from 6 to 2 and the worst case from 19 to 7, drops
/// profiles containing a single-site segment from 30 to 6, removes both of the
/// exact region-vote ties (see [`region_cn_vote`]), and raises agreement with
/// `get_cnvtag` from 51/60 to 54/60 with one disagreement instead of two.
///
/// Note that a global dispersion must *not* be used in place of the
/// per-region one: real CN steps inflate it (0.26-0.64 on the `star68` and
/// `star5_star68` samples), which suppresses the very breakpoints the HMM
/// exists to find.
///
/// A lag-1 first-difference estimator must not be used either. Adjacent
/// diagnostic sites are covered by the same reads, so their D6/D7 split noise
/// is shared: the lag-1 innovation SD is 0.053 against a marginal spread of
/// 0.237 on the cn2 samples, and the estimate lands on the floor for all 61
/// profiles with a d67 profile — a constant wearing an estimator's clothes.
pub fn estimate_sd_per_copy(raw_d6_cn: &[f64], snp_d6: &[usize], snp_d7: &[usize]) -> f64 {
    let n = raw_d6_cn.len().min(snp_d6.len()).min(snp_d7.len());
    let site_ok =
        |i: usize| raw_d6_cn[i].is_finite() && snp_d6[i] + snp_d7[i] >= DEFAULT_MIN_DEPTH;

    let mut region_sds: Vec<f64> = Vec::new();
    for &(_, region_start, region_end) in REGION_BOUNDS.iter() {
        let hi = region_end.min(n);
        if region_start >= hi {
            continue;
        }
        let values: Vec<f64> = (region_start..hi)
            .filter(|&i| site_ok(i))
            .map(|i| raw_d6_cn[i])
            .collect();
        if values.len() >= MIN_REGION_SITES {
            region_sds.push(robust_sd(&values));
        }
    }

    let all: Vec<f64> = (0..n).filter(|&i| site_ok(i)).map(|i| raw_d6_cn[i]).collect();
    let sd_site = if region_sds.is_empty() {
        robust_sd(&all)
    } else {
        crate::stats::median(&region_sds)
    };
    let ref_cn = crate::stats::median(&all).max(1.0);

    (sd_site / ref_cn.sqrt()).max(MIN_SD_PER_COPY)
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
    let mut result = HashMap::new();

    for &(name, region_start, region_end) in REGION_BOUNDS.iter() {
        if let Some(cn) = region_cn_vote(segments, region_start, region_end) {
            result.insert(name.to_string(), cn);
        }
    }

    result
}

/// Overlap-weighted vote for the dominant CN of one region, `[start, end)`.
///
/// Each segment contributes the number of sites it shares with the region.
/// Votes accumulate in a `BTreeMap`, so the outcome does not depend on hash
/// seeding, and an exact tie between two or more CN states returns `None`:
/// a tie carries no information about the region's CN, and resolving it
/// arbitrarily used to pick between structurally different CNV tags at random
/// (review finding #16). `None` is also returned when no segment overlaps.
///
/// A missing region leaves `get_region_cn` without that key, which makes
/// `classify_segments` return `None` — no CNV tag rather than a coin flip.
fn region_cn_vote(
    segments: &[CnSegment],
    region_start: usize,
    region_end: usize,
) -> Option<u32> {
    let mut votes: BTreeMap<u32, usize> = BTreeMap::new();
    for seg in segments {
        let overlap_start = seg.start_idx.max(region_start);
        let overlap_end = (seg.end_idx + 1).min(region_end);
        if overlap_start < overlap_end {
            *votes.entry(seg.cn).or_insert(0) += overlap_end - overlap_start;
        }
    }

    let best = votes.values().copied().max()?;
    let mut winners = votes.iter().filter(|&(_, &count)| count == best);
    let (&cn, _) = winners.next()?;
    if winners.next().is_some() {
        return None; // exact tie — region CN unknown
    }
    Some(cn)
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
        let params = HmmParams::from_sd_per_copy(total_cn, 0.3);

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

    // -----------------------------------------------------------------------
    // Review finding #16: the region vote must be deterministic, and an exact
    // tie must report an unknown region CN rather than a hash-order winner.
    // -----------------------------------------------------------------------

    /// The reviewer's reproduction case: segments [0..5 cn2, 6..116 cn3] split
    /// the 6-site exon9 region exactly 3/3, so that region's vote ties.
    fn tie_segments() -> Vec<CnSegment> {
        vec![
            CnSegment { start_idx: 0, end_idx: 5, cn: 2, n_valid: 6, mean_obs: 2.0 },
            CnSegment { start_idx: 6, end_idx: 116, cn: 3, n_valid: 111, mean_obs: 3.0 },
        ]
    }

    #[test]
    fn test_classify_segments_tie_is_deterministic_and_none() {
        let segments = tie_segments();
        let results: Vec<Option<String>> = (0..200)
            .map(|_| classify_segments(&segments, 5, Some(2)))
            .collect();
        let mut distinct: Vec<&Option<String>> = Vec::new();
        for r in &results {
            if !distinct.contains(&r) {
                distinct.push(r);
            }
        }
        assert_eq!(
            distinct.len(), 1,
            "classify_segments is non-deterministic on a tie: saw {:?} over 200 calls",
            distinct
        );
        assert_eq!(results[0], None,
            "an exact tie in the exon9 region must leave the region CN unknown");
    }

    #[test]
    fn test_region_cn_vote_tie_returns_none() {
        let segments = tie_segments();
        // exon9 (sites 3..9) is split 3/3 between cn2 and cn3 -> unknown
        assert_eq!(region_cn_vote(&segments, REP_END, EXON9_END), None);
        // rep (sites 0..3) is entirely cn2 -> unambiguous
        assert_eq!(region_cn_vote(&segments, 0, REP_END), Some(2));
        // body (sites 9..40) is entirely cn3 -> unambiguous
        assert_eq!(region_cn_vote(&segments, EXON9_END, INTRON4_BP), Some(3));
        // a region no segment overlaps is unknown
        assert_eq!(region_cn_vote(&segments, 200, 300), None);
        // a one-site majority is still a winner, not a tie
        assert_eq!(region_cn_vote(&segments, 0, 5), Some(2));
    }

    // -----------------------------------------------------------------------
    // Review findings #15 / #81: the emission SD must be estimated from the
    // per-site profile, not taken from the depth GMM, or Viterbi over-segments.
    // -----------------------------------------------------------------------

    /// Deterministic synthetic profile with a single true D6 copy number,
    /// generated from the real noise source: a binomial D6/D7 read split at
    /// per-site depths sweeping 20..80 (so the per-site SD, ~0.3, and its tails
    /// match the 0.19-0.34 measured on the cn2 benchmark samples).
    ///
    /// Returns (raw_d6_cn, snp_d6, snp_d7).
    fn binomial_flat_profile(
        n: usize,
        total_cn: u32,
        true_d6_cn: u32,
    ) -> (Vec<f64>, Vec<usize>, Vec<usize>) {
        let mut state: u64 = 0x2d6_cafe;
        let mut next = || {
            // 32-bit LCG (Numerical Recipes constants) — reproducible noise
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223) & 0xffff_ffff;
            (state >> 8) as f64 / ((1u64 << 24) as f64)
        };
        let p = true_d6_cn as f64 / total_cn as f64;
        let mut raw = Vec::with_capacity(n);
        let mut d6 = Vec::with_capacity(n);
        let mut d7 = Vec::with_capacity(n);
        for i in 0..n {
            let depth = 20 + (i * 7) % 61;
            let k = (0..depth).filter(|_| next() < p).count();
            raw.push(total_cn as f64 * k as f64 / depth as f64);
            d6.push(k);
            d7.push(depth - k);
        }
        (raw, d6, d7)
    }

    #[test]
    fn test_noisy_flat_profile_segments_into_one_segment() {
        // A cn2 sample: total D6+D7 CN = 4, per-site D6 CN ~ 2.
        let total_cn = 4u32;
        let (raw, d6, d7) = binomial_flat_profile(117, total_cn, 2);

        let sd_per_copy = estimate_sd_per_copy(&raw, &d6, &d7);
        let params = HmmParams::from_sd_per_copy(total_cn, sd_per_copy);
        let result = hmm_segment(&raw, &d6, &d7, total_cn, &params);
        assert_eq!(
            result.segments.len(), 1,
            "one true CN must give one segment (sd_per_copy={sd_per_copy:.3}), got {:?}",
            result.segments.iter().map(|s| (s.start_idx, s.end_idx, s.cn)).collect::<Vec<_>>()
        );
        assert_eq!(result.segments[0].cn, 2);
        assert_eq!(
            classify_segments(&result.segments, total_cn, Some(2)),
            Some("cn2".to_string())
        );
    }

    #[test]
    fn test_gmm_derived_sd_over_segments_the_same_profile() {
        // Regression guard for findings #15/#81: data/CYP2D6_gmm.txt carries
        // `d67 sd 2:0.056` and main.rs used 0.056 / sqrt(2) = 0.0396, which
        // shattered this single-CN profile.
        let total_cn = 4u32;
        let (raw, d6, d7) = binomial_flat_profile(117, total_cn, 2);
        let params =
            HmmParams::from_sd_per_copy(total_cn, 0.056 / std::f64::consts::SQRT_2);
        let result = hmm_segment(&raw, &d6, &d7, total_cn, &params);
        assert!(
            result.segments.len() > 5,
            "the GMM-derived SD is expected to over-segment this profile, got {} segments",
            result.segments.len()
        );
    }

    #[test]
    fn test_estimate_sd_per_copy_tracks_the_profile_and_ignores_cn_steps() {
        let total_cn = 4u32;
        let (raw, d6, d7) = binomial_flat_profile(117, total_cn, 2);
        let sd_per_copy = estimate_sd_per_copy(&raw, &d6, &d7);
        // Per-site SD of this profile is ~0.31 at CN 2, so per copy ~0.22.
        assert!(
            (0.15..0.30).contains(&sd_per_copy),
            "expected sd_per_copy near 0.22, got {sd_per_copy}"
        );
        // A real CN step raises the dispersion of the one region it falls in
        // (here intron4, sites 40..74); the median across regions absorbs it.
        let mut stepped = raw.clone();
        for x in stepped.iter_mut().take(55) {
            *x += 1.0;
        }
        let sd_stepped = estimate_sd_per_copy(&stepped, &d6, &d7);
        assert!(
            (sd_stepped - sd_per_copy).abs() < 0.03,
            "a CN step moved the estimate from {sd_per_copy} to {sd_stepped}"
        );
        // A global dispersion would not survive that step — the guard that
        // keeps the HMM able to see breakpoints at all.
        assert!(
            robust_sd(&stepped) > 2.0 * robust_sd(&raw),
            "expected the global dispersion to be inflated by the step"
        );
    }

    #[test]
    fn test_estimate_sd_per_copy_floors_on_degenerate_profile() {
        // A perfectly flat profile has a median absolute first difference of 0;
        // the estimator must floor rather than return a degenerate SD of 0.
        let raw = vec![2.0; 50];
        let d6 = vec![30usize; 50];
        let d7 = vec![30usize; 50];
        assert_eq!(estimate_sd_per_copy(&raw, &d6, &d7), MIN_SD_PER_COPY);

        // All sites below min depth -> no usable differences -> floor.
        let low = vec![1usize; 50];
        assert_eq!(estimate_sd_per_copy(&raw, &low, &low), MIN_SD_PER_COPY);

        // Too short to difference -> floor.
        assert_eq!(estimate_sd_per_copy(&[2.0], &[30], &[30]), MIN_SD_PER_COPY);

        // The floored SD still gives a proper (non-degenerate) emission density.
        let params = HmmParams::from_sd_per_copy(4, MIN_SD_PER_COPY);
        assert!(params.emission_sd.iter().all(|sd| *sd > 0.0 && sd.is_finite()));
    }
}

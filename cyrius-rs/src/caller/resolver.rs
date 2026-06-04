//! Likelihood-based tie-breaker for ;-ambiguous diplotype candidate sets.
//! A candidate joins the reported set iff its posterior is within RESOLVE_MARGIN
//! (near-tie epsilon) of the top candidate's posterior. The call is confident
//! iff the top candidate stands alone — no rival within that gap.

use crate::caller::diplotype_caller::CandidateScore;

/// A candidate joins the reported set iff its posterior is within this gap of
/// the top candidate's posterior (a near-tie epsilon). The call is confident
/// iff the top candidate stands alone (no rival within this gap). Calibrated on
/// the GeT-RM benchmark ONLY; conservative — prefer reporting ambiguity over a
/// confident wrong pick. (Starting value; re-confirm during benchmark calibration.)
pub const RESOLVE_MARGIN: f64 = 0.10;

#[derive(Debug, Clone)]
pub struct ResolverOutcome {
    /// Single call, or ";"-joined near-tied set.
    pub genotype: String,
    pub ambiguous: bool,
    pub top_posterior: f64,
}

/// Decide a genotype from scored candidates. Returns None when nothing was
/// scorable (every score non-finite) — the caller should then keep the
/// cascade's own choice.
pub fn resolve(scores: &[CandidateScore]) -> Option<ResolverOutcome> {
    let mut finite: Vec<&CandidateScore> =
        scores.iter().filter(|c| c.log10_likelihood.is_finite()).collect();
    if finite.is_empty() {
        return None;
    }
    finite.sort_by(|a, b| {
        b.posterior
            .total_cmp(&a.posterior)
            .then_with(|| a.diplotype.cmp(&b.diplotype))
    });
    let top = finite[0];
    let kept: Vec<String> = finite
        .iter()
        .filter(|c| top.posterior - c.posterior <= RESOLVE_MARGIN)
        .map(|c| c.diplotype.clone())
        .collect();
    if kept.len() == 1 {
        Some(ResolverOutcome {
            genotype: top.diplotype.clone(),
            ambiguous: false,
            top_posterior: top.posterior,
        })
    } else {
        Some(ResolverOutcome {
            genotype: kept.join(";"),
            ambiguous: true,
            top_posterior: top.posterior,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cs(name: &str, post: f64) -> CandidateScore {
        CandidateScore { diplotype: name.into(), log10_likelihood: post.ln(), posterior: post }
    }

    #[test]
    fn single_candidate_is_confident() {
        let out = resolve(&[cs("*1/*2", 1.0)]).unwrap();
        assert_eq!(out.genotype, "*1/*2");
        assert!(!out.ambiguous);
    }

    #[test]
    fn clear_winner_picks_top() {
        let out = resolve(&[cs("*1/*2", 0.98), cs("*1/*3", 0.02)]).unwrap();
        assert_eq!(out.genotype, "*1/*2");
        assert!(!out.ambiguous);
    }

    #[test]
    fn near_tie_reports_ambiguous_set() {
        let out = resolve(&[cs("*36/*36+*10", 0.5), cs("*10/*36x2", 0.5)]).unwrap();
        assert!(out.ambiguous);
        assert!(out.genotype.contains(';'));
        assert!(out.genotype.contains("*36/*36+*10") && out.genotype.contains("*10/*36x2"));
    }

    #[test]
    fn all_unscorable_returns_none() {
        assert!(resolve(&[cs("*1/*2", 0.0)]).is_none());
    }

    #[test]
    fn dominant_top_excludes_tail() {
        // 0.89 dominates; 0.09 and 0.02 are NOT within 0.10 of the top -> confident single.
        let out = resolve(&[cs("*1/*2", 0.89), cs("*1/*3", 0.09), cs("*1/*4", 0.02)]).unwrap();
        assert_eq!(out.genotype, "*1/*2");
        assert!(!out.ambiguous);
    }

    #[test]
    fn near_tie_pair_excludes_distant_third() {
        // top two within 0.10 of each other are kept; the 0.01 tail is excluded.
        let out = resolve(&[cs("*1/*2", 0.50), cs("*1/*3", 0.49), cs("*1/*4", 0.01)]).unwrap();
        assert!(out.ambiguous);
        assert!(out.genotype.contains("*1/*2") && out.genotype.contains("*1/*3"));
        assert!(!out.genotype.contains("*1/*4"));
    }
}

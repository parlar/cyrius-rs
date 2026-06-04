//! Likelihood-based tie-breaker for ;-ambiguous diplotype candidate sets.
//! Picks a single call when one candidate's posterior clearly dominates,
//! otherwise reports the near-tied set as honest ambiguity.

use crate::caller::diplotype_caller::CandidateScore;

/// Minimum posterior gap between the top two candidates required for a
/// confident single call. Calibrated on the GeT-RM benchmark ONLY and set
/// conservatively: prefer declaring ambiguity over a confident wrong pick.
/// (Starting value; re-confirm during the benchmark calibration step.)
pub const RESOLVE_MARGIN: f64 = 0.90;

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
    if finite.len() == 1 {
        return Some(ResolverOutcome {
            genotype: top.diplotype.clone(),
            ambiguous: false,
            top_posterior: top.posterior,
        });
    }
    let second = finite[1];
    if top.posterior - second.posterior >= RESOLVE_MARGIN {
        Some(ResolverOutcome {
            genotype: top.diplotype.clone(),
            ambiguous: false,
            top_posterior: top.posterior,
        })
    } else {
        let kept: Vec<String> = finite
            .iter()
            .filter(|c| top.posterior - c.posterior < RESOLVE_MARGIN)
            .map(|c| c.diplotype.clone())
            .collect();
        Some(ResolverOutcome {
            ambiguous: kept.len() > 1,
            genotype: kept.join(";"),
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
        let mut a = cs("*1/*2", 0.0);
        a.log10_likelihood = f64::NEG_INFINITY;
        assert!(resolve(&[a]).is_none());
    }
}

//! Likelihood-based tie-breaker for ;-ambiguous diplotype candidate sets.
//! The top candidate is picked only if its posterior leads every rival by more
//! than RESOLVE_MARGIN AND enough individual reads clearly favour it (the same
//! evidence rule phase_disambiguate uses). Otherwise every candidate is kept:
//! the per-read-independent posterior is too overconfident to narrow the set
//! on its own.

use crate::caller::diplotype_caller::CandidateScore;

/// The top candidate's posterior must lead every rival by more than this.
pub const RESOLVE_MARGIN: f64 = 0.10;

/// Minimum number of decisive reads (reads that clearly favour one candidate)
/// the top candidate needs.
pub const MIN_DECISIVE_READS: usize = 3;

/// The top candidate needs at least this many times as many decisive reads
/// as any rival.
pub const DECISIVE_RATIO: usize = 3;

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
    let rivals = &finite[1..];
    let leads = rivals.iter().all(|c| top.posterior - c.posterior > RESOLVE_MARGIN);
    let rival_reads = rivals.iter().map(|c| c.decisive_reads).max().unwrap_or(0);
    let supported = top.decisive_reads >= MIN_DECISIVE_READS
        && top.decisive_reads >= DECISIVE_RATIO * rival_reads;
    if rivals.is_empty() || (leads && supported) {
        Some(ResolverOutcome {
            genotype: top.diplotype.clone(),
            ambiguous: false,
            top_posterior: top.posterior,
        })
    } else {
        let kept: Vec<String> = finite.iter().map(|c| c.diplotype.clone()).collect();
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

    /// Candidate whose decisive-read count follows its posterior (plenty of
    /// evidence for a clear winner).
    fn cs(name: &str, post: f64) -> CandidateScore {
        csd(name, post, (post * 100.0) as usize)
    }

    fn csd(name: &str, post: f64, decisive_reads: usize) -> CandidateScore {
        CandidateScore {
            diplotype: name.into(),
            log10_likelihood: post.ln(),
            posterior: post,
            decisive_reads,
        }
    }

    #[test]
    fn too_few_decisive_reads_keeps_every_candidate() {
        let out = resolve(&[csd("*1/*21", 0.999, 2), csd("*2/*21", 0.001, 0)]).unwrap();
        assert!(out.ambiguous);
        assert!(out.genotype.contains("*1/*21") && out.genotype.contains("*2/*21"));
    }

    #[test]
    fn lead_under_three_to_one_keeps_every_candidate() {
        let out = resolve(&[csd("*1/*21", 0.99, 5), csd("*2/*21", 0.01, 2)]).unwrap();
        assert!(out.ambiguous);
    }

    #[test]
    fn enough_decisive_reads_picks_top() {
        let out = resolve(&[csd("*1/*21", 0.99, 6), csd("*2/*21", 0.01, 2)]).unwrap();
        assert_eq!(out.genotype, "*1/*21");
        assert!(!out.ambiguous);
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
    fn near_tie_keeps_every_candidate() {
        // Without a decided winner nothing is narrowed: the posterior alone is
        // too overconfident to drop the 0.01 tail.
        let out = resolve(&[cs("*1/*2", 0.50), cs("*1/*3", 0.49), cs("*1/*4", 0.01)]).unwrap();
        assert!(out.ambiguous);
        assert!(out.genotype.contains("*1/*2") && out.genotype.contains("*1/*3"));
        assert!(out.genotype.contains("*1/*4"));
    }
}

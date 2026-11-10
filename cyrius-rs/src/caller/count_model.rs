//! Per-site count model for the candidate ranking: Stage A's chosen model M4
//! (exp-H20 result 0a4b020). Beta-binomial with a per-class error rate, allele
//! bias and overdispersion, fitted on 169 cn2 samples (144 HPRC + 25 GeT-RM).
//! The parameters are the `all_cn2` rows of exp-H20 `results/params.tsv`.

use statrs::function::gamma::ln_gamma;

/// Reported as `Candidate_ranking.Model`.
pub const MODEL_ID: &str = "M4 (exp-H20 0a4b020)";

/// Site class: the `variant_type` column, with the homology file as its own class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiteClass {
    Snp,
    SnpMulti,
    Indel,
    IndelMulti,
    GeneConversion,
    Homology,
}

impl SiteClass {
    pub fn from_variant_type(variant_type: &str, homology: bool) -> Option<SiteClass> {
        if homology {
            return Some(SiteClass::Homology);
        }
        match variant_type {
            "snp" => Some(SiteClass::Snp),
            "snp_multi" => Some(SiteClass::SnpMulti),
            "indel" => Some(SiteClass::Indel),
            "indel_multi" => Some(SiteClass::IndelMulti),
            "gene_conversion" => Some(SiteClass::GeneConversion),
            _ => None,
        }
    }

    /// (error rate, allele bias, overdispersion): exp-H20 `results/params.tsv`, M4, all_cn2.
    pub fn params(self) -> (f64, f64, f64) {
        match self {
            SiteClass::Snp => (0.0004015441770681196, 0.9339356268269258, 0.023072497450790153),
            SiteClass::SnpMulti => (0.0005162070342784222, 0.980262902642141, 0.0027068788553678198),
            SiteClass::Indel => (0.0016441425094134746, 0.754160272105761, 0.31571547734888167),
            SiteClass::IndelMulti => (0.0001, 0.8434766740495846, 0.11761644655153744),
            SiteClass::GeneConversion => (0.0002807889779143907, 0.9301523446404545, 0.0),
            SiteClass::Homology => (0.0006937702284259667, 1.0072875994763135, 0.0033292053680482486),
        }
    }
}

fn ln_choose(n: f64, a: f64) -> f64 {
    ln_gamma(n + 1.0) - ln_gamma(a + 1.0) - ln_gamma(n - a + 1.0)
}

/// Natural-log probability of `alt` alternate reads among `alt + ref_` at a site with `n`
/// copies, `k` of them carrying the variant. Same formula as exp-H20 `models.ll_vec("M4")`:
/// p = (k/n)(1 - eps) + (1 - k/n) eps, biased by beta, beta-binomial with rho (binomial when
/// rho < 1e-6).
pub fn site_ln_likelihood(class: SiteClass, alt: u32, ref_: u32, k: u32, n: u32) -> f64 {
    ln_likelihood_with(class.params(), alt, ref_, k, n)
}

/// `site_ln_likelihood` with explicit (error rate, allele bias, overdispersion).
fn ln_likelihood_with((eps, beta, rho): (f64, f64, f64), alt: u32, ref_: u32, k: u32, n: u32) -> f64 {
    let (a, r) = (alt as f64, ref_ as f64);
    let f = k as f64 / n as f64;
    let p = f * (1.0 - eps) + (1.0 - f) * eps;
    let p = (beta * p / (beta * p + 1.0 - p)).clamp(1e-12, 1.0 - 1e-12);
    if rho < 1e-6 {
        ln_choose(a + r, a) + a * p.ln() + r * (1.0 - p).ln()
    } else {
        let s = (1.0 - rho) / rho;
        let (al, be) = (p * s, (1.0 - p) * s);
        ln_choose(a + r, a) + ln_gamma(a + al) + ln_gamma(r + be) - ln_gamma(a + r + al + be)
            - ln_gamma(al) - ln_gamma(be) + ln_gamma(al + be)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference values from exp-H20 `scripts/models.py` (`ll_vec("M4", ...)`) with the
    /// `all_cn2` parameters, printed by Python before this plan was written.
    #[test]
    fn matches_the_stage_a_python_model() {
        let cases = [
            (SiteClass::Snp, 20, 20, 1, 2, -2.424344032803589),
            (SiteClass::Snp, 0, 30, 0, 2, -0.008584880431968145),
            (SiteClass::Indel, 7, 33, 1, 2, -3.5155621349832145),
            (SiteClass::Homology, 10, 50, 1, 4, -3.1241100424440447),
            (SiteClass::GeneConversion, 15, 15, 1, 2, -1.954377886442014),
            (SiteClass::Snp, 3, 40, 0, 3, -7.284606459391441),
            (SiteClass::IndelMulti, 0, 25, 0, 1, -0.000961107793626681),
            (SiteClass::SnpMulti, 12, 0, 2, 2, -0.006228342306030754),
        ];
        for (class, a, r, k, n, want) in cases {
            let got = site_ln_likelihood(class, a, r, k, n);
            assert!((got - want).abs() < 1e-9, "{class:?} {a}/{r} k={k} n={n}: {got} vs {want}");
        }
    }

    #[test]
    fn is_a_distribution_over_the_alt_count() {
        use SiteClass::*;
        for class in [Snp, SnpMulti, Indel, IndelMulti, GeneConversion, Homology] {
            let total: f64 = (0..=30).map(|a| site_ln_likelihood(class, a, 30 - a, 1, 2).exp()).sum();
            assert!((total - 1.0).abs() < 1e-9, "{class:?}: {total}");
        }
    }

    /// Spec "Tests": with no allele bias (beta = 1), swapping the alt and ref reads and
    /// replacing k by n - k gives the same likelihood, for the binomial (rho = 0) and the
    /// beta-binomial (rho > 0). With the fitted SNP bias (beta 0.934) it does not.
    #[test]
    fn binomial_model_is_symmetric_in_ref_and_alt() {
        for rho in [0.0, 0.023, 0.316] {
            for n in 1..=4u32 {
                for k in 0..=n {
                    for (a, r) in [(0, 30), (3, 40), (20, 20), (7, 33), (12, 0)] {
                        let p = (0.0004, 1.0, rho);
                        let (x, y) = (ln_likelihood_with(p, a, r, k, n), ln_likelihood_with(p, r, a, n - k, n));
                        assert!((x - y).abs() < 1e-9, "rho={rho} n={n} k={k} {a}/{r}: {x} vs {y}");
                    }
                }
            }
        }
        let biased = |a, r, k| site_ln_likelihood(SiteClass::Snp, a, r, k, 2);
        assert!((biased(7, 33, 1) - biased(33, 7, 1)).abs() > 1e-3);
    }

    #[test]
    fn variant_types_map_to_classes() {
        assert_eq!(SiteClass::from_variant_type("snp", false), Some(SiteClass::Snp));
        assert_eq!(SiteClass::from_variant_type("snp_multi", false), Some(SiteClass::SnpMulti));
        assert_eq!(SiteClass::from_variant_type("indel", false), Some(SiteClass::Indel));
        assert_eq!(SiteClass::from_variant_type("indel_multi", false), Some(SiteClass::IndelMulti));
        assert_eq!(SiteClass::from_variant_type("gene_conversion", false), Some(SiteClass::GeneConversion));
        assert_eq!(SiteClass::from_variant_type("snp", true), Some(SiteClass::Homology));
        assert_eq!(SiteClass::from_variant_type("weird", false), None);
    }
}

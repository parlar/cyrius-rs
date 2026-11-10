use crate::stats::python_round3;
use statrs::distribution::{Discrete, Poisson};

const POSTERIOR_CUTOFF_STRINGENT: f64 = 0.9;
const ERROR_RATE: f64 = 0.1;

/// Return the reg1 copy number call at each site based on Poisson likelihood,
/// with a minimum read support cutoff.
/// Returns a Vec of 1 element (confident call) or 4 elements (two most likely scenarios).
pub fn call_reg1_cn(full_cn: Option<u32>, count_reg1: f64, count_reg2: f64, min_read: f64) -> Vec<Option<CnProb>> {
    match full_cn {
        None => vec![None],
        Some(0) => vec![Some(CnProb::Single(0))],
        Some(1) if count_reg1 > min_read => vec![Some(CnProb::Single(1))],
        Some(cn) => {
            let Some(post_prob) = cn_posterior(cn, count_reg1, count_reg2) else {
                return vec![None];
            };

            let mut post_prob_sorted: Vec<(usize, f64)> =
                post_prob.iter().copied().enumerate().collect();
            post_prob_sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

            let top_idx = post_prob_sorted[0].0;
            let top_prob = post_prob_sorted[0].1;

            if top_idx != 0 && count_reg1 <= min_read && count_reg2 >= min_read {
                return vec![Some(CnProb::Single(0))];
            }

            if top_prob >= POSTERIOR_CUTOFF_STRINGENT {
                return vec![Some(CnProb::Single(top_idx as u32))];
            }

            // Output the two most likely scenarios
            vec![Some(CnProb::TwoOptions {
                cn1: top_idx as u32,
                prob1: python_round3(top_prob),
                cn2: post_prob_sorted[1].0 as u32,
                prob2: python_round3(post_prob_sorted[1].1),
            })]
        }
    }
}

/// The normalised posterior over 0..=cn that `call_reg1_cn` decides from: the same expected
/// depths, the same smaller-count Poisson term and the same reversal. `None` where
/// `call_reg1_cn` computes no Poisson term. Used only to report decision provenance.
pub fn cn_posterior(cn: u32, count_reg1: f64, count_reg2: f64) -> Option<Vec<f64>> {
    let nsum = count_reg1 + count_reg2;
    if cn == 0 || nsum == 0.0 {
        return None;
    }
    let count = if count_reg1 <= count_reg2 { count_reg1 as u64 } else { count_reg2 as u64 };
    let mut prob = Vec::with_capacity(cn as usize + 1);
    for i in 0..=cn {
        let depth_expected = if i == 0 {
            (ERROR_RATE / 3.0) * nsum
        } else if i == cn {
            nsum - ERROR_RATE * nsum
        } else {
            nsum * i as f64 / cn as f64
        };
        prob.push(Poisson::new(depth_expected).unwrap().pmf(count));
    }
    let sum: f64 = prob.iter().sum();
    if sum == 0.0 {
        return None;
    }
    let mut post: Vec<f64> = prob.iter().map(|p| p / sum).collect();
    if count_reg2 < count_reg1 {
        post.reverse();
    }
    Some(post)
}

/// Represents a CN probability call — either confident (single) or ambiguous (two options).
#[derive(Debug, Clone, Copy)]
pub enum CnProb {
    Single(u32),
    TwoOptions {
        cn1: u32,
        prob1: f64,
        cn2: u32,
        prob2: f64,
    },
}

/// Wrapper matching Python's list return (1-element or 4-element).
#[derive(Debug, Clone)]
pub struct CnProbResult(pub Option<CnProb>);

impl From<Option<CnProb>> for CnProbResult {
    fn from(val: Option<CnProb>) -> Self {
        CnProbResult(val)
    }
}

/// Filter raw CN calls based on posterior probability cutoff.
/// For gene conversion cases (SNVs between paralogs).
pub fn process_raw_call_gc(
    cn_prob: &[CnProbResult],
    post_cutoff: f64,
    keep_none: bool,
) -> Vec<Option<u32>> {
    let mut filtered = Vec::new();
    for cn_call in cn_prob {
        match cn_call.0 {
            None => {
                if keep_none {
                    filtered.push(None);
                }
            }
            Some(CnProb::Single(val)) => {
                filtered.push(Some(val));
            }
            Some(CnProb::TwoOptions { cn1, prob1, .. }) => {
                if prob1 > post_cutoff {
                    filtered.push(Some(cn1));
                } else if keep_none {
                    filtered.push(None);
                }
            }
        }
    }
    filtered
}

/// Filter raw CN calls based on posterior probability cutoff.
/// For de novo variant calling (non-gene-conversion cases).
pub fn process_raw_call_denovo(
    cn_prob: &[CnProbResult],
    post_cutoff1: f64,
    post_cutoff2: f64,
    list_total_cn: Option<&[u32]>,
    keep_none: bool,
) -> Vec<Option<u32>> {
    let mut filtered = Vec::new();
    for (i, cn_call) in cn_prob.iter().enumerate() {
        match cn_call.0 {
            None => {
                if keep_none {
                    filtered.push(None);
                }
            }
            Some(CnProb::Single(val)) => {
                filtered.push(Some(val));
            }
            Some(CnProb::TwoOptions {
                cn1,
                prob1,
                cn2,
                prob2: _,
            }) => {
                let keep_var = if let Some(total_cn) = list_total_cn {
                    let tcn = total_cn[i];
                    (cn1 > 0 && cn2 > 0) || (cn1 == tcn || cn2 == tcn)
                } else {
                    cn1 > 0 && cn2 > 0
                };

                if prob1 > post_cutoff1 {
                    filtered.push(Some(cn1));
                } else if keep_var {
                    if prob1 > post_cutoff2 {
                        filtered.push(Some(cn1));
                    } else {
                        filtered.push(Some(std::cmp::min(cn1, cn2)));
                    }
                } else if keep_none {
                    filtered.push(None);
                }
            }
        }
    }
    filtered
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Checks that the decision layer of call_reg1_cn is consistent with the shared posterior:
    /// wherever call_reg1_cn returns two options or a stringent single value, they must agree
    /// with cn_posterior. Both sides now use the same code, so this is NOT an independent check
    /// of the posterior itself. Orientation is pinned by cn_posterior_is_indexed_by_region_1_copies,
    /// and pre-refactor behaviour by call_reg1_cn_matches_its_pre_refactor_outputs.
    /// min_read = -1 keeps the min-read rule from firing; cn >= 2 keeps the one-copy shortcut out.
    #[test]
    fn cn_posterior_is_the_posterior_call_reg1_cn_decides_from() {
        let mut checked = 0;
        for cn in 2..=4u32 {
            for c1 in 0..30 {
                for c2 in 0..30 {
                    let (c1, c2) = (c1 as f64, c2 as f64);
                    let post = cn_posterior(cn, c1, c2);
                    match call_reg1_cn(Some(cn), c1, c2, -1.0)[0] {
                        None => assert!(post.is_none()),
                        Some(CnProb::Single(v)) => {
                            let p = post.unwrap();
                            let top = (0..p.len()).fold(0, |b, i| if p[i] > p[b] { i } else { b });
                            assert_eq!(top as u32, v, "cn={cn} {c1}/{c2}");
                            checked += 1;
                        }
                        Some(CnProb::TwoOptions { cn1, prob1, .. }) => {
                            let p = post.unwrap();
                            assert!((p[cn1 as usize] - prob1).abs() <= 0.0005 + 1e-12, "cn={cn} {c1}/{c2}");
                            checked += 1;
                        }
                    }
                }
            }
        }
        assert!(checked > 2000);
    }

    /// The test above can no longer catch a broken reversal: call_reg1_cn now decides from
    /// cn_posterior, so both sides move together. This one pins the orientation independently:
    /// the posterior is indexed by copies of region 1, so swapping the counts reverses it.
    #[test]
    fn cn_posterior_is_indexed_by_region_1_copies() {
        let p = cn_posterior(2, 20.0, 0.0).unwrap();
        assert!(p[2] > 0.9, "{p:?}"); // all reads region 1 -> both copies
        let p = cn_posterior(2, 0.0, 20.0).unwrap();
        assert!(p[0] > 0.9, "{p:?}");
        for (a, b) in [(3.0, 11.0), (14.0, 5.0), (9.0, 2.0)] {
            let mut fwd = cn_posterior(3, a, b).unwrap();
            fwd.reverse();
            assert_eq!(fwd, cn_posterior(3, b, a).unwrap());
        }
    }

    /// Golden outputs of call_reg1_cn(Some(cn), c1, c2, min_read)[0], printed by the code at
    /// commit da59f5c (before call_reg1_cn was refactored onto cn_posterior).
    #[test]
    fn call_reg1_cn_matches_its_pre_refactor_outputs() {
        fn two(cn1: u32, prob1: f64, cn2: u32, prob2: f64) -> Option<CnProb> {
            Some(CnProb::TwoOptions { cn1, prob1, cn2, prob2 })
        }
        fn one(v: u32) -> Option<CnProb> {
            Some(CnProb::Single(v))
        }
        let golden = vec![
            ((2, 5.0, 1.0, 4.0), two(2, 0.485, 1, 0.443)),
            ((2, 20.0, 20.0, 4.0), one(1)),
            ((2, 21.0, 19.0, 4.0), one(1)),
            ((3, 30.0, 2.0, 4.0), one(3)),
            ((2, 3.0, 4.0, 4.0), one(0)),
            ((1, 5.0, 95.0, 4.0), one(1)),
            ((0, 10.0, 10.0, 4.0), one(0)),
            ((2, 0.0, 0.0, 4.0), None),
            ((4, 30.0, 10.0, 4.0), one(3)),
            ((3, 10.0, 5.0, 4.0), two(2, 0.803, 1, 0.173)),
            ((3, 5.0, 10.0, 4.0), two(1, 0.803, 2, 0.173)),
            ((4, 12.0, 8.0, 4.0), two(2, 0.559, 3, 0.324)),
            ((4, 8.0, 12.0, 4.0), two(2, 0.559, 1, 0.324)),
            ((2, 10.0, 10.0, 4.0), two(1, 0.893, 2, 0.107)),
            ((3, 15.0, 15.0, 4.0), two(2, 0.57, 1, 0.383)),
            ((2, 1.0, 30.0, 4.0), one(0)),
            ((2, 30.0, 1.0, 4.0), one(2)),
            ((1, 3.0, 95.0, 4.0), one(0)),
            ((2, 9.0, 11.0, 4.0), one(1)),
            ((3, 8.0, 4.0, 4.0), two(2, 0.738, 1, 0.216)),
            ((4, 2.0, 3.0, 4.0), two(2, 0.333, 1, 0.291)),
            ((2, 6.0, 14.0, 4.0), one(1)),
        ];
        for ((cn, c1, c2, m), want) in golden {
            let got = call_reg1_cn(Some(cn), c1, c2, m);
            assert_eq!(got.len(), 1);
            // CnProb has no PartialEq; compare Debug output, which carries every field.
            assert_eq!(format!("{:?}", got[0]), format!("{:?}", want), "cn={cn} {c1}/{c2} min_read={m}");
        }
    }
}

use statrs::distribution::{Discrete, Poisson};

const POSTERIOR_CUTOFF_STRINGENT: f64 = 0.9;
const ERROR_RATE: f64 = 0.1;

/// Return the reg1 copy number call at each site based on Poisson likelihood,
/// with a minimum read support cutoff.
/// Returns a Vec of 1 element (confident call) or 4 elements (two most likely scenarios).
pub fn call_reg1_cn(full_cn: Option<u32>, count_reg1: f64, count_reg2: f64, min_read: f64) -> Vec<Option<CnProb>> {
    match full_cn {
        None => vec![None.into()],
        Some(0) => vec![Some(CnProb::Single(0)).into()],
        Some(1) if count_reg1 > min_read => vec![Some(CnProb::Single(1)).into()],
        Some(cn) => {
            let nsum = count_reg1 + count_reg2;
            if nsum == 0.0 {
                return vec![None.into()];
            }
            let mut prob = Vec::new();
            for i in 0..=cn {
                let depth_expected = if i == 0 {
                    (ERROR_RATE / 3.0) * nsum
                } else if i == cn {
                    nsum - ERROR_RATE * nsum
                } else {
                    nsum * i as f64 / cn as f64
                };
                // Python: poisson.pmf(int(count_reg1), depthexpected)
                let count = if count_reg1 <= count_reg2 {
                    count_reg1 as u64
                } else {
                    count_reg2 as u64
                };
                let poisson = Poisson::new(depth_expected).unwrap();
                prob.push(poisson.pmf(count));
            }
            let sum_prob: f64 = prob.iter().sum();
            if sum_prob == 0.0 {
                return vec![None.into()];
            }
            let mut post_prob: Vec<f64> = prob.iter().map(|&a| a / sum_prob).collect();
            if count_reg2 < count_reg1 {
                post_prob.reverse();
            }

            let mut post_prob_sorted: Vec<(usize, f64)> =
                post_prob.iter().copied().enumerate().collect();
            post_prob_sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

            let top_idx = post_prob_sorted[0].0;
            let top_prob = post_prob_sorted[0].1;

            if top_idx != 0 && count_reg1 <= min_read && count_reg2 >= min_read {
                return vec![Some(CnProb::Single(0)).into()];
            }

            if top_prob >= POSTERIOR_CUTOFF_STRINGENT {
                return vec![Some(CnProb::Single(top_idx as u32)).into()];
            }

            // Output the two most likely scenarios
            vec![Some(CnProb::TwoOptions {
                cn1: top_idx as u32,
                prob1: (top_prob * 1000.0).round() / 1000.0,
                cn2: post_prob_sorted[1].0 as u32,
                prob2: (post_prob_sorted[1].1 * 1000.0).round() / 1000.0,
            })
            .into()]
        }
    }
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

impl CnProbResult {
    pub fn is_single(&self) -> bool {
        matches!(self.0, None | Some(CnProb::Single(_)))
    }

    pub fn single_value(&self) -> Option<u32> {
        match self.0 {
            Some(CnProb::Single(v)) => Some(v),
            _ => None,
        }
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

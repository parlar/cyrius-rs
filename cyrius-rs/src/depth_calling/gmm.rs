use crate::types::{CnCall, GmmParameter};
use statrs::distribution::{ContinuousCDF, Continuous, Normal};

const SIGMA_CN0: f64 = 0.032;
const PV_CUTOFF: f64 = 1e-3;
const POSTERIOR_CUTOFF: f64 = 0.95;
const DEFAULT_GMM_NSTATE: usize = 11;

/// Gaussian Mixture Model for copy number calling.
pub struct Gmm {
    pub mu_state: Vec<f64>,
    pub sigma_state: Vec<f64>,
    pub prior_state: Vec<f64>,
    pub value_shift: f64,
    pub nstate: usize,
}

impl Gmm {
    pub fn new() -> Self {
        Gmm {
            mu_state: Vec::new(),
            sigma_state: Vec::new(),
            prior_state: Vec::new(),
            value_shift: 1.0,
            nstate: DEFAULT_GMM_NSTATE,
        }
    }

    /// Set GMM parameters from the parameter dictionary for a specific variant ID.
    pub fn set_gmm_par(&mut self, dpar_tmp: &GmmParameter, svid: &str) {
        let gmm_parameter = dpar_tmp
            .get(svid)
            .unwrap_or_else(|| panic!("Variant id {} is not recognized.", svid));

        self.value_shift = gmm_parameter["shift"][0].parse::<f64>().unwrap();

        // Means: modeled as 0, 0.5, 1, 1.5, 2, 2.5, etc.
        let mean_cn2: f64 = gmm_parameter["mean"][0].parse().unwrap();
        let mean_cn3: f64 = gmm_parameter["mean"][1].parse().unwrap();
        self.mu_state = vec![0.0, 0.5, mean_cn2, mean_cn3];
        let mu_width = mean_cn3 - mean_cn2;
        for i in 4..self.nstate {
            self.mu_state.push(1.0 + mu_width * (i as f64 - 2.0));
        }

        // Priors
        let prior_values: Vec<f64> = gmm_parameter["prior"]
            .iter()
            .map(|s| s.parse::<f64>().unwrap())
            .collect();
        let sum_prior: f64 = prior_values.iter().sum();
        if sum_prior >= 1.0 {
            panic!("Sum of priors is larger than 1.");
        }
        self.prior_state = Vec::new();
        for i in 0..self.nstate {
            if i < prior_values.len() {
                self.prior_state.push(prior_values[i]);
            } else {
                let prior_value =
                    (1.0 - sum_prior) / (self.nstate - prior_values.len()) as f64;
                self.prior_state.push(prior_value);
            }
        }

        // Standard deviations
        let sd_cn2: f64 = gmm_parameter["sd"][0].parse().unwrap();
        self.sigma_state = vec![SIGMA_CN0];
        for i in 1..self.nstate {
            let sigma_value = sd_cn2 * (i as f64 / 2.0).sqrt();
            self.sigma_state.push(sigma_value);
        }
    }

    /// Return the final copy number call.
    pub fn gmm_call(&self, val: f64) -> CnCall {
        let val_new = (val / 2.0) / self.value_shift;
        let fcall = self.call_post_prob(val_new, POSTERIOR_CUTOFF);

        let fcall = if let Some(cn) = fcall {
            let (_, gauss_p_value) =
                self.get_gauss_pmf_cdf(val_new, self.mu_state[cn], self.sigma_state[cn]);
            if gauss_p_value < PV_CUTOFF {
                None
            } else {
                Some(cn)
            }
        } else {
            None
        };

        CnCall {
            cn: fcall.map(|c| c as u32),
            depth_value: (val / self.value_shift * 1000.0).round() / 1000.0,
        }
    }

    /// Return the (pdf, p_value) of a gaussian distribution.
    fn get_gauss_pmf_cdf(&self, test_value: f64, gauss_mean: f64, gauss_sd: f64) -> (f64, f64) {
        let test_stats = (test_value - gauss_mean) / gauss_sd;
        let standard_normal = Normal::new(0.0, 1.0).unwrap();
        let pdf = standard_normal.pdf(test_stats) / gauss_sd;
        let cdf_val = standard_normal.cdf(test_stats);
        let p_value = cdf_val.min(1.0 - cdf_val);
        (pdf, p_value)
    }

    /// Return the copy number call based on gaussian mixture model.
    fn call_post_prob(&self, val: f64, post_cutoff: f64) -> Option<usize> {
        let number_state = self.prior_state.len();
        let mut prob = Vec::with_capacity(number_state);
        for i in 0..number_state {
            let (gauss_pmf, _) =
                self.get_gauss_pmf_cdf(val, self.mu_state[i], self.sigma_state[i]);
            prob.push(gauss_pmf * self.prior_state[i]);
        }
        let sum_prob: f64 = prob.iter().sum();
        let post_prob: Vec<f64> = prob.iter().map(|&a| a / sum_prob).collect();
        let max_prob = post_prob
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        if max_prob >= post_cutoff {
            post_prob.iter().position(|&p| p == max_prob)
        } else {
            None
        }
    }
}

//! Quality-aware likelihood model for CN calling.
//!
//! Replaces the fixed ERROR_RATE=0.1 Poisson model with per-read base-quality
//! likelihoods. For each read at a differentiating site, the likelihood of
//! observing the read base under each CN hypothesis uses the actual base quality.
//!
//! This module provides an alternative to `call_reg1_cn` that can be swapped in
//! when the quality_aware feature flag is enabled.

use crate::depth_calling::copy_number_call::CnProb;

const POSTERIOR_CUTOFF_STRINGENT: f64 = 0.9;

/// Round to 3 decimal places using Python's banker's rounding.
fn python_round3(x: f64) -> f64 {
    let scaled = x * 1000.0;
    let rounded = if (scaled - scaled.floor() - 0.5).abs() < 1e-9 {
        let floor = scaled.floor();
        if floor as i64 % 2 == 0 {
            floor
        } else {
            floor + 1.0
        }
    } else {
        scaled.round()
    };
    rounded / 1000.0
}

/// Per-read evidence at a variant site: base quality and whether it supports reg1.
#[derive(Debug, Clone, Copy)]
pub struct ReadEvidence {
    /// Base quality (phred score)
    pub base_quality: u8,
    /// True if the read supports region 1 (e.g., CYP2D6 allele)
    pub supports_reg1: bool,
}

/// Quality-aware CN call using per-read likelihoods.
///
/// Returns `Vec<Option<CnProb>>` — same interface as `call_reg1_cn`.
///
/// Instead of counting reads and using Poisson(count | expected), we compute:
///   P(reads | CN=i) = product over reads of P(read_base | CN=i)
///
/// For a read supporting reg1 with error rate e = 10^(-BQ/10):
///   P(read | CN=i) = (i/cn) * (1-e) + ((cn-i)/cn) * (e/3)   for 0 < i < cn
///   P(read | CN=0) = e/3
///   P(read | CN=cn) = 1-e
///
/// Similarly for reads supporting reg2 (swap i with cn-i).
pub fn call_reg1_cn_quality_aware(
    full_cn: Option<u32>,
    reads: &[ReadEvidence],
    min_read: f64,
) -> Vec<Option<CnProb>> {
    match full_cn {
        None => vec![None],
        Some(0) => vec![Some(CnProb::Single(0))],
        Some(cn) => {
            let count_reg1 = reads.iter().filter(|r| r.supports_reg1).count() as f64;
            let count_reg2 = reads.iter().filter(|r| !r.supports_reg1).count() as f64;

            if reads.is_empty() {
                return vec![None];
            }

            // Special case: CN=1 with enough reg1 reads
            if cn == 1 && count_reg1 > min_read {
                return vec![Some(CnProb::Single(1))];
            }

            // Compute log-likelihoods for each CN hypothesis
            let mut log_likelihoods = vec![0.0f64; (cn + 1) as usize];

            for read in reads {
                let error_rate = 10.0f64.powf(-(read.base_quality as f64) / 10.0);
                let error_rate = error_rate.max(0.001).min(0.5); // clamp to reasonable range

                for i in 0..=cn {
                    let p = if read.supports_reg1 {
                        if i == 0 {
                            error_rate / 3.0
                        } else if i == cn {
                            1.0 - error_rate
                        } else {
                            (i as f64 / cn as f64) * (1.0 - error_rate)
                                + ((cn - i) as f64 / cn as f64) * (error_rate / 3.0)
                        }
                    } else {
                        if i == 0 {
                            1.0 - error_rate
                        } else if i == cn {
                            error_rate / 3.0
                        } else {
                            ((cn - i) as f64 / cn as f64) * (1.0 - error_rate)
                                + (i as f64 / cn as f64) * (error_rate / 3.0)
                        }
                    };

                    log_likelihoods[i as usize] += p.max(1e-300).ln();
                }
            }

            // Convert log-likelihoods to posterior probabilities (uniform prior)
            let max_ll = log_likelihoods
                .iter()
                .cloned()
                .fold(f64::NEG_INFINITY, f64::max);
            let probs: Vec<f64> = log_likelihoods
                .iter()
                .map(|&ll| (ll - max_ll).exp())
                .collect();
            let sum_prob: f64 = probs.iter().sum();

            if sum_prob == 0.0 {
                return vec![None];
            }

            let post_prob: Vec<f64> = probs.iter().map(|&p| p / sum_prob).collect();

            let mut post_prob_sorted: Vec<(usize, f64)> =
                post_prob.iter().copied().enumerate().collect();
            post_prob_sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

            let top_idx = post_prob_sorted[0].0;
            let top_prob = post_prob_sorted[0].1;

            // Min-read filter: if few reg1 reads and many reg2, call CN=0
            if top_idx != 0 && count_reg1 <= min_read && count_reg2 >= min_read {
                return vec![Some(CnProb::Single(0))];
            }

            if top_prob >= POSTERIOR_CUTOFF_STRINGENT {
                return vec![Some(CnProb::Single(top_idx as u32))];
            }

            vec![Some(CnProb::TwoOptions {
                cn1: top_idx as u32,
                prob1: python_round3(top_prob),
                cn2: post_prob_sorted[1].0 as u32,
                prob2: python_round3(post_prob_sorted[1].1),
            })]
        }
    }
}

/// Collect per-read evidence from pileup at a SNP site.
/// Each read contributes at most one evidence entry (reg1 takes priority if both match).
pub fn collect_read_evidence(
    reader: &mut rust_htslib::bam::IndexedReader,
    nchr: &str,
    snp_position: i64,
    reg1_allele: &str,
    reg2_allele: &str,
) -> Vec<ReadEvidence> {
    use rust_htslib::bam::{self, Read as BamRead};

    let mut evidence = Vec::new();

    let tid = match reader.header().tid(nchr.as_bytes()) {
        Some(t) => t,
        None => return evidence,
    };

    reader
        .fetch(bam::FetchDefinition::Region(
            tid as i32,
            snp_position - 1,
            snp_position,
        ))
        .unwrap();

    let mut pileups = reader.pileup();
    pileups.set_max_depth(i32::MAX as u32);

    for pileup_result in pileups {
        let pileup = match pileup_result {
            Ok(p) => p,
            Err(_) => continue,
        };
        let site_position = pileup.pos() as i64 + 1;
        if site_position != snp_position {
            continue;
        }

        for alignment in pileup.alignments() {
            let record = alignment.record();
            if alignment.is_del() || alignment.is_refskip() {
                continue;
            }
            if record.is_secondary() || record.is_supplementary() || record.is_duplicate() {
                continue;
            }

            let qpos = match alignment.qpos() {
                Some(p) => p,
                None => continue,
            };

            let seq = record.seq().as_bytes();
            let quals = record.qual();

            let bq = if qpos < quals.len() { quals[qpos] } else { 20 };

            // Check reg1 first, then reg2. Each read contributes at most once.
            let mut matched = false;
            for allele in reg1_allele.split(',') {
                let end_pos = qpos + allele.len();
                if end_pos <= seq.len() {
                    let read_seq = String::from_utf8_lossy(&seq[qpos..end_pos]);
                    if read_seq == allele {
                        evidence.push(ReadEvidence {
                            base_quality: bq,
                            supports_reg1: true,
                        });
                        matched = true;
                        break;
                    }
                }
            }
            if !matched {
                for allele in reg2_allele.split(',') {
                    let end_pos = qpos + allele.len();
                    if end_pos <= seq.len() {
                        let read_seq = String::from_utf8_lossy(&seq[qpos..end_pos]);
                        if read_seq == allele {
                            evidence.push(ReadEvidence {
                                base_quality: bq,
                                supports_reg1: false,
                            });
                            break;
                        }
                    }
                }
            }
        }
    }

    evidence
}

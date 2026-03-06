use crate::types::SnpLookup;
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::HashSet;

/// Check whether a pileup read passes basic filters.
/// Matches Python passing_read().
pub fn passing_read(record: &bam::Record, is_del: bool, is_refskip: bool) -> bool {
    !is_del
        && !is_refskip
        && !record.is_secondary()
        && !record.is_supplementary()
        && !record.is_duplicate()
}

/// Check whether a read passes more stringent filters.
pub fn passing_read_stringent(record: &bam::Record) -> bool {
    let insert_size = record.insert_size();
    insert_size.unsigned_abs() < 1000
}

/// Return the number of reads supporting region1 and region2, forward and reverse.
/// Matches Python get_reads_by_region().
pub fn get_reads_by_region(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    dsnp: &indexmap::IndexMap<String, String>,
    dindex: &std::collections::HashMap<String, usize>,
    min_mapq: u8,
    min_insert_size: Option<i64>,
    stringent: bool,
) -> (
    Vec<HashSet<String>>,
    Vec<HashSet<String>>,
    Vec<HashSet<String>>,
    Vec<HashSet<String>>,
) {
    let n = dsnp.len();
    let mut lsnp1_forward: Vec<HashSet<String>> = (0..n).map(|_| HashSet::new()).collect();
    let mut lsnp1_reverse: Vec<HashSet<String>> = (0..n).map(|_| HashSet::new()).collect();
    let mut lsnp2_forward: Vec<HashSet<String>> = (0..n).map(|_| HashSet::new()).collect();
    let mut lsnp2_reverse: Vec<HashSet<String>> = (0..n).map(|_| HashSet::new()).collect();

    // Get the target ID for the chromosome
    let tid = match reader.header().tid(nchr.as_bytes()) {
        Some(t) => t,
        None => return (lsnp1_forward, lsnp1_reverse, lsnp2_forward, lsnp2_reverse),
    };

    for (snp_position_ori, alleles) in dsnp.iter() {
        let snp_position: i64 = match snp_position_ori.split('_').next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let dsnp_index = match dindex.get(snp_position_ori) {
            Some(&v) => v,
            None => continue,
        };

        let (reg1_allele, reg2_allele) = match alleles.split_once('_') {
            Some(pair) => pair,
            None => continue,
        };

        // Pileup at this position
        reader
            .fetch(bam::FetchDefinition::Region(
                tid as i32,
                snp_position as i64 - 1,
                snp_position as i64,
            ))
            .unwrap();

        let mut pileups = reader.pileup();
        pileups.set_max_depth(i32::MAX as u32);

        for pileup_result in pileups {
            let pileup = pileup_result.unwrap();
            let site_position = pileup.pos() as i64 + 1;
            if site_position != snp_position {
                continue;
            }

            for alignment in pileup.alignments() {
                let record = alignment.record();
                let is_del = alignment.is_del();
                let is_refskip = alignment.is_refskip();

                if !passing_read(&record, is_del, is_refskip) {
                    continue;
                }
                if record.mapq() < min_mapq {
                    continue;
                }
                if let Some(min_ins) = min_insert_size {
                    if (record.insert_size().unsigned_abs() as i64) < min_ins {
                        continue;
                    }
                }
                if stringent && !passing_read_stringent(&record) {
                    continue;
                }

                let qpos = match alignment.qpos() {
                    Some(p) => p,
                    None => continue,
                };

                let read_name =
                    String::from_utf8_lossy(record.qname()).to_string();
                let seq = record.seq().as_bytes();

                // Check reg1 alleles
                for allele in reg1_allele.split(',') {
                    let end_pos = qpos + allele.len();
                    if end_pos <= seq.len() {
                        let read_seq =
                            String::from_utf8_lossy(&seq[qpos..end_pos]);
                        if read_seq == allele {
                            if record.is_reverse() {
                                lsnp1_reverse[dsnp_index].insert(read_name.clone());
                            } else {
                                lsnp1_forward[dsnp_index].insert(read_name.clone());
                            }
                        }
                    }
                }

                // Check reg2 alleles
                for allele in reg2_allele.split(',') {
                    let end_pos = qpos + allele.len();
                    if end_pos <= seq.len() {
                        let read_seq =
                            String::from_utf8_lossy(&seq[qpos..end_pos]);
                        if read_seq == allele {
                            if record.is_reverse() {
                                lsnp2_reverse[dsnp_index].insert(read_name.clone());
                            } else {
                                lsnp2_forward[dsnp_index].insert(read_name.clone());
                            }
                        }
                    }
                }
            }
        }
    }

    (lsnp1_forward, lsnp1_reverse, lsnp2_forward, lsnp2_reverse)
}

/// Merge sets of reads (forward/reverse, region1/region2).
pub fn merge_reads(lists: &[&Vec<HashSet<String>>]) -> Vec<HashSet<String>> {
    let n = lists[0].len();
    let mut merged = Vec::with_capacity(n);
    for i in 0..n {
        let mut combined = HashSet::new();
        for list in lists {
            combined.extend(list[i].iter().cloned());
        }
        merged.push(combined);
    }
    merged
}

/// Return the fraction of reads supporting region1.
pub fn get_fraction(lsnp1: &[usize], lsnp2: &[usize]) -> Vec<f64> {
    lsnp1
        .iter()
        .zip(lsnp2.iter())
        .map(|(&s1, &s2)| {
            let sum = s1 + s2;
            if sum == 0 {
                0.0
            } else {
                s1 as f64 / sum as f64
            }
        })
        .collect()
}

/// Read-pair crossing statistics from diagnostic SNP positions.
#[derive(Debug, Clone)]
pub struct CrossingStats {
    pub n_consistent: usize,
    pub n_crossing: usize,
    pub crossing_fraction: f64,
}

/// Compute read-pair crossing fraction from per-position D6/D7 read name sets.
/// A read name appearing in both D6 and D7 sets (across different positions)
/// indicates a crossing event (mate pair spans paralog boundary).
fn compute_crossing(d6_reads: &[HashSet<String>], d7_reads: &[HashSet<String>]) -> CrossingStats {
    use std::collections::HashMap;
    // Build: read_name -> (seen_in_d6, seen_in_d7)
    let mut assignments: HashMap<&str, (bool, bool)> = HashMap::new();

    for pos_set in d6_reads {
        for name in pos_set {
            assignments.entry(name.as_str()).or_insert((false, false)).0 = true;
        }
    }
    for pos_set in d7_reads {
        for name in pos_set {
            assignments.entry(name.as_str()).or_insert((false, false)).1 = true;
        }
    }

    let mut n_consistent = 0usize;
    let mut n_crossing = 0usize;
    for (_, (is_d6, is_d7)) in &assignments {
        if *is_d6 && *is_d7 {
            n_crossing += 1;
        } else if *is_d6 || *is_d7 {
            n_consistent += 1;
        }
    }

    let total = n_consistent + n_crossing;
    let crossing_fraction = if total > 0 {
        n_crossing as f64 / total as f64
    } else {
        0.0
    };

    CrossingStats {
        n_consistent,
        n_crossing,
        crossing_fraction,
    }
}

/// Return the number of supporting reads at each position in both region1 and region2,
/// plus read-pair crossing statistics.
pub fn get_supporting_reads(
    reader: &mut bam::IndexedReader,
    snp_db: &SnpLookup,
) -> (Vec<usize>, Vec<usize>, CrossingStats) {
    if snp_db.dsnp1.len() != snp_db.dsnp2.len() {
        log::warn!("SNP database region sizes differ: {} vs {}", snp_db.dsnp1.len(), snp_db.dsnp2.len());
        return (Vec::new(), Vec::new(), CrossingStats { n_consistent: 0, n_crossing: 0, crossing_fraction: 0.0 });
    }

    let (lsnp1_reg1_for, lsnp1_reg1_rev, lsnp2_reg1_for, lsnp2_reg1_rev) =
        get_reads_by_region(reader, &snp_db.nchr, &snp_db.dsnp1, &snp_db.dindex, 0, None, false);
    let (lsnp1_reg2_for, lsnp1_reg2_rev, lsnp2_reg2_for, lsnp2_reg2_rev) =
        get_reads_by_region(reader, &snp_db.nchr, &snp_db.dsnp2, &snp_db.dindex, 0, None, false);

    let lsnp1 = merge_reads(&[
        &lsnp1_reg1_for,
        &lsnp1_reg1_rev,
        &lsnp1_reg2_for,
        &lsnp1_reg2_rev,
    ]);
    let lsnp2 = merge_reads(&[
        &lsnp2_reg1_for,
        &lsnp2_reg1_rev,
        &lsnp2_reg2_for,
        &lsnp2_reg2_rev,
    ]);

    let crossing = compute_crossing(&lsnp1, &lsnp2);

    (
        lsnp1.iter().map(|s| s.len()).collect(),
        lsnp2.iter().map(|s| s.len()).collect(),
        crossing,
    )
}

/// Return the number of supporting reads at each position only in region1.
pub fn get_supporting_reads_single_region(
    reader: &mut bam::IndexedReader,
    snp_db: &SnpLookup,
    min_insert_size: Option<i64>,
) -> (Vec<usize>, Vec<usize>, Vec<usize>, Vec<usize>) {
    let (lsnp1_for, lsnp1_rev, lsnp2_for, lsnp2_rev) = get_reads_by_region(
        reader,
        &snp_db.nchr,
        &snp_db.dsnp1,
        &snp_db.dindex,
        10,
        min_insert_size,
        false,
    );

    let lsnp1 = merge_reads(&[&lsnp1_for, &lsnp1_rev]);
    let lsnp2 = merge_reads(&[&lsnp2_for, &lsnp2_rev]);

    (
        lsnp1.iter().map(|s| s.len()).collect(),
        lsnp2.iter().map(|s| s.len()).collect(),
        lsnp1_for.iter().map(|s| s.len()).collect(),
        lsnp1_rev.iter().map(|s| s.len()).collect(),
    )
}

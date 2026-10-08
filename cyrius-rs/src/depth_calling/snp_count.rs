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

        // A site's alleles can differ in length (indels), and the two lengths
        // that differ belong to DIFFERENT regions -- e.g. a 22 bp reg1 ALT
        // against a 13 bp reg2 REF at g.42126656-42126657insAGTGGGCAC. Gating
        // each allele by its own `qpos + allele.len() <= seq.len()` therefore
        // let a read ending within `max_len` of the site match the SHORTER
        // allele while the longer one was never even tested, so reads near
        // their 3' end were counted only for the shorter allele and the allele
        // fractions were biased (review finding #6). The bound is taken over
        // the union of both regions' alleles, so a read without room for the
        // longest allele at the site is counted for neither.
        let max_allele_len = reg1_allele
            .split(',')
            .chain(reg2_allele.split(','))
            .map(str::len)
            .max()
            .unwrap_or(0);

        // Pileup at this position. A fetch failure (bad region / corrupt index)
        // means no reads here — skip rather than aborting the whole caller.
        if reader
            .fetch(bam::FetchDefinition::Region(
                tid as i32,
                snp_position as i64 - 1,
                snp_position as i64,
            ))
            .is_err()
        {
            continue;
        }

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

                let seq = record.seq().as_bytes();
                // One bound for every allele of the site: either the read has
                // room for all of them or it is counted for none.
                if qpos + max_allele_len > seq.len() {
                    continue;
                }

                let read_name =
                    String::from_utf8_lossy(record.qname()).to_string();

                // Check reg1 alleles
                for allele in reg1_allele.split(',') {
                    let read_seq =
                        String::from_utf8_lossy(&seq[qpos..qpos + allele.len()]);
                    if read_seq == allele {
                        if record.is_reverse() {
                            lsnp1_reverse[dsnp_index].insert(read_name.clone());
                        } else {
                            lsnp1_forward[dsnp_index].insert(read_name.clone());
                        }
                    }
                }

                // Check reg2 alleles
                for allele in reg2_allele.split(',') {
                    let read_seq =
                        String::from_utf8_lossy(&seq[qpos..qpos + allele.len()]);
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

#[cfg(test)]
mod tests {
    use super::*;

    /// The real alleles of `g.42126656-42126657insAGTGGGCAC` from
    /// `data/CYP2D6_target_variant_38.txt`: reg1 (D6 ALT) is 22 bp, reg2 (D7
    /// REF) is 13 bp, and they already differ at offset 12 (`A` vs `C`), so a
    /// read carrying ALT can never be mistaken for REF on its first 13 bases.
    const ALT_22: &str = "TCCAGTGGGCACAGTGGGCACC";
    const REF_13: &str = "TCCAGTGGGCACC";
    /// 1-based site position, matching the `dsnp` key below.
    const SITE_1BASED: i64 = 1001;
    /// Every fixture read is 50M, so there is room for the 22 bp allele only
    /// when the site sits at least 22 bases from the read's 3' end.
    const READ_LEN: usize = 50;

    /// Write an indexed single-site BAM holding three 50M reads.
    ///
    /// * `ref_at_end` — positioned so exactly 13 bases remain from the site:
    ///   enough for the 13 bp REF, not for the 22 bp ALT. This is the biased
    ///   read; it must be counted for neither allele.
    /// * `ref_roomy` — 10 bases into the read, so both alleles fit. Counted
    ///   as REF.
    /// * `alt_roomy` — same offset, carrying ALT. Counted as ALT.
    ///
    /// NOTE: `bam::Record::new()` leaves `flags = 4` (unmapped), and htslib's
    /// `bam_endpos` ignores the CIGAR for an unmapped record, so an unflagged
    /// record spans only `[pos, pos + 1)` and silently drops out of region
    /// queries. The `set_flags(0)` calls below are load-bearing.
    fn write_indel_end_bias_bam(dir: &std::path::Path) -> std::path::PathBuf {
        let mut header = bam::header::Header::new();
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr22");
        sq.push_tag(b"LN", 50_818_468);
        header.push_record(&sq);

        let path = dir.join("indel_end_bias.bam");
        let site_0based = SITE_1BASED - 1;
        {
            let mut writer =
                bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
            let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(
                READ_LEN as u32,
            )]);
            let qual = vec![30u8; READ_LEN];

            let mut emit = |name: &[u8], qpos: usize, allele: &str| {
                let mut seq = vec![b'G'; READ_LEN];
                seq[qpos..qpos + allele.len()].copy_from_slice(allele.as_bytes());
                let mut rec = bam::Record::new();
                rec.set(name, Some(&cigar), &seq, &qual);
                rec.set_tid(0);
                rec.set_pos(site_0based - qpos as i64);
                rec.set_mapq(60);
                rec.set_flags(0); // mapped, primary, forward, not a duplicate
                writer.write(&rec).unwrap();
            };

            // Sorted by start position: qpos 37 starts at 963, qpos 10 at 990.
            emit(b"ref_at_end", READ_LEN - REF_13.len(), REF_13);
            emit(b"ref_roomy", 10, REF_13);
            emit(b"alt_roomy", 10, ALT_22);
        }
        bam::index::build(&path, None::<&std::path::PathBuf>, bam::index::Type::Bai, 1)
            .unwrap();
        path
    }

    fn names(sets: &[&Vec<HashSet<String>>]) -> Vec<String> {
        let mut v: Vec<String> = sets
            .iter()
            .flat_map(|s| s[0].iter().cloned())
            .collect();
        v.sort();
        v
    }

    #[test]
    fn a_read_without_room_for_the_longer_allele_counts_for_neither() {
        // Review finding #6. Each allele used to be gated by its own
        // `qpos + allele.len() <= seq.len()`, so a read ending within
        // `max_len` of the site could match the SHORTER allele and never the
        // longer one. `ref_at_end` has exactly 13 bases left: it matched the
        // 13 bp REF while the 22 bp ALT was never even tested, so reads near
        // their 3' end were counted only for the shorter allele and the
        // allele fractions were biased. The bound must be taken over the
        // union of BOTH regions' alleles, because the two lengths that differ
        // here belong to different regions.
        let dir = std::env::temp_dir()
            .join(format!("cyrius_indel_end_bias_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bam_path = write_indel_end_bias_bam(&dir);
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();

        let key = format!("{SITE_1BASED}_0");
        let mut dsnp: indexmap::IndexMap<String, String> = indexmap::IndexMap::new();
        dsnp.insert(key.clone(), format!("{ALT_22}_{REF_13}"));
        let mut dindex: std::collections::HashMap<String, usize> =
            std::collections::HashMap::new();
        dindex.insert(key, 0);

        let (reg1_fwd, reg1_rev, reg2_fwd, reg2_rev) =
            get_reads_by_region(&mut reader, "chr22", &dsnp, &dindex, 0, None, false);

        let reg1 = names(&[&reg1_fwd, &reg1_rev]);
        let reg2 = names(&[&reg2_fwd, &reg2_rev]);

        // Clean up before asserting so a failure leaves nothing behind.
        drop(reader);
        let _ = std::fs::remove_dir_all(&dir);

        assert_eq!(
            reg1,
            vec!["alt_roomy".to_string()],
            "only the read with room for the 22 bp ALT may support reg1; got {reg1:?}"
        );
        assert_eq!(
            reg2,
            vec!["ref_roomy".to_string()],
            "`ref_at_end` has exactly {} bases left from the site -- enough for \
             the {} bp REF but not for the {} bp ALT -- so it must be counted \
             for NEITHER allele, not for the shorter one; got {reg2:?}",
            REF_13.len(),
            REF_13.len(),
            ALT_22.len()
        );
    }

    #[test]
    fn equal_length_alleles_are_unaffected_by_the_read_end_bound() {
        // Guard against over-correcting: when both regions' alleles are the
        // same length (every site in `data/CYP2D6_SNP_38.txt`), a read whose
        // last base is the site must still be counted.
        let dir = std::env::temp_dir()
            .join(format!("cyrius_indel_end_bias_eq_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bam_path = write_indel_end_bias_bam(&dir);
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();

        // 1 bp alleles: `ref_at_end` carries 'T' (REF_13[0]) at the site.
        let key = format!("{SITE_1BASED}_0");
        let mut dsnp: indexmap::IndexMap<String, String> = indexmap::IndexMap::new();
        dsnp.insert(key.clone(), "G_T".to_string());
        let mut dindex: std::collections::HashMap<String, usize> =
            std::collections::HashMap::new();
        dindex.insert(key, 0);

        let (reg1_fwd, reg1_rev, reg2_fwd, reg2_rev) =
            get_reads_by_region(&mut reader, "chr22", &dsnp, &dindex, 0, None, false);

        let reg1 = names(&[&reg1_fwd, &reg1_rev]);
        let reg2 = names(&[&reg2_fwd, &reg2_rev]);

        drop(reader);
        let _ = std::fs::remove_dir_all(&dir);

        // `alt_roomy` has 'T' too (ALT_22[0]), `ref_roomy` has 'T'.
        assert_eq!(
            reg2,
            vec![
                "alt_roomy".to_string(),
                "ref_at_end".to_string(),
                "ref_roomy".to_string()
            ],
            "all three reads carry 'T' at the site, and a 1 bp allele fits in \
             every one of them; got {reg2:?}"
        );
        assert!(reg1.is_empty(), "no read carries 'G' at the site; got {reg1:?}");
    }
}

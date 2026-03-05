use crate::depth_calling::snp_count::passing_read;
use crate::types::SnpLookup;
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::HashMap;

/// Get haplotypes for each read at a set of target positions (both regions).
pub fn get_haplotypes_from_bam(
    reader: &mut bam::IndexedReader,
    base_db: &SnpLookup,
    target_positions: &[usize],
) -> HashMap<String, String> {
    let dread = get_bases_per_read(reader, base_db, target_positions, None, 0);
    let (base1, base2) = get_base1_base2(base_db, target_positions);
    get_hap_counts(&dread, &base1, &base2, target_positions)
}

/// Get haplotypes for each read at a set of target positions (single region only).
pub fn get_haplotypes_from_bam_single_region(
    reader: &mut bam::IndexedReader,
    base_db: &SnpLookup,
    target_positions: &[usize],
) -> HashMap<String, String> {
    let dread = get_bases_per_read(reader, base_db, target_positions, Some(0), 10);
    let (base1, base2) = get_base1_base2(base_db, target_positions);
    get_hap_counts(&dread, &base1, &base2, target_positions)
}

/// Get bases per read at target positions from pileup data.
/// region: None = both, Some(0) = dsnp1 only, Some(1) = dsnp2 only.
fn get_bases_per_read(
    reader: &mut bam::IndexedReader,
    base_db: &SnpLookup,
    target_positions: &[usize],
    region: Option<usize>,
    min_mapq: u8,
) -> HashMap<String, HashMap<usize, Option<String>>> {
    let mut dread: HashMap<String, HashMap<usize, Option<String>>> = HashMap::new();
    let nchr = &base_db.nchr;
    let dindex = &base_db.dindex;

    let dsnps: Vec<&indexmap::IndexMap<String, String>> = match region {
        Some(0) => vec![&base_db.dsnp1],
        Some(1) => vec![&base_db.dsnp2],
        _ => vec![&base_db.dsnp1, &base_db.dsnp2],
    };

    let tid = match reader.header().tid(nchr.as_bytes()) {
        Some(t) => t,
        None => return dread,
    };

    for dsnp in &dsnps {
        for (snp_position_ori, alleles) in dsnp.iter() {
            let dsnp_index = dindex[snp_position_ori];
            let snp_position: i64 = snp_position_ori.split('_').next().unwrap().parse().unwrap();

            if !target_positions.contains(&dsnp_index) {
                continue;
            }

            reader
                .fetch(bam::FetchDefinition::Region(
                    tid as i32,
                    snp_position - 1,
                    snp_position,
                ))
                .unwrap();

            let mut pileups = reader.pileup();
            pileups.set_max_depth(u32::MAX);

            for pileup_result in pileups {
                let pileup = pileup_result.unwrap();
                let site_position = pileup.pos() as i64 + 1;
                if site_position != snp_position {
                    continue;
                }

                let (reg1_allele, reg2_allele) = alleles.split_once('_').unwrap();

                for alignment in pileup.alignments() {
                    let record = alignment.record();
                    let is_del = alignment.is_del();
                    let is_refskip = alignment.is_refskip();

                    if !passing_read(&record, is_del, is_refskip) || record.mapq() < min_mapq {
                        continue;
                    }

                    let qpos = match alignment.qpos() {
                        Some(p) => p,
                        None => continue,
                    };

                    let read_name = String::from_utf8_lossy(record.qname()).to_string();
                    let seq = record.seq().as_bytes();

                    let min_len = std::cmp::min(reg1_allele.len(), reg2_allele.len());
                    let end_pos = qpos + min_len;
                    if end_pos >= seq.len() {
                        continue;
                    }

                    let hap = String::from_utf8_lossy(&seq[qpos..end_pos]).to_string();

                    let entry = dread.entry(read_name).or_insert_with(|| {
                        let mut m = HashMap::new();
                        for &pos in target_positions {
                            m.insert(pos, None);
                        }
                        m
                    });

                    let current = entry.get(&dsnp_index).cloned().flatten();
                    if current.is_some() && current.as_deref() != Some(&hap) {
                        // Conflicting base — set to None
                        entry.insert(dsnp_index, None);
                    } else {
                        entry.insert(dsnp_index, Some(hap));
                    }
                }
            }
        }
    }
    dread
}

/// Translate bases into haplotype strings (1 = reg1 allele, 2 = reg2 allele, x = unknown).
fn get_hap_counts(
    dread: &HashMap<String, HashMap<usize, Option<String>>>,
    base1: &[Option<String>],
    base2: &[Option<String>],
    target_positions: &[usize],
) -> HashMap<String, String> {
    let mut dread_hap = HashMap::new();
    for (read_name, read_bases) in dread {
        let mut pos_list: Vec<char> = vec!['x'; target_positions.len()];
        for (i, &pos) in target_positions.iter().enumerate() {
            if let Some(base_opt) = read_bases.get(&pos) {
                if let Some(ref base) = base_opt {
                    if let Some(ref b1) = base1[i] {
                        for allele in b1.split(',') {
                            if base == &allele.to_uppercase() {
                                pos_list[i] = '1';
                            }
                        }
                    }
                    if let Some(ref b2) = base2[i] {
                        for allele in b2.split(',') {
                            if base == &allele.to_uppercase() {
                                pos_list[i] = '2';
                            }
                        }
                    }
                }
            }
        }
        dread_hap.insert(read_name.clone(), pos_list.iter().collect());
    }
    dread_hap
}

/// Get expected bases at target positions for each allele/paralog.
fn get_base1_base2(
    base_db: &SnpLookup,
    target_positions: &[usize],
) -> (Vec<Option<String>>, Vec<Option<String>>) {
    let n = base_db.dsnp1.len();
    let mut base1: Vec<Option<String>> = vec![None; n];
    let mut base2: Vec<Option<String>> = vec![None; n];

    for (pos, alleles) in &base_db.dsnp1 {
        let dsnp_index = base_db.dindex[pos];
        if target_positions.contains(&dsnp_index) {
            let index: usize = pos.split('_').nth(1).unwrap().parse().unwrap();
            let (allele1, allele2) = alleles.split_once('_').unwrap();
            base1[index] = Some(allele1.to_string());
            base2[index] = Some(allele2.to_string());
        }
    }

    (
        target_positions.iter().map(|&i| base1[i].clone()).collect(),
        target_positions.iter().map(|&i| base2[i].clone()).collect(),
    )
}

/// Extract haplotypes at certain positions and count them.
pub fn extract_hap(
    dhaplotype: &HashMap<String, String>,
    positions_to_extract: &[usize],
) -> HashMap<String, Vec<u32>> {
    let mut hap_count: HashMap<String, Vec<u32>> = HashMap::new();
    for hap in dhaplotype.values() {
        let hap_chars: Vec<char> = hap.chars().collect();
        let hap_base: Vec<char> = positions_to_extract.iter().map(|&p| hap_chars[p]).collect();
        if hap_base.contains(&'x') {
            continue;
        }
        let key: String = hap_base.iter().collect();
        hap_count.entry(key).or_insert_with(Vec::new).push(1);
    }
    hap_count
}

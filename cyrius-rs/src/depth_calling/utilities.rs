use crate::types::{GmmParameter, Region, RegionDic, SnpLookup};
use indexmap::IndexMap;
use rust_htslib::bam;
use std::collections::HashMap;

/// Parse a BED region file from an embedded string.
/// Returns a RegionDic mapping region_type -> Vec<(Region, gc_content)>.
pub fn parse_region_file(content: &str, genome: &str) -> RegionDic {
    let mut region_dic = RegionDic::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 6 {
            continue;
        }
        let mut nchr = fields[0].to_string();
        if genome == "38" {
            nchr = nchr.replace("chr", "");
        }
        let region_start: i64 = fields[1].parse().unwrap();
        let region_end: i64 = fields[2].parse().unwrap();
        let region_name = fields[3].to_string();
        let region_type = fields[4].to_string();
        let region_gc: f64 = fields[5].parse().unwrap();

        let region: Region = (nchr, region_start, region_end, region_name);
        region_dic
            .entry(region_type)
            .or_insert_with(Vec::new)
            .push((region, region_gc));
    }
    region_dic
}

/// Parse GMM parameter file from an embedded string.
/// Returns variant_id -> parameter_name -> list of values.
pub fn parse_gmm_file(content: &str) -> GmmParameter {
    let mut dpar: GmmParameter = HashMap::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 3 {
            continue;
        }
        let variant_id = fields[0].to_string();
        let param_name = fields[1].to_string();
        let values: Vec<String> = fields[2..]
            .iter()
            .map(|s| {
                // Values are in format "key:value", extract the value part
                if let Some(pos) = s.find(':') {
                    s[pos + 1..].to_string()
                } else {
                    s.to_string()
                }
            })
            .collect();
        dpar.entry(variant_id)
            .or_insert_with(HashMap::new)
            .entry(param_name)
            .or_insert(values);
    }
    dpar
}

/// Open a BAM or CRAM alignment file using rust-htslib.
pub fn open_alignment_file(
    alignment_file: &str,
    reference_fasta: Option<&str>,
) -> Result<bam::IndexedReader, Box<dyn std::error::Error>> {
    open_alignment_file_with_index(alignment_file, reference_fasta, None)
}

/// Open a BAM or CRAM alignment file with optional custom index path.
pub fn open_alignment_file_with_index(
    alignment_file: &str,
    reference_fasta: Option<&str>,
    index_filename: Option<&str>,
) -> Result<bam::IndexedReader, Box<dyn std::error::Error>> {
    let mut reader = if alignment_file.to_lowercase().ends_with("cram") {
        let mut r = match index_filename {
            Some(idx) => bam::IndexedReader::from_path_and_index(alignment_file, idx)?,
            None => bam::IndexedReader::from_path(alignment_file)?,
        };
        if let Some(ref_path) = reference_fasta {
            r.set_reference(ref_path)?;
        }
        r
    } else {
        match index_filename {
            Some(idx) => bam::IndexedReader::from_path_and_index(alignment_file, idx)?,
            None => bam::IndexedReader::from_path(alignment_file)?,
        }
    };

    use rust_htslib::bam::Read as HtslibRead;
    reader.set_threads(1)?;
    Ok(reader)
}

/// Parse SNP position file from an embedded string.
/// Matches Python get_snp_position().
pub fn get_snp_position(content: &str, genome: &str, group: Option<&str>) -> SnpLookup {
    let mut dsnp1 = IndexMap::new();
    let mut dsnp2 = IndexMap::new();
    let mut dindex: HashMap<String, usize> = HashMap::new();
    let mut counter: i32 = -1;
    let mut last_nchr = String::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with('\n') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 5 {
            continue;
        }

        let line_group = fields.last().unwrap();
        if group.is_some() && Some(*line_group) != group {
            continue;
        }

        counter += 1;
        let idx = counter as usize;

        let reg1_name = format!("{}_{}", fields[1], idx);
        let reg2_name = format!("{}_{}", fields[3], idx);
        let reg1_base = fields[2].to_uppercase();
        let reg2_base = fields[4].to_uppercase();

        if *line_group != "-" {
            dsnp1
                .entry(reg1_name.clone())
                .or_insert_with(|| format!("{}_{}", reg1_base, reg2_base));
            dsnp2
                .entry(reg2_name.clone())
                .or_insert_with(|| format!("{}_{}", reg1_base, reg2_base));
        } else {
            dsnp1
                .entry(reg1_name.clone())
                .or_insert_with(|| format!("{}_{}", reg1_base, reverse_complement(&reg2_base)));
            dsnp2
                .entry(reg2_name.clone())
                .or_insert_with(|| format!("{}_{}", reverse_complement(&reg1_base), reg2_base));
        }
        dindex.insert(reg1_name, idx);
        dindex.insert(reg2_name, idx);
        last_nchr = fields[0].to_string();
    }

    if counter == -1 {
        panic!("No valid SNP positions found in file");
    }

    let nchr = if genome == "38" {
        last_nchr.replace("chr", "")
    } else {
        last_nchr
    };

    SnpLookup {
        dsnp1,
        dsnp2,
        nchr,
        dindex,
    }
}

/// Get variant names from a variant file (embedded string).
pub fn get_var_names(content: &str) -> Vec<String> {
    let mut var_list = Vec::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if let Some(var_name) = fields.last() {
            var_list.push(var_name.to_string());
        }
    }
    var_list
}

const COMPLEMENT: &[(char, char)] = &[
    ('A', 'T'),
    ('T', 'A'),
    ('C', 'G'),
    ('G', 'C'),
    ('N', 'N'),
];

fn complement_base(base: char) -> char {
    for &(from, to) in COMPLEMENT {
        if from == base {
            return to;
        }
    }
    base
}

fn reverse_complement(sequence: &str) -> String {
    sequence.chars().rev().map(complement_base).collect()
}

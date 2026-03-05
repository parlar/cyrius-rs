use crate::depth_calling::utilities::open_alignment_file;
use crate::stats;
use crate::types::{NormalizedDepth, Region, RegionDic};
use indexmap::IndexMap;
use rust_htslib::bam::{self, Read as BamRead};

/// Return the number of reads that align to a region.
/// Keep duplicate reads. Keep unmapped reads with mapped mates.
pub fn get_read_count(reader: &mut bam::IndexedReader, region: &Region, mapq_cutoff: u8) -> u64 {
    let tid = match reader.header().tid(region.0.as_bytes()) {
        Some(t) => t,
        None => return 0,
    };
    reader
        .fetch(bam::FetchDefinition::Region(
            tid as i32,
            region.1,
            region.2,
        ))
        .unwrap();

    let mut nreads: u64 = 0;
    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        if result.is_err() {
            continue;
        }
        if record.mapq() >= mapq_cutoff
            && !record.is_secondary()
            && !record.is_supplementary()
            && record.pos() >= region.1
            && record.pos() < region.2
        {
            nreads += 1;
        }
    }
    nreads
}

/// Return the median read length from a set of reads (fetches all reads).
pub fn get_read_length(reader: &mut bam::IndexedReader) -> f64 {
    reader.fetch(bam::FetchDefinition::All).unwrap();
    get_read_length_from_fetched(reader)
}

/// Return the median read length from a region.
fn get_read_length_from_region(reader: &mut bam::IndexedReader, region: &Region) -> f64 {
    let tid = match reader.header().tid(region.0.as_bytes()) {
        Some(t) => t,
        None => return 150.0,
    };
    reader
        .fetch(bam::FetchDefinition::Region(tid as i32, region.1, region.2))
        .unwrap();
    get_read_length_from_fetched(reader)
}

fn get_read_length_from_fetched(reader: &mut bam::IndexedReader) -> f64 {
    let mut read_lengths = Vec::new();
    let mut record = bam::Record::new();
    let mut counter = 0;
    while let Some(result) = reader.read(&mut record) {
        if result.is_err() {
            continue;
        }
        read_lengths.push(record.seq_len() as f64);
        counter += 1;
        if counter > 2000 {
            break;
        }
    }
    if read_lengths.is_empty() {
        return 150.0;
    }
    stats::median(&read_lengths)
}

/// GC correction using LOWESS smoothing.
/// Matches Python gc_correction() with scale_coefficient=0.9.
pub fn gc_correction(counts: &[f64], gc: &[f64], scale_coefficient: f64) -> Vec<f64> {
    let med = stats::median(counts);
    if med == 0.0 {
        return counts.to_vec();
    }
    let y_counts: Vec<f64> = counts.iter().map(|&c| c / med).collect();

    // LOWESS smoothing (frac=2/3, default)
    // The lowess crate sorts by x internally; we need to unsort to match
    // Python's return_sorted=False behavior.
    let x_gc: Vec<f64> = gc.to_vec();
    let frac = 2.0 / 3.0;
    let lowess_result = lowess::lowess_with_fraction(&x_gc, &y_counts, frac).unwrap();

    // Build a map from sorted result back to original order.
    // Must use sort_unstable_by to match the lowess crate's internal sort.
    let mut indices: Vec<usize> = (0..x_gc.len()).collect();
    indices.sort_unstable_by(|&a, &b| x_gc[a].partial_cmp(&x_gc[b]).unwrap());
    let mut value_lowess = vec![0.0; x_gc.len()];
    for (sorted_pos, &orig_idx) in indices.iter().enumerate() {
        value_lowess[orig_idx] = lowess_result.y[sorted_pos];
    }

    let sample_median = stats::median(&y_counts);

    let mut gc_corrected = Vec::with_capacity(y_counts.len());
    for i in 0..y_counts.len() {
        let scale_factor = scale_coefficient * y_counts[i].min(2.0);
        gc_corrected.push(y_counts[i] + scale_factor * (sample_median - value_lowess[i]));
    }
    gc_corrected
}

/// Median correction (simple normalization by median).
pub fn median_correction(counts: &[f64]) -> Vec<f64> {
    let med = stats::median(counts);
    if med == 0.0 {
        return counts.to_vec();
    }
    counts.iter().map(|&c| c / med).collect()
}

/// Normalize depth values and return normalized depth, median depth, and MAD.
pub fn normalize(
    counts_for_normalization: &[f64],
    gc_for_normalization: &[f64],
    region_type_cn: &IndexMap<String, f64>,
    read_length: f64,
    gc_correct: bool,
) -> NormalizedDepth {
    let gc_corrected_depth = gc_correction(counts_for_normalization, gc_for_normalization, 0.9);
    let corrected_depth = if gc_correct {
        gc_corrected_depth.clone()
    } else {
        median_correction(counts_for_normalization)
    };

    let vmedian = stats::median(counts_for_normalization) * read_length;
    let gc_med = stats::median(&gc_corrected_depth);
    let vmad = if gc_med == 0.0 {
        0.0
    } else {
        (stats::mad(&gc_corrected_depth) / gc_med * 1000.0).round() / 1000.0
    };

    let mut norm_count = std::collections::HashMap::new();
    for (i, (region_type, &hap_cn)) in region_type_cn.iter().enumerate() {
        if vmedian == 0.0 {
            norm_count.insert(region_type.clone(), None);
        } else {
            norm_count.insert(
                region_type.clone(),
                Some(2.0 * hap_cn * corrected_depth[i]),
            );
        }
    }

    NormalizedDepth {
        normalized: norm_count,
        mediandepth: vmedian,
        mad: vmad,
    }
}

/// Get normalized depth from a BAM file.
pub fn get_normed_depth(
    bam_path: &str,
    region_dic: &RegionDic,
    n_cores: usize,
    reference: Option<&str>,
) -> NormalizedDepth {
    let mut reader = open_alignment_file(bam_path, reference).unwrap();

    let mut counts_for_normalization = Vec::new();
    let mut gc_for_normalization = Vec::new();
    let mut region_type_cn = IndexMap::new();
    let mut last_region: Option<Region> = None;

    for (region_type, regions) in region_dic {
        if region_type == "norm" {
            continue;
        }
        let mut lcount = Vec::new();
        let mut region_length: Option<i64> = None;
        let mut hap_cn: Option<f64> = None;

        for (region, gc) in regions {
            let nreads = get_read_count(&mut reader, region, 0);
            lcount.push(nreads);
            if region.3.contains("_hapcn") {
                region_length = Some(region.2 - region.1);
                gc_for_normalization.push(*gc);
                let cn_str = region.3.split("hapcn").nth(1).unwrap();
                hap_cn = Some(cn_str.parse::<f64>().unwrap());
            }
            last_region = Some(region.clone());
        }

        let region_length = region_length.expect("Problem with region definition. Length not specified.");
        let hap_cn = hap_cn.expect("Problem with region definition. Length not specified.");
        let count_sum: u64 = lcount.iter().sum();
        counts_for_normalization.push(count_sum as f64 / (hap_cn * region_length as f64));
        region_type_cn.insert(region_type.clone(), hap_cn);
    }

    // Get read length from the last region
    let read_length = if let Some(ref region) = last_region {
        get_read_length_from_region(&mut reader, region)
    } else {
        150.0
    };

    // Process normalization regions
    let norm_regions = match region_dic.get("norm") {
        Some(r) => r.clone(),
        None => Vec::new(),
    };

    if n_cores <= 1 {
        for (region, gc) in &norm_regions {
            let nreads = get_read_count(&mut reader, region, 0);
            let region_length = (region.2 - region.1) as f64;
            let norm_depth = nreads as f64 / region_length;
            counts_for_normalization.push(norm_depth);
            gc_for_normalization.push(*gc);
        }
    } else {
        // Parallel processing using rayon
        use rayon::prelude::*;
        let results: Vec<(f64, f64)> = norm_regions
            .par_iter()
            .map(|(region, gc)| {
                let mut r = open_alignment_file(bam_path, reference).unwrap();
                let nreads = get_read_count(&mut r, region, 0);
                let region_length = (region.2 - region.1) as f64;
                (nreads as f64 / region_length, *gc)
            })
            .collect();
        for (norm_depth, gc) in results {
            counts_for_normalization.push(norm_depth);
            gc_for_normalization.push(gc);
        }
    }

    normalize(
        &counts_for_normalization,
        &gc_for_normalization,
        &region_type_cn,
        read_length,
        true,
    )
}

/// Get normalized depth from a count file.
pub fn get_normed_depth_from_count(
    count_file: &str,
    region_dic: &RegionDic,
    read_length: f64,
) -> NormalizedDepth {
    let count_dic = get_count_from_file(count_file);

    let mut counts_for_normalization = Vec::new();
    let mut gc_for_normalization = Vec::new();
    let mut region_type_cn = IndexMap::new();

    for (region_type, regions) in region_dic {
        if region_type == "norm" {
            continue;
        }
        let mut lcount = Vec::new();
        let mut region_length: Option<i64> = None;
        let mut hap_cn: Option<f64> = None;

        for (region, gc) in regions {
            lcount.push(*count_dic.get(&region.3).unwrap_or(&0) as f64);
            if region.3.contains("_hapcn") {
                region_length = Some(region.2 - region.1);
                gc_for_normalization.push(*gc);
                let cn_str = region.3.split("hapcn").nth(1).unwrap();
                hap_cn = Some(cn_str.parse::<f64>().unwrap());
            }
        }

        let region_length = region_length.expect("Problem with region definition.");
        let hap_cn = hap_cn.expect("Problem with region definition.");
        let count_sum: f64 = lcount.iter().sum();
        counts_for_normalization.push(count_sum / (hap_cn * region_length as f64));
        region_type_cn.insert(region_type.clone(), hap_cn);
    }

    if let Some(norm_regions) = region_dic.get("norm") {
        for (region, gc) in norm_regions {
            let region_length = (region.2 - region.1) as f64;
            let count = *count_dic.get(&region.3).unwrap_or(&0) as f64;
            counts_for_normalization.push(count / region_length);
            gc_for_normalization.push(*gc);
        }
    }

    normalize(
        &counts_for_normalization,
        &gc_for_normalization,
        &region_type_cn,
        read_length,
        true,
    )
}

/// Parse count file.
fn get_count_from_file(count_file: &str) -> std::collections::HashMap<String, u64> {
    let mut count_dic = std::collections::HashMap::new();
    let content = std::fs::read_to_string(count_file).unwrap();
    for line in content.lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() >= 5 {
            let name = fields[3].to_string();
            let count: u64 = fields.last().unwrap().parse().unwrap();
            count_dic.entry(name).or_insert(count);
        }
    }
    count_dic
}

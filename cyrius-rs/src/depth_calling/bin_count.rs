use crate::depth_calling::utilities::open_alignment_file_with_index;
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
        stats::python_round3(stats::mad(&gc_corrected_depth) / gc_med)
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

/// A control region must have at least this fraction of the provisional
/// median (computed over the strictly-positive control depths) to be kept.
/// Zero-depth regions always fall below this (the provisional median is
/// computed only from positive depths, so the threshold is strictly
/// positive whenever any region has reads at all).
const MIN_CONTROL_DEPTH_FRACTION_OF_MEDIAN: f64 = 0.10;

/// At least this fraction of control regions must be usable (see
/// [`select_usable_control_regions`]) or normalization is considered
/// unreliable.
const MIN_USABLE_CONTROL_FRACTION: f64 = 0.8;

/// Exclude control ("norm") regions whose depth is zero or negligible
/// (below 10% of the median of the strictly-positive control depths) from
/// the values that feed the normalization median and GC fit.
///
/// A region with literally zero reads is always excluded: the provisional
/// median is computed only over strictly-positive depths, so the 10%
/// threshold derived from it is strictly positive whenever at least one
/// region has reads, and a zero depth can never meet a positive threshold.
/// If every control region is zero, nothing is usable.
///
/// Returns `Err` if fewer than [`MIN_USABLE_CONTROL_FRACTION`] of the
/// control regions remain usable, since normalization on so few good
/// regions is not reliable.
///
/// Shared by both [`get_normed_depth`] (BAM path) and
/// [`get_normed_depth_from_count`] (count-file path) so the two cannot
/// drift apart.
pub fn select_usable_control_regions(
    control_counts: &[f64],
    control_gc: &[f64],
) -> Result<(Vec<f64>, Vec<f64>), String> {
    debug_assert_eq!(control_counts.len(), control_gc.len());
    let total = control_counts.len();
    if total == 0 {
        return Ok((Vec::new(), Vec::new()));
    }

    let positive_depths: Vec<f64> = control_counts.iter().copied().filter(|&c| c > 0.0).collect();

    let mut usable_counts = Vec::new();
    let mut usable_gc = Vec::new();
    if !positive_depths.is_empty() {
        let threshold = MIN_CONTROL_DEPTH_FRACTION_OF_MEDIAN * stats::median(&positive_depths);
        for (&c, &g) in control_counts.iter().zip(control_gc.iter()) {
            if c >= threshold {
                usable_counts.push(c);
                usable_gc.push(g);
            }
        }
    }
    // else: every control region is zero; none are usable.

    let usable_fraction = usable_counts.len() as f64 / total as f64;
    if usable_fraction < MIN_USABLE_CONTROL_FRACTION {
        return Err(format!(
            "only {} of {} control regions usable ({:.1}% < {:.0}% required); \
             depth normalization is unreliable",
            usable_counts.len(),
            total,
            usable_fraction * 100.0,
            MIN_USABLE_CONTROL_FRACTION * 100.0
        ));
    }

    Ok((usable_counts, usable_gc))
}

/// If `coverage_mad` exceeds `threshold` and `filter` is a clean `"PASS"`,
/// downgrade it to `"LowQ_coverage_MAD"` so the Filter field reflects the
/// coverage-quality concern instead of only logging a warning. Any other
/// filter value (ambiguous calls, fuzzy matches, no-calls, ...) is left
/// untouched, since those already explain why the call is not a clean pass.
pub fn downgrade_filter_for_high_mad(
    filter: Option<String>,
    coverage_mad: f64,
    threshold: f64,
) -> Option<String> {
    if coverage_mad > threshold && filter.as_deref() == Some("PASS") {
        Some("LowQ_coverage_MAD".to_string())
    } else {
        filter
    }
}

/// Get normalized depth from a BAM file. Returns Err if the alignment file
/// (or its index) cannot be opened.
pub fn get_normed_depth(
    bam_path: &str,
    region_dic: &RegionDic,
    n_cores: usize,
    reference: Option<&str>,
    index_name: Option<&str>,
) -> Result<NormalizedDepth, String> {
    let open = || {
        open_alignment_file_with_index(bam_path, reference, index_name)
            .map_err(|e| format!("cannot open alignment file {}: {}", bam_path, e))
    };
    let mut reader = open()?;

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
                let cn_str = region.3.split("hapcn").nth(1).unwrap_or("");
                hap_cn = cn_str.parse::<f64>().ok();
            }
            last_region = Some(region.clone());
        }

        let region_length = region_length
            .expect("Problem with region definition: no '_hapcn' region length found.");
        let hap_cn = hap_cn
            .expect("Problem with region definition: '_hapcn' copy number missing or unparseable.");
        let count_sum: u64 = lcount.iter().sum();
        // Guard against a degenerate region (zero length or zero hap_cn) that
        // would produce inf/NaN and silently corrupt downstream normalization.
        let denom = hap_cn * region_length as f64;
        if denom <= 0.0 {
            panic!(
                "Problem with region definition: non-positive normalization denominator \
                 (hap_cn={}, region_length={}).",
                hap_cn, region_length
            );
        }
        counts_for_normalization.push(count_sum as f64 / denom);
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

    let mut control_counts = Vec::new();
    let mut control_gc = Vec::new();

    if n_cores <= 1 {
        for (region, gc) in &norm_regions {
            let nreads = get_read_count(&mut reader, region, 0);
            let region_length = (region.2 - region.1) as f64;
            let norm_depth = nreads as f64 / region_length;
            control_counts.push(norm_depth);
            control_gc.push(*gc);
        }
    } else {
        // Parallel processing using rayon
        use rayon::prelude::*;
        let results: Result<Vec<(f64, f64)>, String> = norm_regions
            .par_iter()
            .map(|(region, gc)| {
                let mut r = open()?;
                let nreads = get_read_count(&mut r, region, 0);
                let region_length = (region.2 - region.1) as f64;
                Ok((nreads as f64 / region_length, *gc))
            })
            .collect();
        for (norm_depth, gc) in results? {
            control_counts.push(norm_depth);
            control_gc.push(gc);
        }
    }

    // Drop zero/negligible-depth control regions before they can pull the
    // normalization median and GC fit off; Err if too few remain usable.
    let (usable_counts, usable_gc) = select_usable_control_regions(&control_counts, &control_gc)?;
    counts_for_normalization.extend(usable_counts);
    gc_for_normalization.extend(usable_gc);

    Ok(normalize(
        &counts_for_normalization,
        &gc_for_normalization,
        &region_type_cn,
        read_length,
        true,
    ))
}

/// Get normalized depth from a count file. Returns Err if the count file
/// cannot be read or parsed.
pub fn get_normed_depth_from_count(
    count_file: &str,
    region_dic: &RegionDic,
    read_length: f64,
) -> Result<NormalizedDepth, String> {
    let count_dic = get_count_from_file(count_file)?;

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

    let mut control_counts = Vec::new();
    let mut control_gc = Vec::new();
    if let Some(norm_regions) = region_dic.get("norm") {
        for (region, gc) in norm_regions {
            let region_length = (region.2 - region.1) as f64;
            let count = *count_dic.get(&region.3).unwrap_or(&0) as f64;
            control_counts.push(count / region_length);
            control_gc.push(*gc);
        }
    }

    // Drop zero/negligible-depth control regions before they can pull the
    // normalization median and GC fit off; Err if too few remain usable.
    // Shared with get_normed_depth (BAM path) so the two cannot drift.
    let (usable_counts, usable_gc) = select_usable_control_regions(&control_counts, &control_gc)?;
    counts_for_normalization.extend(usable_counts);
    gc_for_normalization.extend(usable_gc);

    Ok(normalize(
        &counts_for_normalization,
        &gc_for_normalization,
        &region_type_cn,
        read_length,
        true,
    ))
}

/// Parse count file.
fn get_count_from_file(count_file: &str) -> Result<std::collections::HashMap<String, u64>, String> {
    let mut count_dic = std::collections::HashMap::new();
    let content = std::fs::read_to_string(count_file)
        .map_err(|e| format!("cannot read count file {}: {}", count_file, e))?;
    for line in content.lines() {
        let trimmed = line.trim_start();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() >= 5 {
            let name = fields[3].to_string();
            let last = fields[fields.len() - 1];
            let count: u64 = last
                .parse()
                .map_err(|_| format!("bad count '{}' in count file {}", last, count_file))?;
            count_dic.entry(name).or_insert(count);
        }
    }
    Ok(count_dic)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_bam_is_an_error_not_a_panic() {
        let r = get_normed_depth("/nonexistent/x.bam", &RegionDic::new(), 1, None, None);
        assert!(r.is_err());
    }

    #[test]
    fn missing_count_file_is_an_error_not_a_panic() {
        let r = get_normed_depth_from_count("/nonexistent/x_count.txt", &RegionDic::new(), 150.0);
        assert!(r.is_err());
    }

    #[test]
    fn custom_index_path_is_used() {
        // A symlinked BAM with no index next to it can only be opened through
        // the custom index. Skipped when the benchmark data is not present.
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let bam = root.join("benchmark/HG00096_pgx_regions.bam");
        let bai = root.join("benchmark/HG00096_pgx_regions.bam.bai");
        if !bam.exists() || !bai.exists() {
            eprintln!("skipping: benchmark BAM not present");
            return;
        }
        let dir = std::env::temp_dir().join(format!("cyrius_idx_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let link = dir.join("sample.bam");
        let _ = std::fs::remove_file(&link);
        std::os::unix::fs::symlink(&bam, &link).unwrap();
        let regions = crate::depth_calling::utilities::parse_region_file(crate::data::REGION_BED, "chr38");
        let r = get_normed_depth(link.to_str().unwrap(), &regions, 1, None, bai.to_str());
        let _ = std::fs::remove_dir_all(&dir);
        assert!(r.is_ok(), "{:?}", r.err());
    }

    #[test]
    fn select_usable_control_regions_excludes_zero_and_near_zero_depths() {
        // 9 healthy regions at depth 20.0, one zeroed out and one at 1.0
        // (below 10% of the median of 20.0, i.e. below 2.0). Both should be
        // dropped; 9/11 = 81.8% usable clears the 80% floor.
        let counts = vec![
            20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 0.0, 1.0,
        ];
        let gc = vec![
            0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.11, 0.12,
        ];
        let (kept_counts, kept_gc) = select_usable_control_regions(&counts, &gc).unwrap();
        assert_eq!(kept_counts.len(), 9);
        assert!(kept_counts.iter().all(|&c| c == 20.0));
        assert_eq!(kept_gc.len(), 9);
        // The zero and near-zero region's GC values must not survive.
        assert!(!kept_gc.contains(&0.11));
        assert!(!kept_gc.contains(&0.12));
        // The zero/near-zero entries no longer drag the median down: the
        // median of the kept values is unaffected by them.
        assert_eq!(stats::median(&kept_counts), 20.0);
    }

    #[test]
    fn select_usable_control_regions_all_zero_is_err() {
        let counts = vec![0.0; 10];
        let gc = vec![0.5; 10];
        let r = select_usable_control_regions(&counts, &gc);
        assert!(r.is_err(), "{:?}", r);
    }

    #[test]
    fn select_usable_control_regions_below_80_percent_usable_is_err() {
        // 7 usable of 10 = 70%, below the 80% floor.
        let mut counts = vec![10.0; 7];
        counts.extend(vec![0.0; 3]);
        let gc = vec![0.5; 10];
        let r = select_usable_control_regions(&counts, &gc);
        assert!(r.is_err(), "{:?}", r);
    }

    #[test]
    fn select_usable_control_regions_at_80_percent_is_ok() {
        // 8 usable of 10 = exactly 80%, at the floor (not below it).
        let mut counts = vec![10.0; 8];
        counts.extend(vec![0.0; 2]);
        let gc = vec![0.5; 10];
        let (kept_counts, _) = select_usable_control_regions(&counts, &gc).unwrap();
        assert_eq!(kept_counts.len(), 8);
    }

    #[test]
    fn select_usable_control_regions_empty_input_is_ok() {
        let (counts, gc) = select_usable_control_regions(&[], &[]).unwrap();
        assert!(counts.is_empty());
        assert!(gc.is_empty());
    }

    #[test]
    fn downgrade_filter_for_high_mad_replaces_pass_above_threshold() {
        let r = downgrade_filter_for_high_mad(Some("PASS".to_string()), 0.27, 0.11);
        assert_eq!(r, Some("LowQ_coverage_MAD".to_string()));
    }

    #[test]
    fn downgrade_filter_for_high_mad_leaves_pass_alone_below_threshold() {
        let r = downgrade_filter_for_high_mad(Some("PASS".to_string()), 0.077, 0.11);
        assert_eq!(r, Some("PASS".to_string()));
    }

    #[test]
    fn downgrade_filter_for_high_mad_leaves_non_pass_filters_alone() {
        // A sample already flagged for another reason keeps that reason;
        // it is not clobbered by the MAD check.
        let r = downgrade_filter_for_high_mad(
            Some("Fuzzy_match".to_string()),
            0.27,
            0.11,
        );
        assert_eq!(r, Some("Fuzzy_match".to_string()));
    }

    #[test]
    fn downgrade_filter_for_high_mad_leaves_none_alone() {
        let r = downgrade_filter_for_high_mad(None, 0.27, 0.11);
        assert_eq!(r, None);
    }

    #[test]
    fn count_file_with_comment_header_and_blank_lines_parses() {
        // Reproduces finding #4: a count file whose header is a '#'-prefixed
        // line with >=5 whitespace-separated fields (e.g.
        // '#chr start end name type gc count') used to be read as a data
        // row, failing to parse its last field ('count') as a u64.
        let dir = std::env::temp_dir().join(format!("cyrius_count_header_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("sample_count.txt");
        std::fs::write(
            &path,
            "#chr\tstart\tend\tname\ttype\tgc\tcount\n\
             \n\
             chr1\t100\t200\tregion_a\tnorm\t0.5\t42\n\
             \n\
             chr1\t300\t400\tregion_b\tnorm\t0.6\t7\n",
        )
        .unwrap();

        let count_dic = get_count_from_file(path.to_str().unwrap())
            .expect("a count file with a '#' header line and blank lines should parse");

        assert_eq!(count_dic.get("region_a"), Some(&42));
        assert_eq!(count_dic.get("region_b"), Some(&7));
        // The comment header's literal field value must not leak in as data.
        assert_eq!(count_dic.get("name"), None);
        assert_eq!(count_dic.len(), 2);

        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_dir(&dir);
    }

    #[test]
    fn missing_norm_region_name_is_excluded_like_the_region_never_existed() {
        // Task 1's select_usable_control_regions (shared by both
        // normalization paths) excludes zero/negligible-depth control
        // regions before they can skew the median. This exercises the full
        // exclusion step end-to-end via get_normed_depth_from_count, across
        // three scenarios for one control region out of ten:
        //   1. explicit_zero: the region is present with an explicit count
        //      of 0.
        //   2. missing: the region's name is listed in the norm RegionDic
        //      but absent from the count file, which reads as 0 via
        //      unwrap_or(&0) in the `norm` loop in
        //      get_normed_depth_from_count -- the same 0 as scenario 1 by
        //      the time it reaches select_usable_control_regions.
        //   3. no_zero_region: the region is not in the norm RegionDic at
        //      all, i.e. only the other 9 (real, positive-depth) control
        //      regions exist. This arm is falsifiable: if the exclusion
        //      step were broken and the zero from scenarios 1/2 survived
        //      into their normalization instead of being dropped, those
        //      would fold a 10th (zero) region into the median and GC fit
        //      that this 9-region scenario never sees, and the equality
        //      assertions below would fail.
        // All three must produce identical results if the zero/missing
        // region is truly excluded before use, as Task 1 established.
        let mut region_dic: RegionDic = IndexMap::new();
        region_dic.insert(
            "CYP2D6".to_string(),
            vec![(
                ("chr1".to_string(), 0, 1000, "CYP2D6_hapcn2".to_string()),
                0.5,
            )],
        );

        let gcs = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70];
        let counts: [u64; 9] = [90, 95, 100, 105, 110, 100, 95, 105, 100];
        let mut norm_regions = Vec::new();
        let mut good_lines = String::new();
        for (i, (&gc, &count)) in gcs.iter().zip(counts.iter()).enumerate() {
            let name = format!("norm_{}", i);
            norm_regions.push((("chr2".to_string(), 0, 1000, name.clone()), gc));
            good_lines.push_str(&format!("chr2\t0\t1000\t{}\tnorm\t{}\t{}\n", name, gc, count));
        }

        // Scenario 3's RegionDic: only the 9 real positive-depth regions,
        // no zero-depth region entry at all.
        let mut region_dic_no_zero_region = region_dic.clone();
        region_dic_no_zero_region.insert("norm".to_string(), norm_regions.clone());

        // Scenarios 1 and 2's RegionDic: the same 9 regions plus a 10th,
        // zero-depth region ("norm_var") that must be excluded.
        norm_regions.push((("chr2".to_string(), 0, 1000, "norm_var".to_string()), 0.5));
        region_dic.insert("norm".to_string(), norm_regions);

        let base = format!("chr1\t0\t1000\tCYP2D6_hapcn2\tCYP2D6\t0.5\t200\n{}", good_lines);

        let dir = std::env::temp_dir().join(format!("cyrius_missing_norm_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        let explicit_zero_path = dir.join("explicit_zero_count.txt");
        std::fs::write(
            &explicit_zero_path,
            format!("{}chr2\t0\t1000\tnorm_var\tnorm\t0.5\t0\n", base),
        )
        .unwrap();

        let missing_path = dir.join("missing_count.txt");
        std::fs::write(&missing_path, &base).unwrap();

        let no_zero_region_path = dir.join("no_zero_region_count.txt");
        std::fs::write(&no_zero_region_path, &base).unwrap();

        let explicit_zero =
            get_normed_depth_from_count(explicit_zero_path.to_str().unwrap(), &region_dic, 150.0)
                .expect("explicit zero count should normalize successfully");
        let missing = get_normed_depth_from_count(missing_path.to_str().unwrap(), &region_dic, 150.0)
            .expect("a missing region name should normalize successfully, the same as an explicit zero");
        let no_zero_region = get_normed_depth_from_count(
            no_zero_region_path.to_str().unwrap(),
            &region_dic_no_zero_region,
            150.0,
        )
        .expect("a RegionDic with no zero-depth region at all should normalize successfully");

        assert_eq!(explicit_zero.mediandepth, missing.mediandepth);
        assert_eq!(explicit_zero.mad, missing.mad);
        assert_eq!(explicit_zero.normalized, missing.normalized);

        assert_eq!(
            explicit_zero.mediandepth, no_zero_region.mediandepth,
            "a region absent from the RegionDic entirely must normalize the same as one \
             explicitly zeroed -- if this fails, the zero region leaked into the median/GC fit"
        );
        assert_eq!(explicit_zero.mad, no_zero_region.mad);
        assert_eq!(explicit_zero.normalized, no_zero_region.normalized);

        let _ = std::fs::remove_file(&explicit_zero_path);
        let _ = std::fs::remove_file(&missing_path);
        let _ = std::fs::remove_file(&no_zero_region_path);
        let _ = std::fs::remove_dir(&dir);
    }
}

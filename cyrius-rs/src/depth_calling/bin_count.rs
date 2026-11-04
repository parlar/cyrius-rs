use crate::depth_calling::utilities::{fetch_region, open_alignment_file_with_index};
use crate::stats;
use crate::types::{NormalizedDepth, Region, RegionDic};
use indexmap::IndexMap;
use rust_htslib::bam::{self, Read as BamRead};

/// Return the number of reads that align to a region.
/// Keep duplicate reads. Keep unmapped reads with mapped mates.
pub fn get_read_count(
    reader: &mut bam::IndexedReader,
    region: &Region,
    mapq_cutoff: u8,
) -> Result<u64, String> {
    let tid = match reader.header().tid(region.0.as_bytes()) {
        Some(t) => t,
        None => return Ok(0),
    };
    let what = format!("{}:{}-{}", region.0, region.1, region.2);
    fetch_region(reader, tid, region.1, region.2, &what)?;

    let mut nreads: u64 = 0;
    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        if let Err(e) = result {
            return Err(format!("BAM read failed in {}: {}", what, e));
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
    Ok(nreads)
}

/// Return the median read length from a set of reads (fetches all reads).
pub fn get_read_length(reader: &mut bam::IndexedReader) -> Result<f64, String> {
    // Unlike a region fetch, FetchDefinition::All cannot fail for a missing
    // contig -- it has no contig to look up. It can still fail on a broken
    // or stale index; like the region fetches, that is an error (E1), not a
    // guessed read length.
    reader
        .fetch(bam::FetchDefinition::All)
        .map_err(|e| format!("BAM fetch failed for all reads: {e}"))?;
    get_read_length_from_fetched(reader)
}

/// Return the median read length from a region.
fn get_read_length_from_region(
    reader: &mut bam::IndexedReader,
    region: &Region,
) -> Result<f64, String> {
    let tid = match reader.header().tid(region.0.as_bytes()) {
        Some(t) => t,
        None => return Ok(150.0),
    };
    let what = format!("{}:{}-{}", region.0, region.1, region.2);
    fetch_region(reader, tid, region.1, region.2, &what)?;
    get_read_length_from_fetched(reader)
}

fn get_read_length_from_fetched(reader: &mut bam::IndexedReader) -> Result<f64, String> {
    let mut read_lengths = Vec::new();
    let mut record = bam::Record::new();
    let mut counter = 0;
    while let Some(result) = reader.read(&mut record) {
        if let Err(e) = result {
            return Err(format!("BAM read failed while estimating read length: {}", e));
        }
        read_lengths.push(record.seq_len() as f64);
        counter += 1;
        if counter > 2000 {
            break;
        }
    }
    if read_lengths.is_empty() {
        return Ok(150.0);
    }
    Ok(stats::median(&read_lengths))
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

/// Normalize depth values and return normalized depth, median depth, and MAD.
/// Depths are always GC-corrected.
pub fn normalize(
    counts_for_normalization: &[f64],
    gc_for_normalization: &[f64],
    region_type_cn: &IndexMap<String, f64>,
    read_length: f64,
) -> NormalizedDepth {
    let gc_corrected_depth = gc_correction(counts_for_normalization, gc_for_normalization, 0.9);

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
                Some(2.0 * hap_cn * gc_corrected_depth[i]),
            );
        }
    }

    NormalizedDepth {
        normalized: norm_count,
        mediandepth: vmedian,
        mad: vmad,
        rep7_mapq0_fraction: None,
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
/// (mapped reads, MAPQ-0 reads) among the reads that start in `region`:
/// primary alignments (not secondary, not supplementary), mapped, 0-based
/// start inside the region; duplicates included, as in `get_read_count`.
/// A missing contig gives (0, 0); a failed fetch is an error (E1).
pub fn get_mapq0_counts(
    reader: &mut bam::IndexedReader,
    region: &Region,
) -> Result<(u64, u64), String> {
    let tid = match reader.header().tid(region.0.as_bytes()) {
        Some(t) => t,
        None => return Ok((0, 0)),
    };
    let what = format!("{}:{}-{}", region.0, region.1, region.2);
    fetch_region(reader, tid, region.1, region.2, &what)?;

    let (mut mapped, mut mapq0) = (0u64, 0u64);
    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        if let Err(e) = result {
            return Err(format!("BAM read failed in {}: {}", what, e));
        }
        if record.is_unmapped()
            || record.is_secondary()
            || record.is_supplementary()
            || record.pos() < region.1
            || record.pos() >= region.2
        {
            continue;
        }
        mapped += 1;
        if record.mapq() == 0 {
            mapq0 += 1;
        }
    }
    Ok((mapped, mapq0))
}

/// `mapq0 / mapped` from `get_mapq0_counts`, or None when no mapped read
/// starts in the region.
pub fn mapq0_fraction((mapped, mapq0): (u64, u64)) -> Option<f64> {
    (mapped > 0).then(|| mapq0 as f64 / mapped as f64)
}

/// cyrius-rs reads copy number from reads that fit CYP2D6 and CYP2D7 equally
/// well, which aligners report with MAPQ 0. An `@PG` line in the header shows
/// a pipeline step that drops them: novoalign without `-r` (GIAB's 3.02.07
/// BAMs: multi-mapped reads are not reported), or `samtools view -q N` with
/// N >= 1. Their calls can be wrong at PASS.
///
/// This header check, not the REP7 MAPQ-0 fraction, drives the filter: in
/// the guard experiment (/home/parlar_ai/cyrius-research/exp-guard, stage 2,
/// result e70d315) the fraction flagged 3 correct calls from pipelines that
/// reassign or rescore ties (10x Long Ranger, DRAGEN 4.4.7, vg giraffe) and
/// caught nothing the header missed. The fraction is still reported.
pub fn header_drops_multimappers(header_text: &str) -> bool {
    let no_repeat_option = regex::Regex::new(r"\s-r\s").unwrap();
    let mapq_filter = regex::Regex::new(r"\s-q\s*[1-9]\d*").unwrap();
    header_text.lines().filter(|l| l.starts_with("@PG\t")).any(|line| {
        let field = |key: &str| {
            line.split('\t')
                .skip(1)
                .find_map(|f| f.strip_prefix(key).and_then(|v| v.strip_prefix(':')))
                .unwrap_or("")
        };
        let name = format!("{} {}", field("ID"), field("PN")).to_lowercase();
        let cl = field("CL");
        (name.contains("novoalign") && !no_repeat_option.is_match(cl))
            || (cl.contains("samtools view") && mapq_filter.is_match(cl))
    })
}

/// Replace a PASS, Fuzzy_match, Unknown_sites_unmatched or Unknown_sites_resolved filter with `"Multimapper_depleted"` when
/// the header shows a step that drops multi-mapped reads
/// ([`header_drops_multimappers`]): the input does not meet the caller's
/// assumption, so the call must not read as a clean pass. The genotype is
/// still reported. Header not read (None), other filters and no-calls are
/// left untouched.
pub fn downgrade_filter_for_multimapper_depletion(
    filter: Option<String>,
    header_drops_multimappers: Option<bool>,
) -> Option<String> {
    match (filter.as_deref(), header_drops_multimappers) {
        (Some("PASS" | "Fuzzy_match" | "Unknown_sites_unmatched" | "Unknown_sites_resolved"), Some(true)) => Some("Multimapper_depleted".to_string()),
        _ => filter,
    }
}

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
    let mut rep7_mapq0_fraction = None;

    for (region_type, regions) in region_dic {
        if region_type == "norm" {
            continue;
        }
        let mut lcount = Vec::new();
        let mut region_length: Option<i64> = None;
        let mut hap_cn: Option<f64> = None;

        for (region, gc) in regions {
            let nreads = get_read_count(&mut reader, region, 0)?;
            lcount.push(nreads);
            if region.3 == "REP7" {
                let counts = get_mapq0_counts(&mut reader, region)?;
                log::info!("REP7: {} mapped reads, {} with MAPQ 0", counts.0, counts.1);
                rep7_mapq0_fraction = mapq0_fraction(counts);
            }
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
        get_read_length_from_region(&mut reader, region)?
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
            let nreads = get_read_count(&mut reader, region, 0)?;
            let region_length = (region.2 - region.1) as f64;
            let norm_depth = nreads as f64 / region_length;
            control_counts.push(norm_depth);
            control_gc.push(*gc);
        }
    } else {
        // Parallel processing using rayon. One reader is opened per *chunk*,
        // not per region: `get_read_count` fetches its region and drains it,
        // so it leaves no reader state behind -- the serial branch above
        // already reuses a single reader for all 50 regions. Each open
        // re-reads the index and, for CRAM, the reference, so chunking cuts
        // them where it can: `strace -e trace=openat` counts 60 opens of the
        // BAM (each paired with one of its index) at every `-t >= 2` with
        // one-reader-per-region, and 12 / 14 / 18 / 23 at
        // `-t` 2 / 4 / 8 / 16 with this chunking. The saving runs out above
        // ~50 threads, where `chunk_len` saturates at 1 and the count returns
        // to 60 (measured identical to the unchunked loop at `-t 50` and
        // `-t 72`); the win is at the thread counts people actually pass.
        //
        // The property that matters for correctness is that the chunks are a
        // *partition*: every control region is counted exactly once, none
        // dropped and none duplicated. `normalize` reduces the (depth, gc)
        // pairs order-independently, so it is invariant under a consistent
        // permutation of them -- the pin test below stays green if the whole
        // sequence is reversed, and that is fine. What it does catch is a
        // broken partition: `par_chunks_exact` (drops the final short chunk)
        // fails it on `mediandepth`, and iterating `norm_regions` instead of
        // `chunk` (counts every region in every chunk) fails it on `mad`.
        // Review finding #10.
        use rayon::prelude::*;
        let chunk_len = norm_regions.len().div_ceil(n_cores).max(1);
        let results: Result<Vec<Vec<(f64, f64)>>, String> = norm_regions
            .par_chunks(chunk_len)
            .map(|chunk| {
                let mut r = open()?;
                chunk
                    .iter()
                    .map(|(region, gc)| {
                        let nreads = get_read_count(&mut r, region, 0)?;
                        let region_length = (region.2 - region.1) as f64;
                        Ok((nreads as f64 / region_length, *gc))
                    })
                    .collect::<Result<Vec<_>, String>>()
            })
            .collect();
        for (norm_depth, gc) in results?.into_iter().flatten() {
            control_counts.push(norm_depth);
            control_gc.push(gc);
        }
    }

    // Drop zero/negligible-depth control regions before they can pull the
    // normalization median and GC fit off; Err if too few remain usable.
    let (usable_counts, usable_gc) = select_usable_control_regions(&control_counts, &control_gc)?;
    counts_for_normalization.extend(usable_counts);
    gc_for_normalization.extend(usable_gc);

    let mut normalized = normalize(
        &counts_for_normalization,
        &gc_for_normalization,
        &region_type_cn,
        read_length,
    );
    normalized.rep7_mapq0_fraction = rep7_mapq0_fraction;
    Ok(normalized)
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

    /// An indexed BAM of 2,000 mapped 100 bp reads on chr22, starting at 0, 2, 4, ...
    fn write_test_bam(path: &std::path::Path) {
        let mut header = bam::header::Header::new();
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr22");
        sq.push_tag(b"LN", 50_818_468);
        header.push_record(&sq);
        {
            let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
            let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(100)]);
            for i in 0..2000i64 {
                let seq: Vec<u8> = (0..100).map(|j| b"ACGT"[((i + j) % 4) as usize]).collect();
                let mut r = bam::Record::new();
                r.set(format!("r{i}").as_bytes(), Some(&cigar), &seq, &[30u8; 100]);
                r.set_tid(0);
                r.set_pos(i * 2);
                r.set_mapq(60);
                r.set_flags(0);
                writer.write(&r).unwrap();
            }
        }
        bam::index::build(path, None::<&std::path::Path>, bam::index::Type::Bai, 1).unwrap();
    }

    /// Flip one byte of the CRC32 of BGZF block `block` (block 0 holds the
    /// header). Offsets and the index stay valid, so a fetch succeeds and the
    /// failure only shows up while the records are read.
    fn corrupt_bgzf_block(path: &std::path::Path, block: usize) {
        let mut bytes = std::fs::read(path).unwrap();
        let mut off = 0usize;
        for _ in 0..block {
            off += u16::from_le_bytes([bytes[off + 16], bytes[off + 17]]) as usize + 1;
        }
        let size = u16::from_le_bytes([bytes[off + 16], bytes[off + 17]]) as usize + 1;
        assert!(size > 28, "block {block} is the EOF marker; write more reads");
        bytes[off + size - 8] ^= 0xFF;
        std::fs::write(path, bytes).unwrap();
    }

    fn test_bam_dir(tag: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("cyrius_bamio_{}_{}", tag, std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn a_corrupt_block_is_an_error_not_a_short_count() {
        let dir = test_bam_dir("corrupt");
        let path = dir.join("corrupt.bam");
        write_test_bam(&path);
        corrupt_bgzf_block(&path, 2);
        let region: Region = ("chr22".to_string(), 0, 5000, "test".to_string());
        let count = get_read_count(&mut bam::IndexedReader::from_path(&path).unwrap(), &region, 0);
        let mapq0 = get_mapq0_counts(&mut bam::IndexedReader::from_path(&path).unwrap(), &region);
        let length = get_read_length(&mut bam::IndexedReader::from_path(&path).unwrap());
        let _ = std::fs::remove_dir_all(&dir);
        assert!(count.is_err(), "get_read_count gave {count:?}");
        assert!(mapq0.is_err(), "get_mapq0_counts gave {mapq0:?}");
        assert!(length.is_err(), "get_read_length gave {length:?}");
    }

    #[test]
    fn the_same_bam_uncorrupted_counts_every_read() {
        // Control: the test above must fail because of the corrupt block.
        let dir = test_bam_dir("intact");
        let path = dir.join("ok.bam");
        write_test_bam(&path);
        let region: Region = ("chr22".to_string(), 0, 5000, "test".to_string());
        let count = get_read_count(&mut bam::IndexedReader::from_path(&path).unwrap(), &region, 0);
        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(count, Ok(2000));
    }

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
    fn parallel_control_regions_are_bit_identical_to_the_serial_path() {
        // The n_cores > 1 branch reuses one reader per chunk instead of
        // opening one per region; this pins that it produces exactly the same
        // normalization as the single-reader serial branch, bit for bit.
        // Skipped when the benchmark data is not present.
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let bam = root.join("benchmark/HG00096_pgx_regions.bam");
        if !bam.exists() {
            eprintln!("skipping: benchmark BAM not present");
            return;
        }
        let bam = bam.to_str().unwrap();
        let regions =
            crate::depth_calling::utilities::parse_region_file(crate::data::REGION_BED, "chr38");
        let serial = get_normed_depth(bam, &regions, 1, None, None).unwrap();
        for n_cores in [2usize, 4, 7] {
            let parallel = get_normed_depth(bam, &regions, n_cores, None, None).unwrap();
            assert_eq!(
                serial.mediandepth.to_bits(),
                parallel.mediandepth.to_bits(),
                "mediandepth differs at n_cores={}",
                n_cores
            );
            assert_eq!(
                serial.mad.to_bits(),
                parallel.mad.to_bits(),
                "mad differs at n_cores={}",
                n_cores
            );
            assert_eq!(
                serial.normalized.len(),
                parallel.normalized.len(),
                "region set differs at n_cores={}",
                n_cores
            );
            for (region, s_depth) in &serial.normalized {
                let p_depth = parallel
                    .normalized
                    .get(region)
                    .unwrap_or_else(|| panic!("region {} missing at n_cores={}", region, n_cores));
                assert_eq!(
                    s_depth.map(f64::to_bits),
                    p_depth.map(f64::to_bits),
                    "region {} differs at n_cores={}",
                    region,
                    n_cores
                );
            }
        }
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

    /// Write a minimal, structurally-valid 2-reference BAI index (BAI magic,
    /// n_ref=2, no bins/intervals for either ref, n_no_coor=0) to `path`.
    /// Pairing it with a real multi-contig BAM (via an explicit index path)
    /// reproduces a genuine stale/corrupt-index condition: the index loads
    /// fine, but it covers far fewer references than the BAM's own header.
    fn write_minimal_two_contig_bai(path: &std::path::Path) {
        let bai_bytes: [u8; 32] = [
            0x42, 0x41, 0x49, 0x01, // "BAI\1"
            0x02, 0x00, 0x00, 0x00, // n_ref = 2
            0x00, 0x00, 0x00, 0x00, // ref0 n_bin = 0
            0x00, 0x00, 0x00, 0x00, // ref0 n_intv = 0
            0x00, 0x00, 0x00, 0x00, // ref1 n_bin = 0
            0x00, 0x00, 0x00, 0x00, // ref1 n_intv = 0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_no_coor = 0
        ];
        std::fs::write(path, bai_bytes).unwrap();
    }

    #[test]
    fn get_read_length_with_broken_index_is_an_error() {
        // get_read_length() fetches FetchDefinition::All, not a region, so it
        // cannot fail the way a missing-contig region fetch does. But a broken
        // index can still make the underlying htslib fetch return Err: pairing
        // a real multi-contig benchmark BAM with a minimal, mismatched
        // 2-reference index (valid BAI, but covering far fewer references
        // than the BAM's own header) makes sam_itr_querys(".") return NULL.
        // Confirmed empirically before writing this test: calling
        // reader.fetch(FetchDefinition::All) directly against this exact
        // reader/index pairing returns Err(Fetch).
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let bam = root.join("benchmark/HG00096_pgx_regions.bam");
        if !bam.exists() {
            eprintln!("skipping: benchmark BAM not present");
            return;
        }
        let dir =
            std::env::temp_dir().join(format!("cyrius_broken_idx_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bai_path = dir.join("mismatched.bai");
        write_minimal_two_contig_bai(&bai_path);

        let mut reader = crate::depth_calling::utilities::open_alignment_file_with_index(
            bam.to_str().unwrap(),
            None,
            Some(bai_path.to_str().unwrap()),
        )
        .expect("opening with an explicit (if mismatched) index path should succeed");

        let read_len = get_read_length(&mut reader);

        let _ = std::fs::remove_file(&bai_path);
        let _ = std::fs::remove_dir(&dir);

        // Task 4 made this degrade to the 150.0 fallback instead of panicking;
        // E1 (see below) makes it an error: a guessed read length on a BAM
        // whose index cannot be used is not a basis for a call.
        assert!(
            read_len.is_err(),
            "a fetch failure (broken/mismatched index) must be an error, not a \
             fallback read length: got {:?}",
            read_len
        );
    }

    // ---- E1: a failed BAM fetch is a no-call, not zero reads ----
    //
    // Escalation E1 of review batch 4, decided by the human 2026-10-07: when
    // htslib cannot position the reader (broken or stale index, CRAM reference
    // that cannot be resolved, invalid interval) the reads were never looked
    // at. Returning zero counts makes that look like evidence, and the sample
    // gets a confident call. It must surface as an Err, which main.rs turns
    // into a no-call. An inverted interval is the input used here because it
    // is the one that reliably drives `fetch` into Err through the public API
    // (measured: Region(tid, 100, 50) errs; an out-of-range tid or a negative
    // start does not, they just return no reads).

    fn benchmark_bam() -> Option<std::path::PathBuf> {
        let bam = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../benchmark/HG00096_pgx_regions.bam");
        if bam.exists() {
            Some(bam)
        } else {
            eprintln!("skipping: benchmark BAM not present");
            None
        }
    }

    #[test]
    fn get_read_count_is_an_error_when_the_fetch_fails() {
        let Some(bam) = benchmark_bam() else { return };
        let mut reader = bam::IndexedReader::from_path(&bam).unwrap();
        let inverted: Region = ("chr22".to_string(), 42_130_000, 42_129_000, "inverted".to_string());
        let err = get_read_count(&mut reader, &inverted, 0)
            .expect_err("a failed fetch must be an error, not Ok(0)");
        assert!(err.contains("chr22:42130000-42129000"), "error must name the region: {err}");
    }

    #[test]
    fn get_normed_depth_is_an_error_when_one_control_region_fetch_fails() {
        // One failed control region out of 50: zero counts would be dropped
        // by `select_usable_control_regions` (49/50 usable) and the sample
        // would be called. It must be an error on both the serial and the
        // chunked parallel path instead.
        let Some(bam) = benchmark_bam() else { return };
        let bam = bam.to_str().unwrap();
        let mut regions =
            crate::depth_calling::utilities::parse_region_file(crate::data::REGION_BED, "chr38");
        let control = &mut regions.get_mut("norm").expect("norm regions")[7].0;
        std::mem::swap(&mut control.1, &mut control.2);
        for n_cores in [1usize, 4] {
            let r = get_normed_depth(bam, &regions, n_cores, None, None);
            assert!(r.is_err(), "n_cores={n_cores}: a failed control-region fetch must be an error");
        }
    }

    // ---- Multi-mapper input guard (research P1 / H4; plan locked in
    // /home/parlar_ai/cyrius-research/exp-guard, commit 2b1a973) ----
    //
    // cyrius-rs reads copy number from reads that fit CYP2D6 and CYP2D7
    // equally well (MAPQ 0). An aligner that drops or reassigns them changes
    // the evidence itself: GIAB's novoalign BAMs have 0 MAPQ-0 reads in REP7,
    // and HG002 2x250 novoalign is called *4/*13 PASS against truth *2/*4.

    /// chr22 BAM with 50 bp reads around a test REP7 of [1000, 2000).
    fn write_mapq_fixture_bam(path: &std::path::Path) {
        let mut header = bam::header::Header::new();
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr22");
        sq.push_tag(b"LN", 50_818_468);
        header.push_record(&sq);
        // (name, 0-based start, MAPQ, flags). Coordinate-sorted.
        let reads: [(&str, i64, u8, u16); 8] = [
            ("starts_before_region", 980, 0, 0),
            ("mapq60_a", 1100, 60, 0),
            ("mapq60_b", 1200, 60, 0),
            ("mapq0", 1300, 0, 0),
            ("mapq0_duplicate", 1400, 0, 0x400),
            ("secondary", 1500, 0, 0x100),
            ("supplementary", 1600, 0, 0x800),
            ("unmapped_placed", 1700, 0, 0x4),
        ];
        {
            let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
            let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(50)]);
            for (name, pos, mapq, flags) in reads {
                let mut rec = bam::Record::new();
                rec.set(name.as_bytes(), Some(&cigar), &[b'A'; 50], &[30u8; 50]);
                rec.set_tid(0);
                rec.set_pos(pos);
                rec.set_mapq(mapq);
                rec.set_flags(flags); // Record::new() leaves flags = 4 (unmapped)
                writer.write(&rec).unwrap();
            }
        }
        bam::index::build(path, None::<&std::path::Path>, bam::index::Type::Bai, 1).unwrap();
    }

    fn mapq_fixture(tag: &str) -> (std::path::PathBuf, std::path::PathBuf) {
        let dir = std::env::temp_dir().join(format!("cyrius_mapq_{tag}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let bam_path = dir.join("mapq.bam");
        write_mapq_fixture_bam(&bam_path);
        (dir, bam_path)
    }

    #[test]
    fn get_mapq0_counts_counts_mapped_primary_reads_starting_in_the_region() {
        // Counted: mapq60_a, mapq60_b, mapq0, mapq0_duplicate (duplicates
        // stay in, as in get_read_count). Not counted: the read starting
        // before the region, secondary, supplementary, and the unmapped read.
        let (dir, bam_path) = mapq_fixture("counts");
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();
        let region: Region = ("chr22".to_string(), 1000, 2000, "REP7".to_string());
        let counts = get_mapq0_counts(&mut reader, &region);
        drop(reader);
        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(counts, Ok((4, 2)), "(mapped, MAPQ 0)");
    }

    #[test]
    fn get_mapq0_counts_is_an_error_when_the_fetch_fails() {
        // E1: a failed read is a no-call, here as everywhere else.
        let (dir, bam_path) = mapq_fixture("fetch_err");
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();
        let inverted: Region = ("chr22".to_string(), 2000, 1000, "REP7".to_string());
        let r = get_mapq0_counts(&mut reader, &inverted);
        drop(reader);
        let _ = std::fs::remove_dir_all(&dir);
        assert!(r.is_err(), "{r:?}");
    }

    #[test]
    fn mapq0_fraction_is_none_without_mapped_reads() {
        assert_eq!(mapq0_fraction((0, 0)), None);
        assert_eq!(mapq0_fraction((4, 2)), Some(0.5));
        assert_eq!(mapq0_fraction((133, 0)), Some(0.0));
    }

    #[test]
    fn multimapper_guard_flags_pass_and_fuzzy_match_when_the_header_shows_a_drop() {
        let flag = |f: Option<&str>, h: Option<bool>| {
            downgrade_filter_for_multimapper_depletion(f.map(str::to_string), h)
        };
        let depleted = Some("Multimapper_depleted".to_string());
        assert_eq!(flag(Some("PASS"), Some(true)), depleted);
        assert_eq!(flag(Some("Fuzzy_match"), Some(true)), depleted);
        assert_eq!(flag(Some("Unknown_sites_unmatched"), Some(true)), depleted);
        assert_eq!(flag(Some("Unknown_sites_resolved"), Some(true)), depleted);
        assert_eq!(flag(Some("PASS"), Some(false)), Some("PASS".to_string()));
        // Header not read: no flag.
        assert_eq!(flag(Some("PASS"), None), Some("PASS".to_string()));
        // Filters that already explain a non-clean call, and no-calls, stay.
        assert_eq!(
            flag(Some("More_than_one_possible_genotype"), Some(true)),
            Some("More_than_one_possible_genotype".to_string())
        );
        assert_eq!(flag(None, Some(true)), None);
    }

    #[test]
    fn header_check_flags_novoalign_without_a_repeat_option_and_samtools_mapq_filters() {
        let pg = |id: &str, pn: &str, cl: &str| format!("@HD\tVN:1.6\n@PG\tID:{id}\tPN:{pn}\tCL:{cl}\n");
        // GIAB novoalign 3.02.07: no -r, so multi-mapped reads are not reported.
        assert!(header_drops_multimappers(&pg(
            "novoalign",
            "novoalign",
            "novoalign -d ref.ndx -f r1.fq.gz r2.fq.gz -F STDFQ --Q2Off -t 400 -o SAM -c 10"
        )));
        // samtools merge renames clashing IDs; the PN still says novoalign.
        assert!(header_drops_multimappers(&pg("novoalign-4EFE4A6B", "novoalign", "novoalign -d ref.ndx")));
        // novoalign with -r reports multi-mapped reads.
        assert!(!header_drops_multimappers(&pg("novoalign", "novoalign", "novoalign -d ref.ndx -r Random -f r1.fq")));
        // A MAPQ filter of 1 or more drops MAPQ-0 reads; -q 0 drops nothing.
        assert!(header_drops_multimappers(&pg("samtools", "samtools", "samtools view -b -q 10 in.bam")));
        assert!(header_drops_multimappers(&pg("samtools", "samtools", "samtools view -b -q1 in.bam")));
        assert!(!header_drops_multimappers(&pg("samtools", "samtools", "samtools view -b -q 0 in.bam")));
        assert!(!header_drops_multimappers(&pg("samtools", "samtools", "samtools sort -o out.bam in.bam")));
        // bwa, region extraction and no @PG at all: no flag.
        assert!(!header_drops_multimappers(&pg("bwa", "bwa", "bwa mem -Y -K 100000000 -t 16 ref.fa r1.fq r2.fq")));
        assert!(!header_drops_multimappers(&pg("samtools", "samtools", "samtools view -b -h -T ref.fa in.cram chr22:1-2")));
        assert!(!header_drops_multimappers("@HD\tVN:1.6\n@SQ\tSN:chr22\tLN:50818468\n"));
    }

    #[test]
    fn header_check_flags_the_novoalign_benchmark_bams_and_not_the_nygc_ones() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../benchmark");
        if !dir.join("HG003_pgx_regions.bam").exists() {
            eprintln!("skipping: benchmark BAM not present");
            return;
        }
        let header = |s: &str| {
            let p = dir.join(format!("{s}_pgx_regions.bam"));
            let r = bam::Reader::from_path(&p).unwrap();
            String::from_utf8_lossy(r.header().as_bytes()).to_string()
        };
        for s in ["HG003", "HG004", "HG006", "HG007"] {
            assert!(header_drops_multimappers(&header(s)), "{s}");
        }
        for s in ["HG00096", "NA12878", "HG00589"] {
            assert!(!header_drops_multimappers(&header(s)), "{s}");
        }
    }

    #[test]
    fn get_normed_depth_reports_the_rep7_mapq0_fraction() {
        // Expected counts are exp-H4's features.tsv (n_rep7_mapped,
        // n_rep7_mapq0): HG00096 971 / 813 (NYGC bwa), HG003 133 / 0
        // (GIAB novoalign).
        let regions =
            crate::depth_calling::utilities::parse_region_file(crate::data::REGION_BED, "chr38");
        for (sample, mapped, mapq0) in [("HG00096", 971.0, 813.0), ("HG003", 133.0, 0.0)] {
            let bam = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join(format!("../benchmark/{sample}_pgx_regions.bam"));
            if !bam.exists() {
                eprintln!("skipping: benchmark BAM not present");
                return;
            }
            let nd = get_normed_depth(bam.to_str().unwrap(), &regions, 1, None, None).unwrap();
            assert_eq!(nd.rep7_mapq0_fraction, Some(mapq0 / mapped), "{sample}");
        }
    }
}

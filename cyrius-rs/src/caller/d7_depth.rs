use rust_htslib::bam::{self, Read as BamRead};

/// CYP2D7 sub-region definitions (GRCh38 coordinates).
/// CYP2D7 is on the reverse strand at chr22:42139000-42144575.
/// Gene direction (5'→3') runs from high coords to low coords.
///
/// Sub-regions:
///   D7 exon9:           chr22:42140061-42140855  (exon9 region)
///   D7 intron4-intron8: chr22:42140855-42142569  (gene body, 3' half)
///   D7 exon2-intron8:   chr22:42140600-42143600  (gene body, broad)
///   D7 5'UTR-intron1:   chr22:42143600-42144575  (promoter/5' end)
///
/// D7_EXON9_* and D7_IN4_IN8_* are derived from `data/CYP2D6_SNP_38.txt`
/// (review finding #14 — they used to be swapped, inverting the *36
/// gene-conversion signal). That file's D7 column anchors each named
/// CYP2D6 feature to a CYP2D7 coordinate via the diagnostic D6/D7 SNP
/// pairs. The exon9-annotated rows (0-indexed sites 7 and 8) place exon9
/// at D7 42140315 (16bp allele, to 42140330) and 42140362 (10bp allele, to
/// 42140371). The boundary on each side is drawn at the midpoint to the
/// nearest neighboring annotated feature. On the 3' (low-coord) side, the
/// nearest neighbor is the last downstream_exon9 row at 42139807 (1bp), so
/// midpoint(42139807, 42140315) = 42140061. On the 5' (high-coord) side,
/// the nearest neighbor is the exon7 row at 42141339 (1bp), so
/// midpoint(42140371, 42141339) = 42140855. D7_IN4_IN8 picks up at that
/// same 42140855 boundary (no gap, no overlap) and runs through the
/// exon7/intron6/exon6/intron5/intron4 rows — intron4's last row ends at
/// 42142406 (1bp) — to the midpoint with the next annotated feature,
/// exon3 at 42142732 (1bp): midpoint(42142406, 42142732) = 42142569. This
/// keeps D7_IN4_IN8 a real gene-body window (bounded by exon7 through
/// intron4), not merely "whatever is left after excluding exon9".
///
/// Additional composite regions for *13 detection:
///   D7 3' half:  chr22:42139000-42142455  (intron4 to 3' end)
///   D7 5' half:  chr22:42142465-42144575  (5'UTR to intron4)
const D7_EXON9_START: i64 = 42140061;
const D7_EXON9_END: i64 = 42140855;
const D7_IN4_IN8_START: i64 = 42140855;
const D7_IN4_IN8_END: i64 = 42142569;
const D7_EX2_IN8_START: i64 = 42140600;
const D7_EX2_IN8_END: i64 = 42143600;
const D7_5PR_IN1_START: i64 = 42143600;
const D7_5PR_IN1_END: i64 = 42144575;
// Composite regions for *13 detection
const D7_3PR_HALF_START: i64 = 42139000;
const D7_3PR_HALF_END: i64 = 42142455;
const D7_5PR_HALF_START: i64 = 42142465;
const D7_5PR_HALF_END: i64 = 42144575;

/// Depth measurement for a single region.
#[derive(Debug, Clone, Copy)]
pub struct RegionDepth {
    pub start: i64,
    pub end: i64,
    pub read_count: u64,
    pub avg_depth: f64,
}

/// CYP2D7 sub-region depth profile.
#[derive(Debug, Clone)]
pub struct D7DepthProfile {
    pub d7_exon9: RegionDepth,
    pub d7_in4_in8: RegionDepth,
    pub d7_ex2_in8: RegionDepth,
    pub d7_5pr_in1: RegionDepth,
    pub d7_3pr_half: RegionDepth,
    pub d7_5pr_half: RegionDepth,
}

/// Hybrid classification from D7 depth signals.
#[derive(Debug, Clone)]
pub struct D7HybridSignals {
    /// D7 exon9 / D7 intron4-intron8 ratio (elevated = *36 gene conversion)
    pub exon9_body_ratio: f64,
    /// D7 3' half / D7 5' half ratio (reduced = *13 gene conversion)
    pub three_prime_five_prime_ratio: f64,
    /// Whether D7 exon9 depth is elevated relative to body (suggests *36)
    pub exon9_elevated: bool,
    /// Whether D7 3' region is depleted relative to 5' (suggests *13)
    pub three_prime_depleted: bool,
}

/// Count reads aligning to a region with minimum mapping quality.
fn count_reads_in_region(
    reader: &mut bam::IndexedReader,
    tid: u32,
    start: i64,
    end: i64,
    min_mapq: u8,
) -> Option<u64> {
    if reader
        .fetch(bam::FetchDefinition::Region(tid as i32, start, end))
        .is_err()
    {
        return None;
    }

    let mut count: u64 = 0;
    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        if result.is_err() {
            return None;
        }
        if record.mapq() >= min_mapq
            && !record.is_secondary()
            && !record.is_supplementary()
            && record.pos() >= start
            && record.pos() < end
        {
            count += 1;
        }
    }
    Some(count)
}

fn measure_region(
    reader: &mut bam::IndexedReader,
    tid: u32,
    start: i64,
    end: i64,
    min_mapq: u8,
    read_length: f64,
) -> Option<RegionDepth> {
    let read_count = count_reads_in_region(reader, tid, start, end, min_mapq)?;
    let region_len = (end - start) as f64;
    let avg_depth = if region_len > 0.0 {
        read_count as f64 * read_length / region_len
    } else {
        0.0
    };
    Some(RegionDepth {
        start,
        end,
        read_count,
        avg_depth,
    })
}

/// Compute depth at all CYP2D7 sub-regions.
pub fn compute_d7_depth(
    bam_path: &str,
    chrom: &str,
    reference_fasta: Option<&str>,
) -> Option<D7DepthProfile> {
    let mut reader = crate::depth_calling::utilities::open_alignment_file_with_index(
        bam_path,
        reference_fasta,
        None,
    )
    .ok()?;

    let tid = reader.header().tid(chrom.as_bytes())?;

    // Estimate read length from D7 region
    reader
        .fetch(bam::FetchDefinition::Region(
            tid as i32,
            D7_IN4_IN8_START,
            D7_5PR_IN1_END,
        ))
        .ok()?;
    let mut lengths = Vec::new();
    let mut record = bam::Record::new();
    let mut n = 0;
    while let Some(result) = reader.read(&mut record) {
        if result.is_err() {
            return None;
        }
        lengths.push(record.seq_len() as f64);
        n += 1;
        if n > 500 {
            break;
        }
    }
    let read_length = if lengths.is_empty() {
        150.0
    } else {
        lengths.iter().sum::<f64>() / lengths.len() as f64
    };

    let min_mapq = 0; // Match bin_count behavior

    let d7_exon9 = measure_region(&mut reader, tid, D7_EXON9_START, D7_EXON9_END, min_mapq, read_length)?;
    let d7_in4_in8 = measure_region(&mut reader, tid, D7_IN4_IN8_START, D7_IN4_IN8_END, min_mapq, read_length)?;
    let d7_ex2_in8 = measure_region(&mut reader, tid, D7_EX2_IN8_START, D7_EX2_IN8_END, min_mapq, read_length)?;
    let d7_5pr_in1 = measure_region(&mut reader, tid, D7_5PR_IN1_START, D7_5PR_IN1_END, min_mapq, read_length)?;
    let d7_3pr_half = measure_region(&mut reader, tid, D7_3PR_HALF_START, D7_3PR_HALF_END, min_mapq, read_length)?;
    let d7_5pr_half = measure_region(&mut reader, tid, D7_5PR_HALF_START, D7_5PR_HALF_END, min_mapq, read_length)?;

    Some(D7DepthProfile {
        d7_exon9,
        d7_in4_in8,
        d7_ex2_in8,
        d7_5pr_in1,
        d7_3pr_half,
        d7_5pr_half,
    })
}

/// Classify hybrid signals from D7 depth profile.
///
/// *36 detection: D7 exon9 is elevated relative to D7 gene body.
///   When *36 gene conversion happens, D6 exon9 → D7 sequence, so reads map to D7 exon9.
///   StellarPGx threshold: 2.5 < (2 × D7_exon9 / D7_in4_in8) < 3.5
///
/// *13 detection: D7 3' half (intron4 → 3'end) is depleted relative to 5' half.
///   When *13 conversion happens, D7 loses its exon2-intron4 region to D6.
///   StellarPGx threshold: 0.45 < (D7_3pr / D7_5pr) < 0.75
pub fn classify_d7_signals(profile: &D7DepthProfile) -> D7HybridSignals {
    // D7 exon9 / body ratio
    let exon9_body_ratio = if profile.d7_in4_in8.avg_depth > 0.0 {
        profile.d7_exon9.avg_depth / profile.d7_in4_in8.avg_depth
    } else {
        1.0
    };

    // 3' / 5' ratio for *13 detection
    let three_prime_five_prime_ratio = if profile.d7_5pr_half.avg_depth > 0.0 {
        profile.d7_3pr_half.avg_depth / profile.d7_5pr_half.avg_depth
    } else {
        1.0
    };

    // *36 detection: D7 exon9 elevated
    // StellarPGx uses: 2.5 < (2 * d7_exon9 / d7_in4_in8) < 3.5
    // Equivalent: 1.25 < ratio < 1.75
    let exon9_elevated = (1.25..1.75).contains(&exon9_body_ratio);

    // *13 detection: D7 3' half depleted.
    // StellarPGx cites 0.45 < ratio < 0.75; we deliberately widen to
    // 0.40..0.80 to tolerate depth noise in this (logging-only) signal.
    let three_prime_depleted = (0.40..0.80).contains(&three_prime_five_prime_ratio);

    D7HybridSignals {
        exon9_body_ratio,
        three_prime_five_prime_ratio,
        exon9_elevated,
        three_prime_depleted,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::depth_calling::utilities::get_snp_position;

    /// Review finding #14: D7_EXON9_* and D7_IN4_IN8_* were swapped, so the
    /// window named "exon9" excluded the actual exon9 diagnostic SNP sites
    /// (CYP2D6_SNP_38.txt data rows 7 and 8, 0-indexed) and the window named
    /// "intron4-intron8" contained them instead — inverting the *36
    /// gene-conversion signal.
    ///
    /// Resolve both sites' D7 positions from the real embedded SNP file at
    /// run time (never hardcode the position) so this test breaks if either
    /// the data file's contents or the window constants in this file move.
    #[test]
    fn exon9_window_contains_exon9_snp_sites_and_body_window_does_not() {
        let snp_db = get_snp_position(crate::data::SNP_FILE, "38", None);

        // dsnp2 keys are formatted "<D7 position>_<0-based row index>" by
        // get_snp_position(). Find the D7 position for a given row index
        // without ever naming the position literal ourselves.
        let d7_pos_for_site = |idx: usize| -> i64 {
            let target = idx.to_string();
            let key = snp_db
                .dsnp2
                .keys()
                .find(|k| k.rsplit('_').next() == Some(target.as_str()))
                .unwrap_or_else(|| panic!("site index {idx} not found in SNP_FILE dsnp2"));
            key.rsplit_once('_')
                .expect("dsnp2 key must contain an underscore")
                .0
                .parse()
                .expect("D7 position prefix must be numeric")
        };

        let site7_d7 = d7_pos_for_site(7);
        let site8_d7 = d7_pos_for_site(8);

        for pos in [site7_d7, site8_d7] {
            assert!(
                (D7_EXON9_START..D7_EXON9_END).contains(&pos),
                "exon9 SNP site D7 pos {pos} is not inside the exon9 window \
                 {D7_EXON9_START}..{D7_EXON9_END}"
            );
            assert!(
                !(D7_IN4_IN8_START..D7_IN4_IN8_END).contains(&pos),
                "exon9 SNP site D7 pos {pos} must not fall inside the \
                 intron4-intron8 window {D7_IN4_IN8_START}..{D7_IN4_IN8_END}"
            );
        }
    }
}

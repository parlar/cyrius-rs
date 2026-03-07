use rust_htslib::bam::{self, Read as BamRead};

/// CYP2D7 sub-region definitions (GRCh38 coordinates).
/// CYP2D7 is on the reverse strand at chr22:42139000-42144575.
/// Gene direction (5'→3') runs from high coords to low coords.
///
/// Sub-regions (from StellarPGx test3.bed):
///   D7 intron4-intron8: chr22:42139500-42140600  (gene body, 3' half)
///   D7 exon9:           chr22:42140600-42142500  (exon9 region)
///   D7 exon2-intron8:   chr22:42140600-42143600  (gene body, broad)
///   D7 5'UTR-intron1:   chr22:42143600-42144575  (promoter/5' end)
///
/// Additional composite regions for *13 detection:
///   D7 3' half:  chr22:42139000-42142455  (intron4 to 3' end)
///   D7 5' half:  chr22:42142465-42144575  (5'UTR to intron4)

const D7_EXON9_START: i64 = 42140600;
const D7_EXON9_END: i64 = 42142500;
const D7_IN4_IN8_START: i64 = 42139500;
const D7_IN4_IN8_END: i64 = 42140600;
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
) -> u64 {
    if reader
        .fetch(bam::FetchDefinition::Region(tid as i32, start, end))
        .is_err()
    {
        return 0;
    }

    let mut count: u64 = 0;
    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        if result.is_err() {
            continue;
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
    count
}

fn measure_region(
    reader: &mut bam::IndexedReader,
    tid: u32,
    start: i64,
    end: i64,
    min_mapq: u8,
    read_length: f64,
) -> RegionDepth {
    let read_count = count_reads_in_region(reader, tid, start, end, min_mapq);
    let region_len = (end - start) as f64;
    let avg_depth = if region_len > 0.0 {
        read_count as f64 * read_length / region_len
    } else {
        0.0
    };
    RegionDepth {
        start,
        end,
        read_count,
        avg_depth,
    }
}

/// Compute depth at all CYP2D7 sub-regions.
pub fn compute_d7_depth(
    bam_path: &str,
    chrom: &str,
    reference_fasta: Option<&str>,
) -> Option<D7DepthProfile> {
    let mut reader = if let Some(ref_path) = reference_fasta {
        let mut r = bam::IndexedReader::from_path(bam_path).ok()?;
        r.set_reference(ref_path).ok()?;
        r
    } else {
        bam::IndexedReader::from_path(bam_path).ok()?
    };

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
            continue;
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

    let d7_exon9 = measure_region(&mut reader, tid, D7_EXON9_START, D7_EXON9_END, min_mapq, read_length);
    let d7_in4_in8 = measure_region(&mut reader, tid, D7_IN4_IN8_START, D7_IN4_IN8_END, min_mapq, read_length);
    let d7_ex2_in8 = measure_region(&mut reader, tid, D7_EX2_IN8_START, D7_EX2_IN8_END, min_mapq, read_length);
    let d7_5pr_in1 = measure_region(&mut reader, tid, D7_5PR_IN1_START, D7_5PR_IN1_END, min_mapq, read_length);
    let d7_3pr_half = measure_region(&mut reader, tid, D7_3PR_HALF_START, D7_3PR_HALF_END, min_mapq, read_length);
    let d7_5pr_half = measure_region(&mut reader, tid, D7_5PR_HALF_START, D7_5PR_HALF_END, min_mapq, read_length);

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

    // *13 detection: D7 3' half depleted
    // StellarPGx uses: 0.45 < ratio < 0.75
    let three_prime_depleted = (0.40..0.80).contains(&three_prime_five_prime_ratio);

    D7HybridSignals {
        exon9_body_ratio,
        three_prime_five_prime_ratio,
        exon9_elevated,
        three_prime_depleted,
    }
}

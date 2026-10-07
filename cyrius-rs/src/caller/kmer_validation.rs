//! K-mer-based structural validation for CYP2D6 calls.
//!
//! Uses alignment-free k-mer counting at diagnostic D6/D7 positions to compute
//! a paralog ratio profile, then runs CBS changepoint detection to independently
//! classify the structural configuration. Compares against the depth-based cnvtag
//! to flag discordant calls.
//!
//! Based on the k-mer approach from pdx_caller.

use rust_htslib::bam::{self, Read as HtslibRead};
use serde::Deserialize;
use std::collections::HashMap;

use crate::caller::changepoint::{detect_changepoints, ChangePointResult};

// ---------------------------------------------------------------------------
// Diagnostic panel types (deserialized from embedded JSON)
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct DiagnosticPanel {
    k: usize,
    #[allow(dead_code)]
    n_diagnostic: usize,
    positions: Vec<DiagnosticPosition>,
}

#[derive(Deserialize)]
struct DiagnosticPosition {
    #[allow(dead_code)]
    aln_pos: usize,
    cyp2d6_kmer: Option<String>,
    cyp2d7_kmer: Option<String>,
    kmer_unique: bool,
}

// ---------------------------------------------------------------------------
// K-mer lookup table
// ---------------------------------------------------------------------------

struct KmerLookup {
    /// kmer bytes -> (position_index, is_cyp2d6)
    table: HashMap<Vec<u8>, (usize, bool)>,
    k: usize,
    n_positions: usize,
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect()
}

fn build_lookup(panel: &DiagnosticPanel) -> KmerLookup {
    let mut table = HashMap::new();
    let mut n_positions = 0;

    for (idx, pos) in panel.positions.iter().enumerate() {
        if !pos.kmer_unique {
            continue;
        }
        n_positions = idx + 1;

        if let Some(ref kmer) = pos.cyp2d6_kmer {
            let fwd = kmer.as_bytes().to_vec();
            let rc = reverse_complement(&fwd);
            table.insert(fwd, (idx, true));
            table.insert(rc, (idx, true));
        }
        if let Some(ref kmer) = pos.cyp2d7_kmer {
            let fwd = kmer.as_bytes().to_vec();
            let rc = reverse_complement(&fwd);
            table.insert(fwd, (idx, false));
            table.insert(rc, (idx, false));
        }
    }

    KmerLookup {
        table,
        k: panel.k,
        n_positions,
    }
}

// ---------------------------------------------------------------------------
// K-mer counting from BAM
// ---------------------------------------------------------------------------

/// Per-position k-mer counts.
pub struct KmerCounts {
    pub d6_counts: Vec<u32>,
    pub d7_counts: Vec<u32>,
    pub total_reads: u64,
    pub matched_reads: u64,
}

impl KmerCounts {
    /// Compute D6/(D6+D7) ratio at each position. None if total < min_depth.
    pub fn ratios(&self, min_depth: u32) -> Vec<Option<f64>> {
        self.d6_counts
            .iter()
            .zip(self.d7_counts.iter())
            .map(|(&d6, &d7)| {
                let total = d6 + d7;
                if total >= min_depth {
                    Some(d6 as f64 / total as f64)
                } else {
                    None
                }
            })
            .collect()
    }

    /// Overall D6/(D6+D7) ratio across all positions, or `None` when no
    /// diagnostic k-mers were counted at all.
    ///
    /// The guard is on the *counts*, not on `total_reads`: the ratio is
    /// undefined exactly when its denominator (D6+D7) is zero, whatever the
    /// cause -- a chromosome missing from the BAM header, *both* region
    /// fetches failing, records without SEQ, or reads that simply carry no
    /// diagnostic k-mer. `matched_reads == 0` is equivalent to this test (a
    /// matched read always increments one of the two counters). Returning
    /// `None` keeps "no data" in the type, so it cannot be mistaken for a
    /// measured 50/50 ratio.
    ///
    /// A *single* failed fetch is deliberately not covered: `validate_with_kmers`
    /// merges a D6-region count and a D7-region count, so if only one fetch
    /// fails the other region's counts survive, the denominator is non-zero
    /// and validation proceeds on half the data behind a `log::warn` alone.
    ///
    /// The *ratio* is equally unmeasurable whether nothing was scanned or many
    /// reads were scanned without a single diagnostic match -- but those two
    /// states are not equally informative, so a caller that gives up here must
    /// still report the real `total_reads`/`matched_reads` rather than zeroing
    /// them (see the zero-count path in `validate_with_kmers`).
    pub fn overall_ratio(&self) -> Option<f64> {
        let total_d6: u64 = self.d6_counts.iter().map(|&c| c as u64).sum();
        let total_d7: u64 = self.d7_counts.iter().map(|&c| c as u64).sum();
        let total = total_d6 + total_d7;
        if total == 0 {
            None
        } else {
            Some(total_d6 as f64 / total as f64)
        }
    }
}

/// Decode a BAM record's 4-bit encoded sequence into ASCII bytes.
fn decode_bam_seq(record: &bam::Record) -> Vec<u8> {
    let len = record.seq_len();
    let mut seq = Vec::with_capacity(len);
    let bam_seq = record.seq();
    for i in 0..len {
        seq.push(bam_seq[i]);
    }
    seq
}

/// Count diagnostic k-mers in all reads from a BAM region.
fn count_kmers_in_region(
    reader: &mut bam::IndexedReader,
    chrom: &str,
    start: i64,
    end: i64,
    lookup: &KmerLookup,
) -> KmerCounts {
    let mut d6_counts = vec![0u32; lookup.n_positions];
    let mut d7_counts = vec![0u32; lookup.n_positions];
    let mut total_reads = 0u64;
    let mut matched_reads = 0u64;

    let tid = {
        let header = reader.header().clone();
        let mut found = None;
        for i in 0..header.target_count() {
            let name = String::from_utf8_lossy(header.tid2name(i));
            if name == chrom {
                found = Some(i as i64);
                break;
            }
        }
        match found {
            Some(t) => t,
            None => {
                log::warn!("kmer_validation: chromosome '{}' not found in BAM header", chrom);
                return KmerCounts {
                    d6_counts,
                    d7_counts,
                    total_reads: 0,
                    matched_reads: 0,
                };
            }
        }
    };

    if reader
        .fetch(bam::FetchDefinition::Region(tid as i32, start, end))
        .is_err()
    {
        log::warn!("kmer_validation: failed to fetch region {}:{}-{}", chrom, start, end);
        return KmerCounts {
            d6_counts,
            d7_counts,
            total_reads: 0,
            matched_reads: 0,
        };
    }

    let k = lookup.k;

    for result in reader.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };

        if record.is_secondary() || record.is_supplementary() || record.is_duplicate() {
            continue;
        }

        total_reads += 1;
        let seq = decode_bam_seq(&record);

        if seq.len() < k {
            continue;
        }

        let mut read_matched = false;
        for window in seq.windows(k) {
            // Skip windows with N
            if window.iter().any(|&b| b == b'N' || b == b'n') {
                continue;
            }
            if let Some(&(pos_idx, is_d6)) = lookup.table.get(window) {
                if is_d6 {
                    d6_counts[pos_idx] += 1;
                } else {
                    d7_counts[pos_idx] += 1;
                }
                read_matched = true;
            }
        }
        if read_matched {
            matched_reads += 1;
        }
    }

    KmerCounts {
        d6_counts,
        d7_counts,
        total_reads,
        matched_reads,
    }
}

// ---------------------------------------------------------------------------
// Structural classification from k-mer ratios
// ---------------------------------------------------------------------------

/// Classification result from k-mer ratio analysis.
#[derive(Debug, Clone)]
pub struct KmerStructuralClass {
    /// Broad structural category: "normal", "deletion", "duplication", "hybrid"
    pub category: String,
    /// Estimated D6 copy number from overall ratio
    pub estimated_d6_cn: f64,
    /// Estimated D7 copy number
    pub estimated_d7_cn: f64,
    /// Whether a hybrid breakpoint was detected by CBS
    pub hybrid_detected: bool,
    /// Changepoint result details
    pub changepoint: ChangePointResult,
}

/// Classify structural configuration from k-mer ratio profile.
///
/// Uses the overall D6/(D6+D7) ratio to estimate D6 copy number, which is
/// the most reliable k-mer signal. CBS changepoint detection is run as
/// supplementary information but does NOT drive the classification, since
/// the 189-position panel has too much inherent ratio variation for reliable
/// changepoint detection on regional BAMs.
///
/// `overall` must be a *measured* ratio: it is taken as an argument rather
/// than read back off `counts` so that this function no longer invents a
/// category out of a default ratio when nothing was counted.
///
/// Nothing here *enforces* that, though -- the parameter is a bare `f64`
/// carrying a prose-only contract, so a caller that passes a fabricated ratio
/// would still get a full category and copy-number estimate back. The
/// obligation is on the caller: `validate_with_kmers` returns early when
/// `KmerCounts::overall_ratio()` is `None`, so this function is never reached
/// with a default ratio on the live path (see `KmerCounts::overall_ratio`).
fn classify_from_ratios(
    counts: &KmerCounts,
    total_cn: u32,
    overall: f64,
) -> KmerStructuralClass {
    // Estimated D6 and D7 copy numbers
    let estimated_d6_cn = total_cn as f64 * overall;
    let estimated_d7_cn = total_cn as f64 - estimated_d6_cn;

    // Filter positions for CBS: require both d6 AND d7 to have reads (min 3 each).
    let min_each = 3u32;
    let d6_filtered: Vec<usize> = counts.d6_counts.iter().zip(counts.d7_counts.iter())
        .map(|(&d6, &d7)| if d6 >= min_each && d7 >= min_each { d6 as usize } else { 0 })
        .collect();
    let d7_filtered: Vec<usize> = counts.d6_counts.iter().zip(counts.d7_counts.iter())
        .map(|(&d6, &d7)| if d6 >= min_each && d7 >= min_each { d7 as usize } else { 0 })
        .collect();

    // Run CBS as supplementary info (NOT used for classification)
    let cp_result = detect_changepoints(&d6_filtered, &d7_filtered, 10, 4.0, 8);

    // Classify purely by D6 copy number from overall ratio
    let d6_round = estimated_d6_cn.round() as u32;
    let category = match d6_round {
        0 => "deletion_homozygous".to_string(),
        1 => "deletion".to_string(),
        2 => "normal".to_string(),
        _ => "duplication".to_string(),
    };

    KmerStructuralClass {
        category,
        estimated_d6_cn,
        estimated_d7_cn,
        hybrid_detected: cp_result.is_hybrid,
        changepoint: cp_result,
    }
}

// ---------------------------------------------------------------------------
// Validation result
// ---------------------------------------------------------------------------

/// Result of k-mer validation.
#[derive(Debug, Clone)]
pub struct KmerValidationResult {
    /// Broad structural category from k-mer analysis
    pub kmer_category: String,
    /// Whether it agrees with the depth-based cnvtag
    pub agrees: bool,
    /// Overall D6/(D6+D7) k-mer ratio
    pub overall_ratio: f64,
    /// Estimated D6 copy number from k-mers
    pub estimated_d6_cn: f64,
    /// Estimated D7 copy number from k-mers
    pub estimated_d7_cn: f64,
    /// Whether a hybrid breakpoint was detected
    pub hybrid_detected: bool,
    /// Number of reads scanned
    pub total_reads: u64,
    /// Number of reads with diagnostic k-mer matches
    pub matched_reads: u64,
    /// Per-position ratio profile (None = insufficient depth)
    pub position_ratios: Vec<Option<f64>>,
    /// Per-position D6 counts
    pub d6_counts: Vec<u32>,
    /// Per-position D7 counts
    pub d7_counts: Vec<u32>,
    /// Changepoint segments summary
    pub segments_summary: String,
}

impl KmerValidationResult {
    /// A no-result placeholder used when validation can't run (e.g. the panel
    /// fails to parse, or the BAM fails to open). Reports no agreement and no
    /// detection, with zero read counts -- which is only truthful when nothing
    /// was scanned. A caller that already scanned reads before giving up must
    /// put the real counts back over this constructor (see the zero-k-mer path
    /// in `validate_with_kmers`).
    fn unavailable(reason: &str) -> Self {
        KmerValidationResult {
            kmer_category: format!("unavailable: {}", reason),
            agrees: false,
            overall_ratio: 0.0,
            estimated_d6_cn: 0.0,
            estimated_d7_cn: 0.0,
            hybrid_detected: false,
            total_reads: 0,
            matched_reads: 0,
            position_ratios: Vec::new(),
            d6_counts: Vec::new(),
            d7_counts: Vec::new(),
            segments_summary: String::new(),
        }
    }
}

/// Compute expected k-mer D6 copy number for a given cnvtag and total_cn.
///
/// K-mers count D6-like SEQUENCE, not gene copies. In normal configurations,
/// D7 copy number is always 2 (diploid), so D6_expected = total_cn - 2.
///
/// For star68 hybrids (D6-5'/D7-3'), the hybrid contributes partial D6 k-mers
/// (roughly half D6, half D7), so we count each *68 as ~0.5 D6 equivalent less.
fn expected_d6_cn(tag: &str, total_cn: u32) -> f64 {
    let base_d6 = total_cn as f64 - 2.0;
    match tag {
        // star68 hybrids: the *68 allele has ~50% D6 and ~50% D7 sequence,
        // contributing ~0.5 D6-equivalent less than a full D6 copy
        "star68" => base_d6 - 0.5,
        "star68_star68" | "star68_star68_star68" | "star68_star68_star68_star68" => {
            let n_star68 = tag.matches("star68").count() as f64;
            base_d6 - n_star68 * 0.5
        }
        "star5_star68" | "star5_star5_star68" => base_d6 - 0.5,
        "star13_star68" | "exon9hyb_star68" | "dup_star68" => base_d6 - 0.5,
        // star13 hybrids: *13 replaces a D7 copy (D7→D6 conversion at 5' end),
        // so D7=1 instead of 2. The *13 contributes ~50% D6 k-mers.
        // Net effect: D6_eff = (total_cn - 1) - 0.5 = base_d6 + 0.5
        "star13" | "star13intron1" => base_d6 + 0.5,
        "dup_star13" | "dup_star13intron1" => base_d6 + 0.5,
        "star13_star13" | "star13intron1_star13intron1" => base_d6 + 1.0,
        // For all other cnvtags: D7 is always 2, so D6 = total_cn - 2
        _ => base_d6,
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Run k-mer validation against a BAM file.
///
/// Counts diagnostic k-mers, builds a paralog ratio profile, runs CBS
/// changepoint detection, classifies the structural type, and compares
/// against the depth-based cnvtag.
pub fn validate_with_kmers(
    bam_path: &str,
    chrom: &str,
    total_cn: u32,
    // Currently unused; retained for API symmetry with the depth-based caller.
    _spacer_cn: Option<u32>,
    depth_cnvtag: &str,
    reference_fasta: Option<&str>,
    index_name: Option<&str>,
) -> KmerValidationResult {
    // Parse the embedded diagnostic panel
    let panel: DiagnosticPanel = match serde_json::from_str(crate::data::DIAGNOSTIC_PANEL) {
        Ok(p) => p,
        Err(e) => return KmerValidationResult::unavailable(&format!("panel parse error: {}", e)),
    };
    let lookup = build_lookup(&panel);

    // Open BAM and count k-mers across D6+D7 region
    // D6: chr22:42123000-42132000, D7: chr22:42135000-42146000
    let mut reader = match crate::depth_calling::utilities::open_alignment_file_with_index(
        bam_path,
        reference_fasta,
        index_name,
    ) {
        Ok(r) => r,
        Err(_) => return KmerValidationResult::unavailable("failed to open BAM"),
    };

    // Count in the D6 region
    let counts_d6 = count_kmers_in_region(
        &mut reader,
        chrom,
        42123000,
        42132000,
        &lookup,
    );

    // Count in the D7 region (reads mapping here may carry D7 k-mers)
    let counts_d7 = count_kmers_in_region(
        &mut reader,
        chrom,
        42135000,
        42146000,
        &lookup,
    );

    // Merge counts from both regions
    let mut merged = KmerCounts {
        d6_counts: vec![0u32; lookup.n_positions],
        d7_counts: vec![0u32; lookup.n_positions],
        total_reads: counts_d6.total_reads + counts_d7.total_reads,
        matched_reads: counts_d6.matched_reads + counts_d7.matched_reads,
    };
    for i in 0..lookup.n_positions {
        merged.d6_counts[i] = counts_d6.d6_counts[i] + counts_d7.d6_counts[i];
        merged.d7_counts[i] = counts_d6.d7_counts[i] + counts_d7.d7_counts[i];
    }

    // No diagnostic k-mers at all means there is nothing to validate against:
    // report it as unavailable rather than letting a default ratio fabricate
    // an estimate that happens to "agree" with the depth call.
    //
    // The scanned/matched read counts are real even here, so they are carried
    // through over `unavailable()`'s zeroes: "80,000 reads scanned, none
    // diagnostic" (wrong reference build, SEQ-less records, a panel/BAM
    // mismatch) must stay distinguishable from "the BAM never opened" without
    // parsing the category prose, since `main` emits both fields verbatim.
    let overall_ratio = match merged.overall_ratio() {
        Some(r) => r,
        None => {
            let reason = format!(
                "no diagnostic k-mers counted ({} reads scanned, {} matched)",
                merged.total_reads, merged.matched_reads
            );
            return KmerValidationResult {
                total_reads: merged.total_reads,
                matched_reads: merged.matched_reads,
                ..KmerValidationResult::unavailable(&reason)
            };
        }
    };

    let position_ratios = merged.ratios(5);

    // Classify structural type using D6 copy number
    let classification = classify_from_ratios(&merged, total_cn, overall_ratio);

    // Compare expected D6 CN (from depth cnvtag) vs observed D6 CN (from k-mers)
    let expected = expected_d6_cn(depth_cnvtag, total_cn);
    let agrees = (classification.estimated_d6_cn - expected).abs() < 0.7;

    // Format segments summary
    let segments_summary = classification
        .changepoint
        .segments
        .iter()
        .enumerate()
        .map(|(i, seg)| {
            format!(
                "seg{}[{}..{}]: ratio={:.3} n={}",
                i, seg.start_idx, seg.end_idx, seg.mean_ratio, seg.n_valid,
            )
        })
        .collect::<Vec<_>>()
        .join(", ");

    KmerValidationResult {
        kmer_category: classification.category,
        agrees,
        overall_ratio,
        estimated_d6_cn: classification.estimated_d6_cn,
        estimated_d7_cn: classification.estimated_d7_cn,
        hybrid_detected: classification.hybrid_detected,
        total_reads: merged.total_reads,
        matched_reads: merged.matched_reads,
        position_ratios,
        d6_counts: merged.d6_counts,
        d7_counts: merged.d7_counts,
        segments_summary,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Write a BAM holding nothing but an `@SQ` header line for `chr22`, plus
    /// its BAI index, and return its path. Any region fetched from it yields
    /// zero reads, which leaves the k-mer counters in exactly the state a
    /// missing contig, a failed region fetch or reads without SEQ produce.
    fn write_empty_indexed_bam(dir: &std::path::Path) -> std::path::PathBuf {
        let mut header = bam::header::Header::new();
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr22");
        sq.push_tag(b"LN", 50_818_468);
        header.push_record(&sq);
        let path = dir.join("empty.bam");
        {
            let _writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
        }
        bam::index::build(
            &path,
            None::<&std::path::PathBuf>,
            bam::index::Type::Bai,
            1,
        )
        .unwrap();
        path
    }

    /// How many unmatchable reads `write_unmatchable_indexed_bam` writes.
    const N_UNMATCHABLE_READS: u64 = 12;

    /// Write an indexed BAM whose reads all land inside the D6 region that
    /// `validate_with_kmers` fetches but whose sequence (poly-A) carries no
    /// diagnostic k-mer. Those reads are therefore *scanned* but never
    /// *matched* -- the state a wrong reference build, a panel/BAM mismatch or
    /// SEQ-less records leave behind, and the one where zeroing `total_reads`
    /// would destroy the only evidence that anything was read at all.
    fn write_unmatchable_indexed_bam(dir: &std::path::Path) -> std::path::PathBuf {
        let mut header = bam::header::Header::new();
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr22");
        sq.push_tag(b"LN", 50_818_468);
        header.push_record(&sq);
        let path = dir.join("unmatchable.bam");
        {
            let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
            let seq = vec![b'A'; 150];
            let qual = vec![30u8; 150];
            let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(150)]);
            for i in 0..N_UNMATCHABLE_READS {
                let mut rec = bam::Record::new();
                rec.set(format!("read{}", i).as_bytes(), Some(&cigar), &seq, &qual);
                rec.set_tid(0);
                // Inside the D6 window validate_with_kmers fetches
                // (chr22:42123000-42132000), ascending so the BAM is sorted.
                rec.set_pos(42_124_000 + i as i64);
                rec.set_mapq(60);
                writer.write(&rec).unwrap();
            }
        }
        bam::index::build(
            &path,
            None::<&std::path::PathBuf>,
            bam::index::Type::Bai,
            1,
        )
        .unwrap();
        path
    }

    /// A scratch directory unique to this process *and* the given tag, so the
    /// tests in this module never collide on it.
    fn scratch_dir(tag: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir()
            .join(format!("cyrius_kmer_{}_{}", tag, std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn counts(d6: Vec<u32>, d7: Vec<u32>) -> KmerCounts {
        let matched = d6.iter().chain(d7.iter()).map(|&c| c as u64).sum();
        KmerCounts {
            d6_counts: d6,
            d7_counts: d7,
            total_reads: matched,
            matched_reads: matched,
        }
    }

    #[test]
    fn overall_ratio_is_none_when_nothing_was_counted() {
        // All-zero counters: the ratio's denominator is zero, so there is no
        // ratio -- not a 50/50 one. Holds whether or not reads were scanned.
        assert_eq!(counts(vec![0; 4], vec![0; 4]).overall_ratio(), None);
        let mut scanned_but_unmatched = counts(vec![0; 4], vec![0; 4]);
        scanned_but_unmatched.total_reads = 1000;
        assert_eq!(scanned_but_unmatched.overall_ratio(), None);
        assert_eq!(counts(Vec::new(), Vec::new()).overall_ratio(), None);
    }

    #[test]
    fn overall_ratio_is_the_measured_fraction_when_counts_exist() {
        assert_eq!(counts(vec![3, 1], vec![3, 3]).overall_ratio(), Some(0.4));
        // A single count on one side is still a measurement, not "no data".
        assert_eq!(counts(vec![0, 1], vec![0, 0]).overall_ratio(), Some(1.0));
        assert_eq!(counts(vec![0, 0], vec![0, 1]).overall_ratio(), Some(0.0));
    }

    #[test]
    fn zero_kmer_counts_report_unavailable_not_agreement() {
        // Finding #50: with no diagnostic k-mers counted, overall_ratio()
        // used to hand back a fabricated 0.5, so for a cn2 sample
        // (total_cn = 4) the estimate came out at exactly 4 * 0.5 = 2.0,
        // which is precisely expected_d6_cn("cn2", 4) -- a validation that
        // "agreed" on no data at all, and called the category "normal".
        let dir = scratch_dir("zero_test");
        let bam_path = write_empty_indexed_bam(&dir);

        let result = validate_with_kmers(
            bam_path.to_str().unwrap(),
            "chr22",
            4,
            None,
            "cn2",
            None,
            None,
        );

        let _ = std::fs::remove_dir_all(&dir);

        assert_eq!(
            result.total_reads, 0,
            "an empty BAM has no reads to scan"
        );
        assert_eq!(
            result.matched_reads, 0,
            "an empty BAM must yield no k-mer-matched reads"
        );
        assert!(
            !result.agrees,
            "with zero diagnostic k-mers counted there is nothing to agree with, \
             but agrees={} category={:?} overall_ratio={} estimated_d6_cn={}",
            result.agrees, result.kmer_category, result.overall_ratio, result.estimated_d6_cn,
        );
        assert!(
            result.kmer_category.starts_with("unavailable"),
            "zero counts must report an 'unavailable' category, got {:?}",
            result.kmer_category
        );
        // Pin the *specific* arm. Every other assertion in this test is also
        // satisfied by the "failed to open BAM" arm, so without this the test
        // is one htslib behaviour change away from passing for the wrong
        // reason.
        assert!(
            result.kmer_category.contains("no diagnostic k-mers counted"),
            "must be the zero-k-mer arm, not some other 'unavailable' reason, got {:?}",
            result.kmer_category
        );
        assert_eq!(
            result.overall_ratio, 0.0,
            "zero counts must not report a measured ratio"
        );
        assert_eq!(
            result.estimated_d6_cn, 0.0,
            "zero counts must not report an estimated D6 copy number"
        );
    }

    #[test]
    fn scanned_but_unmatched_reads_are_still_reported_as_scanned() {
        // The zero-k-mer guard must not destroy the read counts on its way
        // out. "N reads scanned, none diagnostic" (wrong reference build,
        // panel/BAM mismatch, SEQ-less records) and "the BAM never opened"
        // are different failures, and main.rs emits total_reads/matched_reads
        // into the JSON verbatim -- so zeroing them here would leave a machine
        // consumer unable to tell the two apart without parsing the category
        // prose.
        let dir = scratch_dir("unmatched_test");
        let bam_path = write_unmatchable_indexed_bam(&dir);

        let result = validate_with_kmers(
            bam_path.to_str().unwrap(),
            "chr22",
            4,
            None,
            "cn2",
            None,
            None,
        );

        let _ = std::fs::remove_dir_all(&dir);

        assert_eq!(
            result.matched_reads, 0,
            "poly-A reads must not match any diagnostic k-mer, got {} matches",
            result.matched_reads
        );
        assert_eq!(
            result.total_reads, N_UNMATCHABLE_READS,
            "all {} scanned reads must be reported as scanned, got total_reads={}",
            N_UNMATCHABLE_READS, result.total_reads
        );
        // Still no fabricated agreement, and still the zero-k-mer arm.
        assert!(!result.agrees, "nothing was measured, so nothing can agree");
        assert!(
            result.kmer_category.contains("no diagnostic k-mers counted"),
            "expected the zero-k-mer arm, got {:?}",
            result.kmer_category
        );
    }

    #[test]
    fn classify_from_ratios_estimates_copy_number_from_the_ratio_it_is_given() {
        // The non-zero path: real counts, a real ratio, a real category.
        // 30 D6 hits against 10 D7 hits is a measured ratio of 0.75.
        let c = counts(vec![10, 10, 10, 0], vec![5, 5, 0, 0]);
        assert_eq!(c.overall_ratio(), Some(0.75));

        // (a) the live arrangement: the caller passes the measured ratio.
        let cls = classify_from_ratios(&c, 4, 0.75);
        assert_eq!(cls.estimated_d6_cn, 3.0);
        assert_eq!(cls.estimated_d7_cn, 1.0);
        assert_eq!(cls.category, "duplication");

        // (b) and (c): same counts, a different ratio. The estimate must track
        // the argument -- a body that re-derived the ratio off `counts`, or
        // that hardcoded `total_cn * 0.5`, would return (a)'s answer or 2.0
        // here instead.
        let cls = classify_from_ratios(&c, 4, 0.25);
        assert_eq!(cls.estimated_d6_cn, 1.0);
        assert_eq!(cls.estimated_d7_cn, 3.0);
        assert_eq!(cls.category, "deletion");

        let cls = classify_from_ratios(&c, 2, 0.0);
        assert_eq!(cls.estimated_d6_cn, 0.0);
        assert_eq!(cls.estimated_d7_cn, 2.0);
        assert_eq!(cls.category, "deletion_homozygous");

        // (d) the ratio that used to be fabricated from no data is a perfectly
        // ordinary measurement when it really is measured.
        let balanced = counts(vec![8, 8], vec![8, 8]);
        assert_eq!(balanced.overall_ratio(), Some(0.5));
        let cls = classify_from_ratios(&balanced, 4, 0.5);
        assert_eq!(cls.estimated_d6_cn, 2.0);
        assert_eq!(cls.estimated_d7_cn, 2.0);
        assert_eq!(cls.category, "normal");
    }
}

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

    /// Overall D6/(D6+D7) ratio across all positions.
    pub fn overall_ratio(&self) -> f64 {
        let total_d6: u64 = self.d6_counts.iter().map(|&c| c as u64).sum();
        let total_d7: u64 = self.d7_counts.iter().map(|&c| c as u64).sum();
        let total = total_d6 + total_d7;
        if total == 0 {
            0.5
        } else {
            total_d6 as f64 / total as f64
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

    reader
        .fetch(bam::FetchDefinition::Region(tid as i32, start, end))
        .unwrap();

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
fn classify_from_ratios(
    counts: &KmerCounts,
    total_cn: u32,
) -> KmerStructuralClass {
    let overall = counts.overall_ratio();

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
    spacer_cn: Option<u32>,
    depth_cnvtag: &str,
    reference_fasta: Option<&str>,
    index_name: Option<&str>,
) -> KmerValidationResult {
    // Parse the embedded diagnostic panel
    let panel: DiagnosticPanel =
        serde_json::from_str(crate::data::DIAGNOSTIC_PANEL).expect("Failed to parse diagnostic panel JSON");
    let lookup = build_lookup(&panel);

    // Open BAM and count k-mers across D6+D7 region
    // D6: chr22:42123000-42132000, D7: chr22:42135000-42146000
    let mut reader = crate::depth_calling::utilities::open_alignment_file_with_index(
        bam_path,
        reference_fasta,
        index_name,
    )
    .expect("Failed to open BAM for kmer_validation");

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

    let position_ratios = merged.ratios(5);
    let overall_ratio = merged.overall_ratio();

    // Classify structural type using D6 copy number
    let classification = classify_from_ratios(&merged, total_cn);

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

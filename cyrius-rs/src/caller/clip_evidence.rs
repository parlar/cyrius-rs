//! Soft-clip cluster detection for structural breakpoint evidence.
//!
//! Scans BAM reads in the CYP2D6 region for positions where multiple reads
//! have large soft clips on the same side. Clusters indicate structural
//! breakpoints (tandem junctions, deletion boundaries, hybrid transitions).
//!
//! Normal diploid samples show NO clip clusters in the CYP2D6 region.
//! Clip clusters are purely structural signals.

use rust_htslib::bam::{self, Read as BamRead};
use std::collections::HashMap;

/// A cluster of soft-clipped reads at a genomic position.
#[derive(Debug, Clone)]
pub struct ClipCluster {
    /// Genomic position (10bp bin start).
    pub position: i64,
    /// Number of reads with left (leading) soft clips >= min_clip_len.
    pub left_count: u32,
    /// Number of reads with right (trailing) soft clips >= min_clip_len.
    pub right_count: u32,
}

impl ClipCluster {
    /// Total clip count (left + right).
    pub fn total(&self) -> u32 {
        self.left_count + self.right_count
    }

    /// Dominant clip side: true if mostly left clips.
    pub fn is_left_dominant(&self) -> bool {
        self.left_count > self.right_count
    }
}

/// Summary of clip evidence across the CYP2D6 region.
#[derive(Debug, Clone)]
pub struct ClipEvidence {
    /// All clip clusters with >= min_reads.
    pub clusters: Vec<ClipCluster>,
    /// Total number of clipped reads (>= min_clip_len) in the region.
    pub total_clipped_reads: u32,
}

impl ClipEvidence {
    /// Check if there's a breakpoint cluster near a specific position (within tolerance).
    pub fn has_cluster_near(&self, pos: i64, tolerance: i64) -> bool {
        self.clusters
            .iter()
            .any(|c| (c.position - pos).abs() <= tolerance)
    }

    /// Get the strongest cluster (most reads) in a genomic range.
    pub fn strongest_cluster_in_range(&self, start: i64, end: i64) -> Option<&ClipCluster> {
        self.clusters
            .iter()
            .filter(|c| c.position >= start && c.position <= end)
            .max_by_key(|c| c.total())
    }

    /// Count clusters in a genomic range.
    pub fn clusters_in_range(&self, start: i64, end: i64) -> Vec<&ClipCluster> {
        self.clusters
            .iter()
            .filter(|c| c.position >= start && c.position <= end)
            .collect()
    }
}

/// Minimum soft-clip length to consider (bp).
const MIN_CLIP_LEN: u32 = 20;
/// Bin size for clustering nearby clip positions.
const BIN_SIZE: i64 = 10;
/// Minimum reads in a bin to form a cluster.
const MIN_CLUSTER_READS: u32 = 3;

/// Scan a BAM file for soft-clip clusters in the CYP2D6 region.
///
/// Examines reads in `chrom:region_start-region_end` and groups soft clips
/// by position (10bp bins). Returns clusters with >= 3 supporting reads.
pub fn detect_clip_clusters(
    bam_path: &str,
    chrom: &str,
    region_start: i64,
    region_end: i64,
    reference: Option<&str>,
) -> Option<ClipEvidence> {
    let mut reader = if let Some(ref_path) = reference {
        let mut r = bam::IndexedReader::from_path(bam_path).ok()?;
        r.set_reference(ref_path).ok()?;
        r
    } else {
        bam::IndexedReader::from_path(bam_path).ok()?
    };

    let header = reader.header().clone();
    let tid = header
        .tid(chrom.as_bytes())
        .or_else(|| {
            // Try with/without "chr" prefix
            let alt = if chrom.starts_with("chr") {
                &chrom[3..]
            } else {
                &format!("chr{}", chrom)
            };
            header.tid(alt.as_bytes())
        })?;

    reader
        .fetch(bam::FetchDefinition::Region(
            tid as i32,
            region_start,
            region_end,
        ))
        .ok()?;

    let mut left_clips: HashMap<i64, u32> = HashMap::new();
    let mut right_clips: HashMap<i64, u32> = HashMap::new();
    let mut total_clipped = 0u32;

    for result in reader.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };

        if record.is_unmapped()
            || record.is_secondary()
            || record.is_supplementary()
            || record.is_duplicate()
        {
            continue;
        }

        let pos = record.pos(); // 0-based
        let cigar = record.cigar();
        let ops: &[rust_htslib::bam::record::Cigar] = cigar.as_ref();
        if ops.is_empty() {
            continue;
        }

        // Check leading soft clip.
        if let rust_htslib::bam::record::Cigar::SoftClip(len) = ops[0] {
            if len as u32 >= MIN_CLIP_LEN {
                let bin = (pos / BIN_SIZE) * BIN_SIZE;
                *left_clips.entry(bin).or_insert(0) += 1;
                total_clipped += 1;
            }
        }

        // Check trailing soft clip.
        if ops.len() > 1 {
            if let rust_htslib::bam::record::Cigar::SoftClip(len) = ops[ops.len() - 1] {
                if len as u32 >= MIN_CLIP_LEN {
                    // Compute reference end position.
                    let ref_end = record.cigar().end_pos();
                    let bin = (ref_end / BIN_SIZE) * BIN_SIZE;
                    *right_clips.entry(bin).or_insert(0) += 1;
                    total_clipped += 1;
                }
            }
        }
    }

    // Build clusters from bins with >= MIN_CLUSTER_READS.
    let all_bins: std::collections::BTreeSet<i64> = left_clips
        .keys()
        .chain(right_clips.keys())
        .cloned()
        .collect();

    let mut clusters = Vec::new();
    for bin in all_bins {
        let l = left_clips.get(&bin).copied().unwrap_or(0);
        let r = right_clips.get(&bin).copied().unwrap_or(0);
        if l >= MIN_CLUSTER_READS || r >= MIN_CLUSTER_READS {
            clusters.push(ClipCluster {
                position: bin,
                left_count: l,
                right_count: r,
            });
        }
    }

    Some(ClipEvidence {
        clusters,
        total_clipped_reads: total_clipped,
    })
}

/// Known CYP2D6 region boundaries (GRCh38).
pub const D6_GENE_START: i64 = 42126499;
pub const D6_GENE_END: i64 = 42130881;
pub const D6_REP6_END: i64 = 42132032;
pub const SPACER_START: i64 = 42132000;
pub const SPACER_END: i64 = 42140000;
pub const D7_GENE_START: i64 = 42139676;
pub const D7_GENE_END: i64 = 42145745;

/// Interpret clip evidence for CYP2D6 structural calling.
///
/// Returns a summary of structural signals derived from clip clusters.
#[derive(Debug, Clone)]
pub struct StructuralClipSignals {
    /// Breakpoint cluster in the D6 gene body (tandem dup junction).
    pub d6_body_breakpoint: bool,
    /// Breakpoint cluster at REP6/spacer boundary.
    pub rep6_spacer_breakpoint: bool,
    /// Breakpoint cluster in the spacer region.
    pub spacer_breakpoint: bool,
    /// Breakpoint cluster in the D7 region.
    pub d7_breakpoint: bool,
    /// Total number of structural clip clusters.
    pub n_clusters: usize,
    /// Strongest cluster read count.
    pub max_cluster_reads: u32,
}

/// Analyze clip evidence and classify structural signals.
pub fn classify_clip_signals(evidence: &ClipEvidence) -> StructuralClipSignals {
    let d6_body = evidence
        .clusters_in_range(D6_GENE_START, D6_GENE_END)
        .iter()
        .any(|c| c.total() >= MIN_CLUSTER_READS);

    let rep6_spacer = evidence
        .clusters_in_range(D6_GENE_END, SPACER_START + 500)
        .iter()
        .any(|c| c.total() >= MIN_CLUSTER_READS);

    let spacer = evidence
        .clusters_in_range(SPACER_START + 500, SPACER_END - 500)
        .iter()
        .any(|c| c.total() >= MIN_CLUSTER_READS);

    let d7 = evidence
        .clusters_in_range(D7_GENE_START, D7_GENE_END)
        .iter()
        .any(|c| c.total() >= MIN_CLUSTER_READS);

    let max_reads = evidence
        .clusters
        .iter()
        .map(|c| c.total())
        .max()
        .unwrap_or(0);

    StructuralClipSignals {
        d6_body_breakpoint: d6_body,
        rep6_spacer_breakpoint: rep6_spacer,
        spacer_breakpoint: spacer,
        d7_breakpoint: d7,
        n_clusters: evidence.clusters.len(),
        max_cluster_reads: max_reads,
    }
}

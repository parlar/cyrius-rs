//! Haplotype-resolved structural typing for tandem vs duplication disambiguation.
//!
//! When the CNV caller detects an exon9 hybrid signature, this module uses
//! heterozygous star-allele-defining variants to partition reads into haplotypes,
//! then checks per-haplotype D6/D7 ratios at exon9 diagnostic positions.
//!
//! Key insight: structural disambiguation is a phasing problem. By separating
//! reads by haplotype, we convert ambiguous aggregate ratios into clean
//! per-haplotype signals.

use rust_htslib::bam::{self, Read as BamRead};
use std::collections::{HashMap, HashSet};

/// Result of haplotype-resolved structural analysis.
#[derive(Debug, Clone)]
pub struct HaplotypeResolveResult {
    /// Number of het sites used for phasing.
    pub n_het_sites: usize,
    /// Reads assigned to haplotype carrying ALT at phasing sites.
    pub n_reads_hap_alt: usize,
    /// Reads assigned to haplotype carrying REF at phasing sites.
    pub n_reads_hap_ref: usize,
    /// D6 fraction at exon9 for the ALT haplotype (None if insufficient reads).
    pub exon9_d6_frac_hap_alt: Option<f64>,
    /// D6 fraction at exon9 for the REF haplotype (None if insufficient reads).
    pub exon9_d6_frac_hap_ref: Option<f64>,
    /// Suggested structural correction, if any.
    /// e.g., "dup" if exon9hyb should be reclassified as duplication.
    pub suggested_cnv: Option<String>,
}

/// Positions for D6/D7 diagnostic comparison at exon9.
/// These are the multi-base comparison sites from CYP2D6_SNP_38.txt.
struct Exon9DiagSite {
    pos: i64,                // 1-based genomic position
    d6_alleles: Vec<Vec<u8>>, // D6-like sequences (may have multiple)
    d7_allele: Vec<u8>,       // D7-like sequence
}

fn exon9_diagnostic_sites() -> Vec<Exon9DiagSite> {
    vec![
        Exon9DiagSite {
            pos: 42126611,
            d6_alleles: vec![
                b"GTCACCAGGAAAGCAA".to_vec(),
                b"CTCACCAGGAAAGCAA".to_vec(),
            ],
            d7_allele: b"GTCACCAGAAAGCTGA".to_vec(),
        },
        Exon9DiagSite {
            pos: 42126658,
            d6_alleles: vec![b"AGTGGGCACC".to_vec()],
            d7_allele: b"GGCGGCCACG".to_vec(),
        },
    ]
}

/// Heterozygous variant site with reads carrying each allele.
struct HetSite {
    _pos: i64,
    reads_alt: HashSet<String>,
    reads_ref: HashSet<String>,
}

/// Find heterozygous target variant positions and collect read names per allele.
/// Uses the same variant positions as the main caller (var_list).
fn find_het_sites(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    var_list: &[String],
    var_alt: &[usize],
    var_ref: &[usize],
) -> Vec<HetSite> {
    let tid = match reader.header().tid(nchr.as_bytes()) {
        Some(t) => t,
        None => return Vec::new(),
    };

    let mut het_sites = Vec::new();

    for (i, var_name) in var_list.iter().enumerate() {
        if i >= var_alt.len() || i >= var_ref.len() {
            break;
        }

        let alt_count = var_alt[i];
        let ref_count = var_ref[i];
        let total = alt_count + ref_count;
        if total < 10 {
            continue;
        }

        // Check for heterozygosity: allele balance between 0.2 and 0.8
        let af = alt_count as f64 / total as f64;
        if af < 0.2 || af > 0.8 {
            continue;
        }

        // Parse position from variant name
        let pos = match parse_variant_pos(var_name) {
            Some(p) => p,
            None => continue,
        };

        // Parse the expected alt base
        let alt_base = match parse_variant_alt_base(var_name) {
            Some(b) => b,
            None => continue,
        };

        // Collect read names at this position
        let (reads_alt, reads_ref) = collect_reads_at_position(
            reader, tid, pos, alt_base,
        );

        if reads_alt.len() >= 3 && reads_ref.len() >= 3 {
            het_sites.push(HetSite {
                _pos: pos,
                reads_alt,
                reads_ref,
            });
        }
    }

    het_sites
}

/// Parse genomic position from variant name like "g.42127941G>A" -> 42127941.
fn parse_variant_pos(name: &str) -> Option<i64> {
    let s = name.strip_prefix("g.")?;
    let pos_str: String = s.chars().take_while(|c| c.is_ascii_digit()).collect();
    pos_str.parse().ok()
}

/// Parse alt base from variant name like "g.42127941G>A" -> 'A'.
fn parse_variant_alt_base(name: &str) -> Option<u8> {
    let s = name.strip_prefix("g.")?;
    // Find X>Y pattern
    if let Some(arrow_pos) = s.find('>') {
        let alt = s.as_bytes().get(arrow_pos + 1)?;
        if alt.is_ascii_alphabetic() {
            return Some(alt.to_ascii_uppercase());
        }
    }
    None
}

/// Collect read names carrying alt vs ref at a specific position.
fn collect_reads_at_position(
    reader: &mut bam::IndexedReader,
    tid: u32,
    pos_1based: i64,
    alt_base: u8,
) -> (HashSet<String>, HashSet<String>) {
    let mut reads_alt = HashSet::new();
    let mut reads_ref = HashSet::new();

    if reader
        .fetch(bam::FetchDefinition::Region(
            tid as i32,
            pos_1based - 1,
            pos_1based,
        ))
        .is_err()
    {
        return (reads_alt, reads_ref);
    }

    let mut pileups = reader.pileup();
    pileups.set_max_depth(i32::MAX as u32);

    for pileup_result in pileups {
        let pileup = match pileup_result {
            Ok(p) => p,
            Err(_) => continue,
        };
        if pileup.pos() as i64 + 1 != pos_1based {
            continue;
        }

        for alignment in pileup.alignments() {
            let record = alignment.record();
            if record.is_secondary()
                || record.is_supplementary()
                || record.is_duplicate()
                || record.is_unmapped()
            {
                continue;
            }
            if record.mapq() < 10 {
                continue;
            }
            if alignment.is_del() || alignment.is_refskip() {
                continue;
            }

            let qpos = match alignment.qpos() {
                Some(p) => p,
                None => continue,
            };

            let seq = record.seq().as_bytes();
            if qpos >= seq.len() {
                continue;
            }

            let read_name = String::from_utf8_lossy(record.qname()).to_string();
            let observed = seq[qpos].to_ascii_uppercase();

            if observed == alt_base {
                reads_alt.insert(read_name);
            } else {
                reads_ref.insert(read_name);
            }
        }
    }

    (reads_alt, reads_ref)
}

/// Assign reads to haplotypes using multiple het sites via majority voting.
/// Returns (hap_alt_reads, hap_ref_reads) — reads consistently on each haplotype.
fn assign_haplotypes(het_sites: &[HetSite]) -> (HashSet<String>, HashSet<String>) {
    if het_sites.is_empty() {
        return (HashSet::new(), HashSet::new());
    }

    // Count per-read: how many sites support ALT vs REF haplotype
    let mut read_alt_votes: HashMap<String, usize> = HashMap::new();
    let mut read_ref_votes: HashMap<String, usize> = HashMap::new();

    for site in het_sites {
        for name in &site.reads_alt {
            *read_alt_votes.entry(name.clone()).or_insert(0) += 1;
        }
        for name in &site.reads_ref {
            *read_ref_votes.entry(name.clone()).or_insert(0) += 1;
        }
    }

    // Assign reads to haplotype by majority vote
    let all_reads: HashSet<String> = read_alt_votes
        .keys()
        .chain(read_ref_votes.keys())
        .cloned()
        .collect();

    let mut hap_alt = HashSet::new();
    let mut hap_ref = HashSet::new();

    for name in all_reads {
        let alt_v = read_alt_votes.get(&name).copied().unwrap_or(0);
        let ref_v = read_ref_votes.get(&name).copied().unwrap_or(0);
        // Only assign if clearly on one side (no ties)
        if alt_v > ref_v {
            hap_alt.insert(name);
        } else if ref_v > alt_v {
            hap_ref.insert(name);
        }
        // Ties: skip (ambiguous)
    }

    (hap_alt, hap_ref)
}

/// Check D6/D7 signal at exon9 for a set of reads.
/// Returns D6 fraction (0.0 = all D7, 1.0 = all D6).
fn exon9_d6_fraction(
    reader: &mut bam::IndexedReader,
    tid: u32,
    hap_reads: &HashSet<String>,
    diag_sites: &[Exon9DiagSite],
) -> Option<f64> {
    if hap_reads.is_empty() {
        return None;
    }

    let mut d6_count = 0usize;
    let mut d7_count = 0usize;

    for site in diag_sites {
        if reader
            .fetch(bam::FetchDefinition::Region(
                tid as i32,
                site.pos - 1,
                site.pos,
            ))
            .is_err()
        {
            continue;
        }

        let mut pileups = reader.pileup();
        pileups.set_max_depth(i32::MAX as u32);

        for pileup_result in pileups {
            let pileup = match pileup_result {
                Ok(p) => p,
                Err(_) => continue,
            };
            if pileup.pos() as i64 + 1 != site.pos {
                continue;
            }

            for alignment in pileup.alignments() {
                let record = alignment.record();
                if record.is_secondary()
                    || record.is_supplementary()
                    || record.is_duplicate()
                    || record.is_unmapped()
                {
                    continue;
                }
                if alignment.is_del() || alignment.is_refskip() {
                    continue;
                }

                let read_name = String::from_utf8_lossy(record.qname()).to_string();
                if !hap_reads.contains(&read_name) {
                    continue; // Not on this haplotype
                }

                let qpos = match alignment.qpos() {
                    Some(p) => p,
                    None => continue,
                };

                let seq = record.seq().as_bytes();
                let allele_len = site.d7_allele.len(); // All alleles same length
                if qpos + allele_len > seq.len() {
                    continue;
                }

                let read_seq = &seq[qpos..qpos + allele_len];

                // Check D6 match
                let is_d6 = site
                    .d6_alleles
                    .iter()
                    .any(|d6| read_seq.eq_ignore_ascii_case(d6));
                let is_d7 = read_seq.eq_ignore_ascii_case(&site.d7_allele);

                if is_d6 {
                    d6_count += 1;
                } else if is_d7 {
                    d7_count += 1;
                }
            }
        }
    }

    let total = d6_count + d7_count;
    if total < 3 {
        return None; // Insufficient data
    }

    Some(d6_count as f64 / total as f64)
}

/// Check D6/D7 signal at exon9 for ALL reads (no haplotype filtering).
fn exon9_d6_fraction_unfiltered(
    reader: &mut bam::IndexedReader,
    tid: u32,
    diag_sites: &[Exon9DiagSite],
) -> Option<f64> {
    let mut d6_count = 0usize;
    let mut d7_count = 0usize;

    for site in diag_sites {
        if reader
            .fetch(bam::FetchDefinition::Region(
                tid as i32,
                site.pos - 1,
                site.pos,
            ))
            .is_err()
        {
            continue;
        }

        let mut pileups = reader.pileup();
        pileups.set_max_depth(i32::MAX as u32);

        for pileup_result in pileups {
            let pileup = match pileup_result {
                Ok(p) => p,
                Err(_) => continue,
            };
            if pileup.pos() as i64 + 1 != site.pos {
                continue;
            }

            for alignment in pileup.alignments() {
                let record = alignment.record();
                if record.is_secondary()
                    || record.is_supplementary()
                    || record.is_duplicate()
                    || record.is_unmapped()
                {
                    continue;
                }
                if record.mapq() < 10 {
                    continue;
                }
                if alignment.is_del() || alignment.is_refskip() {
                    continue;
                }

                let qpos = match alignment.qpos() {
                    Some(p) => p,
                    None => continue,
                };

                let seq = record.seq().as_bytes();
                let allele_len = site.d7_allele.len();
                if qpos + allele_len > seq.len() {
                    continue;
                }

                let read_seq = &seq[qpos..qpos + allele_len];

                let is_d6 = site
                    .d6_alleles
                    .iter()
                    .any(|d6| read_seq.eq_ignore_ascii_case(d6));
                let is_d7 = read_seq.eq_ignore_ascii_case(&site.d7_allele);

                if is_d6 {
                    d6_count += 1;
                } else if is_d7 {
                    d7_count += 1;
                }
            }
        }
    }

    let total = d6_count + d7_count;
    if total < 3 {
        return None;
    }

    log::info!(
        "exon9_d6_fraction_unfiltered: d6={}, d7={}, frac={:.3}",
        d6_count,
        d7_count,
        d6_count as f64 / total as f64,
    );

    Some(d6_count as f64 / total as f64)
}

/// Suggest structural correction from phased exon9 D6 fractions.
fn suggest_from_phased(
    frac_alt: f64,
    frac_ref: f64,
    _n_alt: usize,
    _n_ref: usize,
    _total_cn: u32,
) -> Option<String> {
    // If BOTH haplotypes show high D6 (no gene conversion on either),
    // then exon9hyb is a false positive — it's really a duplication.
    let (higher, lower) = if frac_alt >= frac_ref {
        (frac_alt, frac_ref)
    } else {
        (frac_ref, frac_alt)
    };

    log::info!(
        "suggest_from_phased: higher={:.3}, lower={:.3}",
        higher,
        lower,
    );

    if higher > 0.85 && lower > 0.85 {
        // Both haplotypes are D6 at exon9 → no gene conversion → duplication
        Some("dup".to_string())
    } else {
        // At least one haplotype shows D7 → confirms exon9hyb (tandem)
        None
    }
}

/// Suggest structural correction from unphased (aggregate) exon9 D6 fraction.
/// Uses CN-adjusted thresholds: for total_cn copies, a duplication has
/// (total_cn - 2) D6 copies + 2 D7-origin copies... but actually for dup
/// ALL copies are D6, while for exon9hyb some are D7.
///
/// Expected D6 fractions:
///   - Pure dup (no gene conversion): D6_frac ≈ 1.0
///   - exon9hyb with 1 tandem: D6_frac ≈ (total_cn - 1) / total_cn
///     e.g., CN=5: 4/5 = 0.80 ... but tandem *36 has D7 at exon9,
///     so actually: for CN=5 with one *36: D6_frac ≈ 4/5 = 0.80
///     Wait — the spacer_cn vs exon9 CN difference is what matters.
///
/// Simpler approach: if D6_frac is very high (>0.90), it's likely dup.
/// If D6_frac shows clear D7 presence (<0.80), it's exon9hyb.
fn suggest_from_unphased(frac: f64, total_cn: u32) -> Option<String> {
    // CN-adjusted midpoint between dup and exon9hyb expectations
    // For dup: all copies D6 → frac ≈ 1.0
    // For exon9hyb: one copy has D7 exon9 → frac ≈ (total_cn - 1) / total_cn
    let exon9hyb_expected = if total_cn > 1 {
        (total_cn - 1) as f64 / total_cn as f64
    } else {
        0.5
    };
    let midpoint = (1.0 + exon9hyb_expected) / 2.0;

    log::info!(
        "suggest_from_unphased: frac={:.3}, total_cn={}, exon9hyb_expected={:.3}, midpoint={:.3}",
        frac,
        total_cn,
        exon9hyb_expected,
        midpoint,
    );

    if frac > midpoint {
        // D6 fraction higher than midpoint → looks like pure dup
        Some("dup".to_string())
    } else {
        // D6 fraction lower → confirms exon9hyb
        None
    }
}

/// Run haplotype-resolved structural analysis.
///
/// Called when cnvtag contains "exon9hyb" (potential tandem *36+*10).
/// Uses het variants to phase reads, then checks per-haplotype D6/D7 at exon9.
///
/// If BOTH haplotypes show high D6 fraction at exon9 (no gene conversion),
/// the "exon9hyb" call is likely a false positive from alignment bias,
/// and the true structure is a simple duplication.
pub fn resolve_structure(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    var_list: &[String],
    var_alt: &[usize],
    var_ref: &[usize],
    total_cn: u32,
) -> HaplotypeResolveResult {
    let tid = match reader.header().tid(nchr.as_bytes()) {
        Some(t) => t,
        None => {
            return HaplotypeResolveResult {
                n_het_sites: 0,
                n_reads_hap_alt: 0,
                n_reads_hap_ref: 0,
                exon9_d6_frac_hap_alt: None,
                exon9_d6_frac_hap_ref: None,
                suggested_cnv: None,
            };
        }
    };

    // Step 1: Find het sites for phasing
    let het_sites = find_het_sites(reader, nchr, var_list, var_alt, var_ref);

    log::info!(
        "haplotype_resolve: found {} het sites for phasing",
        het_sites.len()
    );

    if het_sites.is_empty() {
        return HaplotypeResolveResult {
            n_het_sites: 0,
            n_reads_hap_alt: 0,
            n_reads_hap_ref: 0,
            exon9_d6_frac_hap_alt: None,
            exon9_d6_frac_hap_ref: None,
            suggested_cnv: None,
        };
    }

    // Step 2: Assign reads to haplotypes
    let (hap_alt, hap_ref) = assign_haplotypes(&het_sites);

    log::info!(
        "haplotype_resolve: {} reads on hap_alt, {} on hap_ref",
        hap_alt.len(),
        hap_ref.len()
    );

    if hap_alt.len() < 5 || hap_ref.len() < 5 {
        return HaplotypeResolveResult {
            n_het_sites: het_sites.len(),
            n_reads_hap_alt: hap_alt.len(),
            n_reads_hap_ref: hap_ref.len(),
            exon9_d6_frac_hap_alt: None,
            exon9_d6_frac_hap_ref: None,
            suggested_cnv: None,
        };
    }

    // Step 3: Check D6/D7 at exon9 for each haplotype
    let diag_sites = exon9_diagnostic_sites();
    let d6_frac_alt = exon9_d6_fraction(reader, tid, &hap_alt, &diag_sites);
    let d6_frac_ref = exon9_d6_fraction(reader, tid, &hap_ref, &diag_sites);

    log::info!(
        "haplotype_resolve: exon9 D6 fraction — hap_alt={:?}, hap_ref={:?}",
        d6_frac_alt,
        d6_frac_ref,
    );

    // Step 4: Determine structural suggestion
    let suggested_cnv = match (d6_frac_alt, d6_frac_ref) {
        (Some(frac_alt), Some(frac_ref)) => {
            suggest_from_phased(frac_alt, frac_ref, hap_alt.len(), hap_ref.len(), total_cn)
        }
        _ => {
            // Fallback: phased reads don't reach exon9 (het sites too far).
            // Check ALL reads at exon9 for D6/D7 without phasing.
            // Use CN-adjusted thresholds:
            //   For total_cn=5 with exon9hyb: expect D6_frac = 2/5 = 0.40
            //   For total_cn=5 with dup:      expect D6_frac = 3/5 = 0.60
            log::info!("haplotype_resolve: phased exon9 data insufficient, trying unphased fallback");
            let d6_frac_all = exon9_d6_fraction_unfiltered(reader, tid, &diag_sites);
            if let Some(frac) = d6_frac_all {
                suggest_from_unphased(frac, total_cn)
            } else {
                None
            }
        }
    };

    HaplotypeResolveResult {
        n_het_sites: het_sites.len(),
        n_reads_hap_alt: hap_alt.len(),
        n_reads_hap_ref: hap_ref.len(),
        exon9_d6_frac_hap_alt: d6_frac_alt,
        exon9_d6_frac_hap_ref: d6_frac_ref,
        suggested_cnv,
    }
}

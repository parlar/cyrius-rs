//! Physical phasing disambiguator for ambiguous diplotype calls.
//!
//! When match_star returns "more_than_one_match" for CN=2, this module uses
//! read-level phasing to pick the correct diplotype. It examines which defining
//! variants co-occur on the same reads (cis) vs different reads (trans).

use crate::caller::match_star_allele::main_allele_name;
use crate::depth_calling::snp_count::passing_read;
use crate::types::StarCombinations;
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::{HashMap, HashSet};

/// Variant info from the target variant file, mapping variant names to pileup positions
/// and ALT/REF base patterns.
#[derive(Debug, Clone)]
pub struct VariantLookup {
    /// variant_name -> (pileup_pos_0based, alt_bases, ref_bases)
    entries: HashMap<String, (i64, String, String)>,
}

impl VariantLookup {
    /// Parse the target variant file content into a lookup table.
    /// NOTE: this is a module-local variant-name parser with rules specific to
    /// physical phasing; it is not a general-purpose parser.
    pub fn from_target_file(content: &str) -> Self {
        let mut entries = HashMap::new();
        for line in content.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 7 {
                continue;
            }
            let pos: i64 = match fields[1].parse::<i64>() {
                Ok(p) => p - 1, // convert 1-based to 0-based
                Err(_) => continue,
            };
            let alt_bases = fields[2].to_string();
            // REF may have comma-separated alternatives; take the first
            let ref_bases = fields[4].split(',').next().unwrap_or("").to_string();
            let variant_name = fields[6].to_string();
            entries.insert(variant_name, (pos, alt_bases, ref_bases));
        }
        VariantLookup { entries }
    }

    /// Look up a variant by name as the pileup position and ALT/REF patterns
    /// the phasing reads are compared against.
    fn resolve(&self, name: &str) -> Option<ResolvedVariant> {
        self.entries.get(name).map(|(pos, alt, ref_b)| ResolvedVariant {
            pos_0based: *pos,
            alt_bases: alt.clone(),
            ref_bases: ref_b.clone(),
        })
    }
}

/// MAPQ a fragment needs on at least one of its reads before a both-REF
/// observation is trusted to come from CYP2D6 rather than a CYP2D7 read placed
/// on CYP2D6 at random (MAPQ 0). The floor the resolver's fragment scorer
/// (`diplotype_caller::score_candidates`) already uses.
pub const RR_MIN_MAPQ: u8 = 10;

/// Fragment tallies for one variant pair. A fragment is every read sharing a
/// qname, merged as `check_phasing` merges them; "first" and "second" are the
/// pair's two variants in the order they were passed.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PhaseCounts {
    /// ALT at both sites.
    pub aa: usize,
    /// ALT at the first site, REF at the second.
    pub ar: usize,
    /// REF at the first site, ALT at the second.
    pub ra: usize,
    /// REF at both sites, whatever the MAPQ.
    pub rr_all: usize,
    /// REF at both sites, with at least one read at MAPQ >= `RR_MIN_MAPQ`.
    pub rr_anchored: usize,
}

impl PhaseCounts {
    /// Fragments that put the two variants on the same haplotype: ALT at both
    /// sites, or REF at both from an anchored fragment. Phasing only runs on
    /// cn2 calls, so when the two ALTs are cis the other D6 copy carries
    /// neither and shows up as both-REF fragments; when they are trans, no
    /// D6 copy is REF at both sites.
    pub fn cis(&self) -> usize {
        self.aa + self.rr_anchored
    }

    /// Fragments that put the two variants on different haplotypes.
    pub fn trans(&self) -> usize {
        self.ar + self.ra
    }
}

/// One variant pair as counted for phasing: the two variant names, their
/// 0-based pileup positions and the fragment tallies.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairEvidence {
    pub va: String,
    pub vb: String,
    pub pos_a: i64,
    pub pos_b: i64,
    pub counts: PhaseCounts,
}

impl PairEvidence {
    /// Distance between the two variants in bp.
    pub fn distance_bp(&self) -> u64 {
        (self.pos_a - self.pos_b).unsigned_abs()
    }
}

/// A diplotype option with raw (suballele) star names and clean (main allele) names.
#[derive(Debug, Clone)]
struct Diplotype {
    raw_hap1: String,
    raw_hap2: String,
    clean_hap1: String,
    clean_hap2: String,
}

/// Parse raw candidates (e.g., ["*1_*21.002", "*2.001_*21"]) and clean call
/// (e.g., "*1/*21;*2/*21") into diplotype options.
fn parse_diplotypes(raw_candidates: &[String], clean_call: &str) -> Vec<Diplotype> {
    let clean_parts: Vec<&str> = clean_call.split(';').collect();

    // Build diplotypes from raw candidates, matching each with the clean
    // counterpart that has the same main alleles. The two lists are sorted
    // differently, so pairing by position would mislabel them.
    let mut diplotypes = Vec::new();
    for raw in raw_candidates {
        let raw_parts: Vec<&str> = raw.split('_').collect();
        if raw_parts.len() != 2 {
            continue;
        }
        let main1 = main_allele_name(raw_parts[0]);
        let main2 = main_allele_name(raw_parts[1]);
        let mut raw_mains = [main1, main2];
        raw_mains.sort();
        let matched = clean_parts.iter().find_map(|part| {
            let cp: Vec<&str> = part.split('/').collect();
            if cp.len() != 2 {
                return None;
            }
            let mut clean_mains = [cp[0], cp[1]];
            clean_mains.sort();
            (clean_mains == raw_mains).then(|| (cp[0].to_string(), cp[1].to_string()))
        });
        // Fallback: the stripped suballele names
        let (clean_h1, clean_h2) =
            matched.unwrap_or_else(|| (main1.to_string(), main2.to_string()));
        diplotypes.push(Diplotype {
            raw_hap1: raw_parts[0].to_string(),
            raw_hap2: raw_parts[1].to_string(),
            clean_hap1: clean_h1,
            clean_hap2: clean_h2,
        });
    }
    diplotypes
}

/// Get the set of defining variant names for a star allele (using raw suballele name).
/// Returns variant names like ["g.42128945C>T", "g.42129033G>A"].
fn get_star_variants(star: &str, dstar: &HashMap<String, String>) -> Vec<String> {
    if let Some(var_key) = dstar.get(star) {
        if var_key == "NA" {
            return vec![];
        }
        return var_key.split('_').map(|s| s.to_string()).collect();
    }
    // Anything not in the table (including *1, the reference) has no defining
    // variants for phasing purposes.
    vec![]
}


/// Resolved variant info for phasing: pileup position (0-based) + ALT/REF patterns.
#[derive(Debug, Clone)]
struct ResolvedVariant {
    pos_0based: i64,
    alt_bases: String,
    ref_bases: String,
}

/// Half-open BAM fetch window `[start, end)` covering both variant positions
/// for the single-fetch phasing path.
///
/// `end` is one past the later position: the region is half-open, so stopping
/// at the position itself drops any read whose leftmost mapped base is exactly
/// the later variant.
///
/// `start` is `min(pos) - 1`, one base wider than it needs to be. An htslib
/// region query returns every read *overlapping* `[beg, end)`, not only reads
/// starting inside it — a query at a single position comes back with reads
/// that start ~150 bp earlier — so `start = min(pos)` would already fetch
/// everything covering the earlier variant. The `-1` is a harmless leftover of
/// the very misconception that produced the `end` bug. It is kept deliberately:
/// narrowing a live fetch window is a behaviour change, not a comment fix.
///
/// `start` is clamped at 0 so a variant at position 0 cannot produce a negative
/// coordinate. That clamp is unreachable on real CYP2D6 data, where positions
/// are ~42 Mbp; it guards a bad coordinate from a malformed variant file, not
/// any live case.
///
/// Returns `(start, end)`.
fn single_fetch_window(pos1: i64, pos2: i64) -> (i64, i64) {
    let start = (pos1.min(pos2) - 1).max(0);
    let end = pos1.max(pos2) + 1;
    (start, end)
}

/// The read filters every phasing count applies: primary, not a duplicate,
/// and `passing_read`.
fn phase_read_passes(record: &bam::Record) -> bool {
    !(record.is_secondary() || record.is_supplementary() || record.is_duplicate())
        && passing_read(record, false, false)
}

/// What one fragment (all reads sharing a qname) showed at a pair's two sites.
#[derive(Debug, Default)]
struct Fragment {
    /// ALT (`Some(true)`), REF (`Some(false)`) or not seen at the first site.
    allele1: Option<bool>,
    /// The same at the second site.
    allele2: Option<bool>,
    /// Highest MAPQ over every read of this qname that passed the filters.
    max_mapq: u8,
}

/// Merge one fetched read's observations at both sites into its fragment.
fn observe_read(
    fragments: &mut HashMap<String, Fragment>,
    record: &bam::Record,
    rv1: &ResolvedVariant,
    rv2: &ResolvedVariant,
) {
    let read_name = String::from_utf8_lossy(record.qname()).to_string();
    let seq = record.seq().as_bytes();
    let allele1 = check_allele_pileup(record, &seq, rv1);
    let allele2 = check_allele_pileup(record, &seq, rv2);
    let entry = fragments.entry(read_name).or_default();
    if allele1.is_some() {
        entry.allele1 = allele1;
    }
    if allele2.is_some() {
        entry.allele2 = allele2;
    }
    entry.max_mapq = entry.max_mapq.max(record.mapq());
}

/// Check if two variants at nearby positions co-occur on the same reads (cis)
/// or on different reads (trans).
///
/// When `use_read_pairs` is true and positions are >120bp apart, performs
/// two separate fetches so that mate pairs spanning the distance can both
/// contribute observations. Read pairs share the same qname, so a mate
/// covering pos1 and its partner covering pos2 will be linked.
///
/// Returns the fragment tallies; all zero when the contig is missing. A failed
/// fetch or a BAM read error is an `Err`: partial tallies must not look measured.
fn check_phasing(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    rv1: &ResolvedVariant,
    rv2: &ResolvedVariant,
    use_read_pairs: bool,
) -> Result<PhaseCounts, String> {
    let pos1 = rv1.pos_0based;
    let pos2 = rv2.pos_0based;
    let tid = match reader.header().tid(nchr.as_bytes()) {
        Some(t) => t,
        None => return Ok(PhaseCounts::default()),
    };

    let mut fragments: HashMap<String, Fragment> = HashMap::new();

    let distance = (pos1 - pos2).unsigned_abs();
    let use_two_fetches = use_read_pairs && distance > 120;

    if use_two_fetches {
        // Fetch reads overlapping pos1, then pos2. Use a wide window upstream
        // because reads starting up to ~150bp before a position can still
        // cover it.
        let read_len = 150i64;
        for pos in [pos1, pos2] {
            if reader
                .fetch(bam::FetchDefinition::Region(
                    tid as i32,
                    (pos - read_len).max(0),
                    pos + 1,
                ))
                .is_err()
            {
                return Err(format!("BAM fetch failed at {}:{}", nchr, pos1.min(pos2)));
            }
            for record_result in reader.records() {
                let record = match record_result {
                    Ok(r) => r,
                    Err(e) => return Err(format!("BAM read failed at {}:{}: {}", nchr, pos1.min(pos2), e)),
                };
                if !phase_read_passes(&record) {
                    continue;
                }
                observe_read(&mut fragments, &record, rv1, rv2);
            }
        }
    } else {
        // Single fetch spanning both positions
        let (start, end) = single_fetch_window(pos1, pos2);
        if reader
            .fetch(bam::FetchDefinition::Region(tid as i32, start, end))
            .is_err()
        {
            return Err(format!("BAM fetch failed at {}:{}", nchr, pos1.min(pos2)));
        }

        for record_result in reader.records() {
            let record = match record_result {
                Ok(r) => r,
                Err(e) => return Err(format!("BAM read failed at {}:{}: {}", nchr, pos1.min(pos2), e)),
            };
            if !phase_read_passes(&record) {
                continue;
            }
            observe_read(&mut fragments, &record, rv1, rv2);
        }
    }

    let mut counts = PhaseCounts::default();
    for f in fragments.values() {
        match (f.allele1, f.allele2) {
            (Some(true), Some(true)) => counts.aa += 1,  // both ALT → same haplotype
            (Some(true), Some(false)) => counts.ar += 1, // one ALT, one REF → different haplotypes
            (Some(false), Some(true)) => counts.ra += 1,
            (Some(false), Some(false)) => {
                counts.rr_all += 1;
                if f.max_mapq >= RR_MIN_MAPQ {
                    counts.rr_anchored += 1;
                }
            }
            _ => {}
        }
    }

    Ok(counts)
}

/// Count one variant pair: resolve both names, fetch, tally. `None` when
/// either variant has no lookup entry, `Err` on a BAM fetch or read error.
/// The one counting path both
/// `disambiguate` and `count_pair` use.
fn pair_evidence(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    va: &str,
    vb: &str,
    variant_lookup: &VariantLookup,
    use_read_pairs: bool,
) -> Result<Option<PairEvidence>, String> {
    let rv_a = match variant_lookup.resolve(va) {
        Some(rv) => rv,
        None => {
            log::debug!("phase_disambiguate: no lookup for variant {}", va);
            return Ok(None);
        }
    };
    let rv_b = match variant_lookup.resolve(vb) {
        Some(rv) => rv,
        None => {
            log::debug!("phase_disambiguate: no lookup for variant {}", vb);
            return Ok(None);
        }
    };
    let counts = check_phasing(reader, nchr, &rv_a, &rv_b, use_read_pairs)?;
    Ok(Some(PairEvidence {
        va: va.to_string(),
        vb: vb.to_string(),
        pos_a: rv_a.pos_0based,
        pos_b: rv_b.pos_0based,
        counts,
    }))
}

/// Count one variant pair exactly as `disambiguate` does on its default
/// (single-fetch) path. `Ok(None)` when either variant has no lookup entry;
/// `Err` on a BAM fetch or read error.
/// For audits: the counts come from the shipped code, not a port.
pub fn count_pair(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    va: &str,
    vb: &str,
    variant_lookup: &VariantLookup,
) -> Result<Option<PairEvidence>, String> {
    pair_evidence(reader, nchr, va, vb, variant_lookup, false)
}

/// Count the reads overlapping one variant's position as `(alt, ref, none)`,
/// with the phasing read filters and `check_allele_pileup`. `None` when the
/// variant has no lookup entry or when the fetch or a read fails; all zero when
/// the contig is missing. Counts reads, not fragments.
pub fn count_site(
    reader: &mut bam::IndexedReader,
    nchr: &str,
    variant: &str,
    variant_lookup: &VariantLookup,
) -> Option<(usize, usize, usize)> {
    let rv = variant_lookup.resolve(variant)?;
    let tid = match reader.header().tid(nchr.as_bytes()) {
        Some(t) => t,
        None => return Some((0, 0, 0)),
    };
    if reader
        .fetch(bam::FetchDefinition::Region(
            tid as i32,
            rv.pos_0based,
            rv.pos_0based + 1,
        ))
        .is_err()
    {
        return None;
    }
    let (mut alt, mut ref_count, mut none) = (0, 0, 0);
    for record_result in reader.records() {
        let record = match record_result {
            Ok(r) => r,
            Err(_) => return None,
        };
        if !phase_read_passes(&record) {
            continue;
        }
        let seq = record.seq().as_bytes();
        match check_allele_pileup(&record, &seq, &rv) {
            Some(true) => alt += 1,
            Some(false) => ref_count += 1,
            None => none += 1,
        }
    }
    Some((alt, ref_count, none))
}

/// Check whether a read carries the ALT allele at a pileup position.
/// Reads multiple bases from the read starting at rv.pos_0based and compares
/// against ALT and REF patterns from the target variant file.
///
/// A pattern matches only when the read's bases start with the whole of it
/// (case-insensitive), so a read too short to hold a pattern cannot match it:
/// a read ending inside a repeat is not ALT just because what it holds is a
/// prefix of the ALT. If both patterns match (one is a prefix of the other),
/// the longer one wins.
///
/// Returns Some(true) for ALT, Some(false) for REF, None if the position is
/// not covered or neither pattern matches in full.
fn check_allele_pileup(
    record: &rust_htslib::bam::Record,
    seq: &[u8],
    rv: &ResolvedVariant,
) -> Option<bool> {
    let check_len = rv.alt_bases.len().max(rv.ref_bases.len());
    let bases = read_bases_at_pos(record, seq, rv.pos_0based, check_len)?;
    let bases_upper = bases.to_uppercase();
    let alt_upper = rv.alt_bases.to_uppercase();
    let ref_upper = rv.ref_bases.to_uppercase();

    match (
        bases_upper.starts_with(&alt_upper),
        bases_upper.starts_with(&ref_upper),
    ) {
        (true, true) => Some(alt_upper.len() > ref_upper.len()),
        (true, false) => Some(true),
        (false, true) => Some(false),
        (false, false) => None, // holds neither whole pattern
    }
}

/// Read `n` bases from a read starting at a 0-based genomic position.
/// Walks the CIGAR to find the query offset, then extracts bases.
fn read_bases_at_pos(
    record: &rust_htslib::bam::Record,
    seq: &[u8],
    genomic_pos_0based: i64,
    n: usize,
) -> Option<String> {
    let record_start = record.pos(); // 0-based
    let ref_offset = genomic_pos_0based - record_start;
    if ref_offset < 0 {
        return None;
    }

    let mut ref_pos = 0i64;
    let mut query_pos = 0usize;

    for &cigar_elem in record.cigar().iter() {
        use rust_htslib::bam::record::Cigar::*;
        match cigar_elem {
            Match(len) | Equal(len) | Diff(len) => {
                let len = len as i64;
                if ref_pos + len > ref_offset {
                    let qp = query_pos + (ref_offset - ref_pos) as usize;
                    if qp + n <= seq.len() {
                        let bases: String = seq[qp..qp + n]
                            .iter()
                            .map(|&b| b as char)
                            .collect();
                        return Some(bases);
                    } else if qp < seq.len() {
                        let bases: String = seq[qp..seq.len()]
                            .iter()
                            .map(|&b| b as char)
                            .collect();
                        return Some(bases);
                    }
                    return None;
                }
                ref_pos += len;
                query_pos += len as usize;
            }
            Ins(len) => { query_pos += len as usize; }
            Del(len) => {
                let len = len as i64;
                if ref_pos + len > ref_offset { return None; }
                ref_pos += len;
            }
            SoftClip(len) => { query_pos += len as usize; }
            HardClip(_) | Pad(_) => {}
            RefSkip(len) => { ref_pos += len as i64; }
        }
    }
    None
}

/// What `disambiguate` found.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PhaseOutcome {
    /// The diplotype the phase evidence picked, or `None` if it did not decide.
    pub resolved: Option<String>,
    /// One entry per differentiating pair that had both lookups, in the order
    /// evaluated, whether or not its counts were decisive. Empty when the
    /// step stopped before counting (not two options, no differentiating
    /// pair, or the BAM could not be opened).
    pub pairs: Vec<PairEvidence>,
}

/// Attempt to disambiguate an ambiguous CN=2 diplotype call using physical phasing.
///
/// `raw_candidates` are the raw star allele combinations (with suballeles), e.g.,
/// ["*1_*21.002", "*2.001_*21"]. `clean_call` is the main-allele form like "*1/*21;*2/*21".
///
/// Returns the resolved call (`None` unless disambiguation succeeds) together
/// with the evidence for every pair counted.
pub fn disambiguate(
    clean_call: &str,
    raw_candidates: &[String],
    star_combinations: &StarCombinations,
    variant_lookup: &VariantLookup,
    bam_path: &str,
    nchr: &str,
    reference_fasta: Option<&str>,
    index_name: Option<&str>,
    use_read_pairs: bool,
) -> PhaseOutcome {
    let diplotypes = parse_diplotypes(raw_candidates, clean_call);
    if diplotypes.len() != 2 {
        log::debug!(
            "phase_disambiguate: expected 2 diplotype options, got {} (raw: {:?})",
            diplotypes.len(),
            raw_candidates
        );
        return PhaseOutcome::default();
    }

    let dt1 = &diplotypes[0];
    let dt2 = &diplotypes[1];

    let dstar = &star_combinations.dstar;

    // Get variants for each haplotype using raw (suballele) names
    let vars_dt1_h1 = get_star_variants(&dt1.raw_hap1, dstar);
    let vars_dt1_h2 = get_star_variants(&dt1.raw_hap2, dstar);
    let vars_dt2_h1 = get_star_variants(&dt2.raw_hap1, dstar);
    let vars_dt2_h2 = get_star_variants(&dt2.raw_hap2, dstar);

    log::debug!(
        "phase_disambiguate: dt1={}({})/{}({}) dt2={}({})/{}({})",
        dt1.raw_hap1,
        vars_dt1_h1.join(","),
        dt1.raw_hap2,
        vars_dt1_h2.join(","),
        dt2.raw_hap1,
        vars_dt2_h1.join(","),
        dt2.raw_hap2,
        vars_dt2_h2.join(","),
    );

    // Collect all unique variants
    let all_vars: HashSet<&str> = vars_dt1_h1
        .iter()
        .chain(vars_dt1_h2.iter())
        .chain(vars_dt2_h1.iter())
        .chain(vars_dt2_h2.iter())
        .map(|s| s.as_str())
        .collect();

    if all_vars.len() < 2 {
        log::debug!("phase_disambiguate: too few variants ({}) to phase", all_vars.len());
        return PhaseOutcome::default();
    }

    let dt1_h1_set: HashSet<&str> = vars_dt1_h1.iter().map(|s| s.as_str()).collect();
    let dt1_h2_set: HashSet<&str> = vars_dt1_h2.iter().map(|s| s.as_str()).collect();
    let dt2_h1_set: HashSet<&str> = vars_dt2_h1.iter().map(|s| s.as_str()).collect();
    let dt2_h2_set: HashSet<&str> = vars_dt2_h2.iter().map(|s| s.as_str()).collect();

    // Find pairs of variants where the two options predict different phasing.
    //
    // `all_vars` is a HashSet, and std's RandomState seeds its iteration order
    // per process, so draining it unsorted makes both the order of
    // `phasing_pairs` and the (va, vb) assignment within each pair vary from
    // run to run. Every consumer downstream happens to be order-symmetric
    // today — `are_cis_in_diplotype` is symmetric in a/b, `single_fetch_window`
    // is symmetric, the two-fetch branch merges both fetches so swapping
    // rv1/rv2 only swaps `allele1`/`allele2` (and `ar`/`ra`) and cis/trans
    // survive that swap,
    // and the scores are order-independent sums — so the scores and the
    // resolved call are stable. The debug log is not, and any future
    // order-sensitive tiebreak would flip calls between runs. Sort for a
    // reproducible order.
    let mut all_vars_vec: Vec<&str> = all_vars.into_iter().collect();
    all_vars_vec.sort_unstable();
    let mut phasing_pairs: Vec<(&str, &str, bool)> = Vec::new(); // (var_a, var_b, cis_in_dt1)

    for i in 0..all_vars_vec.len() {
        for j in (i + 1)..all_vars_vec.len() {
            let va = all_vars_vec[i];
            let vb = all_vars_vec[j];

            let cis_dt1 = are_cis_in_diplotype(va, vb, &dt1_h1_set, &dt1_h2_set);
            let cis_dt2 = are_cis_in_diplotype(va, vb, &dt2_h1_set, &dt2_h2_set);

            match (cis_dt1, cis_dt2) {
                (Some(true), Some(false)) => phasing_pairs.push((va, vb, true)),
                (Some(false), Some(true)) => phasing_pairs.push((va, vb, false)),
                _ => {} // Same phasing in both options, or indeterminate
            }
        }
    }

    if phasing_pairs.is_empty() {
        log::debug!("phase_disambiguate: no differentiating variant pairs found");
        return PhaseOutcome::default();
    }

    log::info!(
        "phase_disambiguate: found {} differentiating variant pairs for {}",
        phasing_pairs.len(),
        clean_call
    );

    // Open BAM and check phasing for each pair
    let mut reader = match crate::depth_calling::utilities::open_alignment_file_with_index(
        bam_path,
        reference_fasta,
        index_name,
    ) {
        Ok(r) => r,
        Err(e) => {
            log::warn!("phase_disambiguate: failed to open BAM: {}", e);
            return PhaseOutcome::default();
        }
    };

    let mut dt1_score = 0i32;
    let mut dt2_score = 0i32;
    let mut pairs: Vec<PairEvidence> = Vec::new();

    for &(va, vb, cis_in_dt1) in &phasing_pairs {
        let evidence =
            match pair_evidence(&mut reader, nchr, va, vb, variant_lookup, use_read_pairs) {
                Ok(Some(e)) => e,
                Ok(None) => continue,
                Err(e) => {
                    log::warn!("phase_disambiguate: {e}; phase left unresolved");
                    return PhaseOutcome { resolved: None, pairs };
                }
            };
        let cis_count = evidence.counts.cis();
        let trans_count = evidence.counts.trans();
        let (pos_a, pos_b, counts) = (evidence.pos_a, evidence.pos_b, evidence.counts);
        pairs.push(evidence);

        if cis_count + trans_count == 0 {
            continue;
        }

        log::debug!(
            "phase_disambiguate: {}(pos={}) vs {}(pos={}): cis={}, trans={}, expected_cis_dt1={}, rr_all={}, rr_anchored={}",
            va, pos_a, vb, pos_b, cis_count, trans_count, cis_in_dt1,
            counts.rr_all, counts.rr_anchored
        );

        // Only count a pair's evidence when the imbalance is clear (>2:1 ratio)
        let total = cis_count + trans_count;
        let majority = cis_count.max(trans_count);
        if total < 3 || (majority as f64) / (total as f64) < 0.7 {
            log::debug!(
                "phase_disambiguate: skipping weak pair (cis={}, trans={})",
                cis_count, trans_count
            );
            continue;
        }

        if cis_count > trans_count {
            if cis_in_dt1 {
                dt1_score += (cis_count - trans_count) as i32;
            } else {
                dt2_score += (cis_count - trans_count) as i32;
            }
        } else {
            if cis_in_dt1 {
                dt2_score += (trans_count - cis_count) as i32;
            } else {
                dt1_score += (trans_count - cis_count) as i32;
            }
        }
    }

    log::info!(
        "phase_disambiguate: dt1_score={}, dt2_score={} for {}",
        dt1_score, dt2_score, clean_call
    );

    // Require clear winner: score >= 3 and at least 3x the other
    let min_score = 3;
    let resolved = if dt1_score >= min_score && (dt2_score == 0 || dt1_score >= 3 * dt2_score) {
        let resolved = format!("{}/{}", dt1.clean_hap1, dt1.clean_hap2);
        log::info!("phase_disambiguate: resolved {} -> {}", clean_call, resolved);
        Some(resolved)
    } else if dt2_score >= min_score && (dt1_score == 0 || dt2_score >= 3 * dt1_score) {
        let resolved = format!("{}/{}", dt2.clean_hap1, dt2.clean_hap2);
        log::info!("phase_disambiguate: resolved {} -> {}", clean_call, resolved);
        Some(resolved)
    } else {
        log::info!("phase_disambiguate: inconclusive for {}", clean_call);
        None
    };
    PhaseOutcome { resolved, pairs }
}

/// Check if two variants are on the same haplotype (cis) in a given diplotype.
/// Returns Some(true) for cis, Some(false) for trans, None if indeterminate.
fn are_cis_in_diplotype(
    var_a: &str,
    var_b: &str,
    hap1_vars: &HashSet<&str>,
    hap2_vars: &HashSet<&str>,
) -> Option<bool> {
    let a_on_h1 = hap1_vars.contains(var_a);
    let a_on_h2 = hap2_vars.contains(var_a);
    let b_on_h1 = hap1_vars.contains(var_b);
    let b_on_h2 = hap2_vars.contains(var_b);

    // Homozygous variant (on both haplotypes) is not informative for phasing
    if (a_on_h1 && a_on_h2) || (b_on_h1 && b_on_h2) {
        return None;
    }

    match (a_on_h1, b_on_h1, a_on_h2, b_on_h2) {
        (true, true, false, false) => Some(true),  // Both on hap1
        (false, false, true, true) => Some(true),  // Both on hap2
        (true, false, false, true) => Some(false), // a on hap1, b on hap2
        (false, true, true, false) => Some(false), // a on hap2, b on hap1
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_diplotypes_labels_each_raw_candidate_with_its_own_clean_form() {
        // Raw order (*1 before *10) differs from clean order ("*10_" sorts
        // before "*1_"), so pairing by position would swap the labels.
        let raw = vec!["*1_*4.004".to_string(), "*10_*4.012".to_string()];
        let dts = parse_diplotypes(&raw, "*10/*4;*1/*4");
        assert_eq!(dts.len(), 2);
        let by_raw: HashMap<&str, (&str, &str)> = dts
            .iter()
            .map(|d| (d.raw_hap1.as_str(), (d.clean_hap1.as_str(), d.clean_hap2.as_str())))
            .collect();
        assert_eq!(by_raw["*1"], ("*1", "*4"));
        assert_eq!(by_raw["*10"], ("*10", "*4"));
    }

    #[test]
    fn single_fetch_window_reaches_past_the_later_variant() {
        // The fetch region is half-open, so `end` must be one past the later
        // position; at `end == pos` a read whose leftmost mapped base IS the
        // later variant falls outside the window and is never fetched.
        let (start, end) = single_fetch_window(42_126_611, 42_126_666);
        // Exact equality: `start <= 42_126_610` also passes for -1, so it does
        // not guard the clamp, and `end > 42_126_666` also passes for
        // 99_999_999.
        assert_eq!(
            start, 42_126_610,
            "start must reach exactly one base before the earlier variant, got {start}"
        );
        assert_eq!(
            end, 42_126_667,
            "end must be exactly one past the later variant 42126666 so a read \
             starting exactly there is still fetched, got {end}"
        );
        // Argument order must not matter.
        assert_eq!(single_fetch_window(42_126_666, 42_126_611), (start, end));
    }

    #[test]
    fn single_fetch_window_clamps_a_negative_start() {
        // A variant at position 0 must not produce a negative fetch start.
        let (start, end) = single_fetch_window(0, 7);
        assert_eq!(start, 0, "start must be clamped to 0, got {start}");
        assert_eq!(end, 8, "end must be one past the later variant, got {end}");
    }

    /// 0-based position of the earlier of the fixture's two phasing variants.
    const FIXTURE_POS_EARLY: i64 = 1000;
    /// 0-based position of the later one: 100 bp on, so `distance <= 120` and
    /// `check_phasing` takes the single-fetch path under either flag setting.
    const FIXTURE_POS_LATE: i64 = 1100;
    /// 0-based start of the mate. With CIGAR 50M its last mapped base is
    /// exactly `FIXTURE_POS_EARLY`.
    const FIXTURE_MATE_POS: i64 = 951;

    /// Write an indexed BAM holding one two-record pair for the single-fetch
    /// phasing window.
    ///
    /// Record 1 (the mate) is 50M from `FIXTURE_MATE_POS`, so its last mapped
    /// base is exactly `FIXTURE_POS_EARLY`, where it carries the ALT. Record 2
    /// starts at exactly `FIXTURE_POS_LATE` — the later variant — and carries
    /// the ALT there. They share one qname, so the two observations land in the
    /// same `read_alleles` entry and tally as a single cis read, but only if
    /// the fetch window actually reaches record 2.
    ///
    /// NOTE: `bam::Record::new()` leaves `flags = 4` (unmapped), and htslib's
    /// `bam_endpos` ignores the CIGAR for an unmapped record, so an unflagged
    /// record spans only `[pos, pos + 1)` and silently drops out of region
    /// queries. The `set_flags` calls below are load-bearing.
    fn write_phasing_pair_bam(dir: &std::path::Path) -> std::path::PathBuf {
        let mut header = bam::header::Header::new();
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr22");
        sq.push_tag(b"LN", 50_818_468);
        header.push_record(&sq);
        let path = dir.join("phasing_pair.bam");
        {
            let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
            let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(50)]);
            let qual = vec![30u8; 50];
            let span = FIXTURE_POS_LATE + 50 - FIXTURE_MATE_POS;

            // The mate: ALT 'C' against REF 'G' at FIXTURE_POS_EARLY.
            let mut mate_seq = vec![b'G'; 50];
            mate_seq[(FIXTURE_POS_EARLY - FIXTURE_MATE_POS) as usize] = b'C';
            let mut mate = bam::Record::new();
            mate.set(b"pair1", Some(&cigar), &mate_seq, &qual);
            mate.set_tid(0);
            mate.set_pos(FIXTURE_MATE_POS);
            mate.set_mapq(60);
            mate.set_mtid(0);
            mate.set_mpos(FIXTURE_POS_LATE);
            mate.set_insert_size(span);
            // PAIRED | PROPER_PAIR | MATE_REVERSE — mapped, primary, not a dup.
            mate.set_flags(0x1 | 0x2 | 0x20);
            writer.write(&mate).unwrap();

            // The read starting exactly at the later variant: ALT 'T' against
            // REF 'A' at FIXTURE_POS_LATE, query offset 0.
            let mut read_seq = vec![b'A'; 50];
            read_seq[0] = b'T';
            let mut read = bam::Record::new();
            read.set(b"pair1", Some(&cigar), &read_seq, &qual);
            read.set_tid(0);
            read.set_pos(FIXTURE_POS_LATE);
            read.set_mapq(60);
            read.set_mtid(0);
            read.set_mpos(FIXTURE_MATE_POS);
            read.set_insert_size(-span);
            // PAIRED | PROPER_PAIR | REVERSE
            read.set_flags(0x1 | 0x2 | 0x10);
            writer.write(&read).unwrap();
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

    #[test]
    fn check_phasing_fetches_a_read_starting_exactly_at_the_later_variant() {
        // Finding #52, end to end. The single-fetch region used to be
        // [min(pos) - 1, max(pos)); it is half-open, so the read whose leftmost
        // mapped base is exactly the later variant was never fetched. Its
        // mate's ALT at the earlier variant then sat alone in the entry as
        // (Some(true), None) and was dropped as uninformative, so the pair
        // tallied (0, 0) instead of (1, 0).
        //
        // The two window-arithmetic tests above cannot catch that: reverting
        // the call site to `let end = pos1.max(pos2);` leaves both of them
        // green and only this test red. That is what makes the fix
        // demonstrably non-vacuous.
        let dir = std::env::temp_dir()
            .join(format!("cyrius_phase_window_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bam_path = write_phasing_pair_bam(&dir);
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();

        let rv_early = ResolvedVariant {
            pos_0based: FIXTURE_POS_EARLY,
            alt_bases: "C".to_string(),
            ref_bases: "G".to_string(),
        };
        let rv_late = ResolvedVariant {
            pos_0based: FIXTURE_POS_LATE,
            alt_bases: "T".to_string(),
            ref_bases: "A".to_string(),
        };

        let single = check_phasing(&mut reader, "chr22", &rv_early, &rv_late, false).unwrap();
        let single = (single.cis(), single.trans());
        // 100 bp apart, so --phase-readpair takes the same single-fetch path.
        let readpair = check_phasing(&mut reader, "chr22", &rv_early, &rv_late, true).unwrap();
        let readpair = (readpair.cis(), readpair.trans());

        // Clean up before the assertions so a failure still leaves nothing behind.
        drop(reader);
        let _ = std::fs::remove_dir_all(&dir);

        assert_eq!(
            single,
            (1, 0),
            "the ALT-carrying read starting exactly at {FIXTURE_POS_LATE} plus \
             its mate's ALT at {FIXTURE_POS_EARLY} are one cis observation; got \
             (cis, trans) = {single:?}, which is what a window ending at the \
             later position returns"
        );
        assert_eq!(
            readpair,
            (1, 0),
            "the two variants are only 100 bp apart, so --phase-readpair must \
             take the same single-fetch path; got {readpair:?}"
        );
    }

    /// The fixture's two variants as a target-variant file: ALT C / REF G at
    /// `FIXTURE_POS_EARLY`, ALT T / REF A at `FIXTURE_POS_LATE` (1-based in
    /// the file).
    fn fixture_lookup() -> VariantLookup {
        VariantLookup::from_target_file(&format!(
            "#chr\tpos_CYP2D6\tbase_ALT\tpos_CYP2D7\tbase_REF\tvariant_type\tvariant_name\n\
             chr22\t{}\tC\t0\tG\tsnp\tg.early\n\
             chr22\t{}\tT\t0\tA\tsnp\tg.late\n",
            FIXTURE_POS_EARLY + 1,
            FIXTURE_POS_LATE + 1,
        ))
    }

    /// One mapped 50M record of a proper pair on chr22.
    fn fixture_record(qname: &str, pos: i64, seq: &[u8], mapq: u8, flags: u16) -> bam::Record {
        let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(50)]);
        let mut record = bam::Record::new();
        record.set(qname.as_bytes(), Some(&cigar), seq, &vec![30u8; seq.len()]);
        record.set_tid(0);
        record.set_pos(pos);
        record.set_mapq(mapq);
        record.set_mtid(0);
        // `Record::new()` leaves the unmapped flag set; see
        // `write_phasing_pair_bam`. These flags replace it.
        record.set_flags(flags);
        record
    }

    /// The two records of one fragment over the fixture's two variants, laid
    /// out like `write_phasing_pair_bam`: a mate whose last base is
    /// `FIXTURE_POS_EARLY` and a read whose first base is `FIXTURE_POS_LATE`,
    /// each carrying ALT (`true`) or REF (`false`) there, at the given MAPQs.
    fn fixture_fragment(
        qname: &str,
        early_alt: bool,
        late_alt: bool,
        mate_mapq: u8,
        read_mapq: u8,
    ) -> [bam::Record; 2] {
        let mut mate_seq = vec![b'G'; 50];
        mate_seq[(FIXTURE_POS_EARLY - FIXTURE_MATE_POS) as usize] =
            if early_alt { b'C' } else { b'G' };
        let mut read_seq = vec![b'A'; 50];
        read_seq[0] = if late_alt { b'T' } else { b'A' };
        [
            // PAIRED | PROPER_PAIR | MATE_REVERSE
            fixture_record(qname, FIXTURE_MATE_POS, &mate_seq, mate_mapq, 0x1 | 0x2 | 0x20),
            // PAIRED | PROPER_PAIR | REVERSE
            fixture_record(qname, FIXTURE_POS_LATE, &read_seq, read_mapq, 0x1 | 0x2 | 0x10),
        ]
    }

    /// Write `records` (sorted here by position) to an indexed chr22 BAM.
    fn write_fixture_bam(
        dir: &std::path::Path,
        file_name: &str,
        mut records: Vec<bam::Record>,
    ) -> std::path::PathBuf {
        let mut header = bam::header::Header::new();
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr22");
        sq.push_tag(b"LN", 50_818_468);
        header.push_record(&sq);
        let path = dir.join(file_name);
        records.sort_by_key(|r| r.pos());
        {
            let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
            for record in &records {
                writer.write(record).unwrap();
            }
        }
        bam::index::build(&path, None::<&std::path::PathBuf>, bam::index::Type::Bai, 1).unwrap();
        path
    }

    /// One fragment of every kind, with both-REF fragments either side of the
    /// `RR_MIN_MAPQ` floor.
    fn write_rr_fixture_bam(dir: &std::path::Path) -> std::path::PathBuf {
        let records: Vec<bam::Record> = [
            fixture_fragment("aa", true, true, 60, 60),
            fixture_fragment("ar", true, false, 60, 60),
            // An ALT-carrying fragment needs no MAPQ.
            fixture_fragment("ra", false, true, 0, 0),
            fixture_fragment("rr_q0", false, false, 0, 0),
            fixture_fragment("rr_q9", false, false, 9, 0),
            // One read at exactly the floor anchors the fragment.
            fixture_fragment("rr_q10", false, false, 0, 10),
            fixture_fragment("rr_q60", false, false, 60, 60),
        ]
        .into_iter()
        .flatten()
        .collect();
        write_fixture_bam(dir, "rr.bam", records)
    }

    fn fixture_dir(tag: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir()
            .join(format!("cyrius_phase_{}_{}", tag, std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn check_phasing_reports_a_read_error_instead_of_partial_counts() {
        let dir = fixture_dir("readerr");
        // Enough fragments that the BAM spans several BGZF blocks.
        let mut records = Vec::new();
        for i in 0..1500 {
            records.extend(fixture_fragment(&format!("f{i}"), true, true, 60, 60));
        }
        let path = write_fixture_bam(&dir, "readerr.bam", records);
        let mut bytes = std::fs::read(&path).unwrap();
        let size0 = u16::from_le_bytes([bytes[16], bytes[17]]) as usize + 1;
        let size1 = u16::from_le_bytes([bytes[size0 + 16], bytes[size0 + 17]]) as usize + 1;
        let off2 = size0 + size1;
        let size2 = u16::from_le_bytes([bytes[off2 + 16], bytes[off2 + 17]]) as usize + 1;
        assert!(size2 > 28, "too few records for a third data block");
        bytes[off2 + size2 - 8] ^= 0xFF;
        std::fs::write(&path, bytes).unwrap();
        let lookup = fixture_lookup();
        let early = lookup.resolve("g.early").unwrap();
        let late = lookup.resolve("g.late").unwrap();
        let mut reader = bam::IndexedReader::from_path(&path).unwrap();
        let counts = check_phasing(&mut reader, "chr22", &early, &late, false);
        let _ = std::fs::remove_dir_all(&dir);
        assert!(counts.is_err(), "check_phasing gave {counts:?}");
    }

    #[test]
    fn check_phasing_tallies_each_fragment_kind_and_anchors_rr_by_max_mapq() {
        let dir = fixture_dir("rr_tally");
        let bam_path = write_rr_fixture_bam(&dir);
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();
        let lookup = fixture_lookup();
        let rv_early = lookup.resolve("g.early").unwrap();
        let rv_late = lookup.resolve("g.late").unwrap();
        let counts = check_phasing(&mut reader, "chr22", &rv_early, &rv_late, false).unwrap();
        drop(reader);
        let _ = std::fs::remove_dir_all(&dir);

        assert_eq!(
            counts,
            PhaseCounts { aa: 1, ar: 1, ra: 1, rr_all: 4, rr_anchored: 2 },
            "rr_anchored must count only the both-REF fragments with a read at \
             MAPQ >= {RR_MIN_MAPQ} (rr_q10, rr_q60), not rr_q0 or rr_q9"
        );
    }

    #[test]
    fn anchored_rr_fragments_count_as_cis() {
        let dir = fixture_dir("rr_cis");
        let bam_path = write_rr_fixture_bam(&dir);
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();
        let lookup = fixture_lookup();
        let rv_early = lookup.resolve("g.early").unwrap();
        let rv_late = lookup.resolve("g.late").unwrap();
        let counts = check_phasing(&mut reader, "chr22", &rv_early, &rv_late, false).unwrap();
        drop(reader);
        let _ = std::fs::remove_dir_all(&dir);

        assert_eq!(
            counts.cis(),
            3,
            "cis = the AA fragment + the two both-REF fragments with a read at \
             MAPQ >= {RR_MIN_MAPQ}; rr_q0 and rr_q9 must not count; got {counts:?}"
        );
        assert_eq!(counts.trans(), 2, "trans stays AR + RA; got {counts:?}");
    }

    #[test]
    fn rr_fragments_at_mapq_0_do_not_count() {
        let dir = fixture_dir("rr_mapq0");
        let records: Vec<bam::Record> = ["rr1", "rr2", "rr3"]
            .into_iter()
            .flat_map(|q| fixture_fragment(q, false, false, 0, 0))
            .collect();
        let bam_path = write_fixture_bam(&dir, "rr_mapq0.bam", records);
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();
        let lookup = fixture_lookup();
        let rv_early = lookup.resolve("g.early").unwrap();
        let rv_late = lookup.resolve("g.late").unwrap();
        let counts = check_phasing(&mut reader, "chr22", &rv_early, &rv_late, false).unwrap();
        drop(reader);
        let _ = std::fs::remove_dir_all(&dir);

        assert_eq!(counts.rr_all, 3, "{counts:?}");
        assert_eq!(
            (counts.cis(), counts.trans()),
            (0, 0),
            "both-REF fragments with no read at MAPQ >= {RR_MIN_MAPQ} may be \
             CYP2D7 reads placed on CYP2D6; got {counts:?}"
        );
    }

    /// Two options that differ only in phase: `*1/*A` puts the fixture's two
    /// variants on one haplotype (cis), `*B/*C` on two (trans). `g.nolookup`
    /// makes a second differentiating pair that has no lookup entry.
    fn phase_tie() -> (String, Vec<String>, StarCombinations) {
        let dstar: HashMap<String, String> = [
            ("*A", "g.early_g.late_g.nolookup"),
            ("*B", "g.early"),
            ("*C", "g.late_g.nolookup"),
        ]
        .into_iter()
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect();
        // Only `dstar` is read by `disambiguate`.
        let star_combinations = StarCombinations {
            dhap: HashMap::new(),
            dhap2: HashMap::new(),
            dhap3: HashMap::new(),
            dhap3pair: HashMap::new(),
            dhap4pair: HashMap::new(),
            dhap5pair: HashMap::new(),
            dhap6pair: HashMap::new(),
            dhap_exon9_x2: HashMap::new(),
            dhap_exon9_x3: HashMap::new(),
            dhap_exon9_x4: HashMap::new(),
            dhap_dup_exon9: HashMap::new(),
            dstar,
        };
        let raw = vec!["*1_*A".to_string(), "*B_*C".to_string()];
        ("*1/*A;*B/*C".to_string(), raw, star_combinations)
    }

    /// Run `disambiguate` on `phase_tie` over a fixture BAM of `records`.
    fn disambiguate_fixture(tag: &str, records: Vec<bam::Record>) -> PhaseOutcome {
        let dir = fixture_dir(tag);
        let bam_path = write_fixture_bam(&dir, "tie.bam", records);
        let (call, raw, star_combinations) = phase_tie();
        let outcome = disambiguate(
            &call,
            &raw,
            &star_combinations,
            &fixture_lookup(),
            bam_path.to_str().unwrap(),
            "chr22",
            None,
            None,
            false,
        );
        let _ = std::fs::remove_dir_all(&dir);
        outcome
    }

    #[test]
    fn disambiguate_reports_the_pairs_it_counted_and_an_anchored_rr_decides() {
        // 2 AA + 1 anchored RR = 3 cis, 0 trans: a decisive pair for the cis
        // option. Without the RR fragment it would be 2 : 0 and too weak.
        let records: Vec<bam::Record> = [
            fixture_fragment("aa1", true, true, 60, 60),
            fixture_fragment("aa2", true, true, 60, 60),
            fixture_fragment("rr_q60", false, false, 60, 0),
            fixture_fragment("rr_q0", false, false, 0, 0),
        ]
        .into_iter()
        .flatten()
        .collect();
        let outcome = disambiguate_fixture("disamb_rr", records);

        assert_eq!(outcome.resolved.as_deref(), Some("*1/*A"));
        assert_eq!(
            outcome.pairs,
            vec![PairEvidence {
                va: "g.early".to_string(),
                vb: "g.late".to_string(),
                pos_a: FIXTURE_POS_EARLY,
                pos_b: FIXTURE_POS_LATE,
                counts: PhaseCounts { aa: 2, ar: 0, ra: 0, rr_all: 2, rr_anchored: 1 },
            }],
            "only the pair with both lookups is reported; (g.early, g.nolookup) is not"
        );
    }

    #[test]
    fn disambiguate_reports_its_pairs_when_inconclusive() {
        // 2 AA + an RR fragment at MAPQ 0: 2 cis, too weak to decide.
        let records: Vec<bam::Record> = [
            fixture_fragment("aa1", true, true, 60, 60),
            fixture_fragment("aa2", true, true, 60, 60),
            fixture_fragment("rr_q0", false, false, 0, 0),
        ]
        .into_iter()
        .flatten()
        .collect();
        let outcome = disambiguate_fixture("disamb_weak", records);

        assert_eq!(outcome.resolved, None);
        assert_eq!(outcome.pairs.len(), 1, "{:?}", outcome.pairs);
        assert_eq!(
            (outcome.pairs[0].counts.cis(), outcome.pairs[0].counts.trans()),
            (2, 0)
        );
    }

    #[test]
    fn disambiguate_with_nothing_to_phase_reports_no_pairs() {
        // One candidate only: not two options, so nothing is counted.
        let (_, _, star_combinations) = phase_tie();
        let outcome = disambiguate(
            "*1/*A",
            &["*1_*A".to_string()],
            &star_combinations,
            &fixture_lookup(),
            "/nonexistent.bam",
            "chr22",
            None,
            None,
            false,
        );
        assert_eq!(outcome.resolved, None);
        assert!(outcome.pairs.is_empty());
    }

    #[test]
    fn count_pair_matches_check_phasing_and_needs_both_lookups() {
        let dir = fixture_dir("count_pair");
        let bam_path = write_rr_fixture_bam(&dir);
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();
        let lookup = fixture_lookup();
        let pair = count_pair(&mut reader, "chr22", "g.early", "g.late", &lookup).unwrap();
        let rv_early = lookup.resolve("g.early").unwrap();
        let rv_late = lookup.resolve("g.late").unwrap();
        let direct = check_phasing(&mut reader, "chr22", &rv_early, &rv_late, false).unwrap();
        let missing_a = count_pair(&mut reader, "chr22", "g.absent", "g.late", &lookup).unwrap();
        let missing_b = count_pair(&mut reader, "chr22", "g.early", "g.absent", &lookup).unwrap();
        drop(reader);
        let _ = std::fs::remove_dir_all(&dir);

        let pair = pair.expect("both variants are in the lookup");
        assert_eq!(pair.counts, direct);
        assert_eq!((pair.va.as_str(), pair.vb.as_str()), ("g.early", "g.late"));
        assert_eq!((pair.pos_a, pair.pos_b), (FIXTURE_POS_EARLY, FIXTURE_POS_LATE));
        assert_eq!(pair.distance_bp(), 100);
        assert!(missing_a.is_none() && missing_b.is_none());
    }

    #[test]
    fn count_site_counts_reads_as_alt_ref_none() {
        let dir = fixture_dir("count_site");
        let mut records: Vec<bam::Record> = [
            fixture_fragment("a1", true, false, 60, 60),
            fixture_fragment("a2", true, true, 0, 0),
            fixture_fragment("r1", false, true, 60, 60),
        ]
        .into_iter()
        .flatten()
        .collect();
        // Neither ALT 'C' nor REF 'G' at the early site.
        let mut other_seq = vec![b'G'; 50];
        other_seq[(FIXTURE_POS_EARLY - FIXTURE_MATE_POS) as usize] = b'A';
        records.push(fixture_record("n1", FIXTURE_MATE_POS, &other_seq, 60, 0));
        let bam_path = write_fixture_bam(&dir, "site.bam", records);
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();
        let lookup = fixture_lookup();
        let early = count_site(&mut reader, "chr22", "g.early", &lookup);
        let late = count_site(&mut reader, "chr22", "g.late", &lookup);
        let missing = count_site(&mut reader, "chr22", "g.absent", &lookup);
        drop(reader);
        let _ = std::fs::remove_dir_all(&dir);

        // Only the mates overlap the early site, only the reads the late one.
        assert_eq!(early, Some((2, 1, 1)));
        assert_eq!(late, Some((2, 1, 0)));
        assert_eq!(missing, None);
    }

    /// 0-based position of the pileup tests' variant.
    const PILEUP_POS: i64 = 1000;
    /// 0-based start of the pileup tests' reads: 30 bp of flank before the variant.
    const PILEUP_READ_POS: i64 = PILEUP_POS - 30;

    /// The 18-bp insertion `g.42128936-42128937insGGGGCGAAAGGGGCGAAA` with the
    /// patterns the target variant file gives it at 42128927: lower-case
    /// flanks, and a REF that is a prefix of the ALT for its first 16 bases.
    fn repeat_insertion() -> ResolvedVariant {
        ResolvedVariant {
            pos_0based: PILEUP_POS,
            alt_bases: "tGGGGCGAAAGGGGCGAAAGGGGCGAAAGGGGCGt".to_string(),
            ref_bases: "tGGGGCGAAAGGGGCGt".to_string(),
        }
    }

    /// A read at `PILEUP_READ_POS`: 30 bp of 'C' flank, then `at_variant`,
    /// then `tail_len` bases of 'C', aligned with `cigar`.
    fn pileup_record(at_variant: &str, tail_len: usize, cigar: Vec<bam::record::Cigar>) -> bam::Record {
        let mut seq = vec![b'C'; 30];
        seq.extend_from_slice(at_variant.as_bytes());
        seq.extend(std::iter::repeat(b'C').take(tail_len));
        let mut record = bam::Record::new();
        record.set(
            b"r",
            Some(&bam::record::CigarString(cigar)),
            &seq,
            &vec![30u8; seq.len()],
        );
        record.set_pos(PILEUP_READ_POS);
        record
    }

    fn pileup(record: &bam::Record, rv: &ResolvedVariant) -> Option<bool> {
        let seq = record.seq().as_bytes();
        check_allele_pileup(record, &seq, rv)
    }

    #[test]
    fn check_allele_pileup_calls_a_read_truncated_inside_the_repeat_neither() {
        use bam::record::Cigar::Match;
        // The read ends 20 bases into the repeat. Those 20 bases are a prefix
        // of the ALT but hold neither whole pattern (the REF's closing 't' at
        // offset 16 is an 'A' here), so the read cannot tell the alleles apart.
        let truncated = "TGGGGCGAAAGGGGCGAAAG";
        let record = pileup_record(truncated, 0, vec![Match(50)]);
        assert_eq!(
            pileup(&record, &repeat_insertion()),
            None,
            "a read truncated inside the repeat is not ALT just because it is a prefix of ALT"
        );
    }

    #[test]
    fn check_allele_pileup_calls_a_full_alt_read_alt() {
        use bam::record::Cigar::{Ins, Match};
        let alt = "TGGGGCGAAAGGGGCGAAAGGGGCGAAAGGGGCGT";
        // Aligned the way an aligner shows the insertion: 31M 18I 31M.
        let record = pileup_record(alt, 15, vec![Match(31), Ins(18), Match(30)]);
        assert_eq!(pileup(&record, &repeat_insertion()), Some(true));
    }

    #[test]
    fn check_allele_pileup_calls_a_full_ref_read_ref() {
        use bam::record::Cigar::Match;
        let reference = "TGGGGCGAAAGGGGCGT";
        let record = pileup_record(reference, 33, vec![Match(80)]);
        assert_eq!(pileup(&record, &repeat_insertion()), Some(false));
    }

    #[test]
    fn check_allele_pileup_gives_a_read_matching_both_patterns_the_longer_one() {
        use bam::record::Cigar::Match;
        // "TCCAAG..." starts with both "TCC" and "TCCAA".
        let record = pileup_record("TCCAAG", 14, vec![Match(50)]);
        let longer_ref = ResolvedVariant {
            pos_0based: PILEUP_POS,
            alt_bases: "TCC".to_string(),
            ref_bases: "tccaa".to_string(),
        };
        let longer_alt = ResolvedVariant {
            pos_0based: PILEUP_POS,
            alt_bases: "tccaa".to_string(),
            ref_bases: "TCC".to_string(),
        };
        assert_eq!(pileup(&record, &longer_ref), Some(false), "the longer REF must win");
        assert_eq!(pileup(&record, &longer_alt), Some(true), "the longer ALT must win");
    }
}

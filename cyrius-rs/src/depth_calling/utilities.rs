use crate::types::{GmmParameter, Region, RegionDic, SnpLookup};
use indexmap::IndexMap;
use rust_htslib::bam;
use std::collections::HashMap;

/// Parse a BED region file from an embedded string.
/// Returns a RegionDic mapping region_type -> Vec<(Region, gc_content)>.
pub fn parse_region_file(content: &str, genome: &str) -> RegionDic {
    let mut region_dic = RegionDic::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 6 {
            continue;
        }
        let mut nchr = fields[0].to_string();
        if genome == "38" {
            nchr = nchr.replace("chr", "");
        }
        let region_start: i64 = fields[1].parse().unwrap();
        let region_end: i64 = fields[2].parse().unwrap();
        let region_name = fields[3].to_string();
        let region_type = fields[4].to_string();
        let region_gc: f64 = fields[5].parse().unwrap();

        let region: Region = (nchr, region_start, region_end, region_name);
        region_dic
            .entry(region_type)
            .or_default()
            .push((region, region_gc));
    }
    region_dic
}

/// Parse GMM parameter file from an embedded string.
/// Returns variant_id -> parameter_name -> list of values.
pub fn parse_gmm_file(content: &str) -> GmmParameter {
    let mut dpar: GmmParameter = HashMap::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 3 {
            continue;
        }
        let variant_id = fields[0].to_string();
        let param_name = fields[1].to_string();
        let values: Vec<String> = fields[2..]
            .iter()
            .map(|s| {
                // Values are in format "key:value", extract the value part
                if let Some(pos) = s.find(':') {
                    s[pos + 1..].to_string()
                } else {
                    s.to_string()
                }
            })
            .collect();
        dpar.entry(variant_id)
            .or_default()
            .entry(param_name)
            .or_insert(values);
    }
    dpar
}

/// Position `reader` on `[start, end)` of reference `tid`.
///
/// A failed fetch (a broken or stale index, a CRAM reference that cannot be
/// resolved, an inverted interval) means the reads were never looked at, which
/// is not the same as there being none. Return an Err naming `what` so the
/// caller can turn it into a no-call instead of counting zero reads
/// (escalation E1 of review batch 4).
pub fn fetch_region(
    reader: &mut bam::IndexedReader,
    tid: u32,
    start: i64,
    end: i64,
    what: &str,
) -> Result<(), String> {
    reader
        .fetch(bam::FetchDefinition::Region(tid as i32, start, end))
        .map_err(|e| format!("BAM fetch failed for {what}: {e}"))
}

/// Open a BAM or CRAM alignment file using rust-htslib.
pub fn open_alignment_file(
    alignment_file: &str,
    reference_fasta: Option<&str>,
) -> Result<bam::IndexedReader, Box<dyn std::error::Error>> {
    open_alignment_file_with_index(alignment_file, reference_fasta, None)
}

/// Open a BAM or CRAM alignment file with optional custom index path.
pub fn open_alignment_file_with_index(
    alignment_file: &str,
    reference_fasta: Option<&str>,
    index_filename: Option<&str>,
) -> Result<bam::IndexedReader, Box<dyn std::error::Error>> {
    let mut reader = if alignment_file.to_lowercase().ends_with("cram") {
        let mut r = match index_filename {
            Some(idx) => bam::IndexedReader::from_path_and_index(alignment_file, idx)?,
            None => bam::IndexedReader::from_path(alignment_file)?,
        };
        if let Some(ref_path) = reference_fasta {
            r.set_reference(ref_path)?;
        }
        r
    } else {
        match index_filename {
            Some(idx) => bam::IndexedReader::from_path_and_index(alignment_file, idx)?,
            None => bam::IndexedReader::from_path(alignment_file)?,
        }
    };

    use rust_htslib::bam::Read as HtslibRead;
    reader.set_threads(1)?;
    Ok(reader)
}

/// How many leading contigs the genome autodetect inspects for a "chr" prefix.
/// Sampling only the head of the header is cheap but fallible: a header whose
/// first few contigs are bare while the rest carry "chr" resolves to the "38"
/// style and then misses every chr-prefixed lookup.
const AUTODETECT_CONTIG_SAMPLE: usize = 5;

/// How many header contigs `missing_contig_error` names before summarising the
/// rest as a count.
const CONTIG_PREVIEW: usize = 5;

/// All contig names in an alignment file header, in header order.
pub fn header_contig_names(
    input_file: &str,
    index_name: Option<&str>,
) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    use rust_htslib::bam::Read as BamRead;
    let reader = open_alignment_file_with_index(input_file, None, index_name)?;
    let header = reader.header().clone();
    Ok((0..header.target_count())
        .map(|contig_id| String::from_utf8_lossy(header.tid2name(contig_id)).to_string())
        .collect())
}

/// Whether the header uses "chr"-prefixed contig names, judged the way the
/// genome autodetect has always judged it: "chr" occurring anywhere in any of
/// the first `AUTODETECT_CONTIG_SAMPLE` contig names.
pub fn header_has_chr_prefix(header_contigs: &[String]) -> bool {
    header_contigs
        .iter()
        .take(AUTODETECT_CONTIG_SAMPLE)
        .any(|name| name.contains("chr"))
}

/// The distinct CYP2D6 contig names a run will look up, in the order to report
/// them.
///
/// `Vec::dedup` alone would not do: it drops only *consecutive* duplicates, so
/// a third name equal to the first would survive as `[a, b, a]` and log the
/// same error twice. Comparing against everything already collected keeps the
/// order the caller asked for and stays correct for any number of names.
pub fn contig_check_targets(snp_nchr: &str, var_nchr: &str) -> Vec<String> {
    let mut targets: Vec<String> = Vec::new();
    for nchr in [snp_nchr, var_nchr] {
        if !targets.iter().any(|seen| seen == nchr) {
            targets.push(nchr.to_string());
        }
    }
    targets
}

/// Build the message to log when the CYP2D6 contig for the resolved genome
/// build is absent from the alignment header; `None` when the contig is
/// present.
///
/// Names must match exactly, the way every contig lookup in the pipeline
/// matches them: with a resolved genome of `38` the target is `22`, which a
/// `chr22`-named header does not provide.
///
/// The message names the sample and the file because batch runs fan many
/// samples into a single log stream, and it names the *resolved* genome rather
/// than blaming `-g/--genome`, which may never have been passed: a header
/// whose first contigs are bare and whose rest carry "chr" makes the autodetect
/// resolve `38` on its own. For the same reason the no-reads claim is scoped to
/// the missing contig -- lookups on contigs the header does hold, including
/// most of the depth-normalization control regions, still return reads.
pub fn missing_contig_error(
    sample_id: &str,
    input_file: &str,
    genome: &str,
    target: &str,
    header_contigs: &[String],
) -> Option<String> {
    if header_contigs.iter().any(|contig| contig == target) {
        return None;
    }
    let listing = if header_contigs.is_empty() {
        "0 contig(s); the header names none".to_string()
    } else {
        let preview: Vec<&str> = header_contigs
            .iter()
            .take(CONTIG_PREVIEW)
            .map(|contig| contig.as_str())
            .collect();
        let hidden = header_contigs.len() - preview.len();
        let more = if hidden > 0 {
            format!(" (+{hidden} more)")
        } else {
            String::new()
        };
        format!(
            "{} contig(s), first: [{}]{}",
            header_contigs.len(),
            preview.join(", "),
            more
        )
    };
    Some(format!(
        "Sample {} ({}): CYP2D6 contig '{}' for resolved genome '{}' is not in the alignment \
         header ({}). Every region and SNP lookup on that contig will return no reads, so the \
         result will be an uninformative no-call. The resolved genome comes from -g/--genome, \
         or, when that is left at its default, from autodetecting a 'chr' prefix in the first \
         {} header contigs -- check that it matches the reference the reads were aligned to.",
        sample_id, input_file, target, genome, listing, AUTODETECT_CONTIG_SAMPLE
    ))
}

/// Parse SNP position file from an embedded string.
/// Matches Python get_snp_position().
pub fn get_snp_position(content: &str, genome: &str, group: Option<&str>) -> SnpLookup {
    let mut dsnp1 = IndexMap::new();
    let mut dsnp2 = IndexMap::new();
    let mut dindex: HashMap<String, usize> = HashMap::new();
    let mut counter: i32 = -1;
    let mut last_nchr = String::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with('\n') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 5 {
            continue;
        }

        let line_group = fields.last().unwrap();
        if group.is_some() && Some(*line_group) != group {
            continue;
        }

        counter += 1;
        let idx = counter as usize;

        let reg1_name = format!("{}_{}", fields[1], idx);
        let reg2_name = format!("{}_{}", fields[3], idx);
        let reg1_base = fields[2].to_uppercase();
        let reg2_base = fields[4].to_uppercase();

        if *line_group != "-" {
            dsnp1
                .entry(reg1_name.clone())
                .or_insert_with(|| format!("{}_{}", reg1_base, reg2_base));
            dsnp2
                .entry(reg2_name.clone())
                .or_insert_with(|| format!("{}_{}", reg1_base, reg2_base));
        } else {
            dsnp1
                .entry(reg1_name.clone())
                .or_insert_with(|| format!("{}_{}", reg1_base, reverse_complement(&reg2_base)));
            dsnp2
                .entry(reg2_name.clone())
                .or_insert_with(|| format!("{}_{}", reverse_complement(&reg1_base), reg2_base));
        }
        dindex.insert(reg1_name, idx);
        dindex.insert(reg2_name, idx);
        last_nchr = fields[0].to_string();
    }

    if counter == -1 {
        panic!("No valid SNP positions found in file");
    }

    let nchr = if genome == "38" {
        last_nchr.replace("chr", "")
    } else {
        last_nchr
    };

    SnpLookup {
        dsnp1,
        dsnp2,
        nchr,
        dindex,
    }
}

/// Get variant names from a variant file (embedded string).
pub fn get_var_names(content: &str) -> Vec<String> {
    let mut var_list = Vec::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if let Some(var_name) = fields.last() {
            var_list.push(var_name.to_string());
        }
    }
    var_list
}

const COMPLEMENT: &[(char, char)] = &[
    ('A', 'T'),
    ('T', 'A'),
    ('C', 'G'),
    ('G', 'C'),
    ('N', 'N'),
];

fn complement_base(base: char) -> char {
    for &(from, to) in COMPLEMENT {
        if from == base {
            return to;
        }
    }
    base
}

fn reverse_complement(sequence: &str) -> String {
    sequence.chars().rev().map(complement_base).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Review finding #12: the g.42128445delCTGCC row in
    /// data/CYP2D6_target_variant_38.txt (line 64) was off-by-one on its
    /// D6 position and had its ALT/REF alleles swapped, so the site was
    /// never counted. Pin the corrected row as actually parsed by
    /// get_snp_position() (the loader `main.rs` uses for this file),
    /// through the real embedded data, not a hand-rolled re-parse.
    #[test]
    fn target_variant_delctgcc_site_is_parsed_with_corrected_position_and_alleles() {
        let var_db = get_snp_position(crate::data::TARGET_VARIANT, "38", None);

        // This is the 63rd non-comment data row in the file (0-based index 62).
        let d6_key = "42128443_62";
        let d7_key = "42142135_62";

        assert_eq!(
            var_db.dsnp1.get(d6_key).map(|s| s.as_str()),
            Some("GCT_GCCTGCCT"),
            "corrected D6-side entry for g.42128445delCTGCC not found at {d6_key}; dsnp1 did not \
             parse the fixed position/alleles"
        );
        assert_eq!(
            var_db.dsnp2.get(d7_key).map(|s| s.as_str()),
            Some("GCT_GCCTGCCT"),
            "corrected D7-side entry for g.42128445delCTGCC not found at {d7_key}"
        );

        // The old, mis-specified key must be gone.
        assert!(
            !var_db.dsnp1.contains_key("42128442_62"),
            "stale off-by-one key 42128442_62 should no longer be present"
        );
    }

    // ---- Review finding #63: missing-contig diagnostics ------------------
    // These four moved here verbatim from `main.rs`'s own `mod tests` (the bin
    // target, which `cargo test --lib` skips) so the lib suite covers them,
    // then grew the assertions the review asked for.

    fn contigs(names: &[&str]) -> Vec<String> {
        names.iter().map(|n| n.to_string()).collect()
    }

    #[test]
    fn present_contig_is_not_an_error() {
        assert_eq!(
            missing_contig_error(
                "HG00111",
                "/data/HG00111.bam",
                "chr38",
                "chr22",
                &contigs(&["chr1", "chr2", "chr22"])
            ),
            None
        );
    }

    #[test]
    fn absent_contig_message_names_the_target_and_the_header_contigs() {
        let header = contigs(&["chr1", "chr2", "chr3", "chr4", "chr5", "chr22"]);
        let message = missing_contig_error("HG00111", "/data/HG00111.bam", "38", "22", &header)
            .expect("a target missing from the header must produce an error message");
        // The target has to be named, or the message does not say what to fix.
        assert!(
            message.contains("'22'"),
            "message must name the missing target: {message}"
        );
        // Listing what the header does hold is what makes the message useful.
        for name in ["chr1", "chr2", "chr3", "chr4", "chr5"] {
            assert!(
                message.contains(name),
                "message must list header contig {name}: {message}"
            );
        }
        assert!(
            message.contains("6 contig"),
            "message must report the header contig count: {message}"
        );
    }

    /// `compare_to_baseline.sh` fans 62 samples through `xargs -P 32` into a
    /// single stream, so an error line naming neither the sample nor the file
    /// cannot be acted on. Its sibling diagnostics name the sample; so must
    /// this one.
    #[test]
    fn absent_contig_message_attributes_itself_to_a_sample_and_a_file() {
        let message = missing_contig_error(
            "NA18565",
            "/bench/NA18565_pgx_regions.bam",
            "38",
            "22",
            &contigs(&["chr22"]),
        )
        .expect("a target missing from the header must produce an error message");
        assert!(
            message.contains("NA18565"),
            "message must name the sample it belongs to: {message}"
        );
        assert!(
            message.contains("/bench/NA18565_pgx_regions.bam"),
            "message must name the alignment file: {message}"
        );
    }

    /// The message used to claim that *every* region and SNP lookup returns no
    /// reads, and to blame `-g/--genome`. Both are false in a reachable case:
    /// a header whose first five contigs are bare and whose rest carry "chr"
    /// makes autodetect resolve "38" with no `-g` given at all, and the
    /// chr1-chr22 control regions still return reads (the very next log line
    /// reports some control regions usable, not none). Scope the claim to the
    /// missing contig and name the *resolved* genome.
    #[test]
    fn absent_contig_message_scopes_the_claim_and_names_the_resolved_genome() {
        let mut header = contigs(&["1", "2", "3", "4", "5"]);
        header.extend(contigs(&["chr6", "chr7", "chr12", "chr22"]));
        let message = missing_contig_error("HG00111", "/data/HG00111.bam", "38", "22", &header)
            .expect("'22' is absent from a header whose chromosome 22 is named 'chr22'");
        assert!(
            message.contains("on that contig"),
            "the no-reads claim must be scoped to the missing contig, because lookups \
             on contigs the header does hold still return reads: {message}"
        );
        assert!(
            !message.contains("Every region and SNP lookup will return no reads"),
            "the unscoped claim is falsifiable and must be gone: {message}"
        );
        assert!(
            message.contains("resolved genome '38'"),
            "the message must name the resolved genome, so it explains itself when \
             autodetect rather than -g chose it: {message}"
        );
        assert!(
            message.contains("autodetect"),
            "the message must point at autodetect, which is to blame when -g was \
             never passed: {message}"
        );
    }

    #[test]
    fn a_chr_prefixed_header_does_not_satisfy_the_stripped_target() {
        // `-g 38` strips the "chr" prefix, so "22" must not match "chr22":
        // this is exactly the wrong-build case the error exists to report.
        assert!(missing_contig_error("S", "/x.bam", "38", "22", &contigs(&["chr22"])).is_some());
        // Nor may a longer contig name satisfy a shorter target by substring.
        assert!(missing_contig_error(
            "S",
            "/x.bam",
            "chr38",
            "chr22",
            &contigs(&["chr22_KI270879v1"])
        )
        .is_some());
    }

    #[test]
    fn an_empty_header_is_an_error_and_lists_nothing() {
        let message = missing_contig_error("S", "/x.bam", "chr38", "chr22", &[])
            .expect("an empty header cannot hold the target contig");
        assert!(message.contains("'chr22'"), "{message}");
        assert!(message.contains("0 contig"), "{message}");
        // The test name promises the message lists nothing, so check it: the
        // old rendering still printed an empty "first: []".
        assert!(
            !message.contains("first:"),
            "an empty header has no first contigs to show: {message}"
        );
        assert!(
            !message.contains('['),
            "nothing may be listed for an empty header: {message}"
        );
    }

    /// The preview used to end in a bare ", ..." -- rendered even when exactly
    /// one contig was hidden, and silent about how many.
    #[test]
    fn a_long_header_reports_how_many_contigs_it_hides() {
        let six = contigs(&["c1", "c2", "c3", "c4", "c5", "c6"]);
        let message = missing_contig_error("S", "/x.bam", "38", "22", &six)
            .expect("'22' is absent from this header");
        assert!(
            message.contains("(+1 more)"),
            "one hidden contig must be counted, not elided: {message}"
        );
        assert!(
            !message.contains("c6"),
            "the sixth contig is past the preview: {message}"
        );

        let five = contigs(&["c1", "c2", "c3", "c4", "c5"]);
        let message = missing_contig_error("S", "/x.bam", "38", "22", &five)
            .expect("'22' is absent from this header");
        assert!(
            !message.contains("more)"),
            "nothing is hidden, so nothing is summarised: {message}"
        );
    }

    /// `Vec::dedup` removes only *consecutive* duplicates, so the call site's
    /// `targets.dedup()` was correct only because the vec held exactly two
    /// always-equal names; a third equal to the first would have survived as
    /// `[a, b, a]` and logged the same error twice.
    #[test]
    fn contig_check_targets_collapses_two_equal_names_to_one() {
        assert_eq!(contig_check_targets("22", "22"), vec!["22".to_string()]);
        assert_eq!(
            contig_check_targets("chr22", "chr22"),
            vec!["chr22".to_string()]
        );
    }

    /// Nothing tied the *choice* of checked contigs to a test. Both databases
    /// must be reported when their contigs disagree, in the order checked.
    #[test]
    fn contig_check_targets_keeps_both_names_when_they_differ() {
        assert_eq!(
            contig_check_targets("22", "chr22"),
            vec!["22".to_string(), "chr22".to_string()]
        );
    }

    /// Why one target is all the check needs today: both `nchr` values come
    /// from the last data row of their own file, column 1 is `chr22` in both,
    /// and the same `genome` is applied to both. Assert that against the real
    /// embedded data rather than assuming it.
    #[test]
    fn the_real_snp_and_variant_databases_name_the_same_contig() {
        for genome in ["38", "chr38"] {
            let snp_db = get_snp_position(crate::data::SNP_FILE, genome, None);
            let var_db = get_snp_position(crate::data::TARGET_VARIANT, genome, None);
            assert_eq!(
                contig_check_targets(&snp_db.nchr, &var_db.nchr),
                vec![snp_db.nchr.clone()],
                "genome {genome}: snp_db.nchr={} but var_db.nchr={}",
                snp_db.nchr,
                var_db.nchr
            );
        }
    }

    /// One header read now feeds both the genome autodetect and the
    /// missing-contig check, so the autodetect rule it used to own -- "chr"
    /// anywhere in the first five contig names -- has to be pinned exactly as
    /// it was.
    #[test]
    fn header_has_chr_prefix_judges_only_the_first_five_contigs() {
        assert!(header_has_chr_prefix(&contigs(&["chr1", "chr2"])));
        assert!(header_has_chr_prefix(&contigs(&[
            "1", "2", "3", "4", "chr5"
        ])));
        // The sixth contig is past the sample: this is the hole that lets a
        // header resolve to "38" and then miss every chr-prefixed lookup.
        assert!(!header_has_chr_prefix(&contigs(&[
            "1", "2", "3", "4", "5", "chr6", "chr22"
        ])));
        assert!(!header_has_chr_prefix(&contigs(&["1", "2", "22"])));
        assert!(!header_has_chr_prefix(&[]));
    }
}

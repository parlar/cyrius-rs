use clap::Parser;
use indexmap::IndexMap;
use log;
use std::collections::HashMap;
use std::path::Path;

use cyrius_rs::caller::call_variants::{
    call_cn_snp, call_cn_var, call_cn_var_homo, call_exon9gc, call_var42126938,
    call_var42127526_var42127556, call_var42127803hap, call_var42130655ins_a,
    get_allele_counts_var42128936, get_called_variants, reportable_cnvtag, update_var42128936,
    NOISY_VAR,
};
use cyrius_rs::caller::cnv_hybrid::get_cnvtag;
use cyrius_rs::caller::construct_star_table::{combination_table_for, get_hap_table};
use cyrius_rs::caller::match_star_allele::match_star;
use cyrius_rs::caller::strand_bias_all;
use cyrius_rs::data;
use cyrius_rs::types::FeatureFlags;
use cyrius_rs::depth_calling::bin_count::{get_normed_depth, get_normed_depth_from_count};
use cyrius_rs::depth_calling::gmm::Gmm;
use cyrius_rs::depth_calling::snp_count::{
    get_fraction, get_supporting_reads, get_supporting_reads_single_region,
};
use cyrius_rs::depth_calling::utilities::{
    contig_check_targets, get_snp_position, get_var_names, header_contig_names,
    header_has_chr_prefix, missing_contig_error, open_alignment_file_with_index, parse_gmm_file,
    parse_region_file,
};
use cyrius_rs::phenotype;
use cyrius_rs::types::*;

const TOOL_NAME: &str = "BCyrius";
const TOOL_VERSION: &str = "1.0.5";
const MAD_THRESHOLD: f64 = 0.11;
const EXON9_SITE1: usize = 7;
const EXON9_SITE2: usize = 8;
const HIGH_CN_DEPTH_THRESHOLD: f64 = 7.5;
const HAPLOTYPE_VAR: &[&str] = &[
    "g.42126938C>T",
    "g.42127803C>T",
    "g.42127526C>T_g.42127556T>C",
    "g.42130655-42130656insA",
];

#[derive(Parser, Debug)]
#[command(name = TOOL_NAME, version = TOOL_VERSION, about = "CYP2D6 genotyping tool")]
struct Cli {
    #[arg(short = 'i', long = "input", help = "Input BAM/CRAM file")]
    input: String,

    #[arg(
        short = 'g',
        long = "genome",
        help = "Reference genome style",
        default_value = "autodetect",
        value_parser = ["autodetect", "38", "chr38"]
    )]
    genome: String,

    #[arg(short = 'o', long = "outDir", help = "Output directory")]
    out_dir: String,

    #[arg(long = "id", help = "Sample ID (output file name)")]
    id: Option<String>,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        default_value = "1"
    )]
    threads: usize,

    #[arg(long = "countFilePath", help = "Path to count files")]
    count_file_path: Option<String>,

    #[arg(short = 'r', long = "reference", help = "Path to reference FASTA for CRAM")]
    reference: Option<String>,

    #[arg(long = "population-info", help = "Export population frequencies")]
    population_info: bool,

    #[arg(long = "haplotype-info", help = "Output haplotype information")]
    haplotype_info: bool,

    #[arg(long = "print", help = "Print results on screen")]
    print: bool,

    #[arg(long = "strand-bias-all", help = "Apply strand bias filtering to all variant sites")]
    strand_bias_all: bool,

    #[arg(long = "fuzzy-match", help = "Use fuzzy matching when exact star allele lookup fails")]
    fuzzy_match: bool,

    #[arg(long = "quality-aware", help = "[EXPERIMENTAL/INACTIVE] base-quality-aware likelihoods — module present but not yet wired into the calling path; output is unchanged")]
    quality_aware: bool,

    #[arg(long = "phase-disambiguate", help = "Use physical phasing to disambiguate multiple-match diplotypes")]
    phase_disambiguate: bool,

    #[arg(long = "phase-readpair", help = "Use read-pair (mate-pair) phasing for longer-range disambiguation")]
    phase_readpair: bool,

    #[arg(long = "changepoint-hybrid", help = "Use paralog ratio changepoint detection for hybrid identification")]
    changepoint_hybrid: bool,

    #[arg(long = "het-check", help = "Use hemizygosity detection (zero het sites = deletion allele)")]
    het_check: bool,

    #[arg(long = "spacer-cn-check", help = "Flag calls as uncertain when spacer CN is inconsistent with CNV group")]
    spacer_cn_check: bool,

    #[arg(long = "consistency-check", help = "Run post-hoc consistency checks (conversion map, mismatch rate, allele balance)")]
    consistency_check: bool,

    #[arg(long = "hmm-cnv", help = "Use HMM-based CNV segmentation as fallback when consensus fails")]
    hmm_cnv: bool,

    #[arg(long = "diplotype-caller", help = "Run likelihood-based diplotype caller (D6+D7 mixture model)")]
    diplotype_caller: bool,

    #[arg(long = "clip-evidence", help = "Detect soft-clip clusters as structural breakpoint evidence")]
    clip_evidence: bool,

    #[arg(long = "d7-depth", help = "Compute CYP2D7 sub-region depth profile for hybrid confirmation")]
    d7_depth: bool,

    #[arg(long = "kmer-validation", help = "Use alignment-free k-mer paralog ratio validation")]
    kmer_validation: bool,
}

/// Write a file, exiting with a clear error message on failure.
fn write_file_or_exit(path: &str, content: &str) {
    if let Err(e) = std::fs::write(path, content) {
        eprintln!("Cannot write output file {}: {}", path, e);
        std::process::exit(1);
    }
}

/// Prepare all resources from embedded data files.
///
/// Also returns the alignment header's contig names, or `None` when the header
/// could not be read. The header is read exactly once here, because the genome
/// autodetect and the caller's missing-contig check want the same list; before
/// this the file was opened and its header walked twice, one statement apart.
fn prepare_resource(
    genome_arg: &str,
    input_file: &str,
    index_name: Option<&str>,
) -> (ResourceInfo, Option<Vec<String>>) {
    let header_contigs = match header_contig_names(input_file, index_name) {
        Ok(contigs) => Some(contigs),
        Err(e) => {
            // Preserve the autodetect fallback exactly as it was: the "38"
            // style (no "chr" prefix), announced with the same message at the
            // same level, and only when autodetect is what cannot be done.
            // The main pipeline will surface a proper no-call.
            if genome_arg == "autodetect" {
                log::warn!(
                    "Could not open {} for chromosome autodetect ({}); assuming GRCh38 without 'chr' prefix",
                    input_file, e
                );
            }
            log::warn!(
                "Could not read the header of {} to check for the CYP2D6 contig: {}",
                input_file,
                e
            );
            None
        }
    };

    let genome = if genome_arg == "autodetect" {
        let chr_prefixed = match header_contigs.as_deref() {
            Some(contigs) => header_has_chr_prefix(contigs),
            // No header to judge by; keep the historical default.
            None => false,
        };
        if chr_prefixed {
            "chr38".to_string()
        } else {
            "38".to_string()
        }
    } else {
        genome_arg.to_string()
    };

    let gmm_parameter = parse_gmm_file(data::GMM_PARAMS);
    let region_dic = parse_region_file(data::REGION_BED, &genome);
    let snp_db = get_snp_position(data::SNP_FILE, &genome, None);
    let var_db = get_snp_position(data::TARGET_VARIANT, &genome, None);
    let var_homo_db = get_snp_position(data::TARGET_VARIANT_HOMO, &genome, None);

    let mut haplotype_db = HashMap::new();
    for &variant in HAPLOTYPE_VAR {
        haplotype_db.insert(
            variant.to_string(),
            get_snp_position(data::HAPLOTYPE_FILE, &genome, Some(variant)),
        );
    }

    let mut var_list = get_var_names(data::TARGET_VARIANT);
    var_list.extend(get_var_names(data::TARGET_VARIANT_HOMO));

    // No combination table: `d6_star_caller` builds the one its CNV tag needs (#41).
    let star_combinations = get_hap_table(data::STAR_TABLE, &[]);

    (
        ResourceInfo {
            genome,
            gmm_parameter,
            region_dic,
            snp_db,
            var_db,
            var_homo_db,
            haplotype_db,
            var_list,
            star_combinations,
        },
        header_contigs,
    )
}

/// Main CYP2D6 star allele calling pipeline.
fn d6_star_caller(
    bam_path: &str,
    params: &ResourceInfo,
    threads: usize,
    count_file: Option<&str>,
    reference_fasta: Option<&str>,
    index_name: Option<&str>,
    features: &FeatureFlags,
) -> D6Call {
    // 1. Read counting and normalization
    let normalized_depth = if let Some(cf) = count_file {
        let mut reader = match open_alignment_file_with_index(bam_path, reference_fasta, index_name) {
            Ok(r) => r,
            Err(e) => {
                log::error!("Cannot open alignment file {}: {}", bam_path, e);
                return D6Call::no_call(0.0, 0.0);
            }
        };
        cyrius_rs::depth_calling::bin_count::get_read_length(&mut reader)
            .and_then(|read_length| get_normed_depth_from_count(cf, &params.region_dic, read_length))
    } else {
        get_normed_depth(bam_path, &params.region_dic, threads, reference_fasta, index_name)
    };
    let normalized_depth = match normalized_depth {
        Ok(d) => d,
        Err(e) => {
            log::error!("Depth normalization failed: {}", e);
            return D6Call::no_call(0.0, 0.0);
        }
    };

    // No-call after normalization
    if normalized_depth.normalized.get("d67").and_then(|v| *v).is_none() {
        return D6Call {
                rep7_mapq0_fraction: normalized_depth.rep7_mapq0_fraction,
                ..D6Call::no_call(normalized_depth.mad, normalized_depth.mediandepth)
            };
    }

    // 2. GMM and CN call
    let d67_depth = normalized_depth.normalized["d67"].unwrap();
    let spacer_depth = normalized_depth.normalized.get("spacer").and_then(|v| *v).unwrap_or(0.0);

    let mut gmm_d67 = Gmm::new();
    gmm_d67.set_gmm_par(&params.gmm_parameter, "d67");
    let gcall_d67 = gmm_d67.gmm_call(d67_depth);

    let mut gmm_spacer = Gmm::new();
    gmm_spacer.set_gmm_par(&params.gmm_parameter, "spacer");
    let gcall_spacer = gmm_spacer.gmm_call(spacer_depth);

    let mut high_cn_low_confidence = false;
    let raw_cn_call = if gcall_d67.cn.is_none() && gcall_d67.depth_value > HIGH_CN_DEPTH_THRESHOLD {
        high_cn_low_confidence = true;
        RawCnCall {
            d67_cn: Some(gcall_d67.depth_value.round() as u32),
            d67_depth: gcall_d67.depth_value,
            spacer_cn: gcall_spacer.cn,
            spacer_depth: gcall_spacer.depth_value,
        }
    } else {
        RawCnCall {
            d67_cn: gcall_d67.cn,
            d67_depth: gcall_d67.depth_value,
            spacer_cn: gcall_spacer.cn,
            spacer_depth: gcall_spacer.depth_value,
        }
    };

    // 3. Get allele counts (parallel when threads > 1). Each worker opens its
    // own reader. A failed open or a failed BAM fetch is a no-call: reads that
    // could not be looked at are not the same as zero reads.
    let open = || {
        open_alignment_file_with_index(bam_path, reference_fasta, index_name)
            .map_err(|e| format!("cannot open alignment file {}: {}", bam_path, e))
    };
    let allele_counts: Result<_, String> = if threads > 1 {
        std::thread::scope(|s| {
            let h_snp = s.spawn(|| get_supporting_reads(&mut open()?, &params.snp_db));
            let h_var = s.spawn(|| get_supporting_reads_single_region(&mut open()?, &params.var_db, None));
            let h_ins = s.spawn(|| get_allele_counts_var42128936(&mut open()?, &params.genome));
            let h_homo = s.spawn(|| get_supporting_reads(&mut open()?, &params.var_homo_db));

            let (snp_d6, snp_d7, snp_crossing) = h_snp.join().unwrap()?;
            let (var_alt, var_ref, var_alt_forward, var_alt_reverse) = h_var.join().unwrap()?;
            let (ref_read, long_ins_read, short_ins_read) = h_ins.join().unwrap()?;
            let (var_homo_alt, var_homo_ref, _) = h_homo.join().unwrap()?;
            Ok((snp_d6, snp_d7, snp_crossing, var_alt, var_ref, var_alt_forward, var_alt_reverse,
             ref_read, long_ins_read, short_ins_read, var_homo_alt, var_homo_ref))
        })
    } else {
        (|| {
            let mut reader = open()?;
            let (snp_d6, snp_d7, snp_crossing) = get_supporting_reads(&mut reader, &params.snp_db)?;
            let (var_alt, var_ref, var_alt_forward, var_alt_reverse) =
                get_supporting_reads_single_region(&mut reader, &params.var_db, None)?;
            let (ref_read, long_ins_read, short_ins_read) =
                get_allele_counts_var42128936(&mut reader, &params.genome)?;
            let (var_homo_alt, var_homo_ref, _) = get_supporting_reads(&mut reader, &params.var_homo_db)?;
            Ok((snp_d6, snp_d7, snp_crossing, var_alt, var_ref, var_alt_forward, var_alt_reverse,
             ref_read, long_ins_read, short_ins_read, var_homo_alt, var_homo_ref))
        })()
    };
    let (snp_d6, snp_d7, snp_crossing, mut var_alt, mut var_ref, var_alt_forward, var_alt_reverse,
     ref_read, long_ins_read, short_ins_read, var_homo_alt, var_homo_ref) = match allele_counts {
        Ok(v) => v,
        Err(e) => {
            log::error!("Cannot read allele counts: {}", e);
            return D6Call {
                rep7_mapq0_fraction: normalized_depth.rep7_mapq0_fraction,
                ..D6Call::no_call(normalized_depth.mad, normalized_depth.mediandepth)
            };
        }
    };

    update_var42128936(
        &params.var_list,
        &mut var_alt,
        &mut var_ref,
        ref_read,
        long_ins_read,
        short_ins_read,
    );

    // Build raw count dict
    let mut raw_count = IndexMap::new();
    let non_homo_count = var_alt.len();
    for (i, var_name) in params.var_list.iter().enumerate() {
        if i < non_homo_count {
            if NOISY_VAR.contains(&var_name.as_str()) {
                raw_count.insert(
                    var_name.clone(),
                    format!(
                        "{}({}:{}),{}",
                        var_alt[i], var_alt_forward[i], var_alt_reverse[i], var_ref[i]
                    ),
                );
            } else {
                raw_count.insert(var_name.clone(), format!("{},{}", var_alt[i], var_ref[i]));
            }
        } else {
            let homo_idx = i - non_homo_count;
            raw_count.insert(
                var_name.clone(),
                format!("{},{}", var_homo_alt[homo_idx], var_homo_ref[homo_idx]),
            );
        }
    }

    // No-call due to total CN
    if raw_cn_call.d67_cn.is_none() {
        return D6Call {
            coverage_mad: normalized_depth.mad,
            median_depth: normalized_depth.mediandepth,
            total_cn: None, spacer_cn: raw_cn_call.spacer_cn,
            total_cn_raw: Some(raw_cn_call.d67_depth), spacer_cn_raw: Some(raw_cn_call.spacer_depth),
            variants_called: None, cnv_group: None, genotype: None, filter: None,
            raw_star_allele: None, call_info: None, exon9_cn: None, cnv_consensus: None,
            d67_snp_call: None, d67_snp_raw: None, variant_raw_count: Some(raw_count),
            kmer_validation: None, resolver_posterior: None,
            rep7_mapq0_fraction: normalized_depth.rep7_mapq0_fraction,
        };
    }

    let total_cn = raw_cn_call.d67_cn.unwrap();

    // 4. Call CNV and hybrids
    let d6_fraction = get_fraction(&snp_d6, &snp_d7);
    let raw_d6_cn: Vec<f64> = d6_fraction
        .iter()
        .map(|&a| cyrius_rs::stats::python_round3(total_cn as f64 * a))
        .collect();
    let cn_call_snp = call_cn_snp(total_cn, &snp_d6, &snp_d7, 0.6);

    let exon9gc_call_stringent = if snp_d6.len() > EXON9_SITE2 && snp_d7.len() > EXON9_SITE2 {
        call_exon9gc(
            &snp_d6[EXON9_SITE1..=EXON9_SITE2],
            &snp_d7[EXON9_SITE1..=EXON9_SITE2],
            Some(total_cn),
        )
    } else {
        None
    };

    let (cnvtag, consensus) = get_cnvtag(
        total_cn,
        &raw_d6_cn,
        &cn_call_snp,
        exon9gc_call_stringent,
        raw_cn_call.spacer_cn,
    );

    // Format consensus string
    let consensus_str = format!(
        "{},{},{},{},{}",
        consensus.rep,
        consensus.exon9_and_downstream.map_or("None".to_string(), |v| v.to_string()),
        consensus.exon9_to_intron4.map_or("None".to_string(), |v| v.to_string()),
        consensus.intron4_to_intron1.map_or("None".to_string(), |v| v.to_string()),
        consensus.intron1_upstream.map_or("None".to_string(), |v| v.to_string()),
    );
    let cn_call_snp_str = cn_call_snp
        .iter()
        .map(|c| c.map_or("None".to_string(), |v| v.to_string()))
        .collect::<Vec<_>>()
        .join(",");
    let raw_d6_cn_str = raw_d6_cn
        .iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(",");

    // Feature: changepoint_hybrid — use CBS on D6/D7 ratio to detect hybrids
    let cnvtag = if features.changepoint_hybrid && cnvtag.is_none() {
        let cp_result = cyrius_rs::caller::changepoint::detect_changepoints(
            &snp_d6, &snp_d7, 5, 3.0, 5,
        );
        if cp_result.is_hybrid {
            log::info!(
                "changepoint_hybrid: detected hybrid pattern ({} segments, overall_ratio={:.3})",
                cp_result.segments.len(),
                cp_result.overall_ratio,
            );
            for (i, seg) in cp_result.segments.iter().enumerate() {
                log::info!(
                    "  segment {}: idx {}..{}, mean_ratio={:.3}, n_valid={}",
                    i, seg.start_idx, seg.end_idx, seg.mean_ratio, seg.n_valid,
                );
            }
        }
        let suggested = cyrius_rs::caller::changepoint::suggest_cnv_from_changepoints(
            &cp_result, None, total_cn,
        );
        if let Some(ref tag) = suggested {
            log::info!("changepoint_hybrid: suggesting CNV tag '{}'", tag);
        }
        suggested.or(cnvtag)
    } else {
        cnvtag
    };

    // Feature: hmm_cnv — HMM-based segmentation as fallback when consensus fails
    let cnvtag = if features.hmm_cnv && cnvtag.is_none() {
        // Emission SD is measured on the per-site D6-CN profile the HMM runs
        // over. The depth GMM's `sd` describes normalised depth over the whole
        // D6+D7 region and is 4-10x too small here (review findings #15/#81).
        let sd_per_copy =
            cyrius_rs::caller::hmm_cnv::estimate_sd_per_copy(&raw_d6_cn, &snp_d6, &snp_d7);
        let hmm_params =
            cyrius_rs::caller::hmm_cnv::HmmParams::from_sd_per_copy(total_cn, sd_per_copy);
        let hmm_result = cyrius_rs::caller::hmm_cnv::hmm_segment(
            &raw_d6_cn, &snp_d6, &snp_d7, total_cn, &hmm_params,
        );
        let hmm_tag = cyrius_rs::caller::hmm_cnv::classify_segments(
            &hmm_result.segments, total_cn, raw_cn_call.spacer_cn,
        );
        log::info!(
            "hmm_cnv: profile-estimated sd_per_copy={:.3}, segmented into {} segments",
            sd_per_copy,
            hmm_result.segments.len(),
        );
        for (i, seg) in hmm_result.segments.iter().enumerate() {
            log::info!(
                "  segment {}: idx {}..{}, cn={}, n_valid={}, mean_obs={:.2}",
                i, seg.start_idx, seg.end_idx, seg.cn, seg.n_valid, seg.mean_obs,
            );
        }
        match hmm_tag {
            Some(ref tag) => log::info!("hmm_cnv: suggesting CNV tag '{}'", tag),
            None => log::info!("hmm_cnv: segment pattern did not match any known CNV group"),
        }
        hmm_tag
    } else {
        cnvtag
    };

    // Feature: het_check — detect hemizygosity (zero het sites = deletion)
    if features.het_check {
        let het_result = cyrius_rs::caller::het_check::check_hemizygosity(
            &snp_d6, &snp_d7, 10, 3, 0.1,
        );
        log::info!(
            "het_check: n_het={}, n_adequate={}, het_frac={:.3}, hemizygous={}",
            het_result.n_het_sites,
            het_result.n_adequate_depth,
            het_result.het_fraction,
            het_result.is_hemizygous,
        );
        if het_result.is_hemizygous {
            if let Some(ref tag) = cnvtag {
                log::info!(
                    "het_check: hemizygous sample with CNV tag '{}' — possible deletion",
                    tag
                );
            }
        }
    }

    // Feature: clip_evidence — detect soft-clip clusters as structural breakpoint evidence
    if features.clip_evidence {
        let chrom = &params.var_db.nchr;
        let clip_ev = cyrius_rs::caller::clip_evidence::detect_clip_clusters(
            bam_path, chrom, 42123000, 42146000, reference_fasta,
        );
        if let Some(ref ev) = clip_ev {
            let signals = cyrius_rs::caller::clip_evidence::classify_clip_signals(ev);
            log::info!(
                "clip_evidence: {} clusters, {} total clipped reads, \
                 d6_body={}, rep6_spacer={}, spacer={}, d7={}",
                signals.n_clusters, ev.total_clipped_reads,
                signals.d6_body_breakpoint, signals.rep6_spacer_breakpoint,
                signals.spacer_breakpoint, signals.d7_breakpoint,
            );
            for cluster in &ev.clusters {
                log::info!(
                    "  cluster pos={}: left={} right={}",
                    cluster.position, cluster.left_count, cluster.right_count,
                );
            }
        }
    }

    // Feature: d7_depth — CYP2D7 sub-region depth profiling for hybrid confirmation
    if features.d7_depth {
        let chrom = &params.var_db.nchr;
        if let Some(profile) = cyrius_rs::caller::d7_depth::compute_d7_depth(
            bam_path, chrom, reference_fasta,
        ) {
            let signals = cyrius_rs::caller::d7_depth::classify_d7_signals(&profile);
            log::info!(
                "d7_depth: D7_exon9={:.1}x, D7_in4_in8={:.1}x, D7_ex2_in8={:.1}x, D7_5pr_in1={:.1}x",
                profile.d7_exon9.avg_depth, profile.d7_in4_in8.avg_depth,
                profile.d7_ex2_in8.avg_depth, profile.d7_5pr_in1.avg_depth,
            );
            log::info!(
                "d7_depth: exon9/body_ratio={:.3} (elevated={}), 3pr/5pr_ratio={:.3} (depleted={})",
                signals.exon9_body_ratio, signals.exon9_elevated,
                signals.three_prime_five_prime_ratio, signals.three_prime_depleted,
            );
            if signals.exon9_elevated {
                log::info!("d7_depth: D7 exon9 elevated — consistent with *36 gene conversion");
            }
            if signals.three_prime_depleted {
                log::info!("d7_depth: D7 3' region depleted — consistent with *13 (D7→D6 conversion at exon2-intron4)");
            }
        }
    }

    // No-call due to CNV group. A tag outside the accepted vocabulary is
    // logged, never reported as the CNV_group (escalation E2).
    if reportable_cnvtag(cnvtag.as_deref()).is_none() {
        if let Some(tag) = &cnvtag {
            log::warn!("CNV tag \"{}\" is not in the accepted vocabulary -- no-call", tag);
        }
        return D6Call {
            coverage_mad: normalized_depth.mad,
            median_depth: normalized_depth.mediandepth,
            total_cn: Some(total_cn), spacer_cn: raw_cn_call.spacer_cn,
            total_cn_raw: Some(raw_cn_call.d67_depth), spacer_cn_raw: Some(raw_cn_call.spacer_depth),
            variants_called: None, cnv_group: None, genotype: None, filter: None,
            raw_star_allele: None, call_info: None,
            exon9_cn: exon9gc_call_stringent,
            cnv_consensus: Some(consensus_str),
            d67_snp_call: Some(cn_call_snp_str),
            d67_snp_raw: Some(raw_d6_cn_str),
            variant_raw_count: Some(raw_count),
            kmer_validation: None, resolver_posterior: None,
            rep7_mapq0_fraction: normalized_depth.rep7_mapq0_fraction,
        };
    }

    let cnvtag_str = cnvtag.unwrap();

    // Only the combination table this tag matches against: the other nine cost ~3 s
    // and ~1 GB and are never read (finding #41).
    let star_combinations = get_hap_table(data::STAR_TABLE, combination_table_for(&cnvtag_str).as_slice());

    // Feature: kmer_validation — alignment-free k-mer paralog ratio validation.
    // Validation is decoupled from calling: the result is recorded for review
    // but does NOT modify the genotype, filter, or confidence.
    let kmer_validation_summary: Option<cyrius_rs::types::KmerValidationSummary> =
        if features.kmer_validation {
            let kv = cyrius_rs::caller::kmer_validation::validate_with_kmers(
                bam_path,
                &params.snp_db.nchr,
                total_cn,
                raw_cn_call.spacer_cn,
                &cnvtag_str,
                reference_fasta,
                index_name,
            );
            if kv.agrees {
                log::info!(
                    "kmer_validation: AGREE — depth={} kmer_cat={}, ratio={:.3}, \
                     est_d6={:.2} est_d7={:.2}, hybrid={}, reads={}/{}, {}",
                    cnvtag_str, kv.kmer_category, kv.overall_ratio,
                    kv.estimated_d6_cn, kv.estimated_d7_cn, kv.hybrid_detected,
                    kv.matched_reads, kv.total_reads, kv.segments_summary,
                );
            } else {
                log::warn!(
                    "kmer_validation: DISAGREE — depth={} kmer_cat={}, ratio={:.3}, \
                     est_d6={:.2} est_d7={:.2}, hybrid={}, reads={}/{}, {}",
                    cnvtag_str, kv.kmer_category, kv.overall_ratio,
                    kv.estimated_d6_cn, kv.estimated_d7_cn, kv.hybrid_detected,
                    kv.matched_reads, kv.total_reads, kv.segments_summary,
                );
            }
            // Log per-position ratio profile
            let ratio_strs: Vec<String> = kv.position_ratios.iter().enumerate().map(|(i, r)| {
                match r {
                    Some(v) => format!("{}:{:.2}({}/{})", i, v, kv.d6_counts[i], kv.d7_counts[i]),
                    None => format!("{}:-", i),
                }
            }).collect();
            log::info!("kmer_validation profile: {}", ratio_strs.join(" "));

            Some(cyrius_rs::types::KmerValidationSummary {
                agrees: kv.agrees,
                kmer_category: kv.kmer_category,
                overall_ratio: kv.overall_ratio,
                estimated_d6_cn: kv.estimated_d6_cn,
                estimated_d7_cn: kv.estimated_d7_cn,
                hybrid_detected: kv.hybrid_detected,
                total_reads: kv.total_reads,
                matched_reads: kv.matched_reads,
                segments_summary: kv.segments_summary,
            })
        } else {
            None
        };

    // 5. Call variants
    // Feature: strand_bias_all — filter additional sites before CN calling
    if features.strand_bias_all {
        let filtered = strand_bias_all::apply_strand_bias_all(
            &mut var_alt, &var_alt_forward, &var_alt_reverse, &params.var_list,
        );
        if filtered > 0 {
            log::info!("strand_bias_all: filtered {} additional variant sites", filtered);
        }
    }

    let cn_call_var_homo = call_cn_var_homo(total_cn, &var_homo_alt, &var_homo_ref);
    let cn_call_var = call_cn_var(
        &cnvtag_str,
        &mut var_alt,
        &var_ref,
        &var_alt_forward,
        &var_alt_reverse,
        &params.var_list,
        &params.var_db,
    );

    // Call haplotype variants (parallel when threads > 1)
    let hap_db = &params.haplotype_db;

    let hap_calls: Result<_, String> = if threads > 1 {
        std::thread::scope(|s| {
            let h1 = s.spawn(|| call_var42126938(&mut open()?, total_cn, &hap_db["g.42126938C>T"]));
            let h2 = s.spawn(|| {
                call_var42127526_var42127556(&mut open()?, &cnvtag_str, &hap_db["g.42127526C>T_g.42127556T>C"])
            });
            let h3 = s.spawn(|| call_var42127803hap(&mut open()?, &cnvtag_str, &hap_db["g.42127803C>T"]));
            let h4 = s.spawn(|| {
                call_var42130655ins_a(&mut open()?, total_cn, &hap_db["g.42130655-42130656insA"])
            });

            let (c1, v1, g1) = h1.join().unwrap()?;
            let (c2a, c2b, v2) = h2.join().unwrap()?;
            let v3 = h3.join().unwrap()?;
            let (c4, v4) = h4.join().unwrap()?;
            Ok((c1, v1, g1, c2a, c2b, v2, v3, c4, v4))
        })
    } else {
        (|| {
            let mut reader = open()?;
            let (c1, v1, g1) = call_var42126938(&mut reader, total_cn, &hap_db["g.42126938C>T"])?;
            let (c2a, c2b, v2) = call_var42127526_var42127556(&mut reader, &cnvtag_str, &hap_db["g.42127526C>T_g.42127556T>C"])?;
            let v3 = call_var42127803hap(&mut reader, &cnvtag_str, &hap_db["g.42127803C>T"])?;
            let (c4, v4) = call_var42130655ins_a(&mut reader, total_cn, &hap_db["g.42130655-42130656insA"])?;
            Ok((c1, v1, g1, c2a, c2b, v2, v3, c4, v4))
        })()
    };
    let (site42126938_count, var42126938, var42126938_g_haplotype,
     site42127526_count, site42127556_count, var42127526,
     var42127803_diff_haplotype,
     var42130655_count, var42130655ins_a) = match hap_calls {
        Ok(v) => v,
        Err(e) => {
            log::error!("Cannot read haplotype evidence: {}", e);
            return D6Call {
                rep7_mapq0_fraction: normalized_depth.rep7_mapq0_fraction,
                ..D6Call::no_call(normalized_depth.mad, normalized_depth.mediandepth)
            };
        }
    };

    raw_count.entry("g.42126938C>T".to_string()).or_insert_with(||
        format!("{},{}", site42126938_count[1], site42126938_count[0]),
    );
    raw_count.entry("g.42127526C>T".to_string()).or_insert_with(||
        format!("{},{}", site42127526_count[1], site42127526_count[0]),
    );
    raw_count.entry("g.42127556T>C".to_string()).or_insert_with(||
        format!("{},{}", site42127556_count[1], site42127556_count[0]),
    );
    raw_count.entry("g.42130655-42130656insA".to_string()).or_insert_with(||
        format!("{},{}", var42130655_count[1], var42130655_count[0]),
    );

    // 6. Call star allele
    let mut total_callset = get_called_variants(&params.var_list, &cn_call_var, 0);
    let called_var_homo = get_called_variants(&params.var_list, &cn_call_var_homo, cn_call_var.len());
    total_callset.extend(called_var_homo);
    total_callset.extend(var42126938);
    total_callset.extend(var42127526);
    total_callset.extend(var42130655ins_a);

    let exon9_values = Exon9Values {
        exon9_cn: exon9gc_call_stringent,
        exon9cn_in_consensus: consensus.exon9_and_downstream,
        exon9_raw_site1: raw_d6_cn.get(EXON9_SITE1).copied().unwrap_or(0.0),
        exon9_raw_site2: raw_d6_cn.get(EXON9_SITE2).copied().unwrap_or(0.0),
    };

    let star_called = match_star(
        &mut total_callset,
        &cnvtag_str,
        raw_cn_call.spacer_cn,
        &star_combinations,
        &exon9_values,
        var42126938_g_haplotype,
        var42127803_diff_haplotype,
        features,
    );

    let mut genotype_filter = None;
    let mut final_star_allele_call = None;
    let mut resolver_posterior: Option<f64> = None;
    let is_fuzzy = star_called
        .call_info
        .as_deref()
        .map_or(false, |s| s.starts_with("fuzzy_match"));

    if star_called.call_info.is_some()
        && star_called.call_info.as_deref() != Some("no_match")
        && star_called.call_info.as_deref() != None
    {
        final_star_allele_call = star_called.clean_call.clone();

        // Feature: phase_disambiguate — resolve ambiguous CN=2 calls using physical phasing
        if features.phase_disambiguate
            && star_called.call_info.as_deref() == Some("more_than_one_match")
            && cnvtag_str == "cn2"
        {
            if let (Some(ref call), Some(ref raw_call)) =
                (&final_star_allele_call, &star_called.raw_call)
            {
                if call.contains(';') {
                    let variant_lookup = cyrius_rs::caller::phase_disambiguate::VariantLookup::from_target_file(data::TARGET_VARIANT);
                    if let Some(resolved) = cyrius_rs::caller::phase_disambiguate::disambiguate(
                        call,
                        raw_call,
                        &params.star_combinations,
                        &variant_lookup,
                        bam_path,
                        &params.var_db.nchr,
                        reference_fasta,
                        index_name,
                        features.phase_readpair,
                    ) {
                        log::info!("phase_disambiguate: {} -> {}", call, resolved);
                        final_star_allele_call = Some(resolved);
                    }
                }
            }
        }

        // Feature: likelihood resolver — when the cascade leaves a ;-ambiguous
        // candidate set, score each candidate's allele multiset and either pick
        // confidently or report honest ambiguity. Falls back silently to the
        // cascade's own choice if the BAM can't be scored.
        if let Some(ref call) = final_star_allele_call {
            if call.contains(';') {
                let candidates: Vec<String> =
                    call.split(';').map(|s| s.trim().to_string()).collect();
                let raw_candidates = star_called.raw_call.clone().unwrap_or_default();
                if let Ok(mut r) =
                    open_alignment_file_with_index(bam_path, reference_fasta, index_name)
                {
                    let scores = cyrius_rs::caller::diplotype_caller::score_candidates(
                        &mut r,
                        &params.snp_db.nchr,
                        &params.var_list,
                        data::SNP_FILE,
                        &params.star_combinations,
                        &candidates,
                        &raw_candidates,
                    );
                    for s in &scores {
                        log::info!(
                            "resolver: candidate {} log10_lik={:.2} posterior={:.4} decisive_reads={}",
                            s.diplotype, s.log10_likelihood, s.posterior, s.decisive_reads
                        );
                    }
                    if let Some(outcome) = cyrius_rs::caller::resolver::resolve(&scores) {
                        log::info!(
                            "resolver: {} -> {} (ambiguous={}, top_posterior={:.4})",
                            call, outcome.genotype, outcome.ambiguous, outcome.top_posterior
                        );
                        resolver_posterior = Some(outcome.top_posterior);
                        final_star_allele_call = Some(outcome.genotype);
                    } else {
                        log::info!(
                            "resolver: {} -> not scorable ({} candidates), keeping cascade call",
                            call, candidates.len()
                        );
                    }
                } else {
                    log::warn!("resolver: BAM reopen failed for {}, keeping cascade call", bam_path);
                }
            }
        }

        if let Some(ref call) = final_star_allele_call {
            if call.contains(';') {
                genotype_filter = Some("More_than_one_possible_genotype".to_string());
            } else if !call.contains('/') {
                genotype_filter = Some("Not_assigned_to_haplotypes".to_string());
            } else if high_cn_low_confidence {
                genotype_filter = Some("LowQ_high_CN".to_string());
            } else if is_fuzzy {
                genotype_filter = Some("Fuzzy_match".to_string());
            } else {
                genotype_filter = Some("PASS".to_string());
            }
        }
    }

    // Feature: spacer_cn_check — flag calls where spacer CN is inconsistent with CNV group
    // total_cn - spacer_cn = D6 exon9 copy number. Expected values per CNV group:
    //   cn2=2, cn3=3, cn4=4, cn5=5, cn6=6
    //   star5=1, star5_star5=0
    //   exon9hyb=2, exon9hyb_exon9hyb=2, exon9hyb_exon9hyb_exon9hyb=2
    //   dup_exon9hyb=3, dup_star68=3
    //   star68: skip (variable, 1 or 2)
    //   star13/star13intron1 variants: skip (insufficient data)
    if features.spacer_cn_check {
        if let Some(ref filter) = genotype_filter {
            if filter == "PASS" || filter == "Fuzzy_match" {
                if let Some(sp_cn) = raw_cn_call.spacer_cn {
                    let expected_diff: Option<u32> = match cnvtag_str.as_str() {
                        "cn2" => Some(2),
                        "cn3" => Some(3),
                        "cn4" => Some(4),
                        "cn5" => Some(5),
                        "cn6" => Some(6),
                        "star5" => Some(1),
                        "star5_star5" => Some(0),
                        "exon9hyb" | "exon9hyb_exon9hyb"
                        | "exon9hyb_exon9hyb_exon9hyb"
                        | "exon9hyb_exon9hyb_exon9hyb_exon9hyb" => Some(2),
                        "exon9hyb_star5" => Some(1),
                        "dup_exon9hyb" | "dup_star68" => Some(3),
                        "exon9hyb_star68" => Some(2),
                        "star68_star68" | "star68_star68_star68"
                        | "star68_star68_star68_star68" => Some(2),
                        _ => None, // star68, star13 variants — skip
                    };
                    if let Some(expected) = expected_diff {
                        let observed = total_cn.saturating_sub(sp_cn);
                        if observed != expected {
                            log::info!(
                                "spacer_cn_check: inconsistent — CNV group '{}' expects \
                                 total_cn - spacer_cn = {}, got {} - {} = {}, flagging as uncertain",
                                cnvtag_str, expected, total_cn, sp_cn, observed,
                            );
                            genotype_filter = Some("Spacer_CN_inconsistent".to_string());
                        }
                    }
                }
            }
        }
    }

    // Feature: consistency_check — post-hoc validation of the final call
    if features.consistency_check {
        if let Some(ref filter) = genotype_filter {
            if filter == "PASS" || filter == "Fuzzy_match" {
                let consistency = cyrius_rs::caller::consistency::run_consistency_checks(
                    &snp_d6, &snp_d7,
                    &var_alt, &var_ref, &params.var_list,
                    &cnvtag_str, total_cn,
                    final_star_allele_call.as_deref(),
                    &params.star_combinations,
                    &snp_crossing,
                );

                // Pick the first failing check to use as the filter
                if let Some(ref flag) = consistency.crossing_flag {
                    genotype_filter = Some(flag.clone());
                } else if let Some(ref flag) = consistency.conversion_map_flag {
                    genotype_filter = Some(flag.clone());
                } else if let Some(ref flag) = consistency.mismatch_flag {
                    genotype_filter = Some(flag.clone());
                } else if let Some(ref flag) = consistency.balance_flag {
                    genotype_filter = Some(flag.clone());
                }
            }
        }
    }

    // Feature: diplotype_caller — likelihood-based diplotype calling with D6+D7 mixture model
    if features.diplotype_caller {
        let mut dc_reader =
            open_alignment_file_with_index(bam_path, reference_fasta, index_name).unwrap();
        let dc_result = cyrius_rs::caller::diplotype_caller::call_diplotype(
            &mut dc_reader,
            &params.snp_db.nchr,
            &params.var_list,
            data::SNP_FILE,
            &params.star_combinations,
            total_cn,
            20,
        );
        if !dc_result.top_diplotypes.is_empty() {
            let best = &dc_result.top_diplotypes[0];
            log::info!(
                "diplotype_caller: best={}/{} (posterior={:.4}, ll={:.1}), \
                 pileup_call={:?}, reads={} (filtered={})",
                best.allele_a,
                best.allele_b,
                best.posterior,
                best.log10_likelihood,
                final_star_allele_call,
                dc_result.n_informative_reads,
                dc_result.n_filtered_reads,
            );
        }
    }

    D6Call {
        coverage_mad: normalized_depth.mad,
        median_depth: normalized_depth.mediandepth,
        total_cn: Some(total_cn),
        spacer_cn: raw_cn_call.spacer_cn,
        total_cn_raw: Some(raw_cn_call.d67_depth),
        spacer_cn_raw: Some(raw_cn_call.spacer_depth),
        variants_called: star_called
            .variants_called
            .as_ref()
            .map(|v| v.split_whitespace().map(|s| s.to_string()).collect()),
        cnv_group: Some(cnvtag_str),
        genotype: final_star_allele_call,
        filter: genotype_filter,
        raw_star_allele: star_called.raw_call,
        call_info: star_called.call_info,
        exon9_cn: exon9gc_call_stringent,
        cnv_consensus: Some(consensus_str),
        d67_snp_call: Some(cn_call_snp_str),
        d67_snp_raw: Some(raw_d6_cn_str),
        variant_raw_count: Some(raw_count),
        kmer_validation: kmer_validation_summary,
        resolver_posterior,
        rep7_mapq0_fraction: normalized_depth.rep7_mapq0_fraction,
    }
}

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .init();

    let cli = Cli::parse();

    let (bam_name, index_name) = if cli.input.contains("##idx##") {
        let parts: Vec<&str> = cli.input.splitn(2, "##idx##").collect();
        (parts[0].to_string(), Some(parts[1].to_string()))
    } else {
        (cli.input.clone(), None)
    };

    // Derive the default sample id from the BAM half of "input##idx##index".
    let sample_id = cli.id.unwrap_or_else(|| {
        Path::new(&bam_name)
            .file_stem()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| "sample".to_string())
    });

    if let Err(e) = std::fs::create_dir_all(&cli.out_dir) {
        eprintln!("Cannot create output directory {}: {}", cli.out_dir, e);
        std::process::exit(1);
    }

    let out_json = format!("{}/{}.json", cli.out_dir, sample_id);
    let out_tsv = format!("{}/{}.tsv", cli.out_dir, sample_id);

    let count_file = cli
        .count_file_path
        .as_ref()
        .map(|p| format!("{}/{}_count.txt", p, sample_id));

    if !bam_name.contains("://") && !Path::new(&bam_name).exists() {
        log::error!("Input file for sample {} does not exist.", sample_id);
        // Write empty output files (matching Python behavior)
        write_file_or_exit(&out_json, "{}");
        let header = "Sample\tGenotype\tFilter\tConfidence\tActivity score\tPredicted phenotype\tPhasing_ambiguous\tLong_read_recommended\n";
        write_file_or_exit(&out_tsv, header);
        return;
    }

    // Size rayon's global pool from --threads before any rayon work can start.
    // A global pool can only be built once, and whichever rayon call runs first
    // builds the default all-cores pool if we have not, so this has to happen
    // upstream of every rayon user. There are exactly two -- the control-region
    // loop in `bin_count` and `diplotype_caller`'s `par_chunks` -- and both sit
    // inside `d6_star_caller`, below. Nothing between `Cli::parse()` and here
    // touches rayon either (string splitting, the sample-id default,
    // `create_dir_all`, path formatting, the existence check and, on the
    // early-return path, `write_file_or_exit`), which is why this sits *after*
    // the "input does not exist" return rather than straight after
    // `Cli::parse()`: a no-op run on a missing input used to spawn a pool it
    // never used. Review finding #10.
    let pool_threads = cyrius_rs::threading::init_rayon_global_pool(cli.threads);

    // Prepare data (after existence check, with the BAM half of ##idx##)
    log::info!("Starting {} v{}", TOOL_NAME, TOOL_VERSION);
    let (call_parameters, header_contigs) =
        prepare_resource(&cli.genome, &bam_name, index_name.as_deref());

    // Every region and SNP lookup addresses the reference by contig name, and
    // matching is exact: with `-g 38` the target is `22`, which a `chr22`
    // header does not satisfy. When the CYP2D6 contig is absent, every lookup
    // *on that contig* silently yields nothing and the run ends in a no-call
    // that cannot be told apart from a genuine low-quality one. Lookups on
    // contigs the header does hold still return reads, so the depth path's
    // "control regions usable" error is not a reliable substitute, and on the
    // `--countFilePath` path -- where the BAM is still read for SNP counting
    // -- nothing at all is logged. Say so loudly, naming the sample, because
    // batch runs fan many samples into one stream; the no-call itself stays.
    // `prepare_resource` already read this header for the genome autodetect,
    // so there is no second open.
    //
    // Only `snp_db` and `var_db` are checked. `var_homo_db` and the
    // `region_dic` contigs are derived from the same resolved `genome` and so
    // agree with them today; widening the check would add duplicate noise, not
    // coverage.
    if let Some(contigs) = header_contigs.as_deref() {
        let targets =
            contig_check_targets(&call_parameters.snp_db.nchr, &call_parameters.var_db.nchr);
        for target in targets {
            if let Some(message) = missing_contig_error(
                &sample_id,
                &bam_name,
                &call_parameters.genome,
                &target,
                contigs,
            ) {
                log::error!("{}", message);
            }
        }
    }

    let features = FeatureFlags {
        strand_bias_all: cli.strand_bias_all,
        fuzzy_match: cli.fuzzy_match,
        quality_aware: cli.quality_aware,
        phase_disambiguate: cli.phase_disambiguate,
        phase_readpair: cli.phase_readpair,
        changepoint_hybrid: cli.changepoint_hybrid,
        het_check: cli.het_check,
        spacer_cn_check: cli.spacer_cn_check,
        consistency_check: cli.consistency_check,
        hmm_cnv: cli.hmm_cnv,
        diplotype_caller: cli.diplotype_caller,
        clip_evidence: cli.clip_evidence,
        d7_depth: cli.d7_depth,
        kmer_validation: cli.kmer_validation,
    };

    if features.strand_bias_all || features.fuzzy_match || features.quality_aware
        || features.phase_disambiguate || features.phase_readpair
        || features.changepoint_hybrid || features.het_check || features.spacer_cn_check
        || features.consistency_check || features.hmm_cnv
        || features.diplotype_caller || features.d7_depth
        || features.kmer_validation
    {
        log::info!(
            "Experimental features: strand_bias_all={}, fuzzy_match={}, quality_aware={}, \
             phase_disambiguate={}, phase_readpair={}, changepoint_hybrid={}, het_check={}, \
             spacer_cn_check={}, consistency_check={}, hmm_cnv={}, \
             diplotype_caller={}, d7_depth={}, kmer_validation={}",
            features.strand_bias_all, features.fuzzy_match, features.quality_aware,
            features.phase_disambiguate, features.phase_readpair,
            features.changepoint_hybrid, features.het_check, features.spacer_cn_check,
            features.consistency_check, features.hmm_cnv,
            features.diplotype_caller, features.d7_depth, features.kmer_validation
        );
    }

    if features.quality_aware {
        log::warn!(
            "--quality-aware: module is present but not yet wired into the calling path; \
             outputs are identical to runs without this flag."
        );
    }

    match pool_threads {
        Some(n) => log::info!("Processing sample {} with {} rayon thread(s)", sample_id, n),
        None => log::info!(
            "Processing sample {} (rayon thread pool size unknown)",
            sample_id
        ),
    }

    // Now that --threads really bounds rayon, --diplotype-caller at the default
    // -t 1 is single-threaded where it used to borrow every core by accident:
    // ~2.5x slower in wall time (for ~28% less total CPU). Say so here, before
    // the slow part starts, rather than leaving the user to infer it from the
    // "1 rayon thread(s)" line above. No CLI surface change: --help still reads
    // "Number of threads".
    if features.diplotype_caller && pool_threads == Some(1) {
        log::info!(
            "--diplotype-caller scores diplotypes on rayon's thread pool, which \
             --threads has bounded to 1 here; pass -t N to use more cores (-t 8 \
             is roughly 2.5x faster in wall time on a large sample)."
        );
    }

    let mut cyp2d6_call = d6_star_caller(
        &bam_name,
        &call_parameters,
        cli.threads,
        count_file.as_deref(),
        cli.reference.as_deref(),
        index_name.as_deref(),
        &features,
    );

    if cyp2d6_call.coverage_mad > MAD_THRESHOLD {
        log::warn!(
            "Sample {} has uneven coverage. CN calls may be unreliable.",
            sample_id
        );
        cyp2d6_call.filter = cyrius_rs::depth_calling::bin_count::downgrade_filter_for_high_mad(
            cyp2d6_call.filter,
            cyp2d6_call.coverage_mad,
            MAD_THRESHOLD,
        );
    }

    // Input guard: an aligner that dropped the reads fitting CYP2D6 and CYP2D7
    // equally well (MAPQ 0) changed the evidence, so the call is not a PASS.
    // The header says whether a pipeline step dropped them.
    let header_drops_multimappers =
        open_alignment_file_with_index(&bam_name, cli.reference.as_deref(), index_name.as_deref())
            .ok()
            .map(|reader| {
                cyrius_rs::depth_calling::bin_count::header_drops_multimappers(
                    &String::from_utf8_lossy(rust_htslib::bam::Read::header(&reader).as_bytes()),
                )
            });
    let guarded_filter = cyrius_rs::depth_calling::bin_count::downgrade_filter_for_multimapper_depletion(
        cyp2d6_call.filter.clone(),
        header_drops_multimappers,
    );
    if guarded_filter != cyp2d6_call.filter {
        log::warn!(
            "Sample {}: the BAM header shows a step that drops multi-mapped reads (novoalign \
             without -r, or samtools view -q); REP7 MAPQ-0 fraction {:.3}. The call is not reliable.",
            sample_id,
            cyp2d6_call.rep7_mapq0_fraction.unwrap_or(f64::NAN),
        );
        cyp2d6_call.filter = guarded_filter;
    }

    // Compute confidence score
    let confidence = cyrius_rs::caller::confidence::compute_confidence(
        &cyrius_rs::caller::confidence::ConfidenceInput {
            coverage_mad: cyp2d6_call.coverage_mad,
            median_depth: cyp2d6_call.median_depth,
            total_cn_raw: cyp2d6_call.total_cn_raw.unwrap_or(0.0),
            spacer_cn_raw: cyp2d6_call.spacer_cn_raw.unwrap_or(0.0),
            call_info: cyp2d6_call.call_info.as_deref(),
            filter: cyp2d6_call.filter.as_deref(),
            genotype: cyp2d6_call.genotype.as_deref(),
            cnv_group: cyp2d6_call.cnv_group.as_deref(),
            d67_snp_raw: cyp2d6_call.d67_snp_raw.as_deref(),
            variant_raw_count: cyp2d6_call.variant_raw_count.as_ref(),
            star_combinations: &call_parameters.star_combinations,
        },
    );
    let c = &confidence.components;
    log::info!(
        "confidence: {:.2} ({}) — depth={:.2} cn={:.2} match={:.2} completeness={:.2} specificity={:.2} snp={:.2}",
        confidence.score, confidence.label,
        c.depth_quality, c.cn_quality, c.match_quality,
        c.variant_completeness, c.variant_specificity, c.snp_consistency,
    );

    // Compute phenotype predictions and phasing status (needed for both JSON and TSV)
    let haplotype_functionality =
        phenotype::load_haplotype_functionality(data::HAPLOTYPE_FUNC);
    let sorted_genotype = cyp2d6_call.genotype.as_deref().map(phenotype::sort_diplotypes);
    let predictions = phenotype::match_phenotype(
        sorted_genotype.as_deref(),
        &haplotype_functionality,
    );

    let count_of_diplotypes = predictions.len();

    // Determine phasing ambiguity and clinical equivalence.
    // Per README: long-read reflex fires when alternatives differ in
    // activity score OR predicted phenotype — both axes matter — or when
    // any alternative's effect is unknown.
    let phasing_ambiguous = count_of_diplotypes > 1;
    let (clinical_equivalence, long_read_recommended) = if phasing_ambiguous {
        phenotype::clinical_equivalence(&predictions)
    } else {
        ("n/a".to_string(), false)
    };

    if phasing_ambiguous {
        log::info!(
            "phasing_ambiguous=true clinical_equivalence={} long_read_recommended={}",
            clinical_equivalence, long_read_recommended,
        );
    }

    // Write JSON
    log::info!("Writing to json ({})", out_json);
    let sorted_json_genotype = sorted_genotype.clone();
    // K-mer validation summary (None unless --kmer-validation was passed and
    // a depth-based CNV tag was produced). Surfaces the alignment-free QC
    // verdict alongside the call; does not modify the call itself.
    let kmer_validation_json = cyp2d6_call.kmer_validation.as_ref().map(|kv| {
        serde_json::json!({
            "agrees": kv.agrees,
            "kmer_category": kv.kmer_category,
            "overall_ratio": kv.overall_ratio,
            "estimated_d6_cn": kv.estimated_d6_cn,
            "estimated_d7_cn": kv.estimated_d7_cn,
            "hybrid_detected": kv.hybrid_detected,
            "matched_reads": kv.matched_reads,
            "total_reads": kv.total_reads,
            "segments_summary": kv.segments_summary,
        })
    });

    let mut json_map = serde_json::Map::new();
    let call_json = serde_json::json!({
        "Coverage_MAD": cyp2d6_call.coverage_mad,
        "Median_depth": cyp2d6_call.median_depth,
        "REP7_MAPQ0_fraction": cyp2d6_call.rep7_mapq0_fraction,
        "Header_drops_multimappers": header_drops_multimappers,
        "Total_CN": cyp2d6_call.total_cn,
        "Spacer_CN": cyp2d6_call.spacer_cn,
        "Total_CN_raw": cyp2d6_call.total_cn_raw,
        "Spacer_CN_raw": cyp2d6_call.spacer_cn_raw,
        "Variants_called": cyp2d6_call.variants_called,
        "CNV_group": cyp2d6_call.cnv_group,
        "Genotype": sorted_json_genotype,
        "Filter": cyp2d6_call.filter,
        "Confidence": format!("{:.2}", confidence.score),
        "Confidence_label": confidence.label,
        "Phasing_ambiguous": phasing_ambiguous,
        "Clinical_equivalence": if phasing_ambiguous { &clinical_equivalence } else { "n/a" },
        "Long_read_recommended": long_read_recommended,
        "Raw_star_allele": cyp2d6_call.raw_star_allele,
        "Call_info": cyp2d6_call.call_info,
        "Exon9_CN": cyp2d6_call.exon9_cn,
        "CNV_consensus": cyp2d6_call.cnv_consensus,
        "d67_snp_call": cyp2d6_call.d67_snp_call,
        "d67_snp_raw": cyp2d6_call.d67_snp_raw,
        "Variant_raw_count": cyp2d6_call.variant_raw_count,
        "Kmer_validation": kmer_validation_json,
        "Resolver_posterior": cyp2d6_call.resolver_posterior,
    });
    json_map.insert(sample_id.clone(), call_json);
    let json_output = serde_json::to_string(&json_map).unwrap();
    write_file_or_exit(&out_json, &json_output);

    // Write TSV
    log::info!("Writing to tsv ({})", out_tsv);

    let (activity_scores, predicted_phenotypes) = if !predictions.is_empty() {
        if count_of_diplotypes > 1 {
            (
                predictions
                    .iter()
                    .map(|p| p.total_activity.clone())
                    .collect::<Vec<_>>()
                    .join(";"),
                predictions
                    .iter()
                    .map(|p| p.predicted_phenotype.clone())
                    .collect::<Vec<_>>()
                    .join(";"),
            )
        } else {
            (
                predictions[0].total_activity.clone(),
                predictions[0].predicted_phenotype.clone(),
            )
        }
    } else {
        ("-".to_string(), "-".to_string())
    };

    let genotype_str = sorted_genotype.as_deref().unwrap_or("None");

    let mut tsv_content = String::new();
    tsv_content.push_str("Sample\tGenotype\tFilter\tConfidence\tActivity score\tPredicted phenotype\tPhasing_ambiguous\tLong_read_recommended\n");
    tsv_content.push_str(&format!(
        "{}\t{}\t{}\t{:.2} ({})\t{}\t{}\t{}\t{}\n",
        sample_id,
        genotype_str,
        cyp2d6_call.filter.as_deref().unwrap_or("None"),
        confidence.score,
        confidence.label,
        activity_scores,
        predicted_phenotypes,
        phasing_ambiguous,
        long_read_recommended,
    ));

    // Haplotype info
    if cli.haplotype_info && genotype_str != "None" {
        tsv_content.push_str(
            "\nHaplotype\tActivity value\tFunction\tEvidence strength\tEvidence summary\n",
        );
        for prediction in &predictions {
            for hap in &prediction.haplotype_details {
                tsv_content.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\n",
                    hap.haplotype, hap.activity, hap.function, hap.evidence_strength,
                    hap.evidence_summary,
                ));
            }
        }
    }

    // Population frequencies
    let mut frequency_table = None;
    if cli.population_info && genotype_str != "None" {
        tsv_content.push('\n');
        match phenotype::diplotype_frequencies(data::FREQ_HAP, data::FREQ_DIP, genotype_str) {
            Ok(rows) => {
                // Write header
                if let Some(first_row) = rows.first() {
                    let cols: Vec<&str> = first_row.values.keys().map(|s| s.as_str()).collect();
                    tsv_content.push_str(&format!(
                        "Biogeographic group\t{}\n",
                        cols.join("\t")
                    ));
                }
                for row in &rows {
                    let vals: Vec<&str> = row.values.values().map(|s| s.as_str()).collect();
                    tsv_content.push_str(&format!(
                        "{}\t{}\n",
                        row.biogeographic_group,
                        vals.join("\t")
                    ));
                }
                frequency_table = Some(rows);
            }
            Err(e) => {
                tsv_content.push_str(&e);
                tsv_content.push('\n');
            }
        }
    }

    write_file_or_exit(&out_tsv, &tsv_content);

    // Print results
    if cli.print {
        let column_width = 20;
        println!("\n========================================================");
        println!("                     RESULTS");
        println!("========================================================");
        println!("        ");

        let headers = ["Sample", "Genotype", "Filter", "Activity score", "Predicted phenotype"];
        let values = [
            &sample_id,
            &genotype_str.to_string(),
            &cyp2d6_call
                .filter
                .as_deref()
                .unwrap_or("None")
                .to_string(),
            &activity_scores,
            &predicted_phenotypes,
        ];
        for (h, v) in headers.iter().zip(values.iter()) {
            println!("{:<width$} → {}", h, v, width = column_width);
        }
        println!();

        if cli.haplotype_info && genotype_str != "None" {
            let mut processed_haplotypes = std::collections::HashSet::new();
            for (index, prediction) in predictions.iter().enumerate() {
                println!("-------------- Details for each haplotype --------------");
                if index > 0 {
                    println!("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
                    println!("Alternative diplotype solution:");
                }
                for hap in &prediction.haplotype_details {
                    if processed_haplotypes.contains(&hap.haplotype) {
                        continue;
                    }
                    println!("{}", hap.haplotype);
                    println!(
                        "{:<width$} → {}",
                        "Activity value",
                        hap.activity,
                        width = column_width
                    );
                    println!(
                        "{:<width$} → {}",
                        "Function",
                        hap.function,
                        width = column_width
                    );
                    if hap.evidence_strength != "n/a" {
                        println!(
                            "{:<width$} → {}",
                            "Evidence strength",
                            hap.evidence_strength,
                            width = column_width
                        );
                        println!(
                            "{:<width$} → {}",
                            "Evidence summary",
                            hap.evidence_summary,
                            width = column_width
                        );
                    }
                    println!();
                    processed_haplotypes.insert(hap.haplotype.clone());
                }
            }
        }

        if cli.population_info {
            if let Some(ref rows) = frequency_table {
                println!("---------------- Population frequencies ----------------");
                if let Some(first_row) = rows.first() {
                    let cols: Vec<&str> = first_row.values.keys().map(|s| s.as_str()).collect();
                    print!("{:>35}", "Biogeographic group");
                    for col in &cols {
                        print!("  {:>8}", col);
                    }
                    println!();
                    for row in rows {
                        print!("{:>35}", row.biogeographic_group);
                        for val in row.values.values() {
                            print!("  {:>8}", val);
                        }
                        println!();
                    }
                }
            }
        }
    }
}

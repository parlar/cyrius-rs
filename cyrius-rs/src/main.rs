use clap::Parser;
use indexmap::IndexMap;
use log;
use std::collections::HashMap;
use std::path::Path;

use cyrius_rs::caller::call_variants::{
    call_cn_snp, call_cn_var, call_cn_var_homo, call_exon9gc, call_var42126938,
    call_var42127526_var42127556, call_var42127803hap, call_var42130655ins_a,
    get_allele_counts_var42128936, get_called_variants, update_var42128936, NOISY_VAR,
};
use cyrius_rs::caller::cnv_hybrid::get_cnvtag;
use cyrius_rs::caller::construct_star_table::get_hap_table;
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
    get_snp_position, get_var_names, open_alignment_file_with_index, parse_gmm_file, parse_region_file,
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

const CNV_ACCEPTED: &[&str] = &[
    "star5_star5",
    "star13_star13",
    "star13intron1_star13intron1",
    "star5",
    "star13",
    "star13intron1",
    "star5_star5_star68",
    "star5_star68",
    "cn2",
    "exon9hyb_star5",
    "dup_star13",
    "dup_star13intron1",
    "star13_star68",
    "cn3",
    "exon9hyb",
    "star68",
    "cn4",
    "exon9hyb_exon9hyb",
    "star68_star68",
    "dup_exon9hyb",
    "dup_star68",
    "exon9hyb_star68",
    "cn5",
    "exon9hyb_exon9hyb_exon9hyb",
    "star68_star68_star68",
    "cn6",
    "exon9hyb_exon9hyb_exon9hyb_exon9hyb",
    "star68_star68_star68_star68",
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

    #[arg(long = "quality-aware", help = "Use base-quality-aware likelihoods instead of fixed error rate")]
    quality_aware: bool,
}

/// Check if chromosome names contain "chr" prefix.
fn is_chr_in_chromosome_names(input_file: &str) -> bool {
    use rust_htslib::bam::{self, Read as BamRead};
    let reader = bam::IndexedReader::from_path(input_file).unwrap();
    let header = reader.header().clone();
    for contig_id in 0..5u32.min(header.target_count()) {
        let name = String::from_utf8_lossy(header.tid2name(contig_id));
        if name.contains("chr") {
            return true;
        }
    }
    false
}

/// Prepare all resources from embedded data files.
fn prepare_resource(genome_arg: &str, input_file: &str) -> ResourceInfo {
    let genome = if genome_arg == "autodetect" {
        if is_chr_in_chromosome_names(input_file) {
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

    let star_combinations = get_hap_table(data::STAR_TABLE);

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
    }
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
        let mut reader = open_alignment_file_with_index(bam_path, reference_fasta, index_name).unwrap();
        let read_length =
            cyrius_rs::depth_calling::bin_count::get_read_length(&mut reader);
        get_normed_depth_from_count(cf, &params.region_dic, read_length)
    } else {
        get_normed_depth(bam_path, &params.region_dic, threads, reference_fasta)
    };

    // No-call after normalization
    if normalized_depth.normalized.get("d67").and_then(|v| *v).is_none() {
        return D6Call {
            coverage_mad: normalized_depth.mad,
            median_depth: normalized_depth.mediandepth,
            total_cn: None, spacer_cn: None, total_cn_raw: None, spacer_cn_raw: None,
            variants_called: None, cnv_group: None, genotype: None, filter: None,
            raw_star_allele: None, call_info: None, exon9_cn: None, cnv_consensus: None,
            d67_snp_call: None, d67_snp_raw: None, variant_raw_count: None,
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

    // 3. Get allele counts (parallel when threads > 1)
    let (snp_d6, snp_d7, mut var_alt, mut var_ref, var_alt_forward, var_alt_reverse,
     ref_read, long_ins_read, short_ins_read, var_homo_alt, var_homo_ref) = if threads > 1 {
        std::thread::scope(|s| {
            let h_snp = s.spawn(|| {
                let mut r = open_alignment_file_with_index(bam_path, reference_fasta, index_name).unwrap();
                get_supporting_reads(&mut r, &params.snp_db)
            });
            let h_var = s.spawn(|| {
                let mut r = open_alignment_file_with_index(bam_path, reference_fasta, index_name).unwrap();
                get_supporting_reads_single_region(&mut r, &params.var_db, None)
            });
            let h_ins = s.spawn(|| {
                let mut r = open_alignment_file_with_index(bam_path, reference_fasta, index_name).unwrap();
                get_allele_counts_var42128936(&mut r, &params.genome)
            });
            let h_homo = s.spawn(|| {
                let mut r = open_alignment_file_with_index(bam_path, reference_fasta, index_name).unwrap();
                get_supporting_reads(&mut r, &params.var_homo_db)
            });

            let (snp_d6, snp_d7) = h_snp.join().unwrap();
            let (var_alt, var_ref, var_alt_forward, var_alt_reverse) = h_var.join().unwrap();
            let (ref_read, long_ins_read, short_ins_read) = h_ins.join().unwrap();
            let (var_homo_alt, var_homo_ref) = h_homo.join().unwrap();
            (snp_d6, snp_d7, var_alt, var_ref, var_alt_forward, var_alt_reverse,
             ref_read, long_ins_read, short_ins_read, var_homo_alt, var_homo_ref)
        })
    } else {
        let mut reader = open_alignment_file_with_index(bam_path, reference_fasta, index_name).unwrap();
        let (snp_d6, snp_d7) = get_supporting_reads(&mut reader, &params.snp_db);
        let (var_alt, var_ref, var_alt_forward, var_alt_reverse) =
            get_supporting_reads_single_region(&mut reader, &params.var_db, None);
        let (ref_read, long_ins_read, short_ins_read) =
            get_allele_counts_var42128936(&mut reader, &params.genome);
        let (var_homo_alt, var_homo_ref) = get_supporting_reads(&mut reader, &params.var_homo_db);
        (snp_d6, snp_d7, var_alt, var_ref, var_alt_forward, var_alt_reverse,
         ref_read, long_ins_read, short_ins_read, var_homo_alt, var_homo_ref)
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

    // No-call due to CNV group
    if cnvtag.is_none() || !CNV_ACCEPTED.contains(&cnvtag.as_deref().unwrap_or("")) {
        return D6Call {
            coverage_mad: normalized_depth.mad,
            median_depth: normalized_depth.mediandepth,
            total_cn: Some(total_cn), spacer_cn: raw_cn_call.spacer_cn,
            total_cn_raw: Some(raw_cn_call.d67_depth), spacer_cn_raw: Some(raw_cn_call.spacer_depth),
            variants_called: None, cnv_group: cnvtag, genotype: None, filter: None,
            raw_star_allele: None, call_info: None,
            exon9_cn: exon9gc_call_stringent,
            cnv_consensus: Some(consensus_str),
            d67_snp_call: Some(cn_call_snp_str),
            d67_snp_raw: Some(raw_d6_cn_str),
            variant_raw_count: Some(raw_count),
        };
    }

    let cnvtag_str = cnvtag.unwrap();

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

    let (site42126938_count, var42126938, var42126938_g_haplotype,
     site42127526_count, site42127556_count, var42127526,
     var42127803_diff_haplotype,
     var42130655_count, var42130655ins_a) = if threads > 1 {
        std::thread::scope(|s| {
            let h1 = s.spawn(|| {
                let mut r = open_alignment_file_with_index(bam_path, reference_fasta, index_name).unwrap();
                call_var42126938(&mut r, total_cn, &hap_db["g.42126938C>T"])
            });
            let h2 = s.spawn(|| {
                let mut r = open_alignment_file_with_index(bam_path, reference_fasta, index_name).unwrap();
                call_var42127526_var42127556(&mut r, &cnvtag_str, &hap_db["g.42127526C>T_g.42127556T>C"])
            });
            let h3 = s.spawn(|| {
                let mut r = open_alignment_file_with_index(bam_path, reference_fasta, index_name).unwrap();
                call_var42127803hap(&mut r, &cnvtag_str, &hap_db["g.42127803C>T"])
            });
            let h4 = s.spawn(|| {
                let mut r = open_alignment_file_with_index(bam_path, reference_fasta, index_name).unwrap();
                call_var42130655ins_a(&mut r, total_cn, &hap_db["g.42130655-42130656insA"])
            });

            let (c1, v1, g1) = h1.join().unwrap();
            let (c2a, c2b, v2) = h2.join().unwrap();
            let v3 = h3.join().unwrap();
            let (c4, v4) = h4.join().unwrap();
            (c1, v1, g1, c2a, c2b, v2, v3, c4, v4)
        })
    } else {
        let mut reader = open_alignment_file_with_index(bam_path, reference_fasta, index_name).unwrap();
        let (c1, v1, g1) = call_var42126938(&mut reader, total_cn, &hap_db["g.42126938C>T"]);
        let (c2a, c2b, v2) = call_var42127526_var42127556(&mut reader, &cnvtag_str, &hap_db["g.42127526C>T_g.42127556T>C"]);
        let v3 = call_var42127803hap(&mut reader, &cnvtag_str, &hap_db["g.42127803C>T"]);
        let (c4, v4) = call_var42130655ins_a(&mut reader, total_cn, &hap_db["g.42130655-42130656insA"]);
        (c1, v1, g1, c2a, c2b, v2, v3, c4, v4)
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
        &params.star_combinations,
        &exon9_values,
        var42126938_g_haplotype,
        var42127803_diff_haplotype,
        features,
    );

    let mut genotype_filter = None;
    let mut final_star_allele_call = None;
    let is_fuzzy = star_called
        .call_info
        .as_deref()
        .map_or(false, |s| s.starts_with("fuzzy_match"));

    if star_called.call_info.is_some()
        && star_called.call_info.as_deref() != Some("no_match")
        && star_called.call_info.as_deref() != None
    {
        final_star_allele_call = star_called.clean_call.clone();
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
    }
}

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .init();

    let cli = Cli::parse();

    let sample_id = cli.id.unwrap_or_else(|| {
        Path::new(&cli.input)
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .to_string()
    });

    std::fs::create_dir_all(&cli.out_dir).ok();

    // Prepare data
    log::info!("Starting {} v{}", TOOL_NAME, TOOL_VERSION);
    let call_parameters = prepare_resource(&cli.genome, &cli.input);

    let out_json = format!("{}/{}.json", cli.out_dir, sample_id);
    let out_tsv = format!("{}/{}.tsv", cli.out_dir, sample_id);

    let (bam_name, index_name) = if cli.input.contains("##idx##") {
        let parts: Vec<&str> = cli.input.splitn(2, "##idx##").collect();
        (parts[0].to_string(), Some(parts[1].to_string()))
    } else {
        (cli.input.clone(), None)
    };

    let count_file = cli
        .count_file_path
        .as_ref()
        .map(|p| format!("{}/{}_count.txt", p, sample_id));

    if !bam_name.contains("://") && !Path::new(&bam_name).exists() {
        log::warn!("Input file for sample {} does not exist.", sample_id);
        // Write empty output files (matching Python behavior)
        std::fs::write(&out_json, "{}").unwrap();
        let header = "Sample\tGenotype\tFilter\tActivity score\tPredicted phenotype\n";
        std::fs::write(&out_tsv, header).unwrap();
        return;
    }

    let features = FeatureFlags {
        strand_bias_all: cli.strand_bias_all,
        fuzzy_match: cli.fuzzy_match,
        quality_aware: cli.quality_aware,
    };

    if features.strand_bias_all || features.fuzzy_match || features.quality_aware {
        log::info!(
            "Experimental features: strand_bias_all={}, fuzzy_match={}, quality_aware={}",
            features.strand_bias_all, features.fuzzy_match, features.quality_aware
        );
    }

    log::info!("Processing sample {}", sample_id);

    let cyp2d6_call = d6_star_caller(
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
    }

    // Write JSON
    log::info!("Writing to json ({})", out_json);
    let mut json_map = serde_json::Map::new();
    let call_json = serde_json::json!({
        "Coverage_MAD": cyp2d6_call.coverage_mad,
        "Median_depth": cyp2d6_call.median_depth,
        "Total_CN": cyp2d6_call.total_cn,
        "Spacer_CN": cyp2d6_call.spacer_cn,
        "Total_CN_raw": cyp2d6_call.total_cn_raw,
        "Spacer_CN_raw": cyp2d6_call.spacer_cn_raw,
        "Variants_called": cyp2d6_call.variants_called,
        "CNV_group": cyp2d6_call.cnv_group,
        "Genotype": cyp2d6_call.genotype,
        "Filter": cyp2d6_call.filter,
        "Raw_star_allele": cyp2d6_call.raw_star_allele,
        "Call_info": cyp2d6_call.call_info,
        "Exon9_CN": cyp2d6_call.exon9_cn,
        "CNV_consensus": cyp2d6_call.cnv_consensus,
        "d67_snp_call": cyp2d6_call.d67_snp_call,
        "d67_snp_raw": cyp2d6_call.d67_snp_raw,
        "Variant_raw_count": cyp2d6_call.variant_raw_count,
    });
    json_map.insert(sample_id.clone(), call_json);
    let json_output = serde_json::to_string(&json_map).unwrap();
    std::fs::write(&out_json, &json_output).unwrap();

    // Write TSV
    log::info!("Writing to tsv ({})", out_tsv);
    let haplotype_functionality =
        phenotype::load_haplotype_functionality(data::HAPLOTYPE_FUNC);
    let sorted_genotype = phenotype::sort_genotype(cyp2d6_call.genotype.as_deref());
    let predictions = phenotype::match_phenotype(
        sorted_genotype.as_deref(),
        &haplotype_functionality,
    );

    let count_of_diplotypes = predictions.len();
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
    tsv_content.push_str("Sample\tGenotype\tFilter\tActivity score\tPredicted phenotype\n");
    tsv_content.push_str(&format!(
        "{}\t{}\t{}\t{}\t{}\n",
        sample_id,
        genotype_str,
        cyp2d6_call.filter.as_deref().unwrap_or("None"),
        activity_scores,
        predicted_phenotypes,
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

    std::fs::write(&out_tsv, &tsv_content).unwrap();

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

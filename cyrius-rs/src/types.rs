use crate::caller::phase_disambiguate::PairEvidence;
use indexmap::IndexMap;
use std::collections::HashMap;

/// Genomic region: (chromosome, start, end, name)
pub type Region = (String, i64, i64, String);

/// Region dictionary: region_type -> Vec<(Region, gc_content)>
pub type RegionDic = IndexMap<String, Vec<(Region, f64)>>;

/// GMM parameters: variant_id -> parameter_name -> list of values
pub type GmmParameter = HashMap<String, HashMap<String, Vec<String>>>;

/// SNP lookup database
#[derive(Debug, Clone)]
pub struct SnpLookup {
    pub dsnp1: IndexMap<String, String>,
    pub dsnp2: IndexMap<String, String>,
    pub nchr: String,
    pub dindex: HashMap<String, usize>,
}

/// Copy number call result from GMM
#[derive(Debug, Clone, Copy)]
pub struct CnCall {
    pub cn: Option<u32>,
    pub depth_value: f64,
}

/// Normalized depth result
#[derive(Debug, Clone)]
pub struct NormalizedDepth {
    pub normalized: HashMap<String, Option<f64>>,
    pub mediandepth: f64,
    pub mad: f64,
    /// Fraction of mapped reads starting in REP7 that have MAPQ 0 (None when
    /// not measured, e.g. from a count file, or when REP7 has no mapped reads).
    pub rep7_mapq0_fraction: Option<f64>,
}

/// Raw copy number call (d67 + spacer)
#[derive(Debug, Clone, Copy)]
pub struct RawCnCall {
    pub d67_cn: Option<u32>,
    pub d67_depth: f64,
    pub spacer_cn: Option<u32>,
    pub spacer_depth: f64,
}

/// CN regions for CNVTAG lookup
#[derive(Debug, Clone, Copy)]
pub struct CnRegions {
    pub total_cn: u32,
    pub exon9_and_downstream: u32,
    pub exon9_to_intron1: u32,
    pub intron1_upstream: u32,
}

/// Consensus CN across genomic regions (from cnv_hybrid)
#[derive(Debug, Clone)]
pub struct CnConsensus {
    pub rep: u32,
    pub exon9_and_downstream: Option<u32>,
    pub exon9_to_intron4: Option<u32>,
    pub intron4_to_intron1: Option<u32>,
    pub intron1_upstream: Option<u32>,
}

/// Exon9 values passed to star allele matching
#[derive(Debug, Clone, Copy)]
pub struct Exon9Values {
    pub exon9_cn: Option<u32>,
    pub exon9cn_in_consensus: Option<u32>,
    pub exon9_raw_site1: f64,
    pub exon9_raw_site2: f64,
}

/// Star allele combination tables (11 tables)
#[derive(Debug, Clone)]
pub struct StarCombinations {
    pub dhap: HashMap<String, String>,
    pub dhap2: HashMap<String, Vec<String>>,
    pub dhap3: HashMap<String, Vec<String>>,
    pub dhap3pair: HashMap<String, Vec<String>>,
    pub dhap4pair: HashMap<String, Vec<String>>,
    pub dhap5pair: HashMap<String, Vec<String>>,
    pub dhap6pair: HashMap<String, Vec<String>>,
    pub dhap_exon9_x2: HashMap<String, Vec<String>>,
    pub dhap_exon9_x3: HashMap<String, Vec<String>>,
    pub dhap_exon9_x4: HashMap<String, Vec<String>>,
    pub dhap_dup_exon9: HashMap<String, Vec<String>>,
    /// Reverse map: star allele name -> variant key (e.g., "*21" -> "g.42128945C>T_g.42129033G>A")
    pub dstar: HashMap<String, String>,
}

/// Raw star allele call result
#[derive(Debug, Clone)]
pub struct RawStarCall {
    pub call_info: Option<String>,
    pub candidate: Vec<String>,
    pub star_call: Vec<String>,
}

/// Final star allele call result
#[derive(Debug, Clone)]
pub struct StarCall {
    pub call_info: Option<String>,
    pub variants_called: Option<String>,
    pub raw_call: Option<Vec<String>>,
    pub clean_call: Option<String>,
}

/// Compact, serializable summary of k-mer validation for the output.
///
/// Validation is an orthogonal QC channel that records whether an
/// alignment-free paralog ratio agrees with the depth-derived CNV tag.
/// It does NOT influence the genotype, filter, or confidence — only
/// records the verdict alongside them.
#[derive(Debug, Clone)]
pub struct KmerValidationSummary {
    pub agrees: bool,
    pub kmer_category: String,
    pub overall_ratio: f64,
    pub estimated_d6_cn: f64,
    pub estimated_d7_cn: f64,
    pub hybrid_detected: bool,
    pub total_reads: u64,
    pub matched_reads: u64,
    pub segments_summary: String,
}

/// Full D6 call result (matches Python d6_call namedtuple)
#[derive(Debug, Clone)]
pub struct D6Call {
    pub coverage_mad: f64,
    pub median_depth: f64,
    pub total_cn: Option<u32>,
    pub spacer_cn: Option<u32>,
    pub total_cn_raw: Option<f64>,
    pub spacer_cn_raw: Option<f64>,
    pub variants_called: Option<Vec<String>>,
    pub cnv_group: Option<String>,
    pub genotype: Option<String>,
    pub filter: Option<String>,
    pub raw_star_allele: Option<Vec<String>>,
    pub call_info: Option<String>,
    pub exon9_cn: Option<u32>,
    pub cnv_consensus: Option<String>,
    pub d67_snp_call: Option<String>,
    pub d67_snp_raw: Option<String>,
    pub variant_raw_count: Option<IndexMap<String, String>>,
    pub kmer_validation: Option<KmerValidationSummary>,
    /// Top posterior from the likelihood resolver (None if it didn't run).
    pub resolver_posterior: Option<f64>,
    /// Which step narrowed a `;` call: `"fragment_phase"` or `"resolver"`
    /// (None if neither did).
    pub resolved_by: Option<&'static str>,
    /// Every differentiating pair the phase step counted (None if the phase
    /// step did not run).
    pub phase_pairs: Option<Vec<PairEvidence>>,
    /// See `NormalizedDepth::rep7_mapq0_fraction`.
    pub rep7_mapq0_fraction: Option<f64>,
    /// Target variants with no dosage call (see `call_variants::uncallable_variants`).
    pub uncallable_target_sites: Option<Vec<String>>,
}

impl D6Call {
    /// Build an all-`None` no-call result carrying only the coverage stats.
    /// Used when the caller cannot proceed (e.g. depth normalization failed or
    /// the alignment file could not be read).
    pub fn no_call(coverage_mad: f64, median_depth: f64) -> Self {
        D6Call {
            coverage_mad,
            median_depth,
            total_cn: None,
            spacer_cn: None,
            total_cn_raw: None,
            spacer_cn_raw: None,
            variants_called: None,
            cnv_group: None,
            genotype: None,
            filter: None,
            raw_star_allele: None,
            call_info: None,
            exon9_cn: None,
            cnv_consensus: None,
            d67_snp_call: None,
            d67_snp_raw: None,
            variant_raw_count: None,
            kmer_validation: None,
            resolver_posterior: None,
            resolved_by: None,
            phase_pairs: None,
            rep7_mapq0_fraction: None,
            uncallable_target_sites: None,
        }
    }
}

/// Haplotype functionality info
#[derive(Debug, Clone)]
pub struct HaplotypeFunctionality {
    pub activity: String,
    pub function: String,
    pub evidence_strength: String,
    pub evidence_summary: String,
}

/// Haplotype detail for output
#[derive(Debug, Clone)]
pub struct HaplotypeDetail {
    pub haplotype: String,
    pub activity: String,
    pub function: String,
    pub evidence_strength: String,
    pub evidence_summary: String,
}

/// Prediction result per diplotype
#[derive(Debug, Clone)]
pub struct Prediction {
    pub total_activity: String,
    pub predicted_phenotype: String,
    pub haplotype_details: Vec<HaplotypeDetail>,
}

/// Resource info (all loaded parameters for calling)
#[derive(Debug, Clone)]
pub struct ResourceInfo {
    pub genome: String,
    pub gmm_parameter: GmmParameter,
    pub region_dic: RegionDic,
    pub snp_db: SnpLookup,
    pub var_db: SnpLookup,
    pub var_homo_db: SnpLookup,
    pub haplotype_db: HashMap<String, SnpLookup>,
    pub var_list: Vec<String>,
    pub star_combinations: StarCombinations,
}

/// Population frequency row
#[derive(Debug, Clone)]
pub struct FrequencyRow {
    pub biogeographic_group: String,
    pub values: IndexMap<String, String>,
}

/// Feature flags for toggling experimental improvements.
///
/// Serialized as-is into each sample's JSON record (`Flags`), so every field
/// added here is reported there with its effective value.
#[derive(Debug, Clone, Copy, serde::Serialize)]
pub struct FeatureFlags {
    /// Apply strand bias filtering to ALL variant sites (not just NOISY_VAR)
    pub strand_bias_all: bool,
    /// Use fuzzy matching when exact star allele lookup fails (on by default)
    pub fuzzy_match: bool,
    /// Use physical phasing to disambiguate multiple-match diplotypes (on by default)
    pub phase_disambiguate: bool,
    /// Use read-pair (mate-pair) phasing for longer-range disambiguation
    pub phase_readpair: bool,
    /// Use paralog ratio changepoint detection for hybrid identification
    pub changepoint_hybrid: bool,
    /// Diagnostic, log-only: report *5/*5 hemizygosity; never changes the call
    pub het_check: bool,
    /// Flag calls as uncertain when spacer CN is inconsistent with CNV group
    pub spacer_cn_check: bool,
    /// Run post-hoc consistency checks (conversion map, mismatch rate, allele balance)
    pub consistency_check: bool,
    /// Use HMM-based CNV segmentation as fallback when consensus fails
    pub hmm_cnv: bool,
    /// Detect soft-clip clusters as structural breakpoint evidence
    pub clip_evidence: bool,
    /// Compute CYP2D7 sub-region depth profile for hybrid confirmation
    pub d7_depth: bool,
    /// Use k-mer-based paralog ratio validation
    pub kmer_validation: bool,
}

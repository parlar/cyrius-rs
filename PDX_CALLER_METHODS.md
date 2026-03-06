# PDX Caller Methods Inventory

Reference of all methods in `/home/parlar_ai/dev/pdx_caller/python/pdx_caller/` and their status in cyrius-rs.

## Legend

- PORTED: Implemented in cyrius-rs
- PARTIAL: Partially implemented or simplified version exists
- NOT PORTED: Not yet in cyrius-rs
- N/A: Not applicable for cyrius-rs architecture

---

## 1. Depth & Copy Number

### depth.py — Depth-based CN estimation

| Function | What it does | Status | Notes |
|---|---|---|---|
| `collect_depth()` | Per-base depth from BAM with MAPQ + softclip tracking | PORTED | bin_count.rs |
| `estimate_copy_number()` | Gene/flanking depth ratio → CN. Checks flanking reliability (>5.0 AND >30% of gene) | PORTED | gmm.rs + bin_count.rs |
| `collect_spacer_depth()` | Spacer region CN as tandem indicator | PORTED | GMM spacer call in main.rs |
| `triage_sample()` | Fast-path Z-score triage against cohort baseline | NOT PORTED | Needs cohort statistics |
| GC bias correction | LOWESS-based GC content correction | NOT PORTED | cyrius-rs uses raw depth ratios |

### gmm.rs (cyrius-rs) — Gaussian Mixture Model

| Check | Threshold | Notes |
|---|---|---|
| Posterior cutoff | 0.95 | Returns None if below |
| P-value cutoff | 1e-3 | Gaussian p-value check |
| High CN fallback | depth > 7.5 | Rounds depth, sets LowQ_high_CN filter |

---

## 2. Structural Determination

### structural.py — Allele-agnostic structural calling (Step 1 of v2)

| Function | What it does | Status | Notes |
|---|---|---|---|
| `compute_cn_profile()` | Per-position D6 CN from ratio × total_CN | PORTED | get_fraction() + call_cn_snp() |
| `_call_site_cn()` | Poisson discrete CN at each diagnostic site (posterior ≥ 0.9) | PORTED | copy_number_call.rs call_reg1_cn() |
| `_detect_cn_changepoint()` | Consensus CN per gene region, detect step changes | PORTED | cnv_hybrid.rs get_cnvtag() |
| `determine_structure()` | Combines CN + het + spacer → structural config with confidence (0-1) | PARTIAL | cnv_hybrid.rs lacks confidence scores |
| `identify_hybrid_alleles()` | Maps breakpoint + direction → candidate hybrid alleles | NOT PORTED | Hardcoded in match_star_allele.rs |

Confidence scores in structural.py:
- Homozygous deletion: 0.95 (median_cn < 0.3) / 0.70
- Hemizygous deletion: 0.95 (with het agreement) / 0.85 / 0.60
- Normal diploid: 0.95 (1.7 ≤ cn ≤ 2.3) / 0.75
- Duplication: 0.85 (2.5 ≤ cn ≤ 3.5) / 0.70
- Hybrid: 0.70-0.85 depending on type

---

## 3. Paralog Ratio & Changepoint Detection

### paralog_ratio.py — CBS on D6/(D6+D7) ratio

| Function | What it does | Status | Notes |
|---|---|---|---|
| `compute_paralog_ratios()` | Per-position R = D6/(D6+D7), masks low-depth | PORTED | changepoint.rs compute_ratios() |
| `cbs_segment()` | Recursive binary segmentation via Welch t-test | PORTED | changepoint.rs cbs_recurse() |
| `detect_changepoints()` | Two-pass CBS (t=3.0, retry at 2.0 if imbalanced) | PORTED | changepoint.rs detect_changepoints() |
| `_estimate_cn_category()` | Coarse CN from overall ratio | NOT PORTED | |

Parameters: t_thresh=3.0, min_bins=5, min_depth=5, min_ratio_diff=0.15

---

## 4. Gene Conversion Mapping

### gene_conversion.py — Spatial D6/D7 classification (~250 lines)

| Function | What it does | Status | Notes |
|---|---|---|---|
| `build_conversion_map()` | Sliding window (size=5) classifies regions as D6/D7/ambiguous | NOT PORTED | Uses snp_d6/snp_d7 (available) |
| `compute_conversion_likelihood()` | Scores diplotypes against observed D6/D7 pattern | NOT PORTED | Bayesian component |

Classification thresholds:
- D6-like: ratio > 0.35
- D7-like: ratio < 0.15
- Ambiguous: in between
- `is_chimeric = True` if both D6 and D7 segments present

Hardcoded breakpoints for known hybrids:
- *68: ~900bp from gene start
- *36: ~3800bp from gene start

Could be ported as a consistency check: if cn2 but chimeric pattern → flag uncertain.

---

## 5. Variant Calling & Filtering

### call_variants.rs (cyrius-rs) — Already ported

| Check | Details |
|---|---|
| Strand bias (NOISY_VAR) | Fisher exact, p < 0.05. Forward ≤1 OR reverse ≤1 |
| Strand bias all (feature) | Extended to all non-NOISY sites, p < 0.02 |
| Min read support | CLEAN: 2, NOISY: 7, other: 4 |
| Exon9 GC consensus | All sites agree + posterior ≥ 0.88 |
| Haplotype variant checks | g.42126938, g.42127526/556, g.42127803, g.42130655insA |

### inference.py — Allele prefiltering (NOT PORTED)

| Function | What it does | Status |
|---|---|---|
| `prefilter_alleles()` | Contradiction detection: expects alt but VAF<5%, or expects ref but VAF>95% | NOT PORTED |

- mismatch_rate = n_contradictions / n_checked
- Threshold: ≤0.30 (or ≤0.15 with known structure)
- Always keeps: reference alleles, deletions, hybrids
- **Could be used as post-hoc consistency check on final call**

---

## 6. Bayesian Inference

### inference.py — 13 evidence streams (~1800 lines)

| Evidence Stream | Function | Status | Feasibility |
|---|---|---|---|
| K-mer counts | `compute_kmer_likelihood()` | NOT PORTED | Needs k-mer infrastructure |
| Depth | `compute_depth_likelihood()` | PARTIAL | GMM call exists |
| Changepoint | `compute_changepoint_likelihood()` | PARTIAL | CBS exists, no likelihood |
| Read-pair crossing | `compute_readpair_likelihood()` | NOT PORTED | Medium effort |
| SNV VAF | `compute_snv_likelihood()` | NOT PORTED | **Easy — data available** |
| Spacer CN | `compute_spacer_likelihood()` | PARTIAL | Consistency check exists |
| Phasing | `compute_phasing_likelihood()` | PARTIAL | phase_disambiguate exists |
| Pangenome | `compute_pangenome_likelihood()` | NOT PORTED | Needs vg binary |
| Het count | `compute_het_count_likelihood()` | PARTIAL | het_check exists |
| Allele balance | `compute_allele_balance_likelihood()` | NOT PORTED | **Easy — data available** |
| CN phasing | `compute_cn_phasing_likelihood()` | NOT PORTED | Medium effort |
| Realignment | `compute_realign_allele_likelihood()` | NOT PORTED | Needs haplotype FASTA |
| Gene conversion | `compute_conversion_likelihood()` | NOT PORTED | **Easy — data available** |

### VCF FILTER thresholds (vcf_output.py):
- PASS: posterior ≥ 0.95
- Low_posterior: 0.80 ≤ posterior < 0.95
- Low_confidence: posterior < 0.80

---

## 7. Read-Pair Constraint Graph

### read_pairs.py (~520 lines)

| Function | What it does | Status |
|---|---|---|
| `build_read_pair_graph()` | Count consistent vs crossing mate pairs at diagnostic positions | NOT PORTED |
| `crossing_by_position()` | Per-position crossing fraction (reveals breakpoints) | NOT PORTED |
| `compute_readpair_likelihood()` | Beta-binomial scoring: normal=0.02, hybrid=0.15 crossing | NOT PORTED |

ReadPairEdge: pos1_idx, pos2_idx, paralog assignment, is_consistent flag
ReadPairGraph: edges, n_consistent, n_crossing, crossing_fraction

**Could validate hybrid calls**: if called as hybrid but crossing_fraction ~0.02 → flag uncertain.

---

## 8. Copy-Number Phasing

### cn_phasing.py (~200 lines)

| Function | What it does | Status |
|---|---|---|
| `compute_haplotype_cn()` | Per-haplotype CN breakdown | NOT PORTED |
| `_count_alt_copies()` | Count gene copies carrying specific alt | NOT PORTED |
| `_collect_informative_positions()` | SNVs where haplotypes differ in dosage | NOT PORTED |
| `compute_cn_phasing_likelihood()` | Beta-binomial scoring of expected vs observed VAF per CN assignment | NOT PORTED |

Example: CN=3, *1/*2x2. At *2-defining positions:
- *2 duplicated → ~2/3 alt fraction
- *1 duplicated → ~1/3 alt fraction
Min depth: 15, max contribution cap: 15.0, rho=0.02

---

## 9. De Novo Assembly

### assembly.py (~800 lines)

| Function | What it does | Status |
|---|---|---|
| `extract_gene_reads()` | Paralog-resolve reads via diagnostic k-mer voting | NOT PORTED |
| `build_debruijn_graph()` | K-mer graph (k=31, min_count=3) | NOT PORTED |
| `find_heaviest_path()` | Greedy contig assembly | NOT PORTED |
| `align_contig_to_reference()` | Needleman-Wunsch (match=+2, mismatch=-1, gap=-2) | NOT PORTED |
| `find_novel_variants()` | Jaccard similarity to known alleles, report novel variants | NOT PORTED |

---

## 10. Haplotype Sequence Construction

### haplotypes.py

| Function | What it does | Status |
|---|---|---|
| `build_haplotype()` | Apply defining SNVs to gene reference | NOT PORTED |
| `_build_hybrid_haplotype()` | Chimeric D6/D7 sequence at breakpoint | NOT PORTED |
| `build_all_haplotypes()` | All known allele sequences | NOT PORTED |
| `write_haplotype_fasta()` | Output for Rust realigner (pdx_extract) | NOT PORTED |

---

## 11. Statistical Phasing

### phasing.py

| Function | What it does | Status |
|---|---|---|
| `read_phased_vcf()` | Parse WhatsHap/SHAPEIT/Eagle VCF | NOT PORTED |
| `extract_star_allele_phase_constraints()` | Cross-ref phased SNVs with defining positions | NOT PORTED |
| `compute_phasing_likelihood()` | Score diplotypes against phase blocks | NOT PORTED |

---

## 12. Pangenome Graph Alignment

### pangenome.py

| Function | What it does | Status |
|---|---|---|
| `run_vg_giraffe()` | Extract reads → FASTQ → vg giraffe → GAF | NOT PORTED |
| `parse_gaf()` | Parse Graph Alignment Format | NOT PORTED |
| `extract_pangenome_evidence()` | Per-allele read support from graph paths | NOT PORTED |

Requires: vg binary, prebuilt .giraffe.gbz graph

---

## 13. Novel Allele Function Prediction

### function_prediction.py (~670 lines)

| Function | What it does | Status |
|---|---|---|
| `predict_function_rule_based()` | Priority rules: frameshift→NONE, splice→NONE, missense@critical→NONE, SRS/heme→DECREASED, synonymous→NORMAL | NOT PORTED |
| `predict_function_nearest_neighbor()` | Jaccard similarity to known alleles | NOT PORTED |
| `predict_novel_allele_function()` | Combined prediction with confidence levels | NOT PORTED |

Critical residues: {443, 309, 301, 216, 120}

---

## 14. Calibration

### calibration.py

| Function | What it does | Status |
|---|---|---|
| `estimate_depth_noise()` | Std-dev of CN residuals per class | NOT PORTED |
| `estimate_kmer_overdispersion()` | Beta-binomial MLE on k-mer counts | NOT PORTED |
| `estimate_seq_error_rate()` | Median minor allele fraction at hom sites | NOT PORTED |
| `estimate_changepoint_params()` | Ratio noise per structural class | NOT PORTED |
| `estimate_readpair_params()` | Crossing fraction per class | NOT PORTED |
| `calibrate_parameters()` | Orchestrates all estimations | NOT PORTED |

Default parameters:
- depth_std_cn2: 0.15
- spacer_std: 0.30
- seq_error_rate: 0.01
- kmer_rho: 0.02, snv_rho: 0.03
- cbs_t_threshold: 3.0
- changepoint_sigma: normal=0.10, hybrid=0.15
- expected_crossing: normal=0.02, hybrid=0.15
- hybrid_gene_fraction: 0.30

---

## 15. Rust Binaries in pdx_caller

| Binary | Purpose | Status in cyrius-rs |
|---|---|---|
| `pdx_kmer` | Fast k-mer counting | NOT PORTED |
| `pdx_extract` | Per-allele read realignment | NOT PORTED |
| `pdx_graph` | Graph operations (unclear) | NOT PORTED |

---

## Quick-Win Consistency Checks (use existing cyrius-rs data)

### 1. Allele Mismatch Rate
- Data: var_alt, var_ref, star allele definitions
- Logic: For each defining SNV of called allele, check VAF. If expects alt but VAF<5% → contradiction.
- Flag if mismatch_rate > 0.30
- Source: inference.py prefilter_alleles()

### 2. Allele Balance / VAF Consistency
- Data: var_alt, var_ref
- Logic: For cn2 het sites, VAF should be ~0.5. Compute mean VAF deviation from expected.
- Flag if mean deviation > threshold
- Source: inference.py compute_allele_balance_likelihood()

### 3. Gene Conversion Map Consistency
- Data: snp_d6, snp_d7 (already computed)
- Logic: Sliding window D6/D7 classification. If cn2 but chimeric → flag. If hybrid but uniform → flag.
- Source: gene_conversion.py build_conversion_map()

### 4. Read-Pair Crossing Fraction
- Data: BAM reads (new fetch needed)
- Logic: Count crossing pairs. Normal ~0.02, hybrid ~0.15.
- Validates hybrid vs normal calls.
- Source: read_pairs.py build_read_pair_graph()

### 5. Spacer CN Consistency (IMPLEMENTED)
- Data: total_cn, spacer_cn, cnvtag
- Logic: total_cn - spacer_cn must match expected for CNV group
- Feature flag: --spacer-cn-check

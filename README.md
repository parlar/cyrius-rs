# Cyrius-RS

A Rust implementation of [Cyrius](https://github.com/Illumina/Cyrius), a tool for genotyping CYP2D6 from whole-genome sequencing (WGS) data. Cyrius-RS is a faithful port of the Python version with additional experimental modules for structural variant detection, compiled into a single static binary with no runtime dependencies beyond htslib.

## Benchmark

Evaluated on 53 GeT-RM samples with truth annotations from the 1000 Genomes high-coverage dataset:

| Metric | Value |
|---|---|
| **Exact match** | **49 / 53 (92.5%)** |
| Ambiguous (truth included) | 4 |
| Wrong | 0 |

The 4 ambiguous calls report multiple alternative diplotypes (separated by `;`) because the data does not uniquely determine a single answer. In each case, the correct diplotype is among the reported alternatives. In 3 of 4 cases, the alternatives have identical clinical effect (same activity score and metabolizer phenotype). One case (NA19908) has a clinically relevant difference:

- **HG00589, NA18973**: `*1/*21;*2/*21` — cannot distinguish \*1 from \*2. **Same clinical effect** (activity 1.0, Intermediate Metabolizer)
- **NA18565**: `*36/*36+*10;*10/*36x2` — tandem arrangement ambiguity. **Same clinical effect** (activity 0.25, Poor Metabolizer)
- **NA19908**: `*1/*46;*43/*45` — two different diplotype interpretations. **Different clinical effect** (\*1/\*46 = activity 2.0, Normal Metabolizer vs \*43/\*45 = activity 1.0, Intermediate Metabolizer)

Each call includes a confidence score (0–1) with a qualitative label (HIGH, MEDIUM, LOW). On the 53 truth samples: 40 HIGH, 13 MEDIUM, 0 LOW.

<details>
<summary>Full results (53 samples)</summary>

| Sample | Genotype | Confidence | Truth |
|---|---|---|---|
| HG00111 | \*3/\*3 | 0.96 HIGH | \*3/\*3 |
| HG00337 | \*2x2/\*22 | 0.88 HIGH | \*2x2/\*22 |
| HG00373 | \*2/\*2 | 0.93 HIGH | \*2/\*2 |
| HG00421 | \*2/\*10x2 | 0.92 HIGH | \*2/\*10x2 |
| HG00436 | \*2x2/\*71 | 0.98 HIGH | \*2x2/\*71 |
| HG00463 | \*36+\*10/\*36+\*10 | 0.92 HIGH | \*36+\*10/\*36+\*10 |
| HG00589 | \*1/\*21;\*2/\*21 | 0.77 MEDIUM | \*1/\*21 |
| HG01086 | \*31/\*178 | 0.94 HIGH | \*31/\*178 |
| HG01094 | \*1/\*31 | 0.95 HIGH | \*1/\*31 |
| HG01108 | \*2/\*106 | 0.95 HIGH | \*2/\*106 |
| HG01190 | \*5/\*68+\*4 | 0.61 MEDIUM | \*5/\*68+\*4 |
| HG01680 | \*28/\*59 | 0.89 HIGH | \*28/\*59 |
| HG02373 | \*14/\*36+\*10 | 0.95 HIGH | \*14/\*36+\*10 |
| HG03225 | \*5/\*56 | 0.91 HIGH | \*5/\*56 |
| HG03246 | \*5/\*43 | 0.95 HIGH | \*5/\*43 |
| HG03259 | \*5/\*106 | 0.97 HIGH | \*5/\*106 |
| HG03619 | \*2/\*113 | 0.97 HIGH | \*2/\*113 |
| HG03643 | \*2/\*7 | 0.90 HIGH | \*2/\*7 |
| HG03703 | \*1/\*99 | 0.92 HIGH | \*1/\*99 |
| HG03780 | \*1/\*112 | 0.95 HIGH | \*1/\*112 |
| HG03781 | \*2/\*99 | 0.94 HIGH | \*2/\*99 |
| HG03882 | \*1/\*112 | 0.88 HIGH | \*1/\*112 |
| HG04206 | \*2/\*113 | 0.96 HIGH | \*2/\*113 |
| NA06989 | \*9/\*9 | 0.93 HIGH | \*9/\*9 |
| NA11832 | \*1/\*68+\*4 | 0.94 HIGH | \*1/\*68+\*4 |
| NA12154 | \*33/\*68+\*4 | 0.85 HIGH | \*33/\*68+\*4 |
| NA12878 | \*3/\*68+\*4 | 0.91 HIGH | \*3/\*68+\*4 |
| NA18526 | \*1/\*36+\*36+\*10 | 0.74 MEDIUM | \*1/\*36+\*36+\*10 |
| NA18544 | \*10/\*41 | 0.80 MEDIUM | \*10/\*41 |
| NA18545 | \*5/\*36x2+\*10x2 | 0.96 HIGH | \*5/\*36x2+\*10x2 |
| NA18563 | \*1/\*36+\*10 | 0.79 MEDIUM | \*1/\*36+\*10 |
| NA18564 | \*2/\*36+\*10 | 0.80 MEDIUM | \*2/\*36+\*10 |
| NA18565 | \*36/\*36+\*10;\*10/\*36x2 | 0.87 HIGH | \*10/\*36x2 |
| NA18572 | \*36+\*10/\*41 | 0.80 MEDIUM | \*36+\*10/\*41 |
| NA18617 | \*36+\*10/\*36+\*10 | 0.93 HIGH | \*36+\*10/\*36+\*10 |
| NA18632 | \*36+\*36+\*10/\*52 | 0.85 HIGH | \*36+\*36+\*10/\*52 |
| NA18642 | \*1+\*90/\*36+\*10 | 0.94 HIGH | \*1+\*90/\*36+\*10 |
| NA18959 | \*2/\*36+\*10 | 0.79 MEDIUM | \*2/\*36+\*10 |
| NA18973 | \*1/\*21;\*2/\*21 | 0.80 MEDIUM | \*1/\*21 |
| NA18980 | \*2/\*36+\*10 | 0.80 MEDIUM | \*2/\*36+\*10 |
| NA19143 | \*10/\*45 | 0.95 HIGH | \*2(\*45)/\*10 |
| NA19207 | \*2x2/\*10 | 0.93 HIGH | \*2x2/\*10 |
| NA19317 | \*5/\*5 | 1.00 HIGH | \*5/\*5 |
| NA19777 | \*1/\*82 | 0.96 HIGH | \*1/\*82 |
| NA19785 | \*1/\*13+\*2 | 0.95 HIGH | \*1/\*13+\*2 |
| NA19819 | \*2/\*4x2 | 0.74 MEDIUM | \*2/\*4x2 |
| NA19908 | \*1/\*46;\*43/\*45 | 0.81 MEDIUM | \*1/\*46 |
| NA19917 | \*1/\*40 | 0.96 HIGH | \*1/\*40 |
| NA19920 | \*1/\*4x2 | 0.73 MEDIUM | \*1/\*4x2 |
| NA20289 | \*6/\*11 | 0.94 HIGH | \*6/\*11 |
| NA20803 | \*2/\*22 | 0.91 HIGH | \*2/\*22 |
| NA20875 | \*1/\*111 | 0.89 HIGH | \*1/\*111 |
| NA21105 | \*3/\*111 | 0.96 HIGH | \*3/\*111 |

</details>

To run the benchmark:

```bash
bash benchmark/run_benchmark.sh --strand-bias-all --fuzzy-match --changepoint-hybrid --het-check
```

## How It Works

The caller follows a multi-stage pipeline:

### 1. Depth normalization and copy number calling

Reads are counted across the CYP2D6/CYP2D7 region and normalized against control regions on the same chromosome. A Gaussian mixture model (GMM) assigns the total D6+D7 copy number (normal diploid = 4) and spacer copy number.

### 2. D6/D7 SNP ratio analysis

117 diagnostic positions that differ between CYP2D6 and CYP2D7 are examined. At each site, the fraction of reads supporting the D6 vs D7 allele is used to estimate the per-site D6 copy number. A consensus across all sites determines the CNV structure (e.g., deletion, duplication, exon9 hybrid/gene conversion).

### 3. Variant calling

Diagnostic variants that define star alleles are called from read pileups. This includes SNVs, small indels, and the exon9 gene conversion marker. Strand bias filtering removes artifacts at noisy positions.

### 4. Star allele matching

The observed variant pattern is matched against a precomputed table of all possible star allele combinations for the given CNV group. When no exact match is found, fuzzy matching with edit-distance scoring is used as a fallback.

### 5. Phenotype prediction

Activity scores are assigned per CPIC guidelines, and the predicted metabolizer phenotype (Poor, Intermediate, Normal, Ultrarapid) is reported.

## Features

- Genotypes CYP2D6 star alleles from BAM/CRAM files aligned to GRCh38
- Detects copy number variants (deletions, duplications, hybrids, gene conversions)
- Reports activity scores and predicted metabolizer phenotypes based on CPIC guidelines
- Optionally outputs population frequencies and haplotype functional annotations
- Single binary, no Python/pip required

## Building

Requires Rust 1.70+ and htslib development libraries (for BAM/CRAM support).

```bash
cd cyrius-rs
cargo build --release
```

The binary will be at `cyrius-rs/target/release/cyrius-rs`.

## Usage

```bash
cyrius-rs -i sample.bam -o results/
```

### Core Parameters

| Parameter | Description |
|---|---|
| `-i, --input` | Input BAM or CRAM file (indexed) |
| `-o, --outDir` | Output directory for results |
| `-r, --reference` | Reference FASTA (required for CRAM) |
| `-g, --genome` | Contig naming: `autodetect` (default), `38`, or `chr38` |
| `--id` | Sample ID (default: derived from filename) |
| `-t, --threads` | Number of threads (default: 1) |
| `--print` | Print results to stdout |
| `--countFilePath` | Use precomputed count files instead of BAM |

### Output Options

| Parameter | Description |
|---|---|
| `--population-info` | Include population frequencies in output |
| `--haplotype-info` | Include haplotype functional annotations |

### Recommended Flags

These flags are recommended for best accuracy and are used in the benchmark:

| Flag | Description |
|---|---|
| `--strand-bias-all` | Apply strand bias filtering to all variant sites, not just known noisy ones |
| `--fuzzy-match` | Use edit-distance scoring when no exact star allele match is found |
| `--changepoint-hybrid` | Use changepoint detection as a fallback for CNV group assignment |
| `--het-check` | Detect hemizygosity (zero heterozygous sites indicates a deletion) |

### Experimental Flags

These flags enable additional analysis modules that are under development:

| Flag | Description |
|---|---|
| `--read-phasing` | Score diplotype calls against read-level phasing patterns |
| `--cn-classifier` | Run k-NN classifier on per-region CN profiles for CNV group validation |
| `--quality-aware` | Weight variant calls by base quality scores |
| `--clip-evidence` | Detect soft-clip clusters as structural breakpoint evidence |
| `--hmm-cnv` | Use HMM-based CNV segmentation as fallback when consensus fails |
| `--read-voting` | Run read-level allele voting for independent QC |
| `--diplotype-caller` | Run likelihood-based diplotype caller (D6+D7 mixture model) |
| `--spacer-cn-check` | Flag calls where spacer copy number is inconsistent with CNV group |
| `--consistency-check` | Run post-hoc consistency checks (conversion map, mismatch rate, allele balance) |

## Output

Results are written as both JSON and TSV files. The TSV contains:

```
Sample       Genotype   Filter   Confidence     Activity score   Predicted phenotype
sample123    *1/*4      PASS     0.96 (HIGH)    1.0              Intermediate Metabolizer
```

The confidence score (0–1) is a weighted combination of six quality signals: depth quality, CN rounding quality, match quality, variant completeness, variant specificity, and D6/D7 SNP ratio consistency. Labels: HIGH (≥0.85), MEDIUM (0.50–0.84), LOW (<0.50).

The JSON output includes additional fields: total copy number, spacer copy number, raw depth values, CNV group, called variants, confidence score, and the full variant count table.

## Project Structure

```
cyrius-rs/src/
├── main.rs                  # CLI and main pipeline
├── lib.rs
├── types.rs                 # Shared types (D6Call, FeatureFlags, etc.)
├── data/                    # Embedded data files (SNP definitions, star allele tables)
├── stats.rs                 # Statistical utilities
├── fisher.rs                # Fisher's exact test
├── phenotype.rs             # Activity scores and metabolizer phenotype
├── depth_calling/           # Depth-based copy number analysis
│   ├── bin_count.rs         #   Read counting and normalization
│   ├── gmm.rs               #   Gaussian mixture model
│   └── snp_count.rs         #   D6/D7 SNP ratio counting
├── caller/                  # Star allele calling logic
│   ├── call_variants.rs     #   Variant calling from pileups
│   ├── cnv_hybrid.rs        #   CNV group assignment from SNP ratios
│   ├── match_star_allele.rs #   Star allele pattern matching
│   ├── construct_star_table.rs # Star allele combination tables
│   ├── confidence.rs        #   Multi-signal confidence scoring
│   ├── fuzzy_match.rs       #   Edit-distance fallback matching
│   ├── strand_bias_all.rs   #   Strand bias filtering
│   ├── changepoint.rs       #   Changepoint-based CNV detection
│   ├── het_check.rs         #   Hemizygosity detection
│   ├── cn_classifier.rs     #   k-NN CN profile classifier
│   ├── read_phasing.rs      #   Read-level phasing constraints
│   ├── clip_evidence.rs     #   Soft-clip cluster detection
│   ├── hmm_cnv.rs           #   HMM-based CNV segmentation
│   ├── read_voting.rs       #   Read-level allele voting
│   ├── diplotype_caller.rs  #   Likelihood-based diplotype caller
│   └── phase_disambiguate.rs #  Phase-based disambiguation
└── align/                   # Local sequence alignment engine
    ├── mod.rs               #   LocalIndex (q-gram + Myers bit-parallel)
    ├── myers.rs             #   Myers 1999 semi-global edit distance
    ├── pigeonhole.rs        #   Pigeonhole seed filter
    ├── qgram_index.rs       #   Q-gram index
    └── sequence.rs          #   Base encoding, reverse complement
```

## License

Apache-2.0

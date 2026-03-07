# Cyrius-RS

A Rust implementation of [Cyrius](https://github.com/Illumina/Cyrius), a tool for genotyping CYP2D6 from whole-genome sequencing (WGS) data. Cyrius-RS is a faithful port of the Python version with additional experimental modules for structural variant detection, compiled into a single static binary with no runtime dependencies beyond htslib.

## Benchmark

Evaluated on 62 GeT-RM samples (53 with truth annotations) from the 1000 Genomes high-coverage dataset:

| Metric | Value |
|---|---|
| Exact match | 32 / 53 |
| Order-only difference | 19 / 53 |
| **Total correct** | **51 / 53 (96.2%)** |
| Discordant | 2 |

The 2 remaining discrepancies:

- **NA18565**: Called `*36/*36+*10`, truth `*10/*36x2` — a phasing ambiguity. Both represent the same allele composition (1×\*10 + 2×\*36) but differ in which alleles are in tandem. The Python Cyrius produces the same call. Resolving this requires long-read phasing data.
- **NA19143**: Called `*10/*45`, truth `*2(*45)/*10` — a partial duplication where the depth signal is borderline for CN detection.

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
Sample       Genotype   Filter   Activity score   Predicted phenotype
sample123    *1/*4      PASS     1.0              Intermediate Metabolizer
```

The JSON output includes additional fields: total copy number, spacer copy number, raw depth values, CNV group, called variants, and the full variant count table.

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
│   ├── fuzzy_match.rs       #   Edit-distance fallback matching
│   ├── strand_bias_all.rs   #   Strand bias filtering
│   ├── changepoint.rs       #   Changepoint-based CNV detection
│   ├── het_check.rs         #   Hemizygosity detection
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

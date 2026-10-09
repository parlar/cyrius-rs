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

## Improvements over bcyrius

cyrius-rs is a faithful port of bcyrius — same module boundaries, same data resources, same depth/GMM/SNP/matching flow — with a small number of targeted changes that are active by default:

- **Exon 9 GC threshold relaxed (0.75 → 0.70)** in `call_exon9gc` ([`call_variants.rs:323`](cyrius-rs/src/caller/call_variants.rs#L323)). bcyrius's stricter cutoff misses rare `*10D`-style duplication breakpoints downstream of exon 9 (e.g. HG00421, where the exon 9 conversion probability is 0.729). The relaxed threshold lets the exon 9 signal participate in CNV-tag correction.
- **Spacer CN as a structural discriminator** in star allele matching ([`match_star_allele.rs`](cyrius-rs/src/caller/match_star_allele.rs)). bcyrius uses the spacer region only for depth normalization. cyrius-rs uses the spacer CN as a feature to discriminate tandem and deletion-plus-duplication configurations that produce indistinguishable gene-body signals (e.g. NA18545, `exon9hyb_star5` cases).
- **Tandem notation normalized** (`*36+*36` → `*36x2`) to match Pratt 2019 truth conventions.
- **Ambiguous calls reported, not hidden** — when short reads cannot uniquely determine phase or structure, alternatives are reported as `dip1;dip2` rather than picking one. New JSON fields `Phasing_ambiguous`, `Clinical_equivalence`, and `Long_read_recommended` distinguish clinically-equivalent ambiguity (no action) from clinically-different ambiguity (long-read reflex recommended).
- **Confidence scoring** combines six signals (depth quality, CN rounding quality, match quality, variant completeness, variant specificity, D6/D7 SNP consistency) into a single score with HIGH/MEDIUM/LOW labels.
- **Changepoint-based CNV-tag fallback** ([`changepoint.rs`](cyrius-rs/src/caller/changepoint.rs)). When the standard region-consensus rule cannot assign a CNV tag, CBS with Welch's t-statistic on the D6/(D6+D7) ratio profile is used as an alternative path. Enabled by `--changepoint-hybrid`, which is included in the recommended flag set.
- **Several robustness fixes** layered on faithful port: underflow guard in `cnv_hybrid` (`total_cn − spacer_cn`), insertion-order preservation in `get_consensus` (matches Python `Counter` tie-break), missing-file handling that returns empty output instead of panicking on autodetect, and clinical-equivalence comparison that includes both activity score and phenotype.

A larger collection of opt-in flags adds independent QC channels and prototype analyses — see [Feature flags](#feature-flags) below.

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

<a id="feature-flags"></a>
### Feature flags

Flags are grouped by how they affect the output. Pass them after the core options.

#### Tier 1 — Recommended (affect the call)

These flags change the genotype, filter, or CNV tag. They are the default benchmark configuration:

```bash
--strand-bias-all --fuzzy-match --changepoint-hybrid --het-check
```

| Flag | Effect on output |
|---|---|
| `--strand-bias-all` | Strand-bias filtering on all variant sites (not just the small known-noisy set). Removes more variant-calling artifacts. |
| `--fuzzy-match` | Edit-distance fallback when no exact star allele match is found. Recovers near-matches at the cost of a `Fuzzy_match` filter label. **On by default**; the flag is accepted for compatibility and `--no-fuzzy-match` turns it off. |
| `--phase-disambiguate` | Use read-level phasing to resolve ambiguous CN=2 multi-match calls into a single diplotype. **On by default**; the flag is accepted for compatibility and `--no-phase-disambiguate` turns it off. |
| `--changepoint-hybrid` | Circular binary segmentation on the D6/(D6+D7) ratio profile, used as a CNV-tag fallback when the standard region-consensus rule produces no tag. |
| `--het-check` | Diagnostic logging of D6/D7 hemizygosity (zero het sites). Currently QC-only — kept in the recommended set because it is a useful log signal during review. |

#### Tier 2 — Additional call-affecting flags (off by default)

These also affect the genotype or filter, but are off by default because they are either situation-specific or have not yet been validated across the full benchmark:

| Flag | Effect on output |
|---|---|
| `--hmm-cnv` | HMM segmentation on the per-position CN profile, used as a second CNV-tag fallback after changepoint detection. |
| `--phase-readpair` | Extend phasing to mate-pair distance (no effect with `--no-phase-disambiguate`). |
| `--spacer-cn-check` | Flag calls with `Spacer_CN_inconsistent` when the observed spacer-vs-total-CN difference does not match the CNV group. |
| `--consistency-check` | Post-hoc validation; sets the filter to `Conversion_map_inconsistent`, `Mismatch_rate_high`, etc. if checks fail. |

#### Tier 3 — QC and validation (do not affect the call)

These flags add **independent observations** to the output without modifying the genotype, filter, or confidence. Useful for review, audit, and downstream analytics:

| Flag | What it adds |
|---|---|
| `--kmer-validation` | Alignment-free paralog-dosage check (189 diagnostic 31-mers). Emits a `Kmer_validation` JSON block with AGREE/DISAGREE. |
| `--clip-evidence` | Soft-clip cluster scan in the D6+D7+spacer region; flags candidate structural breakpoints. |
| `--d7-depth` | Per-region CYP2D7 depth profile; logs elevation/depletion signals consistent with `*36` and `*13` hybrids. |

## Output

Results are written as both JSON and TSV files. The TSV contains:

```
Sample       Genotype   Filter   Confidence     Activity score   Predicted phenotype   Phasing_ambiguous   Long_read_recommended
sample123    *1/*4      PASS     0.96 (HIGH)    1.0              Intermediate Metabolizer   false           false
```

The confidence score (0–1) is a weighted combination of six quality signals: depth quality, CN rounding quality, match quality, variant completeness, variant specificity, and D6/D7 SNP ratio consistency. Labels: HIGH (≥0.85), MEDIUM (0.50–0.84), LOW (<0.50).

### Phasing ambiguity and long-read reflex

When the data does not uniquely determine a single diplotype, the output contains multiple alternatives separated by `;`. Two additional fields help clinical labs decide when orthogonal confirmation is needed:

- **Phasing_ambiguous**: `true` when multiple alternative diplotypes are reported
- **Long_read_recommended**: `true` only when the alternatives have **different clinical effects** (different activity score or metabolizer phenotype). When all alternatives are clinically equivalent, long-read confirmation would not change the clinical interpretation.

The JSON output additionally includes `Clinical_equivalence` (`same_effect` or `different_effect`) for ambiguous calls.

The JSON output also includes: total copy number, spacer copy number, raw depth values, CNV group, called variants, confidence score, and the full variant count table. Every record carries its provenance: `Tool_version`, `Tool_commit` (the git commit the binary was built from, `-dirty` if the crate had uncommitted changes) and `Flags` (every feature flag with its effective value). `--version` prints the same version and commit.

When `--kmer-validation` is enabled, the JSON additionally contains a `Kmer_validation` object summarizing an alignment-free paralog-dosage check:

```json
"Kmer_validation": {
  "agrees": true,
  "kmer_category": "normal",
  "overall_ratio": 0.524,
  "estimated_d6_cn": 2.09,
  "estimated_d7_cn": 1.91,
  "hybrid_detected": false,
  "matched_reads": 2098,
  "total_reads": 6055,
  "segments_summary": "seg0[0..83]: ratio=0.485 n=80, ..."
}
```

This is an **independent QC channel** — it does not modify the genotype, filter, or confidence. `agrees: false` means the D6 copy number estimated from 189 diagnostic 31-mers disagrees with the dosage expected from the depth-based CNV tag. See `cyrius-rs/src/caller/kmer_validation.rs` for the underlying method.

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
│   ├── clip_evidence.rs     #   Soft-clip cluster detection
│   ├── hmm_cnv.rs           #   HMM-based CNV segmentation
│   ├── diplotype_caller.rs  #   Likelihood-based diplotype caller
│   └── phase_disambiguate.rs #  Phase-based disambiguation
└── align/                   # Local sequence alignment engine
    ├── mod.rs               #   LocalIndex (q-gram + Myers bit-parallel)
    ├── myers.rs             #   Myers 1999 semi-global edit distance
    ├── pigeonhole.rs        #   Pigeonhole seed filter
    ├── qgram_index.rs       #   Q-gram index
    └── sequence.rs          #   Base encoding, reverse complement
```

## Future work: long-read reflex for phasing ambiguity

Some CYP2D6 diplotype calls are inherently ambiguous from short-read WGS because the data cannot phase distant variants onto haplotypes. In our benchmark, 4 of 53 samples report multiple alternative diplotypes. In 3 cases the alternatives are clinically equivalent (same activity score and metabolizer phenotype). In 1 case (~2%) the ambiguity is clinically relevant and requires orthogonal confirmation.

**Proposed strategy: nanopore amplicon sequencing as a targeted reflex test.**

CYP2D6 is only ~4.4 kb, well within a single nanopore read. A long-range PCR amplicon with one primer anchored in a CYP2D6-unique flanking region ensures specificity against CYP2D7 (97% sequence homology) without requiring both primers to be paralog-specific. Each nanopore read then spans all exonic variants simultaneously, directly resolving haplotype phase.

Practical considerations:

- **Reflex criteria**: The `Long_read_recommended` output flag identifies samples where phasing ambiguity changes the clinical interpretation. Only these samples (~2%) need the amplicon assay.
- **Throughput**: Barcoded nanopore library prep (24–96 samples per flow cell) keeps per-sample cost low at clinical volumes.
- **Turnaround**: Same-day results are feasible (long-range PCR + 4-hour sequencing run + analysis).
- **Duplications and tandems**: Each nanopore read represents a single molecule, so reads from different gene copies are naturally separated. Read count ratios reflect copy number directly.
- **Deletions (*5)**: No amplicon is produced from the deleted allele, visible as reduced read count — consistent with the depth-based deletion call from WGS.
- **Published primers**: CYP2D6-specific long-range PCR primers validated by the Gaedigk lab and the PharmVar consortium are available.

This two-tier approach (WGS-based calling for all samples, nanopore reflex for the small fraction with clinically relevant ambiguity) balances accuracy and cost-effectiveness for high-throughput clinical pharmacogenomics.

## License

Apache-2.0

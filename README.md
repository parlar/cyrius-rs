# Cyrius-RS

A Rust reimplementation of [Cyrius](https://github.com/Illumina/Cyrius), a tool for genotyping CYP2D6 from whole-genome sequencing (WGS) data. It started as a port of the Python version (bcyrius) and has since been substantially modified: it keeps the same depth / GMM / SNP-ratio / star-table flow, and adds fragment phasing, a likelihood tie-breaker, input checks, confidence scoring and opt-in QC modules. It builds as a single binary with no runtime dependencies beyond htslib.

## Benchmark

53 GeT-RM samples from the 1000 Genomes 30x WGS data (NYGC), as regional BAMs. Truth is the GeT-RM consensus genotype (Pratt et al. 2019, J Mol Diagn 21:1034-1052). A call counts as **exact** when it is a single diplotype equal to the truth, and as **ambiguous** when it lists several `;`-separated diplotypes, one of which is the truth.

Measured at commit `a75638a`:

| Metric | Value |
|---|---|
| **Exact match** | **51 / 53 (96.2%)** |
| Ambiguous (truth included) | 2 |
| Wrong | 0 |
| No call | 0 |

These 53 samples were also used while developing the calling rules, so this is not a held-out accuracy estimate.

The 2 ambiguous calls:

- **NA18565**: `*10/*36x2;*36/*36+*10`. Same variants, different tandem arrangement; short reads cannot tell them apart. **Same clinical effect** (activity 0.25, Intermediate Metabolizer, for both).
- **NA19908**: `*1/*46;*43/*45`. The variants that would separate the two options are 1.6-4.1 kb apart, and no read pair covers both (`Phase_pairs` shows 0 fragments for all 4 pairs). \*1/\*46 is activity 2.0, Normal Metabolizer; \*43 has no activity value in the function table, so \*43/\*45 is Indeterminate. The output marks this as `Clinical_equivalence: unknown_effect` with `Long_read_recommended: true`.

HG00589 and NA18973 (`*1/*21`) used to be reported as `*1/*21;*2/*21`. Fragment phasing, now on by default, resolves both (`Resolved_by: fragment_phase`).

Each call includes a confidence score (0-1) with a label (HIGH, MEDIUM, LOW). On the 53 truth samples: 48 HIGH, 5 MEDIUM, 0 LOW.

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
| HG00589 | \*1/\*21 | 0.84 MEDIUM | \*1/\*21 |
| HG01086 | \*31/\*178 | 0.94 HIGH | \*178/\*31 |
| HG01094 | \*1/\*31 | 0.95 HIGH | \*1/\*31 |
| HG01108 | \*2/\*106 | 0.95 HIGH | \*2/\*106 |
| HG01190 | \*5/\*68+\*4 | 0.61 MEDIUM | \*68+\*4/\*5 |
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
| NA12154 | \*33/\*68+\*4 | 0.85 MEDIUM | \*68+\*4/\*33 |
| NA12878 | \*3/\*68+\*4 | 0.91 HIGH | \*3/\*68+\*4 |
| NA18526 | \*1/\*36+\*36+\*10 | 0.87 HIGH | \*1/\*36+\*36+\*10 |
| NA18544 | \*10/\*41 | 0.96 HIGH | \*10/\*41 |
| NA18545 | \*5/\*36x2+\*10x2 | 0.96 HIGH | \*5/\*36x2+\*10x2 |
| NA18563 | \*1/\*36+\*10 | 0.94 HIGH | \*1/\*36+\*10 |
| NA18564 | \*2/\*36+\*10 | 0.96 HIGH | \*2/\*36+\*10 |
| NA18565 | \*10/\*36x2;\*36/\*36+\*10 | 0.87 HIGH | \*10/\*36x2 |
| NA18572 | \*36+\*10/\*41 | 0.95 HIGH | \*36+\*10/\*41 |
| NA18617 | \*36+\*10/\*36+\*10 | 0.93 HIGH | \*36+\*10/\*36+\*10 |
| NA18632 | \*36+\*36+\*10/\*52 | 0.85 MEDIUM | \*36+\*36+\*10/\*52 |
| NA18642 | \*1+\*90/\*36+\*10 | 0.94 HIGH | \*36+\*10/\*1+\*90 |
| NA18959 | \*2/\*36+\*10 | 0.94 HIGH | \*2/\*36+\*10 |
| NA18973 | \*1/\*21 | 0.87 HIGH | \*1/\*21 |
| NA18980 | \*2/\*36+\*10 | 0.95 HIGH | \*2/\*36+\*10 |
| NA19143 | \*10/\*45 | 0.95 HIGH | \*2(\*45)/\*10 |
| NA19207 | \*2x2/\*10 | 0.93 HIGH | \*2x2/\*10 |
| NA19317 | \*5/\*5 | 1.00 HIGH | \*5/\*5 |
| NA19777 | \*1/\*82 | 0.96 HIGH | \*1/\*82 |
| NA19785 | \*1/\*13+\*2 | 0.95 HIGH | \*1/\*13+\*2 |
| NA19819 | \*2/\*4x2 | 0.88 HIGH | \*2/\*4x2 |
| NA19908 | \*1/\*46;\*43/\*45 | 0.81 MEDIUM | \*1/\*46 |
| NA19917 | \*1/\*40 | 0.96 HIGH | \*1/\*40 |
| NA19920 | \*1/\*4x2 | 0.87 HIGH | \*1/\*4x2 |
| NA20289 | \*6/\*11 | 0.94 HIGH | \*6/\*11 |
| NA20803 | \*2/\*22 | 0.91 HIGH | \*2/\*22 |
| NA20875 | \*1/\*111 | 0.89 HIGH | \*1/\*111 |
| NA21105 | \*3/\*111 | 0.96 HIGH | \*3/\*111 |

</details>

To run the benchmark (the harness and truth table are versioned in `bench/`; BAMs are read from `benchmark/` or `BENCH_BAM_DIR`):

```bash
bash bench/run_benchmark.sh --strand-bias-all --fuzzy-match --changepoint-hybrid --het-check
```

On these 62 BAMs (53 with truth, 9 without) the calls are identical with these flags, with no flags, and with `--strand-bias-all` or `--changepoint-hybrid` alone (Genotype, Filter, CNV group and Confidence compared).

### Other sequencing pipelines

The same code (`a75638a`, flags as above) on regional extracts of public BAM/CRAM files made by other pipelines. Truth: GeT-RM for NA12878/HG001 (`*3/*68+*4`) and the 1kGP samples; HG002 `*2/*4` (GIAB Q100 assembly); HG005 `*36+*10/*49` (literature, not GeT-RM).

| Pipeline (samples) | Inputs | Exact | Ambiguous | Wrong | No call |
|---|---|---|---|---|---|
| DRAGEN 4.2.7, graph reference (1kGP 30x, the 53 benchmark individuals) | 53 | 51 | 2 | 0 | 0 |
| DRAGEN 4.0.3, graph reference (GIAB HG001, HG002, HG005; 1kGP NA12878) | 4 | 3 | 0 | 0 | 1 |
| DRAGEN 4.4.7, multigenome (same 4) | 4 | 3 | 0 | 0 | 1 |
| DRAGEN 4.2.4 (Platinum Pedigree NA12878) | 1 | 1 | 0 | 0 | 0 |
| bwa-mem, Google NovaSeq 30x (HG001, HG002, HG005) | 3 | 3 | 0 | 0 | 0 |
| bwakit + bwa-postalt, EBI (NA12878) | 1 | 0 | 0 | 0 | 1 |
| 10x Long Ranger 2.1 / 2.2 (GIAB HG001, HG002, HG005) | 3 | 1 | 0 | 0 | 2 |
| vg giraffe, Complete Genomics T7 (HG002, HG005) | 2 | 1 | 0 | 1 | 0 |
| Ultima UG 100 (NA12878, HG002, HG005 x2) | 4 | 0 | 0 | 3 | 1 |
| novoalign 3.02.07, GIAB 300x (HG001) | 1 | 0 | 0 | 0 | 1 |

Wrong calls are not always flagged. vg giraffe HG005 (`*5/*49`) and Ultima HG002 (`*2/*160`) come out with Filter `Fuzzy_match`, not as no-calls. In the Ultima HG002 file, no read at the main \*4 site carries the variant (0 of 19 reads, against 18 of 44 in the bwa file of the same sample); the cause is not known yet. **Treat Ultima and vg giraffe input as unsupported.**

## How It Works

The caller follows a multi-stage pipeline:

### 1. Depth normalization and copy number calling

Reads are counted across the CYP2D6/CYP2D7 region and normalized against 50 control regions on 22 chromosomes; control regions with zero or near-zero depth are left out. A Gaussian mixture model (GMM) assigns the total D6+D7 copy number (normal diploid = 4) and spacer copy number.

### 2. D6/D7 SNP ratio analysis

117 diagnostic positions that differ between CYP2D6 and CYP2D7 are examined. At each site, the fraction of reads supporting the D6 vs D7 allele is used to estimate the per-site D6 copy number. A consensus across all sites determines the CNV structure (e.g., deletion, duplication, exon9 hybrid/gene conversion).

### 3. Variant calling

Diagnostic variants that define star alleles are called from read pileups. This includes SNVs, small indels, and the exon9 gene conversion marker. Strand bias filtering removes artifacts at noisy positions.

### 4. Star allele matching

The observed variant pattern is matched against a precomputed table of all possible star allele combinations for the given CNV group. When no exact match is found, fuzzy matching with edit-distance scoring is used as a fallback (on by default; Filter `Fuzzy_match`).

### 5. Tie-breaking

When several diplotypes fit (a `;` call), two steps try to pick one:

- **Fragment phasing** (cn2 calls only, on by default). For each pair of variants that the options place differently, it counts read fragments (mates merged) that cover both sites: both-ALT and both-REF fragments say "same copy", one-of-each says "different copies". A both-REF fragment counts only if one of its reads has MAPQ >= 10, so randomly placed CYP2D7 reads are not counted. A pair counts when at least 3 fragments cover it and at least 70% agree; an option wins with a score of at least 3 and at least 3x the other. The counts are reported in `Phase_pairs`.
- **Likelihood resolver.** Scores each remaining option's allele set against the reads and narrows the call only on a clear winner; the top posterior is reported in `Resolver_posterior`.

`Resolved_by` says which step, if any, narrowed the call.

### 6. Phenotype prediction

Activity scores are assigned per CPIC guidelines, and the predicted metabolizer phenotype (Poor, Intermediate, Normal, Ultrarapid) is reported.

## Changes from bcyrius

cyrius-rs keeps bcyrius's module boundaries, data resources and depth/GMM/SNP/matching flow, with these changes active by default:

- **Exon 9 GC threshold relaxed (0.75 → 0.70)** in `call_exon9gc` ([`call_variants.rs:376`](cyrius-rs/src/caller/call_variants.rs#L376)). bcyrius's stricter cutoff misses rare `*10D`-style duplication breakpoints downstream of exon 9 (e.g. HG00421, where the exon 9 conversion probability is 0.729). The relaxed threshold lets the exon 9 signal participate in CNV-tag correction.
- **Spacer CN as a structural discriminator** in star allele matching ([`match_star_allele.rs`](cyrius-rs/src/caller/match_star_allele.rs)). bcyrius uses the spacer region only for depth normalization. cyrius-rs uses the spacer CN as a feature to discriminate tandem and deletion-plus-duplication configurations that produce indistinguishable gene-body signals (e.g. NA18545, `exon9hyb_star5` cases).
- **Tandem notation normalized** (`*36+*36` → `*36x2`) to match Pratt 2019 truth conventions.
- **Fuzzy matching and fragment phasing on by default** (see [How It Works](#how-it-works)); `--no-fuzzy-match` and `--no-phase-disambiguate` turn them off.
- **Ambiguous calls reported, not hidden** — when short reads cannot uniquely determine phase or structure, alternatives are reported as `dip1;dip2` rather than picking one. JSON fields `Phasing_ambiguous`, `Clinical_equivalence`, and `Long_read_recommended` distinguish clinically-equivalent ambiguity (no action) from ambiguity that changes, or may change, the clinical result (long-read reflex recommended).
- **Exact copy counts multiply the single-copy activity** ([`phenotype.rs`](cyrius-rs/src/phenotype.rs)). bcyrius, and earlier cyrius-rs, applied the `*Ax≥3` row's activity to every copy number of three or more, so `*41x4/*10` scored 1.0 (Intermediate) instead of 1.25 (Normal). An exact table row still wins; otherwise the activity is the single-copy activity times the exact count, and only the function label comes from the `≥3` row.
- **Disagreeing mates make a phasing site uninformative** ([`phase_disambiguate.rs`](cyrius-rs/src/caller/phase_disambiguate.rs)). When two reads of one fragment disagree at a site, the fragment no longer reports whichever mate was read last; the site is dropped, so the result no longer depends on BAM record order.
- **Bases below quality 13 are not phasing evidence** (`phase_disambiguate.rs`, and every read that `get_bases_per_read` extracts in `haplotype.rs`). Omitted qualities also leave the site unresolved: in SAM/BAM a quality of `0xFF` means "none recorded", not "a perfect base", and a naive `q < 13` test would let it through. The main dosage path for the 140 target variants (`snp_count.rs`) is unchanged.

  The `haplotype.rs` half of that gate reaches further than the g.42126938 phase rule. `get_bases_per_read` is also the only evidence path for three variants that are *called*, not merely phased — `g.42126938C>T` (`call_var42126938`), `g.42127526C>T` and `g.42127556T>C` (`call_var42127526_var42127556`) — and none of the three is in `CYP2D6_target_variant_38.txt`. On a BAM whose base qualities are omitted, every read is dropped at those positions, so the three are never called; `Uncallable_target_sites` cannot report them, because they are not target sites; and the genotype is still returned `PASS`. That path is fail-open, this release does not detect a quality-stripped input, and the fix is filed as its own plan rather than made here.
- **The targeted cn2 phase rules keep both diplotypes when phase is unknown** ([`match_star_allele.rs`](cyrius-rs/src/caller/match_star_allele.rs)). An absence of evidence used to be read as the opposite phase and silently removed one candidate. The unsupported negative branch of the g.42126938 signal was removed rather than replaced, so that signal can only say "yes" or "unknown"; no affirmative non-G pattern is established.
- **A BAM read error is an error, not a short count.** A decoding failure partway through a region now stops the step (a no-call for mandatory evidence; the prior ambiguity is preserved for optional scoring) instead of returning the count of the prefix that was read.
- **Confidence scoring** combines six signals (depth quality, CN rounding quality, match quality, variant completeness, variant specificity, D6/D7 SNP consistency) into a single score with HIGH/MEDIUM/LOW labels.
- **Input checks** (see [Input checks](#input-checks)): a BAM header showing an aligner step that drops multi-mapped reads turns PASS into `Multimapper_depleted`; a failed BAM read gives a no-call instead of a crash; control regions with no reads are left out of normalization.
- **Provenance in every output record**: tool version, git commit and every effective flag.
- **Changepoint-based CNV-tag fallback** ([`changepoint.rs`](cyrius-rs/src/caller/changepoint.rs)). When the standard region-consensus rule cannot assign a CNV tag, CBS with Welch's t-statistic on the D6/(D6+D7) ratio profile is used as an alternative path. Enabled by `--changepoint-hybrid`, which is included in the recommended flag set.
- **Several robustness fixes**: underflow guard in `cnv_hybrid` (`total_cn − spacer_cn`), insertion-order preservation in `get_consensus` (matches Python `Counter` tie-break), missing-file handling that returns empty output instead of panicking on autodetect, and clinical-equivalence comparison that includes both activity score and phenotype.

A larger collection of opt-in flags adds independent QC channels and prototype analyses — see [Feature flags](#feature-flags) below.

## Features

- Genotypes CYP2D6 star alleles from BAM/CRAM files aligned to GRCh38 (GRCh38 only; `chr22` or `22` contig names are detected automatically)
- Detects copy number variants (deletions, duplications, hybrids, gene conversions)
- Reports activity scores and predicted metabolizer phenotypes based on CPIC guidelines
- Optionally outputs population frequencies and haplotype functional annotations
- Single binary, no Python/pip required

## Input checks

cyrius-rs reads copy number from reads that fit CYP2D6 and CYP2D7 equally well, which aligners report with MAPQ 0. A pipeline step that drops those reads changes the evidence.

- **Header check (sets the Filter).** If an `@PG` line shows novoalign without `-r` (multi-mapped reads not reported), or `samtools view -q N` with N >= 1, a PASS or `Fuzzy_match` filter becomes `Multimapper_depleted`. The genotype is still reported. JSON: `Header_drops_multimappers`.
- **REP7 MAPQ-0 fraction (reported only).** The share of reads starting in REP7 that have MAPQ 0. JSON: `REP7_MAPQ0_fraction`. It is not used for the filter: graph and linked-read pipelines (DRAGEN 4.4.7, vg giraffe, 10x Long Ranger) gave low values on correct calls.
- **Uneven coverage.** `Coverage_MAD` above 0.11 turns PASS into `LowQ_coverage_MAD`.
- **Read failures.** A BAM/CRAM read that fails gives a no-call with an error in the log, not a partial result.

The header check only sees what the header records. A file whose header was stripped or rewritten is not checked.

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

#### On by default

| Flag | Effect on output |
|---|---|
| `--fuzzy-match` | Edit-distance fallback when no exact star allele match is found. Recovers near-matches at the cost of a `Fuzzy_match` filter label. The flag is accepted for compatibility; `--no-fuzzy-match` turns it off. |
| `--phase-disambiguate` | Fragment phasing to resolve ambiguous CN=2 multi-match calls into a single diplotype (see [Tie-breaking](#5-tie-breaking)). The flag is accepted for compatibility; `--no-phase-disambiguate` turns it off. |

#### Recommended set

The benchmark is run with:

```bash
--strand-bias-all --fuzzy-match --changepoint-hybrid --het-check
```

| Flag | Effect on output |
|---|---|
| `--strand-bias-all` | Strand-bias filtering on all variant sites (not just the small known-noisy set). |
| `--changepoint-hybrid` | Circular binary segmentation on the D6/(D6+D7) ratio profile, used as a CNV-tag fallback when the standard region-consensus rule produces no tag. |

Both act on the call (variant filtering and a CNV-tag fallback); on the 62 benchmark BAMs they change nothing. `--het-check` (below) only logs.

#### Off by default (can change the call)

These also affect the genotype or filter, but are off by default because they are either situation-specific or have not yet been validated across the full benchmark:

| Flag | Effect on output |
|---|---|
| `--hmm-cnv` | HMM segmentation on the per-position CN profile, used as a second CNV-tag fallback after changepoint detection. |
| `--phase-readpair` | Extend phasing to mate-pair distance (no effect with `--no-phase-disambiguate`). |
| `--spacer-cn-check` | Flag calls with `Spacer_CN_inconsistent` when the observed spacer-vs-total-CN difference does not match the CNV group. |
| `--consistency-check` | Post-hoc validation; sets the filter to the first failing check, e.g. `Conversion_map_chimeric_in_cn2`, `Allele_mismatch_rate_<N>pct` or `Allele_balance_deviation_<N>pct`. |

#### QC and validation (do not affect the call)

These flags add **independent observations** to the output or the log without modifying the genotype, filter, or confidence:

| Flag | What it adds |
|---|---|
| `--het-check` | Logs D6/D7 hemizygosity (zero het sites), which detects `*5/*5`. Log only. |
| `--kmer-validation` | Alignment-free paralog-dosage check (189 diagnostic 31-mers). Emits a `Kmer_validation` JSON block with AGREE/DISAGREE. |
| `--clip-evidence` | Soft-clip cluster scan in the D6+D7+spacer region; flags candidate structural breakpoints. |
| `--d7-depth` | Per-region CYP2D7 depth profile; logs elevation/depletion signals consistent with `*36` and `*13` hybrids. |

## Output

Results are written as both JSON and TSV files. The TSV contains:

```
Sample       Genotype   Filter   Confidence     Activity score   Predicted phenotype   Phasing_ambiguous   Long_read_recommended
sample123    *1/*4      PASS     0.96 (HIGH)    1.0              Intermediate Metabolizer   false           false
```

### Filter values

| Filter | Meaning |
|---|---|
| `PASS` | Single diplotype, exact star-table match, no QC flag. |
| `Fuzzy_match` | No exact match; the closest pattern by edit distance was used. |
| `More_than_one_possible_genotype` | Several diplotypes fit (`;`-separated). |
| `Not_assigned_to_haplotypes` | The variants matched but could not be split into two haplotypes. |
| `LowQ_high_CN` | The GMM gives no copy number and the depth estimate is above 7.5; the CN is that estimate rounded. |
| `LowQ_coverage_MAD` | Uneven coverage (`Coverage_MAD` > 0.11); replaces PASS. |
| `Multimapper_depleted` | The header shows a step that drops multi-mapped reads; replaces PASS or `Fuzzy_match`. |
| `Spacer_CN_inconsistent`, consistency-check values | Only with `--spacer-cn-check` / `--consistency-check`. |
| `None` (null) | No call. |

The confidence score (0–1) is a weighted combination of six quality signals: depth quality, CN rounding quality, match quality, variant completeness, variant specificity, and D6/D7 SNP ratio consistency. Labels: HIGH (≥0.85), MEDIUM (0.50–0.84), LOW (<0.50).

### Phasing ambiguity and long-read reflex

When the data does not uniquely determine a single diplotype, the output contains multiple alternatives separated by `;`. Two additional fields help clinical labs decide when orthogonal confirmation is needed:

- **Phasing_ambiguous**: `true` when multiple alternative diplotypes are reported
- **Long_read_recommended**: `true` when the alternatives have **different clinical effects** (different activity score or metabolizer phenotype), or when an alternative has no activity value. When all alternatives are clinically equivalent, long-read confirmation would not change the clinical interpretation.

The JSON output additionally includes `Clinical_equivalence` for ambiguous calls: `same_effect`, `different_effect` or `unknown_effect` (`n/a` when the call is not ambiguous).

The JSON output also includes: total copy number, spacer copy number, raw depth values, CNV group, called variants, confidence score, the full variant count table (`Variant_raw_count`, `alt,ref` read counts per site), the tie-breaking fields `Resolved_by`, `Phase_pairs` and `Resolver_posterior`, and the input-check fields `Header_drops_multimappers` and `REP7_MAPQ0_fraction`. Every record carries its provenance: `Tool_version`, `Tool_commit` (the git commit the binary was built from, `-dirty` if the crate had uncommitted changes) and `Flags` (every feature flag with its effective value). `--version` prints the same version and commit.

Two further JSON fields report evidence that the call does not use:

- **`Uncallable_target_sites`**: the target variants that received no dosage call from the main dosage path. A site gets there by one of three routes: no usable reads at the site; too few reads for any copy number to reach its posterior cutoff; or the strand-bias filter (`NOISY_VAR`) having zeroed the alt count. A few sites also have a dedicated caller that can resolve them from elsewhere (for example `g.42130655-42130656insA`, via the haplotype database), so a listed site means "the main dosage path could not call it", not "nothing is known about it". The field errs towards over-reporting missing evidence, which is the safe direction for a warning. It is **reported only**: it does not change the genotype, the `Filter` or the confidence in this release, and the matcher still treats an unobserved site as reference. That is a known defect (review item 1) whose fix is deferred to its own plan. Measured on NA19920 with the 68 reads at the `*4` splice site removed: 5 sites are listed (3 with zero reads, 2 with too few at total CN 5), and the call moves from `*1/*4x2` (0.87 HIGH, activity score 1.0, **Intermediate Metabolizer**) to `*1/*10x2` at 0.95 HIGH, still `PASS` — activity score 1.5, **Normal Metabolizer**. The actionable field changes, and the confidence rises rather than falls: the deferred defect turns an Intermediate Metabolizer into a confident Normal Metabolizer.
- **`Phase_pairs_truncated`**: whether `Phase_pairs` is the complete list. It has three states, and `null` is not `false`:
  - `null`: the phase step did not run.
  - `false`: the step ran to completion. This includes the legitimate case where it ran and had nothing to phase, so an empty `Phase_pairs` with `false` is a real result and not a failure.
  - `true`: the list was cut short, either by a BAM read error partway through or by the file not opening at all.

  Every pair listed is a complete measurement even when the list is truncated: it is the list that is short, not the numbers.

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
cyrius-rs/
├── build.rs                 # Stamps the git commit into the binary
├── data/                    # Data files embedded at build time (SNP definitions, star allele tables, GMM parameters)
├── examples/phase_audit.rs  # Prints the phasing counts for a BAM (audit tool)
└── src/
    ├── main.rs              # CLI and main pipeline
    ├── lib.rs
    ├── types.rs             # Shared types (D6Call, FeatureFlags, etc.)
    ├── data.rs              # Embedded data files
    ├── stats.rs             # Statistical utilities
    ├── fisher.rs            # Fisher's exact test
    ├── phenotype.rs         # Activity scores, metabolizer phenotype, clinical equivalence
    ├── threading.rs         # --threads handling
    ├── depth_calling/       # Depth-based copy number analysis
    │   ├── bin_count.rs     #   Read counting, normalization, input checks
    │   ├── copy_number_call.rs # CN calls from depth
    │   ├── gmm.rs           #   Gaussian mixture model
    │   ├── haplotype.rs     #   Haplotype extraction at linked sites
    │   ├── snp_count.rs     #   D6/D7 SNP ratio counting
    │   └── utilities.rs     #   BAM/CRAM opening, region and contig helpers
    └── caller/              # Star allele calling logic
        ├── call_variants.rs #   Variant calling from pileups
        ├── cnv_hybrid.rs    #   CNV group assignment from SNP ratios
        ├── match_star_allele.rs # Star allele pattern matching
        ├── construct_star_table.rs # Star allele combination tables
        ├── fuzzy_match.rs   #   Edit-distance fallback matching
        ├── phase_disambiguate.rs # Fragment phasing for cn2 ties
        ├── diplotype_caller.rs # Read-likelihood scoring of candidate diplotypes
        ├── resolver.rs      #   Likelihood tie-breaker for ; calls
        ├── confidence.rs    #   Multi-signal confidence scoring
        ├── strand_bias_all.rs # Strand bias filtering
        ├── changepoint.rs   #   Changepoint-based CNV detection
        ├── hmm_cnv.rs       #   HMM-based CNV segmentation
        ├── consistency.rs   #   Post-hoc consistency checks
        ├── het_check.rs     #   Hemizygosity logging
        ├── clip_evidence.rs #   Soft-clip cluster detection
        ├── d7_depth.rs      #   CYP2D7 sub-region depth profile
        └── kmer_validation.rs # Alignment-free k-mer dosage check
```

## Future work: long-read reflex for phasing ambiguity

Some CYP2D6 diplotype calls are inherently ambiguous from short-read WGS because the data cannot phase distant variants onto haplotypes. In our benchmark, 2 of 53 samples report multiple alternative diplotypes. In 1 case (NA18565) the alternatives are clinically equivalent (same activity score and metabolizer phenotype). In the other (NA19908, ~2%) one alternative has no activity value, so the output recommends orthogonal confirmation.

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

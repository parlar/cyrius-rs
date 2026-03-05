# Cyrius-RS

A Rust implementation of [BCyrius](../bcyrius/), a tool for genotyping CYP2D6 from whole-genome sequencing (WGS) data. Cyrius-RS is a faithful port that produces identical output to the Python version, compiled into a single static binary with no runtime dependencies (beyond htslib).

> **Warning:** This project is in early development and has not been validated against real data yet. The code is likely still buggy. Do not use for clinical or production purposes.

## Features

- Genotypes CYP2D6 star alleles from BAM/CRAM files aligned to hg38
- Detects copy number variants (deletions, duplications, hybrids, gene conversions)
- Reports activity scores and predicted metabolizer phenotypes based on CPIC guidelines
- Optionally outputs population frequencies and haplotype functional annotations
- Single binary, no Python/pip required

## Building

Requires Rust 1.70+ and htslib development libraries (for BAM/CRAM support).

```bash
cargo build --release
```

The binary will be at `target/release/cyrius-rs`.

## Usage

```bash
cyrius-rs --input sample.bam --outDir results/
```

### Parameters

| Parameter | Description |
|---|---|
| `-i, --input` | Input BAM or CRAM file (indexed) |
| `-o, --outDir` | Output directory for results |
| `-r, --reference` | Reference FASTA (recommended for CRAM) |
| `-g, --genome` | Contig naming: `autodetect` (default), `38`, or `chr38` |
| `--id` | Sample ID (default: derived from filename) |
| `-t, --threads` | Number of threads (default: 1) |
| `--population-info` | Include population frequencies in output |
| `--haplotype-info` | Include haplotype functional annotations |
| `--print` | Print results to stdout |
| `--countFilePath` | Use precomputed count files instead of BAM |

## Output

Results are written as both JSON and TSV files. The TSV contains:

```
Sample       Genotype   Filter   Activity score   Predicted phenotype
sample123    *1/*4      PASS     1.0              Intermediate Metabolizer
```

When `--haplotype-info` and `--population-info` are used, additional rows with functional annotations and population frequency data are appended.

## License

Apache-2.0

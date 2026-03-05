# CLAUDE.md

## Project Overview

bcyrius is a Python-based pharmacogenomics star allele caller. It performs depth-based copy number calling and variant matching to assign star alleles.

## Key Principle

The Python code in this repository is the source of truth. Any rewrites, ports, or translations to other languages (e.g., Rust) must be **faithful to the existing Python implementation** — preserving the same algorithms, logic, data flow, and numerical behavior. Do not "improve" or restructure the logic during translation; replicate it exactly.

## Project Structure

- `run.py` — Entry point
- `star_caller/` — Core star allele calling logic
- `caller/` — Variant calling, star allele matching, CNV hybrid handling
- `depth_calling/` — Depth-based copy number analysis (binning, GMM, haplotyping, SNP counting)

## Guidelines

- When porting to another language, preserve function boundaries and naming where practical.
- Numerical operations (smoothing, GMM fitting, etc.) must produce equivalent results to the Python version.
- Do not add features, optimizations, or abstractions that diverge from the Python behavior without explicit approval.

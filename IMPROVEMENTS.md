# BCyrius - Bugs and Improvement Notes

## Confirmed Bugs

1. **`call_star68` mutates its input list** (match_star_allele.py:139)
   `var_observed.remove()` modifies the original list. The caller makes a backup copy but passes the original to `call_star68`. The copy should be passed instead.

2. **`bamfile.close()` not called on early returns** (star_caller.py:202-477)
   `d6_star_caller` opens a BAM file at line 202 but has early returns at lines 235, 337, 379 that skip `close()` at line 476. Use a context manager.

3. **Dead code / logic error in haplotype conflict resolution** (haplotype.py:92-94)
   Line 92 sets value to `None` on conflict, but line 94 unconditionally overwrites it with `hap`. Line 94 should be in an `else` branch.

4. **`nchr` undefined if SNP file is empty** (snp_count.py:84)
   `split_line` from the loop body is used after the loop. If the file has no data lines, `NameError` is raised.

5. **`output_per_sample` used outside loop scope** (star_caller.py:908)
   If `final_output` is empty, `output_per_sample` and `sorted_genotype` are undefined when `--print` is used.

6. **`isChrInChromosomeNames` doesn't close the alignment file** (star_caller.py:738-749)

## Likely Bugs

7. **Fisher exact test uses incorrect contingency table** (call_variants.py:198-199)
   The second row uses `ntotal/2` (expected values), not observed counts.

8. **Possible alt/ref swap in `get_supporting_reads_single_region`** (snp_count.py:223-238)
   Forward/reverse counts come from `lsnp1` (region1) but callers name them `var_alt_forward/reverse`. Verify SNP file defines region1 as alt.

## Minor Issues

9. **README documents `---haplotype-info` (triple dash)** but argparse uses `--haplotype-info` (double dash).

10. **GC values appended as strings in `process_counts_and_prepare_for_normalization`** (bin_count.py:286) but as floats elsewhere.

11. **Identical if/else branches in `call_var42130655insA`** (call_variants.py:403-407) — both append the same variant, just different counts.

12. **CRAM detection is case-sensitive** (utilities.py:67) — `.CRAM` won't match.

## Performance

### Bottleneck: `get_hap_table` takes ~25 seconds on every run

The star allele combination table builder in construct_star_table.py runs N^3 = 9.6M iterations (N=213 star alleles) to build `dhap3`. This is the dominant startup cost and the data never changes between runs.

| Table           | Iterations | Time  | Memory  |
|-----------------|-----------|-------|---------|
| `dhap3` (N^3)   | 9.6M      | ~19s  | 244 MB  |
| `dhap_exon9_x3` | 2.9M      | ~7s   | 124 MB  |
| Everything else | ~400K     | <1s   | ~90 MB  |

### Recommended fixes (in priority order)

1. **Compute `dhap3` on demand.** It's only used for 2 rare CNV tags (`exon9hyb`, `exon9hyb_star68`). Instead of precomputing all 1.5M entries, decompose observed variants into 3 star alleles at lookup time via Counter matching over 213 candidates. Eliminates the dominant cost entirely.

2. **Cache tables to disk.** Serialize with `marshal`, keyed by hash of `star_table.txt`. Tested at 13x speedup (25s to 2s), but cache is ~500MB. Less elegant than option 1.

3. **Move `loadHaplotypeFunctionality` outside the per-sample loop** (star_caller.py:833). It re-reads and parses the file on every iteration.

4. **Reduce redundant BAM passes.** The pipeline does ~8 separate pileup/fetch passes over overlapping chr22 regions. Could be consolidated.

## Code Quality Suggestions

- Split star_caller.py into separate modules (CLI, phenotype logic, frequency queries, output formatting)
- Add tests (pure functions like `sortGenotype`, `determinePhenotype`, `call_reg1_cn` are very testable)
- Standardize naming to snake_case (BCyrius additions use camelCase)
- Replace magic numbers with named constants
- Move namedtuple definitions out of functions to module level
- Add batch input support (reuse star combination tables across samples)

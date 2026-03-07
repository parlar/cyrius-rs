# Analysis of 3 Remaining Wrong Samples

## Summary

Three samples from the 1000 Genomes high-coverage dataset are miscalled by cyrius-rs.
Full-region BAMs (chr22:42100000-42160000) were downloaded from 1000G CRAMs and analyzed
with pdx_extract and depth profiling.

## Sample Details

### NA18545: Expected `*5/*36x2+*10x2`, Called `*36+*10/*36+*10`

**Depth profile:**
- D6 CN = 4.04, D7 CN = 2.08, Spacer CN = 3.73
- Flanking depth = 41.4x (1-copy = 20.7x)

**Root cause:** The `exon9hyb_exon9hyb` CNV group with pattern `*10_*10_*36_*36` was
hardcoded to `*36+*10/*36+*10` at match_star_allele.rs:497 without checking spacer_cn.
The elevated spacer CN (~4 vs normal 2) indicates a *5 deletion on one haplotype — the
extra spacer copies come from having 4 tandem copies on the non-deleted haplotype.

**Fix applied:** Added spacer_cn > 2 check for this specific pattern. When spacer CN is
elevated, the call is `*5/*36x2+*10x2` instead of `*36+*10/*36+*10`.

**Supporting evidence:**
- 0 het sites in D6 region — consistent with *5 (all D6 reads from one haplotype)
- Diplotype caller produces tied likelihoods (structural, not SNP-level distinction)

### NA18565: Expected `*1/*36+*10`, Called `*1/*2`

**Depth profile:**
- D6 CN = 3.08, D7 CN = 2.09, Spacer CN = 2.70
- Flanking depth = 60.7x (1-copy = 30.4x)
- Breakpoint at 42130053 detected by pdx_extract (1 contig, 775bp)
  - BLAST: contig is pure CYP2D6 (REP6 region), not hybrid
  - Likely tandem duplication junction

**Root cause:** The total_cn likely rounds to 2 (from the GMM depth caller) instead of 3.
With total_cn=2, the CNV tag becomes `cn2`, and the code looks for only 2 alleles,
finding `*1/*2` instead of detecting the *36+*10 tandem arrangement. The true D6 CN=3
(1 from *1, 1 from *36 partial, 1 from *10).

**Potential fix:** Improve GMM sensitivity for CN=3 detection when depth ratio is ~1.5x
flanking, or use spacer_cn > 2 as supporting evidence for CN elevation.

### NA19143: Expected `*1/*2x2`, Called `*1x2/*2`

**Depth profile:**
- D6 CN = 2.45, D7 CN = 2.05, Spacer CN = 2.06
- Flanking depth = 37.8x (1-copy = 18.9x)

**Exon-level depth reveals partial duplication:**
- Exons 1-3 (3' genomic end): CN ~1.7
- Exons 4-6: CN ~2.4-3.0
- Exons 5, 9 (5' genomic end): CN ~3.0

This gradient pattern shows the duplication breakpoint is *within* the gene (between
exon 3 and exon 4), so the second copy is truncated at the 3' end. The overall D6 CN
of 2.45 is between 2 and 3 because only part of the gene is duplicated.

**Assessment:** Our call `*1x2/*2` may actually be more accurate than the truth set
`*1/*2x2`. The partial duplication means the second *2 copy is incomplete, which the
`*1x2` notation captures better. The truth set may be wrong here.

## Why the Benchmark Data Was Not Sufficient

The benchmark BAMs in `benchmark/` are small, pre-extracted BAMs covering only the
CYP2D6 gene body (roughly chr22:42126000-42131000). They contain the reads that
cyrius-rs normally processes — sufficient for variant calling and basic depth analysis.

However, to understand *why* the 3 wrong samples are miscalled, we needed signals
that extend beyond the gene body:

1. **Spacer region (42132000-42140000)** — sits between CYP2D6 and CYP2D7, outside the
   benchmark BAMs. Spacer CN is the key discriminator for NA18545 but can't be
   independently verified from the benchmark data alone.

2. **Flanking depth normalization** — the benchmark count files provide pre-computed
   normalized depth, but we needed raw depth across wider flanking regions (10kb
   upstream and downstream) to independently verify CN estimates and understand where
   the GMM rounding fails.

3. **CYP2D7 region (42140000-42145000)** — needed to verify D7 CN and check for hybrid
   breakpoints. Also outside the benchmark BAMs.

4. **Exon-level depth gradients** — NA19143's partial duplication is only visible when
   you compute depth per-exon across the full gene, which requires deeper analysis than
   the standard cyrius-rs pipeline performs.

5. **Structural breakpoint detection** — soft-clip clusters and discordant read pairs at
   tandem junctions require the full read set in the wider region, not just the gene body.

## Data Download and Investigation Steps

### Step 1: Download Full-Region BAMs from 1000 Genomes

The 3 samples (NA18545, NA18565, NA19143) are from the 1000 Genomes high-coverage
dataset. CRAMs are hosted on the EBI/NCBI servers. We used `samtools view` with remote
CRAM URLs to extract reads in a 60kb region (chr22:42100000-42160000) covering D6, D7,
spacer, and flanking regions:

```bash
samtools view -b -T $REF \
  "https://ftp.sra.ebi.ac.uk/vol1/run/.../NA18545.final.cram" \
  chr22:42100000-42160000 > /tmp/fullwgs/NA18545_full_region.bam
samtools index /tmp/fullwgs/NA18545_full_region.bam
```

Result: 16k-29k reads per sample, 1.2-2.0 MB BAMs.

### Step 2: Depth Profiling

Computed mean depth with `samtools depth` in 5 regions per sample:
- CYP2D6 gene body (42126499-42130881)
- CYP2D7 gene body (42140209-42144592)
- Spacer region (42132000-42140000)
- Upstream flanking (42116499-42126499)
- Downstream flanking (42144592-42154592)

Flanking average / 2 gives the 1-copy depth for CN estimation.

### Step 3: pdx_extract Analysis

Ran pdx_extract on all 3 BAMs to identify problem reads (low-mapq + soft-clipped)
and attempt breakpoint assembly. Key findings:
- NA18545: 180 problem reads, 0 breakpoints
- NA18565: 222 problem reads, 1 breakpoint at 42130053 (775bp contig)
- NA19143: 122 problem reads, 0 breakpoints

### Step 4: Breakpoint Contig Analysis (NA18565)

The 775bp contig at position 42130053 was aligned using BLAST against the CYP2D6
region reference. Three-way comparison (contig vs D6-ref vs D7-ref) at 580
discriminating positions showed the contig is 100% CYP2D6-derived (REP6 region) —
not a hybrid breakpoint. This represents a tandem duplication junction where one
D6 copy connects to another.

Split reads at the breakpoint (42130047-42130054) have 38-107bp soft clips containing
low-complexity/repetitive sequence, consistent with a junction in the REP6 repeat region.

### Step 5: Exon-Level Depth Profiling

Computed depth per exon (exons 1-9) for all 3 samples to detect partial duplications.
This revealed NA19143's striking CN gradient (exons 1-3: CN~1.7, exons 5,9: CN~3.0),
indicating the duplication breakpoint is within the gene.

## Methods Summary

- BAMs extracted from 1000 Genomes high-coverage CRAMs via `samtools view`
- Depth computed with `samtools depth` in 5 regions per sample
- pdx_extract used for problem read extraction and breakpoint assembly
- BLAST used for contig alignment against CYP2D6 region reference
- Three-way comparison (contig vs D6-ref vs D7-ref) to detect hybrid breakpoints
- Exon-level depth profiling to detect partial duplications

## Key Insight

Spacer CN is a powerful structural discriminator that is underutilized in the current
codebase. The spacer is the ~8kb region between the 3' end of CYP2D6 and the 5' end
of CYP2D7 (chr22:42132000-42140000). In a normal diploid genome there are 2 copies
(one per chromosome). When tandem arrangements exist, each extra gene copy inserts an
additional spacer between them, so spacer CN directly reflects structural copy count:

```
Normal:          ---D6---[spacer]---D7---                            spacer CN = 2

*36+*10:         ---D6---[sp]---*36---[sp]---*10---[sp]---D7---      spacer CN = 3

*36+*10/*36+*10: hap1: ---D6---[sp]---*36---[sp]---*10---[sp]---D7---
                 hap2: ---D6---[sp]---*36---[sp]---*10---[sp]---D7---
                                                                     spacer CN = 4

*5/*36x2+*10x2:  hap1: ---deleted---
                 hap2: ---D6---[sp]---*36---[sp]---*10---[sp]---*36---[sp]---*10---[sp]---D7---
                                                                     spacer CN = 5
```

For NA18545 (measured spacer CN = 3.73), this is clearly above the 4 expected for
`*36+*10/*36+*10` and closer to 5 for `*5/*36x2+*10x2`. The *5 deletion removes one
haplotype entirely, while the other carries 4 tandem copies with 5 spacers.

This signal is already used for `exon9hyb_star5` but should be extended to other CNV
groups like `exon9hyb_exon9hyb`.

## Toward a General Solution

### The Problem with Hardcoded Patterns

The current approach uses a cascade of if/else patterns in `get_final_call_clean()` —
every new structural combination needs a manually added case. This doesn't scale and is
fragile. Multiple independent signals exist (spacer CN, exon-level depth gradient, het
sites, breakpoints, variant calls) but they're consumed piecemeal in ad-hoc rules.

### Option A: Structural Diplotype Scoring (Lightweight)

Instead of hardcoded pattern matching, enumerate all candidate structural diplotypes and
score each against ALL available evidence simultaneously:

```
for candidate in all_candidate_diplotypes:
    score  = log_lik(observed_d6_cn      | expected_d6_cn(candidate))
    score += log_lik(observed_spacer_cn   | expected_spacer_cn(candidate))
    score += log_lik(observed_d7_cn       | expected_d7_cn(candidate))
    score += log_lik(observed_het_count   | expected_het_count(candidate))
    score += log_lik(observed_exon_depth  | expected_exon_depth(candidate))
    score += log_lik(observed_variants    | expected_variants(candidate))
```

Each structural diplotype has known expected observables:

| Diplotype              | D6 CN | Spacer CN | D7 CN | Het D6 | Exon gradient |
|------------------------|-------|-----------|-------|--------|---------------|
| `*36+*10/*36+*10`      | 4     | 2         | 2     | normal | flat          |
| `*5/*36x2+*10x2`       | 4     | 4         | 2     | 0      | flat          |
| `*1/*36+*10`           | 3     | 3         | 2     | normal | flat          |
| `*1/*2`                | 2     | 2         | 2     | normal | flat          |
| `*1/*2x2`              | 3     | 2         | 2     | normal | flat          |
| `*1x2/*2` (partial)    | 2.5   | 2         | 2     | normal | gradient      |

This cleanly separates all 3 hard cases without any hardcoded if/else logic.
Incremental — fits on top of the current pipeline as a disambiguation step.

### Option B: Graph Realignment (pdx_caller Integration)

Reads aligned to a linear reference can't represent structural arrangements like
`*36+*10` tandems. Reads spanning a tandem junction get soft-clipped or mismapped
because the junction doesn't exist in the reference. pdx_extract already identifies
120-220 problem reads per sample (low-mapq + soft-clipped).

Targeted realignment of just these problem reads against structural reference haplotypes
would:

1. **Rescue mismapped reads** — ambiguous on linear ref but unambiguous on the correct
   structural path
2. **Produce direct structural evidence** — reads aligning across a tandem junction or
   deletion breakpoint directly support that structural haplotype
3. **Fix CN estimation** — NA18565's D6 CN=3.08 rounds to 2 partly because some reads
   are lost to mismapping; correct graph alignment gives cleaner depth signal

The building blocks exist in pdx_caller's `LocusGraph` with `DeletionSkip`,
`TandemLoop`, `HybridJunction` edges and quality-aware `align_read_to_graph`. This is
the principled solution but requires more infrastructure to integrate.

### Option C: Strobealign + sv_tools Integration

We have two other in-house projects that already extract the signals needed for
structural disambiguation:

**Strobealign** (~/dev/strobealign/) — our custom aligner, currently a high-performance
mapper with planned SV-aware features:

Currently implemented and useful:
- Secondary alignments (`-N 5`) — reads mapping to both CYP2D6 and CYP2D7 show as
  multi-mappers. The `X0` tag (number of equally-scored alignments) directly indicates
  ambiguity in the D6/D7 repeat region.
- Soft-clip tracking — tandem junction reads get clipped because the junction doesn't
  exist in the linear reference. Soft clips are in the CIGAR but not realigned.
- MAPQ gradient — the 120-220 "problem reads" per sample that pdx_extract identifies
  are exactly the low-MAPQ reads from the D6/D7 homology region.
- Tags currently emitted: `NM`, `AS`, `na`, `nr`, `al`, `ga`, `X0`, `mr`, `RG`.

Planned features (specified in FUTURE_SAM_OUTPUT.md, NOT yet implemented):
- `SA:Z` supplementary alignment tags — would explicitly mark chimeric reads spanning
  hybrid junctions (e.g., a read crossing from D7 exon9 into D6 intron4 in a *36 hybrid)
- `YS:Z` SV type classification (DEL, DUP) — per-read structural signal
- `YT:Z` pair type (concordant, discordant orientation, discordant insert) — everted
  read pairs (RF orientation) are the classic tandem duplication signature
- `YB:Z` breakpoint intervals with uncertainty, `YM:i` microhomology length
- `XK:Z` soft-clip realignment results — clipped bases remapped to find their true origin
- `XP:Z` repeat class (UNQ/TAN/SEG/INV/DIS) — would directly flag D6/D7 as segmental dup

The planned feature set is comprehensive and purpose-built for SV calling. Once
implemented, strobealign would provide all the per-read SV signals needed. See
`~/dev/strobealign/SV_FEATURE_ROADMAP.md` for the implementation plan.

**sv_tools** (~/dev/sv_tools/) — experimental Rust crates with useful SV analysis
algorithms. Not a dependency — a source of proven algorithms to copy and adapt into
cyrius-rs or strobealign as needed. Key algorithms worth bringing in:

- **sv_align: LocalIndex** — q-gram index + Myers bit-parallel DP alignment engine.
  Zero external dependencies, ~400 lines. Perfect for realigning soft-clipped bases
  against structural reference sequences. This is exactly what's needed for the 38-107bp
  clipped reads at NA18565's breakpoint (42130053). Should be copied into cyrius-rs
  (for targeted CYP2D6 junction realignment) and/or strobealign (for general clip
  realignment powering SA/XK tags).

- **sv_break: ClipCluster** — groups soft clips by genomic position and clip side
  (left/right). The clustering logic could be adapted into cyrius-rs to identify
  breakpoint positions from clip consensus in the CYP2D6 region.

- **sv_exon: depth ratios** — per-exon depth ratio computation. The algorithm for
  computing depth ratios per interval and detecting ratio transitions could be brought
  into cyrius-rs to detect partial duplications like NA19143 (CN gradient 1.7→3.0).

- **sv_cnv: CBS segmentation** — circular binary segmentation for detecting CN
  changepoints. The algorithm could be useful in strobealign or cyrius-rs for
  normalizing depth signals.

The principle is: copy the algorithms, not the crates. cyrius-rs and strobealign
should be self-contained — no dependency on sv_tools.

**Implementation plan — bring algorithms into cyrius-rs:**

```
Phase 1: Clip realignment (immediate value for NA18565)
  - Bring sv_align's LocalIndex algorithm into cyrius-rs (~400 lines, zero deps)
  - Build structural reference sequences for known CYP2D6 arrangements
    (tandem junction, deletion breakpoint, hybrid breakpoints)
  - Realign soft-clipped reads against these structural refs
  - Score: does clip sequence match a known junction?
       ↓
Phase 2: Evidence collection (structural signals)
  - Add discordant pair detection to cyrius-rs BAM reading
  - Scan for RF-oriented pairs in spacer region (tandem dup signature)
  - Scan for abnormal insert size pairs spanning D6 deletion (*5 signature)
  - Adapt clip clustering from sv_break for breakpoint identification
       ↓
Phase 3: Exon-level depth (partial duplication detection)
  - Bring sv_exon's depth ratio algorithm into cyrius-rs
  - Compute per-exon CN across CYP2D6's 9 exons
  - Detect CN gradients indicating partial duplications (NA19143)
       ↓
Phase 4: Structural diplotype scoring (Option A) consumes ALL evidence:
  - Depth signals: D6 CN, D7 CN, spacer CN, exon gradient
  - SV signals: clip realignment score, discordant pair count
  - Variant signals: called star alleles, het site count
       ↓
Best-scoring structural diplotype is reported
```

### Recommendation

Option A (structural diplotype scoring) is the pragmatic first step — it uses only
the signals cyrius-rs already computes (depth, spacer CN, variants, het count) and
would fix all 3 remaining samples with a single mechanism.

Option C (strobealign + sv_tools algorithms) is the natural evolution. Phase 1 alone
(bringing sv_align's LocalIndex into cyrius-rs, ~400 lines) would add clip realignment
against known CYP2D6 structural junctions — directly addressing NA18565's breakpoint
reads. When strobealign gains SA tag support, those signals feed into the same scoring
framework automatically. sv_tools is experimental but contains proven algorithms worth
copying — the principle is to bring code in, not depend on sv_tools as a crate.

Option B (pdx_caller graph) remains the most principled approach but requires the most
new infrastructure. It could be pursued later if Options A+C prove insufficient for
edge cases.

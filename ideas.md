BCyrius Improvement Opportunities
1. Replace the Fixed Error Rate with Base-Quality-Aware Likelihoods
The problem (copy_number_call.py):


ERROR_RATE = 0.1  # hardcoded for all sites, all reads, all samples

if i == 0:
    depthexpected = (ERROR_RATE / 3) * float(nsum)
Every base at every position is assumed to have a 10% error rate. A Q30 base has actual error rate 0.001. A Q10 base has 0.1. BCyrius ignores this 100× difference.

The fix: Instead of counting reads and feeding counts into Poisson, compute per-read likelihoods using base quality:


# Per-read likelihood at a differentiating site
for read in pileup:
    bq = read.alignment.query_qualities[read.query_position]
    error = 10 ** (-bq / 10)
    if read_base == d6_allele:
        P_d6 = 1 - error
        P_d7 = error / 3
    elif read_base == d7_allele:
        P_d6 = error / 3
        P_d7 = 1 - error
Then the CN likelihood at each site becomes a product over reads rather than a Poisson on counts. This would sharpen the posterior substantially at high-quality sites while appropriately down-weighting low-quality bases.

Impact: Most significant at low coverage (<20x) where every read matters, and at sites near the end of reads where base quality drops.

2. Learn Per-Site Reliability Weights
The problem: All 109 differentiating sites are treated equally in the consensus voting, but they aren't equally informative. From the data file, the sites span diverse contexts:


42126611  GTCACCAGGAAAGCAA  vs  GTCACCAGAAAGCTGA   exon9 (multi-base)
42127634  c                 vs  a                    exon7 (clean SNP)
42131980  c                 vs  t                    upstream_exon1
Some sites are near indels (the exon9 multi-base difference), some are in repetitive sequence, some may have population-specific polymorphisms that confound the D6/D7 distinction.

The fix: From a training cohort (or from the GC-normalization control regions in CYP2D6_region_38.bed), estimate per-site variance. Sites with high variance across known CN=2 samples are unreliable. Weight them down:


# Instead of Counter-based voting:
consensus = Counter(cn_calls).most_common(1)

# Use weighted consensus:
weighted_cn = sum(w[i] * cn_calls[i] for i in range(n_sites)) / sum(w)
Impact: Reduces false hybrid calls caused by noisy sites, especially in the exon9 region where the multi-base differentiating site is harder to call.

3. Exploit the Physical Phasing That Already Exists
The situation: BCyrius's haplotype.py already does read-level physical phasing — it builds per-read haplotype strings across multiple differentiating sites:


dread_hap[read_name] = "11x2x1"  # this read carries D6 at pos 1,2 and D7 at pos 6
But this information seems underutilized. The star allele matching in match_star_allele.py uses per-variant CN counts, not phased haplotypes.

The opportunity: When BCyrius reports More_than_one_possible_genotype, the physical phasing from reads spanning multiple variant sites could resolve the ambiguity. For example:

Two variants A and B are each present at CN=1
Ambiguous: *2/*41 or *1/*69 (both explain A at CN=1 + B at CN=1)
But if reads carrying A also carry B → they're on the same haplotype → *69 (which has both) on one chromosome
At 30x WGS with 150bp reads, adjacent variants within ~500bp are frequently spanned by the same read or read pair. In the ~6kb CYP2D6 coding region, many star-allele-defining variants are close enough for physical phasing.

Impact: Directly reduces the More_than_one_possible_genotype filter rate — arguably the most clinically frustrating outcome.

4. Use Read Pairs for Structural Variant Detection
The problem: BCyrius detects SVs entirely from depth patterns at differentiating sites. It misses information from:

Discordant read pairs: A pair where one read maps to D6 and the mate maps to D7 signals a hybrid junction
Split reads: A read that is soft-clipped at the hybrid breakpoint directly shows the junction sequence
Insert size anomalies: Pairs spanning a deletion have larger-than-expected insert sizes
StellarPGx uses graphtyper's genotype_sv which does consider split reads and discordant pairs (via the template SV VCFs). BCyrius's depth-only approach can detect that a hybrid exists but can't pinpoint the breakpoint as precisely.

The fix: After the depth-based CNV calling identifies a likely hybrid, scan reads in the transition region for:


# Discordant pairs: one mate in D6, other in D7
for read in bamfile.fetch(d6_region):
    mate_pos = read.next_reference_start
    if d7_start <= mate_pos <= d7_end:
        breakpoint_candidates.append(read.reference_end)
Impact: Better breakpoint resolution for hybrids, especially novel ones not in the current CNV_ACCEPTED list. Combined with the HMM idea, this would give BCyrius two independent lines of evidence for hybrid structure.

5. Adaptive Total CN Calling
The problem: The GMM parameters are pre-computed and loaded from CYP2D6_gmm.txt:


# gmm.py: means are fixed at 0, 0.5, 1, 1.5, 2, 2.5, ...
# sigma scales with sqrt(CN/2)
This assumes a specific depth distribution that may not match the actual sample. Batch effects, library prep differences, and sequencing platform variations shift the depth distribution.

The fix: For cohort processing, fit the GMM on the fly from the sample's own GC-normalized depth distribution across the ~4000 "norm" control regions in CYP2D6_region_38.bed. The norm regions provide the expected per-bin depth for CN=2. Sample-specific fitting would handle:

Different mean depths (20x vs 40x vs 60x)
Different depth variance (PCR-free vs PCR libraries)
Batch-specific GC bias profiles
Impact: More accurate CN calling at the extremes (CN=0 and CN≥5) where the pre-built GMM's tails may not match reality. The BCyrius paper reports LowQ_high_CN filter for CN≥6 — adaptive fitting could push this threshold higher.

6. Strand Bias Filtering for All Variant Sites
The problem: BCyrius applies strand bias testing (Fisher's exact) only to NOISY_VAR sites:


if total_var > 0 and var_list[i] in NOISY_VAR:
    oddsratio, pvalue = fisher_exact(
        [[forward, reverse], [ntotal / 2, ntotal / 2]]
    )
    if pvalue < P_CUTOFF or forward <= 1 or reverse <= 1:
        total_var = 0
Other variant sites, including unlisted ones, get no strand bias check. But alignment artifacts (from paralogs or repetitive sequence) often show strong strand bias because mismapped reads tend to come preferentially from one orientation.

The fix: Apply strand bias testing to all variant sites, with a lenient threshold for CLEAN_VAR and stricter for unknowns. The infrastructure already exists — alt_forward and alt_reverse are already tracked.

Impact: Low effort, catches a specific class of false positive that the current code misses.

7. Fuzzy Star Allele Matching with Scoring
The problem: get_star() does exact dictionary lookup:


var_list = "_".join(sorted(var_observed))
if var_list not in dic:
    match_tag = "no_match"
If a sample has one extra variant or is missing one variant compared to a known star allele, it's no_match. This is brittle for:

Novel suballeles (e.g., a known *4 plus one additional variant)
Low-coverage samples where one variant was missed
Samples with a variant at a noisy site that was erroneously filtered
The fix (borrowing from StellarPGx's get_backgroud_alleles() concept, but better):


# Score each star allele combination by overlap with observed variants
for candidate in dic:
    candidate_vars = set(candidate.split("_"))
    observed = set(var_observed)
    score = len(candidate_vars & observed) - 0.5 * len(candidate_vars - observed)
    # Penalize for unexpected extra variants in sample
    score -= 0.3 * len(observed - candidate_vars)
Report the best match with a confidence indicator. When the top match scores well but isn't exact, report it as likely_match rather than no_match.

Impact: Directly addresses the "no call" rate. Combined with novel variant annotation (from the earlier discussion), a report like "likely *1/*4 with novel missense at 42128945" is clinically actionable.

8. Cohort-Level Outlier Detection
The problem: BCyrius processes each sample independently. When processing a cohort, it misses the opportunity to flag samples that are outliers relative to the batch.

The fix: After individual calling, compare each sample's raw metrics against the cohort:

Depth distribution at differentiating sites
Fraction of sites with ambiguous CN calls (many None values)
Consistency of CN calls across adjacent sites (a sample with many CN "jumps" may have alignment issues)
Flag samples that deviate significantly for manual review.

Impact: Quality control for large-scale studies. Catches systematic issues (wrong reference genome, truncated BAM, contamination) that BCyrius currently can't detect.

Prioritized Summary
Improvement	Impact on Accuracy	Effort	Failure Mode Addressed
Base-quality-aware likelihoods	High	Medium	Low-coverage miscalls, noisy sites
Physical phasing for disambiguation	High	Medium	More_than_one_possible_genotype
Fuzzy matching + background alleles	High	Low	no_match / no-call rate
Per-site reliability weights	Medium	Medium	False hybrid calls from noisy sites
Read pair SV evidence	Medium	High	Novel hybrid detection, breakpoint resolution
Strand bias for all sites	Medium	Low	Alignment artifact false positives
Adaptive GMM fitting	Medium	Medium	High-CN miscalls, batch effects
Cohort outlier detection	Low-Medium	Low	QC for large studies
The three highest-value, most-feasible improvements are: base-quality-aware likelihoods (replaces the weakest assumption in the model), physical phasing (the infrastructure already exists in haplotype.py), and fuzzy matching (low effort, high clinical impact on no-call rate).

---

## BayesPhase Integration for Improved Phasing

Findings from the BayesPhase haplotype phasing project (~/dev/bayesphase), which implements Bayesian cis/trans inference for short-read WGS data. Several components could directly improve cyrius-rs's star allele calling.

### 9. Replace phase_disambiguate with Bayesian Fragment Scoring

The problem: `phase_disambiguate.rs` uses simple read counting — it tallies cis vs trans fragment counts and requires ≥3 total observations with a ≥3x score ratio. All reads are weighted equally regardless of quality.

The fix: BayesPhase computes per-fragment log Bayes factors using a proper error model:

```
P(obs_i, obs_j | cis)  = 0.5 * P(a|alt_i)*P(b|alt_j) + 0.5 * P(a|ref_i)*P(b|ref_j)
P(obs_i, obs_j | trans) = 0.5 * P(a|alt_i)*P(b|ref_j) + 0.5 * P(a|ref_i)*P(b|alt_j)

where P(observed|true) = (1-e) if match, e if mismatch, with e from base quality.
```

Each fragment gets a reliability weight based on:
- Base quality (error probability per observation)
- Mapping quality (MAPQ < 20 triggers low_mapq flag)
- Soft-clip fraction (alignment confidence)

The combined log Bayes factor across all fragments gives a posterior P(cis) that is more principled than a count ratio threshold. In CYP2D6's paralogous context (CYP2D7 is ~97% identical), many reads will have degraded MAPQ — BayesPhase automatically downweights these.

On HG002 chr20, this approach achieves 0.26% switch error rate vs WhatsHap's 0.53% on the same data.

Impact: More accurate phase disambiguation for CN=2 diplotypes. Directly reduces More_than_one_possible_genotype calls. The existing `read_phasing.rs` log-odds scoring would also benefit from the quality-aware error model.

### 10. STR Detection for CYP2D6 Variant Sites

The problem: CYP2D6 contains tandem repeats and the surrounding region is rich in repetitive elements. Heterozygous indels within repeats cause reads from different haplotypes to follow different CIGAR paths, creating phantom SNP mismatches with MAPQ=60. These artifacts produce high-confidence wrong phase calls.

From the BayesPhase benchmark on chr20: 67-68% of all switch errors occur at sites inside tandem repeats. The alignment artifact is systematic — all reads from one haplotype agree on the wrong allele because the aligner chose a different indel placement through the repeat.

The fix: BayesPhase detects STR regions on the fly from the reference FASTA:
- Homopolymers ≥ 6bp
- Dinucleotide repeats ≥ 4 copies
- Trinucleotide+ repeats ≥ 3 copies

Variant pairs where either site falls in a tandem repeat are marked as unresolved (log Bayes factor capped to 0). On chr20 this reduced switch errors from 72 to 23 (bwa-mem2) and 102 to 34 (strobealign).

For cyrius-rs: the 109 differentiating sites between CYP2D6/CYP2D7 could be pre-screened against the reference for STR context. Sites in repeats should be downweighted or excluded from phasing evidence, similar to the per-site reliability weights idea (#2 above).

Impact: Prevents wrong phase calls at STR-affected differentiating sites. Particularly relevant for the exon9 multi-base region which is near repetitive sequence.

### 11. Population Star Allele Frequencies as a Bayesian Prior

The idea: Instead of running SHAPEIT on raw variants (which doesn't work well for CYP2D6 due to structural variation, paralog interference, and CN>2), use population-level star allele frequencies as a prior on diplotype calls.

SHAPEIT and similar statistical phasing tools are poorly suited for CYP2D6 because:
- Reference panels don't represent CYP2D6's structural complexity (deletions, duplications, hybrids)
- CYP2D7 paralog contamination corrupts the variant calls SHAPEIT would operate on
- SHAPEIT assumes diploid genomes; CN>2 is common (~5-10% of individuals)
- The gene is only ~4.4kb — read-backed phasing is more reliable at this scale

But population star allele frequencies are well-characterized and already available in cyrius-rs's definition tables. When read evidence is ambiguous between candidate diplotypes, the prior can tip the balance:

```
P(diplotype | reads) ∝ P(reads | diplotype) × P(diplotype)

Example:
  Ambiguous between *1/*4 (freq ~15% in Europeans) and *36/*10 (freq <0.1%)
  → Prior favors *1/*4 by >100x
  → Combined with even weak read evidence, this resolves the ambiguity
```

Implementation: After the current matching step enumerates candidate diplotypes, score each by:
1. Read-level likelihood (from improved Bayesian phasing, #9 above)
2. Population frequency prior (from star allele frequency tables, optionally population-specific)
3. Report posterior probability for each candidate

This is analogous to BayesPhase's population LD prior but operating at the star allele level rather than individual variants — a much better fit for the structured nature of CYP2D6 variation.

Impact: Resolves ambiguous diplotype calls without external tools. Most valuable for rare alleles where the prior strongly distinguishes candidates. Requires population-stratified frequency tables (already partially available from PharmVar and the original Cyrius training data).

### 12. RNA-seq Evidence for Long-Range Phasing

When RNA-seq data is available (CYP2D6 is expressed in liver), spliced reads can span across the ~4.4kb gene body, providing long-range phase information beyond DNA insert size.

BayesPhase supports RNA fragment linkage with:
- Same Bayesian scoring model as DNA fragments
- ASE (allele-specific expression) adjustment: evidence reduced by 70% when allele imbalance >80%
- RNA window extends up to 500kb vs 2kb for DNA

For cyrius-rs, RNA-seq could link variants at opposite ends of CYP2D6 that no DNA read pair can span, resolving diplotype ambiguities where the defining variants are far apart within the gene.

Impact: Niche but high-value for research cohorts with paired DNA/RNA data (e.g., GTEx, TCGA). Not applicable to clinical WGS-only workflows.

### Integration approach

The cleanest path would be to extract BayesPhase's core scoring functions as a shared library crate:
- `fragment_log_bayes_factor()` — per-fragment likelihood ratio
- `fragment_reliability_weight()` — quality-aware weighting
- `score_dna_fragments()` — combined evidence with flag detection
- `detect_str_sites()` — STR detection from reference

cyrius-rs would depend on this crate, replacing the counting logic in `phase_disambiguate.rs`. This keeps the pipeline self-contained (no external tool invocation), which aligns with cyrius-rs's design choice to not port the external WhatsHap/SHAPEIT integration from the Python version.

---

Sources:

BCyrius GitLab Repository
BCyrius Paper (PMC)
Original Cyrius Paper (PMC)
# BCyrius: An updated version of Cyrius

BCyrius (Cyrius B) is a new version of Cyrius, a tool for genotyping CYP2D6 from whole-genome sequencing data. The original version of Cyrius was developed by Illumina and is available on [GitHub](https://github.com/Illumina/Cyrius), however, it has not received updates since mid-2021. BCyrius includes fixes and updates to the software to allow genotyping of all currently known star alleles. In addition to providing genotypes, it also outputs the activity score and the predicted phenotype based on the latest CPIC data, as well as population frequencies for the detected diplotype and haplotypes. BCyrius works on samples aligned to the hg38 reference genome.

## Running the program

Example of running the program with the minimum parameters:
```bash
python3 star_caller.py --input INPUT_FILE_PATH --outDir OUTPUT_DIRECTORY
```
**Additional Parameters:**

| Parameter         | Description/Example                                                                                                                                         |
|------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--reference`      | For CRAM input, it is recommended to provide the path to the reference FASTA file, e.g. `--reference hg38.fasta`.      |
| `--genome`     | Contig naming for the aligned sample is determined automatically by default (e.g., whether it is "chr22" or "22"). If this does not work or if you prefer to set it manually, specify `--genome 38` when contigs do not contain "chr" and `--genome chr38` if they do. |
| `--id`     | Specify sample ID. By default, the file name is used (e.g., sample123.bam ID is sample123 and therefore results are saved into sample123.tsv and sample123.json files). To override this, use "--id". For example, `--id 123` creates 123.tsv and 123.json. |
| `--threads`  | Use "--threads" to allocate more threads to the tool, where X is the number of threads, e.g. `--threads 2`.                                                                    |
| `--population` | Use `--population` to output population frequencies for the determined diplotype as well as for each haplotype (see an example of the output below).                      |
| `--print`      | Besides saving the results into a TSV file, you can also use `--print` to display the results on the screen/STDOUT after the program has finished.                 |

## Example output
Example of an outputted TSV file. First two rows (header and the result line) is outputted in all cases. The population data is outputted only when `--population` parameter is used, e.g. `python3 star_caller.py --input Test_Sample_123.bam --outDir results/ --population`.
```
Sample                               Genotype    Filter      Activity score      Predicted phenotype
Test_Sample_123                      *1/*10      PASS        1.25                CYP2D6 Normal Metabolizer

Biogeographic group                  *1/*10      *1          *10
African American/Afro-Caribbean      1.612%      21.116%     3.816%
American                             1.452%      50.08%      1.45%
Central/South Asian                  4.433%      29.324%     7.559%
East Asian                           22.09%      25.779%     42.845%
European                             0.896%      28.501%     1.571%
Latino                               1.917%      36.455%     2.629%
Near Eastern                         3.368%      24.876%     6.769%
Oceanian                             7.05%       61.702%     5.713%
Sub-Saharan African                  0.401%      4.116%      4.869%
```

### Interpreting the output
| Fields in tsv      | Explanation                                                    |
|:-------------------|:---------------------------------------------------------------|
| Sample             | Sample name                                                    |
| Genotype           | Genotype call                                                  |   
| Filter             | Filters on the genotype call                                   |   
| Activity score     | Activity score from the CPIC table                             |   
| Predicted phenotype| Phenotype from the CPIC table                                  |   

A genotype of "None" indicates a no-call.  
There are currently four possible values for the Filter column:  
-`PASS`: a passing, confident call.   
-`More_than_one_possible_genotype`: In rare cases, Cyrius reports two possible genotypes for which it cannot distinguish one from the other. These are different sets of star alleles that result in the same set of variants that cannot be phased with short reads, e.g. \*1/\*46 and \*43/\*45. The two possible genotypes are reported together, separated by a semicolon.   
-`Not_assigned_to_haplotypes`: In a very small portion of samples with more than two copies of CYP2D6, Cyrius calls a set of star alleles but they can be assigned to haplotypes in more than one way. Cyrius reports the star alleles joined by underscores. For example, \*1_\*2_\*68 is reported and the actual genotype could be \*1+\*68/\*2, \*2+\*68/\*1 or \*1+\*2/\*68.  
-`LowQ_high_CN`: In rare cases, at high copy number (>=6 copies of CYP2D6), Cyrius uses less strict approximation in calling copy numbers to account for higher noise in depth and thus the genotype call could be lower confidence than usual.     
  
A .json file is also produced that contains more information about each sample.  

| Fields in json    | Explanation                                                    |
|:------------------|:---------------------------------------------------------------|
| Coverage_MAD      | Median absolute deviation of depth, measure of sample quality  |
| Median_depth      | Sample median depth                                            |
| Total_CN          | Total copy number of CYP2D6+CYP2D7                             |
| Total_CN_raw      | Raw normalized depth of CYP2D6+CYP2D7                          |
| Spacer_CN         | Copy number of CYP2D7 spacer region                            |
| Spacer_CN_raw     | Raw normalized depth of CYP2D7 spacer region                   |
| Variants_called   | Targeted variants called in CYP2D6                             |
| CNV_group         | An identifier for the sample's CNV/fusion status               |
| Variant_raw_count | Supporting reads for each variant                              |
| Raw_star_allele   | Raw star allele call                                           |
| d67_snp_call      | CYP2D6 copy number call at CYP2D6/7 differentiating sites      |
| d67_snp_raw       | Raw CYP2D6 copy number at CYP2D6/7 differentiating sites       |


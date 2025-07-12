# BCyrius: An upgraded version of Cyrius

BCyrius (Cyrius B) is a new version of Cyrius, a tool for genotyping CYP2D6 from whole-genome sequencing data. The original version of Cyrius was developed by Illumina and is available on [GitHub](https://github.com/Illumina/Cyrius), however, it has not received updates since mid-2021. BCyrius includes fixes and updates to the software to allow genotyping of all currently known star alleles. In addition to providing genotypes, it also outputs the activity score and the predicted phenotype based on the latest CPIC data, as well as population frequencies for the detected diplotype and haplotypes. BCyrius works on samples aligned to the hg38 reference genome.

## Installation
To install BCyrius as a Python package, run the following from the root of the repository:
```bash
pip3 install .
```
This will install all required dependencies and make the bcyrius command-line interface available system-wide.


## Running the program

Once installed, run the tool using the bcyrius CLI:
```bash
bcyrius --input INPUT_FILE_PATH --outDir OUTPUT_DIRECTORY
```
**Additional Parameters:**

| Parameter         | Description/Example                                                                                                                                         |
|------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--reference`      | For CRAM input, it is recommended to provide the path to the reference FASTA file, e.g. `--reference hg38.fasta`.      |
| `--genome`     | Contig naming for the aligned sample is determined automatically by default (e.g., whether it is "chr22" or "22"). If this does not work or if you prefer to set it manually, specify `--genome 38` when contigs do not contain "chr" and `--genome chr38` if they do. |
| `--id`     | Specify sample ID. By default, the file name is used (e.g., sample123.bam ID is sample123 and therefore results are saved into sample123.tsv and sample123.json files). To override this, use "--id". For example, `--id 123` creates 123.tsv and 123.json. |
| `--threads`  | Use "--threads" to allocate more threads to the tool, where X is the number of threads, e.g. `--threads 2`.                                                                    |
| `--population-info` | Use `--population` to output population frequencies for the determined diplotype as well as for each haplotype (see an example of the output below).                      |
| `--haplotype-info` | Use `---haplotype-info` to output information for each haplotype (activity value, function and evidence).                      |
| `--print`      | Besides saving the results into a TSV file, you can also use `--print` to display the results on the screen/STDOUT after the program has finished.                 |

## Legacy usage (not recommended):

If you prefer or need to run the script directly without installation:
```bash
python3 star_caller/star_caller.py --input INPUT_FILE_PATH --outDir OUTPUT_DIRECTORY
```

## Example output
Example of an outputted TSV file. First two rows (header and the result line) is outputted in all cases. The haplotype information is only outputted when `--haplotype-info` parameter is used and population data (last part) only when `--population-info` is used, e.g. to print both use `python3 star_caller.py --input ERR1955391.bam --outDir results/ --population-info --haplotype-info`.
```
Sample     	 Genotype 	 Filter 	 Activity score 	 Predicted phenotype
ERR1955391 	 *4/*35   	 PASS   	 1.0            	 Intermediate Metabolizer

Haplotype 	 Activity value 	 Function        	 Evidence strength 	 Evidence summary
*4        	 0.0            	 No function     	 Strong            	 CYP2D6*4 is assigned no function based on strong evidence. CYP2D6*4 has consistently been described in subjects demonstrating decreased metabolism of various CYP2D6 substrates (2211621, 11266079, 1978251, 1978565). Two in vitro studies found CYP2D6*4 had undetectable protein expression (2211621, 11266079). CYP2D6*4 is defined by a splicing defect resulting in a nonfunctional protein. Therefore, consensus among experts was no function with an activity value of 0 based on strong evidence.
*35       	 1.0            	 Normal function 	 Moderate          	 CYP2D6*35 is assigned normal function based on moderate evidence in heterozygous subjects. CYP2D6*35 was first identified in a subject with CYP2D6*5/*35 genotype demonstrating similar dextromethorphan metabolism compared to CYP2D6 *1/*5 subjects, indicating normal function of the CYP2D6*35 allele (9241659). Similarly, a subject with CYP2D6*4/*35 genotype demonstrated decreased dextromethorphan metabolism compared to wildtype subjects (21833166). Additionally, a study of 396 subjects found the -1584C>G substitution, which is found in CYP2D6*35 and other normal function alleles such as CYP2D6*2, was not observed in any subject demonstrating poor dextromethorphan metabolism, indicating normal function of the CYP2D6*35 allele (12766015). Although two in vitro studies found CYP2D6*35 had decreased enzyme activity compared to wildtype (24647041, 30366777), experts did not find these results convincing as the expression system may not be reflective of in vivo environment as demonstrated by the numerous subjects carrying the CYP2D6*35 allele with substrate metabolism similar to that of wildtype as described previously. Therefore, consensus among experts was normal function with an activity value of 1 based on moderate evidence.

                                	 *4/*35 	 *4      	 *35
African American/Afro-Caribbean 	 0.085% 	 4.811%  	 0.887%
American                        	 0.231% 	 10.194% 	 1.134%
Central/South Asian             	 0.148% 	 8.957%  	 0.827%
East Asian                      	 0.001% 	 0.527%  	 0.05%
European                        	 2.022% 	 18.485% 	 5.468%
Latino                          	 0.641% 	 12.051% 	 2.658%
Near Eastern                    	 0.506% 	 11.406% 	 2.219%
Oceanian                        	 0.013% 	 1.78%   	 0.369%
Sub-Saharan African             	 0.0%   	 2.873%  	 0.0%

```

To output the same information on screen when running the command, use `--print`, which outputs it in the following format:

```
======================================================
                     RESULTS
======================================================
        
Sample               → ERR1955391
Genotype             → *4/*35
Filter               → PASS
Activity score       → 1.0
Predicted phenotype  → Intermediate Metabolizer

--------- Details for each haplotype ---------
*4
Activity value       → 0.0
Function             → No function
Evidence strength    → Strong
Evidence summary     → CYP2D6*4 is assigned no function based on strong evidence. CYP2D6*4 has consistently been described in subjects demonstrating decreased metabolism of various CYP2D6 substrates (2211621, 11266079, 1978251, 1978565). Two in vitro studies found CYP2D6*4 had undetectable protein expression (2211621, 11266079). CYP2D6*4 is defined by a splicing defect resulting in a nonfunctional protein. Therefore, consensus among experts was no function with an activity value of 0 based on strong evidence.

*35
Activity value       → 1.0
Function             → Normal function
Evidence strength    → Moderate
Evidence summary     → CYP2D6*35 is assigned normal function based on moderate evidence in heterozygous subjects. CYP2D6*35 was first identified in a subject with CYP2D6*5/*35 genotype demonstrating similar dextromethorphan metabolism compared to CYP2D6 *1/*5 subjects, indicating normal function of the CYP2D6*35 allele (9241659). Similarly, a subject with CYP2D6*4/*35 genotype demonstrated decreased dextromethorphan metabolism compared to wildtype subjects (21833166). Additionally, a study of 396 subjects found the -1584C>G substitution, which is found in CYP2D6*35 and other normal function alleles such as CYP2D6*2, was not observed in any subject demonstrating poor dextromethorphan metabolism, indicating normal function of the CYP2D6*35 allele (12766015). Although two in vitro studies found CYP2D6*35 had decreased enzyme activity compared to wildtype (24647041, 30366777), experts did not find these results convincing as the expression system may not be reflective of in vivo environment as demonstrated by the numerous subjects carrying the CYP2D6*35 allele with substrate metabolism similar to that of wildtype as described previously. Therefore, consensus among experts was normal function with an activity value of 1 based on moderate evidence.

---------------- Population frequencies ----------------
            Biogeographic group  *4/*35       *4     *35
African American/Afro-Caribbean  0.085%   4.811%  0.887%
                       American  0.231%  10.194%  1.134%
            Central/South Asian  0.148%   8.957%  0.827%
                     East Asian  0.001%   0.527%   0.05%
                       European  2.022%  18.485%  5.468%
                         Latino  0.641%  12.051%  2.658%
                   Near Eastern  0.506%  11.406%  2.219%
                       Oceanian  0.013%    1.78%  0.369%
            Sub-Saharan African    0.0%   2.873%    0.0%
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


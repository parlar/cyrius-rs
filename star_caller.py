#!/usr/bin/env python3
#
# Cyrius: CYP2D6 genotyper
# Copyright (c) 2019-2020 Illumina, Inc.
#
# Author: Xiao Chen <xchen2@illumina.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ----------------------------------------------------------------------------
#
# BCyrius: CYP2D6 genotyper (upgraded version of Cyrius)
# Copyright (c) 2024 Andreas Halman
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import os
import sys
import csv
import argparse
import json
import logging
import datetime
from collections import namedtuple, OrderedDict
import pysam
import pandas as pd
import regex


from depth_calling.snp_count import (
    get_supporting_reads,
    get_supporting_reads_single_region,
    get_fraction,
    get_snp_position,
)
from depth_calling.gmm import Gmm
from depth_calling.utilities import (
    parse_gmm_file,
    parse_region_file,
    open_alignment_file,
)
from depth_calling.bin_count import (
    get_normed_depth,
    get_normed_depth_from_count,
    get_read_length,
)
from caller.call_variants import (
    NOISY_VAR,
    call_cn_snp,
    call_cn_var,
    call_cn_var_homo,
    get_allele_counts_var42128936,
    update_var42128936,
    get_called_variants,
    call_exon9gc,
    call_var42126938,
    call_var42127526_var42127556,
    call_var42127803hap,
    call_var42130655insA,
)
from caller.cnv_hybrid import get_cnvtag
from caller.construct_star_table import get_hap_table
from caller.match_star_allele import match_star

TOOL_NAME = "BCyrius"
TOOL_VERSION = "1.0.2"
MAD_THRESHOLD = 0.11
EXON9_SITE1 = 7
EXON9_SITE2 = 8
HIGH_CN_DEPTH_THRESHOLD = 7.5
HAPLOTYPE_VAR = ["g.42126938C>T", "g.42127803C>T", "g.42127526C>T_g.42127556T>C", "g.42130655-42130656insA"]
resource_info = namedtuple(
    "resource_info",
    "genome gmm_parameter region_dic snp_db var_db var_homo_db haplotype_db var_list star_combinations",
)
exon9_values = namedtuple(
    "exon9_values", "exon9_cn exon9cn_in_consensus exon9_raw_site1 exon9_raw_site2"
)
# Below are the SV configurations that the caller is able to call
CNV_ACCEPTED = [
    "star5_star5",
    "star13_star13",
    "star13intron1_star13intron1",
    "star5",
    "star13",
    "star13intron1",
    "star5_star5_star68",
    "star5_star68",
    "cn2",
    "exon9hyb_star5",
    "dup_star13",
    "dup_star13intron1",
    "star13_star68",
    "cn3",
    "exon9hyb",
    "star68",
    "cn4",
    "exon9hyb_exon9hyb",
    "star68_star68",
    "dup_exon9hyb",
    "dup_star68",
    "exon9hyb_star68",
    "cn5",
    "exon9hyb_exon9hyb_exon9hyb",
    "star68_star68_star68",
    "cn6",
    "exon9hyb_exon9hyb_exon9hyb_exon9hyb",
    "star68_star68_star68_star68",
]


def load_parameters():
    """Return parameters."""
    parser = argparse.ArgumentParser(
        description=f"{TOOL_NAME} v{TOOL_VERSION} - a CYP2D6 genotyping tool"
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Input BAM/CRAM file",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--genome",
        help="Reference genome style, select 'chr38' if contig names start with 'chr' otherwise use '38'. If nothing is selected, it will be detected automatically",
        required=False,
        choices = ["autodetect", "38", "chr38"],
        default="autodetect"
    )
    parser.add_argument("-o", "--outDir", help="Output directory", required=True)
    parser.add_argument("--id", help="Sample ID (output file name)", required=False)
    parser.add_argument(
        "-t",
        "--threads",
        help="Optional, number of threads to use. Default is 1",
        type=int,
        required=False,
        default=1,
    )
    parser.add_argument(
        "--countFilePath", help="Optional path to count files", required=False
    )
    parser.add_argument(
        "-r",
        "--reference",
        help="Optional path to reference fasta file for CRAM",
        required=False,
    )

    parser.add_argument(
        "--population-info",
        help="Export population frequencies information",
        action="store_true",
    )

    parser.add_argument(
        "--haplotype-info",
        help="Output haplotype information",
        action="store_true",
    )

    parser.add_argument(
        "--print",
        help="Print the results on screen",
        action="store_true",
    )

    args = parser.parse_args()

    return args


def d6_star_caller(
    bam, call_parameters, threads, count_file=None, reference_fasta=None, index_name=None
):
    """Return CYP2D6 star allele diplotype calls for each sample."""
    d6_call = namedtuple(
        "d6_call",
        "Coverage_MAD Median_depth Total_CN Spacer_CN Total_CN_raw \
        Spacer_CN_raw Variants_called CNV_group Genotype Filter Raw_star_allele \
        Call_info Exon9_CN CNV_consensus d67_snp_call d67_snp_raw \
        Variant_raw_count",
    )
    # 1. Read counting and normalization
    bamfile = open_alignment_file(bam, reference_fasta, index_filename=index_name)
    if count_file is not None:
        reads = bamfile.fetch()
        read_length = get_read_length(reads)
        normalized_depth = get_normed_depth_from_count(
            count_file, call_parameters.region_dic, read_length
        )
    else:
        normalized_depth = get_normed_depth(
            bam, call_parameters.region_dic, threads, reference=reference_fasta
        )

    # no-call after normalizaton
    if normalized_depth.normalized["d67"] is None:
        sample_call = d6_call(
            normalized_depth.mad,
            normalized_depth.mediandepth,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        return sample_call

    # 2. GMM and CN call
    # There are two regions to call CN based on depth: total CYP2D6+CYP2D7, and CYP2D7 spacer region
    cn_call = namedtuple("cn_call", "d67_cn d67_depth spacer_cn spacer_depth")
    gmm_d67 = Gmm()
    gmm_d67.set_gmm_par(call_parameters.gmm_parameter, "d67")
    gcall_d67 = gmm_d67.gmm_call(normalized_depth.normalized["d67"])
    gmm_spacer = Gmm()
    gmm_spacer.set_gmm_par(call_parameters.gmm_parameter, "spacer")
    gcall_spacer = gmm_spacer.gmm_call(normalized_depth.normalized["spacer"])
    high_cn_low_confidence = False
    if gcall_d67.cn is None and gcall_d67.depth_value > HIGH_CN_DEPTH_THRESHOLD:
        high_cn_low_confidence = True
        raw_cn_call = cn_call(
            int(round(gcall_d67.depth_value)),
            gcall_d67.depth_value,
            gcall_spacer.cn,
            gcall_spacer.depth_value,
        )
    else:
        raw_cn_call = cn_call(
            gcall_d67.cn,
            gcall_d67.depth_value,
            gcall_spacer.cn,
            gcall_spacer.depth_value,
        )

    # 3. Get allele counts at D6/D7 SNP (base difference) sites and target variant sites
    # D6/D7 base difference sites. Get read counts at both D6/D7 positions.
    snp_db = call_parameters.snp_db
    snp_d6, snp_d7 = get_supporting_reads(
        bamfile, snp_db.dsnp1, snp_db.dsnp2, snp_db.nchr, snp_db.dindex
    )
    
    # Variants not in homology regions. Get read counts only at D6 positions.
    var_db = call_parameters.var_db
    var_alt, var_ref, var_alt_forward, var_alt_reverse = get_supporting_reads_single_region(
        bamfile, var_db.dsnp1, var_db.nchr, var_db.dindex
    )
    # Look more carefully for insertions at 42128936 from reads
    var_list = call_parameters.var_list
    ref_read, long_ins_read, short_ins_read = get_allele_counts_var42128936(
        bamfile, call_parameters.genome
    )
    var_alt, var_ref = update_var42128936(
        var_list, var_alt, var_ref, ref_read, long_ins_read, short_ins_read
    )
    # Variants in homology regions. Get read counts at both D6/D7 positions.
    var_homo_db = call_parameters.var_homo_db
    var_homo_alt, var_homo_ref = get_supporting_reads(
        bamfile,
        var_homo_db.dsnp1,
        var_homo_db.dsnp2,
        var_homo_db.nchr,
        var_homo_db.dindex,
    )
    
    # This ordered dictionary is for final reporting.
    raw_count = OrderedDict()
    non_homology_variant_count = len(var_alt)
    for i in range(len(call_parameters.var_list)):
        if i < non_homology_variant_count:
            if var_list[i] in NOISY_VAR:
                raw_count.setdefault(
                    var_list[i],
                    "%i(%i:%i),%i"
                    % (var_alt[i], var_alt_forward[i], var_alt_reverse[i], var_ref[i]),
                )
            else:
                raw_count.setdefault(var_list[i], "%i,%i" % (var_alt[i], var_ref[i]))
        else:
            raw_count.setdefault(
                var_list[i],
                "%i,%i"
                % (
                    var_homo_alt[i - non_homology_variant_count],
                    var_homo_ref[i - non_homology_variant_count],
                ),
            )

    # no-call due to total copy number calling
    if raw_cn_call.d67_cn is None:
        sample_call = d6_call(
            normalized_depth.mad,
            normalized_depth.mediandepth,
            raw_cn_call.d67_cn,
            raw_cn_call.spacer_cn,
            raw_cn_call.d67_depth,
            raw_cn_call.spacer_depth,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            raw_count,
        )
        return sample_call

    # 4. Call CNV and hybrids
    d6_fraction = get_fraction(snp_d6, snp_d7)
    raw_d6_cn = [round(raw_cn_call.d67_cn * a, 3) for a in d6_fraction]
    cn_call_snp = call_cn_snp(raw_cn_call.d67_cn, snp_d6, snp_d7)

    # exon9gc
    exon9gc_call_stringent = call_exon9gc(
        snp_d6[EXON9_SITE1 : EXON9_SITE2 + 1],
        snp_d7[EXON9_SITE1 : EXON9_SITE2 + 1],
        raw_cn_call.d67_cn,
    )
    cnvtag, consensus = get_cnvtag(
        raw_cn_call.d67_cn,
        raw_d6_cn,
        cn_call_snp,
        exon9gc_call_stringent,
        raw_cn_call.spacer_cn,
    )

    # no-call due to CNV group calling
    if cnvtag is None or cnvtag not in CNV_ACCEPTED:
        sample_call = d6_call(
            normalized_depth.mad,
            normalized_depth.mediandepth,
            raw_cn_call.d67_cn,
            raw_cn_call.spacer_cn,
            raw_cn_call.d67_depth,
            raw_cn_call.spacer_depth,
            None,
            cnvtag,
            None,
            None,
            None,
            None,
            exon9gc_call_stringent,
            ",".join(str(a) for a in consensus),
            ",".join(str(a) for a in cn_call_snp),
            ",".join(str(a) for a in raw_d6_cn),
            raw_count,
        )
        return sample_call

    # 5. Call variants
    # homology region
    cn_call_var_homo = call_cn_var_homo(raw_cn_call.d67_cn, var_homo_alt, var_homo_ref)
    # non-homology region
    cn_call_var = call_cn_var(
        cnvtag, var_alt, var_ref, var_alt_forward, var_alt_reverse, var_list, var_db
    )
    # call haplotypes
    haplotype_db = call_parameters.haplotype_db

    site42126938_count, var42126938, var42126938_G_haplotype = call_var42126938(
        bamfile, raw_cn_call.d67_cn, haplotype_db["g.42126938C>T"]
    )
    raw_count.setdefault(
        "g.42126938C>T", "%i,%i" % (site42126938_count[1], site42126938_count[0])
    )

    site42127526_count, site42127556_count, var42127526 = call_var42127526_var42127556(
        bamfile, cnvtag, haplotype_db["g.42127526C>T_g.42127556T>C"]
    )
    raw_count.setdefault(
        "g.42127526C>T", "%i,%i" % (site42127526_count[1], site42127526_count[0])
    )
    raw_count.setdefault(
        "g.42127556T>C", "%i,%i" % (site42127556_count[1], site42127556_count[0])
    )

    var42127803_diff_haplotype = call_var42127803hap(
        bamfile, cnvtag, haplotype_db["g.42127803C>T"]
    )

    var42130655insA_count, var42130655insA = call_var42130655insA(
        bamfile, raw_cn_call.d67_cn, haplotype_db["g.42130655-42130656insA"]
    )
    raw_count.setdefault(
        "g.42130655-42130656insA", "%i,%i" % (var42130655insA_count[1], var42130655insA_count[0])
    )

    # 6. Call star allele
    total_callset = get_called_variants(var_list, cn_call_var)
    called_var_homo = get_called_variants(var_list, cn_call_var_homo, len(cn_call_var))
    total_callset += called_var_homo
    total_callset += var42126938
    total_callset += var42127526
    total_callset += var42130655insA

    star_called = match_star(
        total_callset,
        cnvtag,
        raw_cn_call.spacer_cn,
        call_parameters.star_combinations,
        exon9_values(
            exon9gc_call_stringent,
            consensus.exon9_and_downstream,
            raw_d6_cn[EXON9_SITE1],
            raw_d6_cn[EXON9_SITE2],
        ),
        var42126938_G_haplotype,
        var42127803_diff_haplotype,
    )

    genotype_filter = None
    final_star_allele_call = None
    # no-call due to star allele matching
    if star_called.call_info and star_called.call_info != "no_match":
        final_star_allele_call = star_called.clean_call
        if final_star_allele_call:
            if ";" in final_star_allele_call:
                genotype_filter = "More_than_one_possible_genotype"
            elif "/" not in final_star_allele_call:
                genotype_filter = "Not_assigned_to_haplotypes"
            elif high_cn_low_confidence:
                genotype_filter = "LowQ_high_CN"
            else:
                genotype_filter = "PASS"

    sample_call = d6_call(
        normalized_depth.mad,
        normalized_depth.mediandepth,
        raw_cn_call.d67_cn,
        raw_cn_call.spacer_cn,
        raw_cn_call.d67_depth,
        raw_cn_call.spacer_depth,
        star_called.variants_called.split(),
        cnvtag,
        final_star_allele_call,
        genotype_filter,
        star_called.raw_call,
        star_called.call_info,
        exon9gc_call_stringent,
        ",".join(str(a) for a in consensus),
        ",".join(str(a) for a in cn_call_snp),
        ",".join(str(a) for a in raw_d6_cn),
        raw_count,
    )
    bamfile.close()
    return sample_call


def prepare_resource(datadir, parameters):
    region_file = os.path.join(datadir, "CYP2D6_region_38.bed")
    snp_file = os.path.join(datadir, "CYP2D6_SNP_38.txt")
    gmm_file = os.path.join(datadir, "CYP2D6_gmm.txt")
    star_table = os.path.join(datadir, "star_table.txt")
    variant_file = os.path.join(datadir, "CYP2D6_target_variant_38.txt")
    variant_homology_file = os.path.join(
        datadir, "CYP2D6_target_variant_homology_region_38.txt"
    )
    haplotype_file = os.path.join(datadir, "CYP2D6_haplotype_38.txt")
    star_combinations = get_hap_table(star_table)

    if parameters.genome == "autodetect":
        genome = "chr38" if isChrInChromosomeNames(parameters.input) else "38"
    else:
        genome = parameters.genome

    for required_file in [
        region_file,
        snp_file,
        variant_file,
        variant_homology_file,
        haplotype_file,
        gmm_file,
    ]:
        if os.path.exists(required_file) == 0:
            raise Exception("File %s not found." % required_file)

    snp_db = get_snp_position(snp_file, genome)
    var_db = get_snp_position(variant_file, genome)
    var_homo_db = get_snp_position(variant_homology_file, genome)
    haplotype_db = {}
    for variant in HAPLOTYPE_VAR:
        haplotype_db.setdefault(variant, get_snp_position(haplotype_file, genome, variant))
    var_list = []
    with open(variant_file) as f:
        for line in f:
            if line[0] != "#":
                var_name = line.split()[-1]
                var_list.append(var_name)
    with open(variant_homology_file) as f:
        for line in f:
            if line[0] != "#":
                var_name = line.split()[-1]
                var_list.append(var_name)
    gmm_parameter = parse_gmm_file(gmm_file)
    region_dic = parse_region_file(region_file, genome)
    call_parameters = resource_info(
        genome,
        gmm_parameter,
        region_dic,
        snp_db,
        var_db,
        var_homo_db,
        haplotype_db,
        var_list,
        star_combinations,
    )
    return call_parameters

# Load the CYP2D6 haplotype functionality from the file
def loadHaplotypeFunctionality(file_path):
    haplotype_functionality = {}
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            haplotype = row["Haplotype"].replace("≥", "")
            activity = row["Activity"].replace("≥", "")
            function = row["Function"]
            evidence_strength = row["EvidenceStrength"]
            evidence_summary = row["EvidenceSummary"]
            haplotype_functionality[haplotype] = {
                "activity": activity,
                "function": function,
                "evidence_strength": evidence_strength,
                "evidence_summary": evidence_summary
            }
    return haplotype_functionality

def compareGenotypes(query, haplotype_functionality):
    # Matches haplotype with activity value, functionality, evidence strength, and summary
    haplotypes = []
    diplotypes = query.split("/")
    for haplotype_pair in diplotypes:
        subhaplotypes = haplotype_pair.split("+")
        haplotypes.extend(subhaplotypes)

    haplotype_details = []
    has_na = False  # Flag to track if any haplotype results in "n/a"

    for haplotype in haplotypes:
        haplotype = haplotype.strip()
        if haplotype not in haplotype_functionality:
            has_na = True
            haplotype_details.append({
                "haplotype": haplotype,
                "activity": "n/a",
                "function": "n/a",
                "evidence_strength": "n/a",
                "evidence_summary": "n/a"
            })
        else:
            details = haplotype_functionality[haplotype]
            activity_value = details["activity"]
            if activity_value == "n/a":
                has_na = True
            haplotype_details.append({
                "haplotype": haplotype,
                "activity": activity_value,
                "function": details["function"],
                "evidence_strength": details["evidence_strength"],
                "evidence_summary": details["evidence_summary"]
            })

    return haplotype_details, has_na

def calculateTotalActivityScore(haplotype_details):
    total_activity = 0.0
    for haplotype_info in haplotype_details:
        activity_value = haplotype_info["activity"]
        if activity_value == "n/a":
            return "n/a"  # If any haplotype has "n/a", total activity is "n/a"
        total_activity += float(activity_value.replace("≥", ""))  # Clean the activity value
    return total_activity


def sortGenotype(genotype):
    # Sorts genotype so smaller integer will be the first element
    if not genotype:
        return genotype

    genotype_split = genotype.split("/")
    if len(genotype_split) == 1:
        return genotype

    a1_star = genotype_split[0]
    a2_star = genotype_split[1]

    try:
        if "+" in a1_star:
            a1_int = int(a1_star.split("*")[1].split("+")[0])
        else:
            a1_int = int(a1_star.split("*")[1].split("x")[0])
    except (ValueError, IndexError):
        a1_int = float('inf')  # Assign a very large value or some default

    try:
        if "+" in a2_star:
            a2_int = int(a2_star.split("*")[1].split("+")[0])
        else:
            a2_int = int(a2_star.split("*")[1].split("x")[0])
    except (ValueError, IndexError):
        a2_int = float('inf')  # Assign a very large value or some default

    sorted_genotype = a2_star + "/" + a1_star if a1_int > a2_int else genotype

    return sorted_genotype

def matchPhenotype(phenotypes, genotype, haplotype_functionality):
    predicted = []
    if genotype:
        # Split into possible diplotypes
        diplotypes = genotype.split(";")
        
        for diplotype in diplotypes:
            sorted_diplotype = sortGenotype(diplotype.strip())
            haplotype_details, has_na = compareGenotypes(sorted_diplotype, haplotype_functionality)
            if has_na:
                total_activity = "n/a"
            else:
                total_activity = calculateTotalActivityScore(haplotype_details)
            predicted_phenotype = determinePhenotype(total_activity)
            predicted.append({
                "total_activity": total_activity,
                "predicted_phenotype": predicted_phenotype,
                "haplotype_details": haplotype_details  # Attach haplotype information
            })

    if not predicted:
        return [{'total_activity': 'n/a', 'predicted_phenotype': 'n/a', 'haplotype_details': []}]

    return predicted

def determinePhenotype(activity_score):
    # Function to assign the phenotype based on the activity score
    if activity_score == "n/a":
        return "Indeterminate"
    if activity_score == 0:
        return "Poor Metabolizer"
    elif 0.25 <= activity_score <= 1:
        return "Intermediate Metabolizer"
    elif 1.25 <= activity_score <= 2.25:
        return "Normal Metabolizer"
    elif activity_score > 2.25:
        return "Ultrarapid Metabolizer"
    return "Unknown"

def get_haplotype_percentages(file_path, haplotype_name):
    # Load the tab-delimited file into a DataFrame
    df = pd.read_csv(file_path, delimiter='\t')

    row = df[df['CYP2D6 allele'] == haplotype_name]

    if row.empty:
        return None  # Return None if haplotype is not found

    # Convert the values to percentages and round to one decimal place
    percentages = (row.iloc[0, 1:].astype(float) * 100).round(3)

    # Convert to string and add '%' symbol, but leave empty if NaN
    percentages = percentages.apply(lambda x: f"{x}%" if not pd.isna(x) else "")

    return percentages

def get_diplotype_percentages(file_path, diplotype_name):
    # Load the tab-delimited file into a DataFrame
    df = pd.read_csv(file_path, delimiter='\t')

    # Find the row corresponding to the diplotype name
    row = df[df['CYP2D6 allele'] == diplotype_name]

    if row.empty:
        return None  # Return None if diplotype is not found

    # Convert the values to percentages and round to one decimal place
    percentages = (row.iloc[0, 1:].astype(float) * 100).round(3)

    # Convert to string and add '%' symbol, but leave empty if NaN
    percentages = percentages.apply(lambda x: f"{x}%" if not pd.isna(x) else "")

    return percentages

def diplotype_frequencies(haplotype_file_path, diplotype_file_path, diplotype_string):
    diplotypes = diplotype_string.split(';')
    
    # Dictionary to store frequencies for diplotypes and haplotypes
    frequencies = {}

    for diplotype in diplotypes:
        haplotypes = diplotype.split('/')

        # Get frequencies for the diplotype itself
        frequencies[diplotype] = get_diplotype_percentages(diplotype_file_path, diplotype)
        if frequencies[diplotype] is None:
            return f"Diplotype {diplotype} not found in the file."
        
        # Get frequencies for each haplotype
        for haplotype in haplotypes:
            frequencies[haplotype] = get_haplotype_percentages(haplotype_file_path, haplotype)
            if frequencies[haplotype] is None:
                return f"Haplotype {haplotype} not found in the file."

    frequency_table = pd.DataFrame(frequencies)
    frequency_table.index.name = 'Biogeographic group'
    frequency_table.reset_index(inplace=True)

    return frequency_table

def isChrInChromosomeNames(input_file):
    chr_in_names = False
    input_file_read = pysam.AlignmentFile(input_file, 'rb')

    for contig_id in range(5):  # Check only the first 5 chromosomes to speed it up
        contig_name = input_file_read.get_reference_name(contig_id)

        if "chr" in contig_name:
            chr_in_names = True
            break

    return chr_in_names


def main():
    parameters = load_parameters()
    inputfile = parameters.input
    outdir = parameters.outDir
    sample_id = parameters.id
    reference_fasta = parameters.reference
    threads = parameters.threads
    path_count_file = parameters.countFilePath
    export_population_frequencies = parameters.population_info
    print_results = parameters.print
    output_haplotype_info = parameters.haplotype_info
    logging.basicConfig(level=logging.DEBUG)

    if not sample_id:
        sample_id = os.path.splitext(os.path.basename(inputfile))[0]

    if os.path.exists(outdir) == 0:
        os.makedirs(outdir)

    datadir = os.path.join(os.path.dirname(__file__), "data")

    # Prepare data files
    call_parameters = prepare_resource(datadir, parameters)

    haplotypes_func_file_path = os.path.join(datadir, 'CYP2D6_haplotypes_functionality.txt')
    haplotypes_freq_file_path = os.path.join(datadir, 'CYP2D6_frequency_table_haplotypes.txt')
    diplotypes_freq_file_path = os.path.join(datadir, 'CYP2D6_frequency_table_diplotypes.txt')
    
    out_json = os.path.join(outdir, sample_id + ".json")
    out_tsv = os.path.join(outdir, sample_id + ".tsv")
    final_output = {}

    bam_name = inputfile
    index_name = None
    if '##idx##' in bam_name:
        bam_name, index_name = bam_name.split('##idx##')

    sample_file_name = os.path.basename(bam_name)

    count_file = None
    if path_count_file is not None:
        count_file = os.path.join(path_count_file, sample_id + "_count.txt")
    if "://" not in bam_name and os.path.exists(bam_name) == 0:
        logging.warning("Input file for sample %s does not exist.", sample_id)
    else:
        logging.info(
            "Starting %s v%s", TOOL_NAME, TOOL_VERSION
        )
        logging.info(
            "Processing sample %s (%s) at %s", sample_id, sample_file_name, datetime.datetime.now()
        )
        cyp2d6_call = d6_star_caller(
            bam_name, call_parameters, threads, count_file, reference_fasta, index_name=index_name
        )._asdict()
        # Use normalized coverage MAD across stable regions
        # as a sample QC measure.
        if cyp2d6_call["Coverage_MAD"] > MAD_THRESHOLD:
            logging.warning(
                "Sample %s has uneven coverage. CN calls may be unreliable.",
                sample_id,
            )
        final_output.setdefault(sample_id, cyp2d6_call)


    # Write to json
    logging.info("Writing to json (%s) at %s", out_json, datetime.datetime.now())
    with open(out_json, "w") as json_output:
        json.dump(final_output, json_output)

    # Write to tsv
    logging.info("Writing to tsv (%s) at %s", out_tsv, datetime.datetime.now())
    header = ["Sample", "Genotype", "Filter", "Activity score", "Predicted phenotype"]
    frequency_table = None

    with open(out_tsv, "w") as tsv_output:
        # Write the main header
        tsv_output.write("\t".join(header) + "\n")
        
        for sample_id in final_output:
            final_call = final_output[sample_id]
            sorted_genotype = sortGenotype(final_call["Genotype"])
            haplotype_functionality = loadHaplotypeFunctionality(haplotypes_func_file_path)
            phenotypes = []
            predictions = matchPhenotype(phenotypes, sorted_genotype, haplotype_functionality)

            count_of_diplotypes = sum(1 for diplo in predictions if isinstance(diplo, dict)) 

            # Variables to hold the final outputs
            if predictions[0]:
                if count_of_diplotypes > 1:
                    activity_scores = ";".join([str(p["total_activity"]) for p in predictions])
                    predicted_phenotypes = ";".join([p["predicted_phenotype"] for p in predictions])
                else:
                    activity_scores = predictions[0]["total_activity"]
                    predicted_phenotypes = predictions[0]["predicted_phenotype"]
            else:
                activity_scores = "-"
                predicted_phenotypes = "-"

            haplotype_info_per_solution = []

            for prediction in predictions:
                haplotype_info = []
                for haplotype_detail in prediction["haplotype_details"]:
                    haplotype_info.append({
                        "Haplotype": haplotype_detail["haplotype"],
                        "ActivityValue": haplotype_detail["activity"],
                        "Function": haplotype_detail["function"],
                        "EvidenceStrength": haplotype_detail["evidence_strength"],
                        "EvidenceSummary": haplotype_detail["evidence_summary"]
                    })
                haplotype_info_per_solution.append(haplotype_info)

            # Writing output_per_sample to the file
            output_per_sample = [
                sample_id,
                sorted_genotype,
                final_call["Filter"],
                activity_scores,
                predicted_phenotypes
            ]
            
            tsv_output.write("\t".join(str(a) for a in output_per_sample) + "\n")
            
            # Check if haplotype_info is requested and write haplotype information
            if output_haplotype_info and sorted_genotype != "None":
                # Write header for haplotype details
                tsv_output.write("\nHaplotype\tActivity value\tFunction\tEvidence strength\tEvidence summary\n")
                
                for solution_info in haplotype_info_per_solution:
                    for haplotype in solution_info:
                        # Write each haplotype detail as a row
                        haplotype_row = [
                            haplotype["Haplotype"],
                            haplotype["ActivityValue"],
                            haplotype["Function"],
                            haplotype["EvidenceStrength"],
                            haplotype["EvidenceSummary"]
                        ]
                        tsv_output.write("\t".join(str(h) for h in haplotype_row) + "\n")
                
            # Write population frequencies if requested and sorted_genotype is available
            if export_population_frequencies and sorted_genotype != "None":
                tsv_output.write("\n")
                frequency_table = diplotype_frequencies(haplotypes_freq_file_path, diplotypes_freq_file_path, sorted_genotype)
                if not isinstance(frequency_table, str):
                    frequency_table.to_csv(tsv_output, sep='\t', index=False)

    if print_results:
        column_width = 20
        print("""
========================================================
                     RESULTS
========================================================
        """)

        for header_item, result_item in zip(header, output_per_sample):
            print(f"{str(header_item):<{column_width}} → {str(result_item)}")

        print("")
        processed_haplotypes = set()

        if output_haplotype_info and sorted_genotype != "None":
            # Example: How you might print or process the haplotype info
            for index, solution_info in enumerate(haplotype_info_per_solution):
                print(f"-------------- Details for each haplotype --------------")
                if index > 0:
                    print(f"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                    print(f"Alternative diplotype solution:")

                for haplotype in solution_info:
                    if haplotype['Haplotype'] in processed_haplotypes:
                        continue  # Skip this haplotype if already handled

                    # If it's a new haplotype, print the details
                    print(f"{haplotype['Haplotype']}")
                    print(f"{'Activity value':<{column_width}} → {haplotype['ActivityValue']}")
                    print(f"{'Function':<{column_width}} → {haplotype['Function']}")
                    if haplotype['EvidenceStrength'] != "n/a":
                        print(f"{'Evidence strength':<{column_width}} → {haplotype['EvidenceStrength']}")
                        print(f"{'Evidence summary':<{column_width}} → {haplotype['EvidenceSummary']}")
                    print("")  # For spacing between haplotypes

                    # Add the haplotype to the set of processed ones
                    processed_haplotypes.add(haplotype['Haplotype'])

        if export_population_frequencies and not isinstance(frequency_table, str):
            print(f"---------------- Population frequencies ----------------")
            print(frequency_table.to_string(index=False))




if __name__ == "__main__":
    main()

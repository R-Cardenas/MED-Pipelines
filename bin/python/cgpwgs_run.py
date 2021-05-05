#!/usr/bin/env python
# RPC 161220

import glob
import argparse
import os

# Python script that will run the cgpwgs
# Script was created as cgpwgs (and not cgpwxs) does not like wildcards
# Script will take normal and tumor names and create the whole filename for input

#############
# Functions #
#############

def remove_empties_list(list1):
    x = list(filter(None, list1))
    y = str(x).strip('[').strip(']').strip("'") # remnants removed from list conversion to string
    return y


################################################
# set args parse to allow input from terminal  #
################################################

parser = argparse.ArgumentParser()

parser.add_argument('--tumor', required=True)
parser.add_argument('--normal', required=True)
parser.add_argument('--assembly', required=True)
parser.add_argument('--home', required=True)
parser.add_argument('--reference', required=True)
parser.add_argument('--annotation', required=True)
parser.add_argument('--snv_indel', required=True)
parser.add_argument('--cgp_QC', required=True)
parser.add_argument('--cgp_CNV', required=True)

args = parser.parse_args()
tumor1 = args.tumor
normal1 = args.normal
assembly1 = args.assembly
home1 = args.home
ref1 = args.reference
ann1 = args.annotation
snv_indel1 = args.snv_indel
cgp_QC1 = args.cgp_QC
cgp_CNV1 = args.cgp_CNV

################################################
# use input and read file location     #
################################################

# Possible file locations for bams
locations = ["/input/","/output/BAM/merged_lanes/"]


# Gather file locations
normal_bam = remove_empties_list([glob.glob(home1 + e + normal1 + "*.bam") for e in locations])
normal_bai = remove_empties_list([glob.glob(home1+  e + normal1 + "*.bai") for e in locations])
tumor_bam = remove_empties_list([glob.glob(home1 + e + tumor1 + "*.bam") for e in locations])
tumor_bai = remove_empties_list([glob.glob(home1 + e + tumor1 + "*.bai") for e in locations])


#####################################################
# The singularity / cgpwxs command with file inputs #
#####################################################



cmd1 = fr"""singularity exec --cleanenv \
--home {home1} \
--bind /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes:/var/spool/ref:ro \
--bind {home1}/input:/var/spool/data:ro \
/gpfs/afm/cg_pipelines/Pipelines/singularity/images/cgpwgs_2.1.0.img \
ds-cgpwgs.pl \
-reference {ref1} \
-annot {ann1} \
-snv_indel {snv_indel1} \
-tumour {tumor_bam} \
-tidx {tumor_bai} \
-normal {normal_bam} \
-nidx {normal_bai} \
-exclude NC_007605,hs37d5,GL% \
-outdir {home1}/output/cgpwxs/{tumor1}_vs_{normal1} \
-sp "Homo sapiens" \
-assembly "GRCh37"
"""



cmd2 = fr"""singularity exec --cleanenv \
--home {home1} \
--bind /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes:/var/spool/ref:ro \
--bind {home1}/input:/var/spool/data:ro \
/gpfs/afm/cg_pipelines/Pipelines/singularity/images/cgpwgs_2.1.0.img \
 ds-cgpwgs.pl \
-r {ref1} \
-a {ann1}  \
-si {snv_indel1} \
-cs  {cgp_CNV1} \
-qc {cgp_QC1} \
-pl 3.65 -pu 1.0 \
-e 'MT,GL%,hs37d5,NC_007605' \
-tumour {tumor_bam} \
-tidx {tumor_bai} \
-normal {normal_bam} \
-nidx {normal_bai} \
-outdir {home1}/output/cgpwgs22/{tumor1}_vs_{normal1} \
-sp "Homo sapiens" \
-assembly "GRCh37"
"""
print(cmd2)
os.system(cmd2)

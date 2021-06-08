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

# function will remove the empties from a list, as well as the unwanted punctuation

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
parser.add_argument('--home', required=True)
parser.add_argument('--reference', required=True)
parser.add_argument('--annotation', required=True)
parser.add_argument('--snv_indel', required=True)


args = parser.parse_args()
tumor1 = args.tumor
normal1 = args.normal
home1 = args.home
ref1 = args.reference
ann1 = args.annotation
snv_indel1 = args.snv_indel

################################################
# use input and read file location     #
################################################

# Possible file locations for bams
locations = ["/input/","/output/BAM/merge/RMD/"]


# Gather file locations
normal_bam = os.path.basename(remove_empties_list([glob.glob(home1 + e + normal1 + "*.bam") for e in locations]))
normal_bai = os.path.basename(remove_empties_list([glob.glob(home1+  e + normal1 + "*.bai") for e in locations]))
tumor_bam = os.path.basename(remove_empties_list([glob.glob(home1 + e + tumor1 + "*.bam") for e in locations]))
tumor_bai = os.path.basename(remove_empties_list([glob.glob(home1 + e + tumor1 + "*.bai") for e in locations]))


#####################################################
# The singularity / cgpwxs command with file inputs #
#####################################################



cmd2 = fr"""
	rm -fr {home1}/output/cgpwxs/{tumor1}_vs_{normal1} # during retry if file exist causes fail
    mkdir -p {home1}/output/cgpwxs/{tumor1}_vs_{normal1}

	singularity exec --cleanenv \
    --home {home1} \
    --bind /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes:/var/spool/mail:ro \
    --bind {home1}/output/BAM/merge/RMD:/var/spool/data:ro \
	/gpfs/afm/cg_pipelines/Pipelines/singularity/images/cgpwxs_3.1.6.img \
    ds-cgpwxs.pl \
    -reference {ref1} \
    -annot {ann1} \
    -snv_indel {snv_indel1} \
    -tumour /var/spool/data/{tumor_bam} \
    -tidx /var/spool/data/{tumor_bai} \
    -normal /var/spool/data/{normal_bam} \
    -nidx /var/spool/data/{normal_bai} \
    -exclude NC_007605,hs37d5,GL% \
    -outdir {home1}/output/cgpwxs/{tumor1}_vs_{normal1} \
    -sp "Human" \
    -assembly "GRCh37"
	"""

print(cmd2)
os.system(cmd2)

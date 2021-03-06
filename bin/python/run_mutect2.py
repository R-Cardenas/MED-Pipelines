#!/usr/bin/env python
# RPC 161220

import glob
import argparse
import os
import subprocess

# Script will take the tumor and normal values
# Find files using glob and run the script

#Functions:
def strip(x):
    y = str(x).strip('[').strip(']').strip("'") # removes square brackets and "'" from string
    return y


################################################
# set args parse to allow input from terminal  #
################################################

parser = argparse.ArgumentParser()

parser.add_argument('--tumor', required=True)
parser.add_argument('--normal', required=True)
parser.add_argument('--germline', required=True)
parser.add_argument('--PoN', required=True)
parser.add_argument('--ref', required=True)
parser.add_argument('--intervals', required=True)

args = parser.parse_args()
tumor1 = args.tumor
normal1 = args.normal
germline1 = args.germline
pon1 = args.PoN
ref1 = args.ref
intervals1 = args.intervals


##########################
# Obtain globs for files  #
##########################

tumor2 = tumor1.split("-")
print(tumor2)
tumor3 = tumor2[0] + '-' + tumor2[1] # take sample name (eg sample1) and descriptor (eg normal or tumor )
tumor_path = strip(glob.glob(tumor3 + '*')) # find path of bam file using sample-tumorID (e.g. S01-TUMOR*)


normal2 = normal1.split("-")
normal3 = normal2[0] + '-' + normal2[1] # take sample name (eg sample1) and descriptor (eg normal or normal )
normal_path = strip(glob.glob(normal3 + '*'))

# check samples are the same
if normal[0] == tumor[0]:
    pass
else:
    SyntaxError("Sample names for Tumor and Normal samples are different please revise")
# create output name
final_name = tumor[0] + '-unfiltered-GATK.vcf.gz'

##########################
# Mutect2 command to run #
##########################

cmd = fr"""
mkdir -p tmp
gatk BuildBamIndex \
-I {normal_path} \
-O {normal_path}.bai \
--TMP_DIR tmp

mkdir -p tmp
gatk BuildBamIndex \
-I {tumor_path} \
-O {tumor_path}.bai \
--TMP_DIR tmp

gatk Mutect2 \
-R {ref1} \
-I {normal_path} \
-I {tumor_path} \
-normal {normal_path} \
--germline-resource {germline1} \
--panel-of-normals {pon1} \
-O ${final_name} \
-L {intervals1}
"""

############
### RUN ###
###########
print(cmd)
p = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
print(p.stdout)

# RPC 210920
# This will merge sample bam files from different variant callers - SNPs only!
# !/usr/bin/env python
import re
import argparse
import pathlib
import os

# Test input below:
files = 'sample2-fam2-GATK.vcf.gz sample2-fam2-freebayes.vcf.gz sample1-fam1-freebayes.vcf.gz sample1-fam1-GATK.vcf.gz'

#############
# ARG PARSE #
#############

# Split the input from nextflow (space delim)
# parser = argparse.ArgumentParser()
# parser.add_argument('--bam', required=True)
#
# args = parser.parse_args()
# files = args.bam
files2 = files.split(" ")


###############
# SAMPLE NAME #
###############

# Extract samples names from each input and create unique list
bam_samples = list()
for f in files2:
    Alist = f.split("-")
    sample = Alist[0]
    bam_samples.append(sample)

unique = set(bam_samples)


###############
# BED TOOLS   #
###############

# This will search for files with the same samples name
# Determine how many there are
# Bedtools only accepts 2 files at a time (-a and -b flags)
# If >2 the proccess will have to be repeated with the first intersected file.

for f in unique:

    # --- checks number of files --
    regex = f + "-"
    regex2 = re.compile(fr'{regex}') # searched for 'samplename-' hyphen needed for end of samples or S2 and S22 would be mixed (e.g.).
    selected_files = list(filter(regex2.search, files2)) #searches how many files with same sample name
    print(regex)
    print(selected_files)
    if len(selected_files) == 0:
        print('no files detected')
    elif len(selected_files) == 1:
        print('only one file found using sample names - need 2 minimum for intersection')
    elif len(selected_files) == 2:
        print('selected files = 2')

        # We need to reconstruct the name again for merge_family_indel.py - see line 71
        outputname = selected_files[0].split("-")
        outputname2 = outputname[0] + "-" + outputname[1] + "-" + "caller.merged.indels.vcf"

        # --- will regex specifically GATK or freebayes
        regex_freebayes = f + "(.*)freebayes.vcf"
        regex2_freebayes = re.compile(fr'{regex_freebayes}') # searched for 'samplename-' hyphen needed for end of samples or S2 and S22 would be mixed (e.g.).
        freebayes_files = list(filter(regex2_freebayes.search, selected_files)) #searches how many files with same sample name
        print(freebayes_files)
        regex_GATK = f + "(.*)GATK.vcf"
        regex2_GATK = re.compile(fr'{regex_GATK}') # searched for 'samplename-' hyphen needed for end of samples or S2 and S22 would be mixed (e.g.).
        GATK_files = list(filter(regex2_GATK.search, selected_files)) #searches how many files with same sample name
        print(GATK_files)
        cmd = f"""
        bedtools intersect \
        -a {GATK_files[0]} \
        -b {freebayes_files[0]} \
        -f 0.10 -F 0.10 -wa -header \
        > {outputname2} """
        print(cmd)
        os.system(cmd)

    elif len(selected_files) > 2:
        SyntaxError("Why is there more than 2 files. Should just be freebayes and GATK")


print('Merging of vcf caller files has finished')

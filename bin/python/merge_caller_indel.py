# RPC 210920
# This will merge sample bam files from different variant callers - SNPs only!
# !/usr/bin/env python
import re
import argparse
import pathlib
import os

# Test input below:
files = 'sample1-fam1-GATK.vcf.gz sample1-fam1-freebayes.vcf.gz sample2-fam2-freebayes.vcf.gz sample1-fam1-freebabes.vcf.gz'

#############
# ARG PARSE #
#############

# Split the input from nextflow (space delim)
parser = argparse.ArgumentParser()
parser.add_argument('--bam', required=True)

args = parser.parse_args()
files = args.bam
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
    regex = f + "-"
    regex2 = re.compile(fr'{regex}') # searched for 'samplename-' hyphen needed for end of samples or S2 and S22 would be mixed (e.g.).
    selected_files = list(filter(regex2.search, files2)) #searches how many files with same sample name

    if len(selected_files) == 0:
        print('no files detected')
    elif len(selected_files) == 1:
        print('only one file found using sample names - need 2 minimum for intersection')
    elif len(selected_files) == 2:
        print('selected files = 2')

        # We need to reconstruct the name again for merge_family_indel.py - see line 71
        outputname = selected_files[0].split("-")
        outputname2 = outputname[0] + "-" + outputname[1] + "-" + "caller.merged.indels.vcf"

        cmd = f"""
        bedtools intersect \
        -a {selected_files[1]} \
        -b {selected_files[0]} \
        -f 0.10 -F 0.10 -wa -header \
        > {outputname2} """
        print(cmd)
        os.system(cmd)

    elif len(selected_files) > 2:
        print('selected files > 2')

        # We need to reconstruct the name again for merge_family_indel.py - see line 71
        outputname = selected_files[0].split("-")
        outputname2 = outputname[0] + "-" + outputname[1] + "-" + "caller.merged.indels.vcf"

        file = pathlib.Path(f'{outputname2}')

        def loop1():
            if file.exists():
                print('file exists')
                for z in range(2,len(selected_files)): # notice that '-a' is now the first intersection line 59-67
                    cmd = f"""
                    bedtools intersect \
                    -a {outputname2} \
                    -b {selected_files[z]} \
                    -f 0.10 -F 0.10 -wa -header \
                    > {outputname2}.tmp

                    mv {outputname2}.tmp {outputname2}""" # bedtools needs another name otherwise results in empty file
                    print(cmd)
                    os.system(cmd)
            else:
                cmd = f"""
                bedtools intersect \
                -a {selected_files[1]} \
                -b {selected_files[0]} \
                -f 0.10 -F 0.10 -wa -header \
                > {outputname2} """
                print(cmd)
                os.system(cmd)
                loop1()


        loop1()


print('Merging of vcf caller files has finished')

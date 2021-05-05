# RPC 210920
# This will merge sample bam files from different variant callers - SNPs only!
# !/usr/bin/env python
import os
import argparse
import subprocess
import re
import pathlib

#files = 'sample1-fam1-GATK.vcf.gz sample1-fam1-freebayes.vcf.gz'

#############
# ARG PARSE #
#############

# Split the input from nextflow (space delim)
parser = argparse.ArgumentParser()
parser.add_argument('--bam', required=True)

args = parser.parse_args()
files = args.bam
files2 = files.split(" ")


# Extract samples names from each input and create unique list
bam_samples = list()
for f in files2:
    Alist = f.split("-")
    sample = Alist[1]
    bam_samples.append(sample)

unique = set(bam_samples)
print(unique)

# For loop submits each sample to be lane merged
for i in unique:

    ### Extract full sample name
    regex = i + "-"
    regex2 = re.compile(fr'{regex}') # searched for 'samplename-' hyphen needed for end of samples or S2 and S22 would be mixed (e.g.).
    selected_files = list(filter(regex2.search, files2)) #searches how many files with same sample name

    outputname = selected_files[0].split("-")
    outputname2 = outputname[1] + "-familymerged"

    # count how many files to merge
    count_number = str(len(selected_files))

    if len(selected_files) < 2:
        raise SyntaxError('ERROR there is less than 2 samples in family')
    else:
        # Samples for intersect
        selected_files2 = ' '.join(selected_files)

        script2 = f"""bcftools isec -c indels -n {count_number} -o {outputname2} -p {outputname2} {selected_files2}"""
        print(script2)
        os.system(script2)

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
    sample = Alist[0]
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
    outputname2 = outputname[0] + "-" + outputname[1] + "-callermerged"

    # count how many files to merge
    sample_wild = i + '*.gz'
    cmd_count = 'ls -l ' + sample_wild + ' | wc -l'
    count = subprocess.run([cmd_count], stdout=subprocess.PIPE, shell = True)
    count_number = str(int(count.stdout))

    if int(count.stdout) < 2:
        raise SyntaxError('ERROR there is less than 2 "caller" samples')
    else:
        script2 = 'bcftools isec -c indels -n +' + count_number + ' -o ' + sample + '.caller.snps.merged.vcf -p ' + outputname2 + ' ' + sample_wild
        print(script2)
        os.system(script2)

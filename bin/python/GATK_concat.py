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
parser.add_argument('--vcf', required=True)
parser.add_argument('--vcf2', required=True)

args = parser.parse_args()
filesA = args.vcf
filesB = args.vcf2
files = filesA + filesB

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
    outputname2 = outputname[0] + "-" + outputname[1] + "-GATK.vcf"

    # count how many files to merge
    sample_wild = i + '*.gz'
    cmd_count = 'ls -l ' + sample_wild + ' | wc -l'
    count = subprocess.run([cmd_count], stdout=subprocess.PIPE, shell = True)
    count_number = str(int(count.stdout))

    if int(count.stdout) < 2:
        raise SyntaxError('ERROR there is less than 2 "vcf" samples')
    else:
        script2 = f'bcftools concat -a -O v -o {outputname2} ' + sample_wild
        logscript = fr"echo '{script2}' > {outputname2}_GATK_concat.log"
        print(script2)
        print(logscript)
        os.system(script2)
        os.system(logscript)

# This will count the number of files produced to ensure it is the same number as input.
cmd_count2 = 'ls -l *-GATK.vcf | wc -l'
count2 = subprocess.run([cmd_count2], stdout=subprocess.PIPE, shell = True)

cmd_count3 = 'ls -l *indels_filtered.vcf.gz | wc -l'
count3 = subprocess.run([cmd_count3], stdout=subprocess.PIPE, shell = True)

if int(count2.stdout) != int(count3.stdout):
    raise SyntaxError('Number of input and output files are different. There is a problem')


print("completed everything")

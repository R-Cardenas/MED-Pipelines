# RPC 031120
# This will merge fastq files from different lanes into one using linux 'cat'
# This was built for the nf-core RNA-seq pipeline.
# !/usr/bin/env python
import os
import argparse

#files = "18TB0035-PrimaryTumourMucosa-L002-R2.fastq.gz 18TB0035-PrimaryMucosa-L002-R2.fastq.gz 17TB0788-IsolatedCrypts-L002-R1.fastq.gz"
#############
# ARG PARSE #
#############

# Split the input from nextflow (space delim)
parser = argparse.ArgumentParser()
parser.add_argument('--fq', required=True)

args = parser.parse_args()
files = args.fq
files2 = files.split(" ")

###################################################################
# Will extract unique samples names (ie values 0 and 1 of fq file)
###################################################################

fq_samples = list()
for f in files2:
    filelist = f.split("-")
    samples = filelist[0] + "-" + filelist[1] + "-" + filelist[3] # remove lane from name
    fq_samples.append(samples) # append into new lisgt

fq_unique = set(fq_samples) # remove duplicates

###################################################################
# Performs concatination on samples using bash concat
###################################################################

for f in fq_unique:
    filelist = f.split("-")

    # create wild card using samples names and read number
    wildcard_fq = filelist[0] + "-" + filelist[1] + "*" + filelist[2]

    # create the file output name
    output_name = filelist[0] + "-" + filelist[1] + "-merged-" + filelist[2]

    # command for concatonating files in bash
    cmd = "cat " + wildcard_fq + " > " + output_name

    print(cmd)
    os.system(cmd)


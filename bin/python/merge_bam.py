# RPC 170920
# This will merge sample bam files from different lanes into one using samtools
# !/usr/bin/env python
import os
import argparse
import time

parser = argparse.ArgumentParser()
files = parser.add_argument('--bam', required=True)
#files = "S1208-S12.bam S1208-S13.bam S1205-12.bam"
files2 = files.split(" ")

bam_samples = list()
for f in files2:
    list = f.split("-")
    sample = list[0]
    bam_samples.append(sample)

unique = set(bam_samples)

# For loop submits each sample to be lane merged
for i in unique:
    list = i.split("-")
    sample = list[0]
    script = '"samtools merge ' + sample + '_merged.bam ' + sample + '*bam -@ 5"'
    script2 = '#SBATCH -p compute-16-24 --mem 10G --job-name=bam_merge_' + sample + ' -e %j.err ' + script
    print(script2)
    os.system(script2)
    time.sleep(180)

# Checks that the file exists (ie samtools output)
def loop():
    for i in unique:
        list = i.split("-")
        sample = list[0]
        file = sample + "_merged.bam"
        if os.path.isfile(file):
            print ("File exist")
        else:
            print('file does not exist yet')
            time.sleep(360)
            loop()

# Check that the file is not still increasing in size (ie being processed)
# This is performed as after job is submitted by SBATCH, nextflow will continue and move output files to next step whether or not the have been completed.
def size():
    for i in unique:
        list = i.split("-")
        sample = list[0]
        file = sample + "_merged.bam"
        size1 = os.stat(file).st_size
        time.sleep(180)
        size2 = os.stat(file).st_size
        if size1 == size2: # compares same files two minutes apart. If they are the same size it is presumed processing has finished.
            print("files are the same size")
            os.system("echo 'size1, size2 " + size1 + size2 + file + "' >> size.log")
        else:
            size()

loop()
size()
time.sleep(500)

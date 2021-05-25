# RPC 170920
# This will merge sample bam files from different lanes into one using samtools
# !/usr/bin/env python
import os
import argparse
import time
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--vcf', required=True)
parser.add_argument('--output', required=True)

args = parser.parse_args()
files = args.vcf
output = args.output

if ".gz" in files:
    script = f"""zcat {files} | vcf-annotate --fill-type | grep -oP 'TYPE=\w+' | sort | uniq -c > {output}"""
else:
    script = f"""cat {files} | vcf-annotate --fill-type | grep -oP 'TYPE=\w+' | sort | uniq -c > {output}"""



p = subprocess.run(script, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
print(p)


time.sleep(10)
if os.path.getsize(output) == 0:
    raise SyntaxError("ERROR: File is empty")
else:
    print('Good: File is not empty')

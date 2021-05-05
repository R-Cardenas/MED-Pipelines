# This will be run by the run_pipeline.py if the master_user_config.yaml
# select rna-seq. This script will select merge or nomerge pipeline
# and will select the correct genome
import os
# RPC 080720
# aft19qdu@uea.ac.uk
import glob
from bin.python.data_yaml import data
import subprocess
import re

#################
## MERGE LANES ##
#################

# Select merge or no merge nextflow
if data['merged_lanes'] == 'no':
    variant_nf = "nextflow run nfcore-rnaseq-nomerge.nf" 
    cmd = "cp RNAseq/expression/nfcore-rnaseq-nomerge.nf ."
    p = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
    print(p)

elif data['merged_lanes'] == 'yes':
    variant_nf = "nextflow run nfcore-rnaseq-merge.nf"
    cmd = "cp RNAseq/expression/nfcore-rnaseq-merge.nf ."
    p = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
    print(p)

else:
    print('dna_exome.py - line10')
    raise SyntaxError("dna_exome.py: Incorrect 'merged_lanes' input. Please revise")

########################
## Output Dir         ##
########################

output_dir = data['output_dir']
if output_dir == "":
    raise SyntaxError("YAML output_dir value is empty")
else:
    output_dir2 = "--outdir " + output_dir
    replace_string = f"sed -i 's:--outdir .*:{output_dir2}:g' nfcore-rnaseq*.nf" # replace rna-seq nf with YAML input
    p = subprocess.run(replace_string, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
    print(p)
    print('updated output_dir...')

########################
## Input Dir         ##
########################

input_dir = data['fastq_dir']
if input_dir == "":
    input_dir2 = 'params.fq = "$baseDir/input/*{fq,fastq}.gz"'
else:
    input_dir2 = "params.fq = " + input_dir + "/*{fq,fastq}.gz"
    input_dir3 = re.sub('//','/',input_dir2) # will replce double back slash to single (may be introduced by user in yaml)
    replace_string = f"sed -i 's:$params.fq =.*:{input_dir3}:g' nfcore-rnaseq*.nf" # replace rna-seq nf with YAML input
    p = subprocess.run(replace_string, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
    print(p)
    print('updated input_dir...')


###################
## Genome Config ##
###################

# Select config files for hg19 or hg38
if data['genome_assembly'] == 'hg38':
    rnaseq_genome = "--genome GRCh38 \\"
    replace_string = f"sed -i 's/--genome .*/{rnaseq_genome}/g' nfcore-rnaseq*.nf'"
    os.system(replace_string)
    print('Updated genome to hg38')

elif data['genome_assembly'] == 'hg19':
    rnaseq_genome = "--genome GRCh37 \\"
    replace_string = f"sed -i 's/--genome .*/{rnaseq_genome}/g' nfcore-rnaseq*.nf'"
    os.system(replace_string)
    print('Updated genome to hg19')
else:
    print('dna_exome.py - line24')
    raise SyntaxError("dna_exome.py: Incorrect 'genome_assembly' input. Please revise")

print('RNAseq module completed')

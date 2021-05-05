# DESCRIPTION
# RPC 080720
# aft19qdu@uea.ac.uk
import os
import sys
import subprocess
import glob
import time
import re

from bin.python.data_yaml import data
##################
## HOUSEKEEPING ##
##################


# Remove historical bash scripts
rm_existing = "rm -fr run_selected_pipeline.sh"

SBATCH_header = """#!/bin/bash
#SBATCH -p compute-16-64
#SBATCH --output=test_%j.log
#SBATCH -e test-%j.out
#SBATCH -e test-%j.err
#SBATCH -t 168:00:00
#SBATCH --mem=2G

module add nextflow
module add singularity
module add python/anaconda/2020.07
module add samtools
"""

new_sh = f"echo '{SBATCH_header}' >> run_selected_pipeline.sh"
p = subprocess.run(rm_existing, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
print(p)
p = subprocess.run(new_sh, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
print(p)


# Remove historical variables from config files
rm_projectname = f"""for f in $(find . -name '*config'); do sed -i '/env.projectname/d' $f; done"""
rm_build = f"""for f in $(find . -name '*config'); do sed -i '/env.build/d' $f; done"""
rm_bait = f"""for f in $(find . -name '*config'); do sed -i '/env.bait_interval/d' $f; done"""
rm_target = f"""for f in $(find . -name '*config'); do sed -i '/env.target_interval/d' $f; done"""
p = subprocess.run(rm_projectname, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
print(p)
p = subprocess.run(rm_build, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
print(p)
p = subprocess.run(rm_bait, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
print(p)
p = subprocess.run(rm_target, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
print(p)


## Remove the nextflow files that may have been copied to baseDir
cp_nf = 'rm -fr *.nf'
p = subprocess.run(cp_nf, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
print(p)

#####################
### PROJECT NAME ####
#####################


projectname = data['projectname'].lower() # extract the projects name
projectname2 = 'env.projectname = "' + projectname + '"'
print("Project name is: " + projectname2)
## Below line requires python >3.6
add_projectname = f"""for f in $(find . -name '*config'); do echo 'env.projectname={projectname2}' >> $f; done""" # add projectname to all config files
p = subprocess.run(add_projectname, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)



genome_assembly = data['genome_assembly'].lower() # extract the projects name
genome_assembly2 = 'env.projectname = "' + genome_assembly + '"'
print("Project name is: " + genome_assembly2)
## Below line requires python >3.6
add_build = f"""for f in $(find . -name '*config'); do echo 'env.build={genome_assembly2}' >> $f; done""" # add projectname to all config files
os.system(add_build)



#################################################
### EXOME VS WGS VS RNA-seq --- MAPPING ONLY ####
#################################################

# If BAMs only is yes. The cgpMAP module/ python scripts are not called.
if data['bams_only'].lower() == 'yes':
    print("Using bams only, ommitting cgpMAP module")
elif data['bams_only'].lower() == 'no':
    # The next chunk identifies data type and then imports the python script specific for the datatype e.g. exome
    if data['samples'].lower() == 'dna-exome':
        print('Loading DNA exome module \n')
        from DNAseq.Exome.cgpmap.dna_exome import *
    elif data['samples'].lower() == 'dna-wgs':
        print('Loading DNA WGS module \n')
        from DNAseq.WGS.dna_wgs import *
    elif data['samples'].lower() == 'rna-seq':
        print('Loading RNA-seq module \n')
        from RNAseq.expression.rnaseq import *
    else:
        raise SyntaxError('incorrect "samples" values input in master_user_config.yaml')
else:
    raise SyntaxError('Unrecognisable input for bams_only in master_user_config.yaml please revise')

###########################
### SOMATIC VS GERMLINE ###
###########################

## You need to put in WGS

if data['samples'].lower() == 'dna-exome' and data['variant'] == 'germline':
    variant_nf = "nextflow run freebayes_individual.nf -c DNAseq/Exome/cgpmap/cgpmap_hg38.config & nextflow run haplotypecaller_individual.nf -c DNAseq/Exome/cgpmap/cgpmap_hg38.config"
    cmd1 = "cp DNAseq/Exome/germline/freebayes/freebayes_individual.nf ."
    cmd2 = "cp DNAseq/Exome/germline/gatk/haplotypecaller_individual.nf ."
    p = subprocess.run(cmd1, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
    print(p)
    p = subprocess.run(cmd1, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
    print(p)

elif data['samples'].lower() == 'dna-exome' and data['variant'] == 'somatic':
    variant_nf = "nextflow run cgpwxs_v0.1.nf -c DNAseq/Exome/cgpmap/cgpmap_hg38.config & nextflow run mutect2_individual.nf -c DNAseq/Exome/cgpmap/cgpmap_hg38.config"
    cmd1 = "cp DNAseq/Exome/somatic/cgpwxs_v0.1.nf ."
    cmd2 = "cp DNAseq/Exome/somatic/mutect2_individual.nf ."
    p = subprocess.run(cmd1, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
    print(p)
    p = subprocess.run(cmd1, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
    print(p)
elif data['samples'].lower() == 'rna-seq':
    from RNAseq.expression.rnaseq import variant_nf
else:
    raise SyntaxError('incorrect "variant" values input in master_user_config.yaml')


### Concat nextflow run and config into bash script in repo home Dir
cgpmap_bash2 = "echo '" + variant_nf + "' >> run_selected_pipeline.sh"
os.system(cgpmap_bash2)


###################################
### SOMATIC ONLY - SAMPLE NAMES ###
##################################

if data['matched_identical'].lower() == 'yes':
    print('Tumor/Normal samples haved identical names. Extracting from bam files name...')
    from DNAseq.Exome.somatic.cgpwxs.tumor_normal_var import *
elif data['matched_identical'].lower() == 'no':
    print('Tumor/Normal samples do NOT identical names. Extracting from YAML file')
    from DNAseq.Exome.somatic.cgpwxs.samples_parse import *
else:
    pass # is not essential for germline / rna-seq

##########################
### Input & output Dir ### and out dir? - not for rna-seq
########################## is different

# Input fastqs
fastq_input = data['fastq_dir'] + "/*{1,2}.{fq,fastq}.gz"
fastq_input = re.sub('//','/',fastq_input)
remove_inputDir = f"""find . -name "*.nf" -exec sed -i '/params.fq = /d' {{}} \;"""

p = subprocess.run(remove_inputDir, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
print(p)


# Output dir (not for rna-seq - see rnaseq.py)
output_base = data['output_dir']
remove_outputDir = f"""find . -name "*.nf" -exec sed -i '/params.outputdir = /d' {{}} \;"""

p = subprocess.run(remove_outputDir, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
print(p)

# Add input and output dirs back in with YAML input:
input1 = f"""params.fq = \"{fastq_input}\""""
output1 =  f"""params.outputdir = \"{output_base}\""""

nf_files = glob.glob("*.nf")
for f in nf_files:
    os.system(f"""sed -i '1s|^|\\n|' {f}""")
    os.system(f"""sed -i '1s|^|{input1}|' {f}""")
    os.system(f"""sed -i '1s|^|\\n|' {f}""")
    os.system(f"""sed -i '1s|^|{output1}|' {f}""")

############
### RUN ###
###########


run_master = ["sbatch run_selected_pipeline.sh"]
p = subprocess.run(run_master, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
print("Job submitted - details below:")
print(p.stdout)

#for i in progressbar(range(100)):
time.sleep(0.07)

## Remove the nextflow files that were copied to keep clean
#cp_nf = 'rm -fr *.nf'
#ßos.system(cp_nf)
print("Pipeline has been sucessfully submitted")

exit()

# This will be run by the run_pipeline.py if the master_user_config.yaml
# This script will extract tumor normal samples names from the yaml file
# and add them to XXX file for nextflow to use in cgp pipelines.
# NOTE samples are to be given in order.

import os
# RPC 080720
# aft19qdu@uea.ac.uk
import glob

##################
## Housekeeping ##
##################

# delete tumor/normal environment variables in all config files in pipelines_clean
rm_tumor = f"""for f in $(find . -name '*config'); do sed -i '/params.normal/d' $f; done"""
rm_normal = f"""for f in $(find . -name '*config'); do sed -i '/params.tumor /d' $f; done"""
os.system(rm_tumor)
os.system(rm_normal)
print("Historical tumor/normal samples deleted from configs")

#############
## samples ##
#############

files = sorted(glob.glob('input/*.{bam,fq.gz,fastq.gz,cram}'))
raise SyntaxError("input Dir is empty") if len(files) == 0 else print('files identified in input folder')

tumor_final = list()
normal_final = list()

for f in files:
    Alist = f.replace('input/','') # remove path
    Alist = Alist.replace('.bam','') # remove file ext
    Alist = Alist.replace('.fq.gz','') # remove file ext
    Alist = Alist.replace('.fastq.gz','') # remove file ext
    Alist = Alist.replace('.cram','') # remove file ext
    Alist = Alist.split("-")

    sample = Alist[0].upper() # extract samples name
    family = Alist[1].upper() # extract tumor/sample label

    if family == "TUMOUR" or family == "TUMOR":

        full_name = sample + '-' + family
        tumor_final.append(full_name)

    elif family == "NORMAL":
        full_name = sample + '-' + family
        normal_final.append(full_name)
    else:
        print("ERROR: Cannot find TUMOUR or NORMAL labels on samples")
        print("Please visit site below for proper file naming strategy:")
        print("https://github.com/R-Cardenas/pipelines_clean")
        raise SyntaxError("")

###################
## Check Numbers ##
###################


if len(tumor_final) == len(normal_final):
    print("Tumour and Normal number of samples are the same. Good.")
else:
    raise SyntaxError("Tumour and Normal number of samples are NOT the same. Please check sample numbers and names.")


########################
## Check sample names ##
########################

for i,j in zip(tumor_final,normal_final):

    Alist = i.split("-")
    tumor = Alist[0].upper()

    Alist = j.split("-")
    normal = Alist[0].upper()

    if tumor == normal:
        text = "Sample normal and tumor have matched sample names - tumour: " + i + ' normal: ' + j
        print(text)
    else:
        text = "Sample normal and tumor have NOT matched sample names - tumour: " + i + ' normal: ' + j + '\n PLEASE RENAME or use YAML input.'
        raise SyntaxError(text)

########################
## Add to config file ##
########################


tumor_final_str = str(tumor_final).replace("'",'"').strip("'") # nextflow only accepts double quotations
normal_final_str = str(normal_final).replace("'",'"').strip("'")


tumor_cmd = f"""for f in $(find . -name '*config'); do echo 'params.tumor={tumor_final_str}' >> $f; done""" # add tumor sample names to all config files
normal_cmd = f"""for f in $(find . -name '*config'); do echo 'params.normal={normal_final_str}' >> $f; done"""

os.system(tumor_cmd)
os.system(normal_cmd)

print("Samples have successfully been added to the config files")

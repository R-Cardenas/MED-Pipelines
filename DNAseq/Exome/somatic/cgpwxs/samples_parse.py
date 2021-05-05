# This will be run by the run_pipeline.py if the master_user_config.yaml
# This script will extract tumor normal samples names from the yaml file
# and add them to XXX file for nextflow to use in cgp pipelines.
# NOTE samples are to be given in order.

import os
# RPC 080720
# aft19qdu@uea.ac.uk
from bin.python.data_yaml import data
import glob


##################
## Housekeeping ##
##################

#  delete tumor/normal environment variables
rm_tumor = f"""for f in $(find . -name '*config'); do sed -i '/params.normal/d' $f; done"""
rm_normal = f"""for f in $(find . -name '*config'); do sed -i '/params.tumor /d' $f; done"""
os.system(rm_tumor)
os.system(rm_normal)
print("Historical tumor/normal samples deleted from configs")

#################
## Extract     ##
#################
# Extract data from yaml

tumor = data['Tumour_samples']
normal = data['Normal_samples']


#Split by comma
tumor2 = tumor.replace(" ","").replace("'",'"').split(",")
normal2 = normal.replace(" ","").replace("'",'"').split(",")


###################
## Check Numbers ##
###################
# Check lengths are the same

if len(tumor2) == len(normal2):
    print("Tumor and Normal sample numbers are equal. Good.")
    print("WARNING - Ensure that sample matches are in the same order.")
else:
    raise SyntaxError("Tumour and Normal number of samples are NOT the same. Please check sample numbers and names.")


########################
## Add to config file ##
########################

tumor_final = list()
normal_final = list()

for i,j in zip(tumor2,normal2):

    tumor_final.append(i)
    normal_final.append(j)


tumor_final_str = str(tumor_final).replace(" ","").replace("'",'"')
normal_final_str = str(normal_final).replace(" ","").replace("'",'"')


tumor_cmd = f"""for f in $(find . -name '*config'); do echo 'params.tumor={tumor_final_str}' >> $f; done""" # add tumor sample names to all config files
normal_cmd = f"""for f in $(find . -name '*config'); do echo 'params.normal={normal_final_str}' >> $f; done"""

os.system(tumor_cmd)
os.system(normal_cmd)

print("Samples have successfully been added to the config files")

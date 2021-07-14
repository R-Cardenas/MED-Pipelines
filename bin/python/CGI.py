import os
import glob
import requests
import re
import time
from pathlib import Path
# The waits are necessary to allow the CGI API to prrocess the samples
# as the program allows a set number. In addition the write function
# in the python script will try to write a 'running' job giving an empty file.


files = glob.glob('*.vcf.gz') #contains whole path

print(files)
if len(files) == 0:
    syntaxError("No files detected")
else:
    print("files found")
    for f in files:
        print(f)

# ------------- Send files to cancer genomic interpreter ------------------------

for i in files:
    headers = {'Authorization': 'r.cardenas@uea.ac.uk de60bc15d5b92e56b8b2'}
    payload = {'cancer_type': 'COCA', 'title': 'COCA_mutect2' + str(i)} # i names the title using filename
    r = requests.post('https://www.cancergenomeinterpreter.org/api/v1',
                        headers=headers,
                        files={
                            'mutations': open(str(i), 'rb') # j path of files
                        },
                        data=payload,)
    time.sleep(5)
    print("VCF files submitted")
time.sleep(30)

# ------------------------------ obtain identifiers -----------------------------

headers = {'Authorization': 'r.cardenas@uea.ac.uk de60bc15d5b92e56b8b2'}
ID = requests.get('https://www.cancergenomeinterpreter.org/api/v1', headers=headers)
ID = ID.json()
print(ID)

# -------------------------------- obtain file names ---------------------------
for x in ID:
    info = requests.get('https://www.cancergenomeinterpreter.org/api/v1/' + x, headers=headers)
    info = info.json()
    title = (info['metadata'])['title']
    print(title)

    #Download the files and name correctly
    payload = {'action': 'download'}
    r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/' + x, headers=headers, params=payload)
    with open(title + 'CGI.zip', 'wb') as fd:
        fd.write(r._content)
    time.sleep(10)

    my_file = Path(title + 'CGI.zip')

    if my_file.is_file():
        print("file exists")
    else:
        syntaxError("file does not exist!!")

# ----------------------------------Delete the jobs after it has completed and been downloaded ------
for x in ID:
    headers = {'Authorization': 'r.cardenas@uea.ac.uk de60bc15d5b92e56b8b2'}
    r = requests.delete('https://www.cancergenomeinterpreter.org/api/v1/' + x, headers=headers)
    r.json()

print("CGI Completed")

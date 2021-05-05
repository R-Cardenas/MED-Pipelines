# cgpMAP pipeline

This pipeline utilises the sanger cgpMAP pipeline (https://github.com/cancerit/dockstore-cgpmap), in addition to some QC tools. This manual describes how this works and how to use it on its own without the master-YAML.

There are two versions of this pipeline, 'dna-exome-merge' and 'dna-exome-no-merge' - the merging factor pertains to samples run on different lanes that require to be merged once the BAMs have been formed. This results in an extra 'process' step in the 'merge' pipeline that will merge BAMs before the QC steps, but all other processes are identical in both pipelines. Below the different inputs/processes are described, only the most complicated



## Input

Below the params.read1 variable specifies the fastq directory which are read as pairs of files.
The files are piped into two channels where the fastq files are processed separately. The *$baseDir* is specifies the working directory that nextflow is run from - similar to the $PWD linux command line utility.


```go
params.read1 = "$baseDir/input/*{1,2}.fq.gz"

read1_ch = Channel .fromFilePairs( params.read1 )

read1_ch.into { read2_ch; read3_ch }

```

## Trim_galore

This is the first process. The files are copied into the work dir as usually the files are kept on different drives and sym/hardlinks tend to fail. The destination of the files is set by the *storeDir* variable.

```go
stageInMode = 'copy' // trim_galore doesnt like sym/hardlinks.
storeDir "$baseDir/output/cgpMAP/trim_galore"
```

The script in this process uses trim-galore and then renames the files replacing '_val' extension with '.trim.'. This is done as nextflow allows the automatic removal of file extensions using the '.simpleName' or '.baseName' function, this helps prevent the accumulation of unnecessary file extensions. The copied files (ie not the trimmed files) are deleted to preserve space.

```bash
trim_galore --paired --fastqc --illumina \
--basename ${reads[0].simpleName} \
${reads[0]} ${reads[1]}

# rename val with trim
mv ${reads[0].simpleName}_val_1.fq.gz ${reads[0].simpleName}_1.trim.fq.gz
mv ${reads[0].simpleName}_val_2.fq.gz ${reads[0].simpleName}_2.trim.fq.gz

# delete copied files
rm -fr ${reads[0]} # remove the copied files to prevent memory loss
rm -fr ${reads[1]}
```

## FastQ
This process uses the fastq headers to extract platform_unit (PU), library preparation (LB) and lane information to form a YAML file to inform cgpMAP to annontate the BAM correctly.

First extract the header using fqtools - but only using a header containing no 'N' in the sequencing primer field. This was done to ensure the sequencing primer for read1 and read2 are indentical.

The second part is an if statement - if the txt file produced is empty then repeat the extraction, but excluding the regex grep. Due to differences in fastq header this may occur, however for the majority this will work.

(For more information about fastq headers see: https://angus.readthedocs.io/en/2017/Read_group_info.html)

```bash
# Read1
# Extract header
fqtools -d header ${read1} | grep ":[C,A,T,G]*[+][C,A,T,G]" | head -1 > ${read1.simpleName}.txt

### Counts lines in file1 and will repeat if it is empty
words=`wc -l ${read1.simpleName}.txt  | awk '{print \$1}'`;
if [ \$words -eq 0 ]
then
fqtools -d header ${read1} | head -1 > ${read1.simpleName}.txt
else
echo 'alls good!'
fi
```

Next a python script is used to form a YAML file in the correct format for cgpMAP (described here: https://github.com/cancerit/PCAP-core/wiki/File-Formats-groupinfo.yaml)

This script assumes the platform unit is Illumina, and will add the sample name obtained from the fastq filename (E.g. S0102-P-TUMOR-L01-R1.fastq.gz sample name is S0102-P).

If the txt file is still empty due to a malformed header, python will insert generic inputs into the YAML. As this may affect the QC of BAMs, a warning is produced by the python script and is moved to the logs folder.

```bash
python $baseDir/nextflow_pipelines/bin/python/fastq2config_cgpmap.py \
--fq1 ${read1.simpleName}.txt --fq2 ${read2.simpleName}.txt \
--n1 ${read1} --n2 ${read2} --o ${read1}.yaml

cp *WARNING* $baseDir/logs 2>/dev/null || :
```

An example of the YAML file produced is shown below:

```yaml
SM: S0102-P
READGRPS:
  1.fq:
    PL: ILLUMINA
    LB: HGYYCBBXX_NTTCCT
    PU: HGYYCBBXX.6
  2.fq:
    PL: ILLUMINA
    LB: HGYYCBBXX_NTTCCT
    PU: HGYYCBBXX.6
```

## cgpMAP

The cgpMAP process receives inputs from two processes. The trimmed fastq files (read1,read2) from trim_galore and the YAML config file from fastQ.

```go
storeDir "$baseDir/output/cgpMAP/${read1.simpleName}"
input:
val read1 from read5_ch
val read2 from read10_ch
val yaml from yaml_ch.collect()

```

The following step will run cgpMAP. The name is first extracted from the trimmed fastq filename ($name).

Then the cgpMAP tool is called (ds-cgpmap.pl) and inserts the corresponding inputs. The following variables can be found in the /DNAseq/Exome/cgpmap/cgpmap_hg{19,38}.config:

  *cgpmap_genome
  *gpmap_index

For more information how cgpMAP works visit: https://github.com/cancerit/dockstore-cgpmap.

CgpMAP creates a folder with all the output files into a seperate folder. However, nextflow is expecting the files in a particular folder so it can be passed onto the next process. Therefore following the ds-cgpmap.pl the BAM files are moved one level below.

To ensure the correct fastq pairs have been used, the read1 and read2 variable are recorded into a text file and logged in ${projectname}_cgpmap_samples.log.


```sh
name=\$(echo '${read2}' | sed -e 's/.*[/]//' -e 's/_.*//')

ds-cgpmap.pl  \
-outdir $baseDir/output/cgpMAP/${read1.simpleName} \
-r $cgpmap_genome \
-i $cgpmap_index \
-s \$name \
-t 5 \
-g ${read1}.yaml \
${read1} ${read2}

mv $baseDir/output/cgpMAP/${read1.simpleName}/*.bam \
$baseDir/output/cgpMAP/${read1.simpleName}/${read1.simpleName}.bam

echo 'fq1: ${read1} fq2: ${read2} bam_name: ${read1.simpleName}' >> $baseDir/${projectname}_cgpmap_samples.log
ls -l  ${read1} >> $baseDir/logs/symbolic_test_fastq.log
ls -l  ${read2} >> $baseDir/logs/symbolic_test_fastq.log

```

# Merge BAMs

This is relevant for the dna-exome-merge.nf pipeline. Nextflow will wait until all of the files from the previous process has completed and pool them for this process.
The script below describes collecting these files.

```sh
input:
file bam from bam_merge_ch.collect()
```

The next script extracts the sample name from the bam files (using sed), sort them, and only use unique sample names (uniq) in the for loop.

Using the sample name and wildcard samtools is used to merge BAMs from different lanes but containing the same sample name.

```sh
# This will extract the sample name and create a samtools wild card
for f in \$(ls *.bam | sed -e 's/.*[/]//' -e 's/_.*//' | sort | uniq)
do
samtools merge \${f}_merged.bam \${f}*

# documents the samples used for merging into log
printf "`echo \${f}` `echo $(ls \${f}*)`\\n" >> $baseDir/logs/bam_merge_log.txt

done

```

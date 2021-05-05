# UEA BCRE pipelines - Mutect2

This README explains how the mutect2 pipeline works.

## Input

Firstly, as we are using normal and tumor matched samples - they need to be annotated accordingly. The top of the script contains tumor and normal parameters such as below:

```go
// sample channel (VALUES)

params.tumor = ["PD13382a","PD13399a"]
params.normal = ["PD13382b","PD13399b"]

// Define Tumor vs Normal Variables
tumor_ch = Channel .from (params.tumor )
normal_ch = Channel .from (params.normal )


// bam channel (files)

*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "$baseDir/{output/aligned_sorted/*.rename.bam,input/*.bam}"
bam_ch = Channel .fromPath( params.bam )

bam_ch.into { bam2_ch; bam3_ch }
tumor_ch = Channel .from (params.tumor )
normal_ch = Channel .from (params.normal )
```

NOTE: It is important that the samples are in the correct order, ie params.tumor[1] (PD13382a) and params.normal[1] (PD13382b), as this is the order which nextflow will pipe these into the pipeline. These values are used only in the 'mutect2 process'.

## General proccesses

All of the processes except for mutect2 are straight forward. They accept one file and output one or more files to the next process. Below is an example of the first process which exemplifies this:

```go
process BaseRecalibrator {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_haplotypeCaller"
  input:
  file bam from bam2_ch
  output:
  file "${bam.simpleName}.BQSR.bam" into (mutect2_1_ch,pileup_ch)
	file "${bam}.table"
  script:
  """
	mkdir -p tmp
  gatk BaseRecalibrator \
	-I ${bam} \
  -R $genome_fasta \
  --known-sites $GATK_dbsnp138 \
  --known-sites $GATK_1000G \
  --known-sites $GATK_mills \
  -O ${bam}.table \
	--tmp-dir tmp

	gatk ApplyBQSR \
  -R $genome_fasta \
  -I ${bam} \
  --bqsr-recal-file ${bam}.table \
  -O ${bam.simpleName}.BQSR.bam \
	--tmp-dir tmp
	rm -fr tmp

  """
}
```

The BaseRecalibrator proces will receive one file at a time from the bam2_ch and runs the 'script' on it (e.g. mkdir -p tmp ...) and once the script has finished nextflow will check if the output is present, and if so store it in the 'storeDir'. The files have symlink file created and passed into their own channel, in this case its two channels (mutect2_1_ch and pileup_ch). These channels pass the files onto the next process.

## Mutect2 proccess

The script for mutect2 is more advanced.

```go
process mutect2 {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/mutect2"
  input:
  val x from tumor_ch
  val y from normal_ch
  file bam from mutect2_1_ch.collect()
  output:
  file "*.vcf.gz" into filter_vcf_ch
  script:
  """
	python $baseDir/bin/python/run_mutect2.py \
	--tumor '${x}' \
	--normal '${y}' \
	--germline '$Mutect2_germline' \
	--PoN '$Mutect2_PoN' \
	--ref '$genome_fasta'
    """
}
```

We have two values, one from params.tumor and the other from params.normal - variable x and y. However, we are still inheriting the bam files from the previous process (`file bam from mutect2_1_ch.collect()`), in fact we are inheriting all of the bam files - since the collect function is used. Since there is no collect function on the tumor_ch and normal_ch, this process will repeat as it loops through the values in these channel (e.g. see param.normal). Since nextflow uses symlinks, collecting the bam files each time wont take much space.

The values (x,y) are passed onto a python script (run_mutect2.py) where it uses the values to search for the files and run the mutect2 script. This method (using values and collect etc) ensures that the normal and tumor files used are correct. Other methods within nextflow using symlinks resulted in renaming of files silently with no errors. There are probably more elegant solution using nextflow but this was a robust method I knew how.

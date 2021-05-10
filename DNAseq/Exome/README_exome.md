# UEA BCRE pipelines - Exome

<br />

<!-- TABLE OF CONTENTS -->
## Table of Contents

<br />

- [Input parameters](#Input-parameters)
- [Selecting genome](#Selecting-genome)
- [HPC parameters](#HPC-parameters)
- [Singularity paramters](#Singularity-paramters)
- [Usage](#Usage)
- [Susequent pipelines](#Susequent-pipelines)

## Input parameters

In order to run the exome QC pipeline, the nextflow.config file needs to first be edited. The following variables are important:

```
env.projectname = "PROJECTNAME"
env.params_fq = "PATH/TO/FASTQ.GZ/FILE/{1,2}.fq.gz"
env.bait_interval  = "PATH/TO/BAIT/bait.bed"
env.target_interval = "PATH/TO/TARGET/target.bed"
```

## Selecting genome

You also have to choose which genome you would like to use. The paths are already written out, you need to comment out ("//" in nextflow) the genomes you do NOT want to use. Below shows the two options - with GRCh38 chosen to be used.

```
//includeConfig "DNAseq/Exome/cgpmap/cgpmap_hg19.config"
includeConfig "DNAseq/Exome/cgpmap/cgpmap_hg38.config"
```

These config files contain all the necessary PATHs for other required config files, such as Freebayes/gatk_haplotypecaller/mutect2 configs, and thus providing the specific files for future pipelines for specific genomes.

## HPC parameters

We have included slurm parameters that run efficiently on the UEA HPC. These mainly describe how many jobs per process can be submitted in total and at what rate. These are conservative but are likely to prevent overloading the HPC. These paramters are for ALL pipelines.

```
executor {
  $slurm {
    queueSize = 8
    submitRateLimit = '1 / 45sec'
  }
  $local{
    queueSize = 100
    errorStrategy = 'retry'
    maxRetries = 5
    submitRateLimit = '1 / 10sec'
  }
}
```

## Singularity paramters

This enables singularity to be used with nextflow. It sets the home directory to the directory which the nextflow pipeline was started. It also binds the directory where we keep all our genomes (outside the singularity) to a place within the singularity image (/var/spool/mail). Thus when specifiying the genome files in our pipeline they will look like this '/var/spool/mail/genome.fa'.

```
singularity {
  enabled = true
  runOptions = "--home $baseDir --bind /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes/:/var/spool/mail"
}
```
## Usage

After the config files have been updated. We next have to copy the require files to our working directory and run it - as shown below:

```
module add nextflow
module add singularity

cp DNAseq/Exome/cgpmap/dna-exome-merge.nf .

nextflow run dna-exome-merge.nf

```

This will start the pipeline!

As you may have noticed there are two dna-exome pipelines:

- dna-exome-merge.nf
- dna-exome-nomerge.nf

The merge describes if the samples have been sequenced on multiple lanes. If they have (E.g. sample1-family1-L1.fastq.gz,sample1-family1-L2.fastq.gz,sample1-family1-L3.fastq.gz) the bam files after alignment will be merged. The software will recognise that the sample ID (i.e. sample1 from our example) are all the same and that the lane information (i.e. L1,L2,L3) are different. Therefore, it is imperative that the filename structure be kept in the format described in the [home repo](../../README.md).


## Subsequent pipelines

Following the completion of the pipeline, it is good to have a look at the multiQC report (output/multiqc/multiqc_report.html) which will summarise all of the QC analyses performed, in addition to any logs produced (logs/...). Following this you will run the next pipeline depending on the analysis you wish to perform. Below are the links for the subsequent pipelines:

Germline-pipelines

- Freebayes
- GATK HaplotypeCaller
- Variant processing (family analysis)
- Variant processing (individual analysis)

Somatic-pipelines

- cgpWXS
- GATK Mutect2

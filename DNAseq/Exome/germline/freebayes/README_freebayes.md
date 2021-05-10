# UEA BCRE pipelines - Exome Freebayes

<br />

<!-- TABLE OF CONTENTS -->
## Table of Contents

<br />

- [Usage](#Usage)
- [Singularity paramters](#Singularity-paramters)


## Usage

After running the cgpMAP all files are stored in the output directory and this pipeline uses these files and requires no input. Therefore, we next copy the freebayes nextflow script to the working directory and run to start.

```
module add nextflow
module add singularity

cp /DNAseq/Exome/germline/freebayes/freebayes.nf .

nextflow run freebayes.nf

```

This will run the pipeline.

## Subsequent pipelines

Once completed some of the next pipelines below may be run next.

- Variant processing (family analysis) (req. completion of GATK pipeline also)
- Variant processing (individual analysis) (req. completion of GATK pipeline also)

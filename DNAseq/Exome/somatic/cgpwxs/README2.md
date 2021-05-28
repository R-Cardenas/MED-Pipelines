# UEA BCRE pipelines - Exome cgpWXS

<br />

<!-- TABLE OF CONTENTS -->
## Table of Contents

<br />

- [Usage](#Usage)
- [Singularity paramters](#Singularity-paramters)


## Usage

After running the cgpMAP all files are stored in the output directory and this pipeline uses these files and requires no input. Therefore, we next copy the cgpWXS nextflow script to the working directory and run to start.

*(REMEMBER: ensure you have listed your samples names (e.g. fastq filenames) for matched Tumor and Normal samples. These need to be in the same order for N and T lists)*

```
module add nextflow
module add singularity

cp /DNAseq/Exome/somatic/cgpwxs/cgpwxs.nf .

nextflow run cgpwxs.nf

```

This will run the pipeline.

## Subsequent pipelines

Once completed some of the next pipelines below may be run next.

- Needs updating......

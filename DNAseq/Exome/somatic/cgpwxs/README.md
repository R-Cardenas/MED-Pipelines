# UEA BCRE pipelines - cgpWXS


This README explains how the cgpWXS pipeline works. Although the pipeline is quite short as the cgp-series is a pipeline its self, certain parameters are required to ensure it runs.

## Input

Firstly, as we are using normal and tumor matched samples - they need to be annotated accordingly. The top of the script contains tumor and normal parameters such as below:

```go
params.tumor = ["PD13382a","PD13399a"]
params.normal = ["PD13382b","PD13399b"]

// Define Tumor vs Normal Variables
tumor_ch = Channel .from (params.tumor )
normal_ch = Channel .from (params.normal )
```

NOTE: we are not importing the file locations. But the values of the sample names and this will be used later for wildcards, therfore it is imperative the files are in the correct place (i.e. $baseDir/input or .... where? RMD)

It is important that the samples are in the correct order, ie params.tumor[1] (PD13382a) and params.normal[1] (PD13382b), as this is the order which nextflow will pipe these into the pipeline.

## Process - upper

We have seperated the process into upper and lower sections. The script for the upper is as shown below:

```go
process cgpwxs2
	storeDir "$baseDir/output/cgpwxs/${tumor}_vs_${normal}"
	input:
	val tumor from tumor_ch
	val normal from normal_ch
	output:
	file "WXS_${tumor}_vs_${normal}.result.tar.gz"

```

The storeDir (store directory) includes the file names to give each pair of samples their own directory. This is because cgpWXS creates a number of generic filenames that end up overwriting eachother which causes the pipeline to fail.

Within the input field, we have specified them as val (values) instead of files which is what we normally do. With file, nextflow will create a symlink to the file to its own 'work' directory and uses the symlink to for its proccesses. Val on the other hand will not create any symlinks and will use the full path of the file.

We use val as we are not using files but the sample names (e.g. params.tumor). If you look at the lower process script we use wildcards to find the correct files (--tumour or --tidx).

```go

script:
"""
mkdir -p $baseDir/output/cgpwxs/${tumor}_vs_${normal}

singularity exec --cleanenv \
--home $baseDir \
--bind /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes:/var/spool/ref:ro \
--bind $baseDir/input:/var/spool/data:ro \
/gpfs/afm/cg_pipelines/Pipelines/singularity/images/cgpwxs_3.1.6.img \
ds-cgpwxs.pl \
-reference $cgp_ref \
-annot $cgp_annot \
-snv_indel $cgp_snv_indel \
-tumour /var/spool/data/${tumor}*.bam \
-tidx /var/spool/data/${tumor}*.bam.bai \
-normal /var/spool/data/${normal}*.bam \
-nidx /var/spool/data/${normal}*.bam.bai \
-exclude NC_007605,hs37d5,GL% \
-outdir $baseDir/output/cgpwxs/${tumor}_vs_${normal} \
-sp "Homo sapiens" \
-assembly "GRCh37"
"""

```

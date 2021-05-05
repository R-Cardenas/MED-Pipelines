UEA BCRE pipelines - cgpWGS


This README explains how the cgpWGS pipeline works. Although the pipeline is quite short as the cgp-series is a pipeline its self, certain parameters are required to ensure it runs.

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
process cgpwxs2 {
	executor 'slurm'
	clusterOptions '--qos=hmem'
	queue 'hmem-512'
	memory '400 GB'
	storeDir "$baseDir/output/cgpwgs22"
	input:
	val tumor from tumor_ch
	val normal from normal_ch
	output:
	file "${tumor}_vs_${normal}/WXS_${tumor}_vs_${normal}.result.tar.gz"
	script:
	"""
	python $baseDir/cgpwgs_wrapper.py --tumor '${tumor}' --normal '${normal}' --assembly 'hg19' --home '$baseDir' \
	--reference '$cgp_ref' --annotation '$cgp_annot' --snv_indel '$cgp_snv_indel' \
	--cgp_QC '$cgp_QC' --cgp_CNV '$cgp_CNV'
	"""
}

```

The storeDir (store directory) includes the file names to give each pair of samples their own directory. This is because cgpWGS creates a number of generic filenames that end up overwriting eachother which causes the pipeline to fail.

Within the input field, we have specified them as val (values) instead of files which is what we normally do. With file, nextflow will create a symlink to the file to its own 'work' directory and uses the symlink to for its proccesses. Val on the other hand will not create any symlinks and will use the full path of the file.

We use val as we are not using files but the sample names (e.g. params.tumor). If you look at the lower process script we use wildcards to find the correct files (--tumour or --tidx).

This script part of the pipeline will call upon a python3 script which accepts input values and runs bash script below with the inputs. This python script requires glob to find the files using the sample name (e.g. params.normal) and then run the cgpWGS pipeline. For some reason using a bash script and wildcards wouldn't work with this pipeline/container.

```python
cmd2 = fr"""singularity exec --cleanenv \
--home {home1} \
--bind /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes:/var/spool/ref:ro \
--bind {home1}/input:/var/spool/data:ro \
/gpfs/afm/cg_pipelines/Pipelines/singularity/images/cgpwgs_2.1.0.img \
 ds-cgpwgs.pl \
-r {ref1} \
-a {ann1}  \
-si {snv_indel1} \
-cs  {cgp_CNV1} \
-qc {cgp_QC1} \
-pl 3.65 -pu 1.0 \
-e 'MT,GL%,hs37d5,NC_007605' \
-tumour {tumor_bam} \
-tidx {tumor_bai} \
-normal {normal_bam} \
-nidx {normal_bai} \
-outdir {home1}/output/cgpwgs22/{tumor1}_vs_{normal1} \
-sp "Homo sapiens" \
-assembly "GRCh37"
"""
```

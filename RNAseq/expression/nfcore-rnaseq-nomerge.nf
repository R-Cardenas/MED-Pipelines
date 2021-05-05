/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.outputdir = "/gpfs/afm/cg_pipelines/Pipelines/Williams_RNASeq_processing"
params.fq = "/gpfs/afm/cg_pipelines/Pipelines/Troeberg_sulf2/bsg-ftp.well.ox.ac.uk/201214_A00711_0315_AHMFJJDSXY/fastqs2process/*{fq,fastq}.gz"
fq_ch = Channel .fromPath( params.fq )

println """\
	\
	\
	\
         ==================================
         NF - core: RNAseq
				 v0.1
         ===================================



         """
         .stripIndent()

// needs to be tested
process main_nf{
	//errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	//maxRetries 6
	executor 'local'
//	stageInMode = 'copy'
	input:
	file fastq from fq_ch.collect()
	script:
	"""
	module add nextflow
	module add singularity

	nextflow run nf-core/rnaseq -r 1.4.1 -profile singularity \
	--reads './*{1,2}.{fq,fastq}.gz' \
	-c $baseDir/RNAseq/expression/UEA2.config \
	--genome GRCh37 \
	--outdir '$baseDir/output' \
	--max_memory '49.GB' \
	--saveAlignedIntermediates \
	--saveTrimmed \
	--email_on_fail 'aft19qdu@uea.ac.uk' \
	--max_cpus 2 2>&1 | tee > $baseDir/pipeline.${task.attempt}.output
	"""
}

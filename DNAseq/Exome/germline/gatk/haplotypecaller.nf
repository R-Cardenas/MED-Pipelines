
/*
 * create a channel for bam files produced by Pipeline GATK germline (single)_processing pipeline
 */
params.bam = "$baseDir{/output/BAM/merge/RMD/*.rmd.bam,/input/*.bam}"
bam_ch = Channel .fromPath( params.bam )


println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         GATK Germline - HaplotypeCaller
         Single samples
         v0.2
         ===================================



         """
         .stripIndent()

process BaseRecalibrator {
	//errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	//maxRetries 6
	executor 'slurm'
	memory { 7.GB * task.attempt }
  storeDir "$baseDir/output/GATK_haplotypeCaller"
  input:
  file bam from bam_ch
  output:
  file "${bam.simpleName}.BQSR.bam" into haplotype_bam_ch
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
	--create-output-bam-index true \
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


process haplotypeCaller {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  executor 'slurm'
	clusterOptions '--qos=hmem'
	queue 'hmem-512'
	memory { 30.GB * task.attempt }
  storeDir "$baseDir/output/GATK_haplotypeCaller"
  input:
	file bam from haplotype_bam_ch
  output:
  file "${bam.simpleName}.g.vcf.gz" into haplotype2_ch
	file "${bam}.bai"
  script:
  """
	mkdir -p tmp
	gatk BuildBamIndex \
	-I ${bam} \
	-O ${bam}.bai \
	--TMP_DIR tmp

  gatk HaplotypeCaller \
  -R $genome_fasta \
  -I ${bam} \
	--read-index ${bam}.bai \
  -O ${bam.simpleName}.g.vcf.gz \
  --create-output-variant-index true \
	--intervals ${target_interval} \
  -ERC NONE \
	--tmp-dir tmp
	rm -fr tmp
  """
}


process CNNscoreVariants {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	executor 'slurm'
	memory { 7.GB * task.attempt }
  storeDir "$baseDir/output/GATK_haplotypeCaller"
  input:
  file vcf from haplotype2_ch
  output:
  file "${vcf.simpleName}.CNN.g.vcf.gz" into filterVCF_ch
  script:
  """
	mkdir -p tmp
	gatk IndexFeatureFile \
		-F ${vcf} \
		--tmp-dir tmp

  gatk CNNScoreVariants \
  -V ${vcf} \
  -R $genome_fasta \
	-O ${vcf.simpleName}.CNN.g.vcf.gz \
	--tmp-dir tmp
	"""
}

process FilterVariantTranches {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	executor 'slurm'
	memory { 7.GB * task.attempt }
  storeDir "$baseDir/output/GATK_haplotypeCaller"
	input:
	file vcf from filterVCF_ch
	output:
	file "${vcf.simpleName}-GATK.vcf" into zip_ch
	script:
	"""
	mkdir -p tmp
	gatk IndexFeatureFile \
		-F ${vcf} \
		--tmp-dir tmp

	gatk FilterVariantTranches \
	-V ${vcf} \
	--info-key CNN_1D \
	--resource $GATK_dbsnp138 \
  --resource $GATK_1000G \
  --resource $GATK_mills \
	--resource $GATK_hapmap \
	--snp-tranche 99.95 \
	--indel-tranche 99.4 \
	-O ${vcf.simpleName}-GATK.vcf \
	--tmp-dir tmp
	"""
}

// needs samttools
process zip {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	executor 'local'
  storeDir "$baseDir/output/VCF_collect"
	input:
	file zip from zip_ch
	output:
	file "${zip}.gz" into merge2_ch
	file "${zip}.gz.csi" into merge3_ch
	script:
	"""
	bgzip ${zip}
	bcftools index -c ${zip}.gz


	echo ' Pipeline GATK germline (single) v0.2 completed
	 Project: $projectname
	 Time: ${nextflow.timestamp}
	 BaseRecalibrator - completed
	 (Reference: $genome_fasta
		known sites: $GATK_dbsnp138
		known sites: $GATK_1000G
		known sites: $GATK_mills
		target_bed: ${target_interval})

	 applyBaseRecalibrator - completed
	 bam index  - completed
	 haplotypeCaller - completed

	' >> $baseDir/logs/${projectname}_log.txt

	#mail -s "GATK germline (single) successful" aft19qdu@uea.ac.uk < $baseDir/logs/${projectname}_log.txt
	"""
}


workflow.onError {
	process finish_error{
		script:
		"""
		echo 'Pipeline GATK germline (single) v0.2 FAILED
		Project: $projectname
		Time: ${nextflow.timestamp}

		Error:
		${workflow.errorMessage}' >> $baseDir/logs/${projectname}_GATKerror.txt

	  mail -s "Pipeline GATK germline (single) unsuccessful - $projectname" aft19qdu@uea.ac.uk < $baseDir/logs/${projectname}_GATKerror.txt
	  """
	}
}

/// This works completely 180620

/*
 * create a channel for fastq pairs
 */
params.bam = "$baseDir{/output/BAM/merge/RMD/*.rmd.bam,/input/*.bam}"

Channel
	.fromPath( params.bam )
	.set {bam_ch}


println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         Freebayes - Germline
         v0.2
         ===================================



         """
         .stripIndent()


process Freebayes {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 7
  executor 'slurm'
	clusterOptions '--qos=hmem'
	queue 'hmem-512'
	memory { 30.GB * task.attempt }
  storeDir "$baseDir/output/freebayes"
  input:
  file bam from bam_ch
  output:
  file "${bam.simpleName}.raw.vcf" into (vcf_ch,count1_ch) // inset $project name from config file
  script:
  """
  freebayes ${bam} -f $genome_fasta > ${bam.simpleName}.raw.vcf
  """
}

process vcf_filter{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	memory { 7.GB * task.attempt }
  storeDir "$baseDir/output/freebayes"
  input:
  file vcf from vcf_ch
  output:
  file "${vcf.simpleName}-freebayes.vcf" into (zip2_ch, count2_ch)
  script:
  """
  bcftools filter -i 'QUAL>5 & INFO/DP>5 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1' ${vcf} > ${vcf.simpleName}-freebayes.vcf
  """
}

process count_variants{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	memory { 7.GB * task.attempt }
  storeDir "$baseDir/output/variant_counts/freebayes"
  input:
	file vcf from count1_ch
  file vcf2 from count2_ch
  output:
  file "*.count"
  script:
  """
	python $baseDir/bin/python/vcf_count.py --vcf $vcf --output ${vcf}.unfiltered.count
	python $baseDir/bin/python/vcf_count.py --vcf $vcf2 --output ${vcf2}.filtered.count
  """
}


process zip {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect"
	input:
	file zip from zip2_ch
	output:
	file "${zip}.gz" into merge_ch
	file "${zip}.gz.csi" into csi_ch
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
		""" echo 'Pipeline Freebayes FAILED
		Project: $projectname
		Time: ${nextflow.timestamp}

		Error:
		${workflow.errorMessage}' {input} >> $baseDir/logs/${projectname}_error.txt

	  #mail -s "cgpMAP successful" {input} < $baseDir/logs/${projectname}_error.txt
	  """
	}
}

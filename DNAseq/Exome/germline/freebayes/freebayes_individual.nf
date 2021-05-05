/// This works completely 180620

/*
 * create a channel for fastq pairs
 */
params.bam = "$baseDir{/output/BAM/merged_lanes/*.rmd.bam,/input/*.bam}"

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
	maxRetries 2
  storeDir "$baseDir/output/freebayes"
  input:
  file bam from bam_ch
  output:
  file "${bam.simpleName}.raw.vcf" into vcf_ch // inset $project name from config file
  script:
  """
  freebayes ${bam} -f $genome_fasta > ${bam.simpleName}.raw.vcf
  """
}

process vcf_filter {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/freebayes"
  input:
  file vcf from vcf_ch
  output:
  file "${vcf.simpleName}.filtered.vcf" into zip2_ch
  script:
  """
	bcftools filter -i 'QUAL>5 & INFO/DP>5 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1' \
  ${vcf} > ${vcf.simpleName}.filtered.vcf
  """
}


process zip {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect"
	input:
	file zip from zip2_ch
	output:
	file "${zip}.freebayes.vcf.gz" into merge_ch
	file "${zip}.freebayes.vcf.gz.tbi" into csi_ch
	script:
	"""
	bgzip ${zip}
	mv ${zip}.gz ${zip}.freebayes.vcf.gz
	bcftools index -t ${zip}.freebayes.vcf.gz
	"""
}


workflow.onComplete {
	// create log files and record output
	process finish{
		script:
		""" echo ' Pipeline Freebayes germline v0.2 completed
		Project: $projectname
		Time: ${nextflow.timestamp}

	 Freebayes - completed
	 (Reference - $genome_fasta)
	 VCFTools - completed
	 (Parameters: "bcftools filter -i 'QUAL>5 & INFO/DP>5 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1' )
	 Bgzip - completed
	 BcfTools Index - completed
	 ' >> $baseDir/${projectname}_log.txt

	 mail -s "Freebayes successful" < $baseDir/${projectname}_log.txt
	 """
	}
}

workflow.onError {
	process finish_error{
		script:
		""" echo 'Pipeline Freebayes FAILED
		Project: $projectname
		Time: ${nextflow.timestamp}

		Error:
		${workflow.errorMessage}' {input} >> $baseDir/${projectname}_error.txt

	  mail -s "cgpMAP successful" {input} < $baseDir/${projectname}_error.txt
	  """
	}
}

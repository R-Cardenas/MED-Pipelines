params.tumor = ["PD13382a.cram","PD13389a.cram","PD13399a.cram"]
params.normal = ["PD13382b.cram","PD13389b.cram","PD13399b.cram"]

/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "$baseDir/{output/aligned_sorted/*.rename.bam,input/*.bam}"
bam_ch = Channel .fromPath( params.bam )

bam_ch.into { bam2_ch; bam3_ch }
tumor_ch = Channel .from (params.tumor )
normal_ch = Channel .from (params.normal )

println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         GATK Somatic: Mutect2
				 Single samples
				 v0.1
         ===================================



         """
         .stripIndent()


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

process mutect2 {
	//errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	//maxRetries 6
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

process pileup_summary{
  storeDir "$baseDir/output/mutect2/mutect2"
	input:
	file bam from pileup_ch.flatten()
	output:
	file "${bam.simpleName}.getpileupsummaries.table" into contamination_ch
	script:
	"""
	gatk GetPileupSummaries \
	-I ${bam} \
	-V $Mutect2_germline \
	-L $Mutect2_germline \
	-O ${bam.simpleName}.getpileupsummaries.table
	"""
}

process calculate_contamination{
	storeDir "$baseDir/output/mutect2/mutect2"
	input:
	file table from contamination_ch
	output:
	file "${table.simpleName}_calculatecontamination.table" into filter_vcf2_ch
	script:
	"""
	gatk CalculateContamination \
	-I ${table} \
	-O ${table.simpleName}_calculatecontamination.table
	"""
}


process filter_vcf {
	storeDir "$baseDir/output/mutect2/mutect2/filtered_vcf"
	input:
	file vcf from filter_vcf_ch
	file table from filter_vcf2_ch
	output:
	file "${vcf.simpleName}.filtered.vcf" into zip_ch
	script:
	"""
	gatk FilterMutectCalls \
	-R $genome_fasta \
	-V ${vcf}  \
	--contamination-table ${table} \
	-O ${vcf.simpleName}.filtered.vcf
	"""
}

// needs samttools
process zip {
  storeDir "$baseDir/output/VCF_collect"
	input:
	file zip from zip_ch
	output:
	file "${zip}.gz" into merge_ch
	script:
	"""
	bgzip ${zip}
	"""
}


workflow.onComplete {
	// create log files and record output
	process finish{
		script:
		""" echo ' Pipeline GATK germline (cohort) v0.3 completed
		Project: $projectname
		Time: ${nextflow.timestamp}
		BaseRecalibrator - completed
		(Reference: $genome_fasta
		 known sites: $GATK_dbsnp138
		 known sites: $GATK_1000G
		 known sites: $GATK_mills)

		applyBaseRecalibrator - completed
		bam index  - completed
		Mutect2- completed
		(germline-resource: $Mutect2_germline \
		 panel-of-normal: $Mutect2_PoN )

		Pile up summary - completed
		Contamination calculation - completed
		Filter Mutect VCF - completed


	 ' >> $baseDir/${projectname}_log.txt

	 mail -s "GATK germline (cohort) successful" aft19qdu@uea.ac.uk < $baseDir/${projectname}_log.txt
	 """
	}
}

workflow.onError {
	process finish_error{
		script:
		"""
		echo 'Pipeline cgpmap v0.3 FAILED
		Project: $projectname
		Time: ${nextflow.timestamp}

		Error:
		${workflow.errorMessage}' >> $baseDir/${projectname}_error.txt

	  mail -s "cgpMAP successful" aft19qdu@uea.ac.uk < $baseDir/${projectname}_error.txt
	  """
	}
}

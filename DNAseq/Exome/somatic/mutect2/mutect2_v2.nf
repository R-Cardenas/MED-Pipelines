/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "$baseDir{/output/BAM/merge/RMD/*.rmd.bam,/input/*.bam}"
bam_ch = Channel .fromPath( params.bam )

bam_ch.into { bam2_ch; bam3_ch }
tumor_ch = Channel .from (params_tumor )
normal_ch = Channel .from (params_normal )

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
	maxRetries 7
  executor 'slurm'
	clusterOptions '--qos=hmem'
	queue 'hmem-512'
	memory { 30.GB * task.attempt }
  storeDir "$baseDir/output/mutect2"
  input:
  file bam from bam2_ch
  output:
  file "${bam.simpleName}.BQSR.bam" into (mutect2_1_ch,pileup_ch)
	file "*.bai" into bai_ch
	file "${bam}.table"
  script:
  """
	# cgpmap only add ID (RGID) and Sample (RGSM)
	# we need to add some values. Null should be the same as no value
	gatk AddOrReplaceReadGroups \
	-I ${bam} \
	-O ${bam.simpleName}.replace.head.bam \
	-PU ILLUMINA \
	-SM ${bam.simpleName}.BQSR.bam \
	-LB Exome \
	-PL ILLUMINA

	mkdir -p tmp

	gatk BuildBamIndex \
	-I ${bam} \

  gatk BaseRecalibrator \
	-I ${bam.simpleName}.replace.head.bam \
  -R $genome_fasta \
  --known-sites $GATK_dbsnp138 \
  --known-sites $GATK_1000G \
  --known-sites $GATK_mills \
  -O ${bam}.table \
	--tmp-dir tmp

	gatk ApplyBQSR \
  -R $genome_fasta \
  -I ${bam.simpleName}.replace.head.bam \
  --bqsr-recal-file ${bam}.table \
  -O ${bam.simpleName}.BQSR.bam \
	--tmp-dir tmp
	rm -fr tmp ${bam.simpleName}.replace.head.bam

  """
}

process mutect2 {
	// errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	// maxRetries 3
  executor 'slurm'
	clusterOptions '--qos=hmem'
	queue 'hmem-512'
	memory 85.GB
  storeDir "$baseDir/output/mutect2"
  input:
  val x from tumor_ch
  val y from normal_ch
  file bam from mutect2_1_ch.collect()
  output:
  file "${x}_v_${y}.unfiltered.GATK.vcf.gz" into count1_ch
	file "${x}_v_${y}.filtered.GATK.vcf.gz" into (count2_ch,zip_ch)
  script:
  """
	python3 $baseDir/bin/python/run_mutect2.py \
	--tumor '${x}' \
	--normal '${y}' \
	--intervals ${target_interval} \
	--germline '$Mutect2_germline' \
	--PoN '$Mutect2_PoN' \
	--ref '$genome_fasta'

	gatk FilterMutectCalls \
	-R $genome_fasta \
	-V ${x}_v_${y}.unfiltered.GATK.vcf.gz \
	-O ${x}_v_${y}.filtered.GATK.vcf.gz

    """
}

// process count_variants{
// 	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
// 	maxRetries 6
// 	executor 'local'
//   storeDir "$baseDir/output/variant_counts/GATK"
//   input:
// 	file vcf from count1_ch
//   file vcf2 from count2_ch
//   output:
//   file "*.count"
//   script:
//   """
// 	python $baseDir/bin/python/vcf_count.py --vcf $vcf --output ${vcf}.unfiltered.count
// 	python $baseDir/bin/python/vcf_count.py --vcf $vcf2 --output ${vcf2}.filtered.count
//   """
// }
//
// // needs samttools
// process zip {
// 	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
// 	maxRetries 6
// 	memory { 7.GB * task.attempt }
//   storeDir "$baseDir/output/VCF_collect"
// 	input:
// 	file zip from zip_ch
// 	output:
// 	file "${zip}.gz" into merge_ch
// 	script:
// 	"""
// 	bgzip ${zip}
//
// 	echo ' Pipeline GATK germline (cohort) v0.3 completed
// 	Project: $projectname
// 	Time: ${nextflow.timestamp}
// 	BaseRecalibrator - completed
// 	(Reference: $genome_fasta
// 	 known sites: $GATK_dbsnp138
// 	 known sites: $GATK_1000G
// 	 known sites: $GATK_mills)
//
// 	applyBaseRecalibrator - completed
// 	bam index  - completed
// 	Mutect2- completed
// 	(germline-resource: $Mutect2_germline \
// 	 panel-of-normal: $Mutect2_PoN )
//
// 	Pile up summary - completed
// 	Contamination calculation - completed
// 	Filter Mutect VCF - completed
//
//
//  ' >> $baseDir/logs/${projectname}_log.txt
//
//  #mail -s "GATK germline (cohort) successful" aft19qdu@uea.ac.uk < $baseDir/${projectname}_log.txt
// 	"""
// }
//
// workflow.onError {
// 	process finish_error{
// 		script:
// 		"""
// 		echo 'Pipeline cgpmap v0.3 FAILED
// 		Project: $projectname
// 		Time: ${nextflow.timestamp}
//
// 		Error:
// 		${workflow.errorMessage}' >> $baseDir/${projectname}_error.txt
//
// 	  #mail -s "cgpMAP successful" aft19qdu@uea.ac.uk < $baseDir/${projectname}_error.txt
// 	  """
// 	}
// }

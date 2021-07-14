
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
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 3
  executor 'slurm'
	clusterOptions '--qos=hmem'
	queue 'hmem-512'
	memory { 150.GB * task.attempt }
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
	maxRetries 3
  executor 'slurm'
	clusterOptions '--qos=hmem'
	queue 'hmem-512'
	memory { 150.GB * task.attempt }
  storeDir "$baseDir/output/GATK_haplotypeCaller"
  input:
	file bam from haplotype_bam_ch
  output:
  file "${bam.simpleName}.g.vcf.gz" into (haplotype2_ch,count1_ch)
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
	maxRetries 3
  executor 'slurm'
	clusterOptions '--qos=hmem'
	queue 'hmem-512'
	memory { 150.GB * task.attempt }
  storeDir "$baseDir/output/GATK_haplotypeCaller"
  input:
  file vcf from haplotype2_ch
  output:
  file "${vcf.simpleName}.indels_filtered.vcf.gz" into merge_ch
	file "${vcf.simpleName}.snps_filtered.vcf.gz" into merge2_ch
	file "${vcf.simpleName}*filtered.vcf.gz.tbi" into merge22_ch
  script:
  """


	gatk SortVcf -I ${vcf} -O sorted.${vcf}
	gatk IndexFeatureFile -F sorted.${vcf}

	gatk SelectVariants \
    -V sorted.${vcf} \
    -select-type SNP \
    -O snps.vcf.gz

	gatk SelectVariants \
    -V sorted.${vcf} \
    -select-type INDEL \
    -O indels.vcf.gz

	gatk VariantFiltration \
    -V snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${vcf.simpleName}.snps_filteredd.vcf.gz

		gatk SelectVariants \
		-V ${vcf.simpleName}.snps_filteredd.vcf.gz \
		--exclude-filtered true \
		-O ${vcf.simpleName}.snps_filtered.vcf.gz

	gatk VariantFiltration \
    -V indels.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O ${vcf.simpleName}.indels_filteredd.vcf.gz

	gatk SelectVariants \
		-V ${vcf.simpleName}.indels_filteredd.vcf.gz \
		--exclude-filtered true \
		-O ${vcf.simpleName}.indels_filtered.vcf.gz

	"""
}

process merge_caller_indels {
	//errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	//maxRetries 3
  executor 'slurm'
	memory { 7.GB * task.attempt }
  storeDir "$baseDir/output/GATK_haplotypeCaller"
  input:
  file vcf1 from merge_ch.collect()
	file vcf2 from merge2_ch.collect()
	file tbi from merge22_ch.collect()
  output:
  file "*-GATK.vcf" into (zip_ch,count2_ch)
	file "*concat.log"
  script:
  """
	python $baseDir/bin/python/GATK_concat.py --vcf "${vcf1}" --vcf2 "${vcf2}"
	"""
}



process count_variants{
	//errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	//maxRetries 6
	executor 'local'
  storeDir "$baseDir/output/variant_counts/GATK"
  input:
	file vcf from count1_ch.flatten()
  file vcf2 from count2_ch.flatten()
  output:
  file "{${vcf}.unfiltered.count,${vcf2}.filtered.tsv}"
  script:
  """
	python $baseDir/bin/python/vcf_count.py --vcf $vcf --output ${vcf}.unfiltered.count
	python $baseDir/bin/python/vcf_count.py --vcf $vcf2 --output ${vcf2}.filtered.tsv
  """
}


// needs samttools
process zip {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	executor 'local'
  storeDir "$baseDir/output/VCF_collect"
	input:
	file zip from zip_ch.flatten()
	output:
	file "${zip}.gz"
	file "${zip}.gz.csi"
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

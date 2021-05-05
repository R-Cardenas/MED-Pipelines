// This works on STAR aligned and picard remove duplicates bam files

params.tumor = ["17TB0788-Tumoroids-merged-R1Aligned.sortedByCoord.out.markDups","17TB0815-Tumoroids-merged-R1Aligned.sortedByCoord.out.markDups","18TB0008-Tumouroids-merged-R1Aligned.sortedByCoord.out.markDups","18TB0035-Tumouroids-merged-R1Aligned.sortedByCoord.out.markDups","18TB0071-Tumouroids-merged-R1Aligned.sortedByCoord.out.markDups","18TB0132-Tumouroids-merged-R1Aligned.sortedByCoord.out.markDups"]
params.normal = ["17TB0788-Organoids-merged-R1Aligned.sortedByCoord.out.markDups","17TB0815-Organoids-merged-R1Aligned.sortedByCoord.out.markDups","18TB0008-Organoids-merged-R1Aligned.sortedByCoord.out.markDups","18TB0035-Organoids-merged-R1Aligned.sortedByCoord.out.markDups","18TB0071-Organoids-merged-R1Aligned.sortedByCoord.out.markDups","18TB0132-Organoids-merged-R1Aligned.sortedByCoord.out.markDups"]

/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
//params.bam = "$baseDir/{output/aligned_sorted/*.rename.bam,input/*.bam}"
params.bam ="/gpfs/afm/cg_pipelines/Pipelines/Williams_RNAseq/output/markDuplicates/*bam"
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
         GATK RNAseq: Mutect2
				 Single samples
				 v0.1
         ===================================



         """
         .stripIndent()


process SplitNCigarReads {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	stageInMode 'copy'
  storeDir "$baseDir/output/SplitNCigarReads"
  input:
  file bam from bam2_ch
  output:
  file "${bam.simpleName}_CIGAR.bam" into bsqr_ch

  script:
  """
	chmod 777 ${bam}
gatk SplitNCigarReads \
      -R $genome_fasta \
      -I ${bam} \
      -O ${bam.simpleName}_CIGAR.bam
rm -fr ${bam}
  """
}

process BaseRecalibrator {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_haplotypeCaller"
  input:
  file bam from bsqr_ch
  output:
  file "${bam.simpleName}.BQSR.bam" into mutect2_1_ch
	file "${bam}.table" into table_ch
  script:
  """
	gatk AddOrReplaceReadGroups \
	-I ${bam} \
	-O ${bam.simpleName}.reheader.bam \
	-LB lib1 \
	-PL ILLUMINA \
	-PU instrument1 \
	-SM ${bam.simpleName}

	mkdir -p tmp
  gatk BaseRecalibrator \
	-I ${bam.simpleName}.reheader.bam \
  -R $genome_fasta \
  --known-sites $GATK_dbsnp138 \
  --known-sites $GATK_1000G \
  --known-sites $GATK_mills \
  -O ${bam}.table \
	--tmp-dir tmp

	gatk ApplyBQSR \
  -R $genome_fasta \
  -I ${bam.simpleName}.reheader.bam \
  --bqsr-recal-file ${bam}.table \
  -O ${bam.simpleName}.BQSR.bam \
	--tmp-dir tmp
	rm -fr tmp

	rm -fr ${bam.simpleName}.reheader.bam
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
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/mutect2/mutect2"
	input:
	file table from table_ch
	output:
	file "${table.simpleName}_covariates.pdf"
	script:
	"""
  gatk AnalyzeCovariates \
  -bqsr ${table} \
  -plots ${table.simpleName}_covariates.pdf
	"""
}

process filter_vcf {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	storeDir "$baseDir/output/mutect2/mutect2/filtered_vcf"
	input:
	file vcf from filter_vcf_ch.flatten()
	output:
	file "${vcf.simpleName}.filtered.vcf" into zip_ch
	script:
	"""
	gatk FilterMutectCalls \
	-R $genome_fasta \
	-V ${vcf}  \
	-O ${vcf.simpleName}.filtered.vcf
	"""
}

// needs samttools
process zip {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
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

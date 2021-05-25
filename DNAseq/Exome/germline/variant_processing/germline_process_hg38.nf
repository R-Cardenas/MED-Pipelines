params.vcf= "$baseDir/output/VCF_collect/*.vcf.gz"

vcf_ch = Channel. fromPath (params.vcf)
vcf_ch.into { vcf1_ch; vcf2_ch }


 // use pipeline bundle 1 has python3 installed
process merge_caller_indels {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/merge_vcf/indels/caller"
  input:
  file vcf from vcf2_ch.collect()
  output:
  file "*caller.merged.indels.vcf.gz" into (fam_ch,count1_ch)
  script:
  """
  for file in *.vcf*; do
    bcftools index \$file
  done

  python $baseDir/bin/python/merge_caller_indel.py --bam '$vcf'

  for file in *.vcf; do
    bgzip \$file
  done

  """
}

process count_variants{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	executor 'local'
  storeDir "$baseDir/output/variant_counts/indels/caller_merge"
  input:
	file vcf from count1_ch.flatten()
  output:
  file "${vcf}.CallerMerge.indels.count" into sumCounts1_ch
  script:
  """
  python $baseDir/bin/python/vcf_count.py --vcf $vcf --output ${vcf}.CallerMerge.indels.count1

  grep -v 'snp' ${vcf}.CallerMerge.indels.count1 > ${vcf}.CallerMerge.indels.count
  """
}

process merge_family_indels {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/merge_vcf/indels/family"
  input:
  file vcf from fam_ch.collect()
  output:
  file "*family.merged.indels.vcf.gz" into indels_filter_ch
  script:
  """

  for file in *.vcf*; do
    bcftools index \$file
  done

  python $baseDir/bin/python/merge_family_indel.py --bam '$vcf'

  for file in *.vcf; do
    bgzip \$file
  done
  """
}

process indels_filter {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  memory { 7.GB * task.attempt }
  executor 'slurm'
  storeDir "$baseDir/output/VCF_collect/merge_vcf/indels/family/filter"
  input:
  file vcf from indels_filter_ch.flatten()
  output:
  file "${vcf.simpleName}.indels.vcf.gz" into (indels_sort_ch,count2_ch)
  script:
  """
  bcftools view -O z -v indels ${vcf} > ${vcf.simpleName}.indels.vcf.gz
  """
}

process count_variants2{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	executor 'local'
  storeDir "$baseDir/output/variant_counts/indels/family_merge"
  input:
	file vcf from count2_ch.flatten()
  output:
  file "${vcf}.FamMerge.indels.count" into sumCounts2_ch
  script:
  """
	python $baseDir/bin/python/vcf_count.py --vcf $vcf --output ${vcf}.FamMerge.indels.count1

  grep -v 'snp' ${vcf}.FamMerge.indels.count1 > ${vcf}.FamMerge.indels.count
  """
}

process indels_sort {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/merge_vcf/indels/family/filter"
  input:
  file vcf from indels_sort_ch
  output:
  file "${vcf.simpleName}.indels.family.merged.sorted.vcf.gz" into vep_ch
  script:
  """
  mkdir -p tmp
  bcftools sort -O z -o ${vcf.simpleName}.indels.family.merged.sorted.vcf.gz ${vcf} -T tmp
  rm -fr tmp
  """
}

process VEP {
  executor 'slurm'
  memory { 7.GB * task.attempt }
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/VEP"
  input:
  file vcf from vep_ch.flatten()
  output:
  file "${vcf.baseName}_VEP.vcf"
  script:
  """
  vep -i ${vcf} \
  --dir_cache /var/spool/mail/VEP_hg38/.vep \
  -o ${vcf.baseName}_VEP.vcf \
  --cache \
  --offline \
  --species homo_sapiens \
  --assembly GRCh38 \
  -af_1kg \
  --fork 5 \
  --offline \
  --fasta /var/spool/mail/cgpwgs_ref/GRCh38/core_ref_GRCh38_hla_decoy_ebv/genome.fa \
  --custom /var/spool/mail/conservation/hg38.phastCons7way.bw,Conservation,bigwig,exact \
  --show_ref_allele \
  --use_given_ref \
  --verbose \
  --check_ref \
  --dont_skip \
  --af_gnomad \
  --vcf \
  --pubmed \
  --check_existing \
  --everything \
  --force_overwrite
  """
}


// and VEP_fasta added
// use pipeline bundle 1 has python3 installed
process merge_caller_snps {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/merge_vcf/snps/caller"
  input:
  file vcf from vcf1_ch.collect()
  output:
  file "*callermerged.vcf.gz" into (snps_fam_ch,count3_ch)
  script:
  """
  for file in *.vcf*; do
    bcftools index \$file
  done

  python $baseDir/bin/python/merge_caller_snps.py --bam '$vcf'

  mcp '*/0000.vcf' '#1.vcf'

  for file in *.vcf; do
    bgzip \$file
  done

  # The mcp function renames the file based on the path (where the * is )
  # this will take 0000.vcf which uses GATK vcf info fields.
  """
}

process count_variants3{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	executor 'local'
  storeDir "$baseDir/output/variant_counts/snps/caller_merge"
  input:
	file vcf from count3_ch.flatten()
  output:
  file "${vcf}.CallerMerge.snps.count" into sumCounts3_ch
  script:
  """
	python $baseDir/bin/python/vcf_count.py --vcf $vcf --output ${vcf}.CallerMerge.snps.count1

  grep 'snp' ${vcf}.CallerMerge.snps.count1 > ${vcf}.CallerMerge.snps.count
  """
}


process merge_fam_snps{
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/merge_vcf/snps/family"
  input:
  file vcf from snps_fam_ch.collect()
  output:
  file "*familymerged.vcf.gz" into snps_filter_ch
  script:
  """
  for file in *.vcf.gz; do
    bcftools index \$file
  done

  python $baseDir/bin/python/merge_family_snps.py --bam '$vcf'

  mcp '*/0000.vcf' '#1.vcf'

  for file in *.vcf; do
    bgzip \$file
  done

  # The mcp function renames the file based on the path (where the * is )
  # this will take 0000.vcf which uses GATK vcf info fields.
  """
}


// snps includes snps and MNPs
// rationale for mnps in this section and not indels
// mnps can be covered but whole reads length so variant callers shouldnt struggle as much as with indels
process snps_filter {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/merge_vcf/snps/family/filtered"
  input:
  file vcf from snps_filter_ch.flatten()
  output:
  file "${vcf.simpleName}.snps.vcf.gz" into (snps_sort_ch,count4_ch)
  script:
  """
  bcftools view -O z -V indels ${vcf} > ${vcf.simpleName}.snps.vcf.gz
  """
}

process count_variants4{
	executor 'local'
	memory { 7.GB * task.attempt }
  storeDir "$baseDir/output/variant_counts/snps/family_merge"
  input:
	file vcf from count4_ch.flatten()
  output:
  file "${vcf}.FamMerge.snps.count" into sumCounts4_ch
  script:
  """
	python $baseDir/bin/python/vcf_count.py --vcf $vcf --output ${vcf}.FamMerge.snps.count1

  grep 'snp' ${vcf}.FamMerge.snps.count1 > ${vcf}.FamMerge.snps.count
  """
}


process snps_sort {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/merge_vcf/snps/family/filtered"
  input:
  file vcf from snps_sort_ch
  output:
  file "${vcf.simpleName}.snps.family.merged.sorted.vcf.gz" into vep2_ch
  script:
  """
  mkdir -p tmp
  bcftools sort -O z -o ${vcf.simpleName}.snps.family.merged.sorted.vcf.gz ${vcf} -T tmp
  rm -fr tmp
  """
}

process VEP2 {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  executor 'slurm'
	clusterOptions '--qos=hmem'
	queue 'hmem-512'
	memory { 30.GB * task.attempt }
  storeDir "$baseDir/output/VCF_collect/VEP"
  input:
  file vcf from vep2_ch.flatten()
  output:
  file "${vcf.baseName}_VEP.vcf"
  script:
  """
  vep -i ${vcf} \
  --dir_cache /var/spool/mail/VEP_hg38/.vep \
  -o ${vcf.baseName}_VEP.vcf \
  --cache \
  --offline \
  --species homo_sapiens \
  --assembly GRCh38 \
  -af_1kg \
  --fork 15 \
  --offline \
  --fasta /var/spool/mail/cgpwgs_ref/GRCh38/core_ref_GRCh38_hla_decoy_ebv/genome.fa \
  --custom /var/spool/mail/conservation/hg38.phastCons7way.bw,Conservation,bigwig,exact \
  --show_ref_allele \
  --use_given_ref \
  --verbose \
  --check_ref \
  --dont_skip \
  --af_gnomad \
  --vcf \
  --pubmed \
  --check_existing \
  --everything \
  --force_overwrite
  """
}


process Rtables {
	executor 'local'
  storeDir "$baseDir/output/variant_counts/data"
  input:
  file indelC from sumCounts1_ch.collect()
  file indelF from sumCounts2_ch.collect()
  file snpC from sumCounts3_ch.collect()
  file snpF from sumCounts4_ch.collect()
  output:
  file "CallerOverlapCount.csv" into r1_ch
  file "FamilyOverlapCount.csv" into r2_ch

  script:
  """
  Rscript $baseDir/bin/R/variantCounts.R -v "$indelC" -f NA -w "*.CallerMerge.indels.count" -c NA -o Indels.CallerOverlap.csv -O "Caller Overlap"
  Rscript $baseDir/bin/R/variantCounts.R -v "$indelF" -f NA -w "*.FamMerge.indels.count" -c NA -o Indels.FamilyOverlap.csv -O "Group/Family Overlap"

  Rscript $baseDir/bin/R/variantCounts.R -v "$snpC" -f NA -w "*.CallerMerge.snps.count" -c NA -o Snps.CallerOverlap.csv -O "Caller Overlap"
  Rscript $baseDir/bin/R/variantCounts.R -v "$snpF" -f NA -w "*.FamMerge.snps.count" -c NA -o Snps.FamilyOverlap.csv -O "Group/Family Overlap"

  # Merge the SNPs and Indels
  tail -n +2 Indels.CallerOverlap.csv > Indels.CallerOverlap.csv
  cat Snps.CallerOverlap.csv Indels.CallerOverlap.csv > CallerOverlapCount.csv

  tail -n +2 Indels.FamilyOverlap.csv > Indels.FamilyOverlap.csv
  cat Snps.FamilyOverlap.csv Indels.FamilyOverlap.csv > FamilyOverlapCount.csv
  """
}

params.vcf= "$baseDir/output/VCF_collect/*.vcf.gz"

vcf_ch = Channel. fromPath (params.vcf)
vcf_ch.into { vcf1_ch; vcf2_ch; vcf3_ch }

//This is a fake PED of brentP
process ped {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	executor 'local'
  storeDir "$baseDir/output/VCF_collect/merge_vcf/indels/caller"
  input:
  file vcf from vcf3_ch.collect()
  output:
  file "project.PED" into (ped_ch,ped2_ch)
  script:
  """
  module add python/anaconda/2020.07
  python $baseDir/bin/python/PED_file.py --bam '$vcf'
  """
}

 // use pipeline bundle 1 has python3 installed
process merge_caller_indels {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	executor 'local'
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
	maxRetries 6
	executor 'local'
  storeDir "$baseDir/output/variant_counts/indels/caller_merge"
  input:
	file vcf from count1_ch.flatten()
  output:
  file "${vcf}.CallerMerge.indels.count" into sumCounts1_ch
  script:
  """
  bcftools view -O z -v indels ${vcf} > ${vcf.simpleName}.BCFTOOLS.TEMP.vcf.gz #removes non-indels for the count

  python $baseDir/bin/python/vcf_count.py --vcf ${vcf.simpleName}.BCFTOOLS.TEMP.vcf.gz --output ${vcf}.CallerMerge.indels.count1

  grep -v 'snp' ${vcf}.CallerMerge.indels.count1 > ${vcf}.CallerMerge.indels.count
  """
}

process merge_family_indels {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	executor 'local'
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
	executor 'local'
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
	maxRetries 6
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
	executor 'local'
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
  storeDir "$baseDir/output/VCF_collect/VEP/vcf"
  input:
  file vcf from vep_ch.flatten()
  output:
  file "${vcf.baseName}_VEP.vcf" into vepIndel_ch
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

process slivar2 {
  //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	//maxRetries 6
  storeDir "$baseDir/output/VCF_collect/VEP/slivar"
  input:
  file vcf from vepIndel_ch
  file ped from ped2_ch
  output:
  file "${vcf.simpleName}.indels.slivar.vcf" into vepIndel2_ch
  script:
  """
  slivar expr --js /var/spool/mail/slivar/slivar-functions.js \
  -g /var/spool/mail/slivar/gnomad.hg38.v2.zip \
  -g /var/spool/mail/slivar/topmed.hg38.dbsnp.151.zip \
  --info 'INFO.impactful && INFO.gnomad_popmax_af < 0.1 && variant.FILTER == "PASS" && variant.ALT[0] != "*" && INFO.topmed_af < 0.1'\
  --ped $ped \
  --pass-only \
  --vcf $vcf \
  --out-vcf ${vcf.simpleName}.indels.slivar.vcf
  """
}


process vcf2tsv {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/VEP/tsv/indels"
  input:
  file vcf from vepIndel2_ch
  output:
  file "${vcf.simpleName}.tsv" into tsvF_ch
  script:
  """
  Rscript $baseDir/bin/R/vcf2tsv.R -v $vcf -o ${vcf.simpleName}.tsv
  """
}


process tsvFilter {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/VEP/tsv/indels"
  input:
  file vcf from tsvF_ch.collect()
  output:
  file "ALL_data_indels.csv" into md2_ch
  script:
  """
  Rscript $baseDir/bin/R/filterVariantIndels.R -a 0.05
  """
}


// and VEP_fasta added
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
	maxRetries 6
	executor 'local'
	memory { 7.GB * task.attempt }
  storeDir "$baseDir/output/variant_counts/snps/caller_merge"
  input:
	file vcf from count3_ch.flatten()
  output:
  file "${vcf}.CallerMerge.snps.count" into sumCounts3_ch
  script:
  """
  bcftools view -O z -V indels ${vcf} > ${vcf.simpleName}.BCFTOOLS.TMP.snps.vcf.gz

	python $baseDir/bin/python/vcf_count.py --vcf ${vcf.simpleName}.BCFTOOLS.TMP.snps.vcf.gz --output ${vcf}.CallerMerge.snps.count1

  grep 'snp' ${vcf}.CallerMerge.snps.count1 > ${vcf}.CallerMerge.snps.count
  """
}


process merge_fam_snps{
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	executor 'local'
	memory { 7.GB * task.attempt }
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
	executor 'local'
	memory { 7.GB * task.attempt }
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
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
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
  storeDir "$baseDir/output/VCF_collect/VEP/vcf"
  input:
  file vcf from vep2_ch.flatten()
  output:
  file "${vcf.baseName}_VEP.vcf" into vepSnp_ch
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

process slivar {
  //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	//maxRetries 6
  storeDir "$baseDir/output/VCF_collect/VEP/slivar"
  input:
  file vcf from vepSnp_ch
  file ped from ped_ch
  output:
  file "${vcf.simpleName}.snps.slivar.vcf" into vepSnp2_ch
  script:
  """
  slivar expr --js /var/spool/mail/slivar/slivar-functions.js \
  -g /var/spool/mail/slivar/gnomad.hg38.v2.zip \
  -g /var/spool/mail/slivar/topmed.hg38.dbsnp.151.zip \
  --info 'INFO.impactful && INFO.gnomad_popmax_af < 0.1 && variant.FILTER == "PASS" && variant.ALT[0] != "*" && INFO.topmed_af < 0.1'\
  --ped $ped \
  --pass-only \
  --vcf $vcf \
  --out-vcf ${vcf.simpleName}.snps.slivar.vcf
  """
}


process vcf2tsv2 {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/VEP/tsv/snps"
  input:
  file vcf from vepSnp2_ch
  output:
  file "${vcf.simpleName}.tsv" into tsvF_ch2
  script:
  """
  Rscript $baseDir/bin/R/vcf2tsv.R -v $vcf -o ${vcf.simpleName}.tsv
  """
}

process tsvFilter2 {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/VEP/tsv/snps"
  input:
  file vcf from tsvF_ch2.collect()
  output:
  file "ALL_data_snvs.csv" into md1_ch
  script:
  """
  Rscript $baseDir/bin/R/filterVariantSNPs.R -a 0.05
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
  tail -n +2 Indels.CallerOverlap.csv > Indels.CallerOverlap2.csv
  cat Snps.CallerOverlap.csv Indels.CallerOverlap2.csv > CallerOverlapCount.csv

  tail -n +2 Indels.FamilyOverlap.csv > Indels.FamilyOverlap2.csv
  cat Snps.FamilyOverlap.csv Indels.FamilyOverlap2.csv > FamilyOverlapCount.csv
  """
}


process markdown {
//  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
//	maxRetries 6
  input:
  file vcf from md1_ch.collect()
  file vcf2 from md2_ch.collect()
  val count1 from r1_ch
  val count2 from r2_ch
  output:
  val "$baseDir/output/markdown/${projectname}_report.html"
  script:
  """
  mkdir -p $baseDir/output/markdown

  Rscript -e 'rmarkdown::render("$baseDir/bin/R/germlineExomeTemplate.Rmd",
  output_file="$baseDir/output/markdown/${projectname}_report.html",
  params = list(genome = "$cgpmap_genome",
  index = "$cgpmap_index",
  bait = "$bait_interval",
  target = "$target_interval",
  GATK_db = "$GATK_dbsnp138",
  GATK_1000G = "$GATK_1000G",
  GATK_mills = "$GATK_mills",
  FB_filter = "$baseDir/output/variant_counts/freebayes/",
  GATK_filter = "$baseDir/output/variant_counts/GATK/",
  CallerOverlap = "$count1",
  FamilyOverlap = "$count2"))'
  """
}

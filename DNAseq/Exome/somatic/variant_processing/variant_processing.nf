params.vcf= "$baseDir/output/VCF_collect/*.vcf.gz"

params.vcf= "$baseDir/output/VCF_collect/*.vcf.gz"

vcf_ch = Channel. fromPath (params.vcf)
vcf_ch.into { vcf1_ch; vcf2_ch }


 // use pipeline bundle 1 has python3 installed
process merge_caller_indels {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	executor 'local'
  storeDir "$baseDir/output/VCF_collect/merge_vcf/indels/caller"
  input:
  file vcf from vcf2_ch.collect()
  output:
  file "*caller.merged.indels.vcf.gz" into indels_filter_ch
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
  file "${vcf.simpleName}.indels.family.merged.sorted.vcf.gz" into cgi1_ch
  script:
  """
  mkdir -p tmp
  bcftools sort -O z -o ${vcf.simpleName}.indels.family.merged.sorted.vcf.gz ${vcf} -T tmp
  rm -fr tmp
  """
}


process merge_caller_snps {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/merge_vcf/snps/caller"
  input:
  file vcf from vcf1_ch.collect()
  output:
  file "*callermerged.vcf.gz" into snps_filter_ch
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
  file "${vcf.simpleName}.snps.family.merged.sorted.vcf.gz" into cgi2_ch
  script:
  """
  mkdir -p tmp
  bcftools sort -O z -o ${vcf.simpleName}.snps.family.merged.sorted.vcf.gz ${vcf} -T tmp
  rm -fr tmp
  """
}

process CGI {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/CGI"
  input:
  file vcf from cgi2_ch.collect()
  file vcf2 from cgi1_ch.collect()
  output:
  file "*.zip" into unzip_ch
  script:
  """
  module add python/anaconda/2020.07
  python $baseDir/bin/python/CGI.py
  """
}

//rename the ann and drug files using nextflow
process unzip {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/CGI"
  input:
  file zip from unzip_ch
  output:
  file "*.zip" into unzip_ch.flatten()
  script:
  """
  unzip ${zip}
  mv xx/xx yy
  """
}

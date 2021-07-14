params.vcf= "$baseDir/output/VCF_collect/*.vcf.gz"

vcf_ch = Channel. fromPath (params.vcf)
vcf_ch.into { vcf1_ch; vcf2_ch; vcf3_ch }

//This is a fake PED of brentP
process ped {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	executor 'local'
  storeDir "$baseDir/output_TRAPD/VCF_collect/merge_vcf/indels/caller"
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
  storeDir "$baseDir/output_TRAPD/VCF_collect/merge_vcf/indels/caller"
  input:
  file vcf from vcf2_ch.collect()
  output:
  file "*caller.merged.indels.vcf.gz" into fam_ch
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
  storeDir "$baseDir/output_TRAPD/VCF_collect/merge_vcf/indels/family/filter"
  input:
  file vcf from fam_ch.flatten()
  output:
  file "${vcf.simpleName}.indels.vcf.gz" into vep_ch
  script:
  """
  bcftools view -O z -v indels ${vcf} > ${vcf.simpleName}.indels.vcf.gz
  """
}


process VEP {
  executor 'slurm'
  memory { 7.GB * task.attempt }
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output_TRAPD/VCF_collect/VEP/vcf"
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
  storeDir "$baseDir/output_TRAPD/VCF_collect/VEP/slivar"
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


process merge_caller_snps {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output_TRAPD/VCF_collect/merge_vcf/snps/caller"
  input:
  file vcf from vcf1_ch.collect()
  output:
  file "*callermerged.vcf.gz" into snps_fam_ch
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

process snps_filter {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	executor 'local'
	memory { 7.GB * task.attempt }
  storeDir "$baseDir/output_TRAPD/VCF_collect/merge_vcf/snps/family/filtered"
  input:
  file vcf from snps_fam_ch.flatten()
  output:
  file "${vcf.simpleName}.snps.vcf.gz" into vep2_ch
  script:
  """
  bcftools view -O z -V indels ${vcf} > ${vcf.simpleName}.snps.vcf.gz
  """
}

process VEP2 {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  executor 'slurm'
	clusterOptions '--qos=hmem'
	queue 'hmem-512'
	memory { 30.GB * task.attempt }
  storeDir "$baseDir/output_TRAPD/VCF_collect/VEP/vcf"
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
  storeDir "$baseDir/output_TRAPD/VCF_collect/VEP/slivar"
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

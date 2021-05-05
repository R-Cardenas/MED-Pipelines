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
  storeDir "$baseDir/output/VCF_collect/merge_vcf/indels/filter"
  input:
  file vcf from indels_filter_ch.flatten()
  output:
  file "${vcf.simpleName}.indels.vcf.gz" into indels_sort_ch
  script:
  """
  bcftools view -O z -v indels ${vcf} > ${vcf.simpleName}.indels.vcf.gz
  """
}

process indels_sort {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/merge_vcf/indels/family_sorted"
  input:
  file vcf from indels_sort_ch
  output:
  file "${vcf.simpleName}.indels.family.merged.vcf.gz" into vep_ch
  script:
  """
  mkdir -p tmp
  bcftools sort -O z -o ${vcf.simpleName}.indels.family.merged.vcf.gz ${vcf} -T tmp
  rm -fr tmp

  bcftools view -g het ${vcf.simpleName}.indels.family.merged.vcf.gz > ${vcf.simpleName}_HETS.indels.family.merged.vcf.gz
  """
}

process VEP {
  executor 'slurm'
  cpus 5
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/split_vcf/VEP2"
  input:
  file vcf from vep_ch.flatten()
  output:
  file "${vcf.baseName}_VEP.vcf" into vep32_ch
  script:
  """
  /ensembl-vep/vep -i ${vcf} \
  --dir /var/spool/mail/VEP_hg38/.vep \
  -o ${vcf.baseName}_VEP.vcf \
  --cache homo_sapiens \
  --sift b \
  --polyphen b \
  --af_gnomad \
  --fasta /var/spool/mail/cgpwgs_ref/GRCh38/core_ref_GRCh38_hla_decoy_ebv/genome.fa \
  --show_ref_allele \
  --symbol \
  --verbose \
  --vcf \
  --force_overwrite \
  --fields "Uploaded_variation,Location,Allele,REF_ALLELE,Consequence,Feature_type,Feature,CDS_position,gnomAD_NFE_AF,SIFT,PolyPhen,Amino_acids,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,DOMAINS,ZYG,CLIN_SIG,SYMBOL"
  """
}


process vep_header2 {
  storeDir "$baseDir/output/VCF_collect/split_vcf/VEP2"
  input:
  file txt from vep32_ch
  output:
  file "${txt.baseName}_noheader.txt"
  script:
  """
  sed -i 's/#Uploaded_variation/Uploaded_variation/g' ${txt}
  awk '!/\\#/' ${txt} > ${txt.baseName}_noheader.txt

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

process merge_fam_snps{
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/merge_vcf/snps/family/scratch"
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

  mcp '*/0001.vcf' '#1.vcf'

  for file in *.vcf; do
    bgzip \$file
  done

  # The mcp function renames the file based on the path (where the * is )
  # this will take 0001.vcf which uses freebayes vcf info fields.
  """
}


// snps includes snps and MNPs
// rationale for mnps in this section and not indels
// mnps can be covered but whole reads length so variant callers shouldnt struggle as much as with indels
process snps_filter {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/merge_vcf/snps/family/scratch"
  input:
  file vcf from snps_filter_ch.flatten()
  output:
  file "${vcf.simpleName}.snps.vcf.gz" into snps_sort_ch
  script:
  """
  bcftools view -O z -V indels ${vcf} > ${vcf.simpleName}.snps.vcf.gz
  """
}

process snps_sort {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/merge_vcf/snps/family"
  input:
  file vcf from snps_sort_ch
  output:
  file "${vcf.simpleName}.snps.family.merged.vcf.gz" into vep2_ch
  script:
  """
  mkdir -p tmp
  bcftools sort -O z -o ${vcf.simpleName}.snps.family.merged.vcf.gz ${vcf} -T tmp
  rm -fr tmp

  bcftools view -g het ${vcf.simpleName}.snps.family.merged.vcf.gz > ${vcf.simpleName}_HET.snps.family.merged.vcf
  """
}

process VEP2 {
  executor 'slurm'
  cpus 5
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/split_vcf/VEP2"
  input:
  file vcf from vep2_ch.flatten()
  output:
  file "${vcf.baseName}_VEP.vcf" into vep3_ch
  script:
  """
  /ensembl-vep/vep -i ${vcf} \
  --dir /var/spool/mail/VEP_hg38/.vep \
  -o ${vcf.baseName}_VEP.vcf \
  --cache homo_sapiens \
  --sift b \
  --polyphen b \
  --af_gnomad \
  --fasta /var/spool/mail/cgpwgs_ref/GRCh38/core_ref_GRCh38_hla_decoy_ebv/genome.fa \
  --show_ref_allele \
  --symbol \
  --verbose \
  --vcf \
  --force_overwrite \
  --fields "Uploaded_variation,Location,Allele,REF_ALLELE,Consequence,Feature_type,Feature,CDS_position,gnomAD_NFE_AF,SIFT,PolyPhen,Amino_acids,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,DOMAINS,ZYG,CLIN_SIG,SYMBOL"
  """
}

process vep_header22 {
  storeDir "$baseDir/output/VCF_collect/split_vcf/VEP2"
  input:
  file txt from vep3_ch
  output:
  file "${txt.baseName}_noheader.txt"
  script:
  """
  sed -i 's/#Uploaded_variation/Uploaded_variation/g' ${txt}
  awk '!/\\#/' ${txt} > ${txt.baseName}_noheader.txt

  """
}

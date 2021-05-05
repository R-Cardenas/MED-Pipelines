params.vcf= "$baseDir/output/VCF_collect/indels_overlap/split_VCF/*.vcf"

vcf2_ch = Channel. fromPath (params.vcf)




process het_extract{


  storeDir "$baseDir/output/VCF_collect/indels_overlap/split_VCF/heterozygous_filter"
  input:
  file vcf from vcf2_ch
  output:
  file "${vcf.baseName}_hets.vcf" into vep_ch
  script:
  """
  vcfkeepinfo ${vcf} "AF" | vcfbreakmulti | vcffilter -f "AF = 0.5" > ${vcf.baseName}_hets1.vcf

  vcfkeepgeno ${vcf.baseName}_hets1.vcf "GL" > ${vcf.baseName}_hets2.vcf
  """
}

process VEP2{

  storeDir "$baseDir/output/VCF_collect/indels_overlap/split_VCF/heterozygous_filter/VEP"
  input:
  file vcf from vep_ch
  output:
  file "${vcf.baseName}_VEP.txt" into vep2_ch
  script:
  """
  /ensembl-vep/vep -i ${vcf} \
  --dir /var/spool/mail/VEP_hg38/.vep \
  -o ${vcf.baseName}_VEP.txt \
  --offline \
  --fasta $VEP_fasta \
  --fork 5 \
  --cache homo_sapiens \
  --sift p \
  --polyphen p \
  --variant_class \
  --af_gnomad \
  --no_stats \
  --tab \
  --show_ref_allele \
  --symbol \
  --verbose \
  --domains \
  --regulatory
  """
}

process vep_header {
  storeDir "$baseDir/output/VCF_collect/indels_overlap/split_VCF/heterozygous_filterVEP"
  input:
  file txt from vep2_ch
  output:
  file "${txt.baseName}_noheader.txt"
  script:
  """
  sed -i 's/#Uploaded_variation/Uploaded_variation/g' ${txt}
  awk '!/\\#/' ${txt} > ${txt.baseName}_noheader.txt
  """
}

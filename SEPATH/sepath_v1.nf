// This script is a rewrite of Abe's SEPATH (https://github.com/UEA-Cancer-Genetics-Lab/sepath_tool_UEA/blob/master/bin/fastq_filtering/all_SEPATH.snake)
//but using Nextflow

params.bam = "$baseDir/input/*.bam"
bam_ch = Channel .fromPath( params.bam )


// ask abe which is the big file
// need to put in contam and kraken db

process bam_filter {
  storeDir "$baseDir/output/bam_filter"
  executor 'slurm'
  input:
  file bam from bam_ch
  output:
  file "${bam.simpleName}_R1.fq.gz" into fq1_ch
  file "${bam.simpleName}_R2.fq.gz" into fq2_ch
  script:
  """
  python $baseDir/bin/python/bam_filter.py --input_bam $bam
  """
}


process trimmomatic {
  storeDir "$baseDir/output/trimmomatic"
  executor 'slurm'
  input:
  file fq1 from fq1_ch
  file fq2 from fq2_ch
  output:
  file "${fq1.simpleName}.paired.fastq" into fq1pair_ch
  file "${fq1.simpleName}.unpaired.fastq" into fq1unpair_ch
  file "${fq2.simpleName}.paired.fastq" into fq2pair_ch
  file "${fq2.simpleName}.unpaired.fastq" into fq2unpair_ch
  script:
  """
  trimmomatic PE $fq1 $fq2 ${fq1.simpleName}.paired.fastq ${fq1.simpleName}.unpaired.fastq \
  ${fq2.simpleName}.paired.fastq ${fq2.simpleName}.unpaired.fastq \
  SLIDINGWINDOW:4:20 MINLEN:35
  """
}

process bbduk {
  storeDir "$baseDir/output/bbduk"
  clusterOptions '--qos=hmem --exclusive'
  queue 'hmem-512'
  executor 'slurm'
  memory '400 GB'
  input:
  file fq1 from fq1unpair_ch
  file fq2 from fq2unpair_ch
  file fq11 from fq1pair_ch
  file fq22 from fq2pair_ch
  output:
  file "${fq1.simpleName}.unpaired.cat.bbduk.fastq" into f1unpair_ch
  file "${fq11.simpleName}.out.R1.fastq" into f1pair_ch
  file "${fq11.simpleName}.out.R2.fastq" into f2pair_ch
  script:
  """
  cat ${fq1} ${fq2} > ${fq1.simpleName}_unpaired_cat.fastq

  # single human depletion
  bbduk.sh \
  in1=${fq1.simpleName}_unpaired_cat.fastq out=${fq1.simpleName}.unpaired.cat.bbduk.fastq k=30 \
  -Xmx230g \
  ref=/var/spool/mail/bbduk_db.fa \
  mcf=0.5

  # PE bbduk
  bbduk.sh in1=$fq11 in2=$fq22 \
  out1=${fq11.simpleName}.out.R1.fastq \
  out2=${fq11.simpleName}.out.R2.fastq \
  outs=${fq11.simpleName}.out.R3.fastq \
  k=30 \
  -Xmx230g \
  ref=/var/spool/mail/bbduk_db.fa \
  mcf=0.5
  """
}

// ask abe what is it R1 and R2 we put into kraken
process kraken_paired {
  storeDir "$baseDir/output/kraken"
  clusterOptions '--qos=hmem --exclusive'
  queue 'hmem-512'
  executor 'slurm'
  memory '400 GB'
  input:
  file f1_unpair from f1unpair_ch
  file f1_pair from f1pair_ch
  file f2_pair from f2pair_ch
  output:
  file "${f1_pair.simpleName}.merge.kraken" into thresh_ch
  script:
  """
#_unpaired
kraken \
--preload \
--db /var/spool/kraken \
--threads 10 --fastq-input \
${f1_unpair} \
--output ${f1_pair.simpleName}.unpair.kraken

#paired

kraken --preload \
--db /var/spool/kraken \
--threads 10 \
--fastq-input \
--paired \
${f1_pair} \
${f2_pair} \
--output \
${f1_pair.simpleName}.pair.kraken

sleep 2m

cat ${f1_pair.simpleName}.pair.kraken ${f1_pair.simpleName}.unpair.kraken > ${f1_pair.simpleName}.merge.kraken

rm -fr *.unpair.kraken *.pair.kraken
"""
}

// kraken threshold on the concatentated files
// 20% of kmers has to meet the taxonomy
process threshold {
  clusterOptions '--qos=hmem --exclusive'
  queue 'hmem-512'
  executor 'slurm'
  memory '400 GB'
  storeDir "$baseDir/output/bam_filter"
  input:
  file kraken from thresh_ch
  output:
  file "${kraken.simpleName}.theshold.kraken​​" into report_ch
  script:
  """
  kraken-filter --db /var/spool/kraken --threshold 0.2 ${kraken} > ${kraken.simpleName}.theshold.kraken​​
   """
}

//produces readable file
process report {
  storeDir "$baseDir/output/bam_filter"
  input:
  file kraken from report_ch
  output:
  file "${kraken.simpleName}.report.final" into fq3_ch
  script:
  """
  kraken-report --db /var/spool/kraken ${kraken} > ${kraken.simpleName}.report.final
  """
}

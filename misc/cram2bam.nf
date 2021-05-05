
params.bam = "$baseDir/input/*.cram"
bam_ch = Channel .fromPath( params.bam )


process BaseRecalibrator {
  storeDir "$baseDir/output/bam"
  input:
  file cram from bam_ch
  output:
  file "${cram.simpleName}.bam"
  script:
  """
  samtools view -b -T /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes/hg19_GRCh37d5/GRCh37d5/core_ref_GRCh37d5/genome.fa -o ${cram.simpleName}.bam ${cram}
  """
}
  

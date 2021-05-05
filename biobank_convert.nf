/// This works completely 180620

/*
 * create a channel for fastq pairs
 */
params.bam = "/gpfs/afm/cg_pipelines/Datasets/PCAWG_test_WGS/PD13382a.cram"

Channel
	.fromPath( params.bam )
	.set {bam_ch}


process split {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 2
  cpus 10
	executor 'slurm'
  storeDir "$baseDir/output/split"
  input:
  file bam from bam_ch
  output:
  file "*.bam" into bam2_ch
  script:
  """
  module add samtools
  samtools split --threads 10 \
  --reference /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes/biobank/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  $bam
  """
}

process sort {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 2
	executor 'slurm'
  storeDir "$baseDir/output/sort"
  input:
  file bam from bam2_ch.flatten()
  output:
  file "${bam.baseName}.sorted.bam" into bam4_ch
  script:
  """
	module add samtools
  samtools sort $bam -o ${bam.baseName}.sorted.bam
  rm -fr $bam
  """
}

process fq {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  storeDir "$baseDir/output/fastq"
	executor 'slurm'
  input:
  file bam from bam4_ch
  output:
  file "*.gz"
  script:
  """
  module add samtools
	samtools fastq -@ 8 $bam \
	-1 ${bam.baseName}_1.fq.gz \
	-2 ${bam.baseName}_2.fq.gz
  rm -fr $bam
  """
}

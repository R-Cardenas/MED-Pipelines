/*
 * create a channel for fastq pairss
 */

read1_ch = Channel .fromFilePairs( params_fq )
read1_ch.into { read2_ch; read3_ch }

println """
         ==================================
         E X O M E - N F   P I P E L I N E
         N O  M E R G E

         Mapping and Bam Processing
         v0.1
         ===================================



         """
         .stripIndent()


process trim_galore{
	afterScript "rm -fr ${reads[0]} ${reads[1]}"
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 7
	stageInMode = 'copy' // trim_galore doesnt like sym/hardlinks.
  storeDir "$baseDir/output/trim_galore"
	input:
	tuple val(read2), file(reads) from read2_ch
	output:
	file "${reads[0].simpleName}_val_1.fq.gz" into (read5_ch, read7_ch)
	file "${reads[0].simpleName}_val_2.fq.gz" into (read10_ch, read12_ch)
	file("*.html") optional true
	script: {
	"""
	mkdir -p $baseDir/logs

	trim_galore --paired --fastqc --illumina \
	--basename ${reads[0].simpleName} \
	${reads[0]} ${reads[1]}

	# delete copied files
	rm -fr ${reads[0]} # remove the copied files to prevent memory loss
	rm -fr ${reads[1]}
	"""
}
}


process fqtools{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 7
  storeDir "$baseDir/output/trim_galore/fqtools"
	input:
	file read1 from read7_ch
	file read2 from read12_ch
	output:
	file "${read1}.yaml" into yaml_ch
	file("fqtools_WARNING_?.txt") optional true
	script:
	"""
	# Read1
	# Extract header
	fqtools -d header ${read1} | grep ":[C,A,T,G]*[+][C,A,T,G]" | head -1 > ${read1.simpleName}.txt

	### Counts lines in file1 and will repeat if it is empty
  words=`wc -l ${read1.simpleName}.txt  | awk '{print \$1}'`;
	if [ \$words -eq 0 ]
	then
	fqtools -d header ${read1} | head -1 > ${read1.simpleName}.txt
	else
	echo 'alls good!'
	fi

	# Read2
	# Extract header

  fqtools -d header ${read2} | grep ":[C,A,T,G]*[+][C,A,T,G]" | head -1 > ${read2.simpleName}.txt

	### Counts lines in file2 and will repeat if it is empty
	words=`wc -l ${read2.simpleName}.txt  | awk '{print \$1}'`;
	if [ \$words -eq 0 ]
	then
	fqtools -d header ${read2} | head -1 > ${read2.simpleName}.txt
	else
	echo 'alls good!'
	fi


	python $baseDir/bin/python/fastq2config_cgpmap.py \
	--fq1 ${read1.simpleName}.txt --fq2 ${read2.simpleName}.txt \
	--n1 ${read1} --n2 ${read2} --o ${read1}.yaml

	cp *WARNING* $baseDir/logs 2>/dev/null || :
	"""
}

process cgpMAP {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 3
  executor 'slurm'
	clusterOptions '--qos=hmem'
	queue 'hmem-512'
	memory { 130.GB * task.attempt }
	storeDir "$baseDir/output/cgpMAP/${read1.simpleName}"
  input:
	val read1 from read5_ch
	val read2 from read10_ch
	file yaml from yaml_ch.collect()
  output:
  file "${read1.simpleName}.bam" into cgp_ch
  script:
  """

	rm -fr $baseDir/output/cgpMAP/${read1.simpleName} # delete and remake path in case of retry pipeline
	mkdir -p $baseDir/output/cgpMAP/${read1.simpleName}

  name=\$(echo '${read2}' | sed -e 's/.*[/]//' -e 's/-.*//')
  BaseName=\$(basename ${read1})

  ds-cgpmap.pl  \
  -outdir $baseDir/output/cgpMAP/${read1.simpleName} \
  -r $cgpmap_genome \
  -i $cgpmap_index \
  -t 10 \
  -s \$name \
	-g \$BaseName.yaml \
  ${read1} ${read2}

	mv $baseDir/output/cgpMAP/${read1.simpleName}/*.bam \
	$baseDir/output/cgpMAP/${read1.simpleName}/${read1.simpleName}.bam

	echo 'fq1: ${read1} fq2: ${read2} bam_name: ${read1.simpleName}' >> $baseDir/logs/${projectname}_cgpmap_samples.log
	ls -l  ${read1} >> $baseDir/logs/symbolic_test_fastq.log
	ls -l  ${read2} >> $baseDir/logs/symbolic_test_fastq.log
  """
}


process sam_sort {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 7
  storeDir "$baseDir/output/BAM/sorted"
  memory { 6.GB * task.attempt }
  input:
  file bam from cgp_ch
  output:
  file "${bam}.sorted.bam" into bam_merge_ch
  script:
  """
	# create the tmp file as picard can create large temp files

	mkdir -p tmp
  picard SortSam I=${bam} O=${bam}.sorted.bam SORT_ORDER=coordinate TMP_DIR=tmp
	rm -fr tmp
  """
}

// dont forget to add singularity with python3 installed
process bam_merge {
  executor 'slurm'
  memory { 20.GB * task.attempt }
  storeDir "$baseDir/output/BAM/merge"
  input:
  file bam from bam_merge_ch.collect()
  output:
	file "*-merged.bam" into dup_ch
  script:
  """
  module add python/anaconda/2020.07/3.8
  module add samtools
	python $baseDir/bin/python/merge_bam_v2.py --bam '${bam}'
  """
}

process picard_pcr_removal {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 7
	storeDir "$baseDir/output/BAM/merge/RMD"
  input:
  file bam from dup_ch.flatten()
  output:
  file "${bam.simpleName}.rmd.bam" into (index1_ch, hs_ch, bam10_ch, bam11_ch, bam12_ch)
	file "${bam.simpleName}.log"
  script:
  """
	mkdir -p tmp
  picard MarkDuplicates I=${bam} O=${bam.simpleName}.rmd.bam M=${bam.simpleName}.log TMP_DIR=tmp
	TMP_DIR=tmp
	rm -fr tmp
  """
}

process bam_index {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 7
  memory { 6.GB * task.attempt }
	storeDir "$baseDir/output/BAM/merge/RMD"
  input:
  file bam from index1_ch
  output:
  file "${bam}.bai" into (index_3ch, index_4ch)

  script:
  """
	mkdir -p tmp
  picard BuildBamIndex \
	I=${bam} \
	O=${bam}.bai \
	TMP_DIR=tmp
	rm -fr tmp
  """
}


process hybrid_stats {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 7
  memory { 6.GB * task.attempt }
	storeDir "$baseDir/output/BAM/hybrid_stats"
  input:
  file bam from hs_ch
  output:
  file "${bam.simpleName}_hs_metrics.txt"
  script:
  """

	# create interval files from BED supplied
	picard BedToIntervalList \
	-I ${bait_interval} \
  --CREATE_INDEX true \
	-O ${bait_interval}.interval \
	-SD $genome_fasta

	picard BedToIntervalList \
	-I ${target_interval} \
  --CREATE_INDEX true \
	-O ${target_interval}.interval \
	-SD $genome_fasta


	# perform hybrid stats
	mkdir -p tmp
  picard CollectHsMetrics I=${bam} O=${bam.simpleName}_hs_metrics.txt \
  R=$genome_fasta \
  BAIT_INTERVALS=${bait_interval}.interval \
  TARGET_INTERVALS=${target_interval}.interval \
	TMP_DIR=tmp
	rm -fr tmp
  """
}

process alignment_stats{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 7
  memory { 6.GB * task.attempt }
	storeDir "$baseDir/output/BAM/alignment_stats"
	input:
	file bam from bam10_ch
	output:
	file "${bam.simpleName}_align_stats.txt" into verify_ch
	script:
	"""
	mkdir -p tmp
	picard CollectAlignmentSummaryMetrics \
  R=$genome_fasta \
	I=${bam} \
	O=${bam.simpleName}_align_stats.txt \
	TMP_DIR=tmp
	rm -fr tmp
	"""
}

process verifybamid{
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 7
	cpus 3
  executor 'slurm'
	clusterOptions '--qos=hmem'
	queue 'hmem-512'
	memory '80 GB'
//	stageInMode = 'copy' // somalier doesnt like sym/hardlinks.
	storeDir "$baseDir/output/BAM/verifyBamID"
	input:
	file bam from bam11_ch
	file idx from index_3ch.collect()
	output:
	file "${bam.simpleName}.log"
	file "${bam.simpleName}.selfSM" into con_ch
	script:
	"""
  verifyBamID --vcf $verifybamid \
	--bam ${bam} \
	--out ${bam.simpleName} \
	--maxDepth 1000 \
	--precise \
	--verbose \
	--ignoreRG \
	-ignoreOverlapPair

	rm -fr *.bam
	"""
}

process somalier{
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 7
  executor 'slurm'
	clusterOptions '--qos=hmem'
	queue 'hmem-512'
	memory { 50.GB * task.attempt }
  storeDir "$baseDir/output/BAM/somalier"
  input:
  file bam from bam12_ch.collect()
	file idx from index_4ch.collect()
  output:
	file "*.html" into som_ch
  script:
  """
	python $baseDir/bin/python/PED_file.py --bam '$bam'

	mkdir -p bin

	for f in *.bam; do
    /somalier:v0.2.11/somalier extract -d extracted/ --sites $somalier -f $genome_fasta \$f
	done

	/somalier:v0.2.11/somalier relate --ped project.PED  extracted/*.somalier

	wget -P ancestry_files https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv

	wget https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz
	tar -xzf 1kg.somalier.tar.gz

	/somalier:v0.2.11/somalier ancestry --labels ancestry_files/ancestry-labels-1kg.tsv \
	1kg-somalier/*.somalier ++ extracted/*.somalier

	rm -fr *.bam
  """
}

process multiqc{
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  maxRetries 7
	executor 'slurm'
	storeDir "$baseDir/output/multiqc"
	input:
	file somalier from som_ch
	file con from con_ch.collect()
	output:
	file "*.html"
	script:
	"""
	multiqc $baseDir

  # Creating final log files
  echo ' Pipeline cgpmap v0.3 completed
  Project: $projectname
  Time: ${nextflow.timestamp}

 trimmomatic - completed
 cgpmap - completed
 (reference = $cgpmap_genome)
 (index = $cgpmap_index)
 samsort - completed
 picard pcr removal - completed
 picard rename bam - completed
 samtools bam index - completed
 picard collect insert size - completed
 picard collect HS stats - completed
 (target_intervals = $target_interval)
 (bait_intervals = $bait_interval)
 merge bams - completed
 ' >> $baseDir/logs/${projectname}_cgpmap_log.txt

 mail -s "cgpMAP successful" aft19qdu@uea.ac.uk < $baseDir/logs/${projectname}_cgpmap_log.txt
 mail -s "cgpMAP info" aft19qdu@uea.ac.uk < $baseDir/logs/${projectname}_cgpmap_log.txt

	"""
}

workflow.onError{
	// create log files and record output
	process finish{
		script:
		""" echo ' Pipeline cgpmap v0.3 FAILED
		Project: $projectname
		Time: ${nextflow.timestamp}
	 ' >> $baseDir/logs/cgpMAP.${projectname}.failed.log

	 mail -s "cgpMAP FAIL $projectname" aft19qdu@uea.ac.uk < $baseDir/logs/cgpMAP.${projectname}.failed.log
	 """
	}
}

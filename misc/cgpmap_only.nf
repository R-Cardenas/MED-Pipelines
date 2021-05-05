/*
 * create a channel for fastq pairss
 */

// Input Reads
params.read1 = "$baseDir/input/*{1,2}.fq.gz"


read1_ch = Channel .fromFilePairs( params.read1 )



println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         N O  M E R G E

         Mapping and Bam Processing
         v0.1
         ===================================



         """
         .stripIndent()


process cgpMAP {
	stageInMode = "copy" // trim_galore doesnt like sym/hardlinks.
  input:
	tuple val(read2), file(reads) from read1_ch
  output:
  file "*.bam" into cgp_ch
  script:
  """
	ds-cgpmap.pl  \
  -r $cgpmap_genome \
  -i $cgpmap_index \
	-s ${read2} \
	-t 5 \
	${reads[0]} ${reads[1]}
  """
}

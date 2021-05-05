
// Define Tumor vs Normal Variables
tumor_ch = Channel .from (params.tumor )
normal_ch = Channel .from (params.normal )

println """\
	\
	\
	\
         ==================================
         E X O M E - N F   S O M A T I C
         cgpwxs_3.1.6.img
				 Tumour matched
				 v0.2
         ===================================
         """
         .stripIndent()

process cgpwxs2 {
	storeDir "$baseDir/output/cgpwxs"
	input:
	val tumor from tumor_ch
	val normal from normal_ch
	output:
	file "WXS_${tumor}_vs_${normal}.result.tar.gz"
	script:
	"""
  mkdir -p $baseDir/output/cgpwxs

	singularity exec --cleanenv \
	--home $baseDir \
	--bind /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes:/var/spool/ref:ro \
	--bind $baseDir/input:/var/spool/data:ro \
	/gpfs/afm/cg_pipelines/Pipelines/singularity/images/cgpwxs_3.1.6.img \
ds-cgpwxs.pl \
-reference /var/spool/ref/PCAWG_37_ICGC_Sanger_genome/core_PCAWG_hg19.tar.gz \
-annot $cgp_annot \
-snv_indel $cgp_snv_indel \
-tumour /var/spool/data/${tumor}*.cram \
-tidx /var/spool/data/${tumor}*.cram.crai \
-normal /var/spool/data/${normal}*.cram \
-nidx /var/spool/data/${normal}*.cram.crai \
-exclude NC_007605,hs37d5,GL% \
-outdir $baseDir/output/hg38_decoy/cgpwxs \
-sp "Homo sapiens" \
-assembly ${build}
	"""
}

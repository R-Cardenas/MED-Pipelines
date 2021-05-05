
params.tumor = ["PD13382a.cram","PD13399a.cram"]
params.normal = ["PD13382b.cram","PD13399b.cram"]

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
	storeDir "$baseDir/output/cgpwxs/${tumor}_vs_${normal}"
	input:
	val tumor from tumor_ch
	val normal from normal_ch
	output:
	file "WXS_${tumor}_vs_${normal}.result.tar.gz"
	script:
	"""
  mkdir -p $baseDir/output/cgpwxs/${tumor}_vs_${normal}

	singularity exec --cleanenv \
	--home $baseDir \
	--bind /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes:/var/spool/ref:ro \
	--bind $baseDir/input:/var/spool/data:ro \
	/gpfs/afm/cg_pipelines/Pipelines/singularity/images/cgpwxs_3.1.6.img \
ds-cgpwxs.pl \
-reference $cgp_ref \
-annot $cgp_annot \
-snv_indel $cgp_snv_indel \
-tumour /var/spool/data/${tumor}*.bam \
-tidx /var/spool/data/${tumor}*.bam.bai \
-normal /var/spool/data/${normal}*.bam \
-nidx /var/spool/data/${normal}*.bam.bai \
-exclude NC_007605,hs37d5,GL% \
-outdir $baseDir/output/cgpwxs/${tumor}_vs_${normal} \
-sp "Homo sapiens" \
-assembly "GRCh37"
	"""
}

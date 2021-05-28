
// Define Tumor vs Normal Variables
tumor_ch = Channel .from (params_tumor )
normal_ch = Channel .from (params_normal )

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
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 7
  executor 'slurm'
	clusterOptions '--qos=hmem'
	queue 'hmem-512'
	memory { 30.GB * task.attempt }
	storeDir "$baseDir/output/cgpwxs/${tumor}_vs_${normal}"
	input:
	val tumor from tumor_ch
	val normal from normal_ch
	output:
	file "WXS_${tumor}_vs_${normal}.result.tar.gz" into untar_ch
	script:
	"""
	rm -fr $baseDir/output/cgpwxs/${tumor}_vs_${normal} # during retry if file exist causes fail
  mkdir -p $baseDir/output/cgpwxs/${tumor}_vs_${normal}

	singularity exec --cleanenv \
	--home $baseDir \
	--bind /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes:/var/spool/ref:ro \
	--bind $baseDir/output/BAM/merge/RMD/:/var/spool/data:ro \
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

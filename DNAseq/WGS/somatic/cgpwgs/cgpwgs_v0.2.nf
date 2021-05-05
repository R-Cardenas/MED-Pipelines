// This uses config for cgpWXS (/DNASeq/Exome/somatic/cgpwxs/...) 10/12/20

params.tumor = ["PD13382a","PD13389a","PD13399a"]
params.normal = ["PD13382b","PD13389b","PD13399b"]


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
	executor 'slurm'
	clusterOptions '--qos=hmem'
	queue 'hmem-512'
	memory '400 GB'
	storeDir "$baseDir/output/cgpwgs22"
	input:
	val tumor from tumor_ch
	val normal from normal_ch
	output:
	file "${tumor}_vs_${normal}/WXS_${tumor}_vs_${normal}.result.tar.gz"
	script:
	"""
	python $baseDir/bin/python/cgpwgs_run.py --tumor '${tumor}' --normal '${normal}' --assembly 'hg19' --home '$baseDir' \
	--reference '$cgp_ref' --annotation '$cgp_annot' --snv_indel '$cgp_snv_indel' \
	--cgp_QC '$cgp_QC' --cgp_CNV '$cgp_CNV'
	"""
}

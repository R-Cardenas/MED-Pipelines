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
	storeDir "$baseDir/output/hg38_decoy/cgpwgs/${tumor}_vs_${normal}"
	input:
	val tumor from tumor_ch
	val normal from normal_ch
	output:
	file "${tumor}_vs_${normal}/WXS_${tumor}_vs_${normal}.result.tar.gz"
	script:
	"""
mkdir -p $baseDir/output/cgpwgs

singularity exec --cleanenv \
--home $baseDir \
--bind /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes:/var/spool/ref:ro \
--bind $baseDir/input:/var/spool/data:ro \
/gpfs/afm/cg_pipelines/Pipelines/singularity/images/cgpwgs_2.1.0.img \
ds-cgpwgs.pl \
-reference $cgp_ref \
-annot $cgp_annot \
-snv_indel $cgp_snv_indel \
-cnv_sv $cgp_CNV \
-qcset $cgp_QC \
-tumour /var/spool/data/${tumor}*.bam \
-tidx /var/spool/data/${tumor}*.bam.bai \
-normal /var/spool/data/${normal}*.bam \
-nidx /var/spool/data/${normal}*.bam.bai \
-exclude NC_007605,hs37d5,GL% \
-outdir $baseDir/output/cgpwgs/${tumor}_vs_${normal} \
-sp "Homo sapiens" \
-assembly "GRCh37"
	"""
}

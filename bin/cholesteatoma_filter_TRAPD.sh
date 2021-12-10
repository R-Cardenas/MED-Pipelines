
# first filter the vep fields to replicate what we do in R
# This needs to be done in two round as vep data is in the CSQ

# prior to TRAPD.nf you need this
bcftools norm -m -any -f /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes/cgpwgs_ref/GRCh38/core_ref_GRCh38_hla_decoy_ebv/genome.fa ALL_chole_b1_b2.vcf > ALL_chole_b1_b2_noGQ.vcf



singularity shell --home $PWD /gpfs/afm/cg_pipelines/Pipelines/singularity/images/vep.simg

filter_vep --input_file ALL_chole_b1_b2_noGQ.vcf --format vcf --output_file ALL_chole_b1_b2_noGQ.vcf_VEP_filtered2.vcf \
--filter "Conservation >= 0.1" \
--filter "gnomad_popmax_af != -1 and gnomad_popmax_af != -1" \
--filter "SIFT is not tolerated and PolyPhen is not benign" \
--filter "CDS_position != NA"
#--filter "EUR_AF <= 0.01" \
#--filter "gnomAD_NFE_AF <= 0.01" \
#--filter "EUR_AF != NA" \
#--filter "gnomAD_NFE_AF != NA" \


#I have checked using bcftools and there is no entries less than 0 (ie -1)
#bcftools view -i 'gnomad_popmax_af<0' choleB2-all-merge_filter1.vcf.gz

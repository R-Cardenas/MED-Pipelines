singularity shell --home $PWD --bind /local/aft19qdu/cgpwxs/genomes:/var/spool/mail dockstore-cgpmap-3.0.0.simg

singularity shell --home $PWD --bind /local/aft19qdu/cgpwxs/genomes:/var/spool/mail ../cgpmap_FQ/dockstore-cgpmap-3.0.0.sim

ds-cgpmap.pl \
 -r /var/spool/mail/core_ref_GRCh37d5.tar.gz \
 -i /var/spool/mail/bwa_idx_GRCh37d5.tar.gz \
 -s S1B \
 -o /local/aft19qdu/cgpwxs/cgpmap_TRIMMED/output/1B \
 input/S1B_UKDHE190339-A28-A61_HGV2GDSXX_L3_1_1.trim.fq.gz input/S1B_UKDHE190339-A28-A61_HGV2GDSXX_L3_1_2.trim.fq.gz


  singularity shell --home $PWD --bind /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes/:/var/spool/mail /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes/cgpmap_cgpwxs_NEW_2021/dockstore-cgpwxs_3.1.7.sif
ds-cgpmap.pl \
 -r /var/spool/mail/core_ref_GRCh37d5.tar.gz \
 -i /var/spool/mail/bwa_idx_GRCh37d5.tar.gz \
 -s S1A \
 -o /local/aft19qdu/cgpwxs/cgpmap_fastp/1A \
S1A_fastp_1.fq.gz S1A_fastp_2.fq.gz

ds-cgpmap.pl \
 -r /var/spool/mail/core_ref_GRCh37d5.tar.gz \
 -i /var/spool/mail/bwa_idx_GRCh37d5.tar.gz \
 -s S1A \
 -o /local/aft19qdu/cgpwxs/cgpmap_fastp/2A \
S2A_fastp_1.fq.gz S2A_fastp_2.fq.gz

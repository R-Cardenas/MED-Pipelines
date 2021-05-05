/// This works completely 180620

 /*
  * create a channel for fastq pairs
  */
params.bam = "input/*.vcf"

 Channel
         .fromPath( params.bam )
         .set {bam_ch}


 println """\
         \
         \
         \
          ==================================
          E X O M E - N F   P I P E L I N E
          Freebayes - Germline
          v0.2
          ===================================



          """
          .stripIndent()




 process leftalign {
         errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
         maxRetries 6
   storeDir "$baseDir/output/step01"
         executor 'slurm'
         memory '5 GB'
   input:
   file bam from bam_ch
   output:
   file "${bam.simpleName}.norm.vcf" into vcf_ch
   script:
   """
         bcftools norm -m -any -f /var/spool/mail/cgpwgs_ref/GRCh38/core_ref_GRCh38_hla_decoy_ebv/genome.fa ${bam} > ${bam.simpleName}.norm.vcf
   """
}

process liftover {
              stageInMode = 'copy'
              //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
              //maxRetries 6
        storeDir "$baseDir/output/liftover"
              executor 'slurm'
              memory '5 GB'
        input:
        file bam from vcf_ch
        output:
        file "${bam.simpleName}.vcf.gz" into vep_ch
              file "${bam.simpleName}.reject.vcf.gz"
        script:
        """


      gatk LiftoverVcf -I ${bam} -O ${bam.simpleName}.vcf.gz -R $baseDir/genome/hg19.fa \
      --CHAIN $baseDir/hg38ToHg19.over.chain --REJECT ${bam.simpleName}.reject.vcf.gz
        """
}

process VEP{
        executor 'slurm'
        memory '15 GB'
        cpus 5
//        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
//              maxRetries 6
        storeDir "$baseDir/output/VEP"
        input:
        file vcf from vep_ch
        output:
        file "${vcf.baseName}_VEP.vcf"
        script:
        """
        /ensembl-vep/vep -i ${vcf} \
        --dir /var/spool/mail/VEP_hg19/.vep \
        -o ${vcf.baseName}_VEP.vcf \
        --cache homo_sapiens \
        --sift b \
        --everything \
        --polyphen b \
        --af_gnomad \
        --show_ref_allele \
        --symbol \
        --verbose \
        --vcf \
        --cache \
        --fork 4 --offline \
        --format vcf \
        --use_given_ref \
        --port 3337 \
        --force_overwrite \
        --fasta $baseDir/genome/hg19.fa \
        --plugin CADD,ExAC,dbNSFP,GERP
        #--fields "Uploaded_variation,Location,Allele,REF_ALLELE,Consequence,Feature_type,Feature,CDS_position,gnomAD_NFE_AF,SIFT,PolyPhen,Amino_acids,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,DOMAINS,ZYG,CLIN_SIG,SYMBOL,CADD,Conservation"
        """
}

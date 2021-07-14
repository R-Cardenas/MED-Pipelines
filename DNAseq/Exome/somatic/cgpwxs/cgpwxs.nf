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
  clusterOptions '--qos=hmem --exclusive'
  queue 'hmem-512'
  storeDir "$baseDir/output/cgpwxs/${tumor}_vs_${normal}"
  input:
	val tumor from tumor_ch
  val normal from normal_ch
  output:
  file "WXS_${tumor}_vs_${normal}.result.tar.gz" into untar_ch
  script:
        """
	module add python/anaconda/2020.07
	rm -fr $baseDir/output/cgpwxs/${tumor}_vs_${normal} # remove leftover files if retried
	mkdir -p $baseDir/output/cgpwxs/${tumor}_vs_${normal}

  python $baseDir/bin/python/run_cgpwxs.py \
    --home $baseDir \
    --tumor $tumor \
    --normal $normal \
    --annotation $cgp_annot \
    --reference $cgpmap_genome \
    --snv_indel $cgp_snv_indel

        """
}

process untar {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  maxRetries 7
  executor 'local'
  input:
	file tar from untar_ch
  output:
  file "XXX" into untar_ch
  script:
  """
	tar -xzf ${tar}
  """
}

# Extract vcf file and its associated INFO output as Tab-sep
# This works only with VEP annotated VCF files.

library(VariantAnnotation)
library(tidyverse)
library(ensemblVEP)
library(Repitools)
library("fs")
library(optparse)

#------------ set argparse -----------------

option_list = list(
  make_option(c("-v", "--vcf"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="dataset file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# set filtered variable
vcf <- as.character(opt$vcf)
outputname <- as.character(opt$output)

#------------ create function -----------------

vcf2tab <- function(vcf_path){
  print(vcf_path)
  vcf <- readVcf(vcf_path, "hg38")
  vcf2 <- vcf %>%
    unique()

  # convert VEP CSQ to GRange
  csq <- parseCSQToGRanges(vcf2) %>%
    unique() %>%
    as.data.frame() %>%
    rownames_to_column()


  # Extract the VCF INFO field out
  vcf_info <- info(vcf2)
  vcf_info2 <- cbind(vcf_info@rownames,as.data.frame(vcf_info)) %>%
    dplyr::rename(rowname = "vcf_info@rownames")

  # Merge the FORMAT and INFO
  vcf_df <- inner_join(vcf_info2,csq,by="rowname")

  # Extract filename from variable and add to df
  filename <- path_file(vcf_path)
  filename2 <- paste0(filename,".tsv")
  vcf_df1 <- vcf_df %>%
    mutate(filenames = filename)

  vcf_df2 <- vcf_df1 %>%
    as.data.frame() %>%
    unnest(cols = c(AC,AF.x,AO))


  write_tsv(vcf_df2,file = outputname)
}

#------------ fapply function  -----------------

vcf2tab(vcf)

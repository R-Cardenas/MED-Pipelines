#!/usr/bin/env Rscript
# Script for germlinePipeline RPC 200521
# Will exctract the vcf count file produced by vcf-annotate (vcftools)


library("optparse")
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)

# specify our desired options
# by default ArgumentParser will add an help option


option_list = list(
  make_option(c("-v", "--vcf"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-f", "--filtered"), type="character", default="NA",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-w", "--wildcard"), type="character", default="*.count",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-c", "--variantcaller"), type="character", default="NA",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-o", "--outputname"), type="character", default="NA",
              help="output file name [default= %default]", metavar="character"), # This is for the name of the FILE
  make_option(c("-O", "--outputname2"), type="character", default="NA",
              help="output file name [default= %default]", metavar="character")  # This is adding information into the table
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# set filtered variable
filtered <- as.character(opt$filtered)
wildcard <- as.character(opt$wildcard)
variantcaller <- as.character(opt$variantcaller)
listFiles <- as.character(opt$vcf)
print(paste0("The following files found in Dir: \n   ", listFiles))


# input files are within the current Dir
files <- list.files(path=".", pattern = wildcard,full.names = T)
print(files)

# --------------- Data processing  ---------------

df <- data_frame(filenames = files) %>%
  mutate(file_contents = map(filenames,
                             ~read_table2(.,col_names=FALSE))) %>%
  unnest() %>%
  mutate(filenames = basename(filenames)) %>%
  mutate(variantCaller = variantcaller) %>%
  mutate(filtered = filtered) %>%
  rename(n = X1) %>%
  rename(type2 = X2) %>%
  mutate(type2 = str_extract(type2,"(?<=TYPE=).*")) %>%
  mutate(file = opt$outputname2)

outputname <- as.character(opt$outputname)
write_csv(df,outputname,col_names=T)
print("Process finished.")

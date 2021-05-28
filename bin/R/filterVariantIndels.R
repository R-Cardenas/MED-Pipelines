## This will filter the snps files and then save them
# Part of the MED-genomics varianProcessing pipeline
library(tidyverse)
library(optparse)

#------------ set argparse -----------------

option_list = list(
  make_option(c("-a", "--af"), type="double", default=0.05,
              help="dataset file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# set filtered variable
af <- as.double(opt$af)

#------------ sFilter and save -----------------

files_indels <- list.files(path=".", pattern = '*.tsv')


file_path_indels_ours <- data_frame(filenames = files_indels) %>%
  mutate(file_conetents = map(filenames,
                              ~read_tsv(.) %>%
                                data.frame %>%
                                type.convert)) %>%
  unnest() %>%
  mutate(samples = str_extract(filenames, "^.{10}")) %>%
  mutate(gnomAD_NFE_AF=replace(gnomAD_NFE_AF, gnomAD_NFE_AF=="-", NA)) %>%
  type.convert() %>%
  filter(gnomAD_NFE_AF <= af)

write_csv(file_path_indels_ours,path = "ALL_data_indels.csv",col_names = TRUE)

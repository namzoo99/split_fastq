library(optparse)
library(tidyverse)
library(stringr)
library(reshape2)


option_list <- list(
  make_option(c("--raw_manifest"), type="character", default=NULL, help="legacy_file_sheet"),
  make_option(c("--previous_output_path"), type="character", default=NULL, help="output of 'read_fq_header' and 'split_fq'"),
  make_option(c("--output_path"), type="character", default=NULL, help="output_path")
)

opt = parse_args(OptionParser(option_list=option_list))

################################################################################
# Define input
################################################################################
# params.dna_legacy_file_sheet
raw_fq_manifest <- read_csv(opt$raw_manifest) %>% select(c(-RGID, -RGPU))

# output of process 'read_fq_header'
fq_header <- list.files(opt$previous_output_path, recursive = TRUE, full.names = TRUE, pattern = "*_RGPU_fq*") %>% 
  map_df(~read_tsv(file=., col_names = "fq_header", id="file_name")) %>% unique

# --split_fq_publish_dir
split_fq_publish_dir <- opt$output_path



################################################################################
# Make manifest file according to RGPU/RGID
################################################################################
# output of process 'read_fq_header'
# Check whether the input FQs are multiplexed, non-multiplexed, or both.
header <- fq_header %>% 
  mutate(file_name = sapply(str_split(file_name, "/"), '[', 13),
         aliquot_barcode = sapply(str_split(file_name, "_"), '[', 1))

is_fqs_multiplexing <- unique(header) %>% group_by(aliquot_barcode) %>% summarise(N=n())

# Use the same manifest name as the raw manifest file
name_of_manifest <- opt$raw_manifest %>% str_split("/") %>% .[[1]] %>% str_subset(".csv")

# Write a new manifest file that reflects the status of multiplexing
all_fq_ino <- header %>% 
  mutate(RGPU = sapply(str_split(fq_header, ":"), function(x){paste0(x[3], ".", x[4])}),
         RGID = sapply(str_split(fq_header, ":"), function(x){paste0(str_sub(x[3], 1, 5), ".", x[4])})) %>% 
  select(c(-file_name, -fq_header)) %>% unique

raw_fq_manifest %>%
          left_join(all_fq_ino, by=c("aliquot_barcode")) %>% relocate(RGID,RGPL,RGPU,RGLB,RGDT,RGCN,FQ1,FQ2,action, .after = last_col()) %>% 
          write_csv(paste0(split_fq_publish_dir, Sys.Date(), "-", name_of_manifest), na="")



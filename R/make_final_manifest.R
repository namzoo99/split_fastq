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

# output of process 'split_fq'
splitted_fq_path <- list.files(opt$previous_output_path, recursive = TRUE, full.names = TRUE, pattern = "*_splitted_fq_path.txt") %>% 
  map_df(~read_tsv(file=., col_names = "tmppath")) %>% unique

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
if(header[header$aliquot_barcode %in% is_fqs_multiplexing$aliquot_barcode,] %>% group_by(aliquot_barcode) %>% summarise(N=n()) >= 2){

      if(length(unique(is_fqs_multiplexing$N)) == 1 & length(unique(header)) == 1){
        
      all_fq_ino <- header %>% 
  mutate(RGPU = sapply(str_split(fq_header, ":"), function(x){paste0(x[3], ".", x[4])}),
         RGID = sapply(str_split(fq_header, ":"), function(x){paste0(str_sub(x[3], 1, 5), ".", x[4])})) %>% 
  select(c(-file_name, -fq_header)) %>% unique

      raw_fq_manifest %>%
          left_join(all_fq_ino, by=c("aliquot_barcode")) %>% relocate(RGID,RGPL,RGPU,RGLB,RGDT,RGCN,FQ1,FQ2,action, .after = last_col()) %>% 
          write_csv(paste0(split_fq_publish_dir, Sys.Date(), "-", name_of_manifest), na="")
    
      }else if(length(unique(is_fqs_multiplexing$N)) == 2 & length(unique(header)) == 2){
        
        splittd_fq <- splitted_fq_path[str_detect(splitted_fq_path$tmppath, "/nxf"),] %>% 
  mutate(splitted_fqs = sapply(str_split(tmppath, "/"), '[[', 6),
         aliquot_barcode = sapply(str_split(splitted_fqs, "_"), '[[', 1)) %>% 
  mutate(RGPU = sapply(str_split(splitted_fqs, "_"), function(x){paste0(x[2], ".", x[3])}),
         RGID = sapply(str_split(splitted_fqs, "_"), function(x){paste0(str_sub(x[2], 1, 5), ".", x[3])})) %>% 
  mutate(FQ = paste0(opt$previous_output_path, "/gdc_fastq_splitter/", splitted_fqs)) %>% 
  select(c(-tmppath, -splitted_fqs))

  multiplexed_fq_info <- left_join(splittd_fq[str_detect(splittd_fq$FQ, "_R1.fq.gz"),],
                          splittd_fq[str_detect(splittd_fq$FQ, "_R2.fq.gz"),], by=c("aliquot_barcode", "RGPU", "RGID"),
                          suffix = c("1", "2"))

        raw_fq_manifest %>%
          left_join(multiplexed_fq_info, by=c("aliquot_barcode")) %>% relocate(RGID,RGPL,RGPU,RGLB,RGDT,RGCN,FQ1,FQ2,action, .after = last_col()) %>% 
          write_csv(paste0(split_fq_publish_dir, Sys.Date(), "-", name_of_manifest), na="")
    
      }else{

        all_fq_ino <- header %>% 
  mutate(RGPU = sapply(str_split(fq_header, ":"), function(x){paste0(x[3], ".", x[4])}),
         RGID = sapply(str_split(fq_header, ":"), function(x){paste0(str_sub(x[3], 1, 5), ".", x[4])})) %>% 
  select(c(-file_name, -fq_header)) %>% unique
        
        splittd_fq <- splitted_fq_path[str_detect(splitted_fq_path$tmppath, "/nxf"),] %>% 
  mutate(splitted_fqs = sapply(str_split(tmppath, "/"), '[[', 6),
         aliquot_barcode = sapply(str_split(splitted_fqs, "_"), '[[', 1)) %>% 
  mutate(RGPU = sapply(str_split(splitted_fqs, "_"), function(x){paste0(x[2], ".", x[3])}),
         RGID = sapply(str_split(splitted_fqs, "_"), function(x){paste0(str_sub(x[2], 1, 5), ".", x[3])})) %>% 
  mutate(FQ = paste0(opt$previous_output_path, "/gdc_fastq_splitter/", splitted_fqs)) %>% 
  select(c(-tmppath, -splitted_fqs))

  multiplexed_fq_info <- left_join(splittd_fq[str_detect(splittd_fq$FQ, "_R1.fq.gz"),],
                          splittd_fq[str_detect(splittd_fq$FQ, "_R2.fq.gz"),], by=c("aliquot_barcode", "RGPU", "RGID"),
                          suffix = c("1", "2"))

        rbind(raw_fq_manifest[!raw_fq_manifest$aliquot_barcode %in% multiplexed_fq_info$aliquot_barcode,] %>%
                left_join(all_fq_ino, by=c("aliquot_barcode")) %>% relocate(RGID,RGPL,RGPU,RGLB,RGDT,RGCN,FQ1,FQ2,action, .after = last_col()),
              raw_fq_manifest[raw_fq_manifest$aliquot_barcode %in% multiplexed_fq_info$aliquot_barcode,] %>% select(c(-FQ1, -FQ2)) %>%
                left_join(multiplexed_fq_info, by=c("aliquot_barcode")) %>% relocate(RGID,RGPL,RGPU,RGLB,RGDT,RGCN,FQ1,FQ2,action, .after = last_col())) %>% 
          write_csv(paste0(split_fq_publish_dir, Sys.Date(), "-", name_of_manifest), na="")
  
      }

    }else{
  
  file.create("multiplexed_FQ_not_splitted")
  
}

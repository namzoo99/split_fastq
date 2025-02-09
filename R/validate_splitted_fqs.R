library(tidyverse)
library(optparse)


option_list <- list(
  make_option(c("--previous_output_path"), type="character", default=NULL, help="output of 'read_fq_header' and 'split_fq'")
)

opt = parse_args(OptionParser(option_list=option_list))


################################################################################
# Validate whether the FASTQ files were split correctly
# 1. Compare RGPU value
################################################################################
# Define input files
rawFQ_header <- list.files(opt$previous_output_path, 
                    pattern="RGPU_fq1.txt", recursive = TRUE, full.names = TRUE) %>% 
  map_df(~read_tsv(., id="file_name", col_names = "header"))

splittedFQ_header <- list.files("./", 
                    pattern="fq_RGPU.txt", recursive = TRUE, full.names = TRUE) %>% 
  map_df(~read_tsv(., id="file_name", col_names = "header"))

#
check_condition <- function(condition, filename) {
  tryCatch({
    stopifnot(condition)
  }, error = function(e) {
    writeLines("Error: Condition failed", filename)
    message("Error detected. Log file created: ", filename)
  })
}

# Extract RGPU value from 'splittedFQ_header'
splittedFQ_header2 <- splittedFQ_header %>% 
  mutate(file_name = str_replace(file_name, paste0(c("_R1.fq_RGPU.txt", "_R2.fq_RGPU.txt"), collapse = "|"), ""),
         RGPU = sapply(str_split(file_name, "_"), function(x){paste0(x[2],"_", x[3])})) %>% unique %>% data.frame()

# Check whether both R1 and R2 have same RGPU for each RGPU
#stopifnot(nrow(splittedFQ_header) == nrow(splittedFQ_header2)*2)
check_condition(nrow(splittedFQ_header) == nrow(splittedFQ_header2)*2, "failed.flag")

# Verify whether the RGPU value in the 'splitted fq file name' matches the actual 'RGPU value in splitted file'
#stopifnot(apply(splittedFQ_header2, 1, function(x){str_detect(x[1], x[3])}))
check_condition(apply(splittedFQ_header2, 1, function(x){str_detect(x[1], x[3])}), "failed.flag")

# Check whether the RGPU values in the 'raw fq file' match those in the 'splitted fq files'
# stopifnot(nrow(rawFQ_header) ==
#             left_join(rawFQ_header %>% mutate(RGPU = sapply(str_split(header, ":"), function(x){paste0(x[3], "_", x[4])})),
#                       splittedFQ_header2, by="RGPU") %>% nrow)

check_condition(nrow(rawFQ_header) ==
                  left_join(rawFQ_header %>% mutate(RGPU = sapply(str_split(header, ":"), function(x){paste0(x[3], "_", x[4])})),
                            splittedFQ_header2, by="RGPU") %>% nrow, "failed.flag")

################################################################################
# Validate whether the FASTQ files were split correctly
# 2. Compare # of reads
################################################################################
# Define input files
rawFQ_n_read <- list.files(opt$previous_output_path, 
                    pattern="_fq1_n_of_read.txt", recursive = TRUE, full.names = TRUE) %>% 
  map_df(~read_tsv(., id="file_name", col_names = "N"))

splittedFQ_n_read  <- list.files("./", 
                         pattern="fq_n_read.txt", recursive = TRUE, full.names = TRUE) %>% 
  map_df(~read_tsv(., id="file_name", col_names = "N"))


# Edit 'splittedFQ_n_read'
splittedFQ_n_read2 <- splittedFQ_n_read %>% 
  mutate(file_name = str_replace(file_name, paste0(c("_R1.fq_n_read.txt", "_R2.fq_n_read.txt"), collapse = "|"), "")) %>% unique


# Check whether the read counts match between R1 and R2 files with the same RGPU value
#stopifnot(nrow(splittedFQ_n_read) == nrow(splittedFQ_n_read2)*2)
check_condition(nrow(splittedFQ_n_read) == nrow(splittedFQ_n_read2)*2, "failed.flag")

# Check whether the total read count of the 'splitted fq files' matches the read count of the 'raw fq file'
#stopifnot(apply(splittedFQ_n_read2 %>% select(N), 2, sum) == rawFQ_n_read$N)
check_condition(apply(splittedFQ_n_read2 %>% select(N), 2, sum) == rawFQ_n_read$N, "failed.flag")

#
file.create("success.flag")

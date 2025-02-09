# split_fastq
## This workflow inserts RGID and RGPU values into the manifest, as well as splits the FASTQ file if it contains multiple RGIDs.

## I have only tested it on WGS & RNA-seq (Illumina-TruSeq) data
If you are using other types of sequencing data, please validate it further

## BEFORE RUNNING, PLEASE READ [README](https://github.com/kmhernan/gdc-fastq-splitter) in 'gdc-fastq-splitter' github
_gdc-fastq-splitter currently only support non-interleaved fastq files_

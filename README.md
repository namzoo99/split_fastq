# split_fastq
## This workflow inserts RGID and RGPU values into the manifest, as well as splits the FASTQ file if it contains multiple RGIDs.

## I have only tested it on WGS & RNA-seq (Illumina-TruSeq) data
If you are using other types of sequencing data, please validate it further

## BEFORE RUNNING, PLEASE READ [README](https://github.com/kmhernan/gdc-fastq-splitter) in 'gdc-fastq-splitter' github
_gdc-fastq-splitter currently only support non-interleaved fastq files_

## Example input manifest file
you dont need to fill up *RGID*, *RGPU* columns.
| aliquot_barcode                          | source_barcode | sample_barcode                      | patient_barcode       | sample_type | tumor_or_normal | sequence_type | gender | RGID | RGPL      | RGPU | RGLB       | RGDT | RGCN | FQ1                                                                 | FQ2                                                                 | action |
|------------------------------------------|---------------|-------------------------------------|----------------------|-------------|----------------|---------------|--------|------|-----------|------|------------|------|------|---------------------------------------------------------------------|---------------------------------------------------------------------|--------|
| KUBLCA-SNU16-0001-TPX-B01-5ON806         | SNU16         | KUBLCA-SNU16-0001-TPX-B01          | KUBLCA-SNU16-0001    | TP          | tumor          | CRS           | XX     |      | ILLUMINA  |      | TN2403D0100 |      | CBM  | SNU17_circle-seq_Exo_1.fastq.gz | SNU17_circle-seq_Exo/SNU17_circle-seq_Exo_2.fastq.gz | run    |
| DRUGBR-BT474-0001-TPX-B02-WGS-6GC324     | BT474         | DRUGBR-BT474-0001-TPX-B02          | DRUGBR-BT474-0001    | TP          | tumor          | WGS           | XX     |      | ILLUMINA  |      | TN2409D2511 |      | CBM  | Sample_BT474-R/BT474-R_1.fq.gz | BT474-R_2.fq.gz | run    |
| DRUGBR-LK2XX-0001-TPX-A01-1MA008         | LK2XX         | DRUGBR-LK2XX-0001-TPX-A01          | DRUGBR-LK2XX-0001    | TP          | tumor          | WGS           | XY     |      | ILLUMINA  |      | TN2403D0820 |      | CBM  | Sample_LK-2/LK-2_1.fq.gz |LK-2_2.fq.gz | run    |
| DRUGBR-UMS47-0001-TPX-B02-WGS-9WZ921     | UMS47         | DRUGBR-UMS47-0001-TPX-B02          | DRUGBR-UMS47-0001    | TP          | tumor          | WGS           | XY     |      | ILLUMINA  |      | TN2410D2673 |      | CBM  | UMSCC47-R_1.fq.gz | UMSCC47-R_2.fq.gz | run    |
| DRUGBR-BT474-0001-TPX-B01-WGS-9XY448     | BT474         | DRUGBR-BT474-0001-TPX-B01          | DRUGBR-BT474-0001    | TP          | tumor          | WGS           | XX     |      | ILLUMINA  |      | TN2410D2674 |      | CBM  | BT474-P_1.fq.gz | BT474-P_2.fq.gz | run    |

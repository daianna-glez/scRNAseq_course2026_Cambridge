#!/bin/bash

cellranger count \
  --id ETV6_RUNX1_rep1\
  --transcriptome "Data/references/GRCh38_chr21_index"\
  --fastqs "Data/references/fastq" \
  --sample "SRR9264343"\
  --create-bam "false"\
  --localcores "7" \
  --localmem "24"

#!/bin/bash

# change to directory with the FASTA and GTF files
cd Data/references/

# run mkref
cellranger mkref\
  --fasta "Homo_sapiens.GRCh38.dna.chromosome.21.fa" \
  --genes "gencode.v41.primary_assembly.annotation.chr21.gtf" \
  --genome "GRCh38_chr21_index" \
  --nthreads "7"


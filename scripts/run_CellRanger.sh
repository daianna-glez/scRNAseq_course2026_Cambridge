
## Run CellRanger

# ______________________________________________________________________________
#                         1.  Index the reference genome
# ______________________________________________________________________________
## Current directory
pwd

## Files/subdirs in current directory
ls

ls - lh
# total 8            <-    total disk space used by the listed files

# Permissions   num files    user      group   size  last modified date  file name

# drwxr-xr-x@          14 daianna      staff   448B          5 Apr 16:11 CellRanger_Outputs
# drwxr-xr-x@           8 daianna      staff   256B         15 Apr 23:17 fastq
# drwxr-xr-x@           6 daianna      staff   192B         15 Apr 19:00 references
# -rw-r--r--@           1 daianna      staff   341B          5 Apr 02:42 sample_sheet.tsv
# drwxr-xr-x@           4 daianna      staff   128B          5 Apr 02:43 SRR9264343sub


## Files in a specific subdirectory
ls Data/

## FASTA human reference for chr 21 and corresponding tx anno GFT file from Ensembl
ls Data/references

# ______________________________________________________________________________
## Explore FASTA file
less Data/references/Homo_sapiens.GRCh38.dna.chromosome.21.fa

# We find the header containing the name of the chr
# Followed by letters containing the nucleotides.
# At the start of human chr21 we don't know the sequence and have Ns
# If you scroll down you'll eventually get nts

# >chr21 21
# NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
# NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
# NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

## Exit less pressing "q" for "quit"

## To look at nts
tail -n1000 Data/references/Homo_sapiens.GRCh38.dna.chromosome.21.fa


## Explore GTF file containing tx anno
less Data/references/gencode.v41.primary_assembly.annotation.chr21.gtf

# - Metadata commented with hashes
# - Then each line for a gene/tx, giving its position in that chr and feature data (ids, version, type)

##description: evidence-based annotation of the human genome (GRCh38), version 41 (Ensembl 107)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2022-05-12

# chr21   HAVANA  gene    5086740 5087232 .       +       .       gene_id "ENSG00000289905"; gene_version "1"; gene_type "lncRNA"; gene_name "ENSG00000289905"; level 2;
# chr21   HAVANA  transcript      5086740 5087232 .       +       .       gene_id "ENSG00000289905"; gene_version "1"; transcript_id "ENST00000701545"; transcript_version "1"; gene_type "lncRNA"; gene_name "ENSG00000289905"; transcript_type "lncRNA"; transcript_name "ENST00000701545"; level 2; tag "basic"; tag "Ensembl_canonical"; tag "TAGENE";

## Look first line
head -1
head -n1


## Print one specific line by number
# sed is an editor used to process text line by line
# sed prints every line by default
# -n tells it: “only print what I explicitly ask for”

sed -n '6p' Data/references/gencode.v41.primary_assembly.annotation.chr21.gtf
sed -n '7p' Data/references/gencode.v41.primary_assembly.annotation.chr21.gtf

## Number of annotated features
wc -l Data/references/gencode.v41.primary_assembly.annotation.chr21.gtf


# ______________________________________________________________________________
## Prepare Reference for chr21 from FASTA and GFT

# 1. Create a shell script named "01_prepare_reference.sh"

touch 01_prepare_reference.sh    ## Creates empty file
ls
less 01_prepare_reference.sh
q


## Add text (commands) to .sh file
nano 01_prepare_reference.sh

# #!/bin/bash          👉 “Run this is a
#                          script in Bash language, interpret the text based on it.”
#
# # change to directory where the FASTA and GTF files are
# cd Data/references/
#
# # run mkref
# cellranger mkref\                 👉 Use slash at the end of command to tell the
#                                      operating system this is part of the same line
#   --fasta "YOUR CODE HERE" \
#   --genes "YOUR CODE HERE" \
#   --genome "YOUR CODE HERE" \
#   --nthreads "YOUR CODE HERE"   👉 Tells how many CPUs to use for indexing to paralelize
#                                    the work (depends on the number of CPUs in your machine)

## Run script
bash 01_prepare_reference.sh


## Outputs in Data/references/GRCh38_chr21_index
ls Data/references/GRCh38_chr21_index/

## Look at files and their size
ls -lh Data/references/GRCh38_chr21_index/fasta



# ______________________________________________________________________________
#                    2.  Read alignment and counting
# ______________________________________________________________________________

## Check index folder out from reference indexing
ls -lh Data/references/GRCh38_chr21_index/fasta

## We'll need fastq files to align their reads
ls Data/fastq/
## Dont read gz files! these are binary

## To see them
gunzip Data/fastq/SRR9264343_S1_L001_I1_001.fastq.gz
gunzip Data/fastq/SRR9264343_S1_L001_R1_001.fastq.gz
gunzip Data/fastq/SRR9264343_S1_L001_R2_001.fastq.gz

## Same number of lines each file (4 * number of reads)
wc -l Data/fastq/SRR9264343_S1_L001_I1_001.fastq
wc -l Data/fastq/SRR9264343_S1_L001_R1_001.fastq
wc -l Data/fastq/SRR9264343_S1_L001_R2_001.fastq

head Data/fastq/SRR9264343_S1_L001_I1_001.fastq
head Data/fastq/SRR9264343_S1_L001_R1_001.fastq
head Data/fastq/SRR9264343_S1_L001_R1_001.fastq

nano 02_cellranger_count.sh

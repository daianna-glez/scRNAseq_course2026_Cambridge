
################################################################################
#                       4. Quality Control and Filtering
################################################################################

#-------------------------------------------------------------------------------
#                4.1 Libraries, parallelisation and seed setting
#-------------------------------------------------------------------------------
# Load the libraries we will need for this practical
library(Seurat)
library(sctransform)
library(tidyverse)

# Read 10x data matrices for 1 sample only
etv6_runx1_1_data <- Read10X(
    data.dir = "Data/CellRanger_Outputs/SRR9264343/outs/filtered_feature_bc_matrix/"
)

dim(etv6_runx1_1_data)
# [1] 38606  3153
class(etv6_runx1_1_data)

# Create Seurat object from the loaded matrices
etv6_runx1_1 <- CreateSeuratObject(
    counts = etv6_runx1_1_data,
    project = "ETV6_RUNX1_1"
)

# Look at the structure of the Seurat object
etv6_runx1_1

# The metadata can be accessed using [[]]
head(etv6_runx1_1[[]])
head(etv6_runx1_1@meta.data)
## Per cell, total UMIs and genes
#                      orig.ident nCount_RNA nFeature_RNA
# AAACCTGAGACTTTCG-1 ETV6_RUNX1_1       8354         2935
# AAACCTGGTCTTCAAG-1 ETV6_RUNX1_1      14974         4341
# AAACCTGGTGCAACTT-1 ETV6_RUNX1_1       1566          931
# AAACCTGGTGTTGAGG-1 ETV6_RUNX1_1      10468         3636
# AAACCTGTCCCAAGTA-1 ETV6_RUNX1_1      10437         3340
# AAACCTGTCGAATGCT-1 ETV6_RUNX1_1       2453         1392

## Assays tab
etv6_runx1_1@assays
# $RNA
# Assay (v5) data with 38606 features for 3153 cells
# First 10 features:
#     DDX11L2, MIR1302-2HG, FAM138A, ENSG00000290826, OR4F5, ENSG00000238009,
# ENSG00000239945, ENSG00000239906, ENSG00000241860, ENSG00000241599
# Layers:
#     counts

# Set the default assay to RNA (raw counts)
DefaultAssay(etv6_runx1_1) <- "RNA"

# Check the active assay was set correctly
DefaultAssay(etv6_runx1_1)
etv6_runx1_1@active.assay

# Pull out the raw counts matrix for the RNA assay for exploratory analysis and QC
raw_counts <- etv6_runx1_1[["RNA"]]$counts
raw_counts <- etv6_runx1_1@assays$RNA$counts

#-------------------------------------------------------------------------------
#                     4.2 Properties of Single cell data
#-------------------------------------------------------------------------------

# Calculate the number of genes detected per cell
genes_per_cell <- colSums(raw_counts > 0)

# Plot the distribution of genes detected per cell
plot(density(genes_per_cell),
     main = "",
     xlab = "Genes per cell")

# Calculate the total UMIs per gene
total_umis_per_gene <- rowSums(raw_counts)

# Calculate the proportion of cells expressing each gene
proportion_cells_expressing <- rowMeans(raw_counts > 0)

# Plot the relationship
plot(
    x = total_umis_per_gene,
    y = proportion_cells_expressing,
    log = "x",
    xlab = "Total UMIs per gene",
    ylab = "Proportion of cells expressing the gene"
)


# Calculate the relative expression of each gene in each cell
rel_expression <- t(t(raw_counts) / colSums(raw_counts)) * 100

# Get the 20 most expressed genes
most_expressed <- sort(rowSums(rel_expression), decreasing = TRUE)[20:1]

# Plot the distribution of relative expression for the top 20 genes
plot_data <- as.matrix(t(rel_expression[names(most_expressed), ]))

boxplot(plot_data, cex = 0.1, las = 1,
        xlab = "% total count per cell",
        horizontal = TRUE)


#-------------------------------------------------------------------------------
#                          4.3 Quality Control
#-------------------------------------------------------------------------------

## 4.3.1 Detected Genes
# How many genes are detected in the whole dataset?
table(rowSums(raw_counts) > 0)

## 4.3.2 Mitochondrial Gene Content
# Identify mitochondrial genes
mito_genes <- grep(pattern = "^MT-",
                   x = rownames(raw_counts), value = TRUE)
length(mito_genes)

# Calculate the percentage of UMIs mapped to mitochondrial genes and add this as a column to the metadata
etv6_runx1_1[["percent.mt"]] <- PercentageFeatureSet(etv6_runx1_1,
                                                     pattern = "^MT-")

# Check the metadata to see the new column
head(etv6_runx1_1[[]])
#                      orig.ident nCount_RNA nFeature_RNA percent.mt
# AAACCTGAGACTTTCG-1 ETV6_RUNX1_1       8354         2935   3.507302
# AAACCTGGTCTTCAAG-1 ETV6_RUNX1_1      14974         4341   3.806598
# AAACCTGGTGCAACTT-1 ETV6_RUNX1_1       1566          931  27.522350
# AAACCTGGTGTTGAGG-1 ETV6_RUNX1_1      10468         3636   4.088651
# AAACCTGTCCCAAGTA-1 ETV6_RUNX1_1      10437         3340   5.049344
# AAACCTGTCGAATGCT-1 ETV6_RUNX1_1       2453         1392   3.668977


# 4.3.3 Visualize QC Metrics

# Distribution of number of genes detected per cell
VlnPlot(etv6_runx1_1, features = c("nFeature_RNA"), layer = "counts")

# Distribution of total UMI counts per cell
VlnPlot(etv6_runx1_1, features = c("nCount_RNA"), layer = "counts")

# Distribution of percentage of mitochondrial genes per cell
VlnPlot(etv6_runx1_1, features = c("percent.mt"), layer = "counts")

# Distribution of all three QC metrics together
VlnPlot(etv6_runx1_1,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, layer = "counts")

#-------------------------------------------------------------------------------
#                         4.4 Loading Multiple Samples
#-------------------------------------------------------------------------------

# Read our sample information sheet
sampleinfo <- read_tsv("Data/sample_sheet.tsv")

# Create a named vector of directories for each sample
# We work with a subset of samples for demonstration purposes
samples <- sampleinfo$Sample[c(1, 5, 7, 9)]
list_of_files <- file.path("Data/CellRanger_Outputs/",
                           samples,
                           "/outs/filtered_feature_bc_matrix")
names(list_of_files) <- sampleinfo$SampleName[c(1, 5, 7, 9)]

# Read in the cellranger matrices for all samples
expression_matrix <- Read10X(data.dir = list_of_files)
dim(expression_matrix)
# [1] 38606 14826
expression_matrix[1:4, 1:4]

# Create a new seurat object with the combined expression matrix
multi_seurat_object <- CreateSeuratObject(counts = expression_matrix)
multi_seurat_object

# Pull out the metadata and add sample group and sample name information
temp_metadata <- multi_seurat_object[[]] %>%
    # Add the cell barcodes as a column so we can pull them back in at the end
    rownames_to_column("Cell") %>%
    # Extract the sample group information from the origin identifier
    # The pattern "-.*" matches the first dash and everything after it
    mutate(SampleGroup = str_remove(orig.ident, "-.*")) %>%
    # Add the sample name as a column
    mutate(SampleName = orig.ident) %>%
    # Put the cell barcodes back as rownames
    column_to_rownames("Cell")

# Add the modified metadata back to the Seurat object
multi_seurat_object[[]] <- temp_metadata

## 4.4.2 Remove undetected genes
# Get the raw counts to work with
multi_raw_counts <- multi_seurat_object[["RNA"]]$counts

# How many genes are detected in the whole dataset?
table(rowSums(multi_raw_counts) > 0)

# Filter out genes that are not detected in any cell
filtered_multi_seurat_object <- subset(
    multi_seurat_object,
    features = rownames(multi_seurat_object)[rowSums(multi_raw_counts) > 0]
)

# Check how many genes are detected after filtering
filtered_multi_seurat_object


## 4.4.3 Visualize QC Metrics for Multiple Samples

# Add the percentage of mitochondrial genes to the metadata
filtered_multi_seurat_object[["percent.mt"]] <-
    PercentageFeatureSet(filtered_multi_seurat_object, pattern = "^MT-")

# Plot the distribution of QC metrics for the merged samples
VlnPlot(filtered_multi_seurat_object,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        layer = "counts"
)

## 4.4.4 Filtering out low quality droplets






################################################################################
#                    5. Normalisation and Feature selection
################################################################################

# Load the libraries we will need for this practical
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(tidyverse)
library(patchwork)
use("SparseArray", c("rowMeans", "rowVars", "colMeans", "colVars"))

# set default ggplot theme
theme_set(theme_classic())


# Read preprocessed and filtered data object
seurat_object <- readRDS("RObjects/Filtered.500.rds")
seurat_object

# An object of class Seurat
# 29786 features across 5500 samples within 1 assay
# Active assay: RNA (29786 features, 0 variable features)
# 1 layer present: counts


# Subset seurat object to include only PBMMC-1 sample for demonstration purposes
seurat_pbmmc1 <- subset(seurat_object,
                        subset = SampleName == "PBMMC-1")
seurat_pbmmc1
# An object of class Seurat
# 22644 features across 500 samples within 1 assay
# Active assay: RNA (22644 features, 0 variable features) <- already filtered!!!!
# 1 layer present: counts

# Remove genes that are not expressed in any of the 500 cells in this sample

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#                                !!! Warning !!!
#_______________________________________________________________________________
# This would change if done across all samples or a different subset
# If done with all samples together, you'll be less stric as less genes won't be 0 expressed
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

pbmmc1_expressed_genes <- rowSums(seurat_pbmmc1[["RNA"]]$counts) > 0
# pbmmc1_expressed_genes <- rowSums(GetAssayData(seurat_pbmmc1, assay = "RNA", layer = "counts")) > 0

seurat_pbmmc1 <- seurat_pbmmc1[pbmmc1_expressed_genes, ]
dim(seurat_pbmmc1)
# [1] 22644   500


## Why we normalize:

## Reason 1: Differences in total UMIs between cells
#            (due to biology, doublets or empty droplets, or technical effect)

# Plot the total UMI counts across cells
# Get the raw counts data
seurat_pbmmc1[["RNA"]]$counts %>%
    # calculate the total UMI counts for each cell (column sums)
    colSums() %>%
    # convert the named vector to a data frame with cell names and UMI counts
    enframe(name = "cell", value = "total_umis") %>%
    # reorder the cells by their total UMI counts for plotting
    dplyr::mutate(Cell = fct_reorder(cell, total_umis)) %>%
    # make a bar plot of the total UMI counts per cell
    ggplot(aes(x = Cell, y = total_umis)) +
    geom_col() +
    labs(x = "Cell",
         y = "Total cell UMI counts",
         title = "PBMMC-1: Before Normalization") +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(color = "firebrick")
    )


## Reason 2: Mean-variance relationship across genes
# Get the raw counts data from the object
raw_cts <- seurat_pbmmc1[["RNA"]]$counts

# Calculate summary statistics for each gene
gene_raw_stats <- tibble(gene = rownames(raw_cts),
                         mean = rowMeans(raw_cts),
                         variance = rowVars(raw_cts))

# Plot the mean-variance relationship for the raw counts data
ggplot(gene_raw_stats, aes(x = mean, y = variance)) +
    geom_point(alpha = 0.5) +
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(intercept = 0, slope = 1,
                linewidth = 1, colour = "red") +
    labs(x = "Mean expression (log scale)",
         y = "Variance (log scale)",
         title = "PBMMC-1: Raw Counts - Mean vs Variance")


# ______________________________________________________________________________
#                       1. Shifted-log transformation
# ______________________________________________________________________________

# The default normalisation method (shifted-log transformation)
# 1. counts of each gene are divided by total cell UMIs
# 2. multiplied by a scale factor (default 10,000)
# 3. log-transformed with a pseudocount of 1

seurat_pbmmc1 <- NormalizeData(seurat_pbmmc1)

# A new layer named `data` is added to the RNA assay
seurat_pbmmc1

# Get the log-normalised data from the object
lognorm_cts <- seurat_pbmmc1[["RNA"]]$data
gene_lognorm_stats <- tibble(gene = rownames(lognorm_cts),
                             mean = rowMeans(lognorm_cts),
                             variance = rowVars(lognorm_cts))

# Plot the mean-variance relationship for the log-normalised data
# improved for highly expressed genes
# but still heteroscedastic for lowly expressed genes
ggplot(gene_lognorm_stats, aes(x = mean, y = variance)) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1,
                linewidth = 1, colour = "red") +
    labs(x = "Mean expression",
         y = "Variance",
         title = "PBMMC-1: Log Normalized - Mean vs Variance")


# ______________________________________________________________________________
#                             2. Identification of HVGs
# ______________________________________________________________________________
## Identify HVGs using Seurat's variance stabilising transform (vst) method

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#                             !!! Warning !!!
## This is not normalizing, it assumes normalization has already been done.
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# Extract the 3000 most variable genes
# using the variance stabilising transformation method
seurat_pbmmc1 <- FindVariableFeatures(
    seurat_pbmmc1,
    selection.method = "vst",
    nfeatures = 3000
)

VariableFeaturePlot(seurat_pbmmc1)

# Get the variable genes selected by the vst method
hvgs_vst <- VariableFeatures(seurat_pbmmc1)
head(hvgs_vst)


# ______________________________________________________________________________
#                       3. Sctransform normalisation
# ______________________________________________________________________________

# Run sctransform normalisation
# regressing out the percentage of mitochondrial gene expression
seurat_pbmmc1 <- SCTransform(
    seurat_pbmmc1,
    vars.to.regress = "percent.mt",
    verbose = FALSE
)

# A new assay named `SCT` is added to the Seurat object
# 'counts' layer contains the corrected UMI counts (no residuals, “observed counts adjusted to a common sequencing depth”)
# 'data' layer contains the log-transformed counts data, with a pseudocount of 1 added
# 'scale.data' layer contains the variance-stabilized Pearson residuals of the model fit (“observed − expected” = residuals)
seurat_pbmmc1

# An object of class Seurat
# 37494 features across 500 samples within 2 assays
# Active assay: SCT (14850 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA


# Get the variable genes selected by sctransform
# which are ranked by residual variance
hvgs_sct <- VariableFeatures(seurat_pbmmc1)
head(hvgs_sct)


# Compare the variable genes selected
# by vst on the raw data and from sctransform residual variance
length(intersect(hvgs_vst, hvgs_sct))

# Compare the mean-variance relationship for the sctransform-normalised data
# Get summary statistics for the Pearson residuals from the sctransform normalisation
vst_cts <- seurat_pbmmc1[["SCT"]]$scale.data
gene_vst_stats <- tibble(gene = rownames(vst_cts),
                         mean = rowMeans(vst_cts),
                         variance = rowVars(vst_cts))

# Plot the mean-variance relationship for the sctransform-normalised data
ggplot(gene_vst_stats, aes(x = mean, y = variance)) +
    geom_point(alpha = 0.5) +
    labs(x = "Mean expression",
         y = "Variance",
         title = "PBMMC-1: sctransform Normalized - Mean vs Variance")

# Compare the total UMI counts per cell for the raw counts and sctransform-normalised data

# get the total UMIs for raw UMIs
raw_cts_total <- seurat_pbmmc1[["RNA"]]$counts %>%
    colSums() %>%
    enframe(name = "cell", value = "total_counts")

# get the total UMIs for sctransform-normalised counts
sct_cts_total <- seurat_pbmmc1[["SCT"]]$counts %>%
    colSums() %>%
    enframe(name = "cell", value = "total_counts")

# Plot the total UMI counts per cell for the PBMMC-1 sample before normalisation
p_raw <- raw_cts_total %>%
    # reorder the cells by their total UMI counts for plotting
    mutate(Cell = fct_reorder(cell, total_counts)) %>%
    # make a bar plot of the total UMI counts per cell
    ggplot(aes(x = Cell, y = total_counts)) +
    geom_col() +
    labs(x = "Cell",
         y = "Total cell UMI counts",
         title = "Before Normalization") +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(color = "firebrick")
    )

p_sct <- sct_cts_total %>%
    # reorder the cells by their total UMI counts for plotting
    mutate(Cell = fct_reorder(cell, total_counts)) %>%
    # make a bar plot of the total UMI counts per cell
    ggplot(aes(x = Cell, y = total_counts)) +
    geom_col() +
    labs(x = "Cell",
         y = "Total cell UMI counts",
         title = "SCTransform") +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(color = "firebrick")
    )

# combine the two plots side by side
p_raw + p_sct


################################################################################
#                 7. Dataset Integration and Batch Correction
################################################################################

# Load the libraries we will need for this practical
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(tidyverse)
library(patchwork)


## First we must examine our uncorrected data and decide if we think a correction is needed
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#             7.1 Examining our uncorrected data (clustering)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# 7.1.1 Re-process the data by sample group

# Load the data: 500 cells per sample across 4 groups
seurat_object <- readRDS("RObjects/Filtered.500.rds")
seurat_object

# Split the data by sample group
seurat_object[["RNA"]] <- split(seurat_object[["RNA"]],
                                f = seurat_object$SampleGroup) # f specifies the variable that we
#                                                                want to split by

## Adds new layers, one for each group
seurat_object
# An object of class Seurat
# 29786 features across 5500 samples within 1 assay
# Active assay: RNA (29786 features, 0 variable features)
# 4 layers present: counts.ETV6RUNX1, counts.HHD, counts.PRET, counts.PBMMC

# Normalize and stabilice variance in each sample group separately
seurat_object <- SCTransform(
    seurat_object,
    assay = "RNA",
    vars.to.regress = "percent.mt",
    verbose = FALSE
)

## Note active assay has changed to SCT with 3 layers:
#  - counts: corrected UMI counts (no residuals, “observed counts adjusted to a common sequencing depth”)
#  - data: log-transformed counts data, with a pseudocount of 1 added
#  - scale.data: variance-stabilized Pearson residuals of the model fit (“observed − expected” = residuals)

seurat_object
# An object of class Seurat
# 50798 features across 5500 samples within 2 assays
# Active assay: SCT (21012 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA

## Run PCA on selected HVGs
seurat_object <- RunPCA(seurat_object,
                        features = VariableFeatures(object = seurat_object))

## Added slot for reduction
seurat_object
# An object of class Seurat
# 50798 features across 5500 samples within 2 assays
# Active assay: SCT (21012 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 1 dimensional reduction calculated: pca


## Project data on low dimensional space based on PC1-15
seurat_object <- RunUMAP(seurat_object,
                         reduction = "pca",
                         dims = 1:15)

## Additional reduction slot for UMAP
seurat_object
# An object of class Seurat
# 50798 features across 5500 samples within 2 assays
# Active assay: SCT (21012 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

# Visualise the uncorrected data: see how cells group by sample
uncorrected_plot <- DimPlot(seurat_object,
                            reduction = "umap",
                            group.by = "SampleName") +
    ggtitle("Uncorrected data")


# Run clustering using PCs to see how cells group in the uncorrected data:

## Finds similar cells to every cell
seurat_object <- FindNeighbors(seurat_object, reduction = "pca",
                               dims = 1:15)

## Define clusters based on cell neighbours
seurat_object <- FindClusters(seurat_object,
                              cluster.name = "uncorrected_clusters")

## Plot cluster lables on UMAP
uncorrected_clusters_plot <- DimPlot(seurat_object,
                                     reduction = "umap",
                                     group.by ="uncorrected_clusters") +
    ggtitle("Uncorrected data clusters")

uncorrected_plot + uncorrected_clusters_plot


## Note the large overlap between cluster labels and sample donor!
# Make a table of the clusters and samples
table(seurat_object$uncorrected_clusters, seurat_object$orig.ident)

## Some clusters are almost entirely made of cells from 1 donor,
## and that donor has cells on that cluster only = almost entirely correlated


## But differences across conditions of the samples are expected
## (disease vs control, diff tissue, etc)!

## There could be some true cell types specific to one batch/sample
## (look expression of markers and add replicates !!!!)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#             7.2 Integrate Data and Batch correction - Harmony
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Seurat allows us to integrate our data using its IntegrateLayers()
# This function takes the layers of the RNA assay that we split into and integrates them
# using a specified batch correction method.


# Run Harmony integration:
harmony_object <- IntegrateLayers(object = seurat_object,
                                  method = HarmonyIntegration,
                                  orig.reduction = "pca",
                                  new.reduction = "harmony")
harmony_object
# An object of class Seurat
# 50798 features across 5500 samples within 2 assays
# Active assay: SCT (21012 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, umap, harmony


## Visualize clusters in UMAP using integrated data (not original PCA embedding)
harmony_object <- RunUMAP(harmony_object,
                          reduction = "harmony",
                          dims = 1:15)

# Run clustering
harmony_object <- FindNeighbors(harmony_object, reduction = "harmony",
                                dims = 1:15)
harmony_object <- FindClusters(harmony_object,
                               cluster.name = "harmony_clusters")
# make the plots
harmony_plot <- DimPlot(harmony_object,
                        reduction = "umap",
                        group.by = "SampleName") +
    ggtitle("Harmony integrated data")
harmony_clusters_plot <- DimPlot(harmony_object,
                                 reduction = "umap",
                                 group.by ="harmony_clusters") +
    ggtitle("Harmony data clusters")

harmony_plot + harmony_clusters_plot


## Pars:
#  - theta: diversity clustering penalty parameter,
#           controls batch diversity within clusters, higher values -> more diverse clusters
#  - lambda: controls the strenght of correction,
#            higher values leading to stronger corrections
#  - sigma: controls the bandwidth of the Gaussian kernel used in the correction

# Run Harmony integration with adjusted theta
harmony_object_theta <- IntegrateLayers(object = seurat_object,
                                        method = HarmonyIntegration,
                                        orig.reduction = "pca",
                                        new.reduction = "harmony_theta_01",
                                        theta = 0.1)
## UMAP on harmonized data
harmony_object_theta <- RunUMAP(harmony_object_theta,
                                reduction = "harmony_theta_01",
                                dims = 1:15)

## Clustering
harmony_object_theta <- FindNeighbors(harmony_object_theta,
                                      reduction = "harmony_theta_01",
                                      dims = 1:15)

harmony_object_theta <- FindClusters(harmony_object_theta,
                                     cluster.name = "harmony_theta_01_clusters")

## Plot
harmony_object_theta_plot <- DimPlot(harmony_object_theta,
                                     reduction = "umap",
                                     group.by = "SampleName") +
    ggtitle("Harmony integrated data with theta = 0.1")

harmony_object_theta_clusters_plot <- DimPlot(harmony_object_theta,
                                              reduction = "umap",
                                              group.by ="harmony_theta_01_clusters") +
    ggtitle("Harmony integrated data clusters with theta = 0.1")

harmony_object_theta_plot + harmony_object_theta_clusters_plot

## Note mixture of batches with larger theta as we penalise more by overrepresentation
# of batch in a cluster so prob of a cell with batch b in cluster k is decreased


###################
## Exercise 1: Run integration at the sample level and compare results
###################

## 1. Starting with our filtered data, split the object by SampleName.
# Read in our filtered data
seurat_object_sample <- readRDS("RObjects/Filtered.500.rds")

# Split the data by sample name
seurat_object_sample[["RNA"]] <- split(seurat_object_sample[["RNA"]],
                                       f = seurat_object_sample$SampleName)
# Re-process, Seurat will treat each sample separately
seurat_object_sample <- SCTransform(
    seurat_object_sample,
    assay = "RNA",
    vars.to.regress = "percent.mt",
    verbose = FALSE
)
seurat_object_sample <- RunPCA(seurat_object_sample,
                               features = VariableFeatures(object = seurat_object_sample))








# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                           7.3 Correction Diagnostics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Uncorrected clusters
head(seurat_object$uncorrected_clusters)
# ETV6RUNX1-1_AAAGATGCAAGCTGAG-1 ETV6RUNX1-1_AAAGATGTCTCCCTGA-1 ETV6RUNX1-1_AAAGCAAAGTCCAGGA-1
# 2                              1                              2
# ETV6RUNX1-1_AAAGCAATCTGCTGCT-1 ETV6RUNX1-1_AAATGCCAGAACAATC-1 ETV6RUNX1-1_AACCATGAGAGCAATT-1
# 2                              2                              2
# Levels: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19

## Batches
seurat_object$SampleName %>% head()

# Make a barplot of the clusters and samples
data.frame(Cluster = seurat_object$uncorrected_clusters,
           Sample = seurat_object$SampleName) %>%
    ggplot(aes(x = Cluster)) +
    geom_bar(aes(fill = Sample), position = "fill") +
    labs(title = "Uncorrected data") +
    scale_y_continuous(labels = scales::percent)


# Make a barplot of the clusters and samples after batch correction
data.frame(Cluster = harmony_object_theta$harmony_theta_01_clusters,
           Sample = harmony_object_theta$SampleName) %>%
    ggplot(aes(x = Cluster)) +
    geom_bar(aes(fill = Sample), position = "fill") +
    labs(title = "Harmony Corrected data") +
    scale_y_continuous(labels = scales::percent)

## Note:
# 1. fewer clusters
# 2. more representation of donors per cluster
# 3. some clusters are largely made of donors from same condition / control (expected)


## Has the correction over-worked? Plot cell type marker genes to assess cell type identification
# CD79A (B cells)
# CST3 (monocytes)
# CD3D (T cells)
# HBA1 (erythrocytes)

## Not corrected data: cell types spread (because clusters are not cell types but donors)
FeaturePlot(seurat_object,
            features = c("CD79A", "CST3", "CD3D", "HBA1"),
            reduction = "umap")

## Cell types are there! So biological signal has not been removed
FeaturePlot(harmony_object_theta,
            features = c("CD79A", "CST3", "CD3D", "HBA1"),
            reduction = "umap")


## Adjusted Rand Index: prop of cell pairs with same cluster assignment status
## (relative to each other) before and after correction

library(bluster)
# Calculate the adjusted Rand index
# Input is a vector of the clusters
ari.SG <- pairwiseRand(seurat_object$uncorrected_clusters,
                       harmony_object$harmony_clusters,
                       mode = "index")

ari.SG
# [1] 0.3776693 <- within-batch heterogeneity poorly preserved
#                  (2 cells from same donor now in diff clusters)

## Correcting at sample level is worse
ari.SN <- pairwiseRand(seurat_object_sample$uncorrected_clusters,
                       harmony_object_sample$harmony_sample_clusters,
                       mode = "index")

ari.SN
# [1] 0.3106236

ari.T <- pairwiseRand(seurat_object$uncorrected_clusters,
                      harmony_object_theta$harmony_theta_01_clusters,
                      mode = "index")

ari.T
# [1] 0.5711085  <- improved but still not ideal; would try other theta params


## Rejoin RNA data for downstream analysis (on RNA not SCT nor corrected data)
# Rejoin the data
harmony_object_theta <- JoinLayers(harmony_object_theta, assay = "RNA")








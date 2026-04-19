
################################################################################
#                       6. Dimensionality Reduction
################################################################################

# Load the libraries we will need for this practical
library(Seurat)
library(tidyverse)

# set default ggplot theme
theme_set(theme_classic())

# load the preprocessed Seurat object with 500 cells per sample
seurat_object <- readRDS("RObjects/SCT.500.rds")
seurat_object

#-------------------------------------------------------------------------------
#                                    6.1 PCA
#-------------------------------------------------------------------------------

# Run PCA
seurat_object <- RunPCA(seurat_object,
                        features = VariableFeatures(seurat_object))
seurat_object
# An object of class Seurat
# 53387 features across 5500 samples within 2 assays
# Active assay: SCT (23601 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 1 dimensional reduction calculated: pca

# check the new reductions slot
Reductions(seurat_object)
# [1] "pca"

# check the variance explained by each PC
Stdev(seurat_object, reduction = "pca")

# plot the variance explained by each PC
ElbowPlot(seurat_object, ndims = 50)

# plot the PCA
DimPlot(seurat_object,
        reduction = "pca")

# plot PC2 and PC3
DimPlot(seurat_object,
        reduction = "pca",
        dims = c(2, 3))

# plot PC1 and PC2 split by sample group
DimPlot(seurat_object,
        reduction = "pca",
        split.by = "SampleGroup")


#-------------------------------------------------------------------------------
#                                6.2 t-SNE
#-------------------------------------------------------------------------------

# run t-SNE using default options
seurat_object <- RunTSNE(seurat_object,
                         reduction = "pca")

# confirm a new reduction was added to the object
Reductions(seurat_object)

# plot the t-SNE
DimPlot(seurat_object,
        reduction = "tsne")

########################
## Exercise 1:
########################
seurat_object <- RunTSNE(seurat_object,
                         reduction = "pca",
                         dims = 1:10,
                         perplexity = 30,
                         seed = 123,
                         reduction.name = "TSNE_perplex30",
                         reduction.key = "TSNE30_")
DimPlot(seurat_object,
        reduction = "TSNE_perplex30")

## Part A)
# re-run the t-SNE with a different seed
seurat_object <- RunTSNE(seurat_object,
                         reduction = "pca",
                         dims = 1:10,
                         perplexity = 30,
                         seed = 321,
                         reduction.name = "TSNE_perplex30_seed321",
                         reduction.key = "TSNE30seed321_")

# plot the new t-SNE
DimPlot(seurat_object,
        reduction = "TSNE_perplex30_seed321")


## Part B)
# plot the t-SNE faceted by SampleName
DimPlot(seurat_object,
        reduction = "TSNE_perplex30_seed321",
        split.by = "SampleName",
        ncol = 4)

## Part C)
# re-run the t-SNE with perplexity 5
seurat_object <- RunTSNE(seurat_object,
                         reduction = "pca",
                         dims = 1:10,
                         perplexity = 5,
                         seed = 321,
                         reduction.name = "TSNE_perplex5",
                         reduction.key = "TSNE5_")
# re-run the t-SNE with perplexity 500
seurat_object <- RunTSNE(seurat_object,
                         reduction = "pca",
                         dims = 1:10,
                         perplexity = 500,
                         seed = 321,
                         reduction.name = "TSNE_perplex500",
                         reduction.key = "TSNE500_")

# visualise all projections using DimPlot
DimPlot(seurat_object, reduction = "TSNE_perplex5")
DimPlot(seurat_object, reduction = "TSNE_perplex500")


## Part D)
# colour tsne by CD79A expression
FeaturePlot(seurat_object,
            features = "CD79A",
            reduction = "TSNE_perplex30")


#-------------------------------------------------------------------------------
#                                6.3 UMAP
#-------------------------------------------------------------------------------
# run UMAP using the first 10 PCs
seurat_object <- RunUMAP(seurat_object,
                         reduction = "pca",
                         dims = 1:10)

# confirm a new reduction was added to the object
Reductions(seurat_object)
## [1] "pca"                    "tsne"                   "TSNE_perplex30"
##      "TSNE_perplex30_seed321" "TSNE_perplex5"          "TSNE_perplex500"        "umap"

# visualise the UMAP
DimPlot(seurat_object,
        reduction = "umap")

########################
## Exercise 2:
########################

## Part A)
# run the UMAP with 30 neighbours
seurat_object <- RunUMAP(seurat_object,
                         reduction = "pca",
                         dims = 1:10,
                         seed.use = 123,
                         n.neighbors = 30,
                         reduction.name = "UMAP_neighbors30",
                         reduction.key = "UMAPneighbors30_")

## Part B)
# visualise the resulting UMAP projection
DimPlot(seurat_object,
        reduction = "UMAP_neighbors30")

## Part C)
# run the UMAP with 5 neighbours
seurat_object <- RunUMAP(seurat_object,
                         reduction = "pca",
                         dims = 1:10,
                         seed.use = 123,
                         n.neighbors = 5,
                         reduction.name = "UMAP_neighbors5",
                         reduction.key = "UMAPneighbors5_")
# run the UMAP with 500 neighbours
seurat_object <- RunUMAP(seurat_object,
                         reduction = "pca",
                         dims = 1:10,
                         seed.use = 123,
                         n.neighbors = 500,
                         reduction.name = "UMAP_neighbors500",
                         reduction.key = "UMAPneighbors500_")

# visualise all projections using DimPlot
DimPlot(seurat_object, reduction = "UMAP_neighbors5")
DimPlot(seurat_object, reduction = "UMAP_neighbors500")



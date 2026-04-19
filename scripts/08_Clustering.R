

################################################################################
#                              8. Clustering
################################################################################

# Load libraries
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(cluster) # for silhouette function
library(tidyverse)
library(patchwork)

# set default ggplot theme
theme_set(theme_classic())

# Load the integrated data (batch corrected)
seurat_object <- readRDS("RObjects/DI.500.rds")
seurat_object

# An object of class Seurat
# 50798 features across 5500 samples within 2 assays
# Active assay: SCT (21012 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, harmony, umap


# Confirm the default assay is set to SCT
DefaultAssay(seurat_object)
# [1] "SCT"

# construct the neighbour graph, by default it constructs both:
#   - nearest neighbour (NN) graph
#   - shared nearest neighbour (SNN) graph
seurat_object <- FindNeighbors(seurat_object,
                               reduction = "harmony",
                               k.param = 20, # k = 20 is the default
                               dims = 1:15)

# Check the graphs slot
Graphs(seurat_object)
# [1] "SCT_nn"  "SCT_snn"


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                              8.1 Louvain method
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# identify clusters in the nearest neighbour graph (SNN used by default)
seurat_object <- FindClusters(seurat_object,
                              resolution = 0.8, # default is 0.8
                              algorithm = 1, # default is 1 (louvain)
                              cluster.name = "Louvain_k20_res0.8")

# confirm the new column in the meta.data slot
head(seurat_object[[]])

# visualise the clusters on a UMAP plot
DimPlot(seurat_object,
        reduction = "umap",
        group.by = "Louvain_k20_res0.8")


###############
# Exercise 1:
###############

# Cluster with resolution 0.4
seurat_object <- FindClusters(seurat_object,
                              resolution = 0.4,
                              algorithm = 1,
                              cluster.name = "Louvain_k20_res0.4")

# Cluster with resolution 1.6
seurat_object <- FindClusters(seurat_object,
                              resolution = 1.6,
                              algorithm = 1,
                              cluster.name = "Louvain_k20_res1.6")

# Visualise the clusters on a UMAP plot
DimPlot(seurat_object, reduction = "umap",
        group.by = c("Louvain_k20_res0.4",
                     "Louvain_k20_res0.8",
                     "Louvain_k20_res1.6"))

# -----------------------------------------
# Find neighbours with k = 10
seurat_object <- FindNeighbors(seurat_object,
                               reduction = "harmony",
                               k.param = 10,
                               dims = 1:15)
# Cluster with resolution 0.8
seurat_object <- FindClusters(seurat_object,
                              resolution = 0.8,
                              algorithm = 1,
                              cluster.name = "Louvain_k10_res0.8")
# Find neighbours with k = 40
seurat_object <- FindNeighbors(seurat_object,
                               reduction = "harmony",
                               k.param = 40,
                               dims = 1:15)
# Cluster with resolution 0.8
seurat_object <- FindClusters(seurat_object,
                              resolution = 0.8,
                              algorithm = 1,
                              cluster.name = "Louvain_k40_res0.8")

# Visualise the clusters on a UMAP plot
DimPlot(seurat_object, reduction = "umap",
        group.by = c("Louvain_k10_res0.8",
                     "Louvain_k20_res0.8",
                     "Louvain_k40_res0.8"))



#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                              8.2 Leiden method
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

###############
# Exercise 2:
###############

# Find neighbours with k = 20
seurat_object <- FindNeighbors(seurat_object,
                               reduction = "harmony",
                               k.param = 20,
                               dims = 1:15)

# Cluster with resolution 0.8 using the Leiden method
seurat_object <- FindClusters(seurat_object,
                              resolution = 0.8,
                              algorithm = 4,
                              cluster.name = "Leiden_k20_res0.8",
                              random.seed = 123) # set seed (warning without)

# Visualise the clusters on a UMAP plot
DimPlot(seurat_object, reduction = "umap",
        group.by = c("Louvain_k20_res0.8", "Leiden_k20_res0.8"))


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                        8.3 Assessment of clustering
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#### Silhouette width for Louvain clustering ####

# get Harmony data as we used for clustering
harmony_data <- Embeddings(seurat_object, reduction = "harmony")[, 1:15]
dim(harmony_data)
# [1] 5500   15

# calculate euclidean distance matrix for the harmony data
harmony_dist <- dist(harmony_data)
length(harmony_dist) ## upper triangle of matrix
# [1] 15122250

# get our choice of cluster assignments
# because it's a factor, we convert it to integer
clusters <- seurat_object$Louvain_k20_res0.8
clusters <- as.integer(as.character(clusters))

# calculate silhouette width for each cell
silhouette_width <- silhouette(clusters,
                               harmony_dist) %>%
    # convert to data frame for plotting
    as_tibble()

# look at the top 15 rows
head(silhouette_width, n = 15)
# A tibble: 15 × 3
#   cluster neighbor sil_width
#     <dbl>    <dbl>     <dbl>
# 1       7       15   0.312
# 2       6        7   0.0810
# 3      15        7   0.113
# 4       4       15   0.0783
# 5       4        3   0.145
# 6       4        3   0.296
# 7       4        3   0.336
# 8       4        3   0.243
# 9       7        6   0.306
# 10      13       14   0.517
# 11       4        2  -0.0152
# 12       4        2   0.252
# 13       4       15  -0.00521
# 14       4        2   0.0300
# 15       4        3   0.169

dim(silhouette_width)
# [1] 5500    3


# silhouette width distribution coloured by closest cluster
p1 <- silhouette_width %>%
    # create variable for closest cluster based on silhouette width
    mutate(closest_cluster = ifelse(sil_width > 0, cluster, neighbor)) %>%
    ggplot(aes(x = factor(cluster),
               y = sil_width,
               colour = factor(closest_cluster))) +
    ggbeeswarm::geom_quasirandom() +
    labs(x = "Cluster",
         y = "Silhouette width",
         colour = "Closest cluster")

# UMAP plot for comparison
p2 <- DimPlot(seurat_object,
              reduction = "umap",
              group.by = "Louvain_k20_res0.8",
              label = TRUE) +
    NoLegend()

# plot them together
p1 + p2



#### Silhouette width for Leiden clustering ####

# get our choice of cluster assignments
clusters_leiden <- seurat_object$Leiden_k20_res0.8
clusters_leiden <- as.integer(as.character(clusters_leiden))

# calculate silhouette width for each cell
silhouette_width_leiden <- silhouette(clusters_leiden,
                                      harmony_dist) %>%
    # convert to data frame for plotting
    as_tibble()

# silhouette width distribution coloured by closest cluster
p1_leiden <- silhouette_width_leiden %>%
    # create variable for closest cluster based on silhouette width
    mutate(closest_cluster = ifelse(sil_width > 0, cluster, neighbor)) %>%
    ggplot(aes(x = factor(cluster),
               y = sil_width,
               colour = factor(closest_cluster))) +
    ggbeeswarm::geom_quasirandom() +
    labs(x = "Cluster",
         y = "Silhouette width",
         colour = "Closest cluster")

# UMAP plot for comparison
p2_leiden <- DimPlot(seurat_object,
                     reduction = "umap",
                     group.by = "Leiden_k20_res0.8")

# plot them together
p1_leiden + p2_leiden


#### Calculate mean silhouette width for different values of k ####

# function to calculate mean silhouette width for a given clustering
calc_mean_sil <- function(clusters, dist_matrix) {
    silhouette_width <- silhouette(clusters, dist_matrix)
    mean(silhouette_width[, "sil_width"])
}

# function to calculate number of clusters for a given clustering
calc_num_clusters <- function(clusters) {
    length(unique(clusters))
}

# list of k-values we're interested in
k_values <- c(5, 10, 15, 20, 25)

# get the seurat metadata for convenience
seurat_meta <- seurat_object[[]]

# calculate mean silhouette width and number of clusters for each k
# the lapply() function loops through each k value
# and applies the code inside the function
silhouette_stats <- lapply(k_values, function(k) {
    # column name for the clustering with this k value
    cluster_col <- paste0("Louvain_k", k, "_res0.8")

    # get the cluster assignments for this k value
    clusters <- seurat_meta[, cluster_col]

    # convert to integer
    clusters <- as.integer(as.character(clusters))

    # create a data frame to store the results
    tibble(k = k,
           num_clusters = calc_num_clusters(clusters),
           mean_silhouette_width = calc_mean_sil(clusters, harmony_dist))
}) %>%
    # bind the list into a single data frame
    bind_rows()

# look at the results
silhouette_stats



# Add our final choice of clustering to the seurat object
# Create SNN graph with k = 17
seurat_object <- FindNeighbors(seurat_object,
                               reduction = "harmony",
                               k.param = 17,
                               dims = 1:15)

# Cluster cells with resolution = 1 using the Leiden method
seurat_object <- FindClusters(seurat_object,
                              resolution = 1,
                              algorithm = 4,
                              random.seed = 123,
                              cluster.name = "Leiden_k17_res1")

# Visualise the clusters on a UMAP plot
DimPlot(seurat_object,
        reduction = "umap",
        group.by = "Leiden_k17_res1")

# Visualise the expression of a marker gene on the UMAP plot and across clusters
fplot <- FeaturePlot(seurat_object,
                     features = "HBA1",
                     reduction = "umap")

# visualise the distribution of expression across clusters
vplot <- VlnPlot(seurat_object, features = "HBA1",
                 group.by = "Leiden_k17_res1",
                 pt.size = 0) +
    NoLegend()

# plot them together
fplot + vplot




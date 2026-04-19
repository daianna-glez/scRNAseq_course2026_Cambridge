

################################################################################
#                    9. Cluster marker gene identification
################################################################################

# load the libraries
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(tidyverse)
library(patchwork)

# load in the data
seurat_object <- readRDS("RObjects/Clustered.500.rds")
seurat_object

head(seurat_object$leiden_cluster, 2)
# ETV6RUNX1-1_AAAGATGCAAGCTGAG-1 ETV6RUNX1-1_AAAGATGTCTCCCTGA-1
#                              8                              7
# Levels: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16

# visualise UMAP with cluster labels
DimPlot(seurat_object,
        reduction = "umap",
        group.by = "leiden_cluster",
        label = TRUE,
        label.size = 5) +
    NoLegend()

# Visualise known monocyte marker gene CST3
FeaturePlot(seurat_object,
            reduction = "umap",
            features = "CST3")

# Prepare the data for marker gene detection
# (variance-stabilized residuals cannot be directly used,
# counts are recomputed based on estimated NB model coeffs, and same lib size across all cells)

seurat_object <- PrepSCTFindMarkers(seurat_object,
                                    assay = "SCT")
seurat_object


# finding markers of a single cluster
# Compare cells in cluster 2 against all other cells
markers_cluster2 <- FindMarkers(seurat_object,
                                ident.1 = 2,
                                ident.2 = NULL,
                                assay = "SCT",
                                group.by = "leiden_cluster")
head(markers_cluster2)
#                    p_val avg_log2FC pct.1 pct.2     p_val_adj
# IFITM1      0.000000e+00   2.877734 0.905 0.199  0.000000e+00
# SKAP1       0.000000e+00   3.145875 0.755 0.092  0.000000e+00
# TC2N        0.000000e+00   3.787626 0.503 0.036  0.000000e+00
# NELL2       0.000000e+00   5.667427 0.423 0.012  0.000000e+00
# PCED1B-AS1 4.632024e-318   3.495943 0.666 0.099 9.732808e-314
# RPS27      1.291867e-313   1.390475 1.000 0.988 2.714471e-309

## Num genes tested
dim(markers_cluster2)
# [1] 13455     5

# filter for genes that:
#  - are expressed in at least 70% of cells in cluster 2
#  - and have a log fold change greater than 1
markers_cluster2_filtered <- markers_cluster2 %>%
    filter(pct.1 > 0.70 & abs(avg_log2FC) > 1)

head(markers_cluster2_filtered)
#                  p_val avg_log2FC pct.1 pct.2     p_val_adj
# IFITM1    0.000000e+00   2.877734 0.905 0.199  0.000000e+00
# SKAP1     0.000000e+00   3.145875 0.755 0.092  0.000000e+00
# RPS27    1.291867e-313   1.390475 1.000 0.988 2.714471e-309
# RPL21    5.810414e-307   1.433676 1.000 0.992 1.220884e-302
# LEPROTL1 4.938806e-304   3.088247 0.747 0.152 1.037742e-299
# RPS29    5.575920e-281   1.269738 1.000 0.978 1.171612e-276


# visualise the expression of one of the marker genes on the UMAP
fplot <- FeaturePlot(seurat_object,
                     features = "SKAP1",
                     label = TRUE,
                     reduction = "umap")

# As a violin plot
vplot <- VlnPlot(seurat_object, features = "SKAP1",
                 group.by = "leiden_cluster",
                 pt.size = 0) +
    NoLegend()

# SKAP1 is a TCR gene
# It also appears in cluster 11 next to it
fplot + vplot + plot_layout(ncol = 2)


# find markers for all clusters automatically
markers_all <- FindAllMarkers(seurat_object,
                              assay = "SCT",
                              group.by = "leiden_cluster")

head(markers_all)
dim(markers_all)


#################
## Exercise 1:
#################

# using FindMarkers()
markers_cluster13 <- FindMarkers(seurat_object,
                                 ident.1 = 13,
                                 ident.2 = NULL,
                                 assay = "SCT",
                                 group.by = "leiden_cluster")

markers_cluster13_filtered <- markers_cluster13 %>%
    filter(pct.1 > 0.70 & abs(avg_log2FC) > 1)

head(markers_cluster13_filtered)


# Heatmap for top 20 marker genes for cluster 13
top_markers_cluster13 <- markers_cluster13_filtered %>%
    slice_max(n = 20, order_by = abs(avg_log2FC)) %>%
    rownames()

top_markers_cluster13
# [1] "VCAN"    "CSTA"    "S100A8"  "LYZ"     "S100A9"  "FCN1"    "MNDA"    "CFD"
# [9] "TYROBP"  "FCER1G"  "S100A11" "CST3"    "LGALS1"  "TYMP"    "IFI30"   "NAMPT"
# [17] "S100A6"  "CCDC200" "S100A4"  "CTSS"

# heatmap of expression for these genes
# we change default colour scale
# and remove colour legend (clusters)
DoHeatmap(seurat_object,
          features = top_markers_cluster13,
          group.by = "leiden_cluster") +
    scale_fill_viridis_c() +
    guides(colour = "none")

# dot plot of expression for these genes
# we change default colour scale
# and rotate x-axis labels for better visibility
DotPlot(seurat_object,
        features = top_markers_cluster13,
        group.by = "leiden_cluster") +
    scale_colour_viridis_c() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1))


## Manual annotation
# loop through all marker genes
# and extract top-ranked gene names for each cluster
top_markers_all <- lapply(unique(markers_all$cluster), function(x) {
    markers_all %>%
        filter(cluster == x) %>%
        filter(pct.1 > 0.70, abs(avg_log2FC) > 1) %>%
        slice_max(n = 15, order_by = abs(avg_log2FC)) %>%
        pull(gene)
})

# we get several known markers of immune cells
top_markers_all

# cell type specific genes
known_genes <- c(
    "HBA1", # erythrocytes
    "CST3", # monocytes
    "CD3E", # T cells
    "NKG7", # NK T cells
    "CD79A", # B cells
    "MS4A1" # CD20 B cells
)

# violin plot
VlnPlot(seurat_object,
        features = known_genes,
        group.by = "leiden_cluster",
        pt.size = 0) +
    NoLegend() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1))


# add new labels to the Idents column in the metadata
seurat_object$Idents <- recode_values(
    seurat_object$leiden_cluster,
    "1" ~ "B (c1)",
    "2" ~ "T (c2)",
    "3" ~ "B (c3)",
    "4" ~ "B (c4)",
    "5" ~ "T (c5)",
    "6" ~ "B (c6)",
    "7" ~ "B (c7)",
    "8" ~ "B (c8)",
    "9" ~ "CD20+ B (c9)",
    "10" ~ "T (c10)",
    "11" ~ "NK T (c11)",
    "12" ~ "B (c12)",
    "13" ~ "Monocytes (c13)",
    "14" ~ "Erythrocytes (c14)",
    "15" ~ "Erythrocytes (c15)",
    "16" ~ "Unknown (c16)"
)

# set the Idents to the new labels
Idents(seurat_object) <- "Idents"

# visualise UMAP with new labels
DimPlot(seurat_object,
        reduction = "umap",
        group.by = "Idents",
        label = TRUE,
        label.size = 5) +
    NoLegend()





################################################################################
#                    10. Differential Gene Expression
################################################################################

## Comparison of interest: cells from healthy (PBMMC) vs. leukemia (ETV6-RUNX1)

# Load the libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(DESeq2) # for DGE

## This data is only for 4 ETV6RUNX1 and 3 Ctrl PBMMC samples
## Contains all cells per sample (not only 500 for each)
seurat_object <- readRDS("RObjects/Annotated.full.ETV6.PBMMC.rds")
seurat_object

## 30167 cells
dim(seurat_object)
# [1] 28343 30167

## Num of cells per sample
table(seurat_object$orig.ident)
# ETV6RUNX1-1 ETV6RUNX1-2 ETV6RUNX1-3 ETV6RUNX1-4     PBMMC-1     PBMMC-2     PBMMC-3
#        2715        6533        5012        5991         835        4862        4219

## We have 23 manually annotated cell type clusters from Leiden algorithm with k=48
table(seurat_object$Idents)
seurat_object$Idents %>% unique() %>% length()
# [1] 23

# Plot UMAP to visualize cell types across these 7 samples
p1 <- DimPlot(seurat_object,
        reduction = "umap",
        group.by = "Idents",
        label = F)

## Color by donor
p2 <- DimPlot(seurat_object,
        reduction = "umap",
        group.by = "orig.ident",
        label = F)

p1 + p2

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                       0.1 Creating pseudo-bulk samples
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Tabulate number of cells per label + sample combination
table(seurat_object$Idents, seurat_object$SampleName)
table(seurat_object$SampleName, seurat_object$Idents)

## Note some cell type clusters are almost entirely made up of a few donors
## This is plausible given these clusters correspond to cell types
## present in other clusters from other donors (no cell type is donor-specific)

## But could also makes us go back to batch correction


# Aggregate counts by sample and celltype
pseudo_bulk <- AggregateExpression(seurat_object,
                                   group.by = c("Idents", "SampleName"),
                                   assays = "RNA",
                                   return.seurat = TRUE)
# check the new assay
pseudo_bulk

# An object of class Seurat
# 33869 features across 142 samples within 1 assay
# Active assay: RNA (33869 features, 0 variable features)
# 2 layers present: counts, data

## 142 unique samples out of 23x7 = 161 donor-cluster combinations
## because some had 0 cells
length(unique(pseudo_bulk$orig.ident))
# [1] 142

## Check metadata
head(pseudo_bulk@meta.data)

## Add SampleGroup (condition) for DGE
temp_metadata <- pseudo_bulk@meta.data %>%
    mutate(SampleGroup = str_remove(SampleName, "-.*")) %>%
    mutate(Cluster = Idents)

pseudo_bulk@meta.data <- temp_metadata
# check the new metadata
head(pseudo_bulk@meta.data)


## Check SampleGroup levels
unique(pseudo_bulk$SampleGroup)
# [1] "ETV6RUNX1" "PBMMC"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                        10.2 Filtering Pseudosamples
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# filter out pseudosamples with less than 20 cells
# determine which pseudosamples have at least than 20 cells
ps_keep <- table(seurat_object$Idents, seurat_object$SampleName) %>%
    as.data.frame() %>%
    filter(Freq >= 20) %>%
    mutate(PseudoSample = str_c(Var1, "_", Var2)) %>%
    pull(PseudoSample)

# subset our pseudo-bulk object to keep only the pseudosamples with at least 20 cells
## Rownames of samples in pseudobulk are as in ps_keep
head(ps_keep)
pseudo_bulk@assays[["RNA"]]$counts[1:4, 1:4]

pseudo_bulk <- subset(pseudo_bulk, cells = ps_keep)

## 88 pseudosamples left
pseudo_bulk
# An object of class Seurat
# 33869 features across 88 samples within 1 assay
# Active assay: RNA (33869 features, 0 variable features)
# 2 layers present: counts, data


## Filter lowly expressed genes:

# filter out genes with less than 5 counts across all pseudosamples
# determine which genes have at least 5 counts across all pseudosamples
genes_keep <- rowSums(pseudo_bulk@assays$RNA$counts) >= 5
# subset our pseudo-bulk object to keep only the genes with at least 5 counts across
pseudo_bulk <- subset(pseudo_bulk, features = names(genes_keep)[genes_keep])

dim(pseudo_bulk)
# [1] 27930    88



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  10.3 Differential Testing of the whole dataset (no cell type cluster resolution)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Set Idents to SampleGroup (condition) to let FindMarkers() know what's
## the contrast to compare
Idents(pseudo_bulk)
Idents(pseudo_bulk) <- pseudo_bulk$SampleGroup

## We are comapring 43 vs 45 samples
table(pseudo_bulk$SampleGroup)
# ETV6RUNX1     PBMMC
#        43        45

## Num of cells per condition
table(seurat_object$Idents, seurat_object$SampleName, seurat_object$SampleGroup) %>%
         as.data.frame() %>%
         filter(Freq >= 20) %>% group_by(Var3) %>% summarise(count = sum(Freq))
# A tibble: 2 × 2
#   Var3      count
#   <fct>     <int>
# 1 ETV6RUNX1 19955
# 2 PBMMC      9821


# Run DESeq2 on pseudobulk samples
de <- FindMarkers(pseudo_bulk,
                  ident.1 = "ETV6RUNX1",
                  ident.2 = "PBMMC",
                  test.use = "DESeq2")

## DESeq2 will take raw counts!! not SCT counts, nor PCs nor batch-corrected PCs

# converting counts to integer mode <- to fit NB model
# gene-wise dispersion estimates <- initial gene dispersion estimates based on mean gene expression
# mean-dispersion relationship <- compute relationship between mean and dispertion across genes
# final dispersion estimates <- use relationship across genes to improve per-gene dispertion estimate
##                              Then estimate model coeffs

head(de)

#  - p_val: p-value of the Wald test (logFC diff from 0)
#  - avg_log2FC: log2 fold change between the two groups
#  - pct.1: percentage of pseudobulk samples in group 1 (ETV RUNX1) that have non-zero counts for the gene
#  - pct.2: percentage of pseudobulk samples in group 2 (PBMMC ) that have non-zero counts for the gene
#  - p_val_adj: adjusted p-value for multiple testing



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - -
#            10.4 Differential Testing within cell type clusters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - -

## DGE for Cases vs Ctrls within each cluster

## Number of cases and controls in each
table(pseudo_bulk$Cluster, pseudo_bulk$SampleGroup)
#                    ETV6RUNX1 PBMMC
# B (c1)                     2     3
# B (c12)                    3     1
# B (c16)                    2     0    <- not possible to perform DGE
# B (c22)                    0     2
# B (c23)                    1     2
# B (c3)                     1     3
# B (c4)                     3     3
# B (c5)                     4     3
# B (c6)                     4     3
# B (c7)                     2     0

## Number of cells per condition x cluster
table(seurat_object$Idents, seurat_object$SampleName, seurat_object$SampleGroup) %>%
         as.data.frame() %>%
         filter(Freq >= 20) %>% dplyr::group_by(Var1, Var3) %>%
         dplyr::summarise(total_cells = sum(Freq), .groups = "drop") %>% as.data.frame()

#                  Var1      Var3 total_cells
# 1              B (c1) ETV6RUNX1        4450
# 2              B (c1)     PBMMC        1011
# 3             B (c12) ETV6RUNX1        1740
# 4             B (c12)     PBMMC          67
# 5             B (c16) ETV6RUNX1         955
# 6             B (c22)     PBMMC          69
# 7             B (c23) ETV6RUNX1          23
# 8             B (c23)     PBMMC          78
# 9              B (c3) ETV6RUNX1        2183
# 10             B (c3)     PBMMC         175
# 11             B (c4) ETV6RUNX1        1562
# 12             B (c4)     PBMMC         290
# 13             B (c5) ETV6RUNX1        1476
# 14             B (c5)     PBMMC         247
# 15             B (c6) ETV6RUNX1        1136
# 16             B (c6)     PBMMC         353
# 17             B (c7) ETV6RUNX1         539

# Run DESeq2 on pseudobulk samples for cluster 6

## Visual contrast
p1 <- DimPlot(subset(seurat_object, subset = Idents == "B (c6)"),
        reduction = "umap",
        group.by = "orig.ident",
        label = F)

p2 <- DimPlot(subset(seurat_object, subset = Idents == "B (c6)"),
              reduction = "umap",
              group.by = "SampleGroup",
              label = F)
p1 + p2

# PCA plot of pseudo-bulk samples
c6 <- subset(pseudo_bulk, subset = Cluster == "B (c6)")
c6_counts <-  GetAssayData(c6, assay = "RNA", layer = "counts")
pca <- prcomp(t(c6_counts))
pca_df <- data.frame(pca$x, SampleGroup = c6$SampleGroup)
ggplot(pca_df, aes(x = PC1, y = PC2, color = SampleGroup)) +
    geom_point() +
    theme_classic()

# add cluster information to the Idents of our pseudo-bulk object
pseudo_bulk$Cluster_SampleGroup <- str_c(pseudo_bulk$Cluster,
                                         pseudo_bulk$SampleGroup, sep = "_")
Idents(pseudo_bulk) <- pseudo_bulk$Cluster_SampleGroup

de_cluster6 <- FindMarkers(pseudo_bulk,
                           ident.1 = "B (c6)_ETV6RUNX1",
                           ident.2 = "B (c6)_PBMMC",
                           test.use = "DESeq2")
head(de_cluster6)


# MA plot of DE results for cluster 6
# calculate mean expression for each gene
mean_expression <- rowMeans(c6_counts)
# add mean expression to DE results
de_cluster6$mean_expression <- mean_expression[rownames(de_cluster6)]
# plot MA plot
ggplot(de_cluster6, aes(x = log2(mean_expression), y = avg_log2FC)) +
    geom_point() +
    theme_classic() +
    xlab("Log2 mean expression") +
    ylab("Log2 fold change") +
    ggtitle("MA plot for cluster 6")

# p-value histogram for DE results for cluster 6
hist(de_cluster6$p_val, breaks = 50,
     main = "P-value histogram for cluster 6",
     xlab = "P-value")

# filter the DE results for cluster 6
de_cluster6_sig <- de_cluster6 %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) >1)

dim(de_cluster6_sig)

# expression plots for one of the top differentially expressed genes
top_gene <- "MDK"
# counts for the top gene
top_gene_counts <- GetAssayData(pseudo_bulk,
                                assay = "RNA",
                                layer = "counts")[top_gene, ]

# plot expression of the top gene across the clusters and sample groups
top_gene_df <- data.frame(Expression = top_gene_counts,
                          Cluster = pseudo_bulk$Cluster,
                          SampleGroup = pseudo_bulk$SampleGroup,
                          SampleName = pseudo_bulk$SampleName)

top_gene_df %>% head()

ggplot(top_gene_df,
       aes(x = SampleName, y = Expression, color = SampleGroup)) +
    facet_wrap(~ Cluster) +
    geom_point() +
    xlab("Cluster") +
    ylab("Expression of MDK") +
    ggtitle("Expression of MDK across clusters and sample groups") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Visualize top gene in UMAP
FeaturePlot(seurat_object, features = top_gene, split.by = "SampleGroup")


###################
# Exercise 1:
###################
# Run DESeq2 on pseudobulk samples for cluster 10
de_cluster10 <- FindMarkers(pseudo_bulk,
                            ident.1 = "CD20+ B (c10)_ETV6RUNX1",
                            ident.2 = "CD20+ B (c10)_PBMMC",
                            test.use = "DESeq2")

# filter the DE results for cluster 10
de_cluster10_sig <- de_cluster10 %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) >1)

dim(de_cluster10_sig)
head(de_cluster10_sig)


# plot the top differentially expressed gene for cluster 10
top_gene <- rownames(de_cluster10_sig)[1]
top_gene
top_gene_counts <- GetAssayData(pseudo_bulk,
                                assay = "RNA",
                                layer = "counts")[top_gene, ]
# plot expression of the top gene across the clusters and sample groups
top_gene_df <- data.frame(Expression = top_gene_counts,
                          Cluster = pseudo_bulk$Cluster,
                          SampleGroup = pseudo_bulk$SampleGroup,
                          SampleName = pseudo_bulk$SampleName)
ggplot(top_gene_df,
       aes(x = SampleName, y = Expression, color = SampleGroup)) +
    geom_point() +
    xlab("Cluster") +
    ylab("Expression of SOCS2") +
    ggtitle("Expression of SOCS2 across clusters and sample groups") +
    facet_wrap(~ Cluster) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))



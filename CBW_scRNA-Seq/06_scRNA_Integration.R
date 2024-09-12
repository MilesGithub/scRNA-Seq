# Load necessary libraries
library(Seurat)
library(SeuratDisk)
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratWrappers)
install.packages("/User/username/Downloads/satijalab-seurat-wrappers-d28512f.tar.gz",type='source',repos = NULL)
library(patchwork)
library(harmony)
library(rliger)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(ggplot2) # Ensure ggplot2 is loaded

# Source custom functions
source("./custom_seurat_functions.R")

# Download data files
download.file("https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_PBMC/Parent_NGSC3_DI_PBMC_filtered_feature_bc_matrix.h5",
              destfile = "3p_pbmc10k_filt.h5")
download.file("https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_filtered_feature_bc_matrix.h5",
              destfile = "5p_pbmc10k_filt.h5")

# Read data
matrix_3p <- Read10X_h5("D:/projects/CBW_scRNA/3p_pbmc10k_filt.h5", use.names = T)
matrix_5p <- Read10X_h5("./5p_pbmc10k_filt.h5", use.names = T)$`Gene Expression`

# Create Seurat objects
srat_3p <- CreateSeuratObject(matrix_3p, project = "pbmc10k_3p")
srat_5p <- CreateSeuratObject(matrix_5p, project = "pbmc10k_5p")

# Clean up environment
rm(matrix_3p)
rm(matrix_5p)

# Add percentage mitochondrial and ribosomal content
srat_3p[["percent.mt"]] <- PercentageFeatureSet(srat_3p, pattern = "^MT-")
srat_3p[["percent.rbp"]] <- PercentageFeatureSet(srat_3p, pattern = "^RP[SL]")
srat_5p[["percent.mt"]] <- PercentageFeatureSet(srat_5p, pattern = "^MT-")
srat_5p[["percent.rbp"]] <- PercentageFeatureSet(srat_5p, pattern = "^RP[SL]")

# Generate violin plots using ggplot2
VlnPlot(srat_3p, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), ncol = 4) & theme(
  panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
  panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

VlnPlot(srat_5p, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), ncol = 4) & theme(
  panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
  panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

# Check overlapping genes
table(rownames(srat_3p) %in% rownames(srat_5p))

# TRUE 
# 36601 

# Subset the data
srat_3p <- subset(srat_3p, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 15)
srat_5p <- subset(srat_5p, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

# Integrate data
pbmc_list <- list(srat_3p, srat_5p)
for (i in 1:length(pbmc_list)) {
  i<-2
  pbmc_list[[i]] <- NormalizeData(pbmc_list[[i]], verbose = FALSE)
  pbmc_list[[i]] <- FindVariableFeatures(pbmc_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

# Find integration anchors and integrate data
pbmc_anchors <- FindIntegrationAnchors(object.list = pbmc_list, dims = 1:30)
pbmc_seurat <- IntegrateData(anchorset = pbmc_anchors, dims = 1:30)
rm(pbmc_list)
rm(pbmc_anchors)

# Normalize, find variable features, and scale data
DefaultAssay(pbmc_seurat) <- "RNA"
pbmc_seurat <- NormalizeData(pbmc_seurat, verbose = FALSE)
pbmc_seurat <- FindVariableFeatures(pbmc_seurat, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pbmc_seurat <- ScaleData(pbmc_seurat, verbose = FALSE)

# Run PCA and UMAP
pbmc_seurat <- RunPCA(pbmc_seurat, npcs = 30, verbose = FALSE)
pbmc_seurat <- RunUMAP(pbmc_seurat, reduction = "pca", dims = 1:30, verbose = FALSE)

# Plot UMAP with ggplot2 and plot_annotation
DimPlot(pbmc_seurat, reduction = "umap") + plot_annotation(title = "10k 3' PBMC and 10k 5' PBMC cells, before integration")

# Switch assay and re-run analysis
DefaultAssay(pbmc_seurat) <- "integrated"
pbmc_seurat <- ScaleData(pbmc_seurat, verbose = FALSE)
pbmc_seurat <- RunPCA(pbmc_seurat, npcs = 30, verbose = FALSE)
pbmc_seurat <- RunUMAP(pbmc_seurat, reduction = "pca", dims = 1:30, verbose = FALSE)

# Plot integrated UMAP
p1<-DimPlot(pbmc_seurat, reduction = "umap") + plot_annotation(title = "10k 3' PBMC and 10k 5' PBMC cells, after integration (Seurat 3)") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p1

# Split by original identity
p2<-DimPlot(pbmc_seurat, reduction = "umap", split.by = "orig.ident") + NoLegend() +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p2

# Find neighbors and clusters
pbmc_seurat <- FindNeighbors(pbmc_seurat, dims = 1:30, k.param = 10, verbose = FALSE)
pbmc_seurat <- FindClusters(pbmc_seurat, verbose = FALSE)

# Plot clusters
p3<-DimPlot(pbmc_seurat, label = TRUE) + NoLegend() +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p3

# Generate count table
count_table <- table(pbmc_seurat@meta.data$seurat_clusters, pbmc_seurat@meta.data$orig.ident)
count_table

# pbmc10k_3p pbmc10k_5p
# 0           3       2579
# 1           0       2196
# 2        1532          3
# 3        1293          7
# 4           3       1269
# 5        1231          0
# 6        1220          4
# 7        1158          0
# 8           0       1000
# 9           0        632
# 10        546          5
# 11         40        444
# 12        454          0
# 13          0        405
# 14        360          2
# 15        340          2
# 16        257         53
# 17        291         11
# 18          1        292
# 19          0        273
# 20          0        253
# 21          3        218
# 22          3        185
# 23        130         35
# 24         93         61
# 25        115          0
# 26         64         13
# 27         23         43
# 28         15         10
# 29         15          5
# Plot integrated clusters

source('./custom_seurat_functions.R')

plot_integrated_clusters(pbmc_seurat)

# Clean up environment
rm(pbmc_seurat)

# Harmony integration
pbmc_harmony <- merge(srat_3p, srat_5p)
pbmc_harmony <- NormalizeData(pbmc_harmony, verbose = FALSE)
pbmc_harmony <- FindVariableFeatures(pbmc_harmony, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pbmc_harmony <- ScaleData(pbmc_harmony, verbose = FALSE)
pbmc_harmony <- RunPCA(pbmc_harmony, npcs = 30, verbose = FALSE)
pbmc_harmony <- RunUMAP(pbmc_harmony, reduction = "pca", dims = 1:30, verbose = FALSE)

# UMAP and violin plot with ggplot2 and patchwork
p4 <- DimPlot(pbmc_harmony, reduction = "pca", pt.size = 0.1, group.by = "orig.ident") + NoLegend() +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))
p5 <- VlnPlot(pbmc_harmony, features = "PC_1", group.by = "orig.ident", pt.size = 0.1) + NoLegend() +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))
plot_grid(p4, p5)

# Harmony integration
pbmc_harmony <- pbmc_harmony %>% RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(pbmc_harmony, 'harmony')
harmony_embeddings[1:5, 1:5]

# harmony_1  harmony_2   harmony_3   harmony_4 harmony_5
# AAACCCACATAACTCG-1_1   9.455486   2.427522   1.9122153  -1.1408724 -3.139317
# AAACCCACATGTAACC-1_1  -7.490588 -22.043615  -0.2698843   5.5109152 -2.653551
# AAACCCAGTGAGTCAG-1_1  18.486876  -3.422324  -4.4445333   1.7140767  9.570971
# AAACGAACAGTCAGTT-1_1  17.946506 -15.267514 -14.5427858 -33.9352923 14.068186
# AAACGAACATTCGGGC-1_1 -11.331866   2.417476   3.4331800  -0.4000288  1.094145

# UMAP after harmony
p6 <- DimPlot(pbmc_harmony, reduction = "harmony", pt.size = 0.1, group.by = "orig.ident") + NoLegend()  +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))
p7 <- VlnPlot(pbmc_harmony, features = "harmony_1", group.by = "orig.ident", pt.size = 0.1) + NoLegend()  +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))
plot_grid(p6, p7)

pbmc_harmony <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>%
  FindClusters() %>%
  identity()

# Plot UMAP after harmony
pbmc_harmony <- SetIdent(pbmc_harmony, value = "orig.ident")
p8<-DimPlot(pbmc_harmony, reduction = "umap") + plot_annotation(title = "10k 3' PBMC and 10k 5' PBMC cells, after integration (Harmony)")  +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p8

p9<-DimPlot(pbmc_harmony, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, split.by = 'orig.ident') + NoLegend()  +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p9

# Plot clusters after harmony
pbmc_harmony <- SetIdent(pbmc_harmony, value = "seurat_clusters")
p10<-DimPlot(pbmc_harmony, label = TRUE) + NoLegend() +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p10

plot_integrated_clusters(pbmc_harmony)

# Clean up environment
rm(pbmc_harmony)

# LIGER integration
pbmc_liger <- merge(srat_3p, srat_5p)
pbmc_liger <- pbmc_liger %>% rliger::normalize() %>% selectGenes() %>% scaleNotCenter()

# Run iNMF
install.packages('RcppPlanc', repos = c('https://welch-lab.r-universe.dev', 'https://cloud.r-project.org'))
pbmc_liger <- optimizeALS(pbmc_liger, k = 30)
pbmc_liger <- quantileAlignSNF(pbmc_liger)

# Plot iNMF clusters with ggplot2 and patchwork
pbmc_liger <- pbmc_liger %>% liger::runUMAP(reduction = "iNMF", dims = 1:30)
plotByDatasetAndCluster(pbmc_liger, return.plots = TRUE)

# Plot UMAP after LIGER integration with ggplot2 and plot_annotation
p1 <- DimPlot(pbmc_liger, reduction = "iNMF", group.by = "dataset") + NoLegend()
p2 <- DimPlot(pbmc_liger, reduction = "iNMF", group.by = "cluster") + NoLegend()
plot_grid(p1, p2)

plot_integrated_clusters(pbmc_liger)

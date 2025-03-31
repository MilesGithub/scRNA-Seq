# 1. Install and Load Libraries
# ---------------------------------------------------

# List of packages to check and install
packages <- c('S4Vectors', 'SingleCellExperiment', 'patchwork',
              'SummarizedExperiment', 'SeuratData', 'Seurat', 'devtools', 'dplyr', 'ggplot2')

# Function to check if a package is installed
load_packages <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% rownames(installed.packages())) {
      print(paste0("Package already installed: ", pkg))
    } else if (pkg %in% BiocManager::available()) {
      BiocManager::install(pkg)
      print(paste0("Installed from Bioconductor: ", pkg))
    } else {
      install.packages(pkg)
      print(paste0("Installed from CRAN: ", pkg))
    }
  }
  library(pkg, character.only = TRUE)
  print(paste0("Loaded: ", pkg))
}

# Apply function to each package
lapply(packages, load_packages)

# # For the example dataset:
# install.packages("SeuratData") # Provides easy access to example datasets
# SeuratData::InstallData("pbmc3k") # Download the pbmc3k dataset

# Set a seed for reproducibility (optional)
set.seed(1234)

print(paste("Seurat version:", packageVersion("Seurat")))
print(paste("Current Date:", Sys.Date())) # Seurat analysis can be influenced by package versions

# 2. Load Data and Create Seurat Object
# -------------------------------------

# Load the PBMC3k dataset using SeuratData
# This returns a Seurat object directly
data("pbmc3k")

# Assign it to a shorter variable name
pbmc <- pbmc3k

# Seurat object contains the raw count matrix and basic metadata
print("Initial Seurat object:")
print(pbmc)
# An object of class Seurat
# 13714 features across 2700 samples within 1 assay
# Active assay: RNA (13714 features, 0 variable features)

# The raw counts are stored in pbmc[['RNA']]$counts or pbmc@assays$RNA$counts
# Cell-level metadata is in pbmc@meta.data

# 3. Preprocessing and Quality Control (QC)
# -----------------------------------------

# Calculate the percentage of mitochondrial counts for each cell
# Mitochondrial genes in humans typically start with "MT-"
# High percentage can indicate stressed or dying cells
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics before filtering
# nFeature_RNA: number of genes detected per cell
# nCount_RNA: total number of molecules (UMIs) detected per cell
# percent.mt: percentage of counts mapping to mitochondrial genes
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("seurat_qc_violin_before_filtering.png")

# Visualize relationships between metrics
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("seurat_qc_scatter_before_filtering.png")

# Filter cells based on QC metrics
# These thresholds are common starting points but should be adjusted based on the plots
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

print("Seurat object after QC filtering:")
print(pbmc)
# An object of class Seurat
# 13714 features across 2638 samples within 1 assay
# Active assay: RNA (13714 features, 0 variable features)

# 4. Normalization
# ----------------

# Normalize the data to account for differences in library size (sequencing depth)
# "LogNormalize": Normalizes feature counts for each cell by total counts,
# multiplies by a scale factor (default 10,000), and log-transforms the result.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# The normalized data is stored in pbmc[['RNA']]$data or pbmc@assays$RNA$data

# 5. Feature Selection (Identify Highly Variable Genes - HVGs)
# ------------------------------------------------------------

# Identify genes that exhibit high cell-to-cell variation
# Most informative for distinguishing cell types/states
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) # Find top 2000 HVGs

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# Plot variable features with labels for the top 10
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
ggsave("seurat_hvg_plot.png")

# The list of variable features is stored: VariableFeatures(pbmc)

# 6. Scaling the Data
# -------------------

# Scale the data so that mean expression across cells is 0 and variance is 1
# Only for the variable features by default, important for PCA
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes) # Scale all genes for potential heatmap visualization later
# Scaling is typically applied only to Variable Features before PCA,
# but scaling all makes downstream visualization of any gene easier.
# For PCA specifically, RunPCA uses only VariableFeatures by default if ScaleData was run on them.

# Alternatively, scale only variable features and regress out confounders:
# pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc), vars.to.regress = "percent.mt")

# Scaled data is stored in pbmc[['RNA']]$scale.data

# 7. Linear Dimensionality Reduction (PCA)
# ----------------------------------------

# Perform Principal Component Analysis (PCA) on the scaled data (using HVGs by default)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50, verbose = FALSE)

# Examine and visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5) # Show top genes associated with first 5 PCs

# Visualize feature loadings for first few PCs
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
ggsave("seurat_pca_loadings.png")

# Scatter plot of cells based on first two PCs
DimPlot(pbmc, reduction = "pca")
ggsave("seurat_pca_dimplot.png")

# Heatmap exploring sources of heterogeneity in top PCs
# Helps determine which PCs capture biological signal vs. noise/technical artifacts
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
ggsave("seurat_pca_dimheatmap.png")

# Determine the 'dimensionality' of the dataset (how many PCs to use downstream)
# 'Elbow plot': Look for an 'elbow' where the variance explained by PCs drops off
ElbowPlot(pbmc, ndims = 50) # Examine the first 50 PCs
ggsave("seurat_pca_elbowplot.png")
# Based on the elbow plot, choose the number of dimensions (e.g., 10-20 often reasonable for PBMCs)


# 8. Clustering
# -------------

# Cluster cells based on their PCA space representation
# First, construct a K-Nearest Neighbor (KNN) graph based on Euclidean distance in PCA space
# Then, apply community detection algorithm (default is Louvain, option for Leiden)
dims_to_use <- 1:15 # Choose dimensions based on Elbow Plot
pbmc <- FindNeighbors(pbmc, dims = dims_to_use)

# Find clusters using a resolution parameter (higher value -> more clusters)
pbmc <- FindClusters(pbmc, resolution = 0.5) # 0.5 is a common starting point

# Cluster IDs are stored in pbmc@meta.data$seurat_clusters
print("Cluster sizes:")
print(table(pbmc@meta.data$seurat_clusters))

# 9. Non-linear Dimensionality Reduction (UMAP/tSNE)
# --------------------------------------------------

# Run UMAP or t-SNE for visualization (typically in 2D)
# Uses the same PCs as input for clustering
pbmc <- RunUMAP(pbmc, dims = dims_to_use, verbose = FALSE)
# pbmc <- RunTSNE(pbmc, dims = dims_to_use, verbose = FALSE) # Alternative

# Visualize clusters on UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()
ggsave("seurat_umap_clusters.png")

# Can also visualize metadata on UMAP plot
# DimPlot(pbmc, reduction = "umap", group.by = "percent.mt")

# 10. Finding Cluster Marker Genes
# -------------------------------

# Find genes that are differentially expressed between clusters
# Identifies marker genes that define each cluster

# Find markers for every cluster compared to all remaining cells (one-vs-all)
# Default test is Wilcoxon Rank Sum test
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View top markers per cluster
print("Top markers per cluster:")
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>% # Filter for stronger markers if needed
  slice_head(n = 10) %>% # Get top 10 per cluster
  ungroup()
print(top_markers)

# Find markers for a specific cluster (e.g., cluster 2 vs all others)
# cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
# head(cluster2.markers, n = 5)

# Visualize marker gene expression
# Violin plots show expression distribution across clusters
VlnPlot(pbmc, features = c("MS4A1", "CD79A")) # Example B cell markers
ggsave("seurat_marker_vlnplot.png")

# Feature plots show expression overlaid on UMAP/tSNE
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
ggsave("seurat_marker_featureplot.png")

# Heatmap of top markers across clusters
top10_per_cluster <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10_per_cluster$gene) + NoLegend()
ggsave("seurat_marker_heatmap.png")

# 11. Assigning Cell Type Annotations (Conceptual)
# -----------------------------------------------

# Based on the identified marker genes (e.g., from FeaturePlot, VlnPlot, marker table),
# compare them to known markers from literature to assign cell type identities to clusters.
# Example:
# Cluster expressing CD14, LYZ might be CD14+ Monocytes
# Cluster expressing MS4A1 might be B cells
# Cluster expressing IL7R, CD3D might be T cells
# Cluster expressing GNLY, NKG7 might be NK cells

# Create a vector mapping cluster IDs to cell type names (Example)
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet") # This mapping corresponds to the PBMC3k tutorial result
names(new.cluster.ids) <- levels(pbmc) # Get the current cluster IDs (0, 1, 2, ...)

# Add the annotations to the Seurat object metadata
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Visualize annotated clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("seurat_umap_annotated.png")

# Store annotations in metadata if desired
# pbmc$celltype <- Idents(pbmc)

# 12. Saving Results
# ------------------
# Save the final annotated Seurat object
# saveRDS(pbmc, file = "pbmc3k_final.rds")
# print("Analysis complete. Seurat object saved to pbmc3k_final.rds")

# To load later:
# pbmc_loaded <- readRDS("pbmc3k_final.rds")

print("Basic Seurat analysis workflow completed.")
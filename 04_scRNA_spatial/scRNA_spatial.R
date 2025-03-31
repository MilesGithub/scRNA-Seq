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

# # Required for stxBrain dataset
# SeuratData::InstallData("stxBrain")

# Set a seed for reproducibility (optional)
set.seed(1234)

print(paste("Seurat version:", packageVersion("Seurat")))

# 2. Load Spatial Data
# --------------------

# Load the Anterior1 slice of the mouse brain dataset from SeuratData
# This function directly loads a pre-processed Seurat object containing spatial info
data("stxBrain", package = "SeuratData")
# Let's use the 'anterior1' slice from this dataset
brain <- stxBrain$anterior1

# Explore the Seurat object - notice it contains image data
print("Loaded Spatial Seurat object:")
print(brain)
# An object of class Seurat
# 31053 features across 2696 samples within 1 assay
# Active assay: Spatial (31053 features, 0 variable features)
# 1 image present: anterior1

# If loading from 10x Space Ranger output (typical user scenario):
# You would have folders like 'outs/filtered_feature_bc_matrix/' and 'outs/spatial/'
# brain <- Load10X_Spatial(
#   data.dir = "/path/to/your/spaceranger_output/outs/",
#   filename = "filtered_feature_bc_matrix.h5",
#   assay = "Spatial", # specify name for the assay
#   slice = "slice1", # specify name for the image slice
#   filter.matrix = TRUE,
#   to.lower = FALSE
# )

# 3. Quality Control (QC) and Initial Visualization
# -------------------------------------------------

# Visualize spot locations and basic histology
# ImageDimPlot(brain, cols = "grey", axes = TRUE) # Use this if loaded with Load10X_Spatial

# Calculate mitochondrial percentage (mouse mitochondrial genes start with "mt-")
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^mt-")

# Visualize QC metrics using standard violin plots
VlnPlot(brain, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), ncol = 3, pt.size = 0.1) +
  NoLegend()
ggsave("spatial_qc_violin.png")

# Visualize QC metrics spatially on the tissue slice
# This helps identify regions with low quality (e.g., edge effects, tissue damage)
plot1 <- SpatialFeaturePlot(brain, features = "nCount_Spatial", pt.size.factor = 1.6) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(brain, features = "nFeature_Spatial", pt.size.factor = 1.6) + theme(legend.position = "right")
plot3 <- SpatialFeaturePlot(brain, features = "percent.mt", pt.size.factor = 1.6) + theme(legend.position = "right")
wrap_plots(plot1, plot2, plot3, ncol=3)
ggsave("spatial_qc_featureplots.png")

# Filter spots based on QC metrics (example thresholds, adjust based on plots)
# Keep spots with reasonable counts and low mitochondrial percentage
brain <- subset(brain, subset = nFeature_Spatial > 500 & nFeature_Spatial < 7500 & percent.mt < 15)

print("Seurat object after QC filtering:")
print(brain)


# 4. Normalization and Dimensionality Reduction with SCTransform
# -----------------------------------------------------------------------

# Normalization, variance stabilization, and feature selection with SCTransform
# Replaces NormalizeData, FindVariableFeatures, ScaleData steps
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

# Perform PCA on the SCT assay
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)

# Perform UMAP for visualization
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30, verbose = FALSE) # Use significant PCs

# 5. Clustering
# -------------

# Cluster spots based on gene expression similarity in PCA space
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, resolution = 0.8, verbose = FALSE) # Adjust resolution as needed

# 6. Spatial Visualization of Clusters
# ------------------------------------

# Visualize the clusters directly on the tissue slice coordinates
# This is a key advantage of spatial data - seeing cluster distribution
plot1 <- DimPlot(brain, reduction = "umap", label = TRUE) + NoLegend() + ggtitle("UMAP Clusters")
plot2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3, pt.size.factor = 1.6) + ggtitle("Spatial Clusters")
plot1 + plot2
ggsave("spatial_umap_vs_spatial_clusters.png")

# Visualize specific clusters or metadata spatially
# SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 5)), facet.highlight = TRUE, ncol = 3)


# 7. Identification of Spatially Variable Features (SVFs)
# -------------------------------------------------------

# Genes with strong spatial pattern expression, independent of clusters
# Reveals gradients or region-specific markers missed by clustering alone

# Seurat offers multiple methods, 'markvariogram' is default and often works well
# It models spatial variance based on a variogram
brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
                                       selection.method = "markvariogram")

# View the top spatially variable features
top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 6)
print("Top 6 Spatially Variable Features:")
print(top.features)

# Visualize the expression patterns of these top SVFs spatially
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1), pt.size.factor = 1.6)
ggsave("spatial_top_svf_featureplots.png")


# 8. Relating Clusters to Marker Genes and Spatial Location
# ---------------------------------------------------------

# Find differentially expressed genes between the spatial clusters
# Using the SCT assay data
DefaultAssay(brain) <- "SCT" # Ensure we use SCT for markers
de_markers <- FindAllMarkers(brain, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2, verbose = FALSE)

# Get top markers per cluster
top5_markers <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
print("Top 5 markers per cluster:")
print(top5_markers)

# Visualize expression of key marker genes spatially
# Example: Assuming markers suggest neuronal layers or regions based on known biology
SpatialFeaturePlot(brain, features = c("Ttr", "Plp1", "Hpca", "Gja1"), ncol=2, pt.size.factor = 1.6)
ggsave("spatial_marker_examples_featureplot.png")

# Can also visualize using violin plots or heatmaps as in standard scRNA-seq
# VlnPlot(brain, features = c("Ttr", "Plp1"), group.by = "seurat_clusters", pt.size = 0.1)


# 9. Saving Results
# -----------------
# Save the processed Seurat object with spatial information and analysis results
# saveRDS(brain, file = "mouse_brain_anterior1_processed.rds")
# print("Spatial analysis complete. Seurat object saved.")

print("Basic Seurat spatial analysis workflow completed.")
# 1. Install and Load Libraries
# ---------------------------------------------------

# List of packages to check and install
packages <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'S4Vectors', 'SingleCellExperiment','SummarizedExperiment', 
              'dplyr', 'ggplot2', 'Seurat', 'SeuratData', 'Matrix', 'Matrix.utils', 'DESeq2', 'patchwork')

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


# For reproducibility
set.seed(1234)


# 2. Load citeseq_cbmc Data
# -------------------------

# AvailableData()
# InstallData("citeseq_cbmc")

# Load the dataset
data("cbmc", package = "SeuratData")

print("Loaded citeseq_cbmc dataset:")
print(cbmc)
print(paste("Default assay:", DefaultAssay(cbmc)))


# 3. Explore Metadata
# -------------------

# Check the metadata for sample and condition information
print("\nMetadata columns:\n")
print(colnames(cbmc@meta.data))
print("\nUnique Sample Identifiers (orig.ident):\n")
print(unique(cbmc$orig.ident))
print("\nConditions (stim):\n")
print(table(cbmc$stim))
print("\nSamples per condition:\n")
print(table(cbmc$stim, cbmc$orig.ident))

# Identify metadata columns:
# 'orig.ident' likely represents the biological repliprinte (donor sample)
# 'stim' represents the condition (CTRL vs STIM)
sample_id_col <- "orig.ident"
condition_col <- "stim"

# Set default assay to RNA if it isn't already (it usually is for this dataset)
DefaultAssay(cbmc) <- "RNA"


# 4. Standard Preprocessing and QC (RNA Assay)
# --------------------------------------------

# Calculate mitochondrial percentage (using human MT gene pattern "^MT-")
cbmc[["percent.mt"]] <- PercentageFeatureSet(cbmc, pattern = "^MT-", assay="RNA")

# Visualize QC metrics
VlnPlot(cbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1, group.by=condition_col) + NoLegend()

# Filter cells (adjust thresholds based on VlnPlots for CBMCs)
min_features <- 200
max_features <- 4000 # Example threshold, adjust based on data
max_mt_percent <- 10 # Example threshold, adjust based on data

cbmc <- subset(cbmc, subset = nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mt < max_mt_percent)

print("\nDimensions after QC filtering:\n")
print(cbmc)


# 5. Normalization, Scaling, Dimensionality Reduction, Clustering (RNA Assay)
# ---------------------------------------------------------------------------

cbmc <- NormalizeData(cbmc, assay = "RNA")
cbmc <- FindVariableFeatures(cbmc, selection.method = "vst", nfeatures = 2000, assay = "RNA")
cbmc <- ScaleData(cbmc, features = rownames(cbmc), assay = "RNA") # Scale all genes
cbmc <- RunPCA(cbmc, features = VariableFeatures(object = cbmc), npcs = 30, verbose = FALSE, assay = "RNA")
# ElbowPlot(cbmc, ndims=30)

# Choose appropriate number of dimensions
pca_dims <- 1:15 # Adjust based on ElbowPlot/variance explained

cbmc <- FindNeighbors(cbmc, dims = pca_dims, assay = "RNA")
cbmc <- FindClusters(cbmc, resolution = 0.6) # Adjust resolution as needed
cbmc <- RunUMAP(cbmc, dims = pca_dims, assay = "RNA", reduction="pca", reduction.name="umap.rna") # Specify names to avoid conflict if ADT UMAP exists

# Visualize clusters and conditions
DimPlot(cbmc, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE) + NoLegend() + ggtitle("RNA Clusters")
DimPlot(cbmc, reduction = "umap.rna", group.by = condition_col) + ggtitle("Conditions (CTRL vs STIM)")
DimPlot(cbmc, reduction = "umap.rna", group.by = sample_id_col) + ggtitle("Samples (Donors)") + NoLegend()


# 6. Cell Type Annotation (Using Clusters as Proxy)
# -------------------------------------------------

# Assign clusters as proxy cell types
cbmc$cell_type_proxy <- paste0("Cluster", cbmc$seurat_clusters)
print("\nUsing Seurat clusters as proxy cell types:\n")
print(table(cbmc$cell_type_proxy))
# Use marker genes in real analysis to assign biological cell type names to these clusters.


# 7. Pseudobulk Aggregation
# -------------------------

# Select the cell type identifier to use
cell_type_column <- "cell_type_proxy"

# Aggregate counts per *biological sample* per cell type (cluster)
# Extract RNA counts and metadata needed
counts_matrix <- GetAssayData(cbmc, assay = "RNA", slot = "counts")
meta_data <- cbmc@meta.data %>% select(all_of(c(sample_id_col, cell_type_column))) # Use all_of for safety

# Create unique identifiers for aggregation (sample_celltype)
# Use the actual sample ID (orig.ident) and the proxy cell type
meta_data$sample_celltype <- paste(meta_data[[sample_id_col]], meta_data[[cell_type_column]], sep = "_")

# Check number of cells per group
print("\nNumber of cells per pseudobulk sample (Donor_Cluster):\n")
print(table(meta_data$sample_celltype))

# Sum counts using Matrix::tapply
pseudobulk_counts <- Matrix.utils::aggregate.Matrix(
  t(counts_matrix), # Transpose so cells are rows
  groupings = factor(meta_data$sample_celltype),
  fun = "sum"
)

# Transpose back to get genes x samples_celltype format
pseudobulk_counts <- t(pseudobulk_counts)

print("\nDimensions of pseudobulk matrix (Genes x Sample_CellType):\n")
print(dim(pseudobulk_counts))


# 8. Prepare Metadata for Pseudobulk Samples
# ------------------------------------------

# Create colData for DESeq2
pseudobulk_meta <- data.frame(
  sample_celltype = colnames(pseudobulk_counts),
  row.names = colnames(pseudobulk_counts)
) %>%
  mutate(
    # Extract sample_id (orig.ident) and cell_type carefully
    # Assuming sample IDs don't contain underscores before the cell type part
    !!sample_id_col := sub("_(Cluster\\d+)$", "", sample_celltype),
    cell_type = sub("^.*_", "", sample_celltype)
  ) %>%
  # Get the condition info associated with the original sample ID
  left_join(distinct(cbmc@meta.data[, c(sample_id_col, condition_col)]), by = sample_id_col)

# Reorder metadata rows to match count matrix columns if needed
pseudobulk_meta <- pseudobulk_meta[match(colnames(pseudobulk_counts), rownames(pseudobulk_meta)), ]

# Check consistency
stopifnot(all(rownames(pseudobulk_meta) == colnames(pseudobulk_counts)))

print("\nMetadata for pseudobulk samples (colData for DESeq2):\n")
print(head(pseudobulk_meta))
print(table(pseudobulk_meta$cell_type, pseudobulk_meta[[condition_col]])) # Check sample counts per cell type/condition


# 9. Differential Expression Analysis (Example: CD14 Mono proxy, STIM vs CTRL)
# ----------------------------------------------------------------------------

# Identify a cluster likely representing CD14 Monocytes (often express CD14, LYZ highly).
# Check UMAP/feature plots or marker gene analysis to choose cluster.
# Cluster 1 represents CD14 Monocytes for example
# FeaturePlot(cbmc, c("CD14", "LYZ", "MS4A1", "CD3D", "CD8A", "FCGR3A"), reduction="umap.rna")
target_cell_type <- "Cluster1"

# Filter the pseudobulk counts and metadata
counts_subset <- pseudobulk_counts[, pseudobulk_meta$cell_type == target_cell_type]
meta_subset <- pseudobulk_meta[pseudobulk_meta$cell_type == target_cell_type, ]

# Check if subsetting resulted in samples from only one condition
if (length(unique(meta_subset[[condition_col]])) < 2) {
  stop(paste("Error: Target cell type", target_cell_type, "only found in one condition. Cannot perform DE."))
}
# Check if there are enough repliprintes per condition
condition_counts <- table(meta_subset[[condition_col]])
if (any(condition_counts < 2)) {
  warning(paste("Warning: Fewer than 2 repliprintes in at least one condition for", target_cell_type, ". DESeq2 results might be unreliable or fail."))
  print(condition_counts)
}


stopifnot(all(rownames(meta_subset) == colnames(counts_subset)))

print(paste("\nRunning DESeq2 for:", target_cell_type, "(comparing STIM vs CTRL)\n"))

# Create DESeqDataSet object - design includes the condition
dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                              colData = meta_subset,
                              design = ~ stim) # Use the condition column name ('stim')

# Filter low count genes
min_total_count <- 10
keep <- rowSums(counts(dds)) >= min_total_count
dds <- dds[keep,]
print(paste("Genes kept after filtering low counts:", sum(keep), "\n"))

# Check if enough genes remain
if (nrow(dds) < 10) {
  stop("Error: Too few genes remaining after filtering. Check filtering thresholds or cell counts.")
}

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results - compare STIM vs CTRL
res <- results(dds, contrast=c("stim", "STIM", "CTRL")) # factor name, numerator, denominator
res <- na.omit(res)
resOrdered <- res[order(res$padj),]

# 10. View Results
# ----------------

print(paste("\nDESeq2 Results Summary for", target_cell_type, " (STIM vs CTRL):\n"))
summary(res)

print(paste("\nTop differentially expressed genes for", target_cell_type, "(STIM vs CTRL):\n"))
print(head(resOrdered)) # Expect IFN-response genes (e.g., ISG15, IFITM*) if target is Mono/B cell

# Optional: Volcano plot
# EnhancedVolcano::EnhancedVolcano(resOrdered, lab = rownames(resOrdered), x = 'log2FoldChange', y = 'padj', title = paste(target_cell_type, '- STIM vs CTRL'))

# 10. Saving Results
# -----------------
# saveRDS(resOrdered, file = "pseudobulk_DESeq_results.rds")


print("\nPseudobulk analysis example using citeseq_cbmc dataset complete.\n")
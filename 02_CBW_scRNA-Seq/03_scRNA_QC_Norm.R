library(DropletUtils)
library(SingleCellExperiment)
library(scater)
library(Seurat)
library(ggplot2)
library(patchwork)  # For combining plots

# Load the filtered 10X data into a SingleCellExperiment object
sce <- DropletUtils::read10xCounts("D:/projects/CBW_scRNA/scRNA_Data/HumanLiver/filtered_feature_bc_matrix/")

# Identify mitochondrial genes based on gene symbols
is.mt <- grepl("^MT-", rowData(sce)$Symbol)
sum(is.mt)  # Count of mitochondrial genes

# Calculate per-cell and per-feature QC metrics
cellQC <- perCellQCMetrics(sce, subsets = list(Mito = is.mt))
geneQC <- perFeatureQCMetrics(sce)

# Add the QC metrics to the colData and rowData of the SCE object
colData(sce) <- cbind(colData(sce), cellQC)
rowData(sce) <- cbind(rowData(sce), geneQC)

# Filter genes with at least 1% detection
keep.genes <- geneQC$detected > 0.01
sum(keep.genes)  # Count of genes kept after filtering

#[1] 27083

sce <- sce[keep.genes,]

# Plot the distribution of detected genes using ggplot2
p1<-ggplot(as.data.frame(colData(sce)), aes(x = detected)) +
  geom_histogram(bins = 100, fill = "blue", color = "black") +
  labs(x = "Number of Detected Genes", y = "Frequency") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p1

# Filter cells based on detected genes < 500
cell_filter_detect <- sce$detected < 500

# Plot the distribution of detected genes with a vertical line at the cutoff
p2<-ggplot(as.data.frame(colData(sce)), aes(x = detected)) +
  geom_histogram(bins = 100, fill = "blue", color = "black") +
  geom_vline(xintercept = 500, color = "red", linetype = "dashed") +
  labs(x = "Number of Detected Genes", y = "Frequency") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p2

# Filter cells based on mitochondrial content > 30%
cell_filter_MT <- sce$subsets_Mito_percent > 30

# Plot the distribution of mitochondrial content with a vertical line at the cutoff
p3<-ggplot(as.data.frame(colData(sce)), aes(x = subsets_Mito_percent)) +
  geom_histogram(bins = 100, fill = "blue", color = "black") +
  geom_vline(xintercept = 30, color = "red", linetype = "dashed") +
  labs(x = "Mitochondrial Content (%)", y = "Frequency") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p3

# Summarize the number of cells filtered by different criteria
sum(cell_filter_detect)  # Cells with low detected genes
#[1] 4590
sum(cell_filter_MT)  # Cells with high mitochondrial content
#[1] 4524
sum(cell_filter_detect & cell_filter_MT)  # Cells failing both criteria
#[1] 4086
sum(!(cell_filter_detect | cell_filter_MT))  # Cells passing both filters
#[1] 2596
# Mark cells to be discarded based on QC criteria
sce$discard <- (cell_filter_detect | cell_filter_MT)

# Plot QC metrics with discarded cells highlighted using ggplot2
p4<-ggplot(as.data.frame(colData(sce)), aes(x = sum, y = subsets_Mito_percent, color = discard)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "Total Counts", y = "Mitochondrial Content (%)") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))
p4

p5<-ggplot(as.data.frame(colData(sce)), aes(x = sum, y = detected, color = discard)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "Total Counts", y = "Number of Detected Genes") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p5

# Filter out the discarded cells
sce_filtered <- sce[, !sce$discard]

# Recalculate per-feature QC metrics after filtering
geneQC <- perFeatureQCMetrics(sce_filtered)

# Filter genes with at least 1% detection in the filtered data
keep.genes <- geneQC$detected > 0.01
sum(keep.genes)  # Count of genes kept after filtering
#[1] 26186

sce_filtered <- sce_filtered[keep.genes,]

# Create a Seurat object from the filtered 10X data
seur_obj <- CreateSeuratObject(Read10X("D:/projects/CBW_scRNA/scRNA_Data/HumanLiver/filtered_feature_bc_matrix/"))

# Normalize the SCE object and convert it to a Seurat object
sce_filtered <- logNormCounts(sce_filtered)
seur_obj2 <- as.Seurat(sce_filtered)

# Convert a Seurat object back to a SingleCellExperiment object
sce2 <- as.SingleCellExperiment(seur_obj)

# Calculate the percentage of mitochondrial genes in Seurat object
seur_obj[["percent.mt"]] <- PercentageFeatureSet(seur_obj, "^MT-")

# Create violin plots for quality control metrics using ggplot2
VlnPlot(seur_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
  theme_minimal()

# Create scatter plots for QC metrics and combine them using patchwork
plot1 <- FeatureScatter(seur_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  theme(legend.position = "none") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

plot2 <- FeatureScatter(seur_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  theme(legend.position = "none") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

plot1 + plot2  # Combine the scatter plots

# Mark cells that pass QC filters in Seurat object
seur_obj@meta.data$pass.filter <- seur_obj@meta.data$percent.mt < 30 & seur_obj@meta.data$nFeature_RNA > 500

# Create scatter plots with QC filtering highlighted using ggplot2
plot_mt <- FeatureScatter(seur_obj, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "pass.filter") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 30, color = "red", linetype = "dashed") +
  ggtitle("") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

plot_nfeature <- FeatureScatter(seur_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "pass.filter") +
  geom_hline(yintercept = 500, color = "red", linetype = "dashed") +
  ggtitle("") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

plot_mt + plot_nfeature  # Combine the scatter plots

# Subset the Seurat object based on QC filters
seur_obj_filtered <- subset(seur_obj, subset = nFeature_RNA > 500 & percent.mt < 30)

library(scran)
library(M3Drop)  # Note: Pearson residuals is only available in the development version of M3Drop on GitHub
library(Seurat)
library(ggplot2)

# Normalize the SingleCellExperiment object and the Seurat object
sce_filtered <- scuttle::logNormCounts(sce_filtered)
seur_obj_filtered <- NormalizeData(seur_obj_filtered)

# Compute log-normalized raw counts for PCA
assay(sce_filtered, "logcounts_raw") <- log1p(counts(sce_filtered))
sce_filtered <- runPCA(sce_filtered, exprs_values = "logcounts_raw")

# Plot PCA with ggplot2
assay(sce_filtered, "logcounts_raw") <- log2(counts(sce_filtered) + 1)
sce_filtered <- runPCA(sce_filtered, exprs_values = "logcounts_raw")
plotPCA(sce_filtered, colour_by = "sum", size_by = "detected")


# Perform quick clustering and compute size factors
qclust <- quickCluster(sce_filtered, min.size = 30)
sce_filtered <- computeSumFactors(sce_filtered, clusters = qclust)
sce_filtered <- logNormCounts(sce_filtered)

# Fit a negative binomial model and calculate Pearson residuals
fit <- NBumiFitModel(counts(sce_filtered))
pr <- NBumiPearsonResiduals(counts(sce_filtered), fit)
assay(sce_filtered, "pearson") <- pr

# Run PCA on Pearson residuals
sce_filtered <- runPCA(sce_filtered, exprs_values = "pearson")
plotPCA(sce_filtered, size_by = "detected", colour_by = "sum", ncomponents=1:2)

# Plot PCA with ggplot2
p1<-ggplot(pca_data, aes(x = PC1, y = PC2, color = sum, size = detected)) +
  geom_point() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

# Plot the relationship between mean counts and Pearson residual variance
mean_counts <- apply(counts(sce_filtered), 1, mean)
pr_var <- apply(pr, 1, var)
p2<-ggplot(data.frame(mean_counts, pr_var), aes(x = mean_counts, y = pr_var)) +
  geom_point() +
  labs(x = "Mean UMI counts", y = "Pearson Residual Variance") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p2

# Run PCA on the top 500 expressed genes based on Pearson residuals
top_expressed_genes <- tail(sort(mean_counts), 500)
sce_filtered <- runPCA(sce_filtered, exprs_values = "pearson", subset_row = names(top_expressed_genes))
pca_data <- plotPCA(sce_filtered, size_by = "detected", colour_by = "sum", ncomponents = 1:2)

# Plot PCA with ggplot2
p3<-ggplot(pca_data, aes(x = PC1, y = PC2, color = sum, size = detected)) +
  geom_point() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

remotes::install_version("matrixStats", version="1.1.0")


# Normalize data using SCTransform in Seurat
seur_obj_filtered <- SCTransform(seur_obj_filtered)

# Visualize the relationship between raw and SCT normalized gene counts
p4<-ggplot(seur_obj_filtered@meta.data, aes(x = nFeature_RNA, y = nFeature_SCT)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", lwd = 2) +
  labs(x = "Raw nGenes", y = "SCT nGenes") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p5<-ggplot(seur_obj_filtered@meta.data, aes(x = nCount_RNA, y = nCount_SCT)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", lwd = 2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Raw nUMI", y = "SCT nUMI") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

# Compare the relationship between raw and SCT normalized gene counts
p6<-ggplot(seur_obj_filtered@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point() +
  geom_point(aes(x = nCount_SCT, y = nFeature_SCT), color = "blue") +
  labs(x = "nUMI", y = "nGenes") +
  scale_x_log10() +
  scale_y_log10() +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

# Scale data in Seurat object
seur_obj_filtered <- ScaleData(seur_obj_filtered, assay = "RNA")

# Run PCA on the filtered genes and visualize
pca_genes_ensg <- names(top_expressed_genes)
pca_genes_symbol <- rowData(sce_filtered)$Symbol[rowData(sce_filtered)$ID %in% pca_genes_ensg]

seur_obj_filtered <- RunPCA(seur_obj_filtered, assay = "RNA", features = pca_genes_symbol)
p7<-FeaturePlot(seur_obj_filtered, reduction = "pca", features = "nCount_RNA") + theme_minimal()
p7<-p7+  theme(
  panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
  panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p7

# Run PCA with SCT normalized data
seur_obj_filtered <- RunPCA(seur_obj_filtered, assay = "SCT", features = pca_genes_symbol)
FeaturePlot(seur_obj_filtered, reduction = "pca", features = "nCount_RNA") + theme_minimal()

# Identify highly variable genes and visualize with ggplot2
seur_obj_filtered <- FindVariableFeatures(seur_obj_filtered)
top10 <- head(VariableFeatures(seur_obj_filtered), 10)
plot1 <- VariableFeaturePlot(seur_obj_filtered) + theme_minimal()
plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 <- plot1 +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

plot1

# Add a key marker and remove mitochondrial genes from the HVG list
hvgs <- VariableFeatures(seur_obj_filtered)
hvgs <- c(hvgs, "EOMES")  # Add key marker of liver NK cells
hvgs <- hvgs[!grepl("^MT-", hvgs)]  # Remove mitochondrial genes
VariableFeatures(seur_obj_filtered) <- hvgs

# Select HVGs using M3Drop
hvg_test_m3drop <- M3Drop::NBumiFeatureSelectionCombinedDrop(fit, qval.thresh = 0.05, suppress.plot = FALSE)
hvgs_m3drop <- hvg_test_m3drop$Gene

length(hvgs_m3drop)
#[1] 1001

# Select HVGs using scran and visualize the fit
hvg_model <- modelGeneVar(sce_filtered)
is.hvg <- hvg_model$p.value < 0.05
hvgs_scran <- rownames(hvg_model[is.hvg,])
length(hvgs_scran)

#[1] 2882

scran_fit <- metadata(hvg_model)

p8<-ggplot(data.frame(mean = scran_fit$mean, var = scran_fit$var), aes(x = mean, y = var)) +
  geom_point() +
  geom_smooth(method = "loess", color = "dodgerblue") +
  geom_point(data = data.frame(mean = scran_fit$mean[is.hvg], var = scran_fit$var[is.hvg]), color = "red") +
  labs(x = "Mean of log-expression", y = "Variance of log-expression") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p8

# Run PCA on SCT data and visualize with ggplot2
seur_obj_filtered <- RunPCA(seur_obj_filtered, assay = "SCT")
p9<-DimPlot(seur_obj_filtered, reduction = "pca") + theme_minimal()
p9 <- p9 + theme(panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p9

# Run PCA on SCE data using scran HVGs and visualize with ggplot2
sce_filtered <- runPCA(sce_filtered, exprs_values = "logcounts", subset_row = hvgs_scran)
pca_data <- plotPCA(sce_filtered, ncomponents = 1:2)

p10<-ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

# Run t-SNE and visualize
seur_obj_filtered <- RunTSNE(seur_obj_filtered, dims = 1:20, perplexity = 50)
p11<-DimPlot(seur_obj_filtered, reduction = "tsne") + 
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p11


sce_filtered <- scater::runTSNE(sce_filtered, perplexity = 50, dimred = "PCA", n_dimred = 20)
plot_tsne <- plotReducedDim(sce_filtered, dimred = "TSNE")+ 
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p12<-ggplot(plot_tsne, aes(x = TSNE1, y = TSNE2)) +
  geom_point() +
  labs(title = "t-SNE Plot", x = "t-SNE1", y = "t-SNE2") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

plot_tsne

# Run UMAP and visualize
seur_obj_filtered <- RunUMAP(seur_obj_filtered, dims = 1:20)
p13<-DimPlot(seur_obj_filtered, reduction = "umap") + 
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p13


sce_filtered <- scater::runUMAP(sce_filtered, dimred = "PCA", n_dimred = 20)
plot_umap <- plotReducedDim(sce_filtered, dimred = "UMAP")+ 
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

ggplot(plot_umap, aes(x = UMAP1, y = UMAP2)) +
  geom_point() +
  labs(title = "UMAP Plot", x = "UMAP1", y = "UMAP2") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

# View the reduced dimensions of the SCE object
reducedDims(sce_filtered)

# View the reductions in Seurat object
Seurat::Reductions(seur_obj_filtered)

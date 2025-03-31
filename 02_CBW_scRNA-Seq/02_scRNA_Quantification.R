library(Seurat)
library(ggplot2)
library(DropletUtils)
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

# Load the raw data using Seurat's Read10X function
raw_dat <- Seurat::Read10X(data.dir = "D:/projects/CBW_scRNA/scRNA_Data/HumanLiver/raw_feature_bc_matrix/")
dim(raw_dat)  # Display dimensions of the raw data matrix

#[1]  36601 407866

# Summary of total UMI counts per gene and per droplet
summary(rowSums(raw_dat))  # Total UMI per gene

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       1      24    1164     232 5799038 

summary(colSums(raw_dat))  # Total UMI per droplet

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       1      24    1164     232 5799038 

# Filter out droplets with fewer than 2 UMI counts
raw_dat <- raw_dat[, colSums(raw_dat) > 1]
dim(raw_dat)  # Display dimensions of the filtered data matrix

#[1]  36601 174598

# Create a barcode rank plot using ggplot2
br.out <- barcodeRanks(raw_dat)
br.out_df <- data.frame(rank = br.out$rank, total = br.out$total)

p1<-ggplot(br.out_df, aes(x = rank, y = total)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Rank", y = "Total") +
  geom_hline(yintercept = metadata(br.out)$knee, color = "dodgerblue", linetype = "dashed") +
  geom_hline(yintercept = metadata(br.out)$inflection, color = "forestgreen", linetype = "dashed") +
  theme(
  panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
  panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  annotate("text", x = max(br.out$rank), y = metadata(br.out)$knee, 
           label = "knee", color = "dodgerblue", vjust = -1 , hjust = 1) +
  annotate("text", x = max(br.out$rank), y = metadata(br.out)$inflection, 
           label = "inflection", color = "forestgreen", vjust = 2 , hjust = 1)

p1

# Identifying cells based on knee and inflection points
is.cell.knee <- br.out$total > metadata(br.out)$knee
sum(is.cell.knee)  # Number of cells above the knee point

#[1] 5331

is.cell.inflect <- br.out$total > metadata(br.out)$inflection
sum(is.cell.inflect)  # Number of cells above the inflection point

#[1] 6283

# Run emptyDrops to detect empty droplets
e.out <- emptyDrops(raw_dat, lower = 100, niters = 10000, ignore = NULL, retain = 2 * br.out$knee)
head(e.out)  # Display the first few rows of the emptyDrops output

# Total   LogProb    PValue   Limited       FDR
# <integer> <numeric> <numeric> <logical> <numeric>
#   AAACCTGAGAAACGAG-1       309  -510.136  0.973703     FALSE         1
# AAACCTGAGAAGGCCT-1         3        NA        NA        NA        NA
# AAACCTGAGACAAAGG-1       331  -524.053  0.989701     FALSE         1
# AAACCTGAGACATAAC-1       284  -552.614  0.436556     FALSE         1
# AAACCTGAGACCCACC-1         2        NA        NA        NA        NA
# AAACCTGAGACCGGAT-1       254  -387.596  0.999700     FALSE         1

# Filter out NA values in the PValue column
e.out <- e.out[!is.na(e.out$PValue), ]
head(e.out[e.out$Limited, ])  # Display droplets with limited p-values

# Total   LogProb    PValue   Limited        FDR
# <integer> <numeric> <numeric> <logical>  <numeric>
#   AAACCTGAGCACAGGT-1       252  -678.886 9.999e-05      TRUE 0.00177417
# AAACCTGAGCCAGTAG-1       254  -726.443 9.999e-05      TRUE 0.00177417
# AAACCTGAGCGTTCCG-1      3092 -3373.245 9.999e-05      TRUE 0.00177417
# AAACCTGAGCTGAAAT-1      2233 -5462.634 9.999e-05      TRUE 0.00177417
# AAACCTGAGTACGTTC-1      4805 -9909.250 9.999e-05      TRUE 0.00177417
# AAACCTGCAGTTCATG-1      2662 -5979.771 9.999e-05      TRUE 0.00177417


# Identify cells based on FDR threshold
is.cell <- e.out$FDR <= 0.01

# Plot total UMI count vs -log probability using ggplot2
e.out_df <- data.frame(Total = e.out$Total, LogProb = -e.out$LogProb, IsCell = is.cell)


p2<-ggplot(e.out_df, aes(x = Total, y = LogProb, color = IsCell)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "Total UMI count", y = "-Log Probability") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p2

# Get cell barcodes identified as cells
e.cells <- rownames(e.out)[is.cell]

# Create a barcode rank plot with points colored by cell identification
br.out_df$IsCell <- rownames(br.out) %in% e.cells
p3<-ggplot(br.out_df, aes(x = rank, y = total, color = IsCell)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Rank", y = "Total") +
  geom_hline(yintercept = metadata(br.out)$knee, color = "dodgerblue", linetype = "dashed") +
  geom_hline(yintercept = metadata(br.out)$inflection, color = "forestgreen", linetype = "dashed") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) +
  scale_color_manual(values = c("black", "red")) +
  annotate("text", x = max(br.out$rank), y = metadata(br.out)$knee, 
           label = "knee", color = "dodgerblue", vjust = -1 , hjust = 1) +
  annotate("text", x = max(br.out$rank), y = metadata(br.out)$inflection, 
           label = "inflection", color = "forestgreen", vjust = 2 , hjust = 1)

p3

# Further analysis: points for identified cells (red) and non-cells (blue)
br.out_df$IsCellType <- ifelse(br.out_df$IsCell, "Cell", "Non-Cell")
p4<-ggplot(br.out_df, aes(x = rank, y = total, color = IsCellType)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Rank", y = "Total") +
  geom_hline(yintercept = metadata(br.out)$knee, color = "dodgerblue", linetype = "dashed") +
  geom_hline(yintercept = metadata(br.out)$inflection, color = "forestgreen", linetype = "dashed") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) +
  scale_color_manual(values = c("blue", "red"))

p4

# Create Seurat object with the filtered data
filter_dat <- Seurat::Read10X("D:/projects/CBW_scRNA/scRNA_Data/HumanLiver/filtered_feature_bc_matrix/")
seur_obj <- Seurat::CreateSeuratObject(filter_dat, min.cells = 5, min.features = 100)

# Normalize, find variable features, and scale the data
seur_obj <- Seurat::NormalizeData(seur_obj)
seur_obj <- Seurat::FindVariableFeatures(seur_obj)
seur_obj <- Seurat::ScaleData(seur_obj)

# Run PCA, find neighbors, and cluster the data
seur_obj <- Seurat::RunPCA(seur_obj)
seur_obj <- Seurat::FindNeighbors(seur_obj)
seur_obj <- Seurat::FindClusters(seur_obj)

# Display metadata and cluster frequencies
head(seur_obj@meta.data)

# orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 seurat_clusters
# AAACCTGAGCGTTCCG-1 SeuratProject       3092          782               7               7
# AAACCTGAGCTGAAAT-1 SeuratProject       2233         1142               2               2
# AAACCTGAGGGCTTGA-1 SeuratProject       2467          156               0               0
# AAACCTGAGTACGTTC-1 SeuratProject       4803         2128               2               2
# AAACCTGAGTGCCATT-1 SeuratProject       3899          168               1               1
# AAACCTGCAATGGAAT-1 SeuratProject       3443          149               4               4

table(seur_obj@meta.data$seurat_clusters)

# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 1357 1161 1037  714  602  542  301  296  276  249  229  175  140  110   83   47

# Calculate expected number of doublets
type.freq <- table(seur_obj@meta.data$seurat_clusters) / ncol(seur_obj)
homotypic.prop <- sum(type.freq^2)
nEXP <- 0.009 * (ncol(seur_obj) / 1000) * (1 - homotypic.prop) * ncol(seur_obj)
nEXP

# [1] 429.2712

# Sweep parameters and plot results using ggplot2
pN <- 0.25
PCs <- 18
sweep.out <- paramSweep(seur_obj, PCs = 1:PCs)
sweep.stats <- summarizeSweep(sweep.out)
sweep.stats_df <- as.data.frame(sweep.stats)

p5<-ggplot(sweep.stats_df, aes(x = as.numeric(as.character(sweep.stats_df[,2])), y = sweep.stats_df[,3])) +
  geom_point() +
  labs(x = "pN", y = "Sweep Statistic") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p5

# Filter and identify doublets
pK <- 0.02
pN <- 0.2
seur_obj <- doubletFinder(seur_obj, PCs = 1:PCs, pN = pN, pK = pK, nExp = nEXP)
head(seur_obj@meta.data)

# orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 seurat_clusters pANN_0.2_0.02_429.2712
# AAACCTGAGCGTTCCG-1 SeuratProject       3092          782               7               7             0.28961749
# AAACCTGAGCTGAAAT-1 SeuratProject       2233         1142               2               2             0.08743169
# AAACCTGAGGGCTTGA-1 SeuratProject       2467          156               0               0             0.15846995
# AAACCTGAGTACGTTC-1 SeuratProject       4803         2128               2               2             0.16393443
# AAACCTGAGTGCCATT-1 SeuratProject       3899          168               1               1             0.06557377
# AAACCTGCAATGGAAT-1 SeuratProject       3443          149               4               4             0.03825137
# DF.classifications_0.2_0.02_429.2712
# Singlet
# Singlet
# Singlet
# Singlet
# Singlet
# Singlet

# UMAP plot colored by doublet classification using ggplot2
seur_obj <- RunUMAP(seur_obj, dims = 1:18)
p6<-DimPlot(seur_obj, reduction = "umap")

p6<-p6+  theme(
  panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
  panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p6

# Bar plot showing the proportion of doublets in each cluster
tab <- table(seur_obj@meta.data$seurat_clusters, seur_obj@meta.data$DF.classifications_0.2_0.02_429.2712)
tab_df <- as.data.frame.matrix(tab / rowSums(tab))
tab_df$Cluster <- rownames(tab_df)

p7<-ggplot(melt(tab_df, id.vars = "Cluster"), aes(x = Cluster, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(y = "% of cells", x = "Cluster") +
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#dbe3db"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#dbe3db"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p7

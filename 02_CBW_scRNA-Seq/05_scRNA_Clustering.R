require(cluster)

# Load the Seurat object
seur_obj <- readRDS("Liver_QC_Norm_DimReduc_Seurat.rds")

# Find neighbors using k-nearest neighbors (k = 20)
seur_obj <- FindNeighbors(seur_obj, k.param = 20)

# View the first few entries of the SCT_snn graph
head(seur_obj@graphs$SCT_snn[1:5,1:25])

# Perform clustering at resolution 0.8
seur_obj <- FindClusters(seur_obj, resolution=0.8)

# Plot the clusters using ggplot2-based DimPlot
p1<-DimPlot(seur_obj, label=TRUE) + ggtitle("Clusters at Resolution 0.8")+
    theme(
      panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p1

# View the metadata with cluster assignments

# Check the number of cells in each cluster
table(seur_obj@meta.data$SCT_snn_res.0.8)

# 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14 
# 444 426 415 217 184 162 155 135 105  93  74  53  50  42  38 

# Perform clustering at a lower resolution (0.5)
seur_obj <- FindClusters(seur_obj, resolution=0.5)
head(seur_obj@meta.data)

# orig.ident nCount_RNA nFeature_RNA percent.mt pass.filter nCount_SCT nFeature_SCT SCT_snn_res.0.8 seurat_clusters SCT_snn_res.0.5
# AAACCTGAGCTGAAAT-1 SeuratProject       2233         1142   8.150470        TRUE       2216         1124               8               7               7
# AAACCTGAGTACGTTC-1 SeuratProject       4805         2130   7.159209        TRUE       3011         1957               2               1               1
# AAACCTGCAGTTCATG-1 SeuratProject       2662         1277  11.532682        TRUE       2447         1254               1               2               2
# AAACGGGAGATCCGAG-1 SeuratProject       3697         1778   6.681093        TRUE       2724         1706               4               3               3
# AAACGGGCAGCGTTCG-1 SeuratProject       1505          818  13.421927        TRUE       1894          795               6               5               5
# AAACGGGGTCAACTGT-1 SeuratProject       1762          953  10.329171        TRUE       1983          930               1               2               2

# Plot the clusters using ggplot2-based DimPlot
p2<-DimPlot(seur_obj) + ggtitle("Clusters at Resolution 0.5")+
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p2

# Assign the clusters at resolution 0.8 to a new column in metadata
seur_obj@meta.data$SCT_k20_res0.8 <- seur_obj@meta.data$SCT_snn_res.0.8

# Compute a distance matrix based on the PCA embeddings
dist.matrix <- dist(x = Embeddings(object = seur_obj[["pca"]])[, 1:20])

# Compute silhouette scores for the clusters
sil <- silhouette(as.numeric(as.character(seur_obj@meta.data$SCT_k20_res0.8)), dist=dist.matrix)
head(sil)

# cluster neighbor sil_width
# [1,]       8        0 0.4023732
# [2,]       2        8 0.2446068
# [3,]       1        8 0.4301426
# [4,]       4        7 0.2666117
# [5,]       6        0 0.5633484
# [6,]       1        8 0.4318904

# Calculate the mean silhouette width across all clusters
mean_silhouette <- mean(sil[,3])
mean_silhouette

#[1] 0.3042467

# Aggregate the silhouette widths by cluster
cluster_sil_scores <- aggregate(sil[,3], by=list(sil[,1]), mean)
cluster_sil_scores

# Group.1          x
# 1        0 0.06379908
# 2        1 0.28495395
# 3        2 0.20591677
# 4        3 0.54741813
# 5        4 0.26801986
# 6        5 0.24464399
# 7        6 0.57983728
# 8        7 0.35935293
# 9        8 0.33283824
# 10       9 0.29882867
# 11      10 0.50525859
# 12      11 0.54886500
# 13      12 0.47905592
# 14      13 0.78960128
# 15      14 0.55999610

# Compute Adjusted Rand Index (ARI) to compare clustering solutions
require(igraph) 
ARI <- compare(seur_obj@meta.data$SCT_snn_res.0.8, seur_obj@meta.data$SCT_snn_res.0.5, method="adjusted.rand")
ARI

# [1] 0.8547882

# Identify marker genes for each cluster at resolution 0.8
markers <- FindAllMarkers(seur_obj, group.by="SCT_snn_res.0.8", logfc.threshold = -Inf, only.pos=TRUE, max.cells.per.ident=100)
markers <- markers[ markers[,"p_val_adj"] < 0.05, ]
table(markers$cluster)

# 0   1   2   3   4   5   6   7   8   9  10  11  12 
# 18  53  48 368 376 132 107 379  41 283 376 172 362 

# View the top 20 markers for cluster 0
head(markers[markers$cluster == "0",], 10)

# p_val avg_log2FC pct.1 pct.2    p_val_adj cluster   gene
# RPS29  3.458152e-12  0.8374020 1.000 0.969 6.821550e-08       0  RPS29
# RPS27  5.408515e-12  0.6707508 1.000 0.980 1.066884e-07       0  RPS27
# IL32   9.747951e-12  1.7365211 0.827 0.298 1.922881e-07       0   IL32
# CD3D   3.965245e-10  2.3655857 0.606 0.121 7.821842e-06       0   CD3D
# TRAC   5.557934e-10  2.0305163 0.489 0.118 1.096358e-05       0   TRAC
# CCL5   4.572782e-09  1.3730982 0.863 0.404 9.020270e-05       0   CCL5
# CRTAM  2.065110e-08  3.9743291 0.273 0.035 4.073636e-04       0  CRTAM
# TRBC2  3.573150e-08  1.5093004 0.475 0.215 7.048397e-04       0  TRBC2
# INPP4B 1.226146e-07  2.8814278 0.473 0.084 2.418696e-03       0 INPP4B
# CD2    2.094090e-07  1.9427453 0.505 0.162 4.130802e-03       0    CD2

# Define cluster annotations
nclusters <- length(unique(seur_obj@meta.data$SCT_snn_res.0.8))
cluster_annotation <- rep("unannotated", nclusters)
names(cluster_annotation) <- levels(seur_obj@meta.data$SCT_snn_res.0.8)
cluster_annotation["0"] <- "T cell"

# Define genes used for annotation
anno_genes <- list("cluster0"=c("CD3D","CD3G", "TRAC", "CCL5", "GZMK", "RUNX3"))

# Set cluster annotations in the Seurat object
seur_obj = SetIdent(seur_obj, value=seur_obj@meta.data$SCT_snn_res.0.8)
seur_obj = RenameIdents(seur_obj, cluster_annotation)
seur_obj@meta.data$manual_anno <- Idents(seur_obj)
head(seur_obj@meta.data)

# orig.ident nCount_RNA nFeature_RNA percent.mt pass.filter nCount_SCT nFeature_SCT SCT_snn_res.0.8 seurat_clusters SCT_snn_res.0.5 SCT_k20_res0.8 manual_anno
# AAACCTGAGCTGAAAT-1 SeuratProject       2233         1142   8.150470        TRUE       2216         1124               8               7               7              8 unannotated
# AAACCTGAGTACGTTC-1 SeuratProject       4805         2130   7.159209        TRUE       3011         1957               2               1               1              2 unannotated
# AAACCTGCAGTTCATG-1 SeuratProject       2662         1277  11.532682        TRUE       2447         1254               1               2               2              1 unannotated
# AAACGGGAGATCCGAG-1 SeuratProject       3697         1778   6.681093        TRUE       2724         1706               4               3               3              4 unannotated
# AAACGGGCAGCGTTCG-1 SeuratProject       1505          818  13.421927        TRUE       1894          795               6               5               5              6 unannotated
# AAACGGGGTCAACTGT-1 SeuratProject       1762          953  10.329171        TRUE       1983          930               1               2               2              1 unannotated

# Plot a dot plot for the annotation genes

p3<-DotPlot(seur_obj, features=unique(unlist(anno_genes))) + ggtitle("DotPlot for Annotation Genes")+
  theme(
    panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p3

# Load the reference dataset
ref <- readRDS("MacParlandLiverAtlas_celltype_profiles.rds")
head(ref)

# CentralvenousLSECs Cholangiocytes Non-inflammatoryMacrophages CD3abTcells InflamatoryMacrophages NK-likecells gdTcells1 PeriportalHep Portalendothelialcells gdTcells2
# A1BG                 0.000000              0                           0    0.000000                      0     0.000000         0      2.146611              0.0000000 0.0000000
# A2M                  2.201942              0                           0    0.000000                      0     0.000000         0      2.514386              0.0000000 0.0000000
# AADAC                0.000000              0                           0    0.000000                      0     0.000000         0      0.647771              0.0000000 0.0000000
# AC092580.4           0.000000              0                           0    0.000000                      0     0.000000         0      0.000000              0.0000000 0.0000000
# ACAP1                0.000000              0                           0    1.319308                      0     1.578042         0      0.000000              0.0000000 0.7470038
# ACP5                 2.315153              0                           0    0.000000                      0     0.000000         0      0.000000              0.9899168 0.0000000
# PeriportalLSECs interzonalHep MatureBcells Stellatecells AntibodysecretingBcells Erythroidcells PericentralHep UnidentifiedHep
# A1BG              0.000000             0            0             0                       0              0       1.890293        3.168597
# A2M               1.222757             0            0             0                       0              0       0.652065        1.200735
# AADAC             0.000000             0            0             0                       0              0       1.379483        1.742586
# AC092580.4        0.000000             0            0             0                       0              0       0.000000        0.000000
# ACAP1             0.000000             0            0             0                       0              0       0.000000        0.000000
# ACP5              0.000000             0            0             0                       0              0       0.000000        0.000000

# Calculate average expression in the query dataset
query = AverageExpression(seur_obj, group.by = "SCT_snn_res.0.8", assays="RNA")
head(query$RNA)

# MIR1302-2HG .          . . . . . . .          . . . . . . .
# FAM138A     .          . . . . . . .          . . . . . . .
# OR4F5       .          . . . . . . .          . . . . . . .
# AL627309.1  0.01002337 . . . . . . 0.01682355 . . . . . . .
# AL627309.3  .          . . . . . . .          . . . . . . .
# AL627309.2  .          . . . . . . .          . . . . . . .

# Filter the variable genes to match with the reference dataset
genes = VariableFeatures(seur_obj)
genes = genes[genes %in% rownames(ref)]

# Match genes between query and reference
query2 = query$RNA[match(genes, rownames(query$RNA)),]
ref2 = ref[match(genes, rownames(ref)),]
identical(rownames(ref2), rownames(query2))

# Compute correlations between query and reference
cross_cors = cor(as.matrix(query2), as.matrix(ref2))

# Plot the correlation heatmap using ggplot2
library(reshape2)
cross_cors_melted <- melt(cross_cors)
p4<-ggplot(data = cross_cors_melted, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.background = element_rect(fill = "#f2f2f2", colour = "#f2f2f2", size = 0.5, linetype = "solid"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) +
  ggtitle("Correlation Heatmap between Query and Reference")

p4

# Further annotate clusters based on correlation with reference
cluster_annotation2 <- rep("unannotated", nclusters)
names(cluster_annotation2) <- levels(seur_obj@meta.data$SCT_snn_res.0.8)
cluster_annotation2["0"] = "T cell"
cluster_annotation2["1"] = "T cell"
cluster_annotation2["2"] = "T cell"
cluster_annotation2["3"] = "InfMac"
cluster_annotation2["4"] = "NK cell"
cluster_annotation2["5"] = "gdT cell"
cluster_annotation2["6"] = "CD3T cell"
cluster_annotation2["7"] = "B/T cell"
cluster_annotation2["8"] = "B cell"
cluster_annotation2["9"] = "Endothelial"
cluster_annotation2["10"] = "Macrophage"
cluster_annotation2["11"] = "Hepatocyte"
cluster_annotation2["12"] = "AntiB"
cluster_annotation2["13"] = "cvLSEC"
cluster_annotation2["14"] = "Erythroid"
cluster_annotation2["15"] = "Cholangio"

# Set cluster annotations in the Seurat object based on reference
seur_obj = SetIdent(seur_obj, value=seur_obj@meta.data$SCT_snn_res.0.8)
seur_obj = RenameIdents(seur_obj, cluster_annotation2)
seur_obj@meta.data$ref_anno <- Idents(seur_obj)
head(seur_obj@meta.data)

# orig.ident nCount_RNA nFeature_RNA percent.mt pass.filter nCount_SCT nFeature_SCT SCT_snn_res.0.8 seurat_clusters
# AAACCTGAGCTGAAAT-1 SeuratProject       2233         1142   8.150470        TRUE       2216         1124               8               7
# AAACCTGAGTACGTTC-1 SeuratProject       4805         2130   7.159209        TRUE       3011         1957               2               1
# AAACCTGCAGTTCATG-1 SeuratProject       2662         1277  11.532682        TRUE       2447         1254               1               2
# AAACGGGAGATCCGAG-1 SeuratProject       3697         1778   6.681093        TRUE       2724         1706               4               3
# AAACGGGCAGCGTTCG-1 SeuratProject       1505          818  13.421927        TRUE       1894          795               6               5
# AAACGGGGTCAACTGT-1 SeuratProject       1762          953  10.329171        TRUE       1983          930               1               2
# SCT_snn_res.0.5 SCT_k20_res0.8 manual_anno  ref_anno
# AAACCTGAGCTGAAAT-1               7              8 unannotated    B cell
# AAACCTGAGTACGTTC-1               1              2 unannotated    T cell
# AAACCTGCAGTTCATG-1               2              1 unannotated    T cell
# AAACGGGAGATCCGAG-1               3              4 unannotated   NK cell
# AAACGGGCAGCGTTCG-1               5              6 unannotated CD3T cell
# AAACGGGGTCAACTGT-1               2              1 unannotated    T cell

# Load a mouse brain dataset
mouse_brain_obj <- readRDS("MouseBrain_with_Replicates_Annotated.rds")

# Load custom R functions
source("./My_R_Functions.R")

# View metadata of the mouse brain dataset
head(mouse_brain_obj@meta.data)

# orig.ident nCount_RNA nFeature_RNA       cell_barcode num_features feature_call num_umis  Mouse nCount_SCT
# AAACCCACAAGCACCC-1 SeuratProject        174           69 AAACCCACAAGCACCC-1            1       CMO301    10306 Mouse1         77
# AAACGAAAGGGCTTCC-1 SeuratProject         49           29 AAACGAAAGGGCTTCC-1            1       CMO301    17115 Mouse1         49
# AAAGAACAGGGACTGT-1 SeuratProject         34           20 AAAGAACAGGGACTGT-1            1       CMO301    32528 Mouse1         45
# AAAGAACGTTGGCCTG-1 SeuratProject         36           23 AAAGAACGTTGGCCTG-1            1       CMO301    42172 Mouse1         42
# AAAGGGCTCGCGCCAA-1 SeuratProject         34           23 AAAGGGCTCGCGCCAA-1            1       CMO301    31014 Mouse1         40
# AAAGTCCAGCCTCTGG-1 SeuratProject         35           22 AAAGTCCAGCCTCTGG-1            1       CMO301    29205 Mouse1         45
# nFeature_SCT SCT_snn_res.0.8 seurat_clusters SCT_snn_res.0.2 SCT_snn_res.0.1         cell_type
# AAACCCACAAGCACCC-1           51               3               0               0               0 Excitatory Neuron
# AAACGAAAGGGCTTCC-1           29               0               0               0               0 Excitatory Neuron
# AAAGAACAGGGACTGT-1           20               0               0               0               0 Excitatory Neuron
# AAAGAACGTTGGCCTG-1           23              11               4               5               4   Cortical neuron
# AAAGGGCTCGCGCCAA-1           21              12               5               6               5          Pericyte
# AAAGTCCAGCCTCTGG-1           22               1               0               2               0 Excitatory Neuron

table(mouse_brain_obj@meta.data$cell_type)

# Excitatory Neuron Inhibitory Neuron        Astrocytes     Prukinje cell   Cortical neuron          Pericyte       Interneuron   Red Blood Cells 
# 10296              2623              1973               913               802               764               282               205 

# Subset the mouse brain dataset to a specific cell type
type = "Astrocytes"
subset_obj <- mouse_brain_obj[,mouse_brain_obj@meta.data$cell_type==type]

# Create a table of feature calls by mouse
table( subset_obj@meta.data$feature_call, subset_obj@meta.data$Mouse)

# Mouse1 Mouse2 Mouse3 Mouse4
# CMO301     70      0      0      0
# CMO302     68      0      0      0
# CMO303     99      0      0      0
# CMO304      0    193      0      0
# CMO305      0    123      0      0
# CMO306      0    156      0      0
# CMO307      0      0    258      0
# CMO308      0      0    218      0
# CMO309      0      0    249      0
# CMO310      0      0      0    246
# CMO311      0      0      0    197
# CMO312      0      0      0     96

# Create pseudobulks by summing RNA counts for each feature call
pseudobulks <- group_rowmeans(subset_obj@assays$RNA@counts, subset_obj@meta.data$feature_call, type="sum")
head(pseudobulks)

# CMO301 CMO302 CMO303 CMO304 CMO305 CMO306 CMO307 CMO308 CMO309 CMO310 CMO311 CMO312
# Xkr4         0      0      0      0      0      0      0      0      0      0      0      0
# Gm1992       0      0      0      0      0      0      0      0      0      0      0      0
# Gm19938      0      0      0      0      0      1      0      0      0      0      0      0
# Gm37381      0      0      0      0      0      0      0      0      0      0      0      0
# Rp1          0      0      0      0      0      0      0      0      0      0      0      0
# Sox17        0      0      0      1      0      0      1      0      2      1      0      0

# Filter the pseudobulks to remove low-count features
pseudobulks_filt <- pseudobulks[rowSums(pseudobulks > 0) > 3,]

# Create a table of feature calls by mouse
table( subset_obj@meta.data$feature_call, subset_obj@meta.data$Mouse)
colnames(pseudobulks_filt)
bio_condition <- c("Mouse1","Mouse1","Mouse1",
                   "Mouse2","Mouse2","Mouse2",
                   "Mouse3","Mouse3","Mouse3",
                   "Mouse4","Mouse4","Mouse4")

# Create the design matrix for differential expression analysis
design <- model.matrix(~bio_condition)
head(design)

# (Intercept) bio_conditionMouse2 bio_conditionMouse3 bio_conditionMouse4
# 1           1                   0                   0                   0
# 2           1                   0                   0                   0
# 3           1                   0                   0                   0
# 4           1                   1                   0                   0
# 5           1                   1                   0                   0
# 6           1                   1                   0                   0

# Perform differential expression analysis using edgeR
pseudobulk_dge <- DGEList(counts=pseudobulks_filt, group=bio_condition)
pseudobulk_dge <- calcNormFactors(pseudobulk_dge)
pseudobulk_dge <- estimateDisp(pseudobulk_dge, design)
fit <- glmQLFit(pseudobulk_dge)

# Perform the differential expression test for the second coefficient
res_2vs1 <- glmQLFTest(fit, coef=2)

# View the top 20 results
DEG_df<-data.frame(topTags(res_2vs1, 20))

# logFC    logCPM         F       PValue          FDR
# Npy      4.8122799 11.757091 141.37326 1.566787e-32 6.130839e-29
# Plp1     2.1642033 11.705802  43.68663 3.915801e-11 7.661265e-08
# Chodl    3.1287361  9.906157  33.85454 5.996809e-09 7.821838e-06
# Nts      3.4903373  9.324251  23.59695 1.193490e-06 1.167532e-03
# Ptma     1.6834373 10.978359  20.48722 6.025173e-06 4.715300e-03
# Xist    -3.1416328  9.047962  17.08945 3.575272e-05 2.331673e-02
# Cdk4     5.8365538  8.255458  15.93084 6.584895e-05 3.680956e-02
# Syt6     3.0163166  8.731337  14.47326 1.424338e-04 6.702696e-02
# Dcn     -1.2742017 10.572684  14.32416 1.541637e-04 6.702696e-02

DEG_df$log_padj<- -log(as.numeric(DEG_df$FDR))
DEG_df$names<-rownames(DEG_df)

p <- ggplot(DEG_df, aes(x=as.numeric(as.character(logFC)), y=as.numeric(as.character(log_padj)))) +
  geom_point(aes(color = ifelse(log_padj > 3, "#4d4d4d", "#b3b3b3")), alpha = 1.2) +
  scale_color_identity() +
  labs(title="") +
  geom_hline(yintercept=3, linetype="dashed", color = "black") +
  geom_vline(xintercept=0, linetype="dashed", color = "black") +
  xlab("Log2 Fold Change") +
  ylab("-log(FDR adjusted p-value)") +
  geom_text(data = subset(DEG_df, log_padj > 3), aes(label = names), 
            vjust = -0.5, size = 3.5, color = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

p

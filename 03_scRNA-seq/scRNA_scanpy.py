# 1. Load Libraries
# --------------------------
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Settings for plots
sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

print(f"Scanpy version: {sc.__version__}")

# Load the PBMC3k dataset (3k Peripheral Blood Mononuclear Cells)
# This function downloads the data if not already present
# It returns an AnnData object, which is the core data structure in Scanpy
adata = sc.datasets.pbmc3k()

print("Initial AnnData object:")
print(adata)
# AnnData object with n_obs × n_vars = 2700 × 32738
#     var: 'gene_ids'

# 2. Basic Preprocessing and Quality Control (QC)
# ----------------------------------------------
# Show genes with the highest counts before filtering
sc.pl.highest_expr_genes(adata, n_top=20, save='_highest_expr_before_qc.png')

# Basic filtering:
# - Filter out cells with few genes expressed (potential empty droplets or debris)
# - Filter out cells with very high gene counts (potential doublets)
# - Filter out genes expressed in too few cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print("AnnData object after basic filtering:")
print(adata)
# AnnData object with n_obs × n_vars = 2638 × 13714
#     obs: 'n_genes'
#     var: 'n_cells', 'gene_ids'

# Calculate QC metrics, especially % mitochondrial genes
# Mitochondrial genes often indicate cell stress or apoptosis if percentage is high
adata.var['mt'] = adata.var_names.str.startswith('MT-') # Annotate mitochondrial genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Visualize QC metrics
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save='_qc_violin.png')

# Filter cells based on QC metrics
# These thresholds are dataset-specific and often require inspection of the QC plots
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

print("AnnData object after QC filtering:")
print(adata)
# AnnData object with n_obs × n_vars = 2638 × 13714
#     obs: 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt'
#     var: 'gene_ids', 'n_cells', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts'


# 3. Normalization and Feature Selection
# --------------------------------------
# Normalize counts per cell to a total count (e.g., 10,000 reads per cell)
# This corrects for differences in sequencing depth between cells
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform the data
# This stabilizes variance and makes distributions more symmetrical
sc.pp.log1p(adata)

# Identify highly variable genes (HVGs)
# These genes show high biological variability and are informative for downstream analysis
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, save='_hvg.png')

# Subset the AnnData object to keep only HVGs (optional but common for performance)
adata.raw = adata # Store the full data in .raw for later marker gene analysis
adata = adata[:, adata.var.highly_variable]

# Scale the data to unit variance and zero mean
# Also, optionally regress out unwanted sources of variation like total counts or MT%
# Clipping values avoids extreme outliers having too much influence
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)


# 4. Dimensionality Reduction
# ---------------------------
# Principal Component Analysis (PCA)
# Reduces dimensionality while retaining most variance
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, save='_pca_variance.png')
# Visualize PCA
sc.pl.pca(adata, color='total_counts', save='_pca_total_counts.png') # Example coloring

# Compute neighborhood graph (needed for UMAP and clustering)
# Uses PCA results to find nearest neighbors in high-dimensional space
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40) # Use ~40 PCs based on variance plot

# Uniform Manifold Approximation and Projection (UMAP)
# Further reduces dimensionality for visualization (typically to 2D)
sc.tl.umap(adata)
# Visualize UMAP
sc.pl.umap(adata, color=['total_counts', 'pct_counts_mt'], save='_umap_qc.png')


# 5. Clustering
# ---------------
# Cluster cells using the Leiden algorithm (community detection on the neighborhood graph)
sc.tl.leiden(adata, resolution=0.5) # Resolution parameter influences cluster number/granularity
# Visualize clusters on UMAP
sc.pl.umap(adata, color=['leiden'], legend_loc='on data', title='Leiden Clusters', save='_umap_leiden.png')


# 6. Finding Marker Genes
# -----------------------
# Identify genes differentially expressed in each cluster compared to all other cells
# This helps in annotating the cell types
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon') # Wilcoxon rank-sum test is common

# Visualize top marker genes
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_marker_genes.png')

# Show marker genes as a heatmap
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, groupby="leiden", show_gene_labels=True, save='_marker_heatmap.png')

# Show expression of specific known markers on UMAP
# (Example markers for PBMCs; requires prior biological knowledge)
marker_genes = ['IL7R', 'CD14', 'LYZ', 'MS4A1', 'CD8A', 'FCGR3A', 'MS4A7', 'GNLY', 'NKG7', 'FCER1A', 'CST3', 'PPBP']
# Note: Some might have been filtered out or not be highly variable
# Filter for markers present in the scaled data:
marker_genes_present = [gene for gene in marker_genes if gene in adata.var_names]

if marker_genes_present:
    sc.pl.umap(adata, color=marker_genes_present, save='_umap_known_markers.png')
    sc.pl.dotplot(adata, marker_genes_present, groupby='leiden', save='_dotplot_known_markers.png')
    sc.pl.violin(adata, marker_genes_present, groupby='leiden', rotation=90, save='_violin_known_markers.png')
else:
    print("None of the example marker genes are present in the HVG subset.")

# Accessing results:
print("\nCluster sizes:")
print(adata.obs['leiden'].value_counts())

print("\nTop marker genes for cluster 0:")
# Convert structured array to DataFrame for easier viewing
marker_df = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
print(marker_df['0'].head())


# 7. Annotation (Conceptual)
# --------------------------
# Compare the marker genes found for each cluster
# (e.g., using sc.pl.rank_genes_groups plots and the table in adata.uns['rank_genes_groups'])
# with known cell type markers from literature or databases to assign biological labels
# (e.g., "CD4 T cells", "B cells", "Monocytes") to the Leiden clusters.

# Manually creating an annotation dictionary based on observed markers
# (This is illustrative - actual annotation requires careful biological interpretation)
# cluster_annotation = {
#     '0': 'CD4 T cells',
#     '1': 'CD14+ Monocytes',
#     '2': 'B cells',
#     '3': 'CD8 T cells',
#     '4': 'NK cells',
#     '5': 'FCGR3A+ Monocytes',
#     '6': 'Dendritic Cells',
#     '7': 'Megakaryocytes'
# }
# adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotation).astype('category')
# sc.pl.umap(adata, color='cell_type', legend_loc='on data', title='Annotated Cell Types', save='_umap_annotated.png')


# 8. Saving Results
# -----------------
# Save the final AnnData object (contains all data, metadata, and analysis results)
# results_file = 'pbmc3k_processed.h5ad'
# adata.write(results_file)
# print(f"Analysis complete. Results saved to {results_file}")

# To load later:
# adata_loaded = sc.read_h5ad(results_file)

print("\nBasic Scanpy analysis workflow completed.")

# 1. Install and Load Libraries
# ---------------------------------------------------

# List of packages to check and install
packages <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
              'SummarizedExperiment', 'batchelor', 'Matrix.utils','HDF5Array', 'terra', 'ggrastr', 'devtools', 'dplyr', 'ggplot2')

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

devtools::install_github('cole-trapnell-lab/monocle3')

# Install cicero for the dataset
devtools::install_github('cole-trapnell-lab/cicero-release', ref = 'monocle3')

library(monocle3)
library(cicero) # For loading the example dataset

# Set a seed
set.seed(1234)

print(paste("Monocle 3 version:", packageVersion("monocle3")))

# 2. Load Data and Create Monocle 3 Object
# ---------------------------------------------------------

# Load the human skin fibroblast reprogramming dataset from cicero package
# Loads expression_matrix, cell_metadata, and gene_annotation
data(human_skin_fibroblasts)

# Create the cell_data_set (CDS) object - Monocle 3's main data structure
cds <- new_cell_data_set(
  expression_matrix, # The counts matrix
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

print("Initial cell_data_set object:")
# Print summary, show dimensions (genes x cells) and available metadata
print(cds)
# A cell_data_set object from Monocle 3
#   with 59702 genes and 792 cells
# Cell Metadata:
#              Time Day State Size_Factor UMI ...
# GA मस..ATC   0   0    FB          ...  ... ...
# GA मस..TTC   0   0    FB          ...  ... ...
# Gene Metadata:
#                  gene_short_name num_cells_expressed use_for_ordering ...
# ENSG00000241860           LINC0..                 10          FALSE  ...
# ENSG00000279928            AL62..                 1           FALSE  ...


# 3. Preprocessing (Normalization, Feature Selection, PCA)
# --------------------------------------------------------

# Monocle 3 provides a wrapper function for standard preprocessing steps
# Normalizes, log-transforms, selects features, and performs PCA
cds <- preprocess_cds(cds, num_dim = 100) # num_dim = number of PCs

# Visualize the variance explained by PC elbow plot to potentially adjust num_dim
plot_pc_variance_explained(cds)
ggsave("monocle3_pc_variance.png")

# Note: If batch correction is needed (e.g., different experiments),
# use align_cds() function after preprocess_cds.


# 4. Dimensionality Reduction (UMAP)
# ----------------------------------

# Reduce dimensionality using UMAP (standard in Monocle 3)
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# Visualize cells in UMAP space, colored by metadata
plot_cells(cds, color_cells_by = "Time", label_cell_groups = FALSE, cell_size=1) + ggtitle("UMAP by Time")
ggsave("monocle3_umap_by_time.png")

plot_cells(cds, color_cells_by = "State", label_cell_groups = FALSE, cell_size=1) + ggtitle("UMAP by State")
ggsave("monocle3_umap_by_state.png")


# 5. Clustering Cells (Optional but helpful for finding roots)
# ------------------------------------------------------------

# Cluster cells in the UMAP space using Leiden community detection
cds <- cluster_cells(cds, resolution = 0.005)

# Visualize clusters on UMAP
plot_cells(cds, color_cells_by = "partition", group_cells_by="partition", label_groups_by_cluster = FALSE) + ggtitle("UMAP by Partition/Cluster")
ggsave("monocle3_umap_by_partition.png")


# 6. Learning the Trajectory Graph
# --------------------------------

# Learn the principal graph representing the trajectory backbone
cds <- learn_graph(cds, use_partition = TRUE)

# Visualize the learned graph overlaid on the UMAP plot
# Nodes represent centers of cell populations, edges represent transitions
plot_cells(cds,
           color_cells_by = "State",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           graph_label_size = 3,
           cell_size = 1) + ggtitle("Learned Trajectory Graph")
ggsave("monocle3_trajectory_graph.png")


# 7. Ordering Cells in Pseudotime
# -------------------------------

# Order cells along the trajectory graph to calculate pseudotime
# Identify the root(s) of the trajectory (the starting point)
# Requires biological knowledge (e.g., which cells are the earliest time point,
# which express progenitor markers, which cluster corresponds to the start state).

# Method 1: Programmatically get the root node(s) based on earliest time point/state
get_earliest_principal_node <- function(cds, time_bin="0"){
  cell_ids <- which(colData(cds)[, "Time"] == time_bin) # Find cells at Time 0
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

# Find the principal graph node(s) associated with the earliest time point ("0")
root_nodes <- get_earliest_principal_node(cds, time_bin = "0")
print(paste("Identified root node(s):", paste(root_nodes, collapse=", ")))

# Method 2: Manually specify root based on inspecting plots
# e.g., if UMAP cluster 'X' clearly represents the start based on Time/State plots:
# root_nodes <- colnames(cds)[clusters(cds) == X] # Get cell IDs in that cluster
# Or identify node names by clicking on plot_cells(cds, label_principal_points=TRUE)

# Calculate pseudotime, starting from the identified roots
cds <- order_cells(cds, root_principal_nodes = root_nodes)

# Visualize pseudotime on the UMAP plot
# Cells should be ordered from low (blue/purple) to high (yellow) pseudotime
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           label_roots = TRUE, # Show the root node(s) used
           graph_label_size = 3,
           cell_size = 1) + ggtitle("Cells Ordered by Pseudotime")
ggsave("monocle3_pseudotime_umap.png")


# 8. Finding Trajectory-Dependent Genes
# -------------------------------------

# Identify genes whose expression changes significantly as a function of pseudotime
# Uses regression analysis against the pseudotime values along the graph
print("Finding trajectory-dependent genes (this may take a while)...")
# Running graph_test with default settings
graph_test_results <- graph_test(cds, neighbor_graph = "principal_graph", cores = parallel::detectCores())

# Select top significant genes (e.g., based on q-value)
top_genes <- graph_test_results %>%
  arrange(q_value) %>%
  filter(status == "OK") %>%
  head(5)

print("Top 6 trajectory-dependent genes:")
print(select(top_genes, gene_short_name, status, morans_I, q_value))

# Visualize expression of these top genes as a function of pseudotime
plot_genes_in_pseudotime(cds[top_genes$gene_short_name, ],
                         color_cells_by = "State",
                         min_expr = 0.5, ncol=2)
ggsave("monocle3_top_genes_pseudotime_lineplot.png")

# Visualize expression of specific known markers along pseudotime
# fibroblast_markers <- c("COL1A1", "PDGFRA", "VIM") # Example fibroblast markers
# plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name %in% fibroblast_markers, ], ...)

# Plot a heatmap of multiple trajectory-dependent genes ordered by pseudotime
top_20_genes <- graph_test_results %>% arrange(q_value) %>% filter(status == "OK") %>% head(20) %>% pull(gene_short_name)
plot_pseudotime_heatmap(cds[top_20_genes, ], # Use gene short names if available
                        num_clusters = 3, # Group genes by expression pattern
                        cores = parallel::detectCores(),
                        show_rownames = TRUE,
                        return_heatmap = TRUE)
ggsave("monocle3_top_genes_pseudotime_heatmap.png", width=6, height=8)


# 9. Branch Analysis
# -----------------------------
# If the trajectory has branches (representing cell fate decisions):
# - Branch points are often automatically labeled by plot_cells (label_branch_points=TRUE)
# - `graph_test` can also identify genes differentially expressed *between* branches after a specific branch point.


# 10. Saving Results
# -----------------
# Save the processed CDS object containing trajectory and pseudotime info
# saveRDS(cds, file = "fibroblast_reprogramming_monocle3.rds")
# print("Pseudotime analysis complete. Monocle 3 object saved.")

print("Basic Monocle 3 pseudotime analysis workflow completed.")
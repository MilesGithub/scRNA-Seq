# Define a list of required packages
packages <- c("Seurat", "DropletUtils", "SingleCellExperiment", "Matrix", "DoubletFinder")

# Function to install missing packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))


# Load 10X Genomics data into a SingleCellExperiment object
sce <- DropletUtils::read10xCounts("D:/projects/CBW_scRNA/scRNA_Data/Korrapati_MouseEar/filtered_feature_bc_matrix/")
sce

# class: SingleCellExperiment 
# dim: 32285 552 
# metadata(1): Samples
# assays(1): counts
# rownames(32285): ENSMUSG00000051951 ENSMUSG00000089699 ... ENSMUSG00000095019 ENSMUSG00000095041
# rowData names(3): ID Symbol Type
# colnames: NULL
# colData names(2): Sample Barcode
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):

# Read sparse matrix from a Matrix Market format file
sparsematrix <- Matrix::readMM("D:/projects/CBW_scRNA/scRNA_Data/Korrapati_MouseEar/filtered_feature_bc_matrix/matrix.mtx.gz")

# Read row names (gene names) and column names (cell barcodes)
my_rownames <- read.delim("D:/projects/CBW_scRNA/scRNA_Data/Korrapati_MouseEar/filtered_feature_bc_matrix/features.tsv.gz", header=F)
my_colnames <- read.delim("D:/projects/CBW_scRNA/scRNA_Data/Korrapati_MouseEar/filtered_feature_bc_matrix/barcodes.tsv.gz", header=F)

head(my_rownames)
# V1      V2              V3
# 1 ENSMUSG00000051951    Xkr4 Gene Expression
# 2 ENSMUSG00000089699  Gm1992 Gene Expression
# 3 ENSMUSG00000102331 Gm19938 Gene Expression
# 4 ENSMUSG00000102343 Gm37381 Gene Expression
# 5 ENSMUSG00000025900     Rp1 Gene Expression
# 6 ENSMUSG00000025902   Sox17 Gene Expression

rownames(sparsematrix) <- my_rownames[,1]
colnames(sparsematrix) <- my_colnames[,1]

df<-data.frame(sparsematrix)

class(sparsematrix)
dim(sparsematrix)
#32285   552

# Create a SingleCellExperiment object from the sparse matrix
sce2 <- SingleCellExperiment(
  assays = list(counts = sparsematrix),
  rowData = my_rownames
)
sce2

# class: SingleCellExperiment 
# dim: 32285 552 
# metadata(0):
#   assays(1): counts
# rownames(32285): ENSMUSG00000051951 ENSMUSG00000089699 ... ENSMUSG00000095019 ENSMUSG00000095041
# rowData names(3): V1 V2 V3
# colnames(552): AAACCTGAGCTAGTTC-1 AAACGGGAGTCGTTTG-1 ... TTTGCGCCATGTAGTC-1 TTTGTCAAGAGACTTA-1
# colData names(0):
#   reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):

# Add a new column to the column data indicating the sample name
colData(sce)$Sample <- "Korrapati_exons"
head(colData(sce))

# DataFrame with 6 rows and 2 columns
# Sample            Barcode
# <character>        <character>
#   1 Korrapati_exons AAACCTGAGCTAGTTC-1
# 2 Korrapati_exons AAACGGGAGTCGTTTG-1
# 3 Korrapati_exons AAACGGGCAGCTCCGA-1
# 4 Korrapati_exons AAACGGGTCACCACCT-1
# 5 Korrapati_exons AAAGTAGCACATCCGG-1
# 6 Korrapati_exons AAAGTAGTCGGAAATA-1

# Display the counts matrix from the SingleCellExperiment object
counts(sce)
assay(sce, "counts")
#32285 x 552 sparse Matrix of class "dgCMatrix"

# Load data into a Seurat object from 10X Genomics format
sparsematrix <- Seurat::Read10X("D:/projects/CBW_scRNA/scRNA_Data/Korrapati_MouseEar/filtered_feature_bc_matrix/")
seur_obj <- Seurat::CreateSeuratObject(sparsematrix, project="Korrapati_exons")
seur_obj

# An object of class Seurat 
# 32285 features across 552 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 1 layer present: counts

# Access the counts matrix from the RNA assay in the Seurat object
Seurat::GetAssayData(seur_obj, assay="RNA", "counts")[1:5,1:5]
seur_obj@assays$RNA@layers$counts[1:5,1:5]

#5 x 5 sparse Matrix of class "dgCMatrix"

# Display the metadata associated with cells in the Seurat object
head(seur_obj@meta.data)

# orig.ident nCount_RNA nFeature_RNA
# AAACCTGAGCTAGTTC-1 Korrapati_exons      17174         3933
# AAACGGGAGTCGTTTG-1 Korrapati_exons        639          427
# AAACGGGCAGCTCCGA-1 Korrapati_exons      13500         3012
# AAACGGGTCACCACCT-1 Korrapati_exons       5823         2458
# AAAGTAGCACATCCGG-1 Korrapati_exons      10911         2824
# AAAGTAGTCGGAAATA-1 Korrapati_exons       6335         2083

# Create a new assay object for introns and add it to the Seurat object
introns <- Seurat::CreateAssayObject(sparsematrix)
seur_obj@assays$Introns <- introns

# Calculate the total counts and number of features for each cell and add them to column data
colData(sce)$nCount <- colSums(counts(sce))
colData(sce)$nFeatures <- colSums(counts(sce) > 0)

# Filter cells based on the number of counts and display the dimensions of the filtered object
filtered_sce <- sce[,colData(sce)$nCount > 1000]
dim(filtered_sce)

#[1] 32285   476

# Filter out mitochondrial genes and display the dimensions of the filtered object
filtered_sce <- sce[!grepl("^mt-", rowData(sce)$Symbol),]
dim(filtered_sce)

#[1] 32272   552

### comparison of mice single cell and single nuclei data
setwd("/Users/cristalvillalba/Downloads/DRG_integration/meta_atlas/samples/comparison_sn_sc")
library(Seurat)
library(SeuratData)
library(ggplot2)
library(dplyr)



DimPlot(DRG_nonneurons, reduction = "umap")

DimPlot(DRG_nonneurons, reduction = "umap", group.by = "Cell_or_Nuclei")
# Generate the UMAP plots split by "Cell_or_Nuclei"
DimPlot(DRG_nonneurons, reduction = "umap", group.by = "Atlas_annotation", split.by = "Cell_or_Nuclei") +
  theme_minimal() +
  ggtitle("UMAP Split by Cell or Nuclei") +
  theme(legend.position = "right")

# Load the Seurat package
library(Seurat)

# Check the metadata of the Seurat object to confirm the column exists
head(DRG_nonneurons@meta.data)

# Subset the Seurat object for cells labeled as "Cell"
cell_subset <- subset(DRG_nonneurons, subset = Cell_or_Nuclei == "Cell")

# Subset the Seurat object for nuclei labeled as "Nuclei"
nuclei_subset <- subset(DRG_nonneurons, subset = Cell_or_Nuclei == "Nuclei")

# Optional: Check the resulting subsets
print(cell_subset)
print(nuclei_subset)

# Remove reductions from the subset
cell_subset@reductions <- list()  # Removes existing reductions
nuclei_subset@reductions <- list()
# Find variable features for the subset
cell_subset <- FindVariableFeatures(cell_subset)
nuclei_subset <- FindVariableFeatures(nuclei_subset)

# Recalculate PCA for the subset
cell_subset <- RunPCA(cell_subset)
nuclei_subset <- RunPCA(nuclei_subset)

# Recalculate UMAP
cell_subset <- RunUMAP(cell_subset, dims = 1:30)  # Adjust dims based on your PCA
nuclei_subset <- RunUMAP(nuclei_subset, dims = 1:30)

DimPlot(cell_subset, reduction = "umap")

# Reset the graphics device
dev.off()  # Shut down the current graphics device
DimPlot(cell_subset, reduction = "umap", group.by = "Cell_or_Nuclei")  # Try the plot again

DRG_nonneurons <- RunUMAP(DRG_nonneurons, dims = 1:30)
DimPlot(DRG_nonneurons, reduction = "umap", group.by = "Atlas_annotation", split.by = "Cell_or_Nuclei")

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Ensure the metadata contains relevant columns
head(DRG_nonneurons@meta.data)

# Calculate proportions
cell_proportions <- DRG_neurons@meta.data %>%
  group_by(Atlas_annotation, Cell_or_Nuclei) %>%
  summarize(Cell_Count = n()) %>%
  mutate(Proportion = Cell_Count / sum(Cell_Count))

# Add total counts per cluster
cell_proportions <- cell_proportions %>%
  group_by(Atlas_annotation) %>%
  mutate(Total_Count = sum(Cell_Count))

# Create the bar plot
ggplot(cell_proportions, aes(x = Atlas_annotation, y = Proportion, fill = Cell_or_Nuclei)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Cell_Count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(
    title = "Cell Proportion by Cell Type and Group",
    x = "Cell Type (Atlas Annotation)",
    y = "Proportion",
    fill = "Cell or Nuclei"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
DRG_neurons[["percent.mt"]] <- PercentageFeatureSet(DRG_neurons, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(DRG_neurons, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualize QC metrics as violin plots, split by Cell_or_Nuclei
VlnPlot(
  DRG_neurons, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
  ncol = 3, 
  group.by = "Cell_or_Nuclei"
)


# Save or use the subsets for downstream analysis
saveRDS(cell_subset, file = "drg_nonneuron_cells.rds")
saveRDS(nuclei_subset, file = "drg_nonneuron_nuclei.rds")

# Check the cell names in the main object
all_cells <- colnames(DRG_nonneurons)

# Check if cell names in assays are consistent
assay_cells <- colnames(DRG_nonneurons[["RNA"]])
identical(all_cells, assay_cells)  # Should return TRUE

# Check if cell names in reductions are consistent
reduction_cells <- rownames(DRG_nonneurons@reductions$umap@cell.embeddings)
identical(all_cells, reduction_cells)  # Should return TRUE

# Check if cell names in graphs are consistent
graph_cells <- rownames(DRG_nonneurons@graphs$integrated_nn)
identical(all_cells, graph_cells)  # Should return TRUE


### identical returns FALSE
# Check the differences between the main object and the RNA assay
setdiff(all_cells, assay_cells)  # Cells in the main object but not in the RNA assay
setdiff(assay_cells, all_cells)  # Cells in the RNA assay but not in the main object
# Get the intersection of cell names
shared_cells <- intersect(all_cells, assay_cells)

# Subset the Seurat object to include only the shared cells
DRG_nonneurons <- subset(DRG_nonneurons, cells = shared_cells)
# Check again to ensure alignment
all_cells <- colnames(DRG_nonneurons)
assay_cells <- colnames(DRG_nonneurons[["RNA"]])
identical(all_cells, assay_cells)  # Should now return TRUE

### FALSE AGAIN
# Identify cells in the main object but not in the assay
missing_in_assay <- setdiff(all_cells, assay_cells)

# Identify cells in the assay but not in the main object
missing_in_object <- setdiff(assay_cells, all_cells)

# Check if it's a problem of order rather than missing cells
ordered_cells_issue <- setequal(all_cells, assay_cells)  # TRUE means only the order differs
# Reorder the RNA assay to match the order in the main object
DRG_nonneurons[["RNA"]] <- DRG_nonneurons[["RNA"]][, all_cells]
# Subset and reorder the RNA assay to match the main object
RNA_assay <- DRG_nonneurons[["RNA"]]
RNA_assay <- RNA_assay[, all_cells]  # Reorder the assay
# Replace the RNA assay in the Seurat object
DRG_nonneurons[["RNA"]] <- RNA_assay
all_cells <- colnames(DRG_nonneurons)
assay_cells <- colnames(DRG_nonneurons[["RNA"]])
identical(all_cells, assay_cells)  # Should return TRUE
validObject(DRG_nonneurons)  # Should return TRUE

# Extract the RNA assay and metadata
RNA_assay <- DRG_nonneurons[["RNA"]]
metadata <- DRG_nonneurons@meta.data

# Find the intersection of cell names in the metadata and RNA assay
shared_cells <- intersect(colnames(RNA_assay), rownames(metadata))

# Subset both RNA assay and metadata to shared cells
RNA_assay <- RNA_assay[, shared_cells]
metadata <- metadata[shared_cells, ]

# Rebuild the Seurat object
DRG_nonneurons_fixed <- CreateSeuratObject(counts = RNA_assay, meta.data = metadata)

# Check the structure of RNA_assay
str(RNA_assay)

# If RNA_assay is not a matrix or sparse matrix, extract the counts explicitly
RNA_assay <- GetAssayData(DRG_nonneurons[["RNA"]], layer = "counts")
# Check dimensions of RNA_assay
dim(RNA_assay)

# Check dimensions of metadata
dim(metadata)

# Ensure cell names align
identical(colnames(RNA_assay), rownames(metadata))  # Should return TRUE
str(RNA_assay)
class(RNA_assay)
# Convert to a sparse matrix if not already
if (!inherits(RNA_assay, "dgCMatrix")) {
  library(Matrix)
  RNA_assay <- as(RNA_assay, "dgCMatrix")
}
# Check for non-numeric, negative, or missing values
anyNA(RNA_assay)        # Should return FALSE
any(RNA_assay < 0)      # Should return FALSE
any(RNA_assay %% 1 != 0) # Should return FALSE
RNA_assay <- round(RNA_assay)
str(metadata)
class(metadata)
anyNA(DRG_nonneurons@meta.data)  # Should return FALSE

# Check which columns contain NA values
colSums(is.na(DRG_nonneurons@meta.data))

# Check rows with missing values
metadata[is.na(DRG_nonneurons@meta.data), ]
# Remove rows with NA values
metadata <- metadata[complete.cases(DRG_nonneurons@meta.data), ]

# Replace NA values in a specific column with "Unknown"
metadata$Cell_or_Nuclei[is.na(DRG_nonneurons@meta.data)] <- "Unknown"

# Replace NA values in numeric columns with the median
metadata$SomeNumericColumn[is.na(metadata$SomeNumericColumn)] <- median(metadata$SomeNumericColumn, na.rm = TRUE)
anyNA(metadata)  # Should now return FALSE
DRG_nonneurons_fixed <- CreateSeuratObject(counts = RNA_assay, meta.data = metadata)

# Check for duplicate rownames in metadata
sum(duplicated(rownames(metadata)))  # Should return 0

# Identify the duplicate rownames
duplicated_metadata <- rownames(metadata)[duplicated(rownames(metadata))]
print(duplicated_metadata)
rownames(metadata)[duplicated(rownames(metadata))] <- make.unique(rownames(metadata))
# Check for duplicate column names in RNA_assay
sum(duplicated(colnames(RNA_assay)))  # Should return 0

# Identify duplicate column names
duplicated_columns <- colnames(RNA_assay)[duplicated(colnames(RNA_assay))]
print(duplicated_columns)
RNA_assay <- RNA_assay[, !duplicated(colnames(RNA_assay))]
colnames(RNA_assay)[duplicated(colnames(RNA_assay))] <- make.unique(colnames(RNA_assay))
DRG_nonneurons_fixed <- CreateSeuratObject(counts = RNA_assay, meta.data = metadata)
# Check for duplicate rownames in RNA_assay
sum(duplicated(rownames(RNA_assay)))  # Should return 0
# Identify the duplicate rownames
duplicated_genes <- rownames(RNA_assay)[duplicated(rownames(RNA_assay))]
print(duplicated_genes)

# Combine rows with duplicate gene names by summing their counts
RNA_assay <- aggregate(as.matrix(RNA_assay), by = list(rownames(RNA_assay)), FUN = sum)

# Restore rownames after aggregation
rownames(RNA_assay) <- RNA_assay$Group.1
RNA_assay <- RNA_assay[, -1]  # Remove the grouping column
library(Matrix)

# Check if RNA_assay is sparse
if (!inherits(RNA_assay, "dgCMatrix")) {
  RNA_assay <- as(RNA_assay, "dgCMatrix")  # Convert to a sparse matrix
}

library(Matrix)

# Check if RNA_assay is sparse
if (!inherits(RNA_assay, "dgCMatrix")) {
  RNA_assay <- as(RNA_assay, "dgCMatrix")  # Convert to a sparse matrix
}
# Aggregate rows in a sparse matrix while keeping it sparse
RNA_assay <- as(RNA_assay, "dgCMatrix")
gene_names <- rownames(RNA_assay)
unique_gene_names <- unique(gene_names)

# Create an empty sparse matrix to store aggregated data
RNA_assay_aggregated <- Matrix(0, nrow = length(unique_gene_names), ncol = ncol(RNA_assay), sparse = TRUE)
rownames(RNA_assay_aggregated) <- unique_gene_names
colnames(RNA_assay_aggregated) <- colnames(RNA_assay)

# Aggregate counts for duplicate genes
for (gene in unique_gene_names) {
  gene_rows <- which(gene_names == gene)
  if (length(gene_rows) > 1) {
    RNA_assay_aggregated[gene, ] <- rowSums(RNA_assay[gene_rows, ])
  } else {
    RNA_assay_aggregated[gene, ] <- RNA_assay[gene_rows, ]
  }
}

# Replace the RNA_assay with the aggregated version
RNA_assay <- RNA_assay_aggregated

####### analysis DRG neurons single cell x single nuclei
drg_non_neurons_single_cell   <- NormalizeData(drg_non_neurons_single_cell)
# plot variable features with and without labels
drg_non_neurons_single_cell <- FindVariableFeatures(drg_non_neurons_single_cell, selection.method = "vst", nfeatures = 26595)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(drg_non_neurons_single_cell), 10)
plot1 <- VariableFeaturePlot(drg_non_neurons_single_cell)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale the data
drg_non_neurons_single_cell   <- ScaleData(drg_non_neurons_single_cell)

# Perform PCA
drg_non_neurons_single_cell   <- RunPCA(drg_non_neurons_single_cell)
drg_non_neurons_single_cell <- FindNeighbors(drg_non_neurons_single_cell, dims = 1:30)
drg_non_neurons_single_cell <- FindClusters(drg_non_neurons_single_cell, resolution = 0.08)
drg_non_neurons_single_cell <- RunUMAP(drg_non_neurons_single_cell, dims = 1:30)

DimPlot(drg_non_neurons_single_cell, reduction = "umap")

## Change the idents to use atlas annotation
# Rename Idents based on the Atlas_annotation column
Idents(drg_neuron_nuclei) <- drg_neuron_nuclei@meta.data$Atlas_annotation

library(dplyr)
### find markers for single cell - nuclei
drgnuclei.markers <- FindAllMarkers(drg_neuron_nuclei, only.pos = FALSE)
drgnuclei.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(drgnuclei.markers, "drgmarkers4singlenuclei.csv")


#saverds
saveRDS(drg_non_neurons_single_cell, "drg_non_neurons_single_cell_clustered.RDS")



##### other plots  schwann ###
FeaturePlot(object = drg_non_neurons_single_cell_clustered, features = c("Mpz", "Mbp", "Plp1", "Apoe"))
RidgePlot(object = drg_non_neurons_single_cell_clustered, features = c("Apoe"))
VlnPlot(object = drg_non_neurons_single_cell_clustered, features = c("Apoe"))

### satglia ###
FeaturePlot(object = drg_non_neurons_single_cell_clustered, features = c("Fabp7", "Mpz", "Plp1", "Apoe", "Mbp"))
RidgePlot(object = drg_non_neurons_single_cell_clustered, features = c("Fabp7"))
VlnPlot(object = drg_non_neurons_single_cell_clustered, features = c("Fabp7"))

#### SMC ####
FeaturePlot(object = drg_non_neurons_single_cell_clustered, features = c("Tagln", "Tpm1", "Myh11"))
RidgePlot(object = drg_non_neurons_single_cell_clustered, features = c("Myh11"))
VlnPlot(object = drg_non_neurons_single_cell_clustered, features = c("Myh11"))

#### Fibro 1 genes  Apod+/Dcn high fibroblasts ###

FeaturePlot(object = drg_non_neurons_single_cell_clustered, features = c("Col1a2", "Apod", "Dcn", "Vtn"))
RidgePlot(object = drg_non_neurons_single_cell_clustered, features = c("Vtn"))
VlnPlot(object = drg_non_neurons_single_cell_clustered, features = c("Vtn"))

FeaturePlot(object = drg_non_neurons_single_cell_clustered, features = c("Ccl11", "Mgp", "Il33", "Pdgfra"))
RidgePlot(object = drg_non_neurons_single_cell_clustered, features = c("Ccl11"))
VlnPlot(object = drg_non_neurons_single_cell_clustered, features = c("Ccl11"))

FeaturePlot(object = drg_non_neurons_single_cell_clustered, features = c("Calca"))
RidgePlot(object = drg_non_neurons_single_cell_clustered, features = c("Calca"))
VlnPlot(object = drg_non_neurons_single_cell_clustered, features = c("Calca"))

#### METRICS TO EVALUATE BOTH TECHNOLOGIES ######
# Load libraries
library(Seurat)
library(ggplot2)

### Define gene set - neuronal gene markers
gene_set <- list(
  Neuronal = c("Calca", "Tac1", "Trpv1", "Adra2a", "Bmpr1b", "Oprk1", "Smr2", "Sstr2",
                  "Scn11a", "Mrgpra3", "Mrgprb4", "Mrgprd", "Calb1", "Ntrk3", "Grm8","S100a16",
                  "Ntrk2", "Il31ra", "Nppb", "Sst", "Fam19a4", "Th", "Cdh9", "Trpm8")
)

#Ensure genes exist in the dataset:
gene_set[[1]] <- gene_set[[1]][gene_set[[1]] %in% rownames(drg_neuron_cells)]

##Compute Gene Scores Using AddModuleScore()
drg_neuron_cells <- AddModuleScore(drg_neuron_cells, features = gene_set, name = "Neuronal")

VlnPlot(drg_neuron_cells, features = "Neuronal1", group.by = "Atlas_annotation", pt.size = 0.1) 


# Define gene set
gene_set <- list(
  Neuro = c("Calca", "Tac1", "Trpv1", "Adra2a", "Bmpr1b", "Oprk1", "Smr2", "Sstr2")
)

# Filter genes that exist in the dataset
gene_set[[1]] <- gene_set[[1]][gene_set[[1]] %in% rownames(drg_neuron_nuclei)]

# Check if gene list is empty
if (length(gene_set[[1]]) == 0) {
  stop("Error: No genes from the set exist in the dataset.")
}

# Reduce number of control genes
drg_neuron_nuclei <- AddModuleScore(drg_neuron_nuclei, features = gene_set, name = "Neuro", ctrl = 5)

VlnPlot(drg_neuron_nuclei, features = "Neuro1", group.by = "Atlas_annotation", pt.size = 0.1) 

# Extract metadata
meta_data <- drg_neuron_nuclei@meta.data

# Ensure "cell_type" column exists (modify if needed)
if (!"Atlas_annotation" %in% colnames(meta_data)) {
  stop("Your metadata does not contain a 'cell_type' column. Please check your dataset.")
}

# Select relevant columns
gene_score_table <- meta_data[, c("Atlas_annotation", "Neuro1")]

# View first few rows
head(gene_score_table)

library(dplyr)
library(Polychrome)
# Compute mean score per cell type
gene_score_summary <- gene_score_table %>%
  group_by(Atlas_annotation) %>%
  summarise(Avg_Score = mean(Neuro1, na.rm = TRUE))

# View summarized table
print(gene_score_summary)
write.csv(gene_score_summary, "gene_module_scoresneuro1singlenuclei.csv", row.names = FALSE)

library(viridis)

# Generate 18 colors from viridis palette
colors_18 <- viridis(18, option = "plasma")

# Plot with these colors
ggplot(gene_score_summary, aes(x = reorder(Atlas_annotation, -Avg_Score), y = Avg_Score, fill = Atlas_annotation)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Gene Module Scores by Cell Type",
       x = "Cell Type",
       y = "Average Gene Module Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = colors_18)


############## module 2 ##################

# Define gene set
gene_set <- list(
  Neuro2 = c("Scn11a", "Mrgpra3", "Mrgprb4", "Mrgprd", "Calb1", "Ntrk3", "Grm8","S100a16")
)

# Filter genes that exist in the dataset
gene_set[[1]] <- gene_set[[1]][gene_set[[1]] %in% rownames(drg_neuron_nuclei)]

# Check if gene list is empty
if (length(gene_set[[1]]) == 0) {
  stop("Error: No genes from the set exist in the dataset.")
}

# Reduce number of control genes
drg_neuron_nuclei <- AddModuleScore(drg_neuron_nuclei, features = gene_set, name = "Neuro2", ctrl = 5)
VlnPlot(drg_neuron_nuclei, features = "Neuro21", group.by = "Atlas_annotation", pt.size = 0.1) 

# Extract metadata
meta_data <- drg_neuron_nuclei@meta.data

# Ensure "cell_type" column exists (modify if needed)
if (!"Atlas_annotation" %in% colnames(meta_data)) {
  stop("Your metadata does not contain a 'cell_type' column. Please check your dataset.")
}

# Select relevant columns
gene_score_table <- meta_data[, c("Atlas_annotation", "Neuro21")]

# View first few rows
head(gene_score_table)

library(dplyr)
library(Polychrome)
# Compute mean score per cell type
gene_score_summary <- gene_score_table %>%
  group_by(Atlas_annotation) %>%
  summarise(Avg_Score = mean(Neuro21, na.rm = TRUE))

# View summarized table
print(gene_score_summary)
write.csv(gene_score_summary, "gene_module_scoresneuro21singlenuclei.csv", row.names = FALSE)

library(viridis)

# Generate 18 colors from viridis palette
colors_18 <- viridis(18, option = "plasma")

# Plot with these colors
ggplot(gene_score_summary, aes(x = reorder(Atlas_annotation, -Avg_Score), y = Avg_Score, fill = Atlas_annotation)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Gene Module Scores Neuro21",
       x = "Cell Type",
       y = "Average Gene Module Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = colors_18)



### both sn and sc
FeaturePlot(object = drg_neuron_cells, features = c("Calca"))
RidgePlot(object = drg_neuron_cells, features = c("Bmpr1b"), group.by = "Atlas_annotation")
VlnPlot(DRG_neurons, features = "Bmpr1b", group.by = "Atlas_annotation", pt.size = 0.2) 


## Customize feature plot
FeaturePlot(object = drg_neuron_nuclei, features = "Calca", cols = c("lightgray", "blue", "red"))

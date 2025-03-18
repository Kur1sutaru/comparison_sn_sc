setwd("D:/Baylor/Leelab/singlecell/scprednew")
library("scPred")
library("glue")
library("dplyr")
library("Seurat")
library("magrittr")
library("caret")
options(Seurat.object.assay.version = "v5")
## to avoid compatibility errors, install by devtools::install_github(repo="powellgenomicslab/scPred",  ref="9f407b7436f40d44224a5976a94cc6815c6e837f")
#Training the model


reference <- NormalizeData(reference)
reference <- FindVariableFeatures(reference, selection.method = "vst", nfeatures = 26595)
all.genes <- rownames(reference)
reference <- ScaleData(reference, features = all.genes)
reference <- RunPCA(reference, features = VariableFeatures(object = reference)) 
reference <- FindNeighbors(reference, dims = 1:30)
reference <- FindClusters(reference, resolution = 0.5)
reference <- RunUMAP(reference, dims = 1:30)

DimPlot(reference, group.by = "Atlas_annotation", label = TRUE, repel = TRUE)
reference <- getFeatureSpace(reference, "Atlas_annotation")
reference <- trainModel(reference)
get_probabilities(reference) %>% head()
get_scpred(reference)
plot_probabilities(reference)

### After training the model with the reference, lets predict the query - our data set
# Normalize the data
subset_sample4scaleusinggeneset   <- NormalizeData(subset_sample4scaleusinggeneset)

# Find variable features
subset_sample4scaleusinggeneset   <- FindVariableFeatures(subset_sample4scaleusinggeneset)

# Scale the data
subset_sample4scaleusinggeneset   <- ScaleData(subset_sample4scaleusinggeneset)

# Perform PCA
subset_sample4scaleusinggeneset   <- RunPCA(subset_sample4scaleusinggeneset)

# Check available assays
Assays(subset_sample4scaleusinggeneset)

# Make sure "integrated" is the default assay for predictions
subset_sample4scaleusinggeneset   <- scPredict(subset_sample4scaleusinggeneset, threshold = 0.8, reference)
DimPlot(subset_sample4scaleusinggeneset, group.by = "scpred_prediction", reduction = "scpred")
subset_sample4scaleusinggeneset   <- RunUMAP(subset_sample4scaleusinggeneset, reduction = "scpred", dims = 1:30)
DimPlot(subset_sample4scaleusinggeneset, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
# Store the trained scPred model in the 'misc' slot of the reference Seurat object
subset_sample4scaleusinggeneset@misc$scPred_model <- reference@tools$scPred
# Calculate proportions of cell types
cell_proportions <- table(subset_sample4scaleusinggeneset@meta.data$scpred_prediction)
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("CellType", "Count")
cell_proportions_df$Proportion <- cell_proportions_df$Count / sum(cell_proportions_df$Count)
# Load ggplot2 for plotting
library(ggplot2)

# Create the bar plot with cell counts as labels
ggplot(cell_proportions_df, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +  # Add cell count labels above each bar
  theme_minimal() +
  labs(x = "Predicted Cell Type", y = "Proportion", title = "Cell Type Proportions from sample 4 subset") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


saveRDS(subset_sample4scaleusinggeneset, "subset_sample4scaleusinggeneset_annotatedscpred.rds")
saveRDS(bl_68774 , "bl_68774_annotated.rds")
saveRDS(subset_sample4scaleusinggeneset , "subset_sample4scaleusinggeneset_annotated.rds")
saveRDS(subset_sample4scaleusinggeneset, "subset_sample4scaleusinggeneset_annotated.rds")


### Accessing classifiers
crossTab(subset_sample4scaleusinggeneset  , "cell_type", "scpred_prediction")
crossTab(subset_sample4scaleusinggeneset  , "cell_type", "scpred_prediction", output = "prop")
get_classifiers(reference)

# Each model can be normally treated using the caret enviroment. 
# For example, we can plot the performance resamples using the plot.train:
caret::plot.train(get_classifiers(reference)[["Chondrocytes"]])


# Calculate proportions of cell types
cell_proportions <- table(reference@meta.data$cell_type)
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("CellType", "Count")
cell_proportions_df$Proportion <- cell_proportions_df$Count / sum(cell_proportions_df$Count)
# Load ggplot2 for plotting
library(ggplot2)

# Create the bar plot
ggplot(cell_proportions_df, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Predicted Cell Type", y = "Proportion", title = "Cell Type Proportions from scPred") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


### perform degs between wt and ko samples 
wtcalv<- subset_sample2scaleusinggeneset_annotatedscpred
kocalv<- subset_sample4scaleusinggeneset_annotatedscpred

# Ensure that orig.ident is a factor
wtcalv@meta.data$orig.ident <- as.factor(wtcalv@meta.data$orig.ident)

# Modify the levels, replacing 'SeuratProject' with 'WT'
levels(wtcalv@meta.data$orig.ident)[levels(wtcalv@meta.data$orig.ident) == "SeuratProject"] <- "WT"

## now for KO
# Ensure that orig.ident is a factor
kocalv@meta.data$orig.ident <- as.factor(kocalv@meta.data$orig.ident)

# Modify the levels, replacing 'SeuratProject' with 'WT'
levels(kocalv@meta.data$orig.ident)[levels(kocalv@meta.data$orig.ident) == "SeuratProject"] <- "KO"



# Merge the two Seurat objects
calvariadeg <- merge(wtcalv, y = kocalv, add.cell.ids = c("WT", "KO"), project = "calvaria")
# Combine cell type and condition into a new column
calvariadeg$celltype.stim <- paste(calvariadeg$scpred_prediction, calvariadeg$orig.ident, sep = "_")
Idents(calvariadeg) <- "celltype.stim"


# Differential expression analysis for a specific cell type
osteochondroprog.decalv <- FindMarkers(calvariadeg, ident.1 = "Osteochondrogenic progenitors_WT", ident.2 = "Osteochondrogenic progenitors_KO", verbose = FALSE)

# View the top differentially expressed genes
head(transitioning.decalv, n = 10)
write.csv(transitioning.decalv, "transitioningwtxkocalv.csv")

##### volcano and heatmaps
library(EnhancedVolcano)
# Assume 'mono.de' is the differential expression results data frame
EnhancedVolcano(
  transitioning.decalv,
  lab = rownames(transitioning.decalv),
  x = 'avg_log2FC',        # Log fold change column
  y = 'p_val',         # Adjusted p-value column
  title = 'Volcano Plot: Transitioning cells WT vs KO in calvaria',
  xlab = 'Log2 Fold Change',
  ylab = '-Log10 P-value',
  pCutoff = 0.05,          # Significance threshold for adjusted p-value
  FCcutoff = 0.25,         # Fold change threshold
  pointSize = 2.0,
  labSize = 3.0
)

library(ggplot2)

library(ggplot2)
library(ggrepel)

# Add -log10(p-value) to the data
chondro.de$neg_log10_pval <- -log10(chondro.de$p_val_adj)

# Define significance criteria for labeling genes
significant_genes <- chondro.de[chondro.de$p_val_adj < 0.05 & abs(chondro.de$avg_log2FC) > 0.25, ]

# Create the ggplot volcano plot with gene names
library(ggplot2)

# Add -log10(p-value) to the data in chondro.de
chondro.de$neg_log10_pval <- -log10(chondro.de$p_val_adj)

# Define significance criteria for labeling genes
significant_genes <- chondro.de[chondro.de$p_val_adj < 0.05 & abs(chondro.de$avg_log2FC) > 0.25, ]

# Create the ggplot volcano plot with gene names
ggplot(chondro.de, aes(x = avg_log2FC, y = neg_log10_pval)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "blue") +
  geom_text(
    data = significant_genes,
    aes(label = rownames(significant_genes)),
    size = 3,
    vjust = 1,       # Adjust vertical position
    hjust = 1,       # Adjust horizontal position
    check_overlap = TRUE # Avoid overlapping text
  ) +
  labs(
    title = "Volcano Plot: Chondrocytes of Long Bone WT vs KO",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal()

saveRDS(longbonedeg, "longbonewtandkocombined.rds")

# Extract the 'scpred_prediction' column from the metadata
scpred_predictions <- calvariadeg@meta.data$scpred_prediction

# Get the unique cell type predictions
unique_cell_types <- unique(scpred_predictions)

# Display the unique cell types
print(unique_cell_types)



## perform gene module score analysis

# Define the gene set
gene_set <- list(
  OsteogenesisGeneSet = c(
    "Runx2", "Sp7", "Alpl", "Bglap", "Dmp1", "Sox9", "Col1a1", 
    "Col2a1", "Ctsk", "Ly6a", "Lepr", "Grem1", "Acan", "Adipoq",
    "Pparg", "Fabp4", "Mki67", "Ocn", "Cxcl12", "Prrx1", 
    "Cd140a", "Cd105", "Cd200", "Acta2", "Pthrp"
  )
)


# Ensure the gene names are in the correct format
gene_set <- lapply(gene_set, function(genes) {
  genes[genes %in% chondrocytesdegswtxkolb$X]  # Filter genes that exist in your data
})

# Assuming 'chondro' is your Seurat object
calvariadeg <- AddModuleScore(
  object = calvariadeg,
  features = gene_set,
  name = "GeneSetScore"
)

# Check the results
head(calvariadeg@meta.data)

### Visualize the Module Score
VlnPlot(calvariadeg, features = "GeneSetScore1", group.by = "scpred_prediction")
FeaturePlot(calvariadeg, features = "GeneSetScore1")
write.csv(longbonedeg@meta.data$GeneSetScore1, "genesetscoreforgeneset.csv")


### GSEA
library(clusterProfiler)
library(org.Mm.eg.db)

# Assuming 'de_results' is your differential expression result data frame from Seurat
# Select significantly upregulated genes (e.g., adjusted p-value < 0.05)
# Subset significant genes from chondro.de using row names
significant_genes <- rownames(subset(chondro.decalv, p_val < 0.05))
# GO enrichment analysis using mouse gene annotations
go_enrichment <- enrichGO(
  gene = significant_genes,
  OrgDb = org.Mm.eg.db,  # Use the mouse gene annotation database
  keyType = "SYMBOL",    # Assuming your gene identifiers are symbols
  ont = "BP",            # 'BP' for Biological Process, 'MF' for Molecular Function, 'CC' for Cellular Component
  pAdjustMethod = "none",  # Benjamini-Hochberg correction
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# View the top results
head(go_enrichment)
write.csv(go_enrichment, "(chondrocytescalvgsearesults.csv")
barplot(go_enrichment, showCategory = 10, title = "GO Enrichment in Chondrocytes calvaria wtxko")

library(dplyr)
cell_type_counts <- subset_sample4scaleusinggeneset_annotatedscpred@meta.data %>%
  count(scpred_prediction)
print(cell_type_counts)


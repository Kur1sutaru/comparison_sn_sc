library(Seurat)
library(dplyr)
library(Matrix)
library(edgeR)

# Assuming:
# - scpred_prediction = predicted cluster
# - orig.ident = sample or condition (DMM, Exercise, Sham)

# Get all unique clusters
clusters <- unique(dmmbatch$scpred_prediction)

# Create list to hold results
pseudobulk_results <- list()

# Loop through each cluster
for (clust in clusters) {
  message("Processing cluster: ", clust)
  
  # Subset cells for current cluster
  cells <- WhichCells(dmmbatch, expression = scpred_prediction == clust)
  sub <- subset(dmmbatch, cells = cells)
  
  # Extract counts and metadata
  counts <- GetAssayData(sub, slot = "counts")
  metadata <- sub@meta.data
  
  # Create pseudobulk matrix: sum counts by orig.ident (sample)
  pb_counts <- t(rowsum(t(as.matrix(counts)), group = metadata$orig.ident))
  
  # Create group mapping (DMM, Exercise, Sham)
  group_info <- metadata %>%
    distinct(orig.ident, .keep_all = TRUE) %>%
    select(orig.ident, group = orig.ident) %>%
    column_to_rownames("orig.ident")
  
  # Build DGEList object
  dge <- DGEList(counts = pb_counts, group = group_info$group)
  dge <- calcNormFactors(dge)
  
  # Create design matrix
  design <- model.matrix(~0 + group_info$group)
  colnames(design) <- levels(factor(group_info$group))
  
  # Estimate dispersion
  dge <- estimateDisp(dge, design)
  
  # Fit model and do DE test
  fit <- glmQLFit(dge, design)
  
  # Example contrast: DMM vs Sham
  contrast <- makeContrasts(DMMvsSham = DMM - Sham, levels = design)
  res <- glmQLFTest(fit, contrast = contrast)
  top <- topTags(res, n = Inf)$table
  pseudobulk_results[[clust]] <- top
}

# Example: view results for one cluster
pseudobulk_results[["ClusterName"]]

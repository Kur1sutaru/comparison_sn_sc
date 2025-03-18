## ACLT and sham DRG single nuclei
setwd("/Users/cristalvillalba/Downloads/DRG_integration/bcmsamples")
library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(dplyr)
library(ggplot2)

# Convert and load the Seurat object
Convert("MJE_73733_3v4_ACLT_fresh_NB.h5seurat", dest = "h5seurat", overwrite = TRUE)
MJE_73733_3v4_ACLT_fresh_NB <- LoadH5Seurat("MJE_73733_3v4_ACLT_fresh_NB.h5seurat")

# Check metadata and basic stats
MJE_73733_3v4_ACLT_fresh_NB
MJE_73733_3v4_ACLT_fresh_NB <- NormalizeData(MJE_73733_3v4_ACLT_fresh_NB)
MJE_73733_3v4_ACLT_fresh_NB <- FindVariableFeatures(MJE_73733_3v4_ACLT_fresh_NB, selection.method = "vst", nfeatures = 33696)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(MJE_73733_3v4_ACLT_fresh_NB), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(MJE_73733_3v4_ACLT_fresh_NB)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(MJE_73733_3v4_ACLT_fresh_NB)
MJE_73733_3v4_ACLT_fresh_NB <- ScaleData(MJE_73733_3v4_ACLT_fresh_NB, features = all.genes)

MJE_73733_3v4_ACLT_fresh_NB <- RunPCA(MJE_73733_3v4_ACLT_fresh_NB, features = VariableFeatures(object = MJE_73733_3v4_ACLT_fresh_NB))

MJE_73733_3v4_ACLT_fresh_NB <- FindNeighbors(MJE_73733_3v4_ACLT_fresh_NB, dims = 1:30)

#test different resolutions
MJE_73733_3v4_ACLT_fresh_NB <- FindClusters(MJE_73733_3v4_ACLT_fresh_NB, resolution = 1.5)
MJE_73733_3v4_ACLT_fresh_NB <- RunUMAP(MJE_73733_3v4_ACLT_fresh_NB, dims = 1:30)

DimPlot(MJE_73733_3v4_ACLT_fresh_NB)
saveRDS(MJE_73733_3v4_ACLT_fresh_NB, file = "MJE_73733_3v4_ACLT_fresh_NB.rds")


#Gene markers - neuronal
# Calca
# Tac1
# Trpv1
# Adra2a
# Bmpr1b
# Oprk1
# Smr2
# Sstr2
# Scn11a
# Mrgpra3
# Mrgprb4
# Mrgprd
# Calb1
# Ntrk3
# Grm8
# S100a16
# Ntrk2
# Il31ra
# Nppb
# Sst
# Fam19a4
# Th
# Cdh9
# Trpm8

# Gene markers - non neuronal
# Ednrb
# Fabp7
# Cadm2
# Scn7a
# Prx
# Mbp
# Mpz
# Ptgds
# Mgp
# Dcn
# Egfl7
# Cldn5
# Pecam1
# Kcnj8
# Notch3
# Lyz2
# S100a8
# Cd3e
# Ms4a1
# Cd79a
# Cd79b
# Isg15
# Folr2
# Mrc1
# S100a4
# Plac8
# Ccr2
# S100a9
# Naaa
# Xcr1
# Itgax
# Cd69
# Cx3cr1
# Fcrls
# Bcl2a1b
# B cells Ms4a1 Cd79a Cd79b
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Ms4a1", "Cd79a","Cd79b"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Ms4a1", "Cd79a","Cd79b"))
# Calca+Adra2a neurons
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Calca", "Adra2a"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Calca", "Adra2a"))
#neurons Calca+Bmpr1b - Calca Bmpr1b Grm8 S100a16
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Calca", "Bmpr1b","Grm8", "S100a16"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Calca", "Bmpr1b","Grm8", "S100a16"))
# Calca+Oprk1 neurons -  Calca Tac1 Trpv1 Oprk1
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Calca", "Tac1","Trpv1", "Oprk1"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Calca", "Tac1","Trpv1", "Oprk1"))
#Calca+Smr2 neurons  = Calca Smr2 Grm8
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Calca", "Smr2","Grm8", "S100a16"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Calca", "Smr2","Grm8", "S100a16"))

# Calca+Sstr2+ neurons =  Calca Sstr2
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Calca", "Sstr2"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Calca", "Sstr2"))

# Ccr2 macrophages Akopian markers Cx3cr1 Fcrls Bcl2a1b
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Cx3cr1", "Fcrls", "Bcl2a1b"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Cx3cr1", "Fcrls", "Bcl2a1b"))

# Dendritic cells Naaa Xcr1 Itgax
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Naaa", "Xcr1", "Itgax"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Naaa", "Xcr1", "Itgax"))

# Endothelial Egfl7 Cldn5 Pecam1
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Egfl7", "Cldn5", "Pecam1"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Egfl7", "Cldn5", "Pecam1"))

#fibroblast Ptgds Mgp Dcn
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Ptgds", "Mgp", "Dcn"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Ptgds", "Mgp", "Dcn"))
#Immune - general markers Lyz2 S100a8 Cd3e
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Lyz2", "S100a8", "Cd3e"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Lyz2", "S100a8", "Cd3e"))
# M2-like macrophages akopian Cx3cr1 Fcrls Bcl2a1b
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Cx3cr1", "Fcrls", "Bcl2a1b"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Cx3cr1", "Fcrls", "Bcl2a1b"))
#monocytes S100a4 Plac8 Ccr2
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("S100a4", "Plac8", "Ccr2"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("S100a4", "Plac8", "Ccr2"))

#neuronal Mrgpra3+Mrgprb4 = Scn11a Mrgpra3 Mrgprb4
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Scn11a", "Mrgpra3", "Mrgprb4"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Scn11a", "Mrgpra3", "Mrgprb4"))
# Mrgpra3+Trpv1 = Scn11a Mrgpra3 Trpv1
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Scn11a", "Mrgpra3", "Trpv1"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Scn11a", "Mrgpra3", "Trpv1"))
# Mrgprd neurons = Mrgprd Scn11a
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Mrgprd", "Scn11a"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Mrgprd", "Scn11a"))
#neutrophils S100a8 S100a9
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("S100a8", "S100a9"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("S100a8", "S100a9"))
# Ntrk3 high + Ntrk2 = Calb1 Ntrk3 Ntrk2
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Calb1", "Ntrk3", "Ntrk2"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Calb1", "Ntrk3", "Ntrk2"))
# Ntrk3 high + S100a16 = Ntrk3 Grm8 S100a16. Ntrk2
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("S100a16", "Ntrk3", "Ntrk2","Grm8"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("S100a16", "Ntrk3", "Ntrk2", "Grm8"))
#. Ntrk3 low + Ntrk2  Ntrk3 Ntrk2
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Ntrk3", "Ntrk2"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Ntrk3", "Ntrk2"))
# Pericyte = Kcnj8 Notch3
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Kcnj8", "Notch3"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Kcnj8", "Notch3"))
# satglia Ednrb Fabp7
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Ednrb", "Fabp7"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Ednrb", "Fabp7"))
#schwann_M Prx Mbp Mpz
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Prx", "Mbp", "Mpz"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Prx", "Mbp", "Mpz"))
#schwann_N Cadm2 Scn7a
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Cadm2", "Scn7a"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Cadm2", "Scn7a"))
# sst neurons. Il31ra Nppb Sst
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Il31ra", "Nppb", "Sst"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Il31ra", "Nppb", "Sst"))
# T cells - Cd3e Cd69
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Cd3e", "Cd69"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Cd3e", "Cd69"))
# th neurons - Fam19a4 Th Cdh9
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Fam19a4", "Th", "Cdh9"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Fam19a4", "Th", "Cdh9"))
# tlf macrophages akopian markers = Folr2 Mrc1
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Folr2", "Mrc1"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Folr2", "Mrc1"))
# Trpm8 neurons - Trpm8
VlnPlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Trpm8"))
FeaturePlot(MJE_73733_3v4_ACLT_fresh_NB, features = c("Trpm8"))


### sample 73732 analysis high variable genes
DimPlot(MJE_73733_3v4_ACLT_fresh_NB)
all_markers <- FindAllMarkers(object = MJE_73733_3v4_ACLT_fresh_NB, 
                              only.pos = FALSE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)
write.csv(all_markers, "MJE_73733_3v4_ACLT_fresh_NBdegsclustersres1_5.csv")
library(dplyr)
top3 <- all_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 3, wt = avg_log2FC)

DoHeatmap(MJE_73733_3v4_ACLT_fresh_NB, features = top3) + scale_fill_gradientn(colors = c("blue", "white", "red"))

## RENAME CLUSTERS
levels(MJE_73733_3v4_ACLT_fresh_NB)
MJE_73733_3v4_ACLT_fresh_NB
new.cluster.ids <- c("Neuron", "SatGlia", "Fibroblast", "Fibroblast", "Schwann_M", "Endothelial",
                     "Schwann_N", "Neuron", "Neuron", "Endothelial", "Neuron", "Schwann_M",
                     "Neuron", "Endothelial", "Doublets", "Fibroblast", "Macrophage", "Neuron",
                     "Fibroblast", "Neuron", "SatGlia", "Neuron", "Pericyte", "SatGlia",
                     "Neuron", "Schwann_M", "Neuron", "Macrophage", "Doublets",
                     "B cell", "Fibroblast", "Neuron", "Neuron")
names(new.cluster.ids) <- levels(MJE_73733_3v4_ACLT_fresh_NB)  # Assign names to match original clusters
MJE_73733_3v4_ACLT_fresh_NB <- RenameIdents(MJE_73733_3v4_ACLT_fresh_NB, new.cluster.ids)
Idents(MJE_73733_3v4_ACLT_fresh_NB)
DimPlot(MJE_73733_3v4_ACLT_fresh_NB, reduction = "umap", label = TRUE, pt.size = 0.5)
# Assuming 'seurat_obj' is your Seurat object
# and 'seurat_clusters' or 'ident' contains the cluster information

MJE_73733_3v4_ACLT_fresh_NB <- subset(MJE_73733_3v4_ACLT_fresh_NB, idents = "Doublets", invert = TRUE)
DimPlot(MJE_73733_3v4_ACLT_fresh_NB, reduction = "umap", label = TRUE, pt.size = 0.5)



# Load required library
library(ggplot2)

# Extract identity classes (clusters) and count occurrences
cell_counts <- as.data.frame(table(Idents(MJE_73733_3v4_ACLT_fresh_NB)))
colnames(cell_counts) <- c("Cluster", "Count")

# Compute percentages
cell_counts$Percentage <- (cell_counts$Count / sum(cell_counts$Count)) * 100

# Ensure Cluster is a factor (for ordered plotting)
cell_counts$Cluster <- factor(cell_counts$Cluster, levels = cell_counts$Cluster[order(-cell_counts$Count)])

ggplot(cell_counts, aes(x = Cluster, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", round(Percentage, 1), "%)")), vjust = -0.5, size = 4) +
  theme_minimal() +
  labs(title = "Cell Type Proportions",
       x = "Cell Type",
       y = "Number of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

saveRDS(MJE_73733_3v4_ACLT_fresh_NB, "cellkannotMJE_73733_3v4_ACLT_fresh_NB.rds")




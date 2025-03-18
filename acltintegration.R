#aclt integration withR - after with scvi
library(Seurat)
library(dplyr)
library(ggplot2)
seurat1 <- LoadH5Seurat("~/Downloads/DRG_integration/bcmsamples/MJE_73732/MJE_73732_3v4_sham_fresh.h5seurat")
seurat2 <- LoadH5Seurat("~/Downloads/DRG_integration/bcmsamples/MJE_73733/MJE_73733_3v4_ACLT_fresh_NB.h5seurat")
seurat3 <- LoadH5Seurat("~/Downloads/DRG_integration/bcmsamples/MJE_73734/MJE_73734_3v4_sham_fresh_NB.h5seurat")
seurat1@project.name
seurat1@project.name <- "73732_sham"
seurat1@project.name
seurat2@project.name <- "73733_aclt"
seurat2@project.name
seurat3@project.name <- "73734_sham"
seurat3@project.name
seurat_list <- list(seurat1, seurat2, seurat3)
seurat_list <- lapply(seurat_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 33696)
  return(x)
})

anchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                  anchor.features = 33696)





library(SeuratDisk)
library(Seurat)
library(hdf5r)
seurat_obj <- LoadH5Seurat("metanonneurons_sn.h5Seurat")
# Ensure that 'RNA' is the default assay
DefaultAssay(drg_neuron_nuclei) <- "RNA"

# Manually set counts as main expression matrix
drg_neuron_nuclei[["X"]] <- drg_neuron_nuclei[["RNA"]]@data

# Save and Convert again
SaveH5Seurat(Non_neuron_mouse_v3singlenuclei, filename = "Non_neuron_mouse_v3singlenuclei.h5Seurat", overwrite = TRUE)
file.exists("Non_neuron_mouse_v3singlenuclei.h5Seurat")

Convert("Non_neuron_mouse_v3singlenuclei.h5Seurat", dest = "h5ad", overwrite = TRUE)

library(Seurat)
library(SeuratDisk)

library(Seurat)
library(SeuratDisk)
library(hdf5r)

# Convert Assay5 to a standard Seurat Assay
Non_neuron_mouse_v3singlenuclei[["RNA"]] <- as(Non_neuron_mouse_v3singlenuclei[["RNA"]], "Assay")

# Ensure RNA assay is set as default
DefaultAssay(Non_neuron_mouse_v3singlenuclei) <- "RNA"

# Save the converted Seurat object
SaveH5Seurat(Non_neuron_mouse_v3singlenuclei, filename = "drg_non_neuron_sn.h5Seurat", overwrite = TRUE)

# Try converting to h5ad again
Convert("drg_non_neuron_sn.h5Seurat", dest = "h5ad", overwrite = TRUE)














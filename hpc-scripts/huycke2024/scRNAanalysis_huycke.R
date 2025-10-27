library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork)

sample_names <- c("E15.5/16.5 mouse intestine all cells lane 1",
                  "E15.5/16.5 mouse intestine all cells lane 2",
                  "E15.5/16.5 mouse intestine enriched cells",
                  "E13.5/14.5/E18.5 mouse intestine all cells lane 1",
                  "E13.5/14.5/E18.5 mouse intestine all cells lane 2",
                  "E13.5/14.5/E18.5 mouse intestine enriched cells")

data <- Read10X(data.dir = "data")
seurat_obj <- CreateSeuratObject(
  counts = data
)

VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"))

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt")

seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

DimPlot(seurat_obj, reduction = "umap", label = TRUE)
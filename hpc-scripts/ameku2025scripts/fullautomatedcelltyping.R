library(Seurat)
library(Azimuth)
library(dplyr)
library(ggplot2)
library(patchwork)

options(future.globals.maxSize = Inf)
options(Seurat.object.assay.version = "v5")

reference <- readRDS("~/intestinal_reference_azimuth.rds")
reference <- fread()

load("/gpfs/data/zwicklab/raz/xeniumsubsets/epi_obj_annotated.RData")
epi.obj$cell_type_original <- epi.obj$cluster_annotation

epi.obj <- RunAzimuth(
  query = epi.obj,
  reference = "~/intestinal_reference_azimuth",
  umap.name = "umap",
  assay = "RNA",
  verbose = TRUE
)

epi.obj$cell_type_reference <- epi.obj$predicted.cell_type.l1
xenium$cell_type_reference_l2 <- xenium$predicted.cell_type.l2
xenium$reference_prediction_score <- xenium$predicted.celltype.l1.score
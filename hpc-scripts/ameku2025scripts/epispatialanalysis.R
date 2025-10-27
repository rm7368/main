suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(ggplot2)
  library(interp)
  library(rdist)
  library(Matrix)
  library(dplyr)
  library(future)
  library(BPCells)
  library(data.table)
  library(patchwork)
  library(arrow)
})

options(future.globals.maxSize = Inf)
options(Seurat.object.assay.version = "v5")
options(Seurat.checkdots = FALSE)

load("/gpfs/data/zwicklab/raz/xenium.RData")

samples = c("TA1", "TA2", "TA3", "TA4", "TA5", "TA6",
            "TA7", "TA8", "TA9", "TA10", "TA11", "TA12")

for (sample in samples) {
  cells = read_parquet(paste0("/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium/",
                       sample, "/cells.parquet"))
  
  cell_centroid_df = data.frame(
      x = cells$x_centroid,
      y = cells$y_centroid,
      cell = paste0(sample, "-", cells$cell_id),
      stringsAsFactors = FALSE
  )
  cell_centroid_df = cell_centroid_df[cell_centroid_df$cell %in% colnames(xenium),]
  centroid_data = list(
    "centroids" = CreateCentroids(cell_centroid_df)
  )
  
  coords <- CreateFOV(
    coords = centroid_data,
    type = c("centroids"),
    assay = "RNA"
  )
  
  xenium[[sample]] = coords
}

save(xenium, file="/gpfs/data/zwicklab/raz/xenium.RData")
write.csv(xenium@meta.data, "/gpfs/data/zwicklab/raz/nuclei_meta.csv")

image_dir = ("/gpfs/home/rm7368/piezoProject/image_dim_res0.5/")
dir.create(image_dir)

Idents(xenium) <- "cluster_annotation"
xenium@active.ident = factor(xenium$cluster_annotation, levels = levels(xenium@active.ident))
for (i in 1:12) {
  sample = paste0("TA", i)
  ImageDimPlot(xenium, size = 0.7, axes = T, fov = sample, cols = "polychrome")
  ggsave(paste0(image_dir,"res0.5_image_dim_plot_",sample,".png"), height = 30, width = 20)
}

DefaultAssay(xenium) <- "RNA"

image_dir = ("/gpfs/home/rm7368/piezoProject/image_dim_epi_res0.5/")
dir.create(image_dir)

for (i in 1:12) {
  sample = paste0("TA", i)
  ImageDimPlot(xenium, size=0.7, axes = T, fov = sample, cols = "polychrome", group.by = "epi_cluster_annotation")
  ggsave(paste0(image_dir,"epi_res0.5_image_dim_plot_",sample,".png"), height=30, width=20)
}

image_dir = ("/gpfs/home/rm7368/piezoProject/image_feature_Slc5a4a/")
dir.create(image_dir)

for (i in 1:12) {
  sample = paste0("TA", i)
  ImageFeaturePlot(xenium, "Slc5a4a", size = 0.7, axes = T, fov = sample)
  ggsave(paste0(image_dir, "Slc5a4a_image_", sample, ".png"), height = 30, width = 20)
}

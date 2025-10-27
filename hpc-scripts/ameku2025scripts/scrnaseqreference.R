library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork)
library(Azimuth)
library(Rcpp)

sample_names <- c("WT_Virgin_1", "WT_Virgin_2", "WT_Virgin_3", "WT_Lact_7", "WT_Lact_8",
                  "WT_Lact_9", "WT_Lact_10", "WT_Lact_11")
seurat_objects <- list()

dirs.list <- list()
dirs.list <- c("/gpfs/home/rm7368/scseq_analysis/WT_Virgin_1",
               "/gpfs/home/rm7368/scseq_analysis/WT_Virgin_2",
               "/gpfs/home/rm7368/scseq_analysis/WT_Virgin_3",
               "/gpfs/home/rm7368/scseq_analysis/WT_Lact_7",
               "/gpfs/home/rm7368/scseq_analysis/WT_Lact_8",
               "/gpfs/home/rm7368/scseq_analysis/WT_Lact_9",
               "/gpfs/home/rm7368/scseq_analysis/WT_Lact_10",
               "/gpfs/home/rm7368/scseq_analysis/WT_Lact_11"
               )

for (i in 1:8) {
  data <- Read10X(data.dir = dirs.list[i])
  seurat_obj <- CreateSeuratObject(
    counts = data,
    project = sample_names[i],
    min.cells = 3,
    min.features = 200
  )
  seurat_obj$sample <- sample_names[i]
  seurat_obj$orig.ident <- sample_names[i]
  seurat_objects[[i]] <- seurat_obj
}

combined_seurat <- merge(seurat_objects[[1]],
                         y = c(seurat_objects[[2]], seurat_objects[[3]],
                               seurat_objects[[4]], seurat_objects[[5]],
                               seurat_objects[[6]], seurat_objects[[7]],
                               seurat_objects[[8]]),
                               add.cell.ids = sample_names,
                         project = "reference")
combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "mt-")
VlnPlot(combined_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0)

combined_seurat <- subset(combined_seurat,
                          subset = ((percent.mt < 25) & (nFeature_RNA > 500) & (nCount_RNA > 100))
)

combined_seurat <- NormalizeData(combined_seurat)
combined_seurat <- FindVariableFeatures(combined_seurat,
                                        selection.method = "vst",
                                        nfeatures = 2000)
all.genes <- rownames(combined_seurat)
combined_seurat <- ScaleData(combined_seurat, features = all.genes)
combined_seurat <- RunPCA(combined_seurat, features = all.genes)
ElbowPlot(combined_seurat, ndims = 50)

combined_seurat <- FindNeighbors(combined_seurat, dims = 1:30)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.7)
combined_seurat <- RunUMAP(combined_seurat, dims = 1:30)
DimPlot(combined_seurat, reduction = "umap", label=TRUE)
DimPlot(combined_seurat, reduction = "umap", group.by = "sample")


combined_seurat <- JoinLayers(combined_seurat)
cluster_markers <- FindAllMarkers(combined_seurat,
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n=15, wt = avg_log2FC)
print(top_markers)

intestinal_markers <- list(
  "Isthmus progenitor" = c("Fgfbp1", "Mki67"),
  "Bottom-villus EC" = c("Krt19"),
  "Top-villus EC" = c("Ada", "Neat1"),
  "Lgr5+ CBC" = c("Lgr5"),
  "Mid-villus EC" = c("Slc5a1", "Slc2a2"),
  "Goblet" = c("Muc2"),
  "Enteroendocrine" = c("Chgb"),
  "Paneth" = c("Lyz1"),
  "Microfold" = c("Ccl20"),
  "Secretory progenitor" = c("Atoh1")
)

for (cell_type in names(intestinal_markers)) {
  markers <- intestinal_markers[[cell_type]]
  available_markers <- markers[markers %in% rownames(combined_seurat)]
  
  if (length(available_markers) > 0) {
    p <- FeaturePlot(combined_seurat, features = available_markers[1:min(4, length(available_markers))],
                     ncol = 2) + plot_annotation(title = cell_type)
    
    print(p)
  }
}

myplot <- DimPlot(combined_seurat, reduction = "umap", label = TRUE)
ggsave(myplot, file="/gpfs/data/zwicklab/raz/myplot.png")


cell_type_annotations <- c(
  "0" = "Mid/top-villus EC",
  "1" = "Top-villus EC",
  "2" = "Mid-villus EC",
  "3" = "Mid-villus EC",
  "4" = "Bottom-villus EC",
  "5" = "Mid/top-villus EC",
  "6" = "Top-villus EC",
  "7" = "Top-villus EC",
  "8" = "Mid-villus EC",
  "9" = "Bottom-villus EC",
  "10" = "Lgr5+ CBC",
  "11" = "Isthmus progenitor",
  "12" = "Mid-villus EC",
  "13" = "Mid-villus EC",
  "14" = "Mid/top-villus EC",
  "15" = "Secretory progenitor",
  "16" = "Goblet",
  "17" = "Paneth",
  "18" = "Enteroendocrine",
  "19" = "Microfold",
  "20" = "Unknown proliferative cells"
)


combined_seurat$cell_type <- as.character(combined_seurat$seurat_clusters) 
for (cluster in names(cell_type_annotations)) {
  combined_seurat$cell_type[combined_seurat$seurat_clusters == cluster] <- cell_type_annotations[cluster]
}                                          


Idents(combined_seurat) <- "cell_type"

combined_seurat$predicted.celltype.l1 <- combined_seurat$cell_type
combined_seurat$predicted.celltype.l2 <- combined_seurat$cell_type

pca_embeddings <- Embeddings(combined_seurat, reduction = "refDR")
n_dims <- ncol(pca_embeddings)
n_cells <- nrow(pca_embeddings)
annoy_index <- new(RcppAnnoy::AnnoyEuclidean, n_dims)
for (i in 1:n_cells) {
  annoy_index$addItem(i-1, pca_embeddings[i,])
}
annoy_index$build(50)
annoy_index$save("idx.annoy")
saveRDS(combined_seurat,file="intestinal_reference_azimuth.rds")

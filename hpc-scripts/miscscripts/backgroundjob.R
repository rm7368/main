suppressPackageStartupMessages({
  library(Seurat)           
  library(SeuratObject)     
  library(BPCells)          
  library(DESeq2)           
  library(speckle)          
  library(lme4)             
  library(emmeans)          
  library(dplyr)            
  library(purrr)            
  library(jsonlite)         
  library(readxl)           
  library(BiocParallel)     
  library(biomaRt)          
  library(AnnotationDbi)    
  library(glmGamPoi)
  library(presto)
  library(rhdf5)
  library(Matrix)
  library(stringr)
  library(RcppAnnoy)
  library(ggplot2)
  library(future)
  library(arrow)
})
options(future.globals.maxSize = Inf)
options(Seurat.object.assay.version = "v5")
options(Seurat.checkdots = FALSE)
register(SerialParam())

load("/gpfs/data/zwicklab/raz/xenium.RData")
epi.obj = subset(xenium, subset = cluster_full %in% c("0","11", "18","10", "12","4","17","1","19"))
save(epi.obj, file="/gpfs/data/zwicklab/raz/xeniumsubsets/epi_obj.RData")
dim(epi.obj)
DefaultAssay(epi.obj) <- "RNA"
dim(epi.obj)

epi.obj = ScaleData(epi.obj)
epi.obj = RunPCA(epi.obj, npcs=30, features = rownames(epi.obj))
epi.obj = RunUMAP(epi.obj, dims = 1:20, return.model = TRUE, min.dist = 0.15, n.neighbors = 35, spread = 0.4, metric = "cosine")
epi.obj = FindNeighbors(epi.obj, dims = 1:15, k.param=35, prune.SNN = 1/15)
epi.obj = FindClusters(epi.obj, resolution = 0.5, algorithm = 4)
save(epi.obj, file="/gpfs/data/zwicklab/raz/xeniumsubsets/epi_obj_clustered.RData")

DefaultAssay(epi.obj) = "RNA"
Idents(epi.obj) <- "seurat_clusters"

epi.obj$seurat_clusters = factor(epi.obj$seurat_clusters, levels = levels(epi.obj@active.ident))

getMarkers = function(seuratObj, name, dir, nmarkers = 5, height = 10, width = 20){
  markers = FindAllMarkers(seuratObj, only.pos = TRUE)
  markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)
  write.csv(markers, paste0(dir,name, "_markers.csv"))
  markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = nmarkers) %>%
    ungroup() -> top
  return(top)
}

top5 = getMarkers(epi.obj, "epi_res0.5", "/gpfs/data/zwicklab/raz/")


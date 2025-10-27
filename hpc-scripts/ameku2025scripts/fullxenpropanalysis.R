suppressPackageStartupMessages({
  library(openssl)
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
})
options(future.globals.maxSize = Inf)
options(Seurat.object.assay.version = "v5")
options(Seurat.checkdots = FALSE)
register(SerialParam())
#=============================
# You can skip this section if you have a xenium object
#=============================

xenium <- readRDS("/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium/maternalIntestineXenium.rds")
xenium <- subset(xenium, idents = 20, invert = TRUE)
DefaultAssay(xenium) <- "sketch"

read_xenium_h5 <- function(h5_path) {
  h5ls(h5_path)
  
  counts <- h5read(h5_path, "matrix/data")
  indices <- h5read(h5_path, "matrix/indices")
  indptr <- h5read(h5_path, "matrix/indptr")
  shape <- h5read(h5_path, "matrix/shape")
  
  features <- h5read(h5_path, "matrix/features/name")
  barcodes <- h5read(h5_path, "matrix/barcodes")
  
  mat <- sparseMatrix(
    i = indices + 1,
    p = indptr,
    x = as.numeric(counts),
    dims = shape,
    dimnames = list(features, barcodes)
  )
  
  return(mat)
}

sample_dirs <- list.dirs("/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium", recursive = FALSE)
matrices_list <- list()

for (sample in sample_dirs) {
  h5_file <- file.path(sample, "cell_feature_matrix.h5")
  if(file.exists(h5_file)) {
    cat("Loading", sample, "\n")
    matrices_list[[sample]] <- read_xenium_h5(h5_file)
  }
}

table(xenium$sample)
head(colnames(xenium))

merged_matrix <- NULL
for (i in seq_along(matrices_list)) {
  sample_name <- names(matrices_list)[i]
  mat <- matrices_list[[i]]
  
  colnames(mat) <- paste0(sample_name, "-", colnames(mat))
  
  if(is.null(merged_matrix)) {
    merged_matrix <- mat
  } else {
    all_genes <- union(rownames(merged_matrix), rownames(mat))
    
    merged_matrix <- merged_matrix[match(all_genes, rownames(merged_matrix)),]
    mat <- mat[match(all_genes, rownames(mat)), ]
    
    merged_matrix[is.na(merged_matrix)] <- 0
    mat[is.na(mat)] <- 0
    
    merged_matrix <- cbind(merged_matrix, mat)
  }
}

current_names <- colnames(merged_matrix)
fixed_names <- sapply(current_names, function(x) {
  parts <- strsplit(x, "-")[[1]]
  
  path_part <- parts[1]
  sample_name <- basename(path_part)
  
  barcode <- parts[2]
  
  paste(sample_name, barcode, sep = "-")
})

fixed_names <- unname(fixed_names)
colnames(merged_matrix) <- fixed_names
head(colnames(xenium))
head(colnames(merged_matrix))

if(all(grepl("-1$", colnames(xenium))) & !all(grepl("-1$", colnames(merged_matrix)))) {
  colnames(merged_matrix) <- paste0(colnames(merged_matrix), "-1")
} else if (!all(grepl("-1$", colnames(xenium))) & all(grepl("-1$", colnames(merged_matrix)))) {
  colnames(merged_matrix) <- sub("-1$", "", colnames(merged_matrix))
}

common_cells <- intersect(colnames(xenium), colnames(merged_matrix))
save(merged_matrix, file="~/piezoProject/merged_matrix.rds")

names(xenium@assays$RNA@layers)
xenium@assays$RNA@layers$counts <- merged_matrix
if("data" %in% names(xenium@assays$RNA@layers)) {
  if(!identical(xenium@assays$RNA@layers$counts, xenium@assays$RNA@layers$data)) {
    xenium@assays$RNA@layers$data <- LogNormalize(merged_matrix, verbose = TRUE)
  } else {
    xenium@assays$RNA@layers$data <- merged_matrix
  }
}

#========================

load(file="/gpfs/data/zwicklab/raz/xenium.RData")

common_cells <- intersect(colnames(xenium), colnames(merged_matrix))
save(merged_matrix, file="~/piezoProject/merged_matrix.rds")

names(xenium@assays$sketch@layers)
xenium@assays$sketch@layers$counts <- merged_matrix
if("data" %in% names(xenium@assays$sketch@layers)) {
  if(!identical(xenium@assays$sketch@layers$counts, xenium@assays$sketch@layers$data)) {
    xenium@assays$sketch@layers$data <- LogNormalize(merged_matrix, verbose = TRUE)
  } else {
    xenium@assays$sketch@layers$data <- merged_matrix
  }
}

cluster_annotations <- c(
  "0" = "Epithelial",
  "1" = "Epithelial",
  "2" = "Mesenchymal",
  "3" = "Immune", 
  "4" = "Epithelial",
  "5" = "Immune",
  "6" = "Muscle",
  "7" = "Endothelial",
  "8" = "Immune",
  "9" = "Immune",
  "10" = "Epithelial",
  "11" = "Epithelial",
  "12" = "Epithelial",
  "13" = "Enteric glia",
  "14" = "Lymphatics",
  "15" = "Mesothelial",
  "16" = "Fat",
  "17" = "Epithelial",
  "18" = "Epithelial",
  "19" = "Epithelial"
)


pca_embeddings <- Embeddings(xenium, "pca.full")
sketch_cells <- colnames(xenium[["sketch"]])
sketch_clusters <- xenium@meta.data[sketch_cells, "seurat_clusters"]
has_clusters <- !is.na(xenium@meta.data[[cluster_column]])
remaining_cells <- colnames(xenium)[!has_clusters]
sketch_pca <- pca_embeddings[sketch_cells, 1:30]
query_pca <- pca_embeddings[remaining_cells, 1:30]

index <- new(AnnoyAngular, ncol(sketch_pca))
for(i in 1:nrow(sketch_pca)) {
  index$addItem(i-1, sketch_pca[i,])
}
index$build(50)

predicted_clusters <- sapply(1:nrow(query_pca), function(i) {
  if(i %% 50000 == 0) cat("Processed", i, "cells\n")
  
  neighbors <- index$getNNsByVector(query_pca[i,],10)
  neighbor_clusters <- sketch_clusters[neighbors + 1]
  names(sort(table(neighbor_clusters), decreasing=TRUE))[1]
})

full_clusters <- rep(NA, ncol(xenium))
names(full_clusters) <- colnames(xenium)
full_clusters[sketch_cells] <- sketch_clusters
full_clusters[remaining_cells] <- predicted_clusters

xenium@meta.data$seurat_clusters <- full_clusters

xenium@meta.data$cluster_annotation <- cluster_annotations[as.character(xenium@meta.data$seurat_clusters)]

#================================
# DIM PLOTS
#=================================

DimPlot(xenium, group.by = "cluster_full", cols = "polychrome")
DimPlot(xenium, group.by = "cluster_annotation", cols = "polychrome")
DimPlot(xenium, group.by = "condition", cols = "polychrome")
DimPlot(xenium, group.by = "sample", cols = "polychrome")
FeaturePlot(xenium, features="Piezo1")

#===========================
# Cell type proportions - full dataset
#===========================

props = getTransformedProps(xenium$cluster_annotation, xenium$sample, transform="logit")

colnames(props$TransformedProps)
sample = colnames(props$TransformedProps)
condition = c("V", "V", "P", "L", "P", "L", "V", "P", "L", "V", "P", "L")
slide = c("1","4", "4", "4", "1", "1", "2", "2", "2", "3", "3", "3")

design = model.matrix(~0 + condition + slide)
design

anova.res = propeller.anova(prop.list=props,design=design, coef = c(1,2,3),
                            robust = TRUE, trend=FALSE, sort=TRUE)

write.csv(anova.res, file="~/piezoProject/anovaresultsfull.RData")

mycontr = limma::makeContrasts(conditionV-conditionP, levels = design)
ttest_V_P = propeller.ttest(props, design, contrasts = mycontr, robust=TRUE, trend= FALSE, sort = TRUE)
mycontr = limma::makeContrasts(conditionV-conditionL, levels = design)
ttest_V_L = propeller.ttest(props, design, contrasts = mycontr, robust=TRUE, trend= FALSE, sort = TRUE)
mycontr = limma::makeContrasts(conditionP-conditionL, levels = design)
ttest_P_L = propeller.ttest(props, design, contrasts = mycontr, robust=TRUE, trend=FALSE, sort = TRUE)

pvals = c(ttest_V_P$P.Value,ttest_V_L$P.Value, ttest_P_L$P.Value)
p.adj.overall = p.adjust(pvals, method = "fdr")


write.csv(ttest_V_P, file="~/piezoProject/ttest_V_P.RData")
write.csv(ttest_V_L, file="~/piezoProject/ttest_V_L.RData")
write.csv(ttest_P_L, file="~/piezoProject/ttest_P_L.RData")

sig_clusters = rownames(anova.res)[anova.res$FDR < 0.05]

props$Counts = props$Counts[rownames(props$counts) %in% sig_clusters,] 
props$TransformedProps = props$TransformedProps[rownames(props$TransformedProps) %in% sig_clusters,]
props$Proportions = props$Proportions[rownames(props$Proportions) %in% sig_clusters,]

mycontr = limma::makeContrasts(conditionV-conditionP, levels=design)
ttest_V_P = propeller.ttest(props, design, contrasts = mycontr, robust=TRUE, trend=FALSE, 
                            sort=TRUE)

mycontr = limma::makeContrasts(conditionV-conditionL, levels=design)
ttest_V_L = propeller.ttest(props, design, contrasts = mycontr, robust=TRUE, trend=FALSE, 
                            sort=TRUE)

mycontr = limma::makeContrasts(conditionP-conditionL, levels=design)
ttest_P_L = propeller.ttest(props, design, contrasts = mycontr, robust=TRUE, trend=FALSE, 
                            sort=TRUE)

pvals = c(ttest_V_P$P.Value, ttest_V_L$P.Value, ttest_P_L$P.Value)
p.adj.overall = p.adjust(pvals, method = "fdr")
ttest_V_P$FDR_overall = p.adj.overall[1:9]
ttest_V_L$FDR_overall = p.adj.overall[10:18]
ttest_P_L$FDR_overall = p.adj.overall[19:27]

write.csv(ttest_V_P, file="~/piezoProject/ttest_V_P_sigclustersfull.RData")
write.csv(ttest_V_L, file="~/piezoProject/ttest_V_L_sigclustersfull.RData")
write.csv(ttest_P_L, file="~/piezoProject/ttest_P_L_sigclustersfull.RData")


#=====================
# their plots for proportions
# ===========================
proportions = props$Proportions

proportions = reshape2::melt(proportions, id.vars = "clusters")

proportions$condition = ""
proportions$condition[proportions$sample %in% c("TA1","TA4","TA7","TA10")] = "Virgin"
proportions$condition[proportions$sample %in% c("TA2","TA5","TA8","TA11")] = "Pregnant"
proportions$condition[proportions$sample %in% c("TA3","TA6","TA9","TA12")] = "Lactating"

proportions$replicate = ""
proportions$replicate[proportions$sample %in% c("TA1","TA2","TA3")] = 1
proportions$replicate[proportions$sample %in% c("TA4","TA5","TA6")] = 2
proportions$replicate[proportions$sample %in% c("TA7","TA8","TA9")] = 3
proportions$replicate[proportions$sample %in% c("TA10","TA11","TA12")] = 4

proportions$replicate = factor(proportions$replicate, levels = c(1,2,3,4))

proportions$clusters = factor(proportions$clusters, levels = c("Epithelial", "Fat", "Endothelial", "Enteric glia", "Immune", "Lymphatics", "Mesenchymal", "Muscle", "Mesothelial"))
proportions$condition = factor(proportions$condition, levels = c("Virgin", "Pregnant","Lactating"))


pal = DiscretePalette(9, palette = "polychrome")

ggplot(proportions, aes(fill=clusters, y=value, x=replicate)) + 
  geom_bar(position="stack", stat="identity") + facet_wrap(~condition) + scale_fill_manual(values = pal) + theme_classic() +
  theme(text = element_text(size = 20))


save(xenium, file="/gpfs/data/zwicklab/raz/xenium.rds")

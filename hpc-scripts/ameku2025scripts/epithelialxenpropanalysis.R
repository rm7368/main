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
})
options(future.globals.maxSize = Inf)
options(Seurat.object.assay.version = "v5")
options(Seurat.checkdots = FALSE)
register(SerialParam())


epithelial_annotations <- c(
  "0" = "Isthmus progenitor",
  "1" = "Bottom-villus EC",
  "2" = "Top (+mid)-villus EC",
  "3" = "Lgr5+ CBC/Paneth",
  "4" = "Mid (+top)-villus EC",
  "5" = "Goblet",
  "6" = "Mid-villus EC, anterior",
  "7" = "ECM/mesenchymal",
  "8" = "Other EC 1",
  "9" = "Enteroendocrine",
  "10" = "Other EC 2",
  "11" = "Microfold",
  "12" = "Tuft",
  "13" = "Enteroendocrine progenitor",
  "14" = "Other EC 3",
  "15" = "Other EC 4",
  "16" = "Telocyte, top-villus",
  "17" = "Other EC 5",
  "18" = "Telocyte"
)

  
#=======================
# Object subsetting
#=======================

load("/gpfs/data/zwicklab/raz/xenium.RData")
epi.obj = subset(xenium, subset = xenium@meta.data$epi_cluster_annotation != "NA")
Idents(epi.obj) <- epi.obj@meta.data$epi_cluster_annotation
save(epi.obj, file="/gpfs/data/zwicklab/raz/xeniumsubsets/epi_obj.RData")
dim(epi.obj)
DefaultAssay(epi.obj) <- "RNA"
dim(epi.obj)

#=====================

epi.obj = ScaleData(epi.obj)
epi.obj = RunPCA(epi.obj, npcs=30, features = rownames(epi.obj))
epi.obj = RunUMAP(epi.obj, dims = 1:20, return.model = TRUE, min.dist = 0.15, n.neighbors = 35, spread = 0.4, metric = "cosine")
epi.obj = FindNeighbors(epi.obj, dims = 1:20, k.param=20)
epi.obj = FindClusters(epi.obj, resolution = 0.5, algorithm = 4)
save(epi.obj, file="/gpfs/data/zwicklab/raz/xeniumsubsets/epi_obj_clustered.RData")

DimPlot(epi.obj, cols = "polychrome")
ggsave("/gpfs/data/zwicklab/raz/epiUMAP_res0.5.png")

DimPlot(epi.obj, group.by = "condition")
ggsave(paste0("/gpfs/data/zwicklab/raz/epiUMAP_condition.png"))

DimPlot(epi.obj, group.by = "sample", cols = "polychrome")
ggsave(paste0("/gpfs/data/zwicklab/raz/epiUMAP_sample.png"))

clusters_to_keep = c("1","2","3","4","5","6","7",
                     "8","9","10","11","12","13","14",
                     "15","16", "17","18")

epi.obj <- subset(epi.obj, subset = epi.obj$seurat_clusters %in% clusters_to_keep)

DefaultAssay(epi.obj) = "RNA"
Idents(epi.obj) <- "seurat_clusters"

save(epi.obj, file="/gpfs/data/zwicklab/raz/xeniumsubsets/epi_obj_clustered.RData")

res = table(epi.obj$condition, epi.obj$seurat_clusters)
write.csv(res, "/gpfs/data/zwicklab/raz/epi_numbers_res0.5.csv", row.names = F, quote = F)

epi.obj$seurat_clusters = factor(epi.obj$seurat_clusters, levels = levels(epi.obj@active.ident))
p = DimPlot(epi.obj, cols = "polychrome") + ylab("umap_2") + xlab("umap_1")
p + ggtitle("")

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

xenium@meta.data$epi_cluster_annotation <- epithelial_annotations[as.character(xenium@meta.data$epi_res0.5)]

epi.obj@meta.data$cluster_annotation <- epithelial_annotations[as.character(epi.obj@meta.data$seurat_clusters)]

DimPlot(epi.obj, group.by = "cluster_annotation", cols = "polychrome")


save(epi.obj, file="/gpfs/data/zwicklab/raz/xeniumsubsets/epi_obj_annotated.RData")

FeaturePlot(epi.obj, features = "Piezo1")


#======================
# Proportions
#=====================

props = getTransformedProps(epi.obj$epi_cluster_annotation, epi.obj$sample, transform="logit")
colnames(props$TransformedProps)
sample = colnames(props$TransformedProps)
condition = c("V", "V", "P", "L", "P",
              "L", "V", "P", "L", "V",
              "P", "L")
slide = c("1","4","4","4","1","1","2","2","2","3","3","3")
design = model.matrix(~ 0 + condition + slide)
design

anova.res = propeller.anova(prop.list = props, design=design, coef = c(1,2,3),
                            robust = TRUE, trend=FALSE, sort=TRUE)
write.csv(anova.res, "~/piezoProject/epipropeller/anovaresultsepi.csv")

mycontr = limma::makeContrasts(conditionV-conditionP, levels=design)          
ttest_V_P = propeller.ttest(props, design, contrasts = mycontr, robust = TRUE, trend = FALSE,
                            sort = TRUE)          
mycontr = limma::makeContrasts(conditionV-conditionL, levels=design)          
ttest_V_L = propeller.ttest(props, design, contrasts = mycontr, robust = TRUE, trend = FALSE, 
                            sort = TRUE)          

mycontr = limma::makeContrasts(conditionP-conditionL, levels=design)          
ttest_P_L = propeller.ttest(props, design, contrasts = mycontr, robust = TRUE, trend = FALSE, 
                            sort = TRUE)       

pvals = c(ttest_V_P$P.Value, ttest_V_L$P.Value, ttest_P_L$P.Value)
p.adj.overall = p.adjust(pvals, method = "fdr")

ttest_V_P$FDR_overall = p.adj.overall[1:19]
ttest_V_L$FDR_overall = p.adj.overall[20:38]
ttest_P_L$FDR_overall = p.adj.overall[39:57]

write.csv(ttest_V_P, "~/piezoProject/epipropeller/propeller_results_V_P.csv")
write.csv(ttest_V_L, "~/piezoProject/epipropeller/propeller_results_V_L.csv")
write.csv(ttest_P_L, "~/piezoProject/epipropeller/propeller_results_P_L.csv")

sig_clusters = rownames(anova.res)[anova.res$FDR < 0.05]

#=============
# plotting
#=============

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

proportions$clusters = factor(proportions$clusters, levels = unique(epithelial_annotations))
proportions$condition = factor(proportions$condition, levels = c("Virgin", "Pregnant","Lactating"))


pal = DiscretePalette(19, palette = "polychrome")

ggplot(proportions, aes(fill=clusters, y=value, x=replicate)) +
  geom_bar(position="stack", stat="identity") + facet_wrap(~condition) + scale_fill_manual(values=pal) + theme_classic() +
  theme(text = element_text(size = 20))
ggsave("~/piezoProject/epipropeller/cell_type_proportions_epi.pdf")

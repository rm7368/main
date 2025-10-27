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
  library(glmGamPoi)
  library(EnhancedVolcano)
})

# High-performance computing setup (sequential processing)
options(future.globals.maxSize = Inf)
options(Seurat.object.assay.version = "v5")
options(Seurat.checkdots = FALSE)

# BiocParallel setup for DESeq2 (this is separate from future)
register(SerialParam())

load("/gpfs/data/zwicklab/raz/xeniumsubsets/epi_obj_annotated.RData")

DefaultAssay(epi.obj) = "RNA"

pseudo_bulk = AggregateExpression(epi.obj, assays = "RNA", group.by = c("condition", "sample", "epi_cluster_annotation"))

pseudo_bulk_meta = lapply(colnames(pseudo_bulk$RNA),str_split_1, "_")
pseudo_bulk_meta = do.call(rbind.data.frame,pseudo_bulk_meta)
colnames(pseudo_bulk_meta) = c("condition", "sample", "cluster")
rownames(pseudo_bulk_meta) = colnames(pseudo_bulk$RNA)

pseudo_bulk_meta$cluster[pseudo_bulk_meta$cluster == "ECM/mesenchymal"] <- "ECMmesenchymal"

pseudo_bulk_meta$slide = ""
pseudo_bulk_meta$slide[pseudo_bulk_meta$sample %in% c("TA1","TA2","TA3")] = "S1"
pseudo_bulk_meta$slide[pseudo_bulk_meta$sample %in% c("TA4","TA5","TA6")] = "S2"
pseudo_bulk_meta$slide[pseudo_bulk_meta$sample %in% c("TA7","TA8","TA9")] = "S3"
pseudo_bulk_meta$slide[pseudo_bulk_meta$sample %in% c("TA10","TA11","TA12")] = "S4"

run_deseq2 = function(condition1, condition2){
  results = list()
  for (cluster in unique(pseudo_bulk_meta$cluster)){
    
    selected = pseudo_bulk$RNA[,(pseudo_bulk_meta$cluster == cluster) & (pseudo_bulk_meta$condition %in% c(condition1,condition2)) ]
    selected_meta = pseudo_bulk_meta[(pseudo_bulk_meta$cluster == cluster) & (pseudo_bulk_meta$condition %in% c(condition1,condition2)),]
    selected_meta$condition = factor(selected_meta$condition, levels = c(condition1,condition2))
    
    dds = DESeqDataSetFromMatrix(countData = selected,
                                 colData = selected_meta,
                                 design = ~slide + condition)
    
    dds = DESeq(dds)
    res = results(dds)
    
    res$cluster = cluster
    res$gene = rownames(res)
    
    res = res[order(res$pvalue),]
    res$gene = rownames(res)
    
    results[[cluster]] = res
  }
  results = do.call(rbind,results)
  
  # add adjusted p value that takes into account testing was done for multiple clusters.
  results$padj_overall = p.adjust(results$pvalue,method = "BH")
  return(results)
}

results_v_p = run_deseq2("V","P")
results_v_l = run_deseq2("V","L")
results_p_l = run_deseq2("P","L")

write.csv(results_p_l, "/gpfs/home/rm7368/piezoProject/epipseudobulk/results_p_l.csv")
write.csv(results_v_p, "~/piezoProject/epipseudobulk/results_v_p.csv")
write.csv(results_v_l, "~/piezoProject/epipseudobulk/results_v_l.csv")

results_p_l = results_p_l[!(is.na(results_p_l$padj_overall)),]
results_v_p = results_v_p[!(is.na(results_v_p$padj_overall)),]
results_v_l = results_v_l[!(is.na(results_v_l$padj_overall)),]

write.csv(results_p_l[results_p_l$padj_overall < 0.05,], "~/piezoProject/epipseudobulk/results_p_l_sig.csv")
write.csv(results_v_p[results_v_p$padj_overall < 0.05,], "~/piezoProject/epipseudobulk/results_v_p_sig.csv")
write.csv(results_v_l[results_v_l$padj_overall < 0.05,], "~/piezoProject/epipseudobulk/results_v_l_sig.csv")

volcanoDir = "~/piezoProject/epipseduobulk/volcanoplots/"
dir.create(volcanoDir, recursive=TRUE)

makeVolcano = function(results, cluster, name){
  
  if (cluster == "ECMmesenchymal") {
    res = results[results$cluster == "ECM/mesenchymal",]
    name = name
    plot_name = paste0(name, " cluster ECMmesenchymal")
  } else {
    res = results[results$cluster == cluster,]
    name = name
    plot_name = paste0(name," cluster ", cluster)
  }
  
  # Clean the data - remove rows with NA or infinite values
  res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & 
               is.finite(res$log2FoldChange) & is.finite(res$padj), ]
  
  # Check if we still have data after cleaning
  if (nrow(res) == 0) {
    cat("No valid data after cleaning for cluster:", cluster, "in", name, "\n")
    return(NULL)
  }
  
  # Replace any remaining problematic values
  res$padj[res$padj == 0] <- .Machine$double.xmin  # Replace 0 with very small number
  
  sig_genes <- res$gene[res$padj < 0.05]
  if (length(sig_genes) == 0) {
    sig_genes <- character(0)
  }
  
  EnhancedVolcano(res,
                  title = plot_name,
                  subtitle = "",
                  lab = res$gene,
                  x = 'log2FoldChange',
                  y = 'padj',
                  xlab = bquote(~Log[2]~ 'fold change'),
                  ylab = bquote(~-Log[10]~adjusted~italic(P)),
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  selectLab = sig_genes,
                  pointSize = 4.0,
                  labSize = 6.0,
                  legendPosition = 'right',
                  legendLabSize = 12,
                  legendIconSize = 4.0,
                  drawConnectors = ifelse(length(sig_genes) > 0, TRUE, FALSE),
                  widthConnectors = 0.75,
                  max.overlaps = Inf)
  
  safe_cluster_name <- gsub("/", "_", cluster)
  filename <- paste0(volcanoDir, name, " cluster ", safe_cluster_name, ".pdf")
  ggsave(filename=filename, width = 8, height = 6)
}

for (cluster in unique(pseudo_bulk_meta$cluster)){
  makeVolcano(results_p_l, cluster,  "Lactating vs pregnant")
  makeVolcano(results_v_p, cluster,  "Pregnant vs virgin")
  makeVolcano(results_v_l, cluster,  "Lactating vs virgin")
}

for (cluster in unique(pseudo_bulk_meta$cluster)) {
  cat("Cluster:", cluster, "\n")
  cat("  results_p_l rows:", nrow(results_p_l[results_p_l$cluster == cluster,]), "\n")
  cat("  results_v_p rows:", nrow(results_v_p[results_v_p$cluster == cluster,]), "\n")
  cat("  results_v_l rows:", nrow(results_v_l[results_v_l$cluster == cluster,]), "\n")
}

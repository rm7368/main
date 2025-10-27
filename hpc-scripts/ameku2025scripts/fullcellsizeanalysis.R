library(reticulate)
library(arrow)
library(ggplot2)
library(interp)
library(rdist)
library(Matrix)
library(dplyr)
library(future)
library(BPCells)
options(future.globals.maxSize = Inf)
library(data.table)
library(lme4)
library(emmeans)
library(car)
library(lmerTest)

load("/gpfs/data/zwicklab/raz/xenium.RData")

sample_names = unique(xenium$sample)

sample_paths = c("/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium/TA1",
                 "/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium/TA2",
                 "/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium/TA3",
                 "/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium/TA4",
                 "/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium/TA5",
                 "/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium/TA6",
                 "/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium/TA7",
                 "/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium/TA8",
                 "/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium/TA9",
                 "/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium/TA10",
                 "/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium/TA11",
                 "/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium/TA12")

names(sample_paths) <- sample_names

meta = c("V", "P", "L",
         "V", "P", "L",
         "V", "P", "L",
         "V", "P", "L")
names(meta) <- names(sample_paths)

cells = read_parquet(paste0(sample_paths[1], "/cells.parquet"))
cells$id = paste0("TA1-", cells$cell_id)
cells$sample = "TA1"

for (sample in sample_names[!is.name("TA1")]) {
  cellsTemp = read_parquet(paste0(sample_paths[[sample]], "/cells.parquet"))
  cellsTemp$id = paste0(sample, "-", cellsTemp$cell_id)
  cellsTemp$sample = sample
  cells = rbind(cells, cellsTemp)
}


directory = "/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium/"
write.csv(xenium@meta.data,paste0(directory,"nuclei_meta.csv"))

nuclei_meta = read.csv(paste0(directory, "nuclei_meta.csv"))
rownames(nuclei_meta) = nuclei_meta$cluster_annotation

cells_unique <- cells[!duplicated(cells$id), ]
cat("Cells after removing duplicates:", nrow(cells_unique), "\n")

cells_selected <- cells_unique[cells_unique$id %in% nuclei_meta$X, ]
cat("Final cells_selected rows:", nrow(cells_selected), "\n")

cells_selected$cluster = nuclei_meta$cluster_annotation

cell_stats_dir = paste0(directory,"cell_stats/")
dir.create(cell_stats_dir)

cells_selected$condition = ""
cells_selected$condition[cells_selected$sample %in% c("TA1","TA4","TA7","TA10")] = "Virgin"
cells_selected$condition[cells_selected$sample %in% c("TA2","TA5","TA8","TA11")] = "Pregnant"
cells_selected$condition[cells_selected$sample %in% c("TA3","TA6","TA9","TA12")] = "Lactating"

cells_selected$replicate = ""
cells_selected$replicate[cells_selected$sample %in% c("TA1","TA2","TA3")] = 1
cells_selected$replicate[cells_selected$sample %in% c("TA4","TA5","TA6")] = 2
cells_selected$replicate[cells_selected$sample %in% c("TA7","TA8","TA9")] = 3
cells_selected$replicate[cells_selected$sample %in% c("TA10","TA11","TA12")] = 4

cells_selected$condition = factor(cells_selected$condition, levels = c("Virgin","Pregnant","Lactating"))
cells_selected$replicate = factor(cells_selected$replicate)

cells_selected$cell_area_log = log(cells_selected$cell_area)
cells_selected$nucleus_area_log = log(cells_selected$nucleus_area)

#=================
# Full Xenium dataset cell/nucleus areas
#================

p = ggplot(data = cells_selected, aes(x = replicate, y = cell_area, color = condition)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  facet_wrap(~condition) +
  theme_classic() +
  theme(text = element_text(size=20))

print(p)
ggsave(paste0(cell_stats_dir, "cell_area_epi.pdf"))  

p = ggplot(data = cells_selected, aes(x = replicate, y = cell_area_log, color = condition)) + 
  geom_boxplot(width=0.5) + 
  geom_violin(alpha=0.7) +
  facet_wrap(~condition) + 
  theme_classic() +
  theme(text = element_text(size = 20))

print(p)
ggsave(paste0(cell_stats_dir, "cell_area_log_epi.pdf"))

p = ggplot(data = cells_selected, aes(x = replicate, y = nucleus_area, color = condition)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  facet_wrap(~condition) +
  theme_classic() +
  theme(text = element_text(size=20))

print(p)
ggsave(paste0(cell_stats_dir, "nucleus_area_full.pdf"))  

p = ggplot(data = cells_selected, aes(x = replicate, y = nucleus_area_log, color = condition)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  facet_wrap(~condition) +
  theme_classic() +
  theme(text = element_text(size=20))

print(p)
ggsave(paste0(cell_stats_dir, "nucleus_area_log_epi.pdf"))

#=========================
# Clustered dataset cell and nucleus areas
#==========================

for (i in cells_selected$cluster) {
  print(i)
  
  selected = cells_selected[cells_selected$cluster == i,]
  
  p = ggplot(data = cells_selected, aes(x = replicate, y = cell_area, color = condition)) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    facet_wrap(~condition) +
    theme_classic() +
    theme(text = element_text(size=20))
  
  print(p)
  ggsave(paste0(cell_stats_dir, "cell_area_cluster_", i, ".pdf"))  
  
  p = ggplot(data = cells_selected, aes(x = replicate, y = cell_area_log, color = condition)) + 
    geom_boxplot(width=0.5) + 
    geom_violin(alpha=0.7) +
    facet_wrap(~condition) + 
    theme_classic() +
    theme(text = element_text(size = 20))
  
  print(p)
  ggsave(paste0(cell_stats_dir, "cell_area_log_cluster_", i, ".pdf"))
  
  p = ggplot(data = cells_selected, aes(x = replicate, y = nucleus_area, color = condition)) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    facet_wrap(~condition) +
    theme_classic() +
    theme(text = element_text(size=20))
  
  print(p)
  ggsave(paste0(cell_stats_dir, "nucleus_area_cluster_", i, ".pdf"))  
  
  p = ggplot(data = cells_selected, aes(x = replicate, y = nucleus_area_log, color = condition)) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    facet_wrap(~condition) +
    theme_classic() +
    theme(text = element_text(size=20))
  
  print(p)
  ggsave(paste0(cell_stats_dir, "nucleus_area_log_cluster_", i, ".pdf"))
  
}


  
  
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
  library(ggplot2)
  library(future)
  library(ggridges)
  library(ggpubr)
  library(rstatix)
  library(patchwork)
  library(arrow)
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

load("/gpfs/data/zwicklab/raz/xeniumsubsets/epi_obj.RData")

Idents(epi.obj) = "epi_cluster_annotation"

####################
# Create statistics table and dot plots
#####################

results <- data.frame()
all_cell_types = epithelial_annotations

for(i in 1:length(all_cell_types)) {
  ct <- all_cell_types[i]
  
  # Get cells of this type
  cells_of_type <- which(epi.obj$epi_cluster_annotation == ct)
  
  # Only proceed if enough cells
  if(length(cells_of_type) > 10) {
    
    # Get data
    expr_data <- LayerData(epi.obj, assay = "RNA", layer = "data")["Prlr", cells_of_type]
    conditions <- epi.obj$condition[cells_of_type]
    
    # Create dataframe
    df <- data.frame(
      expression = as.numeric(expr_data),
      condition = conditions,
      cell_type = ct
    )
    
    # Test if multiple conditions
    if(length(unique(df$condition)) > 1) {
      
      # Statistical test
      stat_test <- df %>%
        wilcox_test(expression ~ condition) %>%
        add_significance("p") %>%
        mutate(cell_type = ct)
      
      # Add to results
      results <- rbind(results, stat_test)
      
      print(paste("Completed:", ct, "- found", nrow(stat_test), "comparisons"))
    } else {
      print(paste("Skipped (single condition):", ct))
    }
  } else {
    print(paste("Skipped (too few cells):", ct))
  }
}

print(paste("Total comparisons:", nrow(results)))

# Apply multiple testing correction
results$p.adj.global <- p.adjust(results$p, method = "BH")

# Add global significance levels
results$significance_global <- ifelse(results$p.adj.global < 0.001, "***",
                                      ifelse(results$p.adj.global < 0.01, "**",
                                             ifelse(results$p.adj.global < 0.05, "*", "ns")))

# View all results
print("All statistical results:")
print(results)

# Filter for significant results
significant_results <- results %>%
  filter(p.adj.global < 0.05) %>%
  arrange(p.adj.global)

print("Globally significant results after MTC:")
print(significant_results)
# Create summary data for visualizations
summary_data <- data.frame()

for(ct in unique(epi.obj$epi_cluster_annotation)) {
  for(cond in unique(epi.obj$condition)) {
    cells <- which(epi.obj$epi_cluster_annotation == ct & epi.obj$condition == cond)
    
    if(length(cells) > 5) {
      expr <- LayerData(epi.obj, assay = "RNA", layer = "data")["Prlr", cells]
      
      summary_data <- rbind(summary_data, data.frame(
        cell_type = ct,
        condition = cond,
        mean_expr = mean(expr),
        pct_expr = sum(expr > 0) / length(expr) * 100,
        n_cells = length(cells)
      ))
    }
  }
}

print("Summary data created:")
head(summary_data)

print(significant_results)
# Look at the significant results in detail
print("Significant comparisons:")
print(significant_results)

# Which cell types show significant changes?
sig_cell_types <- unique(significant_results$cell_type)
print("Cell types with significant changes:")
print(sig_cell_types)

# What types of comparisons are significant?
print("Comparison breakdown:")
table(paste(significant_results$group1, "vs", significant_results$group2))

library(patchwork)

# Create DotPlots for each condition
p1 <- DotPlot(virgin_obj, features = "Prlr", group.by = "epi_cluster_annotation", 
              cols = c("lightgrey", "red")) + 
  ggtitle("Virgin") + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- DotPlot(pregnant_obj, features = "Prlr", group.by = "epi_cluster_annotation",
              cols = c("lightgrey", "red")) + 
  ggtitle("Pregnant") + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- DotPlot(lactating_obj, features = "Prlr", group.by = "epi_cluster_annotation",
              cols = c("lightgrey", "red")) + 
  ggtitle("Lactating") + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine DotPlots

sig_cell_types_df <- data.frame(
  cell_type = sig_cell_types,
  has_significance = TRUE
)

summary_data_annotated <- merge(summary_data, sig_cell_types_df, all.x = TRUE)
summary_data_annotated$has_significance[is.na(summary_data_annotated$has_significance)] <- FALSE
summary_data_annotated$condition <- factor(summary_data_annotated$condition, 
                                           levels = c("V", "P", "L"))


# Clean combined DotPlot (standard approach)
combined_dotplot <- ggplot(summary_data_annotated, aes(x = condition, y = cell_type)) +
  geom_point(aes(size = pct_expr, color = mean_expr)) +
  scale_color_gradient(low = "lightgrey", high = "red", name = "Mean\nExpression") +
  scale_size_continuous(range = c(1, 10), name = "% Expressing") +
  scale_alpha_manual(values = c(0.5, 1.0), name = "Significant", guide = "none") +
  theme_bw() +
  labs(x = "Condition", y = "Cell Type", title = "Prlr Expression Across Conditions") +
  theme(axis.text.x = element_text(angle = 0),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))

print(combined_dotplot)
ggsave(paste0("/gpfs/data/zwicklab/raz/xeniumctplots/", features, " Expression_CellTypes.pdf"),
       height = 12, width = 7)

# Create a clean statistical summary table
stats_table <- significant_results %>%
  select(cell_type, group1, group2, p, p.adj.global, significance_global) %>%
  arrange(p.adj.global) %>%
  mutate(
    comparison = paste(group1, "vs", group2),
    p_value = round(p, 4),
    adj_p_value = round(p.adj.global, 4)
  ) %>%
  select(cell_type, comparison, p_value, adj_p_value, significance_global)

# Rename columns for clarity
colnames(stats_table) <- c("Cell_Type", "Comparison", "P_Value", "Adjusted_P_Value", "Significance")

print("Statistical Summary Table:")
print(stats_table)

# Save the table
write.csv(stats_table, paste0("/gpfs/data/zwicklab/raz/xeniumctplots/", features, "_ctexpression_statistical_comparisons.csv"), row.names = FALSE)




####################
# Add FOVs to epi.obj
####################

base_path <- "/gpfs/data/zwicklab/raz/maternalIntestine/maternalIntestineXenium"
add_spatial_to_sample <- function(sample_name, base_path, seurat_obj) {
  sample_path <- file.path(base_path, sample_name)
  cells_data <- read_parquet(file.path(sample_path, "cells.parquet"))
  coords_df <- data.frame(
    x = cells_data$x_centroid,
    y = cells_data$y_centroid,
    cell = paste0(sample_name, "-", cells_data$cell_id),
    stringsAsFactors = FALSE
  )
  coords_df = coords_df[coords_df$cell %in% colnames(epi.obj),]
  centroid_data = list(
    "centroids" = CreateCentroids(coords_df)
  )
  fov <- CreateFOV(
    coords = centroid_data,
    assay = "RNA",
    key = paste0(sample_name, "_"),
  )
  
  return(fov)
}

fov_TA2 <- add_spatial_to_sample("TA2", base_path, epi.obj)
epi.obj[["TA2"]] <- fov_TA2

for (i in range(7,8)) {
  fov <- add_spatial_to_sample(paste0("TA",i), base_path, epi.obj)
  epi.obj[[paste0("TA",i)]] <- fov
}

save(epi.obj, file="/gpfs/data/zwicklab/raz/xeniumsubsets/epi_obj.RData")



####################
# Spatial feature plots
####################

receptors = c("Esr1", "Esr2", "Esrrg", "Ar", "Pgr", "Prlr")

for (gene in receptors) {
  for (sample in unique(epi.obj@meta.data$sample)) {
    dir.create(paste0("/gpfs/data/zwicklab/raz/xeniumctplots/", sample, "/virgin_gutroll_centroid_plots"), recursive = TRUE)
    dir.create(paste0("/gpfs/data/zwicklab/raz/xeniumctplots/", sample, "/pregnant_gutroll_centroid_plots"), recursive = TRUE)
    dir.create(paste0("/gpfs/data/zwicklab/raz/xeniumctplots/", sample, "/lactating_gutroll_centroid_plots"), recursive = TRUE)
    if (sample %in% c("TA1", "TA4", "TA7", "TA10")) {
      p <- ImageFeaturePlot(epi.obj, features = gene, fov = sample)
      p
      ggsave(filename=paste0("/gpfs/data/zwicklab/raz/", sample, "/virgin_gutroll_centroid_plots/", gene, "_expression_centroids.pdf"), height = 20, width = 20)
    } else if (sample %in% c("TA2", "TA5", "TA8", "TA11")) {
      p <- ImageFeaturePlot(epi.obj, features = gene, fov = sample)
      p
      ggsave(filename=paste0("/gpfs/data/zwicklab/raz/", sample, "/pregnant_gutroll_centroid_plots/", gene, "_expression_centroids.pdf"), height = 20, width = 20)
    } else {
      p <- ImageFeaturePlot(epi.obj, features = gene, fov = sample)
      p
      ggsave(filename=paste0("/gpfs/data/zwicklab/raz/", sample, "/lactating_gutroll_centroid_plots/", gene, "_expression_centroids.pdf"), height = 20, width = 20)
    }
  }
}


###########DEBUG AREA#########################
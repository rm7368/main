library(Seurat)
library(caret)
seurat_obj <- readRDS("/gpfs/home/rm7368/AllCells_1.11.22.decompressed")

set.seed(42)

# Domain distribution
table(seurat_obj$MetabolicDomains)

# Cell type distribution
table(seurat_obj$cell.type)

# Cross-tabulation: which cell types are in which domains?
table(seurat_obj$cell.type, seurat_obj$MetabolicDomains)

# Visualize
DimPlot(seurat_obj, group.by = "MetabolicDomains")
DimPlot(seurat_obj, group.by = "cell.type", label = TRUE)
DimPlot(seurat_obj, group.by = "MetabolicDomains", split.by = "cell.type", ncol = 3)

# Get gene names for SCT transform
hvg_3000 <- VariableFeatures(seurat_obj)
length(hvg_3000)

# Subset RNA assay for the 3000 most variable genes
expr_matrix <- as.matrix(seurat_obj@assays$RNA@data[hvg_3000,])
dim(expr_matrix)
range(expr_matrix)
summary(as.vector(expr_matrix[1:100,1:100]))

# Extract metadata
metadata <- data.frame(
  cell_id = colnames(expr_matrix),
  mouse = seurat_obj$mouse,
  segment = seurat_obj$segment,
  cell_type = seurat_obj$cell.type,
  domain = seurat_obj$MetabolicDomains,
  stringsAsFactors = FALSE
)

# First, log-normalize raw RNA counts
normalize_counts <- function(counts_matrix) {
  lib_size <- colSums(counts_matrix)
  cpm <- t(t(counts_matrix)/lib_size * 10000)
  log_norm <- log1p(cpm)
  return(log_norm)
}

expr_matrix_norm <- normalize_counts(expr_matrix)

# Verify it worked
cat("Before normalization:\n")
print(summary(as.vector(expr_matrix[1:1000, 1:100])))
cat("\nAfter normalization:\n")
print(summary(as.vector(expr_matrix_norm[1:1000, 1:100])))

# Check max value (should be reasonable now)
cat("\nMax value:", max(expr_matrix_norm), "\n")

# Visualize
hist(as.vector(expr_matrix_norm[, 1:1000]), breaks = 100, 
     main = "Log-normalized expression distribution",
     xlab = "log1p(CPM)")

# Update to use normalized data
expr_matrix <- expr_matrix_norm

# Clean up
rm(expr_matrix_norm)

# Caret stratification

metadata$strata <- paste(metadata$cell_type, metadata$domain, sep = "_")
train_idx_random <- createDataPartition(metadata$strata, p = 0.8, list = FALSE)
train_data_random <- list(
  expr = expr_matrix[, train_idx_random],
  cell_type = metadata$cell_type[train_idx_random],
  domain = metadata$domain[train_idx_random],
  mouse = metadata$mouse[train_idx_random]
)

test_data_random <- list(
  expr = expr_matrix[,-train_idx_random],
  cell_type = metadata$cell_type[-train_idx_random],
  domain = metadata$domain[-train_idx_random],
  mouse = metadata$mouse[-train_idx_random]
)

# One-mouse-out splits
train_data_mouseA <- list(
  expr = expr_matrix[,metadata$mouse == "A"],
  cell_type = metadata$cell_type[metadata$mouse == "A"],
  domain = metadata$domain[metadata$mouse == "A"]
)

test_data_mouseB <- list(
  expr = expr_matrix[,metadata$mouse == "B"],
  cell_type = metadata$cell_type[metadata$mouse == "B"],
  domain = metadata$domain[metadata$mouse == "B"]
)

train_data_mouseB <- list(
  expr = expr_matrix[, metadata$mouse == "B"],
  cell_type = metadata$cell_type[metadata$mouse == "B"],
  domain = metadata$domain[metadata$mouse == "B"]
)

test_data_mouseA <- list(
  expr = expr_matrix[, metadata$mouse == "A"],
  cell_type = metadata$cell_type[metadata$mouse == "A"],
  domain = metadata$domain[metadata$mouse == "A"]
)

# Verify splits
cat("Strategy A (Random split):\n")
cat("Train:", ncol(train_data_random$expr), "cells\n")
cat("Test:", ncol(test_data_random$expr), "cells\n")
table(train_data_random$cell_type, train_data_random$domain)

cat("\nStrategy B (Mouse A train):\n")
cat("Train:", ncol(train_data_mouseA$expr), "cells\n")
cat("Test:", ncol(test_data_mouseB$expr), "cells\n")

# Save everything for later use
save(train_data_random, test_data_random,
     train_data_mouseA, test_data_mouseA,
     train_data_mouseB, test_data_mouseB,
     metadata, expr_matrix,
     file = "intestinal_domain_data.RData")

cat("\nData saved! Ready for modeling.\n")

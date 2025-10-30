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

# ============Enterocyte subsetting==========
enterocyte_types <- c("enterocyte", "enterocyte.mature")
enterocyte_mask <- metadata$cell_type %in% enterocyte_types
print(table(metadata$cell_type))
print(table(metadata$domain[enterocyte_mask]))
print(table(metadata$mouse[enterocyte_mask]))

# ============Enterocyte filtering=========

metadata_enterocyte <- metadata[enterocyte_mask,]
expr_matrix_enterocyte <- expr_matrix[,enterocyte_mask]
print(table(metadata_enterocyte$cell_type))

# =============One-mouse-out splits==============
train_data_enterocyte_mouseA <- list(
  expr = expr_matrix_enterocyte[,metadata_enterocyte$mouse == "A"],
  cell_type = metadata_enterocyte$cell_type[metadata_enterocyte$mouse == "A"],
  domain = metadata_enterocyte$domain[metadata_enterocyte$mouse == "A"]
)

test_data_enterocyte_mouseB <- list(
  expr = expr_matrix_enterocyte[,metadata_enterocyte$mouse == "B"],
  cell_type = metadata_enterocyte$cell_type[metadata_enterocyte$mouse == "B"],
  domain = metadata_enterocyte$domain[metadata_enterocyte$mouse == "B"]
)

# =============Verify splits=============
cat("\n===== MOUSE-BASED SPLIT =====\n")
cat("Training cells (Mouse A):", ncol(train_data_enterocyte_mouseA$expr), "\n")
cat("Test cells (Mouse B):", ncol(test_data_enterocyte_mouseB$expr), "\n")

cat("\nTraining domain distribution (Mouse A):\n")
print(table(train_data_enterocyte_mouseA$domain))

cat("\nTest domain distribution (Mouse B):\n")
print(table(test_data_enterocyte_mouseB$domain))

cat("\nTraining enterocyte subtypes (Mouse A):\n")
print(table(train_data_enterocyte_mouseA$cell_type))

cat("\nTest enterocyte subtypes (Mouse B):\n")
print(table(test_data_enterocyte_mouseB$cell_type))

# Check class balance
cat("\n===== CLASS BALANCE CHECK =====\n")
cat("Training proportions:\n")
print(prop.table(table(train_data_enterocyte_mouseA$domain)))
cat("\nTest proportions:\n")
print(prop.table(table(test_data_enterocyte_mouseB$domain)))

# ===== SAVE R DATA =====

save(train_data_enterocyte_mouseA, test_data_enterocyte_mouseB,
     metadata_enterocyte, expr_matrix_enterocyte,
     file = "enterocyte_domain_data_mouseAB.RData")

cat("\n\nData saved to enterocyte_domain_data_mouseAB.RData\n")


#==========Export training data=============
write.csv(t(train_data_enterocyte_mouseA$expr),
          "pydata/train_mouseA_expr.csv",
          row.names = FALSE)
train_labels <- data.frame(
  domain = train_data_enterocyte_mouseA$domain,
  cell_type = train_data_enterocyte_mouseA$cell_type,
  stringsAsFactors = FALSE
)
write.csv(train_labels, "pydata/train_mouseA_labels.csv", row.names = FALSE)
write.csv(t(test_data_enterocyte_mouseB$expr),
          "pydata/test_mouseB_expr.csv",
          row.names = FALSE)
test_labels <- data.frame(
  domain = test_data_enterocyte_mouseB$domain,
  cell_type = test_data_enterocyte_mouseB$cell_type,
  stringsAsFactors = FALSE
)
write.csv(test_labels, "pydata/test_mouseB_labels.csv", row.names = FALSE)
write.csv(data.frame(gene = rownames(train_data_enterocyte_mouseA$expr)),
          "pydata/gene_names.csv", row.names = FALSE)
domain_mapping <- data.frame(
  domain = c("Domain_A", "Domain_B", "Domain_C", "Domain_D", "Domain_E"),
  label = 0:4,
  stringsAsFactors = FALSE
)
write.csv(domain_mapping, "pydata/domain_mapping.csv", row.names = FALSE)
files <- list.files("pydata", pattern = "*.csv")
for(f in files){
  size_mb <- file.info(file.path("pydata",f))$size / 1024^2
  cat(sprintf("%-40s %.2f MB\n", f, size_mb))
}

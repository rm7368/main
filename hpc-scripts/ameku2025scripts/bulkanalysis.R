library(DESeq2)
library(BiocParallel)
library(RUVSeq)
library(tidyverse)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(pheatmap)
library(data.table)
library(jsonlite)
library(GEOquery)
library(ggpubr)
library(introdataviz)

register(SerialParam())

COUNT_DIR <- "~/rnaseq_analysis/data/raw_counts"
METADATA_FILE <- "~/rnaseq_analysis/data/metadata/gse.csv"
GENE_SETS_DIR <- "~/rnaseq_analysis/data/gene_sets"
OUTPUT_DIR <- "~/rnaseq_analysis/results"

MIN_COUNT <- 10
MIN_SAMPLES_FRACTION <- 0.1
PADJ_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 1.0
RUV_K <- 3

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE)

# Load metadata
metadata <- read.csv(METADATA_FILE, row.names = 1)

# Fix the jejunum misspelling
metadata$tissue.ch1 <- gsub("Jeunum", "Jejunum", metadata$tissue.ch1)
print("Fixed jejunum spelling in metadata")

# Check your metadata
print("Metadata structure:")
str(metadata)

# Remove outlier

metadata <- metadata[-14,]
View(metadata)

print("Sample distribution:")
table(metadata$mating.status.ch1, metadata$Sex.ch1, metadata$tissue.ch1)

# Set factor levels (important for comparisons)
metadata$mating.status <- factor(metadata$mating.status.ch1, 
                                 levels = c("Virgin", "Lactating"))
metadata$sex <- factor(metadata$Sex.ch1,
                       levels = c("Female", "Male"))
metadata$region <- factor(metadata$tissue.ch1, 
                          levels = c("Duodenum", "Jejunum", "Ileum"))

# Create combined condition factor for easier comparisons
# Note: Males can only be virgin, so we have 3 conditions total
metadata$condition <- factor(
  paste(metadata$sex, metadata$mating.status, sep = "_"),
  levels = c("Female_Virgin", "Male_Virgin", "Female_Lactating")
)

print("Combined conditions:")
table(metadata$condition, metadata$region)

# Load count data from TSV files
tsv_files <- list.files(COUNT_DIR, pattern = "\\.tsv$", full.names = TRUE)
print(sprintf("Found %d TSV files", length(tsv_files)))

# Read the first file to check format
first_file <- fread(tsv_files[1])
print("First file structure:")
head(first_file)

# Load all count files
# Assuming first column is gene ID and second is counts
count_list <- list()
for (i in seq_along(tsv_files)) {
  file <- tsv_files[i]
  sample_name <- gsub("\\.tsv$", "", basename(file))
  
  data <- fread(file, data.table = FALSE)
  count_vector <- data[[2]]  # Second column has counts
  names(count_vector) <- data[[1]]  # First column has gene names
  
  count_list[[sample_name]] <- count_vector
  
  if (i %% 10 == 0) {
    cat(sprintf("Loaded %d/%d files\n", i, length(tsv_files)))
  }
}

# Combine into matrix
count_matrix <- do.call(cbind, count_list)
print(sprintf("Count matrix dimensions: %d genes x %d samples", 
              nrow(count_matrix), ncol(count_matrix)))

# Check what sample names we have
print("Sample names in count matrix:")
print(colnames(count_matrix))

print("Sample names in metadata:")
print(metadata$title)

# Check if they match
matching <- metadata$title %in% colnames(count_matrix)
if (!all(matching)) {
  print("WARNING: Some metadata samples not found in count matrix:")
  print(metadata$title[!matching])
}

# Also check the reverse
in_metadata <- colnames(count_matrix) %in% metadata$title
if (!all(in_metadata)) {
  print("WARNING: Some count matrix samples not found in metadata:")
  print(colnames(count_matrix)[!in_metadata])
}

# Try to match them - this assumes the order might be the same
if (ncol(count_matrix) == nrow(metadata)) {
  print("Same number of samples. Checking if order matches...")
  
  # Option 1: If samples are in the same order but names don't match
  # colnames(count_matrix) <- metadata$title
  
  # Option 2: Try to match based on partial names
  # This is useful if one has extensions like .tsv and the other doesn't
}

# Ensure column order matches metadata
# Only do this if we can match the samples
if (all(metadata$title %in% colnames(count_matrix))) {
  count_matrix <- count_matrix[, metadata$title]
} else {
  stop("Cannot match samples between count matrix and metadata. Please check sample names.")
}

# Quick look at the data
print("Count matrix preview:")
count_matrix[1:5, 1:5]

# ==============================================================================
# Create DESeq2 Object, initial QC
# ==============================================================================

# First, let's check our data before creating DESeq2 object
print("Checking data before DESeq2 creation:")
print(sprintf("Count matrix dimensions: %d x %d", nrow(count_matrix), ncol(count_matrix)))
print(sprintf("Metadata dimensions: %d x %d", nrow(metadata), ncol(metadata)))

# Check for NULL or duplicate row/column names
print("Checking for NULL names:")
print(sprintf("NULL rownames in count_matrix: %s", is.null(rownames(count_matrix))))
print(sprintf("NULL colnames in count_matrix: %s", is.null(colnames(count_matrix))))

# Check for duplicates
if (!is.null(rownames(count_matrix))) {
  dup_rows <- duplicated(rownames(count_matrix))
  if (any(dup_rows)) {
    print("WARNING: Duplicate row names found:")
    print(head(rownames(count_matrix)[dup_rows]))
    # Fix by making unique
    rownames(count_matrix) <- make.unique(rownames(count_matrix))
  }
}

if (!is.null(colnames(count_matrix))) {
  dup_cols <- duplicated(colnames(count_matrix))
  if (any(dup_cols)) {
    print("WARNING: Duplicate column names found:")
    print(colnames(count_matrix)[dup_cols])
    # Fix by making unique
    colnames(count_matrix) <- make.unique(colnames(count_matrix))
  }
}

# Check that metadata rownames match count matrix colnames
print("Checking sample name alignment:")
print(sprintf("Metadata rownames match count matrix colnames: %s", 
              identical(rownames(metadata), colnames(count_matrix))))

if (!identical(rownames(metadata), colnames(count_matrix))) {
  print("Sample name mismatch!")
  print("First few metadata rownames:")
  print(head(rownames(metadata)))
  print("First few count matrix colnames:")
  print(head(colnames(count_matrix)))
  
  # Try to fix by setting metadata rownames to match
  if (length(rownames(metadata)) == ncol(count_matrix)) {
    print("Attempting to align by position...")
    rownames(metadata) <- colnames(count_matrix)
  }
}

# Ensure count matrix is numeric
if (!is.numeric(count_matrix)) {
  print("Converting count matrix to numeric...")
  count_matrix <- as.matrix(count_matrix)
  mode(count_matrix) <- "numeric"
}

# Remove any rows with all zeros before creating DESeq2 object
zero_rows <- rowSums(count_matrix) == 0
if (any(zero_rows)) {
  print(sprintf("Removing %d rows with all zeros", sum(zero_rows)))
  count_matrix <- count_matrix[!zero_rows, ]
}

# Create DESeq2 dataset
print("Creating DESeq2 dataset...")
dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),  # Ensure integers
  colData = metadata,
  design = ~ region + condition + region:condition
)

print(sprintf("Initial dataset: %d genes x %d samples", nrow(dds), ncol(dds)))

# Pre-filtering to remove low count genes
min_samples <- ceiling(ncol(dds) * MIN_SAMPLES_FRACTION)
keep <- rowSums(counts(dds) >= MIN_COUNT) >= min_samples
dds <- dds[keep,]

print(sprintf("After filtering: %d genes (removed %d low-count genes)", 
              nrow(dds), sum(!keep)))

# Calculate some basic statistics
print("Library sizes (total counts per sample):")
summary(colSums(counts(dds)))

# Variance stabilizing transformation for visualization
vsd <- vst(dds, blind = FALSE)

# PCA plot - see how samples cluster
pca_data <- plotPCA(vsd, intgroup = c("condition", "region"), returnData = TRUE)
percentVar <- attr(pca_data, "percentVar")

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = region)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  scale_color_manual(values = c("Female_Virgin" = "#E31A1C", 
                                "Male_Virgin" = "#1F78B4", 
                                "Female_Lactating" = "#33A02C")) +
  theme_minimal() +
  ggtitle("PCA Before Batch Correction")

print(pca_plot)
ggsave(file.path(OUTPUT_DIR, "plots", "pca_before_correction.pdf"), 
       pca_plot, width = 10, height = 8)

# Sample distance heatmap
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- paste(vsd$condition, vsd$region, sep = "_")
colnames(sample_dist_matrix) <- NULL

pdf(file.path(OUTPUT_DIR, "plots", "sample_distances.pdf"), width = 12, height = 10)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

# Check for batch effects
print("Do we see condition or region effects in PC1 or PC2?")
print("This will help us understand if RUVseq correction is needed")

# Find empirical control genes (least variable genes)
gene_vars <- apply(counts(dds), 1, var)
empirical_controls <- order(gene_vars)[1:1000]  # Use 1000 least variable genes

print(sprintf("Using %d empirical control genes for RUVseq", length(empirical_controls)))

# Apply RUVg correction
counts_matrix <- counts(dds)
ruv_result <- RUVg(counts_matrix, 
                   cIdx = empirical_controls,
                   k = RUV_K)

# Look at the RUV factors
print("RUV factors (first few samples):")
head(ruv_result$W)

# Create new DESeq2 object with RUV factors
ruv_metadata <- cbind(colData(dds), ruv_result$W)
dds_ruv <- DESeqDataSetFromMatrix(
  countData = ruv_result$normalizedCounts,
  colData = ruv_metadata,
  design = ~ W_1 + W_2 + W_3 + region + condition + region:condition
)

# Check PCA after batch correction
vsd_ruv <- vst(dds_ruv, blind = FALSE)
pca_data_ruv <- plotPCA(vsd_ruv, intgroup = c("condition", "region"), returnData = TRUE)

pca_plot_ruv <- ggplot(pca_data_ruv, aes(x = PC1, y = PC2, color = condition, shape = region)) +
  geom_point(size = 4) +
  theme_minimal() +
  scale_color_manual(values = c("Female_Virgin" = "#E31A1C", 
                                "Male_Virgin" = "#1F78B4", 
                                "Female_Lactating" = "#33A02C")) +
  ggtitle("PCA After RUVseq Batch Correction")

print(pca_plot_ruv)
ggsave(file.path(OUTPUT_DIR, "plots", "pca_after_correction.pdf"), 
       pca_plot_ruv, width = 10, height = 8)

#######################

# Run DESeq2 analysis with glmGamPoi for speed
# This is the main statistical analysis step
print("Running DESeq2 analysis (this may take 15-30 minutes)...")
system.time({
  dds_ruv <- DESeq(dds_ruv, 
                   fitType = "glmGamPoi",  # Much faster than standard
                   test = "LRT",
                   reduced = ~ W_1 + W_2 + W_3 + region)
})

# Check dispersion estimates
plotDispEsts(dds_ruv)

# Extract results for specific comparisons
mod_mat <- model.matrix(design(dds_ruv), colData(dds_ruv))
lfi <- colMeans(mod_mat[dds_ruv$condition == "Female_Lactating" & dds_ruv$region == "Ileum",])
vfi <- colMeans(mod_mat[dds_ruv$condition == "Female_Virgin" & dds_ruv$region == "Ileum",])
lfj <- colMeans(mod_mat[dds_ruv$condition == "Female_Lactating" & dds_ruv$region == "Jejunum",])
vfj <- colMeans(mod_mat[dds_ruv$condition == "Female_Virgin" & dds_ruv$region == "Jejunum",])
lfd <- colMeans(mod_mat[dds_ruv$condition == "Female_Lactating" & dds_ruv$region == "Duodenum",])
vfd <- colMeans(mod_mat[dds_ruv$condition == "Female_Virgin" & dds_ruv$region == "Duodenum",])

# Comparisons: Covariate Female Lactating Regions vs Female Virgin Regions
res_lfi_vfi <- results(dds_ruv, contrast = (lfi - vfi),
                       alpha = PADJ_CUTOFF)
res_lfj_vfj <- results(dds_ruv, contrast = (lfj - vfj),
                       alpha = PADJ_CUTOFF)
res_lfd_vfd <- results(dds_ruv, contrast = (lfd - vfd),
                       alpha = PADJ_CUTOFF)


# Comparison 1: Female Lactating vs Female Virgin
res_jejunum_duodenum <- results(dds_ruv, 
                                contrast = c("region", "Jejunum", "Duodenum"),
                                alpha = PADJ_CUTOFF)

print("Female Lactating vs Female Virgin results:")
summary(res_jejunum_duodenum)

# Comparison 2: Male Virgin vs Female Virgin  
res_ileum_duodenum <- results(dds_ruv,
                           contrast = c("region", "Ileum", "Duodenum"),
                           alpha = PADJ_CUTOFF)

print("Male Virgin vs Female Virgin results:")
summary(res_ileum_duodenum)

res_ileum_jejunum <- results(dds_ruv,
                              contrast = c("region", "Ileum", "Jejunum"),
                              alpha = PADJ_CUTOFF)

print("Male Virgin vs Female Virgin results:")
summary(res_ileum_jejunum)

# Save results tables
write.csv(as.data.frame(res_jejunum_duodenum), 
          file.path(OUTPUT_DIR, "tables", "jejunum_vs_duodenum.csv"))
write.csv(as.data.frame(res_ileum_jejunum), 
          file.path(OUTPUT_DIR, "tables", "ileum_vs_jejunum.csv"))
write.csv(as.data.frame(res_ileum_duodenum), 
          file.path(OUTPUT_DIR, "tables", "ileum_vs"))

summary(res_lfi_vfi)
summary(res_lfd_vfd)
summary(res_lfj_vfj)

write.csv(as.data.frame(res_lfi_vfi), 
          file.path(OUTPUT_DIR, "tables", "LFileum_vs_VFileum.csv"))
write.csv(as.data.frame(res_lfj_vfj), 
          file.path(OUTPUT_DIR, "tables", "LFjejunum_vs_VFjejunum.csv"))
write.csv(as.data.frame(res_lfd_vfd), 
          file.path(OUTPUT_DIR, "tables", "LFduodenum_vs_VFduodenum.csv"))

################################

# Female only

dds_ruv_female <- dds_ruv[,dds_ruv$condition %in% c("Female_Virgin", "Female_Lactating")]
vsd_ruv_female <- vst(dds_ruv_female, blind=TRUE)
select <- order(rowVars(assay(vsd_ruv_female)), decreasing=TRUE)[1:50]
mat <- assay(vsd_ruv_female)[select,]
mat <- t(scale(t(mat)))
annotation_df <- data.frame(
  Condition = colData(dds_ruv_female)$condition,
  Region = colData(dds_ruv_female)$region,
  row.names = colnames(mat)
)
pheatmap(mat,
         annotation_col = annotation_df,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D2",
         fontsize_col = 8)

# ==============================================================================
# STEP 6: VISUALIZATION OF RESULTS
# ==============================================================================

# Volcano plot for female lactating vs female virgin
volcano1 <- EnhancedVolcano(res_lfd_vfd,
                            lab = rownames(res_lfd_vfd),
                            x = 'log2FoldChange',
                            y = 'padj',
                            title = 'Female Lactating Duodenum vs Female Virgin Duodenum',
                            pCutoff = PADJ_CUTOFF,
                            FCcutoff = LOG2FC_CUTOFF,
                            pointSize = 1.0,
                            labSize = 3.0,
                            colAlpha = 0.7,
                            drawConnectors = TRUE,
                            widthConnectors = 0.5,
                            max.overlaps = 50)

print(volcano1)
ggsave(file.path(OUTPUT_DIR, "plots", "volcano_LFDuodenum_vs_VFDuodenum.pdf"), 
       volcano1, width = 12, height = 10)

# MA plot
plotMA(res_lactating_virgin, ylim = c(-50, 30000), xlim = c(0, 15000))

# Heatmap of top 50 variable genes
top_var_genes <- head(order(rowVars(assay(vsd_ruv)), decreasing = TRUE), 50)
mat <- assay(vsd_ruv)[top_var_genes,]
mat <- t(scale(t(mat)))  # Z-score normalization

# Create annotation for heatmap
annotation_col <- data.frame(
  Condition = colData(dds_ruv)$condition,
  Region = colData(dds_ruv)$region,
  row.names = colnames(mat)
)

ann_colors <- list(
  Condition = c(Female_Virgin = "#E31A1C", 
                Male_Virgin = "#1F78B4", 
                Female_Lactating = "#33A02C"),
  Region = c(Duodenum = "#FF7F00", 
             Jejunum = "#6A3D9A", 
             Ileum = "#CAB2D6")
)

pdf(file.path(OUTPUT_DIR, "plots", "top50_genes_heatmap.pdf"), width = 14, height = 10)
pheatmap(mat,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D2",
         fontsize_col = 8)
dev.off()

#####################

# Get significant genes from female lactating vs female virgin comparison
sig_genes_up <- rownames(res_lactating_virgin)[
  !is.na(res_lactating_virgin$padj) & 
    res_lactating_virgin$padj < PADJ_CUTOFF & 
    res_lactating_virgin$log2FoldChange > LOG2FC_CUTOFF
]

sig_genes_down <- rownames(res_lactating_virgin)[
  !is.na(res_lactating_virgin$padj) & 
    res_lactating_virgin$padj < PADJ_CUTOFF & 
    res_lactating_virgin$log2FoldChange < -LOG2FC_CUTOFF
]

print(sprintf("Significant genes (Female Lactating vs Virgin): %d up, %d down", 
              length(sig_genes_up), length(sig_genes_down)))

# Convert gene symbols to Entrez IDs for enrichment analysis
gene_list <- c(sig_genes_up, sig_genes_down)

# This might need adjustment based on your gene ID format
entrez_ids <- bitr(sig_genes_up, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Mm.eg.db)

# GO enrichment analysis
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       ont = "BP",  # Biological Process
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

print("Top enriched GO terms:")
head(go_results)

# Visualize GO results
dotplot(go_results, showCategory = 20)
ggsave(file.path(OUTPUT_DIR, "plots", "GO_dotplot.pdf"), width = 10, height = 8)

# KEGG pathway analysis
kegg_results <- enrichKEGG(gene = entrez_ids$ENTREZID,
                           organism = 'mmu',  # Mouse
                           pvalueCutoff = 0.05)

if (nrow(kegg_results) > 0) {
  dotplot(kegg_results, showCategory = 20)
  ggsave(file.path(OUTPUT_DIR, "plots", "KEGG_dotplot.pdf"), width = 10, height = 8)
}

######################

# Get list of all JSONL files in the gene sets directory
jsonl_files <- list.files(GENE_SETS_DIR, pattern = "\\.jsonl$", full.names = TRUE)
print(sprintf("Found %d gene set JSONL files:", length(jsonl_files)))
print(basename(jsonl_files))

# Load and parse all gene sets
all_gene_sets <- list()

for (jsonl_file in jsonl_files) {
  set_name <- gsub("\\.jsonl$", "", basename(jsonl_file))
  print(sprintf("\nLoading gene set: %s", set_name))
  
  # Read JSONL file line by line
  lines <- readLines(jsonl_file)
  print(sprintf("  Found %d entries", length(lines)))
  
  # Parse each line as JSON
  gene_data <- list()
  for (i in seq_along(lines)) {
    if (nchar(trimws(lines[i])) > 0) {  # Skip empty lines
      gene_data[[i]] <- fromJSON(lines[i])
    }
  }
  
  # Remove any NULL entries
  gene_data <- gene_data[!sapply(gene_data, is.null)]
  
  # Show structure of first entry to understand the format
  if (length(gene_data) > 0) {
    print("  Structure of first gene entry:")
    str(gene_data[[1]], max.level = 2)
  }
  
  # Extract gene symbols - NCBI gene JSONL typically has these structures
  gene_symbols <- character()
  gene_descriptions <- character()
  
  for (entry in gene_data) {
    # Extract symbol - try multiple possible locations
    symbol <- NULL
    if (!is.null(entry$symbol)) {
      symbol <- entry$symbol
    } else if (!is.null(entry$gene$symbol)) {
      symbol <- entry$gene$symbol
    } else if (!is.null(entry$gene_symbol)) {
      symbol <- entry$gene_symbol
    }
    
    # Extract description
    description <- NULL
    if (!is.null(entry$description)) {
      description <- entry$description
    } else if (!is.null(entry$gene$description)) {
      description <- entry$gene$description
    } else if (!is.null(entry$name)) {
      description <- entry$name
    } else {
      description <- ""
    }
    
    if (!is.null(symbol)) {
      gene_symbols <- c(gene_symbols, symbol)
      gene_descriptions <- c(gene_descriptions, description)
    }
  }
  
  print(sprintf("  Successfully extracted %d gene symbols", length(gene_symbols)))
  
  if (length(gene_symbols) > 0) {
    print("  First 5 genes:")
    print(head(gene_symbols, 5))
  } else {
    warning(sprintf("No gene symbols found in %s", jsonl_file))
    next
  }
  
  # Store the gene set
  all_gene_sets[[set_name]] <- list(
    symbols = gene_symbols,
    descriptions = gene_descriptions,
    size = length(gene_symbols)
  )
  
  print(sprintf("  Loaded %d genes", length(gene_symbols)))
}

# Summary of all gene sets
gene_set_summary <- data.frame(
  GeneSet = names(all_gene_sets),
  NumGenes = sapply(all_gene_sets, function(x) x$size),
  row.names = NULL
)
print("\nGene set summary:")
print(gene_set_summary)

# Check overlap between gene sets
if (length(all_gene_sets) > 1) {
  print("\nOverlap between gene sets:")
  for (i in 1:(length(all_gene_sets)-1)) {
    for (j in (i+1):length(all_gene_sets)) {
      set1_name <- names(all_gene_sets)[i]
      set2_name <- names(all_gene_sets)[j]
      overlap <- length(intersect(all_gene_sets[[i]]$symbols, 
                                  all_gene_sets[[j]]$symbols))
      if (overlap > 0) {
        print(sprintf("  %s ∩ %s: %d genes", set1_name, set2_name, overlap))
      }
    }
  }
}

# Analyze each gene set against your DE results
gene_set_results <- list()

for (set_name in names(all_gene_sets)) {
  print(sprintf("\nAnalyzing gene set: %s", set_name))
  
  gene_set <- all_gene_sets[[set_name]]$symbols
  
  # Check overlap with significant genes
  overlap_up <- intersect(sig_genes_up, gene_set)
  overlap_down <- intersect(sig_genes_down, gene_set)
  overlap_all <- c(overlap_up, overlap_down)
  
  print(sprintf("  Overlap: %d up-regulated, %d down-regulated (total: %d/%d)", 
                length(overlap_up), length(overlap_down), 
                length(overlap_all), length(gene_set)))
  
  # Store results
  gene_set_results[[set_name]] <- list(
    total_genes = length(gene_set),
    overlap_up = overlap_up,
    overlap_down = overlap_down,
    overlap_all = overlap_all,
    percent_overlap = 100 * length(overlap_all) / length(gene_set)
  )
  
  # Create visualization if there's meaningful overlap
  if (length(overlap_all) >= 5) {  # At least 5 genes
    # Extract expression values for overlapping genes
    overlap_expr <- assay(vsd_ruv)[overlap_all, ]
    
    # Create heatmap for this gene set
    pdf(file.path(OUTPUT_DIR, "plots", sprintf("%s_heatmap.pdf", set_name)), 
        width = 12, height = max(8, length(overlap_all) * 0.3))
    
    # Add row annotation showing up/down regulation
    row_annotation <- data.frame(
      Direction = ifelse(overlap_all %in% overlap_up, "Up", "Down"),
      row.names = overlap_all
    )
    
    row_colors <- list(Direction = c(Up = "red", Down = "blue"))
    
    pheatmap(overlap_expr,
             annotation_col = annotation_col,
             annotation_row = row_annotation,
             annotation_colors = c(ann_colors, row_colors),
             scale = "row",
             show_rownames = TRUE,
             show_colnames = FALSE,
             main = sprintf("%s Gene Set Expression", set_name),
             fontsize_row = min(10, 300/length(overlap_all)))  # Adjust font size
    dev.off()
    
    # Create a volcano plot highlighting this gene set
    volcano_data <- as.data.frame(res_lactating_virgin) %>%
      mutate(
        gene = rownames(.),
        in_set = gene %in% gene_set,
        sig_in_set = gene %in% overlap_all
      )
    
    volcano_geneset <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = sig_in_set), alpha = 0.6, size = 2) +
      scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red"),
                         labels = c("Other genes", sprintf("%s genes", set_name))) +
      geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed") +
      geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed") +
      theme_minimal() +
      labs(title = sprintf("Volcano Plot - %s Gene Set (Female Lactating vs Virgin)", set_name),
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-value",
           color = "Gene Set") +
      theme(legend.position = "bottom")
    
    # Add labels for top genes in the set
    top_genes_in_set <- volcano_data %>%
      filter(sig_in_set) %>%
      arrange(padj) %>%
      head(10)
    
    if (nrow(top_genes_in_set) > 0) {
      volcano_geneset <- volcano_geneset +
        ggrepel::geom_text_repel(
          data = top_genes_in_set,
          aes(label = gene),
          size = 3,
          max.overlaps = 20
        )
    }
    
    print(volcano_geneset)
    ggsave(file.path(OUTPUT_DIR, "plots", sprintf("%s_volcano.pdf", set_name)), 
           volcano_geneset, width = 10, height = 8)
  }
}

# Create summary table of all gene set results
gene_set_summary_table <- do.call(rbind, lapply(names(gene_set_results), function(x) {
  data.frame(
    GeneSet = x,
    TotalGenes = gene_set_results[[x]]$total_genes,
    UpRegulated = length(gene_set_results[[x]]$overlap_up),
    DownRegulated = length(gene_set_results[[x]]$overlap_down),
    TotalOverlap = length(gene_set_results[[x]]$overlap_all),
    PercentOverlap = round(gene_set_results[[x]]$percent_overlap, 1)
  )
}))

print("\nGene set enrichment summary:")
print(gene_set_summary_table)

write.csv(gene_set_summary_table, 
          file.path(OUTPUT_DIR, "tables", "gene_set_summary.csv"), 
          row.names = FALSE)

# Create a combined visualization showing all gene sets
if (length(all_gene_sets) > 1) {
  # Prepare data for upset plot (shows overlaps between sets and DE genes)
  library(UpSetR)
  
  # Create a matrix of gene membership
  all_genes_in_sets <- unique(unlist(lapply(all_gene_sets, function(x) x$symbols)))
  
  # Create binary matrix
  upset_matrix <- matrix(0, 
                         nrow = length(all_genes_in_sets), 
                         ncol = length(all_gene_sets) + 2)
  rownames(upset_matrix) <- all_genes_in_sets
  colnames(upset_matrix) <- c(names(all_gene_sets), "Sig_Up", "Sig_Down")
  
  # Fill in the matrix
  for (i in seq_along(all_gene_sets)) {
    upset_matrix[all_genes_in_sets %in% all_gene_sets[[i]]$symbols, i] <- 1
  }
  upset_matrix[all_genes_in_sets %in% sig_genes_up, "Sig_Up"] <- 1
  upset_matrix[all_genes_in_sets %in% sig_genes_down, "Sig_Down"] <- 1
  
  # Create upset plot
  pdf(file.path(OUTPUT_DIR, "plots", "gene_set_overlap_upset.pdf"), 
      width = 12, height = 8)
  upset(as.data.frame(upset_matrix), 
        order.by = "freq",
        nsets = ncol(upset_matrix),
        nintersects = 20,
        main.bar.color = "black",
        sets.bar.color = c(rep("darkblue", length(all_gene_sets)), "red", "blue"))
  dev.off()
}

# Perform enrichment analysis using your custom gene sets
# Convert to format for clusterProfiler
custom_term2gene <- do.call(rbind, lapply(names(all_gene_sets), function(set_name) {
  data.frame(
    term = set_name,
    gene = all_gene_sets[[set_name]]$symbols,
    stringsAsFactors = FALSE
  )
}))

custom_term2name <- data.frame(
  term = names(all_gene_sets),
  name = names(all_gene_sets),  # Use set names as descriptions
  stringsAsFactors = FALSE
)

# Run enrichment analysis with your custom gene sets
custom_enrichment <- enricher(
  gene = c(sig_genes_up, sig_genes_down),
  TERM2GENE = custom_term2gene,
  TERM2NAME = custom_term2name,
  pvalueCutoff = 1,  # Show all results
  qvalueCutoff = 1
)

if (nrow(custom_enrichment) > 0) {
  print("\nCustom gene set enrichment results:")
  print(as.data.frame(custom_enrichment))
  
  # Visualize enrichment
  dotplot(custom_enrichment, showCategory = length(all_gene_sets)) +
    ggtitle("Custom Gene Set Enrichment")
  ggsave(file.path(OUTPUT_DIR, "plots", "custom_geneset_enrichment.pdf"), 
         width = 10, height = 6)
}

# Collect all comparisons
all_comparisons <- list(
  "Female_Lactating_vs_Female_Virgin" = res_lactating_virgin,
  "Male_Virgin_vs_Female_Virgin" = res_male_female
)

# Additional comparisons you might want:
# Female Lactating vs Male Virgin (if biologically relevant)
res_lactating_male <- results(dds_ruv,
                              contrast = c("condition", "Female_Lactating", "Male_Virgin"),
                              alpha = PADJ_CUTOFF)
all_comparisons[["Female_Lactating_vs_Male_Virgin"]] <- res_lactating_male

# Region comparisons
res_duodenum_jejunum <- results(dds_ruv,
                                contrast = c("region", "Duodenum", "Jejunum"),
                                alpha = PADJ_CUTOFF)
res_jejunum_ileum <- results(dds_ruv,
                             contrast = c("region", "Jejunum", "Ileum"),
                             alpha = PADJ_CUTOFF)
all_comparisons[["Duodenum_vs_Jejunum"]] <- res_duodenum_jejunum
all_comparisons[["Jejunum_vs_Ileum"]] <- res_jejunum_ileum

# Summary statistics
summary_stats <- data.frame(
  Comparison = names(all_comparisons),
  Total_Genes = sapply(all_comparisons, nrow),
  Significant_Genes = sapply(all_comparisons, function(x) 
    sum(!is.na(x$padj) & x$padj < PADJ_CUTOFF)),
  Up_Regulated = sapply(all_comparisons, function(x) 
    sum(!is.na(x$padj) & x$padj < PADJ_CUTOFF & x$log2FoldChange > LOG2FC_CUTOFF)),
  Down_Regulated = sapply(all_comparisons, function(x) 
    sum(!is.na(x$padj) & x$padj < PADJ_CUTOFF & x$log2FoldChange < -LOG2FC_CUTOFF))
)

print("Summary of all comparisons:")
print(summary_stats)

write.csv(summary_stats, file.path(OUTPUT_DIR, "summary_statistics.csv"), row.names = FALSE)

# Save the final DESeq2 object for later use
saveRDS(dds_ruv, file.path(OUTPUT_DIR, "final_deseq2_object.rds"))

# ==============================================================================
# STEP 10: ADVANCED VISUALIZATIONS (OPTIONAL)
# ==============================================================================

# Gene-concept network plot for GO results
if (nrow(go_results) > 0) {
  # Add fold change information
  foldchanges <- res_lactating_virgin$log2FoldChange
  names(foldchanges) <- rownames(res_lactating_virgin)
  
  # Calculate pairwise term similarity
  go_sim <- pairwise_termsim(go_results)
  
  # Create concept network
  cnet <- cnetplot(go_sim, 
                   showCategory = 10,
                   foldChange = foldchanges,
                   cex_gene = 0.8,
                   cex_category = 1.2)
  
  print(cnet)
  ggsave(file.path(OUTPUT_DIR, "plots", "GO_concept_network.pdf"), 
         cnet, width = 14, height = 12)
}

# Interactive volcano plot using plotly
library(plotly)

volcano_data <- as.data.frame(res_lactating_virgin) %>%
  mutate(
    gene = rownames(.),
    sig = ifelse(!is.na(padj) & padj < PADJ_CUTOFF & abs(log2FoldChange) > LOG2FC_CUTOFF,
                 "Significant", "Not Significant"),
    hover_text = paste("Gene:", gene, 
                       "<br>Log2FC:", round(log2FoldChange, 3),
                       "<br>Adj P-value:", format(padj, scientific = TRUE, digits = 3))
  )

interactive_volcano <- plot_ly(
  data = volcano_data,
  x = ~log2FoldChange,
  y = ~-log10(padj),
  color = ~sig,
  colors = c("Not Significant" = "gray", "Significant" = "red"),
  type = "scatter",
  mode = "markers",
  hovertext = ~hover_text,
  hoverinfo = "text"
) %>%
  layout(
    title = "Interactive Volcano Plot",
    xaxis = list(title = "Log2 Fold Change"),
    yaxis = list(title = "-Log10 Adjusted P-value")
  )

# Save interactive plot
htmlwidgets::saveWidget(interactive_volcano, 
                        file.path(OUTPUT_DIR, "plots", "interactive_volcano.html"))

print("Analysis complete! Check the results directory for all outputs.")
print(sprintf("Results saved in: %s", OUTPUT_DIR))

# ==============================================================================
# HELPFUL COMMANDS FOR EXPLORING YOUR RESULTS
# ==============================================================================

# To explore specific genes:
plotCounts(dds_ruv, gene = "Piezo1", intgroup = c("condition", "region"), returnData = TRUE)

# To get normalized counts for specific genes:
# normalized_counts <- counts(dds_ruv, normalized = TRUE)
# normalized_counts["YourGene", ]

# To export normalized counts for all genes:
write.csv(counts(dds_ruv, normalized = TRUE), 
          file.path(OUTPUT_DIR, "tables", "normalized_counts.csv"))

# To look at specific contrasts between regions:
res_duodenum_jejunum <- results(dds_ruv, 
                                 contrast = c("region", "duodenum", "jejunum"))

# Memory usage check:
gc()
object.size(dds_ruv) / 1024^2  # Size in GB


#########violin###############

print("\n=== Creating violin plots for specific genes ===")

# Define your genes of interest here
# You can modify this list with your specific genes
genes_of_interest <- c("Pgr") # REPLACE WITH YOUR GENES

# Or if you want to use genes from your NCBI gene sets:
# genes_of_interest <- unique(unlist(lapply(all_gene_sets, function(x) x$symbols)))[1:10]

# Get normalized counts
normalized_counts <- counts(dds_ruv, normalized = TRUE)

# Check which genes are actually in the dataset
genes_found <- genes_of_interest[genes_of_interest %in% rownames(normalized_counts)]
genes_not_found <- setdiff(genes_of_interest, genes_found)

# Extract normalized counts for genes of interest
gene_counts <- as.data.frame(t(normalized_counts[genes_found, , drop = FALSE]))
gene_counts$sample <- rownames(gene_counts)

# Add metadata
gene_counts <- merge(gene_counts, 
                     data.frame(sample = rownames(colData(dds_ruv)),
                                condition = colData(dds_ruv)$condition,
                                region = colData(dds_ruv)$region,
                                sex = colData(dds_ruv)$sex),
                     by = "sample")

# Convert to long format for ggplot
gene_counts_long <- pivot_longer(gene_counts, 
                                 cols = all_of(genes_found),
                                 names_to = "gene",
                                 values_to = "normalized_count")

# Add log2 transformed counts (adding pseudocount of 1)
gene_counts_long$log2_count <- log2(gene_counts_long$normalized_count + 1)
gene_counts_long <- gene_counts_long %>%
  filter(condition == "Female_Lactating"| condition == "Female_Virgin")

# Define comparisons for condition (adjust based on your specific comparisons of interest)
condition_comparisons <- list(
  c("Female_Virgin", "Female_Lactating")
)

# Create violin plot with statistical significance for conditions
violin_all <- ggplot(gene_counts_long, aes(x = condition, y = log2_count, fill = condition)) +
  geom_violin(alpha = 0.8, scale = "width") +
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.05, alpha = 0.6, size = 2) +
  stat_compare_means(comparisons = condition_comparisons,
                     method = "t.test",  # or "wilcox.test" for non-parametric
                     label = "p.signif",
                     step.increase = 0.1,
                     bracket.size = 0.6,
                     size = 3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Female_Virgin" = "#E31A1C",
                               "Female_Lactating" = "#33A02C")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(title = "Gene Expression Across Conditions",
       x = "Condition",
       y = "Log2(Normalized Count + 1)",
       fill = "Condition")

print(violin_all)
# ggsave(file.path(OUTPUT_DIR, "plots", "violin_genes_by_condition.pdf"), 
#       violin_all, width = 12, height = 3 * ceiling(length(genes_found)/2))


# Define comparisons for region
region_comparisons <- list(
  c("Duodenum", "Jejunum"),
  c("Duodenum", "Ileum"),
  c("Jejunum", "Ileum")
)

# Create violin plot by region with statistical significance
violin_region <- ggplot(gene_counts_long, aes(x = region, y = log2_count, fill = region)) +
  geom_violin(alpha = 0.8, scale = "width") +
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.05, alpha = 0.6, size = 2) +
  stat_compare_means(comparisons = region_comparisons,
                     method = "t.test",
                     label = "p.signif",
                     step.increase = 0.1,
                     bracket.size = 0.6,
                     size = 3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Duodenum" = "#FF7F00", 
                               "Jejunum" = "#6A3D9A", 
                               "Ileum" = "#CAB2D6")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(title = "Gene Expression Across Intestinal Regions",
       x = "Region",
       y = "Log2(Normalized Count + 1)",
       fill = "Region")

print(violin_region)
# ggsave(file.path(OUTPUT_DIR, "plots", "violin_genes_by_region.pdf"), 
#       violin_region, width = 12, height = 3 * ceiling(length(genes_found)/2))

# Split violin plots
violin_region <- ggplot(gene_counts_long, aes(x = region, y = log2_count, fill = condition)) +
  introdataviz::geom_split_violin(alpha = 0.8, scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.5, fatten = NULL, outlier.shape = NA) +
  geom_jitter(width = 0.05, alpha = 0.6, size = 2) +
  stat_compare_means(comparisons = region_comparisons,
                     method = "t.test",
                     label = "p.signif",
                     step.increase = 0.1,
                     bracket.size = 0.6,
                     size = 3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Female_Lactating" = "#FF7F00", 
                               "Female_Virgin" = "#6A3D9A")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(title = "Gene Expression Across Intestinal Regions",
       x = "Region",
       y = "Log2(Normalized Count + 1)",
       fill = "Condition")
print(violin_region)

print("\n=== Analysis complete! ===")
print(sprintf("All results saved to: %s", OUTPUT_DIR))

############raincloud plots#############################
rain_height <- .1

rain <- ggplot(gene_counts_long, aes(x = region, y = log2_count, fill = condition)) +
    # clouds
    introdataviz::geom_flat_violin(trim=FALSE, alpha = 0.4,
                                   position = position_nudge(x = rain_height+.05)) +
    # rain
    geom_point(aes(colour = condition), size = 2, alpha = .5, show.legend = FALSE, 
               position = position_jitter(width = rain_height, height = 0)) +
    # boxplots
    geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE, 
                 outlier.shape = NA,
                 position = position_nudge(x = -rain_height*2)) +
    # mean and SE point in the cloud
    stat_summary(fun.data = mean_cl_normal, mapping = aes(color = condition), show.legend = FALSE,
                 position = position_nudge(x = rain_height * 3)) +
    stat_compare_means(comparisons = region_comparisons,
                       method = "t.test",
                       label = "p.signif",
                       step.increase = 0.1,
                       bracket.size = 0.6,
                       size = 3) +
    coord_flip() +
    facet_wrap(~ gene, scales = "free_y", nrow = 2) +
    # custom colours and theme
    scale_fill_brewer(palette = "Dark2", name = "Condition") +
    scale_colour_brewer(palette = "Dark2") +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          legend.position = c(0.2, 0.1),
          legend.background = element_rect(fill = "white", color = "white"))
print(rain)

########### intestinal region volcano plots #############

volcano1 <- EnhancedVolcano(res_jejunum_duodenum,
                            lab = rownames(res_jejunum_duodenum),
                            x = 'log2FoldChange',
                            y = 'padj',
                            title = 'Jejunum vs Duodenum',
                            pCutoff = PADJ_CUTOFF,
                            FCcutoff = LOG2FC_CUTOFF,
                            pointSize = 1.0,
                            labSize = 3.0,
                            colAlpha = 0.7,
                            drawConnectors = TRUE,
                            widthConnectors = 0.5,
                            max.overlaps = 80)

print(volcano1)
ggsave(file.path(OUTPUT_DIR, "plots", "volcano_jejunum_vs_duodenum.pdf"), 
       volcano1, width = 12, height = 10)


volcano2 <- EnhancedVolcano(res_ileum_duodenum,
                            lab = rownames(res_ileum_duodenum),
                            x = 'log2FoldChange',
                            y = 'padj',
                            title = 'Ileum vs Duodenum',
                            pCutoff = PADJ_CUTOFF,
                            FCcutoff = LOG2FC_CUTOFF,
                            pointSize = 1.0,
                            labSize = 3.0,
                            colAlpha = 0.7,
                            drawConnectors = TRUE,
                            widthConnectors = 0.5,
                            max.overlaps = 80)

print(volcano2)
ggsave(file.path(OUTPUT_DIR, "plots", "volcano_ileum_vs_duodenum.pdf"), 
       volcano2, width = 12, height = 10)

volcano3 <- EnhancedVolcano(res_ileum_jejunum,
                            lab = rownames(res_ileum_jejunum),
                            x = 'log2FoldChange',
                            y = 'padj',
                            title = 'Ileum vs Jejunum',
                            pCutoff = PADJ_CUTOFF,
                            FCcutoff = LOG2FC_CUTOFF,
                            pointSize = 1.0,
                            labSize = 3.0,
                            colAlpha = 0.7,
                            drawConnectors = TRUE,
                            widthConnectors = 0.5,
                            max.overlaps = 80)

print(volcano3)
ggsave(file.path(OUTPUT_DIR, "plots", "volcano_ileum_vs_jejunum.pdf"), 
       volcano3, width = 12, height = 10)

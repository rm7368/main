# ================================================================================
# Comprehensive Single-Cell RNA Sequencing Analysis Pipeline
# ================================================================================
# This pipeline processes 10x Genomics single-cell RNA-seq data following
# the methods described in your study. Optimized for high-performance computing.
# 
# Study design: 17 samples (3 lactating KO, 11 lactating WT, 3 virgin WT mice)
# ================================================================================

# Load required libraries
# ================================================================================
library(Seurat)          # v4.3.0+ for single-cell analysis
library(SeuratObject)    # Seurat object handling
library(dplyr)           # Data manipulation
library(ggplot2)         # Plotting
library(patchwork)       # Combining plots
library(Matrix)          # Sparse matrix operations
library(future)          # Parallel processing
library(DESeq2)          # Pseudobulk differential expression
library(clusterProfiler) # Gene ontology analysis
library(enrichplot)      # GO visualization
library(org.Mm.eg.db)    # Mouse gene annotations
library(RColorBrewer)    # Color palettes
library(viridis)         # Color schemes
library(pheatmap)        # Heatmaps
library(cowplot)         # Plot arrangements
library(stringr)         # String manipulation
library(readr)           # File reading
library(tidyr)           # Data tidying
library(GEOquery)

# Set up parallel processing for HPC
# ================================================================================
options(future.globals.maxSize = Inf)  # 50GB memory limit

# Set random seed for reproducibility
set.seed(42)

# Define sample information
# ================================================================================
# Create sample metadata based on your study design
# Create sample metadata based on your study design with batch effects
# Batch 1: KO Lactating (3) + WT Lactating (6) = 9 samples
# Batch 2: WT Lactating (5) + WT Virgin (3) = 8 samples

sample_info <- data.frame(
  sample_id = c(paste0("KO_Lact_", 1:3),      # Batch 1
                paste0("WT_Lact_", 1:6),      # Batch 1
                paste0("WT_Lact_", 7:11),     # Batch 2
                paste0("WT_Virgin_", 1:3)),   # Batch 2
  condition = c(rep("KO_Lactating", 3),
                rep("WT_Lactating", 11),
                rep("WT_Virgin", 3)),
  genotype = c(rep("KO", 3),
               rep("WT", 14)),
  state = c(rep("Lactating", 14),
            rep("Virgin", 3)),
  batch = c(rep("Batch1", 9),                 # First instrument
            rep("Batch2", 8)),                # Second instrument
  instrument = c(rep("NextSeq 2000", 9),
                 rep("Illumina HiSeq 2500", 8)),
  stringsAsFactors = FALSE
)

# Display sample information
print("Sample Information with Batch Effects:")
print(sample_info)
# Set file paths
# ================================================================================
# Adjust these paths to match your directory structure
base_path <- "~/scseq_analysis"  # Change this to your data directory
metadata_file_batch1 <- file.path(base_path, "GSE247927-GPL30172_series_matrix.txt")  # First batch metadata
metadata_file_batch2 <- file.path(base_path, "GSE247927-GPL17021_series_matrix.txt")  # Second batch metadata

# Function to read series matrix files
# ================================================================================
gse1 <- getGEO(filename="~/scseq_analysis/GSE247927-GPL30172_series_matrix.txt")
gse2 <- getGEO(filename="~/scseq_analysis/GSE247927-GPL17021_series_matrix.txt")

# Read metadata files
batch1_metadata <- (gse1)
batch2_metadata <- (gse2)

# Display metadata structure if successfully loaded
if(!is.null(batch1_metadata)) {
  cat("Batch 1 metadata structure:\n")
  str(batch1_metadata)
}

if(!is.null(batch2_metadata)) {
  cat("Batch 2 metadata structure:\n")
  str(batch2_metadata)
}

# Function to read 10x data with error handling
# ================================================================================
read_10x_safe <- function(data_dir, sample_name) {
  tryCatch({
    # Read 10x data
    data <- Read10X(data.dir = data_dir)
    
    # Create Seurat object with initial filtering
    integrated_seurat <- CreateSeuratObject(
      counts = data,
      project = sample_name,
      min.cells = 3,      # Filter genes present in at least 3 cells
      min.features = 100  # Filter cells with at least 100 features
    )
    
    # Add sample metadata
    integrated_seurat$sample <- sample_name
    integrated_seurat$orig.ident <- sample_name
    
    return(integrated_seurat)
  }, error = function(e) {
    warning(paste("Failed to read sample", sample_name, ":", e$message))
    return(NULL)
  })
}

# Function to read 10x data
# ================================================================================
read_10x_safe <- function(data_dir, sample_name) {
  tryCatch({
    # Read 10x data
    data <- Read10X(data.dir = data_dir)
    
    # Handle case where Read10X returns a list (multiple feature types)
    if(is.list(data)) {
      if("Gene Expression" %in% names(data)) {
        data <- data[["Gene Expression"]]
      } else {
        # Use the first element if "Gene Expression" not found
        data <- data[[1]]
      }
    }
    
    # Create Seurat object with initial filtering
    integrated_seurat <- CreateSeuratObject(
      counts = data,
      project = sample_name,
      min.cells = 3,      # Filter genes present in at least 3 cells
      min.features = 100  # Filter cells with at least 100 features
    )
    
    # Add sample metadata
    integrated_seurat$sample <- sample_name
    integrated_seurat$orig.ident <- sample_name
    
    return(integrated_seurat)
  }, error = function(e) {
    warning(paste("Failed to read sample", sample_name, ":", e$message))
    return(NULL)
  })
}

# Load and create Seurat objects for all samples
# ================================================================================
seurat_list <- list()

for(i in 1:nrow(sample_info)) {
  sample_id <- sample_info$sample_id[i]
  sample_path <- file.path(base_path, sample_id)
  
  cat("Loading sample:", sample_id, "\n")
  
  # Read 10x data
  integrated_seurat <- read_10x_safe(sample_path, sample_id)
  
  if(!is.null(integrated_seurat)) {
    # Add metadata from sample_info
    integrated_seurat$condition <- sample_info$condition[i]
    integrated_seurat$genotype <- sample_info$genotype[i]
    integrated_seurat$state <- sample_info$state[i]
    integrated_seurat$batch <- sample_info$batch[i]
    integrated_seurat$instrument <- sample_info$instrument[i]
    
    seurat_list[[sample_id]] <- integrated_seurat
    cat("Successfully loaded", ncol(integrated_seurat), "cells for", sample_id, "\n")
  }
}

# Remove any NULL objects
seurat_list <- seurat_list[!sapply(seurat_list, is.null)]
cat("Successfully loaded", length(seurat_list), "samples\n")


# Quality Control Metrics Calculation
# ================================================================================
# Calculate QC metrics for each sample
for(i in names(seurat_list)) {
  # Calculate mitochondrial gene percentage
  seurat_list[[i]][["percent.mt"]] <- PercentageFeatureSet(seurat_list[[i]], pattern = "^mt-")
  
  # Calculate total UMI and feature counts
  seurat_list[[i]][["log10GenesPerUMI"]] <- log10(seurat_list[[i]]$nFeature_RNA) / log10(seurat_list[[i]]$nCount_RNA)
}

# Comprehensive Quality Control Visualization
# ================================================================================
# Function to create QC plots
create_qc_plots <- function(integrated_seurat, sample_name) {
  # Violin plots for QC metrics
  p1 <- VlnPlot(integrated_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, pt.size = 0.1) + 
    plot_annotation(title = paste("QC Metrics -", sample_name))
  
  # Scatter plots for QC relationships
  p2 <- FeatureScatter(integrated_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    geom_hline(yintercept = 25, color = "red", linetype = "dashed") +
    ggtitle("UMI vs Mitochondrial %")
  
  p3 <- FeatureScatter(integrated_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    ggtitle("UMI vs Feature Count")
  
  return(list(violin = p1, scatter1 = p2, scatter2 = p3))
}

# Generate QC plots for each sample
qc_plots <- lapply(names(seurat_list), function(x) {
  create_qc_plots(seurat_list[[x]], x)
})
names(qc_plots) <- names(seurat_list)

# Display QC summary statistics
cat("\nQuality Control Summary:\n")
for(sample in names(seurat_list)) {
  obj <- seurat_list[[sample]]
  cat("\nSample:", sample, "\n")
  cat("Cells:", ncol(obj), "\n")
  cat("Genes:", nrow(obj), "\n")
  cat("Median UMI:", median(obj$nCount_RNA), "\n")
  cat("Median Genes:", median(obj$nFeature_RNA), "\n")
  cat("Median MT%:", round(median(obj$percent.mt), 2), "\n")
  cat("Batch:", unique(obj$batch), "\n")
  cat("Instrument:", unique(obj$instrument), "\n")
}


# Quality Control Filtering
# ================================================================================
# Apply filtering criteria based on your methods
filter_cells <- function(integrated_seurat) {
  # Initial cell count
  initial_cells <- ncol(integrated_seurat)
  
  # Apply filters: UMI > 500, Features > 100, MT% < 25%
  integrated_seurat <- subset(integrated_seurat, 
                       subset = nCount_RNA > 500 & 
                         nFeature_RNA > 100 & 
                         nFeature_RNA < 8000 & # Remove potential doublets
                         percent.mt < 25)
  
  # Final cell count
  final_cells <- ncol(integrated_seurat)
  
  cat("Filtered", initial_cells - final_cells, "cells, retained", final_cells, "cells\n")
  
  return(integrated_seurat)
}

# Apply filtering to all samples
seurat_list <- lapply(seurat_list, filter_cells)

# Remove samples with too few cells (< 100)
min_cells <- 100
seurat_list <- seurat_list[sapply(seurat_list, ncol) >= min_cells]

cat("Retained", length(seurat_list), "samples after quality filtering\n")

# Merge all samples into a single Seurat object
# ================================================================================
cat("Merging all samples...\n")
combined_seurat <- merge(seurat_list[[1]], 
                         y = seurat_list[2:length(seurat_list)], 
                         add.cell.ids = names(seurat_list),
                         project = "IntestinalSC")

cat("Combined object contains", ncol(combined_seurat), "cells and", nrow(combined_seurat), "genes\n")

# SCTransform Normalization and Integration
# ================================================================================
cat("Performing SCTransform normalization...\n")

# Split object by sample for integration
sample_list <- SplitObject(combined_seurat, split.by = "orig.ident")

# SCTransform normalization for each sample with batch-specific regression
sample_list <- lapply(sample_list, function(x) {
  # Get batch information for this sample
  sample_batch <- unique(x$batch)
  
  # Determine which variables to regress based on variation within sample
  vars_to_regress <- c("percent.mt", "nCount_RNA")
  
  # Only include batch if there's variation within this sample (shouldn't be the case, but safety check)
  if(length(unique(x$batch)) > 1) {
    vars_to_regress <- c(vars_to_regress, "batch")
  }
  
  # Perform SCTransform with enhanced parameters including batch effects
  x <- SCTransform(x, 
                   method = "glmGamPoi",  # Faster method
                   vars.to.regress = vars_to_regress,  # Only regress variables with variation
                   verbose = FALSE,
                   return.only.var.genes = FALSE,
                   variable.features.n = 3000)
  
  cat("Completed SCTransform for sample in", sample_batch, "\n")
  return(x)
})

# Run PCA

sample_list <- lapply(sample_list, function(x) {
  x <- RunPCA(x, npcs = 50, verbose = TRUE)
})

# Select integration features with batch consideration
cat("Selecting integration features with batch awareness...\n")
features <- SelectIntegrationFeatures(object.list = sample_list, nfeatures = 3000)

# Prepare for integration
sample_list <- PrepSCTIntegration(object.list = sample_list, anchor.features = features)

# Enhanced integration with batch correction
cat("Finding integration anchors with batch correction...\n")
anchors <- FindIntegrationAnchors(object.list = sample_list, 
                                  normalization.method = "SCT",
                                  anchor.features = features,
                                  reduction = "rpca",  # More robust for batch effects
                                  k.anchor = 20,
                                  reference = c(1, 10))  # Use samples from both batches as reference

# Integrate data with batch correction
cat("Integrating data across batches...\n")
integrated_seurat <- IntegrateData(anchorset = anchors, 
                                   normalization.method = "SCT",
                                   k.weight = 100)

# Add batch information to integrated object - DEBUGGING AND FIXED METHOD
cat("Debugging metadata assignment...\n")

# Check what orig.ident values we have
cat("Unique orig.ident values in integrated object:\n")
print(unique(integrated_seurat$orig.ident))

cat("Sample IDs in sample_info:\n")
print(sample_info$sample_id)

# Create mapping from sample to batch/instrument
sample_to_batch <- setNames(sample_info$batch, sample_info$sample_id)
sample_to_instrument <- setNames(sample_info$instrument, sample_info$sample_id)

cat("Sample to batch mapping:\n")
print(sample_to_batch)

# Check if all orig.ident values have matches
missing_samples <- setdiff(unique(integrated_seurat$orig.ident), names(sample_to_batch))
if(length(missing_samples) > 0) {
  cat("WARNING: These samples in integrated object don't have batch info:\n")
  print(missing_samples)
}

tryCatch({
  integrated_seurat$batch <- sample_to_batch[integrated_seurat$orig.ident]
  integrated_seurat$instrument <- sample_to_instrument[integrated_seurat$orig.ident]
  cat("Successfully added batch and instrument information to integrated object\n")
}, error = function(e) {
  cat("Error adding metadata:", e$message, "\n")
  
  # Alternative approach - add directly to meta.data
  cat("Trying alternative approach...\n")
  
  # Get current metadata
  current_meta <- integrated_seurat@meta.data
  
  # Add new columns
  current_meta$batch <- sample_to_batch[current_meta$orig.ident]
  current_meta$instrument <- sample_to_instrument[current_meta$orig.ident]
  
  # Replace metadata
  integrated_seurat@meta.data <- current_meta
  
  cat("Successfully added metadata using alternative method\n")
})

# Clean up memory
rm(sample_list, anchors)
gc()


# Dimensionality Reduction and Clustering
# ================================================================================
cat("Performing dimensionality reduction...\n")

# Set default assay to integrated
DefaultAssay(integrated_seurat) <- "integrated"

# Run PCA
integrated_seurat <- RunPCA(integrated_seurat, 
                            npcs = 50, 
                            verbose = TRUE)

# Determine optimal number of PCs
# ElbowPlot for PC selection
elbow_plot <- ElbowPlot(integrated_seurat, ndims = 50)
print(elbow_plot)

# Use first 30 PCs based on typical scRNA-seq practices
pcs_to_use <- 30

# Run UMAP
integrated_seurat <- RunUMAP(integrated_seurat, 
                             reduction = "pca", 
                             dims = 1:pcs_to_use,
                             n.neighbors = 30,
                             min.dist = 0.3,
                             spread = 1)

# Find neighbors and clusters
integrated_seurat <- FindNeighbors(integrated_seurat, 
                                   reduction = "pca", 
                                   dims = 1:pcs_to_use,
                                   k.param = 20)

# Find clusters with resolution 0.6 as specified in methods
integrated_seurat <- FindClusters(integrated_seurat, 
                                  resolution = 0.6,
                                  algorithm = 1)  # Louvain algorithm

cat("Identified", length(unique(integrated_seurat$seurat_clusters)), "clusters\n")

save(integrated_seurat, file="/gpfs/home/rm7368/scseq_analysis/integrated_seurat.RData")

# Enhanced Visualization of Integration and Clustering with Batch Assessment
# ================================================================================
# Create comprehensive visualization plots including batch effects
create_integration_plots <- function(integrated_seurat) {
  # UMAP colored by sample
  p1 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "orig.ident", 
                pt.size = 0.5, raster = TRUE) +
    ggtitle("Integration by Sample") +
    theme(legend.position = "bottom")
  
  # UMAP colored by condition
  p2 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "condition", 
                pt.size = 0.5, raster = TRUE) +
    ggtitle("Integration by Condition") +
    theme(legend.position = "bottom")
  
  # UMAP colored by clusters
  p3 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "seurat_clusters", 
                label = TRUE, pt.size = 0.5, raster = TRUE) +
    ggtitle("Clusters") +
    theme(legend.position = "none")
  
  # UMAP colored by batch - CRITICAL for assessing batch correction
  p4 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "batch", 
                pt.size = 0.5, raster = TRUE) +
    ggtitle("Integration by Batch (Instrument)") +
    theme(legend.position = "bottom")
  
  # UMAP colored by genotype
  p5 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "genotype", 
                pt.size = 0.5, raster = TRUE) +
    ggtitle("Integration by Genotype") +
    theme(legend.position = "bottom")
  
  # Split by batch to assess integration quality
  p6 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "condition", 
                split.by = "batch", pt.size = 0.3, raster = TRUE) +
    ggtitle("Conditions Split by Batch") +
    theme(legend.position = "bottom")
  
  return(list(sample = p1, condition = p2, clusters = p3, 
              batch = p4, genotype = p5, batch_split = p6))
}

# Generate integration plots
integration_plots <- create_integration_plots(integrated_seurat)

# Display plots with special attention to batch effects
print(integration_plots$sample)
print(integration_plots$condition)
print(integration_plots$clusters)
print(integration_plots$batch)  # CRITICAL: Check if batches are well-mixed
print(integration_plots$genotype)
print(integration_plots$batch_split)  # CRITICAL: Check condition representation across batches

# Cell Type Annotation Based on Marker Genes
# ================================================================================
# Define marker genes based on your methods
marker_genes <- list(
  "Lgr5_CBC" = "Lgr5",                    # Lgr5+ CBCs
  "Paneth" = "Lyz1",                      # Lyz1+ Paneth cells
  "Isthmus_Prog" = c("Fgfbp1", "Mki67"), # Fgfbp1+Mki67+ isthmus progenitors
  "Goblet" = "Muc2",                      # Muc2+ goblet cells
  "Enteroendocrine" = "Chgb",             # Chgb+ enteroendocrine cells
  "Secretory progenitors" = "Atoh1",             # Atoh1+ secretory progenitors
  "Microfold" = "Ccl20",                  # Ccl20+ microfold cells
  "Top-villus EC" = c("Ada", "Neat1"), # Top-villus enterocytes
  "Mid-villus EC" = c("Slc5a1", "Slc2a2"), # Mid-villus enterocytes
  "Bottom-villus EC" = "Krt19"        # Bottom-villus enterocytes
)

# Function to calculate marker gene expression scores
calculate_marker_scores <- function(integrated_seurat, markers) {
  DefaultAssay(integrated_seurat) <- "SCT"
  
  for(cell_type in names(markers)) {
    genes <- markers[[cell_type]]
    # Check if genes exist in the dataset
    genes_present <- genes[genes %in% rownames(integrated_seurat)]
    
    if(length(genes_present) > 0) {
      if(length(genes_present) == 1) {
        # Single gene - use expression directly
        integrated_seurat[[paste0(cell_type, "_score")]] <- integrated_seurat@assays$SCT@data[genes_present, ]
      } else {
        # Multiple genes - use AddModuleScore
        integrated_seurat <- AddModuleScore(integrated_seurat, 
                                     features = list(genes_present),
                                     name = paste0(cell_type, "_score"),
                                     assay = "SCT")
      }
      cat("Added score for", cell_type, "using genes:", paste(genes_present, collapse = ", "), "\n")
    } else {
      cat("Warning: No genes found for", cell_type, "\n")
    }
  }
  
  return(integrated_seurat)
}

# Calculate marker scores
integrated_seurat <- calculate_marker_scores(integrated_seurat, marker_genes)

# Find cluster identity
#===============================================================
colnames(integrated_seurat@meta.data)
DimPlot(integrated_seurat, reduction = "umap", group.by = "integrated_snn_res.0.6", label = TRUE)
FeaturePlot(integrated_seurat, reduction = "umap", features = c("Ada", "Neat1"))

# Keep SCT as default
DefaultAssay(integrated_seurat) <- "SCT"

# Annotate cell types based on marker expression
# ================================================================================
# Function to annotate clusters based on marker expression

  cluster_annotations["0"] <- "Top-villus EC"
  cluster_annotations["1"] <- "Mid/top-villus EC"
  cluster_annotations["2"] <- "Mid-villus EC"
  cluster_annotations["3"] <- "Mid-villus EC"
  cluster_annotations["4"] <- "Bottom-villus EC"
  cluster_annotations["5"] <- "Mid-villus EC"
  cluster_annotations["6"] <- "Mid/top-villus EC*"
  cluster_annotations["7"] <- "Bottom-villus EC"
  cluster_annotations["8"] <- "Mid-villus EC"
  cluster_annotations["9"] <- "Isthmus progenitors"
  cluster_annotations["10"] <- "Bottom-villus EC"
  cluster_annotations["11"] <- "Lgr5+ CBC"
  cluster_annotations["12"] <- "Top-villus EC*"
  cluster_annotations["13"] <- "Mid/top-villus EC"
  cluster_annotations["14"] <- "Goblet"
  cluster_annotations["15"] <- "Bottom-villus EC"
  cluster_annotations["16"] <- "Top-villus EC"
  cluster_annotations["17"] <- "Mid-villus EC"
  cluster_annotations["18"] <- "Secretory progenitors"
  cluster_annotations["19"] <- "Mid/top-villus EC"
  cluster_annotations["20"] <- "Paneth"
  cluster_annotations["21"] <- "Enteroendocrine"
  cluster_annotations["22"] <- "Microfold"
  
integrated_seurat$cell_type <- "zero"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 0] <- "Top-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 1] <- "Mid/top-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 2] <- "Mid-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 3] <- "Mid-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 4] <- "Bottom-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 5] <- "Mid-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 6] <- "Mid/top-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 7] <- "Bottom-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 8] <- "Mid-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 9] <- "Isthmus progenitors"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 10] <- "Bottom-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 11] <- "Lgr5+ CBC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 12] <- "Top-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 13] <- "Mid/top-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 14] <- "Goblet"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 15] <- "Bottom-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 16] <- "Top-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 17] <- "Mid-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 18] <- "Secretory progenitors"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 19] <- "Mid/top-villus EC"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 20] <- "Paneth"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 21] <- "Enteroendocrine"
integrated_seurat$cell_type[integrated_seurat$integrated_snn_res.0.6 == 22] <- "Microfold"




# Pseudobulk DE
#==================================

raw_counts <- GetAssayData(integrated_seurat, layer = "counts", assay="RNA")
meta_df <- data.frame(
  cell_id = rownames(integrated_seurat@meta.data),
  condition = integrated_seurat@meta.data$condition,
  sample = integrated_seurat@meta.data$sample,
  cell_type = integrated_seurat@meta.data$cell_type
)
table(meta_df$condition, meta_df$sample)

# Aggregate counts by sample and cell type
pseudobulk_counts <- list()

for(ct in unique(meta_df$cell_type)) {
  # Get cells for this cell type
  cells_ct <- meta_df[meta_df$cell_type == ct, ]
  
  # Get counts for these cells
  counts_ct <- raw_counts[, cells_ct$cell_id]
  
  # Aggregate by sample
  samples <- unique(cells_ct$sample)
  sample_counts <- matrix(0, nrow = nrow(counts_ct), ncol = length(samples))
  rownames(sample_counts) <- rownames(counts_ct)
  colnames(sample_counts) <- samples
  
  for(s in samples) {
    cells_sample <- cells_ct[cells_ct$sample == s, "cell_id"]
    if(length(cells_sample) > 1) {
      sample_counts[, s] <- Matrix::rowSums(counts_ct[, cells_sample])
    } else {
      sample_counts[, s] <- counts_ct[, cells_sample]
    }
  }
  
  pseudobulk_counts[[ct]] <- sample_counts
}

# Check one cell type
names(pseudobulk_counts)
dim(pseudobulk_counts[[1]])

# Create sample metadata
sample_meta <- unique(meta_df[, c("sample", "condition")])
rownames(sample_meta) <- sample_meta$sample

# Update the DESeq2 function for WT comparison
pseudobulk_results_wt <- list()

for(i in 1:length(pseudobulk_counts)) {
  celltype <- names(pseudobulk_counts)[i]
  cat("Processing pseudobulk for:", celltype, "\n")
  
  # Get samples for WT conditions only
  wt_samples <- sample_meta[sample_meta$condition %in% c("WT_Lactating", "WT_Virgin"), ]
  counts_wt <- pseudobulk_counts[[i]][, rownames(wt_samples)]
  res <- run_deseq_celltype(counts_wt, celltype, c("WT_Lactating", "WT_Virgin"))
  pseudobulk_results_wt[[celltype]] <- res
}

length(pseudobulk_results_wt)

# Check significant genes (FDR < 0.05) across all cell types
sig_counts_wt <- sapply(pseudobulk_results_wt, function(res) {
  res_df <- as.data.frame(res)
  sum(res_df$padj < 0.05 & !is.na(res_df$padj))
})

print(sig_counts_wt)

# Also check nominal significance (p < 0.05)
nominal_sig_counts_wt <- sapply(pseudobulk_results_wt, function(res) {
  res_df <- as.data.frame(res)
  sum(res_df$pvalue < 0.05 & !is.na(res_df$pvalue))
})

print(nominal_sig_counts_wt)

# Look at the most significant results across cell types
for(celltype in names(pseudobulk_results_wt)) {
  res_df <- as.data.frame(pseudobulk_results_wt[[celltype]])
  res_df <- res_df[order(res_df$pvalue), ]
  cat("\nTop genes in", celltype, ":\n")
  print(head(res_df[, c("log2FoldChange", "pvalue", "padj")], 3))
}

# Create directory for results
dir.create("sc_pseudobulk_results", showWarnings = FALSE)

# Function to create volcano plot
make_volcano_plot <- function(deseq_result, celltype_name, padj_threshold = 0.05, fc_threshold = 1, n_labels = 50) {
  res_df <- as.data.frame(deseq_result)
  res_df <- res_df[!is.na(res_df$padj), ]
  res_df$gene <- rownames(res_df) 
  
  # Add significance categories
  res_df$significance <- "Not Significant"
  res_df$significance[res_df$padj < padj_threshold & res_df$log2FoldChange > fc_threshold] <- "Up in Lactating"
  res_df$significance[res_df$padj < padj_threshold & res_df$log2FoldChange < -fc_threshold] <- "Up in Virgin"
  
  top_genes <- res_df[order(res_df$pvalue), ][1:n_labels, ]
  
  # Create plot
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("Up in Lactating" = "red", 
                                  "Up in Virgin" = "blue", 
                                  "Not Significant" = "grey")) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", alpha = 0.5) +
    geom_text_repel(data = top_genes, 
                    aes(label = gene), 
                    size = 3, 
                    max.overlaps = 40,
                    color = "black") +
    labs(title = paste("Volcano Plot:", celltype_name),
         subtitle = "WT_Lactating vs WT_Virgin",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Save plot and CSV before returning
  clean_name <- gsub("[^A-Za-z0-9_]", "_", celltype_name)
  filename <- paste0("sc_pseudobulk_results/volcano_", clean_name, ".pdf")
  ggsave(filename, p, width = 10, height = 8, dpi = 300)
  cat("Saved plot:", filename, "\n")
  
  # Save CSV
  res_df_export <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  res_df_export <- res_df_export[order(res_df_export$pvalue), ]
  csv_filename <- paste0("sc_pseudobulk_results/", clean_name, "_WT_Lactating_vs_WT_Virgin.csv")
  write.csv(res_df_export, csv_filename, row.names = FALSE)
  cat("Saved CSV:", csv_filename, "\n")
  
  return(p)
}

# Create volcano plots for all cell types
volcano_plots <- list()
for(celltype in names(pseudobulk_results_wt)) {
  volcano_plots[[celltype]] <- make_volcano_plot(pseudobulk_results_wt[[celltype]], celltype)
}

# Sample aggregation DE
#====================================

# First, create the whole-sample aggregated counts matrix
whole_sample_counts <- matrix(0, nrow = nrow(raw_counts), ncol = length(unique(meta_df$sample)))
rownames(whole_sample_counts) <- rownames(raw_counts)
colnames(whole_sample_counts) <- unique(meta_df$sample)

for(s in unique(meta_df$sample)) {
  cells_sample <- meta_df[meta_df$sample == s, "cell_id"]
  if(length(cells_sample) > 1) {
    whole_sample_counts[, s] <- Matrix::rowSums(raw_counts[, cells_sample])
  } else {
    whole_sample_counts[, s] <- raw_counts[, cells_sample]
  }
}

# Now get WT samples only
wt_sample_meta <- sample_meta[sample_meta$condition %in% c("WT_Lactating", "WT_Virgin"), ]
wt_counts <- whole_sample_counts[, rownames(wt_sample_meta)]

# Check dimensions
dim(wt_counts)
table(wt_sample_meta$condition)

# Run DESeq2 on whole-sample aggregated data
keep_genes <- rowSums(wt_counts >= 10) >= 3
wt_counts_filtered <- wt_counts[keep_genes, ]

# Create DESeq2 object for whole-sample analysis
dds_whole <- DESeqDataSetFromMatrix(
  countData = wt_counts_filtered,
  colData = wt_sample_meta,
  design = ~ condition
)

# Run DESeq2
dds_whole <- DESeq(dds_whole)

# Get results
whole_sample_result <- results(dds_whole, contrast = c("condition", "WT_Lactating", "WT_Virgin"))

summary(whole_sample_result)

# Create volcano plot function for whole-sample
make_volcano_plot_whole <- function(deseq_result, analysis_name, padj_threshold = 0.05, fc_threshold = 1, n_labels = 50) {
  res_df <- as.data.frame(deseq_result)
  res_df <- res_df[!is.na(res_df$padj), ]
  res_df$gene <- rownames(res_df) 
  
  # Add significance categories
  res_df$significance <- "Not Significant"
  res_df$significance[res_df$padj < padj_threshold & res_df$log2FoldChange > fc_threshold] <- "Up in Lactating"
  res_df$significance[res_df$padj < padj_threshold & res_df$log2FoldChange < -fc_threshold] <- "Up in Virgin"
  
  top_genes <- res_df[order(res_df$pvalue), ][1:n_labels, ]
  
  # Create plot
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("Up in Lactating" = "red", 
                                  "Up in Virgin" = "blue", 
                                  "Not Significant" = "grey")) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", alpha = 0.5) +
    geom_text_repel(data = top_genes, 
                    aes(label = gene), 
                    size = 3, 
                    max.overlaps = 40,
                    color = "black") +
    labs(title = paste("Volcano Plot:", analysis_name),
         subtitle = "WT_Lactating vs WT_Virgin",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Save plot and CSV
  filename <- paste0("sc_pseudobulk_results/volcano_", analysis_name, ".pdf")
  ggsave(filename, p, width = 10, height = 8, dpi = 300)
  cat("Saved plot:", filename, "\n")
  
  # Save CSV
  res_df_export <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  res_df_export <- res_df_export[order(res_df_export$pvalue), ]
  csv_filename <- paste0("sc_pseudobulk_results/", analysis_name, "_WT_Lactating_vs_WT_Virgin.csv")
  write.csv(res_df_export, csv_filename, row.names = FALSE)
  cat("Saved CSV:", csv_filename, "\n")
  
  return(p)
}

# Create whole-sample volcano plot
whole_sample_plot <- make_volcano_plot_whole(whole_sample_result, "Whole_Sample_Analysis")
print(whole_sample_plot)

library(Seurat)
library(biomaRt)
library(dplyr)

# Function to convert ENSMUSG IDs to gene symbols
convert_ensmusg_to_symbols <- function(seurat_obj) {
  
  # Get all unique ENSMUSG IDs from the Seurat object
  all_genes <- rownames(seurat_obj)
  ensmusg_ids <- unique(all_genes[grepl("^ENSMUSG", all_genes)])
  
  cat("Found", length(ensmusg_ids), "ENSMUSG IDs to convert\n")
  
  # Connect to biomaRt
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  # Query biomaRt for gene symbols
  # Split into chunks if there are many genes (biomaRt has query limits)
  chunk_size <- 500
  n_chunks <- ceiling(length(ensmusg_ids) / chunk_size)
  
  conversion_df <- data.frame()
  
  for (i in 1:n_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, length(ensmusg_ids))
    chunk_ids <- ensmusg_ids[start_idx:end_idx]
    
    cat("Processing chunk", i, "of", n_chunks, "\n")
    
    # Get conversions for this chunk
    chunk_result <- getBM(
      attributes = c("ensembl_gene_id", "mgi_symbol"),
      filters = "ensembl_gene_id",
      values = chunk_ids,
      mart = mart
    )
    
    conversion_df <- rbind(conversion_df, chunk_result)
  }
  
  # Create conversion mapping
  conversion_map <- setNames(conversion_df$mgi_symbol, conversion_df$ensembl_gene_id)
  
  # Remove empty symbols
  conversion_map <- conversion_map[conversion_map != ""]
  
  cat("Successfully mapped", length(conversion_map), "genes\n")
  
  # Get current gene names
  old_names <- rownames(seurat_obj)
  new_names <- old_names
  
  # Replace ENSMUSG IDs with symbols where available
  for (i in seq_along(old_names)) {
    if (old_names[i] %in% names(conversion_map)) {
      new_names[i] <- conversion_map[old_names[i]]
    }
  }
  
  # Check for duplicates after conversion
  if (any(duplicated(new_names))) {
    cat("Warning: Found", sum(duplicated(new_names)), "duplicate gene symbols after conversion\n")
    
    # Handle duplicates by keeping ENSMUSG ID for duplicates
    dup_symbols <- unique(new_names[duplicated(new_names)])
    for (symbol in dup_symbols) {
      indices <- which(new_names == symbol)
      # Keep the first occurrence as symbol, revert others to ENSMUSG
      if (length(indices) > 1) {
        new_names[indices[-1]] <- old_names[indices[-1]]
      }
    }
  }
  
  # For Seurat v5, we need to create a new object with updated feature names
  # First, extract all the data
  assay_list <- list()
  assay_names <- names(seurat_obj@assays)
  
  for (assay_name in assay_names) {
    cat("Processing assay:", assay_name, "\n")
    
    # Get the current assay
    assay <- seurat_obj@assays[[assay_name]]
    
    if (inherits(assay, "Assay5")) {
      # Extract layers
      layer_names <- names(assay@layers)
      layers_data <- list()
      
      for (layer in layer_names) {
        if (!is.null(assay@layers[[layer]])) {
          mat <- assay@layers[[layer]]
          rownames(mat) <- new_names
          layers_data[[layer]] <- mat
        }
      }
      
      # Create new v5 assay with updated feature names
      if ("counts" %in% names(layers_data)) {
        new_assay <- CreateAssay5Object(counts = layers_data[["counts"]])
        
        # Add other layers
        for (layer in setdiff(names(layers_data), "counts")) {
          new_assay@layers[[layer]] <- layers_data[[layer]]
        }
      } else if ("data" %in% names(layers_data)) {
        new_assay <- CreateAssay5Object(data = layers_data[["data"]])
        
        # Add other layers
        for (layer in setdiff(names(layers_data), "data")) {
          new_assay@layers[[layer]] <- layers_data[[layer]]
        }
      }
      
      # Copy over metadata if it exists
      if (!is.null(assay@meta.data) && nrow(assay@meta.data) > 0) {
        meta_data <- assay@meta.data
        rownames(meta_data) <- new_names
        new_assay@meta.data <- meta_data
      }
      
      # Copy over variable features if they exist
      var_features <- VariableFeatures(assay)
      if (length(var_features) > 0) {
        # Update variable feature names
        var_features_updated <- var_features
        for (i in seq_along(var_features)) {
          if (var_features[i] %in% old_names) {
            idx <- which(old_names == var_features[i])
            var_features_updated[i] <- new_names[idx]
          }
        }
        VariableFeatures(new_assay) <- var_features_updated
      }
      
      assay_list[[assay_name]] <- new_assay
      
    } else {
      # Handle older assay versions
      counts <- GetAssayData(assay, slot = "counts")
      rownames(counts) <- new_names
      
      data <- GetAssayData(assay, slot = "data")
      rownames(data) <- new_names
      
      new_assay <- CreateAssayObject(counts = counts, data = data)
      
      if (nrow(GetAssayData(assay, slot = "scale.data")) > 0) {
        scale_data <- GetAssayData(assay, slot = "scale.data")
        rownames(scale_data) <- new_names
        SetAssayData(new_assay, slot = "scale.data", new.data = scale_data)
      }
      
      assay_list[[assay_name]] <- new_assay
    }
  }
  
  # Create new Seurat object with updated assays
  new_seurat <- CreateSeuratObject(
    counts = assay_list[[1]], 
    assay = names(assay_list)[1],
    meta.data = seurat_obj@meta.data
  )
  
  # Add additional assays if present
  if (length(assay_list) > 1) {
    for (i in 2:length(assay_list)) {
      new_seurat[[names(assay_list)[i]]] <- assay_list[[i]]
    }
  }
  
  # Set default assay to match original
  DefaultAssay(new_seurat) <- DefaultAssay(seurat_obj)
  
  # Copy over reductions
  for (reduc_name in names(seurat_obj@reductions)) {
    new_seurat@reductions[[reduc_name]] <- seurat_obj@reductions[[reduc_name]]
  }
  
  # Copy over graphs
  for (graph_name in names(seurat_obj@graphs)) {
    new_seurat@graphs[[graph_name]] <- seurat_obj@graphs[[graph_name]]
  }
  
  # Copy over images if they exist (for spatial data)
  if (length(seurat_obj@images) > 0) {
    for (image_name in names(seurat_obj@images)) {
      new_seurat@images[[image_name]] <- seurat_obj@images[[image_name]]
    }
  }
  
  # Copy over any project info
  new_seurat@project.name <- seurat_obj@project.name
  
  # Report conversion statistics
  n_converted <- sum(old_names != new_names)
  n_unchanged <- sum(old_names == new_names)
  
  cat("\nConversion complete:\n")
  cat("- Converted:", n_converted, "genes\n")
  cat("- Unchanged:", n_unchanged, "genes\n")
  cat("- Failed to map:", length(ensmusg_ids) - length(conversion_map), "ENSMUSG IDs\n")
  
  return(new_seurat)
}

merged_xenium <- convert_ensmusg_to_symbols(merged_xenium)

# Alternative: If biomaRt is slow or unavailable, use a local annotation package
convert_with_annotationdbi <- function(seurat_obj) {
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  
  # Get all unique ENSMUSG IDs
  all_genes <- rownames(seurat_obj)
  ensmusg_ids <- unique(all_genes[grepl("^ENSMUSG", all_genes)])
  
  cat("Found", length(ensmusg_ids), "ENSMUSG IDs to convert\n")
  
  # Map ENSMUSG to gene symbols
  symbols <- mapIds(org.Mm.eg.db,
                    keys = ensmusg_ids,
                    column = "SYMBOL",
                    keytype = "ENSEMBL",
                    multiVals = "first")
  
  # Create conversion map
  conversion_map <- symbols[!is.na(symbols)]
  
  cat("Successfully mapped", length(conversion_map), "genes\n")
  
  # Apply same conversion logic as above
  old_names <- rownames(seurat_obj)
  new_names <- old_names
  
  for (i in seq_along(old_names)) {
    if (old_names[i] %in% names(conversion_map)) {
      new_names[i] <- conversion_map[old_names[i]]
    }
  }
  
  # Handle duplicates
  if (any(duplicated(new_names))) {
    cat("Warning: Found", sum(duplicated(new_names)), "duplicate gene symbols after conversion\n")
    
    dup_symbols <- unique(new_names[duplicated(new_names)])
    for (symbol in dup_symbols) {
      indices <- which(new_names == symbol)
      if (length(indices) > 1) {
        new_names[indices[-1]] <- old_names[indices[-1]]
      }
    }
  }
  
  # Use the same reconstruction logic as biomaRt version
  # ... (rest of the code is the same as in convert_ensmusg_to_symbols from line 65 onwards)
  
  return(new_seurat)
}
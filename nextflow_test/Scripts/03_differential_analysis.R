#!/usr/bin/env Rscript
################################################################################
# 03_differential_analysis.R
#
# Differential Analysis Pipeline:
# 1. DA - Differential Abundance (MiloR)
# 2. DS - Differential State (Seurat + Muscat)
################################################################################

################################################################################
# SETUP
################################################################################

# Load core libraries
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# Load optional libraries
safe_lib <- function(pkgs) {
  for (p in pkgs) {
    suppressWarnings(suppressMessages(
      require(p, character.only = TRUE, quietly = TRUE)
    ))
  }
}

safe_lib(c(
  "plotly", "htmlwidgets", "pheatmap", "RColorBrewer", "tidyr",
  "muscat", "miloR", "SingleCellExperiment", "edgeR", 
  "patchwork", "UpSetR", "scater", "tibble"
))

# Source utility functions
# Source bootstrap utility to get script directory
source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))), "utils", "bootstrap.R"))
script_dir <- get_script_dir()
source(file.path(script_dir, "utils", "common.R"))
source(file.path(script_dir, "utils", "performance.R"))
source(file.path(script_dir, "utils", "03_differential_plots.R"))
source(file.path(script_dir, "utils", "pvalue_helpers.R"))

################################################################################
# UTILITY FUNCTIONS
################################################################################

# Note: All plot functions have been moved to utils/03_differential_plots.R

################################################################################
# DIFFERENTIAL ABUNDANCE (DA) - MILOR
################################################################################

#' Calculate optimal MiloR parameters with iterative validation
#' 
#' Uses empirically recommended parameters and validates neighborhood sizes.
#' Iteratively adjusts k and prop if neighborhood size distribution is suboptimal.
#' Target: peak mean between 50-100 and average > 5 x N_samples
#' 
#' @param seu Seurat object
#' @param annotation_col Cell type annotation column name
#' @param default_k Default k neighbors (from buildGraph)
#' @param default_prop Default proportion of cells for neighborhoods (0.1-0.2)
#' @param k_graph Number of PCA dimensions to use
#' @param max_iterations Maximum number of parameter adjustment iterations
#' @param log Logging function
#' @return List with recommended parameters
calculate_optimal_milo_params <- function(seu, annotation_col = "cell.type", 
                                         default_k = 10, default_prop = 0.1, 
                                         k_graph = 30, max_iterations = 3, 
                                         log = cat) {
  
  log("=== MiloR Parameter Selection (Using Empirically Recommended Values) ===")
  
  # Check PCA availability
  if (!"pca" %in% names(seu@reductions)) {
    return(list(k = default_k, prop = default_prop, d = k_graph, optimized = FALSE))
  }
  
  # Get dimensions from step 02 (actual PCA dimensions in the object)
  d_optimal <- min(k_graph, ncol(seu@reductions$pca))
  log(paste("Using d =", d_optimal, "PCs (from", ncol(seu@reductions$pca), "available in Seurat object)"))
  
  # Get cell count and adjust prop based on dataset size
  n_cells <- ncol(seu)
  log(paste("Dataset size:", n_cells, "cells"))
  
  # Adjust prop according to MiloR recommendations:
  # prop=0.1 for <30k cells, prop=0.05 for >=30k cells
  if (n_cells < 30000) {
    prop_optimal <- 0.1
    log("Using prop=0.1 (dataset < 30k cells)")
  } else {
    prop_optimal <- 0.05
    log("Using prop=0.05 (dataset >= 30k cells, faster computation)")
  }
  
  # Get cell type info for logging
  cell_type_table <- if (annotation_col %in% colnames(seu@meta.data)) table(seu@meta.data[[annotation_col]][!is.na(seu@meta.data[[annotation_col]])]) else NULL
  smallest_cluster_size <- if (!is.null(cell_type_table)) min(cell_type_table) else NA
  if (!is.na(smallest_cluster_size)) {
    log(paste("Cell types:", length(cell_type_table), "| Smallest:", smallest_cluster_size, "cells"))
    
    # Warning if smallest cluster is too small
    if (smallest_cluster_size < 10) {
      log(paste("WARNING: Very small cell type detected (", smallest_cluster_size, "cells). This may cause issues with neighbourhood construction."))
      log("Consider filtering out rare cell types before DA analysis.")
    }
  }
  
  # Use k from defaults (should match the k used in graph building)
  k_optimal <- default_k
  
  # Determine number of samples for validation
  n_samples <- if ("sample" %in% colnames(seu@meta.data)) {
    length(unique(seu@meta.data$sample))
  } else if ("orig.ident" %in% colnames(seu@meta.data)) {
    length(unique(seu@meta.data$orig.ident))
  } else {
    NA
  }
  
  min_avg_nhood_size <- if (!is.na(n_samples)) 5 * n_samples else 50
  
  log(paste("Initial parameters: k =", k_optimal, ", prop =", prop_optimal, ", d =", d_optimal))
  if (!is.na(n_samples)) {
    log(paste("Number of samples:", n_samples, "| Minimum avg neighborhood size:", min_avg_nhood_size))
  }
  log("Expected neighbourhood size distribution: peak mean between 50-100")
  log("Rule of thumb: Average neighbourhood size should be > 5 x N_samples")
  

  
  return(list(
    k = k_optimal,
    prop = prop_optimal,
    d = d_optimal,
    min_avg_nhood_size = min_avg_nhood_size,
    n_samples = n_samples,
    optimized = TRUE
  ))
}

#' Validate MiloR neighborhood parameters
#' 
#' Checks quality of neighborhoods by analyzing purity and cell type coverage
#' 
#' @param milo Milo object
#' @param seu Seurat object
#' @param annotation_col Cell type annotation column
#' @param log Logging function
#' @return List with validation results
validate_milor_neighborhoods <- function(milo, seu, annotation_col = "cell.type", log = cat) {
  
  log("=== Validating MiloR Neighborhood Parameters ===")
  
  # Check if annotation column exists
  if (!annotation_col %in% colnames(seu@meta.data)) return(list(validation_passed = FALSE, reason = "missing_annotation"))
  
  # Extract neighborhoods and annotations
  nhoods <- nhoods(milo)
  cell_annotations <- seu@meta.data[[annotation_col]]
  names(cell_annotations) <- colnames(seu)
  
  n_nhoods <- ncol(nhoods)
  log(paste("Analyzing", n_nhoods, "neighborhoods..."))
  
  # Calculate purity for each neighborhood
  purity_scores <- sapply(1:n_nhoods, function(i) {
    nhood_cells <- rownames(nhoods)[nhoods[, i] > 0]
    if (length(nhood_cells) == 0) return(NA)
    
    nhood_annot <- cell_annotations[nhood_cells]
    nhood_annot <- nhood_annot[!is.na(nhood_annot)]
    if (length(nhood_annot) == 0) return(NA)
    
    annot_table <- table(nhood_annot)
    return(max(annot_table) / length(nhood_annot))
  })
  
  # Remove NAs
  purity_scores <- purity_scores[!is.na(purity_scores)]
  
  # Calculate summary statistics
  mean_purity <- mean(purity_scores)
  median_purity <- median(purity_scores)
  pct_high_purity <- sum(purity_scores >= 0.8) / length(purity_scores) * 100
  
  log(paste("Neighborhood Purity Statistics:"))
  log(paste("  Mean purity:", round(mean_purity, 3)))
  log(paste("  Median purity:", round(median_purity, 3)))
  log(paste("  % neighborhoods with purity >= 80%:", round(pct_high_purity, 1), "%"))
  
  # Assign majority annotation to each neighborhood
  nhood_annotations <- sapply(1:n_nhoods, function(i) {
    nhood_cells <- rownames(nhoods)[nhoods[, i] > 0]
    if (length(nhood_cells) == 0) return(NA)
    
    nhood_annot <- cell_annotations[nhood_cells][!is.na(cell_annotations[nhood_cells])]
    if (length(nhood_annot) == 0) return(NA)
    
    names(which.max(table(nhood_annot)))
  })
  
  # Check coverage: are all cell types represented?
  unique_cell_types <- unique(cell_annotations[!is.na(cell_annotations)])
  covered_cell_types <- unique(nhood_annotations[!is.na(nhood_annotations)])
  missing_cell_types <- setdiff(unique_cell_types, covered_cell_types)
  
  coverage_pct <- length(covered_cell_types) / length(unique_cell_types) * 100
  
  log(paste("Cell Type Coverage:"))
  log(paste("  Total annotated cell types:", length(unique_cell_types)))
  log(paste("  Cell types represented in neighborhoods:", length(covered_cell_types)))
  log(paste("  Coverage:", round(coverage_pct, 1), "%"))
  
  if (length(missing_cell_types) > 0) {
    log(paste("  Missing cell types:", paste(missing_cell_types, collapse=", ")))
  }
  
  # Check smallest cluster sizes
  cell_type_counts <- table(cell_annotations)
  smallest_cluster_size <- min(cell_type_counts)
  smallest_cluster_name <- names(cell_type_counts)[which.min(cell_type_counts)]
  
  log(paste("Smallest annotated cluster:", smallest_cluster_name, "with", smallest_cluster_size, "cells"))
  
  validation_passed <- mean_purity >= 0.7 && coverage_pct >= 90
  log(paste(ifelse(validation_passed, "✓ PASSED", "✗ FAILED"), "- Purity:", round(mean_purity, 2), "Coverage:", round(coverage_pct, 1), "%"))
  
  return(list(
    validation_passed = validation_passed,
    mean_purity = mean_purity,
    median_purity = median_purity,
    pct_high_purity = pct_high_purity,
    purity_scores = purity_scores,
    coverage_pct = coverage_pct,
    covered_cell_types = covered_cell_types,
    missing_cell_types = missing_cell_types,
    smallest_cluster_name = smallest_cluster_name,
    smallest_cluster_size = smallest_cluster_size
  ))
}

#' Optimize MiloR parameters with validation
#' 
#' Determines optimal parameters and validates them before use
#' 
#' @param seu Seurat object
#' @param milo_params List of MiloR parameters (may be updated)
#' @param analysis_params Analysis parameters from previous steps
#' @param log Logging function
#' @return List with optimized parameters and optimization result
optimize_and_validate_milo_params <- function(seu, milo_params, analysis_params = NULL, log = cat) {
  
  # Find annotation column
  annotation_col <- c("cell.type", "celltype", "seurat_clusters")[
    c("cell.type", "celltype", "seurat_clusters") %in% colnames(seu@meta.data)
  ][1]
  
  if (is.na(annotation_col) || annotation_col == "seurat_clusters") {
    log("Warning: No cell.type annotation found, using seurat_clusters for optimization")
  }
  
  param_optimization <- if (!is.na(annotation_col)) {
    log("\n--- Optimizing MiloR Parameters Based on Data ---")
    
    # Use k_neighbors from step 02 if available
    k_graph_value <- if (!is.null(analysis_params) && !is.null(analysis_params$k_neighbors)) {
      analysis_params$k_neighbors
    } else {
      milo_params$k %||% 30
    }
    
    log(paste("Using k_graph =", k_graph_value, "from", 
              if (!is.null(analysis_params$k_neighbors)) "step 02" else "defaults"))
    
    # Calculate optimal parameters
    opt <- calculate_optimal_milo_params(
      seu = seu,
      annotation_col = annotation_col,
      default_k = milo_params$k %||% 10,
      default_prop = milo_params$prop %||% 0.1,
      k_graph = k_graph_value,
      log = log
    )
    
    if (opt$optimized) {
      log(paste("\nUsing optimized parameters instead of defaults:"))
      log(paste("  k:", milo_params$k %||% 10, "→", opt$k))
      log(paste("  prop:", milo_params$prop %||% 0.1, "→", opt$prop))
      
      # Validate optimized parameters before using them
      if (!is.null(opt$k) && !is.na(opt$k) && is.numeric(opt$k) && opt$k > 0) {
        milo_params$k <- opt$k
      } else {
        log("WARNING: Optimized k is invalid, keeping default")
      }
      
      if (!is.null(opt$prop) && !is.na(opt$prop) && is.numeric(opt$prop) && 
          opt$prop > 0 && opt$prop <= 1) {
        milo_params$prop <- opt$prop
      } else {
        log("WARNING: Optimized prop is invalid, keeping default")
      }
    }
    opt
  } else {
    log("Warning: No annotation available - using default parameters without optimization")
    NULL
  }
  
  return(list(
    milo_params = milo_params,
    optimization_result = param_optimization
  ))
}

# Note: create_neighborhood_validation_plots() moved to utils/03_differential_plots.R

#' Build Milo object from Seurat
#' 
#' Creates Milo object with precomputed KNN graph from Seurat
#' 
#' @param seu Seurat object
#' @param milo_k K neighbors
#' @param milo_d Number of dimensions
#' @param milo_prop Proportion for neighborhoods
#' @param milo_refined Use refined neighborhoods
#' @param log Logging function
#' @return Milo object with neighborhoods
build_milo_object <- function(seu, milo_k, milo_d, milo_prop, milo_refined, log = cat) {
  
  # Extract counts
  count_matrix <- GetAssayData(seu, assay = "RNA", slot = "counts")
  log(paste("Count matrix:", nrow(count_matrix), "genes x", ncol(count_matrix), "cells"))
  
  # Check PCA availability
  has_pca <- "pca" %in% names(seu@reductions)
  
  # Create SingleCellExperiment with reductions
  sce_temp <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = count_matrix))
  
  if (has_pca) {
    n_pcs <- ncol(seu@reductions$pca)
    milo_dims <- min(milo_d, n_pcs)
    log(paste0("Using ", milo_dims, " PCs from step 02"))
    SingleCellExperiment::reducedDim(sce_temp, "PCA") <- Embeddings(seu, "pca")[, 1:milo_dims]
    if ("umap" %in% names(seu@reductions)) {
      SingleCellExperiment::reducedDim(sce_temp, "UMAP") <- Embeddings(seu, "umap")
    }
  } else {
    milo_dims <- milo_d
  }
  
  # Initialize Milo object
  milo <- Milo(sce_temp)
  
  # Use precomputed Seurat KNN graph if available
  if ("RNA_nn" %in% names(seu@graphs)) {
    log("Using precomputed Seurat KNN graph (from FindNeighbors in step 02)")
    knn_graph <- seu@graphs$RNA_nn
    milo <- buildFromAdjacency(milo, knn_graph, k = milo_k)
    log(paste0("Graph built with k=", milo_k, " from Seurat"))
  } else {
    log("No precomputed graph found, building new graph")
    milo <- buildGraph(milo, k = milo_k, d = milo_dims)
  }
  
  # Build neighborhoods
  if (has_pca) {
    milo <- makeNhoods(milo, prop = milo_prop, k = milo_k, d = milo_dims, 
                      reduced_dim = "PCA", refined = milo_refined)
  } else {
    milo <- makeNhoods(milo, prop = milo_prop, k = milo_k, d = milo_d, 
                      refined = milo_refined)
  }
  
  # Validate neighborhood sizes
  nhood_sizes <- colSums(nhoods(milo))
  mean_nhood_size <- mean(nhood_sizes)
  median_nhood_size <- median(nhood_sizes)
  
  log(paste("Neighborhood size statistics:"))
  log(paste("  Mean:", round(mean_nhood_size, 1)))
  log(paste("  Median:", round(median_nhood_size, 1)))
  log(paste("  Min:", min(nhood_sizes)))
  log(paste("  Max:", max(nhood_sizes)))
  
  # Check if sizes are in recommended range (50-100 peak)
  if (mean_nhood_size < 50) {
    log("WARNING: Mean neighborhood size < 50. Consider increasing k or prop.")
  } else if (mean_nhood_size > 100) {
    log("WARNING: Mean neighborhood size > 100. Consider decreasing k or prop.")
  } else {
    log("✓ Neighborhood sizes are in recommended range (50-100)")
  }
  
  milo <- buildNhoodGraph(milo)
  
  return(milo)
}

#' Setup sample metadata and design for MiloR
#' 
#' Creates sample metadata and design matrix for differential testing
#' 
#' @param seu Seurat object
#' @param cond_field Condition field name
#' @param log Logging function
#' @return List with design_df, sample_metadata, and design matrix
setup_milo_design <- function(seu, cond_field, log = cat) {
  
  # Determine sample column
  sample_col <- if ("sample" %in% colnames(seu@meta.data)) {
    "sample"
  } else if ("orig.ident" %in% colnames(seu@meta.data)) {
    "orig.ident"
  } else {
    seu$milo_sample <- paste0(
      seu@meta.data[[cond_field]], 
      "_rep", 
      ave(seq_len(ncol(seu)), seu@meta.data[[cond_field]], FUN = seq_along)
    )
    "milo_sample"
  }
  
  # Sanitize condition names to be valid R names
  sanitized_conditions <- make.names(seu@meta.data[[cond_field]])
  
  # Create mapping from original to sanitized names for later use
  condition_name_mapping <- unique(data.frame(
    original = seu@meta.data[[cond_field]],
    sanitized = sanitized_conditions
  ))
  
  # Design dataframe (keep original condition labels for reporting)
  design_df <- data.frame(
    sample = seu@meta.data[[sample_col]],
    condition = sanitized_conditions,
    condition_orig = seu@meta.data[[cond_field]]
  )
  rownames(design_df) <- colnames(seu)
  
  # Sample metadata
  sample_metadata <- design_df %>% distinct(sample, condition, condition_orig)
  rownames(sample_metadata) <- sample_metadata$sample
  
  # Design matrix - use ~0 + condition for easier contrast specification
  design <- model.matrix(~0 + condition, data = sample_metadata)
  
  log(paste("Design: ", nrow(sample_metadata), "samples,", 
           length(unique(sample_metadata$condition)), "conditions"))
  
  return(list(
    design_df = design_df,
    sample_metadata = sample_metadata,
    design = design,
    sample_col = sample_col,
    condition_name_mapping = condition_name_mapping
  ))
}

#' Perform MiloR Differential Abundance Analysis
#' 
#' Main wrapper function for MiloR DA analysis
#' 
#' @param seu Seurat object
#' @param cond_field Condition field name
#' @param milo_params MiloR parameters
#' @param analysis_params Analysis parameters from previous steps
#' @param log Logging function
#' @return List with milo object, DA results, and metadata
perform_milor_analysis <- function(seu, cond_field, milo_params, analysis_params = NULL, log = cat) {
  
  if (!"miloR" %in% rownames(installed.packages())) {
    log("MiloR not installed, skipping DA analysis")
    return(list(milo = NULL, da = NULL, sample_metadata = NULL))
  }
  
  log("Starting MiloR DA analysis...")
  
  # Optimize and validate parameters
  opt_result <- optimize_and_validate_milo_params(
    seu = seu,
    milo_params = milo_params,
    analysis_params = analysis_params,
    log = log
  )
  
  milo_params <- opt_result$milo_params
  param_optimization <- opt_result$optimization_result
  
  # Extract parameters (now potentially optimized)
  milo_use_pca <- milo_params$use_pca %||% TRUE
  milo_k <- milo_params$k %||% 10
  milo_d <- milo_params$d %||% 30
  milo_prop <- milo_params$prop %||% 0.1
  milo_refined <- milo_params$refined %||% TRUE
  milo_fdr <- milo_params$fdr %||% 0.1
  milo_fdr_weighting <- milo_params$fdr_weighting %||% "none"
  
  log("=== Final MiloR parameters after optimization and validation ===")
  log(paste("  - k (neighbors):", milo_k))
  log(paste("  - prop:", milo_prop))
  log(paste("  - d (dimensions):", milo_d))
  log(paste("  - Using PCA:", milo_use_pca))
  log(paste("  - refined:", milo_refined))
  log(paste("  - FDR threshold:", milo_fdr))
  log(paste("  - FDR weighting:", milo_fdr_weighting))
  log("===============================================================")
  
  # Build Milo object
  log("Building Milo object...")
  milo <- build_milo_object(
    seu = seu,
    milo_k = milo_k,
    milo_d = milo_d,
    milo_prop = milo_prop,
    milo_refined = milo_refined,
    log = log
  )
  
  # Validate neighborhoods
  annotation_col <- c("cell.type", "celltype", "seurat_clusters")[
    c("cell.type", "celltype", "seurat_clusters") %in% colnames(seu@meta.data)
  ][1]
  
  validation_results <- if (!is.na(annotation_col)) {
    validate_milor_neighborhoods(milo, seu, annotation_col, log)
  } else {
    NULL
  }
  
  # Setup design
  log("Setting up design matrix...")
  design_setup <- setup_milo_design(seu, cond_field, log)
  
  # Count cells
  milo <- countCells(milo, meta.data = design_setup$design_df, sample = "sample")
  
  # Test neighborhoods with pairwise contrasts
  log("Testing neighborhoods (pairwise comparisons)...")
  
  # Get all unique conditions (original labels for reporting)
  conditions <- unique(design_setup$sample_metadata$condition_orig)
  log(paste("Found", length(conditions), "conditions:", paste(conditions, collapse=", ")))
  
  # Perform pairwise DA tests
  da_results <- list()
  
  for (i in 1:(length(conditions)-1)) {
    for (j in (i+1):length(conditions)) {
      cond1 <- conditions[i]
      cond2 <- conditions[j]
      contrast_name <- paste0(cond1, "_vs_", cond2)
      
      log(paste("  Testing:", contrast_name))
      
      # Create contrast string for model.contrasts parameter
      # Format: "condition2 - condition1" to match column names in design matrix
      cond1_sanitized <- design_setup$condition_name_mapping$sanitized[
        match(cond1, design_setup$condition_name_mapping$original)
      ]
      cond2_sanitized <- design_setup$condition_name_mapping$sanitized[
        match(cond2, design_setup$condition_name_mapping$original)
      ]
      contrast_str <- paste0("condition", cond2_sanitized, " - condition", cond1_sanitized)
      
      # Test this contrast using model.contrasts
      da_pair <- testNhoods(milo, design = design_setup$design, 
                           design.df = design_setup$sample_metadata,
                           model.contrasts = contrast_str,
                           fdr.weighting = milo_fdr_weighting)
      
      da_results[[contrast_name]] <- da_pair
    }
  }
  
  return(list(
    milo = milo,
    da = da_results,  # Now a list of pairwise comparisons
    sample_metadata = design_setup$sample_metadata,
    milo_fdr = milo_fdr,
    validation_results = validation_results,
    conditions = conditions
  ))
}

# Note: create_milor_plots() moved to utils/03_differential_plots.R

# Note: create_celltype_proportion_plots() moved to utils/03_differential_plots.R

# REMOVED: Orphaned cell type proportion code (previously lines 600-741)
# This code was outside any function and caused "object 'celltype_field' not found" error
# All cell type proportion plotting functionality is now in utils/03_differential_plots.R

################################################################################
# DIFFERENTIAL STATE (DS) - SEURAT
################################################################################

#' Get cell type field from Seurat metadata
#' 
#' Determines which cell type annotation to use
#' 
#' @param seu Seurat object
#' @return Cell type field name
get_celltype_field <- function(seu) {
  celltype_field <- if ("data_driven_celltype" %in% colnames(seu@meta.data)) {
    "data_driven_celltype"
  } else if ("individual_celltype" %in% colnames(seu@meta.data)) {
    "individual_celltype"
  } else {
    "seurat_clusters"
  }
  return(celltype_field)
}

#' Perform Seurat-based Differential State Analysis
#' 
#' Tests for differential gene expression between conditions within each cluster
#' 
#' @param seu Seurat object
#' @param cond_field Condition field name
#' @param milo Milo object (unused, for compatibility)
#' @param da DA results (unused, for compatibility)
#' @param seurat_params Seurat parameters
#' @param log Logging function
#' @return Data frame with DS results
perform_seurat_ds_analysis <- function(seu, cond_field, milo = NULL, da = NULL, seurat_params, log = cat) {
  
  log("Starting Seurat DE analysis...")
  
  # Extract parameters
  seurat_test_use <- seurat_params$test_use %||% "wilcox"
  seurat_logfc_threshold <- seurat_params$logfc_threshold %||% CFG_LOG2FC_THRESHOLD
  seurat_min_pct <- seurat_params$min_pct %||% CFG_MIN_PCT
  seurat_plot_fdr <- seurat_params$plot_fdr %||% CFG_FDR_THRESHOLD
  seurat_plot_logfc <- seurat_params$plot_logfc %||% CFG_LOG2FC_THRESHOLD
  
  log("Seurat DE parameters:")
  log("  - test:", seurat_test_use)
  log("  - logfc.threshold:", seurat_logfc_threshold)
  log("  - min.pct:", seurat_min_pct)
  
  # Get conditions and cell types
  conditions <- unique(seu@meta.data[[cond_field]])
  celltype_field <- get_celltype_field(seu)
  clusters <- unique(seu@meta.data[[celltype_field]])
  clusters <- clusters[!is.na(clusters)]
  
  log(paste("Analyzing", length(conditions), "conditions across", length(clusters), "cell types"))
  log(paste("Using cell type field:", celltype_field))
  
  # Storage for results
  all_markers <- list()
  
  # Perform pairwise comparisons
  if (length(conditions) >= 2) {
    condition_pairs <- combn(conditions, 2, simplify = FALSE)
    
    for (pair in condition_pairs) {
      cond1 <- pair[1]
      cond2 <- pair[2]
      comp_name <- paste(cond1, "vs", cond2, sep="_")
      
      log(paste("Testing:", comp_name))
      
      for (cluster in clusters) {
        cluster_cells <- colnames(seu)[seu@meta.data[[celltype_field]] == cluster]
        cells_cond1 <- colnames(seu)[seu@meta.data[[cond_field]] == cond1 & seu@meta.data[[celltype_field]] == cluster]
        cells_cond2 <- colnames(seu)[seu@meta.data[[cond_field]] == cond2 & seu@meta.data[[celltype_field]] == cluster]
        
        if (length(cells_cond1) >= CFG_MIN_CELLS_PSEUDOBULK && length(cells_cond2) >= CFG_MIN_CELLS_PSEUDOBULK) {
            markers <- FindMarkers(
              seu,
              ident.1 = cells_cond1,
              ident.2 = cells_cond2,
              test.use = seurat_test_use,
              logfc.threshold = seurat_logfc_threshold,
              min.pct = seurat_min_pct
            )
            
            if (nrow(markers) > 0) {
              markers$gene <- rownames(markers)
              markers$cluster <- cluster
              markers$comparison <- comp_name
              markers$condition1 <- cond1
              markers$condition2 <- cond2
              markers$analysis_type <- "all_cells"
              markers$significant <- markers$p_val_adj < seurat_plot_fdr & abs(markers$avg_log2FC) > seurat_plot_logfc
              all_markers[[paste(cluster, comp_name, sep="_")]] <- markers
              log(paste("  Cluster", cluster, ":", nrow(markers), "DE genes"))
            }
        }
      }
    }
  }
  
  # Combine results
  if (length(all_markers) > 0) {
    ds_tbl <- bind_rows(all_markers)
    log(paste("Total DE genes found:", nrow(ds_tbl)))
    return(ds_tbl)
  } else {
    log("No DE results found")
    return(data.frame())
  }
}

################################################################################
# DIFFERENTIAL STATE (DS) - MUSCAT
################################################################################

#' Setup SingleCellExperiment for Muscat
#' 
#' Creates SCE with proper sample_id, group_id, and cluster_id
#' 
#' @param seu Seurat object
#' @param cond_field Condition field name
#' @param log Logging function
#' @return SingleCellExperiment object
setup_muscat_sce <- function(seu, cond_field, log = cat) {
  
  # Create SCE
  count_matrix <- GetAssayData(seu, assay = "RNA", slot = "counts")
  sce <- SingleCellExperiment(assays = list(counts = count_matrix), colData = seu@meta.data)
  col_data <- colData(sce)
  
  # Debug logging
  log("SEURAT METADATA COLUMNS:")
  log(paste("  Available columns:", paste(colnames(col_data), collapse=", ")))
  
  if ("sample" %in% colnames(col_data)) {
    log(paste("  'sample' column values:", paste(unique(col_data$sample), collapse=", ")))
  }
  if ("orig.ident" %in% colnames(col_data)) {
    log(paste("  'orig.ident' column values:", paste(unique(col_data$orig.ident), collapse=", ")))
  }
  log(paste("  Condition field '", cond_field, "' values:", paste(unique(col_data[[cond_field]]), collapse=", ")))
  
  # Set cluster_id (prefer annotated cell types)
  sce$cluster_id <- if ("singler_main" %in% colnames(col_data) && !all(grepl("^Cluster_", col_data$singler_main))) {
    col_data$singler_main
  } else if ("cell_type" %in% colnames(col_data) && !all(grepl("^Cluster_", col_data$cell_type))) {
    col_data$cell_type
  } else {
    col_data$seurat_clusters
  }
  
  # Set sample_id
  sce$sample_id <- if ("sample" %in% colnames(col_data)) {
    col_data$sample
  } else if ("orig.ident" %in% colnames(col_data)) {
    col_data$orig.ident
  } else {
    paste0(sce[[cond_field]], "_rep", ave(seq_len(ncol(sce)), sce[[cond_field]], FUN = seq_along))
  }
  
  sce$group_id <- sce[[cond_field]]
  
  # Debug: Log what muscat sees
  log("MUSCAT SAMPLE INTERPRETATION:")
  log(paste("  Unique sample_id values:", paste(unique(sce$sample_id), collapse=", ")))
  log(paste("  Unique group_id (condition) values:", paste(unique(sce$group_id), collapse=", ")))
  log(paste("  Total samples:", length(unique(sce$sample_id))))
  log(paste("  Samples per condition:"))
  
  sample_per_cond <- table(sce$group_id, sce$sample_id)
  print(sample_per_cond)
  
  for (cond in unique(sce$group_id)) {
    samples_in_cond <- unique(sce$sample_id[sce$group_id == cond])
    log(paste("    -", cond, ":", length(samples_in_cond), "samples -", paste(samples_in_cond, collapse=", ")))
  }
  
  return(sce)
}

#' Filter clusters for Muscat analysis
#' 
#' Keeps only clusters with sufficient cells across samples
#' 
#' @param sce SingleCellExperiment object
#' @param min_cells_per_sample Minimum cells per sample
#' @param min_samples Minimum number of samples
#' @param log Logging function
#' @return Filtered SCE object or NULL
filter_muscat_clusters <- function(sce, min_cells_per_sample = 10, min_samples = 2, log = cat) {
  
  cluster_sample_counts <- table(sce$cluster_id, sce$sample_id)
  
  valid_clusters <- rownames(cluster_sample_counts)[
    rowSums(cluster_sample_counts >= min_cells_per_sample) >= min_samples
  ]
  
  log(paste("Valid clusters with sufficient cells:", paste(valid_clusters, collapse=", ")))
  
  if (length(valid_clusters) == 0) {
    log("ERROR: No clusters have sufficient cells in at least 2 samples.")
    return(NULL)
  }
  
  # Filter SCE to valid clusters
  sce_filtered <- sce[, sce$cluster_id %in% valid_clusters]
  log(paste("Filtered SCE:", ncol(sce_filtered), "cells in", length(valid_clusters), "clusters"))
  
  return(sce_filtered)
}

#' Setup design and contrast for Muscat
#' 
#' Creates design matrix and contrast for differential testing
#' 
#' @param sce_filtered Filtered SCE object
#' @param log Logging function
#' @return List with design, contrast, and experiment info
setup_muscat_design <- function(sce_filtered, log = cat) {
  
  # Experiment info
  ei <- data.frame(
    sample_id = unique(sce_filtered$sample_id),
    group_id = sce_filtered$group_id[match(unique(sce_filtered$sample_id), sce_filtered$sample_id)]
  )
  rownames(ei) <- ei$sample_id
  
  # Design matrix
  design <- model.matrix(~ 0 + group_id, data = ei)
  colnames(design) <- make.names(gsub("group_id", "", colnames(design)))
  
  condition_levels <- colnames(design)
  
  if (length(condition_levels) <= 1) {
    log("Only one condition - skipping muscat")
    return(NULL)
  }
  
  # Create contrasts
  contrast_pairs <- combn(condition_levels, 2, simplify = FALSE)
  contrast_strings <- sapply(contrast_pairs, function(pair) paste(pair[1], "-", pair[2]))
  contrast <- limma::makeContrasts(contrasts = contrast_strings, levels = design)
  
  log(paste("Created", length(contrast_strings), "contrasts"))
  
  return(list(
    design = design,
    contrast = contrast,
    ei = ei
  ))
}

#' Perform Muscat-based Differential State Analysis
#' 
#' Uses pseudobulk aggregation and edgeR/limma for differential state testing
#' 
#' @param seu Seurat object
#' @param cond_field Condition field name
#' @param muscat_params Muscat parameters
#' @param log Logging function
#' @return List with muscat results and SCE object
perform_muscat_ds_analysis <- function(seu, cond_field, muscat_params, log = cat) {
  
  if (!"muscat" %in% rownames(installed.packages())) {
    log("muscat not installed, skipping pseudobulk DS")
    return(list(res = NULL, sce = NULL, pb = NULL))
  }
  
  log("Starting muscat pseudobulk DS analysis...")
  
  # Extract parameters
  muscat_fdr <- muscat_params$fdr %||% CFG_FDR_THRESHOLD
  muscat_logfc <- muscat_params$logfc %||% 0
  
  # Setup SCE
  sce <- setup_muscat_sce(seu, cond_field, log)
  
  # Filter clusters
  sce_filtered <- filter_muscat_clusters(
    sce = sce,
    min_cells_per_sample = CFG_MIN_CELLS_PSEUDOBULK,
    min_samples = 2,
    log = log
  )
  
  if (is.null(sce_filtered)) {
    return(list(res = NULL, sce = NULL, pb = NULL, muscat_fdr = muscat_fdr, muscat_logfc = muscat_logfc))
  }
  
  # Pseudobulk aggregation
  pb <- aggregateData(sce_filtered, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
  log("Using sum aggregation (recommended by Junttila et al. benchmark)")
  
  # Setup design and contrast
  design_setup <- setup_muscat_design(sce_filtered, log)
  
  if (is.null(design_setup)) {
    log("Only one condition - skipping muscat")
    res <- NULL
  } else {
    # Add experiment info to pb
    metadata(pb)$experiment_info <- design_setup$ei
    
    # Run differential state analysis
    log("Running pbDS with edgeR method (recommended for high sensitivity)")
    res <- pbDS(pb, design = design_setup$design, contrast = design_setup$contrast, method = "edgeR")
    
    # Check if all results are empty
    has_results <- if (!is.null(res) && !is.null(res$table)) {
      log(paste("muscat returned", length(names(res$table)), "comparisons"))
      any(sapply(names(res$table), function(comparison) {
        log(paste("  Comparison:", comparison, "has", length(names(res$table[[comparison]])), "clusters"))
        any(sapply(names(res$table[[comparison]]), function(cluster) {
          tbl <- res$table[[comparison]][[cluster]]
          if (!is.null(tbl) && is.data.frame(tbl) && nrow(tbl) > 0) {
            log(paste("    Cluster", cluster, ":", nrow(tbl), "genes"))
            if ("FDR" %in% colnames(tbl)) {
              n_sig <- sum(tbl$FDR < muscat_fdr, na.rm = TRUE)
              log(paste("      Significant genes (FDR <", muscat_fdr, "):", n_sig))
              log(paste("      FDR range:", min(tbl$FDR, na.rm=TRUE), "-", max(tbl$FDR, na.rm=TRUE)))
            }
            return(TRUE)
          }
          FALSE
        }))
      }))
    } else FALSE
    
    # If no significant results with current FDR threshold, try more lenient threshold
    if (!has_results && muscat_fdr == 0.05 && !is.null(res) && !is.null(res$table)) {
      muscat_fdr <- 0.1
      has_results <- any(sapply(names(res$table), function(comparison) {
        any(sapply(names(res$table[[comparison]]), function(cluster) {
          tbl <- res$table[[comparison]][[cluster]]
          if (!is.null(tbl) && is.data.frame(tbl) && nrow(tbl) > 0 && "FDR" %in% colnames(tbl)) {
            n_sig_new <- sum(tbl$FDR < muscat_fdr, na.rm = TRUE)
            if (n_sig_new > 0) log(paste("    Cluster", cluster, "with FDR < 0.1:", n_sig_new, "significant genes"))
            return(n_sig_new > 0)
          }
          FALSE
        }))
      }))
    }
  }
  
  return(list(
    res = res,
    sce = sce,
    pb = pb,
    muscat_fdr = muscat_fdr,
    muscat_logfc = muscat_logfc
  ))
}

################################################################################
# VISUALIZATION - MUSCAT PLOTS
################################################################################
# Note: create_muscat_plots() has been moved to utils/03_differential_plots.R

# Perform combined Seurat->Muscat validation analysis
# Strategy: Use Seurat for high sensitivity candidate discovery, then validate with Muscat
# Returns: Data frame with validated genes and validation metrics
perform_combined_ds_validation <- function(seurat_result, muscat_result, plots_dir, dash_dir, data_dir, log = cat) {
  
  log("=== Combined DS Analysis: Seurat Discovery -> Muscat Validation ===")
  log("Strategy: Seurat (Wilcoxon) for sensitive candidate discovery")
  log("          Muscat (Pseudobulk) for robust validation across samples")
  
  # Check if both methods have results
  if (is.null(seurat_result) || nrow(seurat_result) == 0) {
    log("No Seurat results for validation")
    return(NULL)
  }
  
  if (is.null(muscat_result$res) || is.null(muscat_result$res$table)) {
    log("No muscat results for validation")
    return(NULL)
  }
  
  # Extract Seurat candidates (sensitive discovery)
  seurat_candidates <- seurat_result[seurat_result$significant == TRUE, ]
  log(paste("Seurat identified", nrow(seurat_candidates), "candidate genes (high sensitivity)"))
  
  # Build muscat lookup table
  muscat_validated <- data.frame()
  
  for (comparison in names(muscat_result$res$table)) {
    for (cluster in names(muscat_result$res$table[[comparison]])) {
      tbl <- muscat_result$res$table[[comparison]][[cluster]]
      if (!is.null(tbl) && is.data.frame(tbl) && nrow(tbl) > 0) {
        if ("FDR" %in% colnames(tbl) && "logFC" %in% colnames(tbl)) {
          tbl$gene <- rownames(tbl)
          tbl$cluster <- cluster
          tbl$comparison <- comparison
          tbl$muscat_significant <- tbl$FDR < muscat_result$muscat_fdr & 
                                   abs(tbl$logFC) > muscat_result$muscat_logfc
          muscat_validated <- rbind(muscat_validated, tbl)
        }
      }
    }
  }
  
  if (nrow(muscat_validated) == 0) {
    log("No muscat data available for validation")
    return(NULL)
  }
  
  log(paste("muscat tested", length(unique(muscat_validated$gene)), "genes across clusters"))
  
  # Merge Seurat candidates with Muscat validation
  validated_genes <- merge(
    seurat_candidates,
    muscat_validated[, c("gene", "cluster", "comparison", "logFC", "FDR", "muscat_significant")],
    by = c("gene", "cluster", "comparison"),
    suffixes = c("_seurat", "_muscat")
  )
  
  if (nrow(validated_genes) == 0) {
    log("No overlapping genes between Seurat candidates and muscat results")
    return(NULL)
  }
  
  # Classification of genes
  validated_genes$validation_status <- ifelse(
    validated_genes$muscat_significant,
    "Validated",
    "Seurat-only"
  )
  
  # Calculate concordance (same direction)
  validated_genes$direction_concordant <- sign(validated_genes$avg_log2FC) == sign(validated_genes$logFC_muscat)
  
  # Validated = muscat significant AND concordant direction
  validated_genes$final_validated <- validated_genes$muscat_significant & validated_genes$direction_concordant
  
  # Summary statistics
  n_candidates <- nrow(seurat_candidates)
  n_tested <- nrow(validated_genes)
  n_validated <- sum(validated_genes$final_validated, na.rm = TRUE)
  n_discordant <- sum(!validated_genes$direction_concordant, na.rm = TRUE)
  validation_rate <- round(100 * n_validated / n_tested, 1)
  
  log(paste("\nValidation Results:"))
  log(paste("  Seurat candidates:", n_candidates))
  log(paste("  Tested by muscat:", n_tested))
  log(paste("  Validated (significant + concordant):", n_validated, paste0("(", validation_rate, "% of tested)")))
  log(paste("  Seurat-only (not validated):", n_tested - n_validated))
  log(paste("  Direction discordant:", n_discordant))
  
  # Save results
  write.csv(validated_genes, file.path(data_dir, "DS_combined_validated_genes.csv"), row.names = FALSE)
  
  validated_only <- validated_genes[validated_genes$final_validated == TRUE, ]
  write.csv(validated_only, file.path(data_dir, "DS_combined_validated_only.csv"), row.names = FALSE)
  
  # Create summary table
  validation_summary <- data.frame(
    Category = c("Seurat Candidates (Discovery)", "Tested by Muscat", "Validated (Robust)", 
                "Seurat-only (Not Validated)", "Direction Discordant"),
    Count = c(n_candidates, n_tested, n_validated, n_tested - n_validated, n_discordant),
    Percentage = c(100, 
                  round(100 * n_tested / n_candidates, 1),
                  validation_rate,
                  round(100 * (n_tested - n_validated) / n_tested, 1),
                  round(100 * n_discordant / n_tested, 1))
  )
  write.csv(validation_summary, file.path(data_dir, "DS_validation_summary.csv"), row.names = FALSE)
  
  # Visualization: Validation funnel
  p_funnel <- create_validation_funnel_plot(
    validation_summary = validation_summary,
    n_candidates = n_candidates,
    n_validated = n_validated,
    validation_rate = validation_rate,
    plots_dir = plots_dir,
    base_size = CFG_BASE_SIZE,
    title_size = CFG_TITLE_SIZE,
    axis_title_size = CFG_AXIS_TITLE_SIZE
  )
  
  # Scatter plot: Seurat vs Muscat log fold changes for tested genes
  p_scatter <- create_validation_scatter_plot(
    validated_genes = validated_genes,
    n_validated = n_validated,
    n_tested = n_tested,
    plots_dir = plots_dir,
    dash_dir = dash_dir,
    base_size = CFG_BASE_SIZE,
    title_size = CFG_TITLE_SIZE,
    log = log
  )
  
  # Per-cluster validation rates
  p_cluster <- create_validation_by_cluster_plot(
    validated_genes = validated_genes,
    plots_dir = plots_dir,
    data_dir = data_dir,
    base_size = CFG_BASE_SIZE,
    title_size = CFG_TITLE_SIZE,
    axis_title_size = CFG_AXIS_TITLE_SIZE
  )
  
  log("=== Combined DS Validation Completed ===")
  log(paste("Recommendation: Use", n_validated, "validated genes for downstream analysis"))
  log("These genes show robust effects across biological replicates (samples)")
  
  return(list(
    validated_genes = validated_genes,
    validated_only = validated_only,
    summary = validation_summary,
    n_candidates = n_candidates,
    n_validated = n_validated,
    validation_rate = validation_rate
  ))
}

# Wrapper function to run complete differential analysis pipeline
# Executes DA analysis (MiloR) and DS analysis (Seurat + Muscat)
# Returns: List with all analysis results
run_differential_analysis <- function(seu, cond_field, milo_params, seurat_params, muscat_params,
                                     plots_dir, dash_dir, data_dir, analysis_params = NULL, log = cat) {
  
  log("=== Starting Complete Differential Analysis Pipeline ===")
  
  # Create output directories if they don't exist
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(dash_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 0. Cell Type Proportion Analysis
  log("Creating cell type proportion plots...")
  create_celltype_proportion_plots(seu, cond_field, celltype_field = "data_driven_celltype",
                                  plots_dir, data_dir, log)
  
  # 1. DA Analysis - MiloR
  milor_result <- perform_milor_analysis(seu, cond_field, milo_params, analysis_params, log)
  
  if (!is.null(milor_result$da)) {
    if (is.list(milor_result$da)) {
      da_combined <- do.call(rbind, lapply(names(milor_result$da), function(nm) {
        df <- as.data.frame(milor_result$da[[nm]])
        df$contrast <- nm
        df
      }))
      write.csv(da_combined, file.path(data_dir, "DA_milo_results.csv"), row.names = FALSE)
    } else {
      write.csv(as.data.frame(milor_result$da), file.path(data_dir, "DA_milo_results.csv"), row.names = FALSE)
    }
    
    # Save validation results if available
    if (!is.null(milor_result$validation_results)) {
      saveRDS(milor_result$validation_results, 
              file.path(data_dir, "milor_validation_results.rds"))
    }
    
    create_milor_plots(milor_result$milo, milor_result$da, seu, cond_field, 
                      milor_result$sample_metadata, milor_result$milo_fdr, 
                      plots_dir, dash_dir, 
                      validation_results = milor_result$validation_results,
                      log)
  }
  
  # 2.1 DS Analysis - Seurat
  seurat_result <- perform_seurat_ds_analysis(seu, cond_field, milor_result$milo, 
                                              milor_result$da, seurat_params, log)
  
  if (nrow(seurat_result) > 0) {
    write.csv(seurat_result, file.path(data_dir, "DS_seurat_results.csv"), row.names = FALSE)
    
    # Create volcano plots for Seurat DE results
      if ("avg_log2FC" %in% colnames(seurat_result) && "p_val_adj" %in% colnames(seurat_result)) {
        # Overall volcano plot
        seurat_result$significant <- seurat_result$p_val_adj < seurat_params$plot_fdr & 
                                     abs(seurat_result$avg_log2FC) > seurat_params$plot_logfc
        
        # Mark top genes
        seurat_result$label <- ""
        seurat_result$regulation <- "Not Significant"
        sig_seurat <- seurat_result[seurat_result$significant, ]
        if (nrow(sig_seurat) > 0) {
          sig_seurat$score <- abs(sig_seurat$avg_log2FC) * -log10(sig_seurat$p_val_adj)
          sig_seurat <- sig_seurat[order(-sig_seurat$score), ]
          top_seurat <- head(sig_seurat, 12)
          seurat_result$label <- ifelse(seurat_result$gene %in% top_seurat$gene, seurat_result$gene, "")
          seurat_result$regulation <- ifelse(!seurat_result$significant, "Not Significant",
                                            ifelse(seurat_result$avg_log2FC > 0, "Up-regulated", "Down-regulated"))
        }
        
        if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel", repos = "https://cloud.r-project.org/")
        
        # Create main volcano plot using plot function
        p_volcano_seurat <- create_seurat_volcano_plot(
          seurat_result = seurat_result,
          plot_fdr = seurat_params$plot_fdr,
          plot_logfc = seurat_params$plot_logfc,
          plots_dir = plots_dir,
          title = "Seurat DS (2.1): Volcano Plot (All Comparisons)",
          filename = "DS_seurat_volcano_all.png",
          base_size = CFG_BASE_SIZE,
          title_size = CFG_TITLE_SIZE,
          legend_size = CFG_LEGEND_SIZE
        )
        
        # Volcano plot per comparison if multiple comparisons exist
        if ("comparison" %in% colnames(seurat_result)) {
          comparison_plots <- create_comparison_volcano_per_condition(
            seurat_result = seurat_result,
            plot_fdr = seurat_params$plot_fdr,
            plot_logfc = seurat_params$plot_logfc,
            plots_dir = plots_dir,
            base_size = CFG_BASE_SIZE
          )
        }
      }
  }
  
  # 3. DS Analysis - muscat
  muscat_result <- perform_muscat_ds_analysis(seu, cond_field, muscat_params, log)
  
  if (!is.null(muscat_result$res) && !is.null(muscat_result$res$table)) {
    write.csv(muscat_result$res_tbl, file.path(data_dir, "DS_muscat_results.csv"), row.names = FALSE)
    create_muscat_plots(muscat_result, seu, plots_dir, dash_dir, log)
    
    # Create muscat volcano plots
    muscat_fdr <- muscat_result$muscat_fdr
    muscat_logfc <- muscat_result$muscat_logfc
    muscat_volcano_plots <- list()
    all_muscat_data <- list()
    
    for (comparison in names(muscat_result$res$table)) {
      for (cluster in names(muscat_result$res$table[[comparison]])) {
        tbl <- muscat_result$res$table[[comparison]][[cluster]]
        if (!is.null(tbl) && is.data.frame(tbl) && nrow(tbl) > 0 && 
            "logFC" %in% colnames(tbl) && "FDR" %in% colnames(tbl)) {
          tbl$gene <- rownames(tbl)
          tbl$comparison <- comparison
          tbl$cluster <- cluster
          all_muscat_data[[paste(comparison, cluster, sep = "_")]] <- tbl
          
          # Create individual volcano plot using plot function
          p_volcano_muscat <- create_muscat_volcano_plot(
            tbl = tbl,
            muscat_fdr = muscat_fdr,
            muscat_logfc = muscat_logfc,
            cluster = cluster,
            comparison = comparison,
            plots_dir = plots_dir,
            base_size = CFG_BASE_SIZE,
            title_size = CFG_TITLE_SIZE,
            axis_title_size = CFG_AXIS_TITLE_SIZE,
            legend_size = CFG_LEGEND_SIZE
          )
          
          # Store for combined figure
          muscat_volcano_plots[[paste(cluster, comparison, sep = "_")]] <- p_volcano_muscat
        }
      }
    }
    
    # Create combined muscat volcano figure
    if (length(muscat_volcano_plots) > 0) {
      p_combined_muscat <- create_combined_plot_grid(
        plot_list = muscat_volcano_plots,
        ncols = 4,
        title = "muscat Differential State: All Cell Types",
        plots_dir = plots_dir,
        filename = "DS_muscat_volcano_all_combined.png",
        width_per_plot = 8,
        height_per_plot = 7,
        title_size = CFG_TITLE_SIZE
      )
      
      log(paste("Created combined muscat volcano plot with", length(muscat_volcano_plots), "clusters"))
    }
    
    # Overall volcano plot
    if (length(all_muscat_data) > 0) {
      combined_data <- dplyr::bind_rows(all_muscat_data)
      if (nrow(combined_data) > 0 && "logFC" %in% colnames(combined_data) && "FDR" %in% colnames(combined_data)) {
        p_volcano_all <- create_muscat_overall_volcano(
          combined_data = combined_data,
          muscat_fdr = muscat_fdr,
          muscat_logfc = muscat_logfc,
          plots_dir = plots_dir,
          base_size = CFG_BASE_SIZE
        )
        log("Created muscat volcano plots")
      }
    }
  }
  
  log("=== Differential Analysis Pipeline Completed ===")
  
  # 4. Combined validation: Seurat discovery -> Muscat validation
  validation_result <- perform_combined_ds_validation(
    seurat_result, muscat_result, plots_dir, dash_dir, data_dir, log
  )
  
  # 5. Compare Seurat and muscat results (traditional comparison)
  compare_seurat_muscat(seurat_result, muscat_result, plots_dir, dash_dir, data_dir, log)
  
  # Save analysis summary
  log("Saving analysis summary...")
  analysis_summary <- list(
    timestamp = Sys.time(),
    milo_performed = !is.null(milor_result$da),
    seurat_performed = !is.null(seurat_result),
    muscat_performed = !is.null(muscat_result$res),
    validation_performed = !is.null(validation_result)
  )
  saveRDS(analysis_summary, file.path(data_dir, "analysis_summary.rds"))
  log("=== All Analyses Completed Successfully ===")
  
  # End performance tracking
  end_performance_tracking(
    perf_tracker,
    success = TRUE,
    additional_metrics = list(
      n_cells = ncol(seu),
      n_genes = nrow(seu),
      n_comparisons = nrow(analysis_summary$comparisons)
    )
  )
  
  return(list(
    milor = milor_result,
    seurat = seurat_result,
    muscat = muscat_result,
    validated = validation_result
  ))
}

# Compare Seurat and muscat differential expression results
compare_seurat_muscat <- function(seurat_result, muscat_result, plots_dir, dash_dir, data_dir, log = cat) {
  
  log("=== Comparing Seurat and muscat Results ===")
  
  # Check if both have results
  if (is.null(seurat_result) || nrow(seurat_result) == 0) {
    log("No Seurat results to compare")
    return()
  }
  
  if (is.null(muscat_result$res) || is.null(muscat_result$res$table)) {
    log("No muscat results to compare")
    return()
  }
  
    # Extract significant genes from Seurat
    seurat_sig <- seurat_result[seurat_result$significant == TRUE, ]
    seurat_genes <- unique(seurat_sig$gene)
    log(paste("Seurat: Found", length(seurat_genes), "significant genes"))
    
    # Extract significant genes from muscat
    # Track genes at FDR 0.05 separately to show fallback in plot
    muscat_genes <- character(0)
    muscat_genes_at_fdr_005 <- character(0)  # Track original 0.05 threshold
    muscat_data <- list()
    
    for (comparison in names(muscat_result$res$table)) {
      for (cluster in names(muscat_result$res$table[[comparison]])) {
        tbl <- muscat_result$res$table[[comparison]][[cluster]]
        if (!is.null(tbl) && is.data.frame(tbl) && nrow(tbl) > 0) {
          if ("FDR" %in% colnames(tbl) && "logFC" %in% colnames(tbl)) {
            # Count genes at original FDR 0.05 (before any fallback)
            sig_genes_005 <- rownames(tbl)[tbl$FDR < 0.05 & 
                                           abs(tbl$logFC) > muscat_result$muscat_logfc]
            muscat_genes_at_fdr_005 <- c(muscat_genes_at_fdr_005, sig_genes_005)
            
            # Count genes at final FDR (may be 0.05 or 0.1 after fallback)
            sig_genes <- rownames(tbl)[tbl$FDR < muscat_result$muscat_fdr & 
                                       abs(tbl$logFC) > muscat_result$muscat_logfc]
            muscat_genes <- c(muscat_genes, sig_genes)
            
            # Store for later comparison
            tbl$gene <- rownames(tbl)
            tbl$cluster <- cluster
            tbl$comparison <- comparison
            muscat_data[[paste(cluster, comparison, sep = "_")]] <- tbl
          }
        }
      }
    }
    
    muscat_genes <- unique(muscat_genes)
    muscat_genes_at_fdr_005 <- unique(muscat_genes_at_fdr_005)
    log(paste("muscat: Found", length(muscat_genes), "significant genes"))
    
    # Calculate overlaps
    overlap_genes <- intersect(seurat_genes, muscat_genes)
    seurat_only <- setdiff(seurat_genes, muscat_genes)
    muscat_only <- setdiff(muscat_genes, seurat_genes)
    
    log(paste("Overlap:", length(overlap_genes), "genes"))
    log(paste("Seurat only:", length(seurat_only), "genes"))
    log(paste("muscat only:", length(muscat_only), "genes"))
    
    # Create comparison summary
    comparison_summary <- data.frame(
      Method = c("Seurat", "muscat", "Overlap", "Seurat Only", "muscat Only"),
      Count = c(length(seurat_genes), length(muscat_genes), 
               length(overlap_genes), length(seurat_only), length(muscat_only)),
      Percentage = c(100, 100, 
                    round(100 * length(overlap_genes) / length(seurat_genes), 2),
                    round(100 * length(seurat_only) / length(seurat_genes), 2),
                    round(100 * length(muscat_only) / length(muscat_genes), 2))
    )
    
    write.csv(comparison_summary, file.path(data_dir, "method_comparison_summary.csv"), row.names = FALSE)
    
    # UpSet plot
    create_method_comparison_upset(
      seurat_genes = seurat_genes,
      muscat_genes = muscat_genes,
      plots_dir = plots_dir,
      log = log
    )
    
    # Barplot comparison
    p_barplot <- create_method_comparison_barplot(
      comparison_summary = comparison_summary,
      seurat_genes = seurat_genes,
      muscat_genes = muscat_genes,
      overlap_genes = overlap_genes,
      seurat_only = seurat_only,
      muscat_only = muscat_only,
      plots_dir = plots_dir,
      muscat_genes_at_fdr_threshold = length(muscat_genes_at_fdr_005),
      muscat_fdr = muscat_result$muscat_fdr,
      base_size = CFG_BASE_SIZE,
      title_size = CFG_TITLE_SIZE,
      axis_title_size = CFG_AXIS_TITLE_SIZE
    )
    
    # Scatter plot: Compare log fold changes for overlapping genes
    p_scatter <- create_method_comparison_logfc_scatter(
      seurat_sig = seurat_sig,
      overlap_genes = overlap_genes,
      muscat_data = muscat_data,
      plots_dir = plots_dir,
      data_dir = data_dir,
      base_size = CFG_BASE_SIZE,
      title_size = CFG_TITLE_SIZE,
      log = log
    )
    
    log("=== Method Comparison Completed ===")
}

#===============================================================================
# MAIN
#===============================================================================

# Parse command-line arguments (same pattern as 01 and 02)
option_list <- list(
  make_option("--seurat", type="character", default="sce_clustered.rds", help="Clustered Seurat object"),
  make_option("--out", type="character", default="output/", help="Output path"),
  make_option("--config", type="character", default=NULL, help="Path to nextflow.config file"),
  make_option("--milo_fdr", type="numeric", default=NULL),
  make_option("--milo_fdr_weighting", type="character", default=NULL),
  make_option("--seurat_test_use", type="character", default=NULL),
  make_option("--seurat_logfc_threshold", type="numeric", default=NULL),
  make_option("--seurat_min_pct", type="numeric", default=NULL),
  make_option("--seurat_plot_fdr", type="numeric", default=NULL),
  make_option("--seurat_plot_logfc", type="numeric", default=NULL),
  make_option("--muscat_fdr", type="numeric", default=NULL),
  make_option("--muscat_logfc", type="numeric", default=NULL)
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Parse nextflow configuration
nf_params <- parse_nextflow_config(opt$config)

# Load configuration from nextflow.config with fallback values
# Command-line arguments override config file values
CFG_RANDOM_SEED <- as.integer(nf_params$differential.random_seed %||% 42)
CFG_BASE_SIZE <- as.integer(nf_params$differential.base_size %||% 18)
CFG_TITLE_SIZE <- as.integer(nf_params$differential.title_size %||% 24)
CFG_AXIS_TITLE_SIZE <- as.integer(nf_params$differential.axis_title_size %||% 20)
CFG_LEGEND_SIZE <- as.integer(nf_params$differential.legend_size %||% 16)
CFG_MIN_CELLS_PSEUDOBULK <- as.integer(nf_params$differential.min_cells_pseudobulk %||% 10)
CFG_MIN_CELLS_CLUSTER <- as.integer(nf_params$differential.min_cells_cluster %||% 3)
CFG_LOG2FC_THRESHOLD <- as.numeric(nf_params$differential.log2fc_threshold %||% 0.25)
CFG_MIN_PCT <- as.numeric(nf_params$differential.min_pct %||% 0.10)
CFG_FDR_THRESHOLD <- as.numeric(nf_params$differential.fdr_threshold %||% 0.05)
CFG_TOP_GENES_HEATMAP <- as.integer(nf_params$differential.top_genes_heatmap %||% 20)

# Method-specific parameters with command-line overrides
milo_params <- list(
  fdr = opt$milo_fdr %||% nf_params$differential.milo_fdr %||% CFG_FDR_THRESHOLD,
  fdr_weighting = opt$milo_fdr_weighting %||% nf_params$differential.milo_fdr_weighting %||% "none"
)

seurat_params <- list(
  test_use = opt$seurat_test_use %||% nf_params$differential.seurat_test_use %||% "wilcox",
  logfc_threshold = opt$seurat_logfc_threshold %||% nf_params$differential.seurat_logfc_threshold %||% CFG_LOG2FC_THRESHOLD,
  min_pct = opt$seurat_min_pct %||% nf_params$differential.seurat_min_pct %||% CFG_MIN_PCT,
  plot_fdr = opt$seurat_plot_fdr %||% nf_params$differential.seurat_plot_fdr %||% CFG_FDR_THRESHOLD,
  plot_logfc = opt$seurat_plot_logfc %||% nf_params$differential.seurat_plot_logfc %||% CFG_LOG2FC_THRESHOLD
)

muscat_params <- list(
  fdr = opt$muscat_fdr %||% nf_params$differential.muscat_fdr %||% CFG_FDR_THRESHOLD,
  logfc = opt$muscat_logfc %||% nf_params$differential.muscat_logfc %||% 0
)

# Set global random seed
set.seed(CFG_RANDOM_SEED)

# Setup output directories
output_dir <- opt$out
plots_dir <- file.path(output_dir, "plots")
dash_dir <- file.path(output_dir, "interactive")
data_dir <- file.path(output_dir, "data")
for (d in c(output_dir, plots_dir, dash_dir, data_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Start performance tracking
perf_tracker <- start_performance_tracking("04_differential_analysis", data_dir)

# Load Seurat object
cat("Loading clustered Seurat object...\n")
seu <- readRDS(opt$seurat)

# Load analysis parameters
analysis_params <- NULL
seurat_real_path <- Sys.readlink(opt$seurat)
if (seurat_real_path == "") seurat_real_path <- opt$seurat

param_paths <- c(
  file.path(dirname(seurat_real_path), "data", "analysis_parameters.rds"),
  file.path(dirname(opt$seurat), "data", "analysis_parameters.rds")
)

for (param_path in param_paths) {
  if (file.exists(param_path)) {
      analysis_params <- readRDS(param_path)
      cat("Loaded analysis parameters from:", param_path, "\n")
      break
  }
}

# Run the analysis
run_differential_analysis(
  seu = seu,
  cond_field = "condition",
  milo_params = milo_params,
  seurat_params = seurat_params,
  muscat_params = muscat_params,
  plots_dir = plots_dir,
  dash_dir = dash_dir,
  data_dir = data_dir,
  analysis_params = analysis_params,
  log = cat
)

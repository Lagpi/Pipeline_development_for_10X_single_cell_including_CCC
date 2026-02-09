#!/usr/bin/env Rscript
#===============================================================================
# 02_clustering_annotation.R
#
# Single-cell RNA-seq Clustering and Cell Type Annotation Pipeline
# 
# Main Processing Steps:
#   1. Cell cycle scoring (Phase, S.Score, G2M.Score)
#   2. PCA dimensionality reduction with JackStraw validation
#   3. UMAP embedding for 2D visualization
#   4. Multi-resolution graph-based clustering with quality metrics
#   5. Automated cell type annotation using SingleR
#   6. Marker gene identification for each cluster
#   7. Comprehensive visualization and quality assessment
#
# Functions defined:
#   - perform_cell_cycle_scoring(): Cell cycle phase annotation
#   - run_pca(): PCA dimensionality reduction
#   - perform_jackstraw(): Statistical PC selection
#   - perform_pca_analysis(): Complete PCA with JackStraw
#   - create_umap(): UMAP embedding generation
#   - perform_multires_clustering(): Multi-resolution clustering
#   - calculate_clustering_quality(): Quality metrics (silhouette, WCSS/BCSS)
#   - select_optimal_resolution(): Optimal resolution selection
#   - perform_singler_annotation(): SingleR cell type annotation
#   - find_cluster_markers(): Differential marker gene identification
#
# Output files:
#   - clustered_annotated_seurat.rds: Final Seurat object with clusters and annotations
#   - cluster_markers.csv: Marker genes for each cluster
#   - cluster_to_celltype_mapping.csv: Cluster to cell type assignments
#   - clustering_quality_metrics.csv: Quality metrics per resolution
#   - Various visualization plots (PNG + HTML interactive)
#===============================================================================

#===============================================================================
# 0. SETUP & INITIALIZATION
#===============================================================================

#-------------------------------------------------------------------------------
# 0.1. Load Required Libraries
#-------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)  # Command-line argument parsing
  library(Seurat)    # Single-cell RNA-seq analysis
  library(ggplot2)   # Data visualization
  library(patchwork) # Combine multiple plots
  library(dplyr)     # Data manipulation
  library(tidyr)     # Data tidying
  library(plotly)    # Interactive plots
  library(htmlwidgets)  # Save interactive widgets
  library(SingleR)   # Cell type annotation
  library(celldex)   # Reference datasets for SingleR
  library(SingleCellExperiment)  # Data structure for single-cell data
  library(future)    # Parallel processing
  library(future.apply)  # Apply functions in parallel
})

#-------------------------------------------------------------------------------
# 0.2. Load Utility Functions
#-------------------------------------------------------------------------------
# Source bootstrap utility to get script directory
source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))), "utils", "bootstrap.R"))
script_dir <- get_script_dir()
source(file.path(script_dir, "utils", "common.R"))
source(file.path(script_dir, "utils", "logging.R"))
source(file.path(script_dir, "utils", "plotting.R"))
source(file.path(script_dir, "utils", "colors.R"))
source(file.path(script_dir, "utils", "performance.R"))
source(file.path(script_dir, "utils", "02_clustering_plots.R"))

# Define helper operator for NULL coalescing
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (is.character(a) && a == "")) b else a

#===============================================================================
# 1. MODULAR ANALYSIS FUNCTIONS
#===============================================================================

#' Perform cell cycle scoring on Seurat object
#' 
#' Assigns cell cycle phase (G1, S, G2M) to each cell based on expression
#' of canonical cell cycle genes. Converts human gene names to mouse if needed.
#' 
#' @param seu Seurat object with normalized data
#' @param species Species identifier ("mouse" or "human")
#' @return Seurat object with Phase, S.Score, G2M.Score metadata
perform_cell_cycle_scoring <- function(seu, species = "mouse") {
  log("Performing cell cycle scoring...")
  
  # Load cell cycle genes from Seurat package
  # Load built-in cell cycle genes
  cc.genes <- NULL
  tryCatch({
    data("cc.genes", package = "Seurat", envir = environment())
  }, error = function(e) {
    log(paste("Warning: Could not load cc.genes:", e$message))
  })
  
  if (is.null(cc.genes)) {
    log("ERROR: cc.genes not available from Seurat package")
    return(seu)
  }
  
  # Get available features in the object
  available_features <- rownames(seu)
  
  # Convert to appropriate species
  if (tolower(species) == "mouse") {
    # Convert human cell cycle genes to mouse (first letter uppercase, rest lowercase)
    s.genes <- tolower(cc.genes$s.genes)
    s.genes <- paste0(toupper(substring(s.genes, 1, 1)), substring(s.genes, 2))
    g2m.genes <- tolower(cc.genes$g2m.genes)
    g2m.genes <- paste0(toupper(substring(g2m.genes, 1, 1)), substring(g2m.genes, 2))
    
    # Filter to only genes present in the dataset
    s.genes <- s.genes[s.genes %in% available_features]
    g2m.genes <- g2m.genes[g2m.genes %in% available_features]
    
    log(sprintf("Converted and filtered to %d S-phase and %d G2M-phase genes available in mouse dataset", 
                length(s.genes), length(g2m.genes)))
  } else {
    # Use default human genes, but filter to available
    s.genes <- cc.genes$s.genes[cc.genes$s.genes %in% available_features]
    g2m.genes <- cc.genes$g2m.genes[cc.genes$g2m.genes %in% available_features]
    log(sprintf("Using %d S-phase and %d G2M-phase genes for human", 
                length(s.genes), length(g2m.genes)))
  }
  
  # Perform cell cycle scoring (suppress warnings about missing genes)
  # NOTE: Cell cycle effects are NOT regressed out in ScaleData, as recommended
  # by Heumos et al. (2023) - cell cycle information is important for trajectory inference
  seu <- suppressWarnings(
    CellCycleScoring(
      seu, 
      s.features = s.genes, 
      g2m.features = g2m.genes,
      set.ident = FALSE
    )
  )
  
  # Log results
  if ("Phase" %in% colnames(seu@meta.data)) {
    phase_counts <- table(seu$Phase)
    log("Cell cycle phases identified:")
    for (phase in names(phase_counts)) {
      log(sprintf("  %s: %d cells (%.1f%%)", phase, phase_counts[phase], 
                  100 * phase_counts[phase] / ncol(seu)))
    }
    log("Cell cycle scoring completed successfully")
  } else {
    log("Warning: Phase column not created by CellCycleScoring")
    seu$Phase <- "Unknown"
    seu$S.Score <- 0
    seu$G2M.Score <- 0
  }
  
  return(seu)
}

#' Run PCA dimensionality reduction
#' 
#' Performs Principal Component Analysis on scaled data to reduce
#' dimensionality while preserving variance. Cell cycle effects are NOT
#' regressed out (important for trajectory inference, Heumos et al. 2023).
#' 
#' @param seu Seurat object with scaled data
#' @param npcs Number of principal components to compute
#' @return Seurat object with PCA reduction
run_pca <- function(seu, npcs = CFG_PCA_SEARCH_SPACE) {
  log(sprintf("Computing PCA with %d components...", npcs))
  # ScaleData without vars.to.regress - preserving cell cycle signals for trajectory
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = npcs, verbose = FALSE)
  log("PCA completed")
  return(seu)
}

#' Perform JackStraw analysis for PC selection
#' 
#' Uses permutation testing to identify statistically significant PCs.
#' Randomly permutes gene-cell associations to create null distribution,
#' then compares observed PC scores to null.
#' 
#' @param seu Seurat object with PCA computed
#' @param dims Number of PCs to test (max PC number)
#' @param num_replicate Number of permutation replicates
#' @return List with:
#'   - seu: Seurat object with JackStraw results
#'   - significant_pcs: Vector of significant PC indices
#'   - success: Whether JackStraw succeeded
perform_jackstraw <- function(seu, dims = CFG_PCA_SEARCH_SPACE, num_replicate = CFG_JACKSTRAW_REPLICATES) {
  significant_pcs <- NULL
  success <- FALSE
  
  # Ensure dims is a scalar (max PC number to test)
  # If dims is already a vector like 1:30, take the max value
  if (length(dims) > 1) {
    dims <- max(dims)
  }
  
  # Convert to integer to be safe
  dims <- as.integer(dims)
  
  # Validate dims is reasonable
  if (is.na(dims) || dims < 1 || dims > 100) {
    log(sprintf("Invalid dims parameter (%s), using fallback", dims))
    return(list(seu = seu, significant_pcs = NULL, success = FALSE))
  }
  
  # Store old limit before try-catch
  old_limit <- getOption("future.globals.maxSize")
  
  tryCatch({
    log(sprintf("Running JackStraw to identify significant PCs (testing %d:%d)...", CFG_PCA_MIN_DIMS, dims))
    log(sprintf("Using %d replicates with %d cores for parallelization", num_replicate, CFG_JACKSTRAW_CORES))
    
    # Increase memory limit for large Seurat objects (each worker needs a copy)
    # For 78K cells: ~10GB per worker × 4 workers = 40GB buffer
    options(future.globals.maxSize = 50000 * 1024^2)  # 50 GB
    
    # Enable parallel processing with multisession (works better on SLURM than multicore)
    plan(multisession, workers = CFG_JACKSTRAW_CORES)
    
    # Call JackStraw - test full range 1:dims but will filter to min_dims:dims later
    seu <- JackStraw(seu, dims = dims, num.replicate = num_replicate, verbose = TRUE)
    seu <- ScoreJackStraw(seu, dims = 1:dims)
    
    # Reset to sequential processing
    plan(sequential)
    
    js_pvals <- seu[["pca"]]@misc$jackstraw$overall.p.values
    # Only consider PCs >= CFG_PCA_MIN_DIMS (10-50 searchspace per Xiang et al. 2021)
    pcs_in_range <- CFG_PCA_MIN_DIMS:dims
    sig_pcs <- pcs_in_range[js_pvals[pcs_in_range, 2] < 0.05]  # p < 0.05 threshold
    
    if (length(sig_pcs) >= 5) {
      significant_pcs <- sig_pcs
      success <- TRUE
      log(sprintf("JackStraw identified %d significant PCs (p < 0.05) in range %d-%d", 
                  length(sig_pcs), CFG_PCA_MIN_DIMS, dims))
    } else {
      log(sprintf("JackStraw found <%d significant PCs in range %d-%d, will use fallback dimensions",
                  5, CFG_PCA_MIN_DIMS, dims))
    }
  }, error = function(e) {
    log(sprintf("JackStraw failed: %s", e$message))
    log("Will use fallback PC dimensions")
  }, finally = {
    # Always restore limit
    options(future.globals.maxSize = old_limit)
  })
  
  return(list(
    seu = seu,
    significant_pcs = significant_pcs,
    success = success
  ))
}



#' Complete PCA analysis with JackStraw and visualization
#' 
#' Runs PCA, optionally performs JackStraw for automatic PC selection,
#' and falls back to provided dimensions if JackStraw fails. (30)
#' 
#' @param seu Seurat object
#' @param pca_search_space Maximum number of PCs to compute
#' @param pc_dims Fallback PC dimensions if JackStraw fails
#' @return List with:
#'   - seu: Seurat object with PCA
#'   - selected_pcs: Vector of selected PC indices
#'   - jackstraw_used: Whether JackStraw selection was used
perform_pca_analysis <- function(seu, pca_search_space, pc_dims) {
  log(sprintf("Starting PCA analysis (searchspace: %d-%d PCs)", CFG_PCA_MIN_DIMS, pca_search_space))
  
  # Run PCA
  seu <- run_pca(seu, npcs = pca_search_space)
  
  jackstraw_used <- FALSE
  selected_pcs <- pc_dims
  
  # Try JackStraw for automatic PC selection (if enabled)
  if (CFG_JACKSTRAW_ENABLED) {
    log("Attempting JackStraw analysis for PC selection...")
    js_result <- perform_jackstraw(seu, dims = pca_search_space)
    seu <- js_result$seu
    
    if (js_result$success) {
      selected_pcs <- js_result$significant_pcs
      jackstraw_used <- TRUE
      log(paste("JackStraw successful - using", length(selected_pcs), "significant PCs"))
    } else {
      log("JackStraw not successful - using fallback PC dimensions:", paste(range(pc_dims), collapse = "-"))
    }
  } else {
    log("JackStraw disabled - using fallback PC dimensions:", paste(range(pc_dims), collapse = "-"))
  }
  
  return(list(
    seu = seu,
    selected_pcs = selected_pcs,
    jackstraw_used = jackstraw_used
  ))
}

#' Create UMAP embedding
#' 
#' Projects high-dimensional PCA space into 2D using UMAP
#' 
#' @param seu Seurat object with PCA computed
#' @param pc_dims PC dimensions to use for UMAP
#' @param umap_neighbors Number of nearest neighbors for UMAP
#' @param reduction_name Name for UMAP reduction (default: "umap")
#' @param seed Random seed for reproducibility
#' @return Seurat object with UMAP reduction
create_umap <- function(seu, pc_dims, umap_neighbors, reduction_name = "umap", seed = CFG_RANDOM_SEED) {
  log(sprintf("Creating UMAP '%s' using %d PCs...", reduction_name, length(pc_dims)))
  seu <- RunUMAP(seu, dims = pc_dims, n.neighbors = umap_neighbors, 
                 reduction.name = reduction_name, seed.use = seed, verbose = FALSE)
  return(seu)
}

#' Perform multi-resolution clustering
#' 
#' Builds shared nearest neighbor (SNN) graph and performs Louvain clustering
#' 
#' @param seu Seurat object with PCA computed
#' @param pc_dims PC dimensions to use for graph construction
#' @param resolutions Vector of clustering resolutions to test
#' @param seed Random seed for reproducibility
#' @return Seurat object with clustering results at each resolution
perform_multires_clustering <- function(seu, pc_dims, resolutions, seed = CFG_RANDOM_SEED) {
  log("Building nearest neighbor graph...")
  # Set seed before FindNeighbors for SNN graph reproducibility
  set.seed(seed)
  seu <- FindNeighbors(seu, dims = pc_dims, verbose = FALSE)
  
  log("Performing clustering at multiple resolutions...")
  for (res in resolutions) {
    log(sprintf("  Clustering at resolution: %.2f", res))
    # algorithm = 1: Louvain 
    seu <- FindClusters(seu, resolution = res, algorithm = 1, 
                       random.seed = seed, verbose = FALSE)
  }
  
  return(seu)
}

#' Calculate clustering quality metrics
#' 
#' Computes multiple quality metrics for each clustering resolution:
#' - Silhouette score (cluster separation)
#' - WCSS (Within-Cluster Sum of Squares - compactness)
#' - BCSS (Between-Cluster Sum of Squares - separation)
#' - Inertia ratio (BCSS/(BCSS+WCSS))
#' 
#' @param seu Seurat object with clustering results
#' @param resolutions Vector of resolutions to evaluate
#' @param pc_dims PC dimensions used for clustering
#' @param max_cells_for_metrics Maximum cells to sample for efficiency
#' @return Data frame with quality metrics per resolution
calculate_clustering_quality <- function(seu, resolutions, pc_dims, max_cells_for_metrics = CFG_MAX_CELLS_METRICS) {
  log("Calculating clustering quality metrics...")
  
  if (!requireNamespace("cluster", quietly = TRUE)) {
    log("Warning: 'cluster' package not available, skipping quality metrics")
    return(data.frame())
  }
  
  # Sample cells if needed
  n_cells <- ncol(seu)
  idx <- seq_len(n_cells)
  if (n_cells > max_cells_for_metrics) {
    set.seed(CFG_RANDOM_SEED)
    idx <- sort(sample(idx, max_cells_for_metrics))
    log(sprintf("Sampling %d of %d cells for quality metrics", length(idx), n_cells))
  }
  
  pca_coords <- Embeddings(seu, "pca")[idx, pc_dims, drop = FALSE]
  dist_mat <- dist(pca_coords)
  
  quality_metrics <- data.frame()
  
  for (res in resolutions) {
    col <- find_resolution_column(seu, res)
    if (is.null(col)) next
    
    clusters <- seu@meta.data[[col]][idx]
    cluster_int <- as.integer(factor(clusters))
    n_clusters <- length(unique(clusters))
    
    if (n_clusters <= 1) next
    
      # Silhouette score
      sil <- cluster::silhouette(cluster_int, dist_mat)
      avg_silhouette <- mean(sil[, 3])
      
      # WCSS and BCSS
      overall_centroid <- colMeans(pca_coords)
      wcss <- 0
      bcss <- 0
      
      for (k in unique(cluster_int)) {
        cluster_data <- pca_coords[cluster_int == k, , drop = FALSE]
        cluster_size <- nrow(cluster_data)
        
        if (cluster_size > 1) {
          cluster_centroid <- colMeans(cluster_data)
          wcss <- wcss + sum(apply(cluster_data, 1, function(x) sum((x - cluster_centroid)^2)))
          bcss <- bcss + cluster_size * sum((cluster_centroid - overall_centroid)^2)
        }
      }
      
      inertia_ratio <- bcss / (bcss + wcss)
      
      quality_metrics <- rbind(quality_metrics, data.frame(
        Resolution = res,
        NumClusters = n_clusters,
        AvgSilhouette = round(avg_silhouette, 4),
        WCSS = round(wcss, 2),
        BCSS = round(bcss, 2),
        InertiaRatio = round(inertia_ratio, 4)
      ))
      
  }
  
  return(quality_metrics)
}

#' Select optimal clustering resolution
#' 
#' Chooses best resolution using clustering stability analysis 
#' Analyzes tree structure to identify resolution with minimal cluster splitting/merging
#' and maximal stability (high proportion of cells maintaining cluster membership).
#' 
#' @param seu Seurat object with clustering at multiple resolutions
#' @param resolutions Vector of tested resolutions
#' @param quality_metrics Data frame of quality metrics (unused, kept for compatibility)
#' @return List with:
#'   - resolution: Selected optimal resolution
#'   - method: Method used for selection
select_optimal_resolution <- function(seu, resolutions, quality_metrics = NULL) {
  log("Selecting optimal clustering resolution using clustree...")
  
  if (!requireNamespace("clustree", quietly = TRUE)) {
    stop("clustree package is required for resolution selection. Install with: install.packages('clustree')")
  }
  
  # Create clustree and extract tree data
  ct <- clustree::clustree(seu, prefix = paste0(DefaultAssay(seu), "_snn_res."))
  tree_data <- ct$data
  
  if (is.null(tree_data) || nrow(tree_data) == 0) {
    log("Warning: clustree data is empty, using median resolution as fallback")
    stable_res <- resolutions[ceiling(length(resolutions) / 2)]
    return(list(resolution = stable_res, method = "clustree_fallback"))
  }
  
  # Calculate stability metrics for each resolution
  # Lower resolution edges = more stable clustering
  stability_scores <- data.frame()
  
  for (i in seq_along(resolutions)) {
    res <- resolutions[i]
    
    # Count transitions FROM this resolution (outgoing edges)
    res_pattern <- paste0("_", res, "$")
    outgoing <- sum(grepl(res_pattern, tree_data$node))
    
    # Calculate mean proportion of cells staying together (in_prop)
    # High in_prop = stable clusters
    res_edges <- tree_data[grepl(res_pattern, tree_data$node), ]
    mean_in_prop <- if(nrow(res_edges) > 0) mean(res_edges$in_prop, na.rm = TRUE) else 0
    
    # Calculate stability score: prefer fewer splits and higher in_prop
    # Normalize by number of clusters to avoid bias toward low resolutions
    n_clusters <- length(unique(seu@meta.data[[find_resolution_column(seu, res)]]))
    
    if (n_clusters > 1) {
      # Score: high in_prop is good, low edge count relative to clusters is good
      edges_per_cluster <- outgoing / n_clusters
      stability_score <- mean_in_prop - (0.3 * edges_per_cluster)  # Weight transitions less than coherence
    } else {
      stability_score <- -Inf  # Avoid single cluster solutions
    }
    
    stability_scores <- rbind(stability_scores, data.frame(
      Resolution = res,
      NumClusters = n_clusters,
      Outgoing = outgoing,
      MeanInProp = round(mean_in_prop, 3),
      EdgesPerCluster = round(edges_per_cluster, 3),
      StabilityScore = round(stability_score, 3)
    ))
  }
  
  # Select resolution with highest stability score
  best_idx <- which.max(stability_scores$StabilityScore)
  stable_res <- stability_scores$Resolution[best_idx]
  
  log("Clustree stability analysis:")
  for (i in 1:nrow(stability_scores)) {
    log(sprintf("  Res %.2f: %d clusters, %.3f in_prop, %.3f stability%s",
                stability_scores$Resolution[i],
                stability_scores$NumClusters[i],
                stability_scores$MeanInProp[i],
                stability_scores$StabilityScore[i],
                ifelse(i == best_idx, " <- SELECTED", "")))
  }
  
  log(sprintf("Selected resolution %.2f with stability score %.3f", 
              stable_res, stability_scores$StabilityScore[best_idx]))
  
  return(list(
    resolution = stable_res,
    method = "clustree",
    stability_metrics = stability_scores
  ))
}

#' Perform SingleR cell type annotation
#' 
#' Automatically annotates cell types using reference transcriptomic datasets.
#' Uses correlation-based scoring to match cells to reference cell types.
#' Reference must be pre-downloaded (saves time, avoids dependency issues).
#' Uses ref$label.main for annotation as specified in pipeline requirements.
#' 
#' @param seu Seurat object with normalized data
#' @param species Species identifier ("mouse" or "human")
#' @param singler_ref_path Path to pre-downloaded reference (optional)
#' @return Seurat object with singler_main, singler_fine, singler_main_score metadata
perform_singler_annotation <- function(seu, species, singler_ref_path = NULL) {
  log("Running SingleR cell type annotation...")
  
  if (!requireNamespace("SingleR", quietly = TRUE) || !requireNamespace("celldex", quietly = TRUE)) {
    log("Warning: SingleR or celldex not available, skipping annotation")
    return(seu)
  }
  
  sce <- as.SingleCellExperiment(seu)
  
  # Load reference
  if (!is.null(singler_ref_path) && file.exists(singler_ref_path)) {
    log("Loading pre-downloaded SingleR reference from:", singler_ref_path)
    ref <- readRDS(singler_ref_path)
  } else {
    log(sprintf("Downloading %s reference dataset...", species))
    if (species == "mouse") {
      ref <- celldex::MouseRNAseqData()
    } else {
      ref <- celldex::HumanPrimaryCellAtlasData()
    }
  }
  
  # Run annotation
  log("Annotating cells (this may take several minutes)...")
  # Using ref$label.main as specified in pipeline requirements
  pred_main <- SingleR::SingleR(test = sce, ref = ref, labels = ref$label.main)
  pred_fine <- SingleR::SingleR(test = sce, ref = ref, labels = ref$label.fine)
  
  # Add to Seurat object
  seu$singler_main <- pred_main$labels
  seu$singler_fine <- pred_fine$labels
  seu$singler_main_score <- pred_main$scores[cbind(seq_len(nrow(pred_main)), 
                                                     match(pred_main$labels, colnames(pred_main$scores)))]
  
  log("SingleR annotation completed")
  log("Main types identified:", length(unique(pred_main$labels)))
  log("Fine types identified:", length(unique(pred_fine$labels)))
  
  return(seu)
}

#' Find marker genes for all clusters
#' 
#' Identifies differentially expressed genes that distinguish each cluster
#' 
#' @param seu Seurat object with clustering
#' @param resolution_col Metadata column name for cluster assignments
#' @return Data frame with marker genes (gene, cluster, avg_log2FC, p_val_adj, etc.)
find_cluster_markers <- function(seu, resolution_col) {
  log("Finding marker genes for all clusters...")
  
  # Increase memory limit for large objects in parallel processing
  old_limit <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = 8000 * 1024^2)  # 8 GB
  
  Idents(seu) <- resolution_col
  
  # Enable parallel processing for FindMarkers
  plan(multisession, workers = 8)
  
  # Use only highly variable features for faster computation
  features_to_test <- VariableFeatures(seu)
  if (length(features_to_test) == 0) {
    features_to_test <- rownames(seu)
  }
  log(sprintf("Testing %d features for markers", length(features_to_test)))
  
  # Downsample cells per cluster if dataset is large
  n_cells <- ncol(seu)
  max_cells_per_cluster <- if (n_cells > 50000) 500 else if (n_cells > 20000) 1000 else NULL
  
  if (!is.null(max_cells_per_cluster)) {
    log(sprintf("Downsampling to max %d cells per cluster for efficiency", max_cells_per_cluster))
  }
  
  markers <- FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    max.cells.per.ident = max_cells_per_cluster,
    features = features_to_test,
    test.use = "wilcox",  # Wilcoxon test (Pullin & McCarthy 2024), presto auto-used if available
    verbose = FALSE
  )
  
  # Return to sequential processing
  plan(sequential)
  
  # Restore old limit
  options(future.globals.maxSize = old_limit)
  
  # Check if markers dataframe is empty or has no rows
  if (nrow(markers) == 0) {
    log("Warning: No marker genes found. Creating empty dataframe with expected structure.")
    markers <- data.frame(
      cluster = character(0),
      gene = character(0),
      avg_log2FC = numeric(0),
      pct.1 = numeric(0),
      pct.2 = numeric(0),
      p_val_adj = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  log(sprintf("Found %d marker genes across %d clusters", 
              nrow(markers), length(unique(markers$cluster))))
  
  return(markers)
}

#' Helper function to find resolution column name
#' 
#' Searches metadata for clustering column matching given resolution.
#' Checks common prefixes (integrated_snn_res., RNA_snn_res., SCT_snn_res.).
#' 
#' @param seu Seurat object
#' @param res Resolution value to search for
#' @return Column name if found, NULL otherwise
find_resolution_column <- function(seu, res) {
  for (prefix in c("integrated_snn_res.", "RNA_snn_res.", "SCT_snn_res.")) {
    col_name <- paste0(prefix, res)
    if (col_name %in% colnames(seu@meta.data)) {
      return(col_name)
    }
  }
  return(NULL)
}

#' Improved cell type assignment using cell-first approach
#' 
#' Assigns cell types using two-step approach:
#' 1. Individual cell annotations from SingleR (fine labels)
#' 2. Cluster-level consensus based on majority vote within each cluster
#' Combines with top marker genes for validation.
#' 
#' @param seurat_obj Seurat object with clustering and SingleR annotations
#' @param markers Data frame from FindAllMarkers with cluster-specific genes
#' @return List with:
#'   - seurat_obj: Updated Seurat object with cell type assignments
#'   - cluster_mappings: Data frame mapping clusters to cell types
assign_celltypes_improved <- function(seurat_obj, markers) {
  log("Improved cell type assignment (cell-first approach)...")
  
  # Step 1: Assign individual cell types from SingleR (cell-level annotations)
  n_cells <- ncol(seurat_obj)
  
  if ("singler_fine" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$individual_celltype <- seurat_obj$singler_fine
    seurat_obj$individual_confidence <- 1.0  # SingleR confidence
    log("Individual cell type annotations assigned from SingleR fine labels")
  } else {
    seurat_obj$individual_celltype <- paste0("Cluster_", seurat_obj$seurat_clusters)
    seurat_obj$individual_confidence <- 0
    log("Warning: No SingleR annotations available, using cluster IDs")
  }
  
  # Step 2: For each cluster, find consensus cell type using top marker genes
  cluster_mappings <- data.frame()
  all_clusters <- unique(seurat_obj$seurat_clusters)
  
  for (cluster_id in all_clusters) {
    cluster_cells <- seurat_obj$seurat_clusters == cluster_id
    
    if (sum(cluster_cells) == 0) next
    
    # Get individual cell type annotations for this cluster
    cluster_individual_types <- seurat_obj$individual_celltype[cluster_cells]
    cluster_individual_types <- cluster_individual_types[!is.na(cluster_individual_types)]
    
    # Find top 2-3 most common cell types in cluster
    if (length(cluster_individual_types) > 0) {
      type_table <- table(cluster_individual_types)
      # Get top 3 most common types
      top_types <- head(names(sort(type_table, decreasing = TRUE)), 3)
      dominant_type <- top_types[1]
      confidence <- max(type_table) / sum(type_table)
      
      # Get top marker genes for this cluster
      cluster_markers <- markers[markers$cluster == cluster_id, ]
      if (nrow(cluster_markers) > 0) {
        top_markers <- head(cluster_markers$gene, 3)
        avg_logfc <- round(mean(head(cluster_markers$avg_log2FC, 3)), 2)
      } else {
        top_markers <- "No markers found"
        avg_logfc <- 0
      }
      
      cluster_mappings <- rbind(cluster_mappings, data.frame(
        Cluster = cluster_id,
        N_cells = sum(cluster_cells),
        Dominant_Type = dominant_type,
        Confidence = round(confidence, 3),
        Alternative_Types = paste(top_types[-1], collapse = ", "),
        Top_Markers = paste(top_markers, collapse = ", "),
        Avg_LogFC = avg_logfc,
        stringsAsFactors = FALSE
      ))
    } else {
      cluster_mappings <- rbind(cluster_mappings, data.frame(
        Cluster = cluster_id,
        N_cells = sum(cluster_cells),
        Dominant_Type = "Unknown",
        Confidence = 0,
        Alternative_Types = "",
        Top_Markers = "",
        Avg_LogFC = 0,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Step 3: Assign final cluster-based cell types
  seurat_obj$data_driven_celltype <- rep("Unknown", n_cells)
  seurat_obj$celltype_confidence <- rep(0, n_cells)
  
  names(seurat_obj$data_driven_celltype) <- colnames(seurat_obj)
  names(seurat_obj$celltype_confidence) <- colnames(seurat_obj)
  
  for (i in 1:nrow(cluster_mappings)) {
    cluster_id <- cluster_mappings$Cluster[i]
    cluster_cells <- seurat_obj$seurat_clusters == cluster_id
    
    celltype <- cluster_mappings$Dominant_Type[i]
    confidence <- cluster_mappings$Confidence[i]
    
    cell_names_in_cluster <- colnames(seurat_obj)[cluster_cells]
    seurat_obj$data_driven_celltype[cell_names_in_cluster] <- celltype
    seurat_obj$celltype_confidence[cell_names_in_cluster] <- confidence
  }
  
  # Also set cell_type for backward compatibility
  seurat_obj$cell_type <- seurat_obj$data_driven_celltype
  
  log("Cell type assignment completed")
  log("Individual annotations preserved in 'individual_celltype'")
  log("Cluster consensus in 'data_driven_celltype'")
  
  return(list(
    seurat_obj = seurat_obj,
    cluster_celltype_map = cluster_mappings
  ))
}

# Complete PCA, Clustering and Annotation Pipeline
# Integrates: PCA analysis with JackStraw, UMAP embedding, multi-resolution clustering,
# quality metrics, SingleR annotation, and cell type assignment
# Returns: List with annotated Seurat object, selected resolution, quality metrics, and markers
perform_pca_clustering_annotation <- function(seu, 
                                              pca_search_space = 50,
                                              pc_dims = 10:50,  # Default 10-50 searchspace per literature
                                              resolutions = c(0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 1.1, 1.3),
                                              umap_neighbors = 30,
                                              species = "mouse",
                                              singler_ref_path = NULL,
                                              max_cells_for_metrics = 5000,
                                              seed = 42) {
  
  log("=== Starting PCA, Clustering and Annotation Pipeline ===")
  
  # Step 1: PCA Analysis
  log("Step 1: Running PCA analysis...")
  pca_result <- perform_pca_analysis(seu, pca_search_space, pc_dims)
  seu <- pca_result$seu
  selected_pcs <- pca_result$selected_pcs
  log(sprintf("PCA completed - using %d PCs", length(selected_pcs)))
  
  # Step 2: UMAP Embedding
  log("Step 2: Creating UMAP embedding...")
  seu <- create_umap(seu, selected_pcs, umap_neighbors, seed = seed)
  log("UMAP completed")
  
  # Step 3: Multi-resolution Clustering
  log("Step 3: Performing multi-resolution clustering...")
  seu <- perform_multires_clustering(seu, selected_pcs, resolutions, seed = seed)
  log(sprintf("Clustering completed at %d resolutions", length(resolutions)))
  
  # Step 4: Calculate Clustering Quality Metrics
  log("Step 4: Calculating clustering quality metrics...")
  quality_metrics <- calculate_clustering_quality(seu, resolutions, selected_pcs, max_cells_for_metrics)
  
  # Step 5: Select Optimal Resolution
  log("Step 5: Selecting optimal clustering resolution...")
  optimal_res <- select_optimal_resolution(seu, resolutions, quality_metrics)
  selected_resolution <- optimal_res$resolution
  log(sprintf("Selected resolution: %.2f (method: %s)", selected_resolution, optimal_res$method))
  
  # Set identity to selected resolution
  res_col <- find_resolution_column(seu, selected_resolution)
  if (!is.null(res_col)) {
    Idents(seu) <- res_col
    seu$seurat_clusters <- Idents(seu)
  }
  
  # Step 6: Find Cluster Markers
  log("Step 6: Finding cluster marker genes...")
  markers <- find_cluster_markers(seu, res_col)
  
  # Step 7: SingleR Annotation
  log("Step 7: Performing SingleR cell type annotation...")
  seu <- perform_singler_annotation(seu, species, singler_ref_path)
  
  # Step 8: Assign Cell Types (combine SingleR and markers)
  log("Step 8: Assigning final cell types...")
  celltype_result <- assign_celltypes_improved(seu, markers)
  seu <- celltype_result$seurat_obj
  cluster_celltype_map <- celltype_result$cluster_celltype_map
  
  log("=== Clustering and annotation completed ===")
  log(sprintf("Summary: %d cells, %d clusters, %d cell types identified",
              ncol(seu), length(unique(seu$seurat_clusters)), length(unique(seu$cell_type))))
  
  return(list(
    seu = seu,
    selected_resolution = selected_resolution,
    selected_pcs = selected_pcs,
    quality_metrics = quality_metrics,
    markers = markers,
    cluster_celltype_map = cluster_celltype_map,
    jackstraw_used = pca_result$jackstraw_used
  ))
}

#===============================================================================
# 0.5. HELPER AND UTILITY FUNCTIONS
#===============================================================================

#===============================================================================
# 3. COMMAND-LINE ARGUMENT PARSING
#===============================================================================
# Parse command-line arguments and load configuration from nextflow.config
# Command-line args can override config values
#===============================================================================

# Parameters are read directly from nextflow.config
# Command-line args can override config values
option_list <- list(
  make_option(c("-i", "--input_rds"), type="character", default=NULL, help="Path to integrated Seurat RDS object"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Output directory"),
  make_option(c("-r", "--resolutions"), type="character", default=NULL, help="Clustering resolutions (comma-separated)"),
  make_option(c("-p", "--pc_dims"), type="character", default=NULL, help="PC dimensions (e.g., '30' or '1:30' or '1,2,5')"),
  make_option(c("--pca_search_space"), type="integer", default=NULL, help="Number of PCs to test with JackStraw"),
  make_option(c("--umap_neighbors"), type="integer", default=NULL, help="Number of UMAP neighbors"),
  make_option(c("-s", "--species"), type="character", default=NULL, help="Species: mouse or human"),
  make_option(c("--singler_ref"), type="character", default=NULL, help="Path to pre-downloaded SingleR reference"),
  make_option(c("--seed"), type="integer", default=NULL, help="Random seed (0 = random, or specify integer)"),
  make_option(c("--config"), type="character", default=NULL, help="Path to nextflow.config file")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Read parameters from nextflow.config
nf_params <- parse_nextflow_config(opt$config)

# Use config values as defaults, but allow command-line overrides
opt$output_dir <- opt$output_dir %||% nf_params$output_dir %||% "outputs/03_clustering"
opt$resolutions <- opt$resolutions %||% nf_params$resolutions %||% "0.1,0.3,0.5,0.7,0.8,0.9,1.1,1.3"
opt$pc_dims <- opt$pc_dims %||% nf_params$pc_dims %||% "10-50"  # PCA searchspace 10-50 as per literature
opt$pca_search_space <- opt$pca_search_space %||% as.integer(nf_params$pca_search_space) %||% 50L
opt$umap_neighbors <- opt$umap_neighbors %||% as.integer(nf_params$umap_neighbors) %||% 30L
opt$species <- opt$species %||% nf_params$species %||% "mouse"
opt$singler_ref <- opt$singler_ref %||% nf_params$singler_ref
opt$seed <- opt$seed %||% as.integer(nf_params$seed) %||% as.integer(nf_params$clustering.random_seed) %||% 42L

#===============================================================================
# 4. INITIALIZE ANALYSIS PARAMETERS
#===============================================================================
# Load parameters from nextflow.config (params.clustering block)
# Set global configuration variables
#===============================================================================

# Reproducibility
CFG_RANDOM_SEED <- as.integer(nf_params$clustering.random_seed %||% opt$seed %||% 42)

# PCA Parameters
CFG_PCA_MIN_DIMS <- as.integer(nf_params$clustering.pca_min_dims %||% 10)
CFG_PCA_SEARCH_SPACE <- as.integer(nf_params$clustering.pca_search_space %||% opt$pca_search_space %||% 50)
CFG_JACKSTRAW_ENABLED <- as.logical(nf_params$clustering.jackstraw_enabled %||% TRUE)
CFG_JACKSTRAW_REPLICATES <- as.integer(nf_params$clustering.jackstraw_replicates %||% 100)
CFG_JACKSTRAW_CORES <- as.integer(nf_params$clustering.jackstraw_cores %||% 8)

# Marker Gene Parameters
CFG_MARKER_WORKERS <- as.integer(nf_params$clustering.marker_workers %||% 8)
marker_max_val <- nf_params$clustering.marker_max_cells_per_cluster
CFG_MARKER_MAX_CELLS_PER_CLUSTER <- if(!is.null(marker_max_val) && marker_max_val != "null") as.integer(marker_max_val) else NULL

# Clustering Parameters
CFG_K_NEIGHBORS <- as.integer(nf_params$clustering.k_neighbors %||% 20)
CFG_MAX_CELLS_METRICS <- as.integer(nf_params$clustering.max_cells_metrics %||% 5000)
CFG_TOP_LOADING_GENES <- as.integer(nf_params$clustering.top_loading_genes %||% 10)

# Set global random seed
set.seed(CFG_RANDOM_SEED)

# Set future globals size limit to unlimited for large Seurat objects
options(future.globals.maxSize = Inf)

# Parse user-provided parameters
resolutions      <- parse_resolutions(opt$resolutions)
pc_dims          <- parse_dims(opt$pc_dims)
pca_search_space <- opt$pca_search_space
umap_neighbors   <- opt$umap_neighbors
species          <- tolower(opt$species)

# Random seed handling
seed_param <- opt$seed
if (seed_param == 0) {
  analysis_seed <- sample(1:10000, 1)
  set.seed(analysis_seed)
  log("Using random seed:", analysis_seed)
} else {
  analysis_seed <- seed_param
  set.seed(analysis_seed)
  log("Using specified seed:", analysis_seed)
}

log("Parameters loaded from nextflow.config:")
log("  - resolutions:", paste(resolutions, collapse=", "))
log("  - pc_dims:", paste(range(pc_dims), collapse="-"))
log("  - pca_search_space:", pca_search_space)
log("  - umap_neighbors:", umap_neighbors)
log("  - species:", species)
log("  - seed:", analysis_seed)

#===============================================================================
# 5. SETUP OUTPUT DIRECTORIES AND LOGGING
#===============================================================================
# Create directory structure and initialize logging
#===============================================================================

module_dir     <- opt$output_dir
data_dir       <- file.path(module_dir, "data")
plots_dir      <- file.path(module_dir, "plots")
dash_dir       <- file.path(module_dir, "interactive")
rds_dir        <- file.path(module_dir, "rds_objects")
final_obj_path <- file.path(module_dir, "clustered_annotated_seurat.rds")

# Start performance tracking
perf_tracker <- start_performance_tracking("02_clustering_annotation", module_dir)

for (d in c(module_dir, data_dir, plots_dir, dash_dir, rds_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

log_file <- file.path(module_dir, "02_clustering_annotation.log")
setup_logging(log_file)

#===============================================================================
# 6. LOG ANALYSIS PARAMETERS
#===============================================================================
log("==================================================")
log("Starting Clustering and Cell Type Annotation")
log("==================================================")
log("Output directory:", module_dir)
log("Species:", species)
log("PCA search space:", pca_search_space, "PCs")
log("Clustering resolutions:", paste(resolutions, collapse = ", "))
log("PC dims:", paste(min(pc_dims), "-", max(pc_dims)))
log("UMAP neighbors:", umap_neighbors)
log("==================================================")

#===============================================================================
# 7. LOAD INTEGRATED SEURAT OBJECT
#===============================================================================
log("Loading integrated Seurat object...")

seurat_paths <- character(0)
if (!is.null(opt$input_rds) && file.exists(opt$input_rds)) {
  seurat_paths <- c(seurat_paths, opt$input_rds)
} else {
  seurat_paths <- c(
    file.path(opt$output_dir, "seurat_integrated.rds"),
    file.path(opt$output_dir, "data", "seurat_integrated.rds"),
    file.path(dirname(opt$output_dir), "02_qc", "data", "seurat_integrated.rds")
  )
}

seu <- NULL
for (path in seurat_paths) {
  if (file.exists(path)) {
    log(paste("Loading Seurat object from:", path))
    tryCatch({
      seu <- readRDS(path)
      log(paste("Successfully loaded Seurat object with", ncol(seu), "cells"))
      break
    }, error = function(e) {
      log(paste("ERROR loading", path, ":", e$message))
      seu <- NULL
    })
  } else {
    log(paste("Path does not exist:", path))
  }
}

if (is.null(seu)) {
  stop("Could not find or load Seurat object. Tried paths: ", paste(seurat_paths, collapse = ", "))
}

# Set default assay
tryCatch({
  DefaultAssay(seu) <- if ("integrated" %in% names(seu@assays)) "integrated" else DefaultAssay(seu)
  log(paste("Using assay:", DefaultAssay(seu)))
}, error = function(e) {
  log(paste("Warning: Could not set DefaultAssay:", e$message))
})

# Factor ordering for consistent visualization
if ("condition" %in% colnames(seu@meta.data)) {
  seu$condition <- factor(seu$condition)
  log("Condition levels:", paste(levels(seu$condition), collapse = ", "))
}
if ("sample" %in% colnames(seu@meta.data)) {
  seu$sample <- factor(seu$sample)
  log("Sample levels:", paste(levels(seu$sample), collapse = ", "))
}

#===============================================================================
# 8. CELL CYCLE SCORING
#===============================================================================
log("==================================================")
log("Performing Cell Cycle Scoring")
log("==================================================")

seu <- perform_cell_cycle_scoring(seu, species)

log("==================================================")

#===============================================================================
# 9. INITIAL PCA ANALYSIS AND VISUALIZATION
#===============================================================================
log("Running PCA analysis...")

# Run PCA
seu <- run_pca(seu, npcs = pca_search_space)

# Create initial PCA visualization plots (before JackStraw)
pca_plots <- create_pca_plots(seu, pc_dims, jackstraw_used = FALSE)

# Save plots (static + interactive)
save_dual_plot(pca_plots$variance, "pca_variance_analysis", plots_dir, dash_dir, width = 12, height = 6)
save_dual_plot(pca_plots$elbow, "pca_elbow_plot", plots_dir, dash_dir, width = 10, height = 6)

if (!is.null(pca_plots$pca_condition)) {
  save_dual_plot(pca_plots$pca_condition, "pca_condition", plots_dir, dash_dir, width = 10, height = 8)
}

if (!is.null(pca_plots$pca_sample)) {
  save_dual_plot(pca_plots$pca_sample, "pca_sample", plots_dir, dash_dir, width = 12, height = 8)
}

save_dual_plot(pca_plots$loadings, "pca_top_loadings", plots_dir, dash_dir, width = 14, height = 8)

log("PCA visualization completed (static + interactive)")

#===============================================================================
# 10. RUN COMPLETE CLUSTERING AND ANNOTATION PIPELINE
#===============================================================================
# Executes the integrated pipeline:
# PCA (with JackStraw) → UMAP → Multi-res Clustering → SingleR → Markers
#===============================================================================

log("Running complete PCA, Clustering and Annotation pipeline...")

# Run the integrated pipeline function
pipeline_result <- perform_pca_clustering_annotation(
  seu = seu,
  pca_search_space = pca_search_space,
  pc_dims = pc_dims,
  resolutions = resolutions,
  umap_neighbors = umap_neighbors,
  species = species,
  singler_ref_path = opt$singler_ref,
  max_cells_for_metrics = 5000,
  seed = analysis_seed
)

# Extract results
seu <- pipeline_result$seu
selected_resolution <- pipeline_result$selected_resolution
quality_metrics <- pipeline_result$quality_metrics
markers <- pipeline_result$markers
cluster_celltype_map <- pipeline_result$cluster_celltype_map
jackstraw_used <- pipeline_result$jackstraw_used

log(sprintf("Pipeline completed: %d clusters, %d cell types identified",
            length(unique(seu$seurat_clusters)), length(unique(seu$cell_type))))

# Create updated PCA plots with JackStraw information
log("Creating updated PCA plots with JackStraw information...")
pca_plots_final <- create_pca_plots(seu, pipeline_result$selected_pcs, jackstraw_used)
save_dual_plot(pca_plots_final$variance, "pca_variance_analysis_final", plots_dir, dash_dir, width = 12, height = 6)
save_dual_plot(pca_plots_final$elbow, "pca_elbow_plot_final", plots_dir, dash_dir, width = 10, height = 6)
log("Updated PCA plots saved with JackStraw selection info")

#===============================================================================
# 11. VISUALIZATION AND QUALITY ASSESSMENT
#===============================================================================

log("Creating visualization and quality assessment plots...")

#-------------------------------------------
# 11.1. Clustering Quality Metrics Visualization
#-------------------------------------------
log("Creating clustering quality metrics visualization...")

if (!is.null(quality_metrics) && nrow(quality_metrics) > 1) {
    # Quality metrics comparison plot
    # Create quality metrics plots
    quality_plots <- create_quality_metrics_plots(quality_metrics)
    save_dual_plot(quality_plots$quality_metrics, "clustering_quality_metrics", plots_dir, dash_dir, width = 12, height = 6)
    save_dual_plot(quality_plots$quality_vs_clusters, "clustering_quality_vs_clusters", plots_dir, dash_dir, width = 10, height = 6)
    log("Saved clustering quality visualization plots (static + interactive)")
    
}

#-------------------------------------------
# 11.2. Clustering Stability Analysis
#-------------------------------------------
log("Creating clustering stability analysis...")

# Simple stability analysis
stab <- data.frame()
for (res in resolutions) {
  col <- find_resolution_column(seu, res)
  if (is.null(col)) next
  k <- length(unique(seu@meta.data[[col]]))
  stab <- rbind(stab, data.frame(Resolution = res, Clusters = k))
}

# Create stability plot
if (nrow(stab) > 0) {
  p_stab <- create_stability_plot(stab, selected_resolution, length(unique(seu$seurat_clusters)))
  save_dual_plot(p_stab, "cluster_stability", plots_dir, dash_dir, width = 10, height = 6)
  write.csv(stab, file.path(data_dir, "stability_metrics.csv"), row.names = FALSE)
}

#-------------------------------------------
# 11.3. Save Cluster Mapping
#-------------------------------------------
write.csv(cluster_celltype_map, file.path(data_dir, "cluster_to_celltype_mapping.csv"), row.names = FALSE)


#===============================================================================
# 12. UMAP VISUALIZATIONS  
#===============================================================================

log("Creating UMAP visualizations...")

#-------------------------------------------
# 12.1. Basic UMAP Plots by Metadata
#-------------------------------------------

# Create custom color palette for cell types
n_cell_types <- length(unique(seu$cell_type))
log("Number of cell types:", n_cell_types)

if(requireNamespace("ggsci", quietly = TRUE)) {
  cell_type_colors <- ggsci::pal_igv()(min(n_cell_types, 51))
  if(n_cell_types > 51) {
    additional_colors <- viridis::turbo(n_cell_types - 51)
    cell_type_colors <- c(cell_type_colors, additional_colors)
  }
} else {
  cell_type_colors <- viridis::turbo(n_cell_types)
}

# UMAP plots with consistent font sizes
p_cond <- create_umap_plot(seu, group_by = "condition", title = "UMAP: By Condition")
save_dual_plot(p_cond, "umap_condition", plots_dir, dash_dir, width = 12, height = 8)

p_samp <- create_umap_plot(seu, group_by = "sample", title = "UMAP: By Sample")
save_dual_plot(p_samp, "umap_sample", plots_dir, dash_dir, width = 14, height = 8)

p_clus <- create_umap_plot(seu, group_by = "cell_type", title = "UMAP: By Clusters (Cell Type Names)", label = TRUE, repel = TRUE, colors = cell_type_colors)
save_dual_plot(p_clus, "umap_clusters_numeric", plots_dir, dash_dir, width = 14, height = 10)

# Create cluster plot with cell type names
if ("cell_type" %in% colnames(seu@meta.data)) {
  p_clus_named <- create_umap_plot(seu, group_by = "cell_type", 
                                   title = "UMAP: By Clusters (Cell Type Names)",
                                   label = TRUE, repel = TRUE, colors = cell_type_colors)
  save_dual_plot(p_clus_named, "umap_clusters", plots_dir, dash_dir, width = 14, height = 10)
} else {
  # Fallback if cell_type is not available
  save_dual_plot(p_clus, "umap_clusters", plots_dir, dash_dir, width = 12, height = 8)
}

#-------------------------------------------
# 12.2. Cell Cycle Visualizations
#-------------------------------------------

# Cell cycle phase visualization
if ("Phase" %in% colnames(seu@meta.data)) {
  log("Creating cell cycle phase visualizations...")
  
  # Define cell cycle colors
  phase_colors <- c("G1" = "#4DAF4A", "S" = "#377EB8", "G2M" = "#E41A1C", "Unknown" = "grey70")
  
  p_phase <- create_umap_plot(seu, group_by = "Phase", title = "UMAP: Cell Cycle Phase", 
                              colors = phase_colors)
  save_dual_plot(p_phase, "umap_cell_cycle_phase", plots_dir, dash_dir, width = 12, height = 8)
  
  # Cell cycle phase with annotated clusters
  if ("cell_type" %in% colnames(seu@meta.data)) {
    p_phase_labeled <- create_umap_plot_with_labels(seu, group_by = "Phase", 
                                        label_by = "cell_type",
                                        title = "UMAP: Cell Cycle Phase (Cell Type Labels)", 
                                        colors = phase_colors, label_size = 6)
    save_dual_plot(p_phase_labeled, "umap_cell_cycle_phase_annotated", plots_dir, dash_dir, width = 12, height = 8)
  }
  
  # Cell cycle scores
  if ("S.Score" %in% colnames(seu@meta.data)) {
    p_s_score <- create_feature_plot(seu, features = "S.Score", title = "S Phase Score")
    save_dual_plot(p_s_score, "umap_s_phase_score", plots_dir, dash_dir, width = 10, height = 8)
    
    # Additional plots with cell type and cluster annotations
    if ("cell_type" %in% colnames(seu@meta.data)) {
      p_s_celltype <- create_feature_plot(seu, features = "S.Score", 
                                          title = "S Phase Score (Cell Types Labeled)",
                                          label = TRUE, label.size = 3, repel = TRUE, 
                                          group.by = "cell_type")
      save_dual_plot(p_s_celltype, "umap_s_phase_score_celltypes", plots_dir, dash_dir, width = 12, height = 8)
    }
    
    # With cluster labels
    p_s_clusters <- create_feature_plot(seu, features = "S.Score", 
                                        title = "S Phase Score (Clusters Labeled)",
                                        label = TRUE, label.size = 4)
    save_dual_plot(p_s_clusters, "umap_s_phase_score_clusters", plots_dir, dash_dir, width = 12, height = 8)
  }
  
  if ("G2M.Score" %in% colnames(seu@meta.data)) {
    p_g2m_score <- create_feature_plot(seu, features = "G2M.Score", title = "G2M Phase Score")
    save_dual_plot(p_g2m_score, "umap_g2m_phase_score", plots_dir, dash_dir, width = 10, height = 8)
    
    # Additional plots with cell type and cluster annotations
    if ("cell_type" %in% colnames(seu@meta.data)) {
      p_g2m_celltype <- create_feature_plot(seu, features = "G2M.Score", 
                                            title = "G2M Phase Score (Cell Types Labeled)",
                                            label = TRUE, label.size = 3, repel = TRUE,
                                            group.by = "cell_type")
      save_dual_plot(p_g2m_celltype, "umap_g2m_phase_score_celltypes", plots_dir, dash_dir, width = 12, height = 8)
    }
    
    # With cluster labels
    p_g2m_clusters <- create_feature_plot(seu, features = "G2M.Score", 
                                          title = "G2M Phase Score (Clusters Labeled)",
                                          label = TRUE, label.size = 4)
    save_dual_plot(p_g2m_clusters, "umap_g2m_phase_score_clusters", plots_dir, dash_dir, width = 12, height = 8)
  }
  
  # Cell cycle distribution by cluster (using cell_type names)
  if ("cell_type" %in% colnames(seu@meta.data)) {
    phase_by_cluster <- seu@meta.data %>%
      group_by(cell_type, Phase) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(cell_type) %>%
      mutate(percentage = count / sum(count) * 100)
    
    p_phase_by_cluster <- create_cell_cycle_barplot(phase_by_cluster, phase_colors)
    save_dual_plot(p_phase_by_cluster, "cell_cycle_by_cluster", plots_dir, dash_dir, width = 12, height = 7)
    
    # Create individual plots for each cell type
    create_cell_cycle_individual_plots(seu, phase_colors, plots_dir, dash_dir, group_by = "cell_type")
  }
  
  # Cell cycle distribution by numeric cluster (if needed)
  # Removed to focus on cell type names instead of seurat_clusters
  
  log("Cell cycle visualizations saved")
}

#-------------------------------------------
# 12.3. Cell Type UMAP and Annotation Comparison
#-------------------------------------------

p_cell <- create_umap_plot(seu, group_by = "cell_type", title = "UMAP: By Cell Type",
                           label = TRUE, repel = TRUE, colors = cell_type_colors)
save_dual_plot(p_cell, "umap_celltype", plots_dir, dash_dir, width = 14, height = 10)

log("Creating cell type annotation comparison plots...")

# Step 7.3.1: Side-by-side UMAP comparison
if ("individual_celltype" %in% colnames(seu@meta.data) && 
    "data_driven_celltype" %in% colnames(seu@meta.data)) {
  
  # Create consistent color mapping based on data_driven_celltype
  unique_celltypes <- sort(unique(seu$data_driven_celltype))
  n_types <- length(unique_celltypes)
  
  if(requireNamespace("ggsci", quietly = TRUE)) {
    consistent_colors <- ggsci::pal_igv()(min(n_types, 51))
    if(n_types > 51) {
      additional_colors <- viridis::turbo(n_types - 51)
      consistent_colors <- c(consistent_colors, additional_colors)
    }
  } else {
    consistent_colors <- viridis::turbo(n_types)
  }
  names(consistent_colors) <- unique_celltypes
  
  # Map individual_celltype to data_driven_celltype colors
  individual_to_datadriven <- seu@meta.data %>%
    select(individual_celltype, data_driven_celltype) %>%
    distinct()
  individual_colors <- consistent_colors[individual_to_datadriven$data_driven_celltype]
  names(individual_colors) <- individual_to_datadriven$individual_celltype
  
  p_comparison <- create_annotation_comparison_sidebyside(seu, individual_colors, consistent_colors)
  save_dual_plot(p_comparison, "annotation_comparison_sidebyside", plots_dir, dash_dir, 
                width = 20, height = 8)
  
  # Step 7.3.2: Confidence distribution
  if ("celltype_confidence" %in% colnames(seu@meta.data)) {
    p_confidence <- create_confidence_plot(seu@meta.data)
    save_dual_plot(p_confidence, "celltype_confidence_distribution", plots_dir, dash_dir, 
                  width = 10, height = 7)
    
    # Create confidence statistics table by cell type
    confidence_stats <- seu@meta.data %>%
      group_by(cell_type) %>%
      summarise(
        n_cells = n(),
        mean_confidence = mean(celltype_confidence, na.rm = TRUE),
        median_confidence = median(celltype_confidence, na.rm = TRUE),
        sd_confidence = sd(celltype_confidence, na.rm = TRUE),
        min_confidence = min(celltype_confidence, na.rm = TRUE),
        max_confidence = max(celltype_confidence, na.rm = TRUE),
        q25_confidence = quantile(celltype_confidence, 0.25, na.rm = TRUE),
        q75_confidence = quantile(celltype_confidence, 0.75, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_confidence))
    
    # Save confidence statistics table
    write.csv(confidence_stats, 
              file.path(data_dir, "celltype_confidence_statistics.csv"), 
              row.names = FALSE)
    log("Saved celltype confidence statistics to celltype_confidence_statistics.csv")
    
    # Create overall confidence summary
    overall_confidence <- data.frame(
      metric = c("Overall_Mean", "Overall_Median", "Overall_SD", "Overall_Min", "Overall_Max"),
      value = c(
        mean(seu$celltype_confidence, na.rm = TRUE),
        median(seu$celltype_confidence, na.rm = TRUE),
        sd(seu$celltype_confidence, na.rm = TRUE),
        min(seu$celltype_confidence, na.rm = TRUE),
        max(seu$celltype_confidence, na.rm = TRUE)
      )
    )
    
    write.csv(overall_confidence, 
              file.path(data_dir, "celltype_confidence_overall.csv"), 
              row.names = FALSE)
    log("Saved overall confidence summary to celltype_confidence_overall.csv")
  }
  
  # Step 7.3.3: Comparison matrix
  p_comparison_matrix <- create_comparison_matrix_plot(seu)
  save_dual_plot(p_comparison_matrix, "annotation_comparison_matrix", plots_dir, dash_dir, 
                width = 14, height = 12)
  
  # Step 7.3.4: Bar plot comparison
  p_bar_comparison <- create_celltype_barplot_comparison(seu)
  save_dual_plot(p_bar_comparison, "annotation_distribution_comparison", plots_dir, dash_dir, 
                width = 18, height = 10)
  
  log("Cell type annotation comparison plots saved")
}

#-------------------------------------------
# 12.4. Interactive Cluster Plot
#-------------------------------------------
if(requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
  p_interactive <- create_interactive_cluster_plot(seu)
  htmlwidgets::saveWidget(
    p_interactive,
    file = file.path(plots_dir, "interactive_cluster_plot.html"),
    selfcontained = TRUE
  )
  log("Interactive cluster plot saved")
}

#-------------------------------------------
# 12.5. Cell Type Centers Plot
#-------------------------------------------
if(requireNamespace("dplyr", quietly = TRUE)) {
  um <- as.data.frame(Seurat::Embeddings(seu, "umap"))
  colnames(um) <- c("UMAP_1", "UMAP_2")
  df_um <- cbind(seu@meta.data, um)
  
  centers <- df_um %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(
      UMAP1_center = mean(UMAP_1),
      UMAP2_center = mean(UMAP_2),
      n_cells = n(),
      .groups = "drop"
    )
  write.csv(centers, file.path(data_dir, "celltype_centers_coordinates.csv"), row.names = FALSE)
  
  p_cent <- create_celltype_centers_plot(centers)
  save_dual_plot(p_cent, "umap_celltype_centers_sizes", plots_dir, dash_dir, width = 12, height = 10)
}

#===============================================================================
# 13. SAVE FINAL OBJECT AND GENERATE SUMMARY STATISTICS
#===============================================================================
# Save final annotated Seurat object and generate comprehensive summaries:
# - Cell type distribution by condition
# - Cluster distribution by condition
# - Annotation quality metrics
# - JSON summary for automated processing
#===============================================================================

log("Saving final annotated object...")
saveRDS(seu, file = final_obj_path)
log("Final annotated object saved to:", final_obj_path)

# Save parameters for MiloR analysis in step 03
log("Saving MiloR parameters...")
milor_params <- list(
  num_pcs = length(pipeline_result$selected_pcs),
  k_neighbors = CFG_K_NEIGHBORS
)
saveRDS(milor_params, file = file.path(data_dir, "analysis_parameters.rds"))
log("MiloR parameters saved:", sprintf("k=%d, d=%d", milor_params$k_neighbors, milor_params$num_pcs))

# Generate final data summaries
log("Generating summary statistics...")

# Cell type distribution summary
cell_type_distribution <- table(seu$cell_type, seu$condition)
write.csv(cell_type_distribution, file.path(data_dir, "cell_type_distribution_by_condition.csv"))

cluster_distribution <- table(seu$seurat_clusters, seu$condition)
write.csv(cluster_distribution, file.path(data_dir, "cluster_distribution_by_condition.csv"))

# Annotation quality metrics
if ("singler_main" %in% colnames(seu@meta.data)) {
  annotation_summary <- data.frame(
    total_cells = ncol(seu),
    unique_clusters = length(unique(seu$seurat_clusters)), 
    unique_cell_types = length(unique(seu$singler_main)),
    annotation_method = "SingleR + cluster markers",
    stringsAsFactors = FALSE
  )
  write.csv(annotation_summary, file.path(data_dir, "annotation_summary_stats.csv"), row.names = FALSE)
}

# Clustering summary
summary_stats <- data.frame(
  Metric = c("Total Cells", "Total Genes", "Species", "Final Resolution", "Number of Clusters", 
             "Cell Types Identified", "UMAP Neighbors", "Method", "Status"),
  Value = c(ncol(seu), nrow(seu), species, max(resolutions), length(unique(seu$seurat_clusters)),
            length(unique(seu$cell_type)), umap_neighbors, "Seurat + SingleR", "Completed")
)
write.csv(summary_stats, file.path(data_dir, "clustering_summary.csv"), row.names = FALSE)

# JSON summary for automated processing
if (requireNamespace("jsonlite", quietly = TRUE)) {
  json_summary <- list(
    status = "completed",
    total_cells = ncol(seu),
    total_genes = nrow(seu), 
    num_clusters = length(unique(seu$seurat_clusters)),
    num_cell_types = length(unique(seu$cell_type)),
    species = species,
    final_resolution = max(resolutions),
    umap_neighbors = umap_neighbors,
    annotation_method = "improved_SingleR_plus_markers"
  )
  writeLines(jsonlite::toJSON(json_summary, pretty = TRUE, auto_unbox = TRUE),
             file.path(data_dir, "analysis_summary.json"))
}

# Set final identity to cell type names for downstream analysis
if ("data_driven_celltype" %in% colnames(seu@meta.data)) {
  Idents(seu) <- seu$data_driven_celltype
  log("Set active identity to data_driven_celltype for downstream analysis")
} else if ("cell_type" %in% colnames(seu@meta.data)) {
  Idents(seu) <- seu$cell_type
  log("Set active identity to cell_type for downstream analysis")
}

log("==================================================")
log("=== Pipeline completed successfully ===")
log("==================================================")
log(sprintf("Final results: %d cells, %d genes, %d clusters, %d cell types",
            ncol(seu), nrow(seu), length(unique(seu$seurat_clusters)), length(unique(seu$cell_type))))
log("All outputs saved to:", module_dir)
log("Completed successfully.")
cat("Clustering and annotation analysis completed successfully!\n")

# End performance tracking
end_performance_tracking(
  perf_tracker,
  success = TRUE,
  additional_metrics = list(
    n_cells = ncol(seu),
    n_genes = nrow(seu),
    n_clusters = length(unique(seu$seurat_clusters)),
    n_cell_types = length(unique(seu$cell_type))
  )
)

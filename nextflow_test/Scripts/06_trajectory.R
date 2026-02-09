#!/usr/bin/env Rscript
#===============================================================================
# 06_trajectory_inference.R
# Trajectory inference using Slingshot and Monocle3
#===============================================================================

################################################################################
# SETUP
################################################################################
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(RColorBrewer)
})

safe_lib <- function(pkgs){ for (p in pkgs) suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE))) }
safe_lib(c("plotly", "htmlwidgets", "slingshot", "SingleCellExperiment", "tradeSeq", "mgcv", "pheatmap"))

# Get script directory using bootstrap utility
script_dir <- tryCatch({
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    dirname(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    getwd()
  }
}, error = function(e) getwd())

source(file.path(script_dir, "utils", "bootstrap.R"))
script_dir <- get_script_dir()
source(file.path(script_dir, "utils", "common.R"))
source(file.path(script_dir, "utils", "logging.R"))
source(file.path(script_dir, "utils", "plotting.R"))
source(file.path(script_dir, "utils", "colors.R"))
source(file.path(script_dir, "utils", "performance.R"))
source(file.path(script_dir, "utils", "06_trajectory_plots.R"))
source(file.path(script_dir, "utils", "pvalue_helpers.R"))

# Font configuration (override defaults)
CFG_BASE_SIZE <- 18       # Base font size (was 16)
CFG_TITLE_SIZE <- 24      # Title font size (was 20)
CFG_AXIS_TITLE_SIZE <- 20 # Axis title size (was 18)

set.seed(42)

################################################################################
# HELPER FUNCTIONS
################################################################################

# Get preferred dimensional reduction for trajectory inference
# Priority: lmds > pca (with optimal dims) > umap > first available
# 
# @param seu Seurat object with dimensional reductions
# @param analysis_params Optional list with pc_dims from clustering analysis
# @return Matrix of cell embeddings (cells x dimensions)
get_preferred_dimred <- function(seu, analysis_params = NULL) {
  if ("lmds" %in% names(seu@reductions)) {
    return(Embeddings(seu, "lmds"))
  }
  if ("pca" %in% names(seu@reductions)) {
    pdims <- if (!is.null(analysis_params$pc_dims)) analysis_params$pc_dims else 1:2
    pdims <- pdims[1:min(2, length(pdims))]
    return(Embeddings(seu, "pca")[, pdims, drop = FALSE])
  }
  if ("umap" %in% names(seu@reductions)) {
    return(Embeddings(seu, "umap")[, 1:2, drop = FALSE])
  }
  
  # Fallback to first available reduction
  if (length(seu@reductions) > 0) {
    rd_name <- names(seu@reductions)[1]
    return(Embeddings(seu, rd_name)[, 1:2, drop = FALSE])
  }
  
  stop("No dimensional reduction found in Seurat object.")
}

# Get preferred clustering from Seurat metadata
# Priority: spliced_snn_res.0.25 > seurat_clusters > active ident
# 
# @param seu Seurat object with cluster assignments
# @return Factor vector of cluster assignments for each cell
get_preferred_clustering <- function(seu) {
  if ("spliced_snn_res.0.25" %in% colnames(seu@meta.data)) return(seu@meta.data$spliced_snn_res.0.25)
  if ("seurat_clusters" %in% colnames(seu@meta.data)) return(seu@meta.data$seurat_clusters)
  return(as.factor(Idents(seu)))
}

################################################################################
# ARGUMENT PARSING
################################################################################

option_list <- list(
  make_option("--seurat", type="character", default="clustered_annotated_seurat.rds", 
              help="Clustered and annotated Seurat object"),
  make_option("--out", type="character", default="output/",
              help="Output path")
)
opt <- parse_args(OptionParser(option_list = option_list))

base_dir <- file.path(opt$out)
dirs <- setup_analysis_directories(base_dir, "07_trajectory")
data_dir <- dirs$data_dir
plots_dir <- dirs$plots_dir
dash_dir <- dirs$dash_dir
log_dir <- file.path(base_dir, "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(log_dir, "trajectory_analysis.log")
setup_logging(log_file)

# Start performance tracking
perf_tracker <- start_performance_tracking("07_trajectory", base_dir)

start_time <- Sys.time()
cat("=== Trajectory Analysis Started ===", "\n", file = log_file, append = FALSE)
cat("Start time:", as.character(start_time), "\n", file = log_file, append = TRUE)
cat("Base directory:", base_dir, "\n", file = log_file, append = TRUE)

cat("Loading Seurat object\n", file = log_file, append = TRUE)
seu <- readRDS(opt$seurat)
cat("Loaded Seurat object from:", opt$seurat, "\n")
cat("Seurat object:", ncol(seu), "cells,", nrow(seu), "genes\n", file = log_file, append = TRUE)

annotations_file <- file.path(dirname(opt$out), "cluster_annotations.csv")
if (file.exists(annotations_file)) {
  cluster_annotations <- read.csv(annotations_file)
  seu$cell_type <- cluster_annotations$annotation[match(seu$seurat_clusters, cluster_annotations$cluster)]
  cat("Loaded cluster annotations\n", file = log_file, append = TRUE)
} else {
  if (!"cell_type" %in% colnames(seu@meta.data)) {
    seu$cell_type <- paste0("Cluster_", seu$seurat_clusters)
  }
  cat("Using existing annotations\n", file = log_file, append = TRUE)
}

################################################################################
# LOAD UPSTREAM RESULTS
################################################################################

cat("Loading analysis parameters from previous steps...\n")
cat("Loading upstream analysis parameters\n", file = log_file, append = TRUE)

# Load analysis parameters from Step 02 (clustering)
analysis_params <- NULL
param_paths <- c(
  file.path(dirname(opt$seurat), "data", "analysis_parameters.rds"),
  file.path(dirname(dirname(opt$seurat)), "03_clustering", "data", "analysis_parameters.rds"),
  file.path(opt$out, "..", "03_clustering", "data", "analysis_parameters.rds")
)
for (param_path in param_paths) {
  if (file.exists(param_path)) {
    analysis_params <- readRDS(param_path)
    cat("Loaded clustering parameters:", length(analysis_params$pc_dims), "PCs\n")
    cat("  - PCs:", length(analysis_params$pc_dims), "\n", file = log_file, append = TRUE)
    break
  }
}

if (is.null(analysis_params)) {
  analysis_params <- list(
    pc_dims = 1:30,
    num_pcs = 30,
    num_clusters = length(unique(seu$seurat_clusters)),
    jackstraw_used = FALSE
  )
  cat("Using default: 30 PCs\n", file = log_file, append = TRUE)
}

de_genes_all <- c()
da_dir_candidates <- c(
  file.path(dirname(dirname(opt$seurat)), "04_differential"),
  file.path(dirname(opt$seurat), "..", "04_differential"),
  file.path(opt$out, "..", "04_differential")
)

da_dir <- NULL
for (candidate in da_dir_candidates) {
  if (dir.exists(candidate)) {
    da_dir <- candidate
    cat("Found differential analysis results at:", da_dir, "\n")
    cat("Found DA results:", da_dir, "\n", file = log_file, append = TRUE)
    break
  }
}

if (!is.null(da_dir)) {
  ds_file <- file.path(da_dir, "data", "DS_seurat_results_all_cells.csv")
  if (file.exists(ds_file)) {
    ds_results <- read.csv(ds_file, stringsAsFactors = FALSE)
    de_genes_all <- ds_results %>%
      filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
      pull(gene) %>%
      unique()
    cat("Loaded", length(de_genes_all), "condition-responsive DE genes from Step 03\n")
    cat("  - DE genes loaded:", length(de_genes_all), "\n", file = log_file, append = TRUE)
  }
}

# Load Ligand-Receptor pairs from Step 04 (cell communication) - BEST PRACTICE
lr_genes_all <- c()
lr_pairs <- data.frame()
comm_dir_candidates <- c(
  file.path(dirname(dirname(opt$seurat)), "05_cellcomm"),
  file.path(dirname(opt$seurat), "..", "05_cellcomm"),
  file.path(opt$out, "..", "05_cellcomm")
)

comm_dir <- NULL
for (candidate in comm_dir_candidates) {
  if (dir.exists(candidate)) {
    comm_dir <- candidate
    cat("Found cell communication results at:", comm_dir, "\n")
    cat("Found communication results:", comm_dir, "\n", file = log_file, append = TRUE)
    break
  }
}

if (!is.null(comm_dir)) {
  # Load NichNet ligand activities (all cell types)
  nichnet_files <- list.files(file.path(comm_dir, "data"), 
                              pattern = "^nichnet_ligand_activities_.*\\.csv$", 
                              full.names = TRUE)
  
  if (length(nichnet_files) > 0) {
    lr_genes_list <- list()
    for (nf in nichnet_files) {
      tryCatch({
        nichnet_data <- read.csv(nf, stringsAsFactors = FALSE)
        if ("aupr_corrected" %in% colnames(nichnet_data)) {
          top_ligands <- nichnet_data %>%
            filter(aupr_corrected > 0.01) %>%
            pull(test_ligand) %>%
            unique()
          lr_genes_list <- c(lr_genes_list, list(top_ligands))
        }
      }, error = function(e) {
        cat("Warning: Failed to load", basename(nf), "\n")
      })
    }
    lr_genes_all <- unique(unlist(lr_genes_list))
    cat("Loaded", length(lr_genes_all), "ligand genes from Step 04 (NichNet)\n")
    cat("  - L-R genes loaded:", length(lr_genes_all), "\n", file = log_file, append = TRUE)
  }
  
}

cond_field <- NULL
conditions <- NULL
for (candidate in c("condition", "treatment", "group")) {
  if (candidate %in% colnames(seu@meta.data)) {
    cond_field <- candidate
    conditions <- unique(seu@meta.data[[cond_field]])
    cat("Condition field:", cond_field, "\n", file = log_file, append = TRUE)
    cat("Found", length(conditions), "conditions:", paste(conditions, collapse=", "), "\n")
    break
  }
}

cat("Integration complete: cells=", ncol(seu), ", DE genes=", length(de_genes_all), ", L-R genes=", length(lr_genes_all), "\n")

trajectory_results <- list()
trajectory_method_used <- "none"

# If conditions found, analyze separately
if (!is.null(cond_field) && !is.null(conditions) && length(conditions) > 1) {
  cat("Trajectory analysis will be performed SEPARATELY for each condition:\n")
  for (cond in conditions) {
    cat(sprintf("  - %s\n", cond))
  }
} else if (!is.null(cond_field)) {
  cat("Only 1 condition found - analyzing all cells together\n")
}

################################################################################
# SLINGSHOT TRAJECTORY ANALYSIS
# 
# Trajectory inference using annotated and integrated expression values.
# The slingshot function is applied with:
#   - clusterLabels: cluster assignments from Seurat clustering
#   - reducedDim: PCA embedding to preserve global distance structure
#   - Random seed: 42 for reproducibility
#   - Multiple lineages allowed with automatic refinement
#   - Lineage assignment via simultaneous principal curves
#   - Default smoothing parameters selected by cross-validation
################################################################################

if ("slingshot" %in% rownames(installed.packages())) {
  cat("Running Slingshot analysis\n", file = log_file, append = TRUE)
  
  tryCatch({
    safe_lib("slingshot")
    safe_lib("SingleCellExperiment")
    
    # Determine which datasets to analyze
    datasets_to_analyze <- list()
    if (!is.null(cond_field) && !is.null(conditions) && length(conditions) > 1) {
      # Analyze each condition separately
      for (cond in conditions) {
        cells_in_cond <- seu@meta.data[[cond_field]] == cond
        n_cells_cond <- sum(cells_in_cond, na.rm = TRUE)
        if (n_cells_cond > 100) {  # Only if enough cells
          # Use subset() to properly preserve reductions
          subset_expr <- paste0(cond_field, " == '", cond, "'")
          seu_cond <- subset(seu, subset = eval(parse(text = subset_expr)))
          datasets_to_analyze[[as.character(cond)]] <- list(seu=seu_cond, name=as.character(cond))
          cat(sprintf("  - %s: %d cells\n", cond, ncol(seu_cond)), file = log_file, append = TRUE)
        }
      }
    } else {
      # Analyze all cells together
      datasets_to_analyze[["all_cells"]] <- list(seu=seu, name="all_cells")
    }
    
    # Run trajectory analysis for each dataset
    for (dataset_name in names(datasets_to_analyze)) {
      cat(sprintf("\nAnalyzing trajectories for: %s\n", dataset_name), file = log_file, append = TRUE)
      
      tryCatch({
      seu_sub <- datasets_to_analyze[[dataset_name]]$seu
      sce <- as.SingleCellExperiment(seu_sub)
    
    # Check available reductions
    available_reductions <- names(seu_sub@reductions)
    cat("  Available reductions:", paste(available_reductions, collapse=", "), "\n", 
        file = log_file, append = TRUE)
    
    if (length(available_reductions) == 0) {
      stop("No dimensionality reductions found in Seurat object")
    }
    
    # Use PCA embedding to preserve global distance structure
    if ("pca" %in% available_reductions) {
      pca_obj <- seu_sub@reductions$pca
      if (is.null(pca_obj)) {
        stop("PCA reduction is NULL")
      }
      pca_mat <- Embeddings(seu_sub, "pca")
      pc_dims <- analysis_params$pc_dims
      pc_dims <- pc_dims[pc_dims >= 1 & pc_dims <= ncol(pca_mat)]
      if (length(pc_dims) == 0) pc_dims <- 1:min(2, ncol(pca_mat))
      rd <- pca_mat[, pc_dims, drop = FALSE]
      rd_name <- "pca"
      cat("  Using PCA:", ncol(rd), "dimensions\n", file = log_file, append = TRUE)
    } else if ("umap" %in% available_reductions) {
      umap_obj <- seu_sub@reductions$umap
      if (is.null(umap_obj)) {
        stop("UMAP reduction is NULL")
      }
      rd <- Embeddings(seu_sub, "umap")
      rd_name <- "umap"
      cat("  Using UMAP\n", file = log_file, append = TRUE)
    } else {
      red_obj <- seu_sub@reductions[[available_reductions[1]]]
      if (is.null(red_obj)) {
        stop(paste("First reduction", available_reductions[1], "is NULL"))
      }
      rd <- Embeddings(seu_sub, available_reductions[1])
      rd_name <- available_reductions[1]
      cat("  Using", rd_name, "\n", file = log_file, append = TRUE)
    }
    
    # Set cluster assignments from Seurat clustering
    colData(sce)$cluster <- Idents(seu_sub)
    colData(sce)$cell_type <- seu_sub$cell_type
    
    if (!is.null(cond_field)) {
      colData(sce)$condition <- seu_sub@meta.data[[cond_field]]
    }
    
    # Configure parameters
    n_cells <- ncol(sce)
    approx_points_value <- min(150, n_cells)
    
    # Allow multiple lineages with omega parameter
    use_omega <- FALSE
    n_clusters <- length(unique(colData(sce)$cluster))
    if (n_clusters > 10) {
      use_omega <- TRUE
      cat("  - Using omega=TRUE for", n_clusters, "clusters\n")
    }
    
      # Set random seed for reproducibility
      set.seed(42)
      
      # Run slingshot with simultaneous principal curves
      # Multiple lineages allowed, refined automatically
      # Smoothing parameters selected via cross-validation (default)
      run_slingshot <- function(reduced_dim) {
        if (use_omega) {
          slingshot::slingshot(sce, clusterLabels = 'cluster', reducedDim = reduced_dim,
                               approx_points = approx_points_value, omega = TRUE)
        } else {
          slingshot::slingshot(sce, clusterLabels = 'cluster', reducedDim = reduced_dim,
                               approx_points = approx_points_value)
        }
      }

      sce <- tryCatch({
        run_slingshot(rd)
      }, error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("computationally singular", msg, ignore.case = TRUE)) {
          cat("  - Slingshot singular matrix detected. Retrying with fewer PCs/UMAP\n", file = log_file, append = TRUE)
          
          # Try 10 PCs first
          if (rd_name == "pca" && ncol(rd) > 10) {
            tryCatch({
              rd_fallback <- rd[, 1:10, drop = FALSE]
              cat("  - Trying with 10 PCs\n", file = log_file, append = TRUE)
              return(run_slingshot(rd_fallback))
            }, error = function(e2) {
              cat("  - 10 PCs also failed, trying UMAP\n", file = log_file, append = TRUE)
            })
          }
          
          # Fallback to UMAP
          if ("umap" %in% names(seu_sub@reductions)) {
            tryCatch({
              rd_fallback <- Embeddings(seu_sub, "umap")
              rd_name <<- "umap"
              rd <<- rd_fallback
              cat("  - Using UMAP embedding\n", file = log_file, append = TRUE)
              return(run_slingshot(rd_fallback))
            }, error = function(e3) {
              cat("  - UMAP also failed:", conditionMessage(e3), "\n", file = log_file, append = TRUE)
            })
          }
        }
        stop(e)
      })
      
      # Save with dataset name
      sce_file <- file.path(data_dir, sprintf("slingshot_sce_%s.rds", dataset_name))
      saveRDS(sce, sce_file)
      cat(sprintf("  Saved SCE for %s\n", dataset_name), file = log_file, append = TRUE)
      
      curves <- slingshot::slingCurves(sce)
      pseudotime_curves <- slingshot::slingPseudotime(sce)
      curve_weights <- slingshot::slingCurveWeights(sce)
      
      # Save trajectory data
      trajectory_data <- data.frame(
        cell_id = colnames(sce),
        cluster = colData(sce)$cluster,
        cell_type = colData(sce)$cell_type,
        pseudotime_curves,
        curve_weights
      )
      
      # Add coordinates
      trajectory_data[[paste0(rd_name, "_1")]] <- rd[, 1]
    trajectory_data[[paste0(rd_name, "_2")]] <- rd[, 2]
      
      # Save with dataset name
      traj_file <- file.path(data_dir, sprintf("slingshot_trajectory_data_%s.csv", dataset_name))
      write.csv(trajectory_data, traj_file, row.names = FALSE)
      
      trajectory_results[[dataset_name]] <- list(
        sce = sce,
        curves = curves,
        pseudotime = pseudotime_curves,
        weights = curve_weights,
        data = trajectory_data,
        name = dataset_name
      )
      
      cat(sprintf("  %s: %d trajectories completed\n", dataset_name, length(curves)), file = log_file, append = TRUE)
      
      # Also run getLineages for visualization/demonstration
      tryCatch({
        # Determine dimensionality reduction for getLineages
        dimred_for_lineages <- NULL
        if ("lmds" %in% names(seu_sub@reductions)) {
          dimred_for_lineages <- Embeddings(seu_sub, "lmds")
        } else if ("pca" %in% names(seu_sub@reductions)) {
          pdims <- if (!is.null(analysis_params$pc_dims)) analysis_params$pc_dims else 1:2
          pdims <- pdims[1:min(2, length(pdims))]
          dimred_for_lineages <- Embeddings(seu_sub, "pca")[, pdims, drop = FALSE]
        } else if ("umap" %in% names(seu_sub@reductions)) {
          dimred_for_lineages <- Embeddings(seu_sub, "umap")[, 1:2, drop = FALSE]
        } else {
          dimred_for_lineages <- rd[, 1:min(2, ncol(rd)), drop = FALSE]
        }

        # Prioritize cell_type for labeling if available
        clustering_for_lineages <- NULL
        if ("cell_type" %in% colnames(seu_sub@meta.data)) {
          clustering_for_lineages <- as.factor(seu_sub$cell_type)
        } else if ("spliced_snn_res.0.25" %in% colnames(seu_sub@meta.data)) {
          clustering_for_lineages <- seu_sub@meta.data$spliced_snn_res.0.25
        } else if ("seurat_clusters" %in% colnames(seu_sub@meta.data)) {
          clustering_for_lineages <- seu_sub$seurat_clusters
        } else {
          clustering_for_lineages <- as.factor(Idents(seu_sub))
        }

        set.seed(1)
        lineages <- slingshot::getLineages(dimred_for_lineages, clustering_for_lineages)
        
        # Save with dataset name
        safe_cond <- gsub("[^A-Za-z0-9_-]+", "_", dataset_name)
        saveRDS(lineages, file.path(data_dir, sprintf("slingshot_lineages_getLineages_%s.rds", safe_cond)))

        # Create plot for each condition
        create_lineages_plot(dimred_for_lineages, lineages, clustering_for_lineages, plots_dir, cat, condition = dataset_name)
        
        cat(sprintf("  - getLineages visualization saved for %s\n", dataset_name), file = log_file, append = TRUE)
      }, error = function(e) {
        cat(sprintf("  - getLineages visualization failed for %s: %s\n", dataset_name, conditionMessage(e)), file = log_file, append = TRUE)
      })
      
      }, error = function(e) {
        cat(sprintf("  Warning: Analysis failed for %s: %s\n", dataset_name, conditionMessage(e)), file = log_file, append = TRUE)
        cat(sprintf("  Skipping %s and continuing with next condition\n", dataset_name), file = log_file, append = TRUE)
      })
    }  # End dataset loop
    
    if (length(trajectory_results) > 0) {
      trajectory_method_used <- "slingshot"
    }
    
    if (length(trajectory_results) > 1) {
      cat(sprintf("Slingshot analysis completed for %d conditions\n", length(trajectory_results)), file = log_file, append = TRUE)
    } else {
      cat("Slingshot analysis completed\n", file = log_file, append = TRUE)
    }
    
    ################################################################################
    # COMBINED TRAJECTORY ANALYSIS (ALL CONDITIONS TOGETHER) + TRADESEQ
    # Required for statistical comparison between conditions using tradeSeq
    ################################################################################
    
    if (!is.null(cond_field) && !is.null(conditions) && length(conditions) > 1 && 
        "tradeSeq" %in% rownames(installed.packages())) {
      
      cat("\n=== Combined Trajectory Analysis (All Conditions) ===\n", file = log_file, append = TRUE)
      cat("Running Slingshot on combined dataset for tradeSeq differential testing...\n")
      
      tryCatch({
        safe_lib("tradeSeq")
        
        # Use all cells together (not split by condition)
        sce_combined <- as.SingleCellExperiment(seu)
        
        # Use PCA embedding
        if ("pca" %in% names(seu@reductions)) {
          pca_mat <- Embeddings(seu, "pca")
          pc_dims <- analysis_params$pc_dims
          pc_dims <- pc_dims[pc_dims >= 1 & pc_dims <= ncol(pca_mat)]
          if (length(pc_dims) == 0) pc_dims <- 1:min(30, ncol(pca_mat))
          rd_combined <- pca_mat[, pc_dims, drop = FALSE]
          cat("Using PCA:", ncol(rd_combined), "dimensions for combined analysis\n", file = log_file, append = TRUE)
        } else if ("umap" %in% names(seu@reductions)) {
          rd_combined <- Embeddings(seu, "umap")
        } else {
          stop("No PCA or UMAP reduction found")
        }
        
        # Set cluster and condition labels
        colData(sce_combined)$cluster <- Idents(seu)
        colData(sce_combined)$cell_type <- seu$cell_type
        colData(sce_combined)$condition <- seu@meta.data[[cond_field]]
        
        # Configure parameters
        n_cells_combined <- ncol(sce_combined)
        approx_points_combined <- min(150, n_cells_combined)
        use_omega_combined <- length(unique(colData(sce_combined)$cluster)) > 10
        
        set.seed(42)
        
        # Run slingshot on combined data
        cat("Fitting trajectories on combined dataset...\n", file = log_file, append = TRUE)
        if (use_omega_combined) {
          sce_combined <- slingshot::slingshot(sce_combined, clusterLabels = 'cluster', 
                                              reducedDim = rd_combined,
                                              approx_points = approx_points_combined, 
                                              omega = TRUE)
        } else {
          sce_combined <- slingshot::slingshot(sce_combined, clusterLabels = 'cluster', 
                                              reducedDim = rd_combined,
                                              approx_points = approx_points_combined)
        }
        
        # Save combined SCE
        saveRDS(sce_combined, file.path(data_dir, "slingshot_sce_combined_all_conditions.rds"))
        cat("Saved combined SCE with", length(slingshot::slingCurves(sce_combined)), "trajectories\n", file = log_file, append = TRUE)
        
        # Extract pseudotime and weights
        pseudotime_combined <- slingshot::slingPseudotime(sce_combined)
        weights_combined <- slingshot::slingCurveWeights(sce_combined)
        
        # Prepare trajectory data with conditions
        trajectory_data_combined <- data.frame(
          cell_id = colnames(sce_combined),
          cluster = colData(sce_combined)$cluster,
          cell_type = colData(sce_combined)$cell_type,
          condition = colData(sce_combined)$condition,
          pseudotime_combined,
          weights_combined
        )
        
        write.csv(trajectory_data_combined, 
                 file.path(data_dir, "slingshot_trajectory_data_combined_all_conditions.csv"), 
                 row.names = FALSE)
        
        ################################################################################
        # TRADESEQ: DIFFERENTIAL TRAJECTORY TESTING
        ################################################################################
        
        cat("\n=== tradeSeq Differential Analysis ===\n", file = log_file, append = TRUE)
        cat("Fitting GAM models for differential testing between conditions...\n")
        
        # Get expression data (use normalized counts)
        counts_matrix <- as.matrix(GetAssayData(seu, assay = "RNA", slot = "data"))
        
        # Subset to genes with sufficient expression (speed up analysis)
        gene_means <- Matrix::rowMeans(counts_matrix)
        expressed_genes <- names(gene_means[gene_means > 0.01])
        
        # Further subset to HVG or DE genes if available
        if (length(de_genes_all) > 0) {
          # Use DE genes from differential analysis
          test_genes <- intersect(expressed_genes, de_genes_all)
          cat("Using", length(test_genes), "DE genes for tradeSeq testing\n", file = log_file, append = TRUE)
        } else if (length(VariableFeatures(seu)) > 0) {
          # Use highly variable genes
          test_genes <- intersect(expressed_genes, VariableFeatures(seu)[1:min(2000, length(VariableFeatures(seu)))])
          cat("Using", length(test_genes), "HVG genes for tradeSeq testing\n", file = log_file, append = TRUE)
        } else {
          # Use top expressed genes
          test_genes <- names(sort(gene_means, decreasing = TRUE)[1:min(2000, length(expressed_genes))])
          cat("Using", length(test_genes), "highly expressed genes for tradeSeq testing\n", file = log_file, append = TRUE)
        }
        
        counts_subset <- counts_matrix[test_genes, , drop = FALSE]
        
        # Fit GAM models with condition as covariate
        cat("Fitting GAM models (this may take several minutes)...\n", file = log_file, append = TRUE)
        
        gam_models <- tryCatch({
          tradeSeq::fitGAM(
            counts = counts_subset,
            pseudotime = pseudotime_combined,
            cellWeights = weights_combined,
            conditions = factor(colData(sce_combined)$condition),
            nknots = 6,
            verbose = TRUE
          )
        }, error = function(e) {
          cat("Error fitting GAM models:", conditionMessage(e), "\n", file = log_file, append = TRUE)
          return(NULL)
        })
        
        if (!is.null(gam_models)) {
          saveRDS(gam_models, file.path(data_dir, "tradeseq_gam_models.rds"))
          cat("GAM models fitted successfully\n", file = log_file, append = TRUE)
          
          # Test 1: conditionTest - genes with different patterns between conditions
          cat("Running conditionTest...\n", file = log_file, append = TRUE)
          condition_test_results <- tryCatch({
            tradeSeq::conditionTest(gam_models)
          }, error = function(e) {
            cat("conditionTest failed:", conditionMessage(e), "\n", file = log_file, append = TRUE)
            return(NULL)
          })
          
          if (!is.null(condition_test_results)) {
            condition_test_df <- as.data.frame(condition_test_results)
            condition_test_df$gene <- rownames(condition_test_df)
            condition_test_df <- condition_test_df[order(condition_test_df$pvalue), ]
            
            write.csv(condition_test_df, 
                     file.path(data_dir, "tradeseq_conditionTest_results.csv"), 
                     row.names = FALSE)
            
            n_sig <- sum(condition_test_df$pvalue < 0.05, na.rm = TRUE)
            cat(sprintf("conditionTest: %d genes with p < 0.05\n", n_sig), file = log_file, append = TRUE)
          }
          
          # Test 2: patternTest - genes with different expression patterns
          cat("Running patternTest...\n", file = log_file, append = TRUE)
          pattern_test_results <- tryCatch({
            tradeSeq::patternTest(gam_models)
          }, error = function(e) {
            cat("patternTest failed:", conditionMessage(e), "\n", file = log_file, append = TRUE)
            return(NULL)
          })
          
          if (!is.null(pattern_test_results)) {
            pattern_test_df <- as.data.frame(pattern_test_results)
            pattern_test_df$gene <- rownames(pattern_test_df)
            pattern_test_df <- pattern_test_df[order(pattern_test_df$pvalue), ]
            
            write.csv(pattern_test_df, 
                     file.path(data_dir, "tradeseq_patternTest_results.csv"), 
                     row.names = FALSE)
            
            n_sig <- sum(pattern_test_df$pvalue < 0.05, na.rm = TRUE)
            cat(sprintf("patternTest: %d genes with p < 0.05\n", n_sig), file = log_file, append = TRUE)
          }
          
          # Test 3: diffEndTest - genes with different endpoints between conditions
          cat("Running diffEndTest...\n", file = log_file, append = TRUE)
          diffend_test_results <- tryCatch({
            tradeSeq::diffEndTest(gam_models)
          }, error = function(e) {
            cat("diffEndTest failed:", conditionMessage(e), "\n", file = log_file, append = TRUE)
            return(NULL)
          })
          
          if (!is.null(diffend_test_results)) {
            diffend_test_df <- as.data.frame(diffend_test_results)
            diffend_test_df$gene <- rownames(diffend_test_df)
            diffend_test_df <- diffend_test_df[order(diffend_test_df$pvalue), ]
            
            write.csv(diffend_test_df, 
                     file.path(data_dir, "tradeseq_diffEndTest_results.csv"), 
                     row.names = FALSE)
            
            n_sig <- sum(diffend_test_df$pvalue < 0.05, na.rm = TRUE)
            cat(sprintf("diffEndTest: %d genes with p < 0.05\n", n_sig), file = log_file, append = TRUE)
          }
          
          # Store results in trajectory_results
          trajectory_results[["tradeseq"]] <- list(
            models = gam_models,
            condition_test = condition_test_results,
            pattern_test = pattern_test_results,
            diffend_test = diffend_test_results
          )
          
          cat("tradeSeq differential testing completed successfully\n", file = log_file, append = TRUE)
          
        } else {
          cat("GAM model fitting failed - skipping differential tests\n", file = log_file, append = TRUE)
        }
        
        # Store combined trajectory in results
        trajectory_results[["combined"]] <- list(
          sce = sce_combined,
          pseudotime = pseudotime_combined,
          weights = weights_combined,
          data = trajectory_data_combined
        )
        
      }, error = function(e) {
        cat("Combined trajectory analysis failed:", conditionMessage(e), "\n", file = log_file, append = TRUE)
        cat("Continuing with per-condition results only\n", file = log_file, append = TRUE)
      })
      
    } else if (!("tradeSeq" %in% rownames(installed.packages()))) {
      cat("tradeSeq package not installed - skipping differential testing\n", file = log_file, append = TRUE)
      cat("Install with: BiocManager::install('tradeSeq')\n")
    }
    # Demonstrate predict() functionality using combined trajectory if available
    tryCatch({
      # Use combined trajectory if available, otherwise use first condition
      sce_for_demo <- if ("combined" %in% names(trajectory_results)) {
        trajectory_results[["combined"]]$sce
      } else if (length(trajectory_results) > 0) {
        trajectory_results[[1]]$sce
      } else {
        NULL
      }
      
      if (!is.null(sce_for_demo)) {
        # Get the reducedDim used for this SCE
        rd_demo <- NULL
        rd_name_demo <- "pca"
        
        if ("combined" %in% names(trajectory_results)) {
          if (exists("rd_combined")) {
            rd_demo <- rd_combined
            rd_name_demo <- "pca"
          }
        } else if (length(trajectory_results) > 0) {
          # Try to get from the first trajectory result
          first_result <- trajectory_results[[1]]
          if (!is.null(first_result$sce)) {
            available_dims <- reducedDimNames(first_result$sce)
            if (length(available_dims) > 0) {
              rd_demo <- reducedDim(first_result$sce, available_dims[1])
              rd_name_demo <- available_dims[1]
            }
          }
        }
        
        if (is.null(rd_demo)) {
          stop("No reduced dimensions available for projection demo")
        }
        
        # Simulate "new" cells by adding noise to existing coordinates
        newPCA <- rd_demo + matrix(rnorm(length(rd_demo), sd = 0.5), nrow = nrow(rd_demo), ncol = ncol(rd_demo))
        
        pto_original <- slingshot::SlingshotDataSet(sce_for_demo)
        newPTO <- slingshot::predict(pto_original, newPCA)
        
        projection_data <- data.frame(
          cell_id = paste0("projected_", 1:nrow(newPCA)),
          slingPseudotime(newPTO),
          slingCurveWeights(newPTO)
        )
        write.csv(projection_data, file.path(data_dir, "slingshot_projected_cells.csv"), row.names = FALSE)
        
        # Create visualization
        pseudotime_demo <- slingshot::slingPseudotime(sce_for_demo)
        create_projection_demo_plot(rd_demo, rd_name_demo, pto_original, newPCA, newPTO, pseudotime_demo, plots_dir, cat)
        
        cat("Projection demo completed\n", file = log_file, append = TRUE)
      }
    }, error = function(e) {
      cat("  - Warning: Cell projection demo failed:", conditionMessage(e), "\n", file = log_file, append = TRUE)
    })
    
  }, error = function(e) {
    message("Slingshot analysis failed: ", conditionMessage(e))
    cat("Slingshot failed: ", conditionMessage(e), "\n", file = log_file, append = TRUE)
  })
}

################################################################################
# MONOCLE3 TRAJECTORY ANALYSIS
################################################################################

if ("monocle3" %in% rownames(installed.packages())) {
  cat("Starting Monocle3 analysis\n", file = log_file, append = TRUE)
  
  tryCatch({
    safe_lib("monocle3")
    
    gene_metadata <- data.frame(gene_short_name = rownames(seu))
    rownames(gene_metadata) <- rownames(seu)
    
    cds <- monocle3::new_cell_data_set(
      expression_data = GetAssayData(seu, assay = "RNA", slot = "counts"),
      cell_metadata = seu@meta.data,
      gene_metadata = gene_metadata
    )
    
    if ("pca" %in% names(seu@reductions)) {
      pca_mat <- Embeddings(seu, "pca")
      pc_dims <- analysis_params$pc_dims
      pc_dims <- pc_dims[pc_dims >= 1 & pc_dims <= ncol(pca_mat)]
      if (length(pc_dims) == 0) pc_dims <- 1:min(30, ncol(pca_mat))
      pca_embeddings <- pca_mat[, pc_dims, drop = FALSE]
      cds@int_colData@listData$reducedDims@listData$PCA <- pca_embeddings
      cat("Reused PCA:", ncol(pca_embeddings), "dimensions\n", file = log_file, append = TRUE)
      
      cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
    } else {
      num_dims_monocle <- analysis_params$num_pcs
      cds <- monocle3::preprocess_cds(cds, num_dim = num_dims_monocle)
      cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")
    }
    
    cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")
    cds <- monocle3::learn_graph(cds, use_partition = FALSE)
    cds <- monocle3::order_cells(cds)
    
    saveRDS(cds, file.path(data_dir, "monocle3_cds.rds"))
    
    pseudotime_data <- data.frame(
      cell_id = colnames(cds),
      pseudotime = cds@principal_graph_aux@listData$UMAP$pseudotime,
      cluster = cds@clusters@listData$UMAP$clusters,
      partition = cds@clusters@listData$UMAP$partitions
    )
    
    umap_coords <- cds@int_colData@listData$reducedDims@listData$UMAP
    pseudotime_data$UMAP_1 <- umap_coords[, 1]
    pseudotime_data$UMAP_2 <- umap_coords[, 2]
    
    write.csv(pseudotime_data, file.path(data_dir, "monocle3_trajectory_data.csv"), row.names = FALSE)
    
    trajectory_results[["monocle3"]] <- list(
      cds = cds,
      pseudotime_data = pseudotime_data
    )
    
    if (trajectory_method_used == "none") {
      trajectory_method_used <- "monocle3"
    } else {
      trajectory_method_used <- paste(trajectory_method_used, "monocle3", sep = "+")
    }
    cat("Monocle3 analysis completed\n")
    cat("Monocle3 completed successfully\n", file = log_file, append = TRUE)
    
  }, error = function(e) {
    message("Monocle3 analysis failed: ", conditionMessage(e))
    cat("Monocle3 failed: ", conditionMessage(e), "\n", file = log_file, append = TRUE)
  })
}


################################################################################
# GENERATE VISUALIZATIONS
################################################################################

if (trajectory_method_used != "none" && length(trajectory_results) > 0) {
  cat("Generating visualizations\n", file = log_file, append = TRUE)
  
  # Get the trajectory data for visualization - prioritize slingshot over monocle3
  if ("slingshot" %in% names(trajectory_results)) {
    viz_data <- trajectory_results$slingshot$data
    coord_cols <- c("umap_1", "umap_2")
    if (!all(coord_cols %in% colnames(viz_data))) {
      coord_cols <- colnames(viz_data)[grepl("_1$|_2$", colnames(viz_data))][1:2]
    }
  } else if ("monocle3" %in% names(trajectory_results)) {
    viz_data <- trajectory_results$monocle3$pseudotime_data
    coord_cols <- c("UMAP_1", "UMAP_2")
  } else {
    cat("No trajectory results available for visualization\n", file = log_file, append = TRUE)
    viz_data <- NULL
  }
  
  # Skip visualization if no data
  if (is.null(viz_data)) {
    cat("Skipping visualization - no trajectory data available\n", file = log_file, append = TRUE)
  } else {
    # 4.1 Static plots
    tryCatch({
    # Trajectory overview plot
    if ("slingshot" %in% names(trajectory_results)) {
      create_slingshot_trajectory_overview(viz_data, coord_cols, trajectory_results$slingshot$curves, plots_dir, cat)
      cat("  - Trajectory overview completed\n", file = log_file, append = TRUE)
    } else {
      create_generic_trajectory_plot(viz_data, coord_cols, plots_dir, cat)
      cat("  - Generic trajectory plot completed\n", file = log_file, append = TRUE)
    }
    
    # Additional static plots
    if ("slingshot" %in% names(trajectory_results)) {
      create_pseudotime_histograms(viz_data, plots_dir, cat)
      cat("  - Pseudotime histograms completed\n", file = log_file, append = TRUE)
      create_pseudotime_boxplots(viz_data, plots_dir, cat)
      cat("  - Pseudotime boxplots completed\n", file = log_file, append = TRUE)
      create_curve_weights_heatmap(viz_data, plots_dir, cat)
      cat("  - Curve weights heatmap completed\n", file = log_file, append = TRUE)
    } else if ("pseudotime" %in% colnames(viz_data)) {
      # Run a minimal Slingshot workflow: find lineages (prefer upstream lmds),
      # optionally set a start cluster by Dpp4 expression, fit principal curves,
      # and save objects + a simple plot. Heavy downstream testing (tradeSeq, etc.)
      # is intentionally omitted here.
      if ("slingshot" %in% rownames(installed.packages())) {
        cat("Running minimal Slingshot lineage finding...\n")
        cat("Starting minimal Slingshot\n", file = log_file, append = TRUE)

        tryCatch({
          safe_lib("slingshot")

          dimred <- get_preferred_dimred(seu, analysis_params)
          clustering <- get_preferred_clustering(seu)

          # initial lineage finding
          set.seed(1)
          lineages <- slingshot::getLineages(dimred, clustering)

          # optional: pick start cluster by Dpp4/Cd26 expression if present
          start_cluster_id <- NULL
          dpp4_idx <- which(tolower(rownames(seu)) %in% tolower(c("Dpp4", "Cd26")))
          if (length(dpp4_idx) > 0) {
            counts_vec <- if ("spliced" %in% names(seu@assays)) as.numeric(seu@assays$spliced@counts[dpp4_idx[1], ]) else as.numeric(GetAssayData(seu, slot = "counts")[dpp4_idx[1], ])
            start_cluster_id <- dplyr::tibble(expression = counts_vec, cluster_id = clustering) %>%
              dplyr::group_by(cluster_id) %>%
              dplyr::summarise(expression = mean(expression), .groups = "drop") %>%
              dplyr::arrange(desc(expression)) %>%
              dplyr::pull(cluster_id) %>%
              dplyr::first() %>%
              as.character()
          }

          if (!is.null(start_cluster_id)) {
            set.seed(1)
            lineages <- slingshot::getLineages(dimred, clustering, start.clus = start_cluster_id)
          }

          curves <- slingshot::getCurves(lineages)

          saveRDS(lineages, file.path(data_dir, "slingshot_lineages.rds"))
          saveRDS(curves, file.path(data_dir, "slingshot_curves.rds"))

          create_minimal_slingshot_plot(dimred, curves, clustering, plots_dir, cat)
          
          cat("Minimal Slingshot: saved objects\n", file = log_file, append = TRUE)
          trajectory_results[["slingshot"]] <- list(lineages = lineages, curves = curves)
          trajectory_method_used <- "slingshot"

        }, error = function(e) {
          message("Minimal Slingshot failed: ", conditionMessage(e))
          cat("Minimal Slingshot failed: ", conditionMessage(e), "\n", file = log_file, append = TRUE)
        })
      }  # Close if ("slingshot" %in% rownames(installed.packages()))
    }
    
    # Pseudotime comparison by condition
    if (!is.null(cond_field)) {
      if (any(grepl("curve", colnames(viz_data)))) {
        pt_cols <- colnames(viz_data)[grepl("curve", colnames(viz_data)) & !grepl("weight", colnames(viz_data))]
        
        for (pt_col in pt_cols) {
          create_interactive_pseudotime_by_condition(viz_data, cond_field, pt_col, dash_dir, cat)
        }
        
        # Create combined pseudotime plots for all conditions
        create_pseudotime_by_condition_combined(viz_data, cond_field, plots_dir, cat)
        create_pseudotime_histograms_by_condition(viz_data, cond_field, plots_dir, cat)
      }
      
      cat("  - Saved condition-specific trajectory plots\n", file = log_file, append = TRUE)
    }
    
    # Interactive trajectory plot by cell type
    create_interactive_celltype_plot(viz_data, coord_cols, dash_dir, cat)
    cat("  - Interactive cell type plot completed\n", file = log_file, append = TRUE)
    
    # Pseudotime visualization
    if ("slingshot" %in% names(trajectory_results)) {
      create_interactive_pseudotime_slingshot(viz_data, coord_cols, dash_dir, cat)
      cat("  - Interactive pseudotime slingshot completed\n", file = log_file, append = TRUE)
      create_interactive_curve_comparison(viz_data, dash_dir, cat)
      cat("  - Interactive curve comparison completed\n", file = log_file, append = TRUE)
      create_interactive_curve_weights(viz_data, coord_cols, dash_dir, cat)
      cat("  - Interactive curve weights completed\n", file = log_file, append = TRUE)
    } else if ("pseudotime" %in% colnames(viz_data)) {
      create_interactive_pseudotime_generic(viz_data, coord_cols, dash_dir, cat)
      cat("  - Interactive pseudotime generic completed\n", file = log_file, append = TRUE)
    }
    
    # Pseudotime distributions by cell type
    create_interactive_violin_by_celltype(viz_data, trajectory_method_used, dash_dir, cat)
    cat("  - Interactive violin by cell type completed\n", file = log_file, append = TRUE)
    
    # 3D trajectory visualization
    if ("slingshot" %in% names(trajectory_results)) {
      create_interactive_3d_trajectory(viz_data, coord_cols, dash_dir, cat)
      cat("  - Interactive 3D trajectory completed\n", file = log_file, append = TRUE)
    }
    
    # Cell type progression along pseudotime
    if ("slingshot" %in% names(trajectory_results)) {
      create_celltype_progression_plot(viz_data, dash_dir, cat)
      cat("  - Cell type progression plot completed\n", file = log_file, append = TRUE)
    }
    
  }, error = function(e) {
    message("Interactive plot generation failed: ", conditionMessage(e))
    cat("Interactive plot generation failed: ", conditionMessage(e), "\n", file = log_file, append = TRUE)
  })
  
  cat("Visualization generation completed\n")
  cat("Visualization generation completed\n", file = log_file, append = TRUE)
  }  # Close if (!is.null(viz_data))
}

# End performance tracking
end_performance_tracking(
  perf_tracker,
  success = TRUE,
  additional_metrics = list(
    n_cells = ncol(seu)
  )
)

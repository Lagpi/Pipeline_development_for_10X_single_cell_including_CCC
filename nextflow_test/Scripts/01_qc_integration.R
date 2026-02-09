#!/usr/bin/env Rscript
#===============================================================================
# 01_qc_integration.R
# 
# Single-cell RNA-seq QC and Integration Pipeline
# 
# Main Processing Steps:
#   1. Sample-level QC: filtering, SoupX, doublet detection, normalization
#   2. Comprehensive QC visualization (before/after)
#   3. Dataset integration using STACAS (with SingleR annotation)
#   4. Quality metrics calculation
#   5. Final output generation and export
#
# Functions defined:
#   - perform_qc_and_filtering(): Adaptive QC thresholds and filtering
#   - remove_ambient_rna(): SoupX contamination removal
#   - detect_doublets(): scDblFinder doublet detection
#   - normalize_and_select(): Normalization and HVG selection
#   - integrate_stacas(): STACAS integration with SingleR
#   - perform_integration(): Main integration wrapper

# Output files:
#   - seurat_qc.rds: Unnormalized merged object after QC filtering
#   - seurat_normalized.rds: Normalized merged object (all genes, no integration)
#   - seurat_normalized_hvg.rds: Normalized merged object (HVGs only, no integration)
#   - integrated_seurat.rds: Seurat CCA integrated object
#   - integrated_stacas.rds: STACAS integrated object
#   - integrated_primary.rds: Primary integrated object for downstream analysis
#===============================================================================


#===============================================================================
# 0.1. SETUP & INITIALIZATION
#===============================================================================
# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(plotly)
  library(patchwork)
})

# Load utility functions
script_dir <- {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    dirname(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    getwd()
  }
}
source(file.path(script_dir, "utils", "common.R"))
source(file.path(script_dir, "utils", "logging.R"))
source(file.path(script_dir, "utils", "plotting.R"))
source(file.path(script_dir, "utils", "colors.R"))
source(file.path(script_dir, "utils", "performance.R"))
source(file.path(script_dir, "utils", "01_qc_plots.R"))
safe_lib("ggrepel")
safe_lib("pheatmap")
safe_lib("cowplot")

# Parse command line arguments
option_list <- list(
  make_option(c("-s","--samplesheet"), type="character", help="Path to samplesheet TSV"),
  make_option(c("--starsolo_dir"), type="character", help="Root directory of STARsolo outputs"),
  make_option(c("-o","--output_dir"), type="character", default=".", help="Output base directory"),
  make_option(c("--config"), type="character", default=NULL, help="Path to nextflow.config file")
)
opt <- parse_args(OptionParser(option_list=option_list))
if(is.null(opt$samplesheet) || is.null(opt$starsolo_dir)) stop("Please provide --samplesheet and --starsolo_dir")

# Read parameters from nextflow.config
nf_params <- parse_nextflow_config(opt$config)

# Setup directories and logging (all-in-one from common.R)
setup_analysis_env(opt$output_dir, "01_qc_integration")

log("Starting 01_qc integration step")

#===============================================================================
# 0.2. CONFIGURATION PARAMETERS
#===============================================================================
# All parameters read from nextflow.config (params.qc block)
# Fallback to defaults if config not available

# Reproducibility
CFG_RANDOM_SEED <- as.integer(nf_params$qc.random_seed %||% 42)

# QC and Filtering Parameters
CFG_MIN_CELLS_PER_GENE <- as.integer(nf_params$qc.min_cells_per_gene %||% 3)
CFG_MIN_FEATURES_PER_CELL <- as.integer(nf_params$qc.min_features_per_cell %||% 200)
CFG_QC_QUANTILE_LOW <- as.numeric(nf_params$qc.quantile_low %||% 0.05)
CFG_QC_QUANTILE_HIGH <- as.numeric(nf_params$qc.quantile_high %||% 0.95)

# Feature Selection Parameters
CFG_N_VARIABLE_FEATURES <- as.integer(nf_params$qc.n_variable_features %||% 3000)
CFG_TOP_HVG_LABEL <- as.integer(nf_params$qc.top_hvg_label %||% 10)
CFG_TOP_HVG_DOTPLOT <- as.integer(nf_params$qc.top_hvg_dotplot %||% 10)
CFG_HVG_PER_CONDITION <- as.integer(nf_params$qc.hvg_per_condition %||% 3000)

# Integration Parameters
integration_dims <- as.integer(nf_params$qc.integration_dims %||% 30)
CFG_INTEGRATION_DIMS <- 1:integration_dims
CFG_CLUSTERING_RESOLUTION <- as.numeric(nf_params$qc.clustering_resolution %||% 0.5)

# Dimensionality Reduction Parameters
CFG_PCA_NPCS <- as.integer(nf_params$qc.pca_npcs %||% 30)
umap_dims <- as.integer(nf_params$qc.umap_dims %||% 30)
tsne_dims <- as.integer(nf_params$qc.tsne_dims %||% 30)
CFG_UMAP_DIMS <- 1:umap_dims
CFG_TSNE_DIMS <- 1:tsne_dims

# Plotting Parameters
CFG_MAX_CELLS_INTERACTIVE <- as.integer(nf_params$qc.max_cells_interactive %||% 30000)
CFG_PLOT_WIDTH <- as.numeric(nf_params$qc.plot_width %||% 12)
CFG_PLOT_HEIGHT <- as.numeric(nf_params$qc.plot_height %||% 6)
CFG_BASE_FONT_SIZE <- as.integer(nf_params$qc.base_font_size %||% 18)
CFG_TITLE_FONT_SIZE <- as.integer(nf_params$qc.title_font_size %||% 22)
CFG_AXIS_TITLE_SIZE <- as.integer(nf_params$qc.axis_title_size %||% 18)

# Silhouette Score Parameters
CFG_SILHOUETTE_SUBSAMPLE <- as.integer(nf_params$qc.silhouette_subsample %||% 5000)
CFG_SILHOUETTE_THRESHOLD <- as.numeric(nf_params$qc.silhouette_threshold %||% 0.3)

# Set global random seed for reproducibility
set.seed(CFG_RANDOM_SEED)

# Increase memory limit for large dataset integration 
options(future.globals.maxSize = 8 * 1024^3)  # 8 GB for large datasets

# Start performance tracking
perf_tracker <- start_performance_tracking("01_qc_integration", base_dir)

#===============================================================================
# 0.3. DATA LOADING FUNCTIONS
#===============================================================================

#' Load 10X Genomics count matrix
#' 
#' Reads gzipped 10X format files (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
#' and returns a sparse count matrix. Handles both single and multi-assay outputs.
#' 
#' @param matrix_path Path to directory containing 10X files
#' @param sample_name Sample identifier for logging
#' @return Sparse matrix (genes x cells) or NULL if loading fails
load_sample_matrix <- function(matrix_path, sample_name) {
  required_files_gz <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  
  if(!all(file.exists(file.path(matrix_path, required_files_gz)))) {
    log("Matrix not found for", sample_name, "at:", matrix_path)
    return(NULL)
  }
  
  counts <- Read10X(matrix_path)
  
  # Extract gene expression matrix if multiple assays are present
  if (is.list(counts)) {
    counts <- if("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
  }
  
  return(counts)
}

#===============================================================================
# 1. QC AND PREPROCESSING FUNCTIONS
#===============================================================================

#' Perform quality control and cell filtering
#' 
#' Creates Seurat object, calculates comprehensive QC metrics (mitochondrial %,
#' ribosomal %), determines adaptive filtering thresholds based on quantiles,
#' and filters low-quality cells. Uses data-driven adaptive thresholds rather
#' than fixed cutoffs.
#' 
#' @param counts Raw count matrix (genes x cells)
#' @param sample_name Sample identifier for metadata
#' @param condition Experimental condition for metadata
#' @param plots_dir Directory for plots (currently unused)
#' @param data_dir Directory for data outputs (currently unused)
#' @return List with:
#'   - obj: Filtered Seurat object
#'   - orig: Original unfiltered Seurat object
#'   - thresholds: Data frame with filtering thresholds
#'   - stats: List of filtering statistics
perform_qc_and_filtering <- function(counts, sample_name, condition, plots_dir, data_dir) {
  # 1. Create Seurat object with basic filtering
  seu <- CreateSeuratObject(counts=counts, project=sample_name, 
                           min.cells=CFG_MIN_CELLS_PER_GENE, 
                           min.features=CFG_MIN_FEATURES_PER_CELL)
  if(ncol(seu) == 0) return(NULL)
  
  # 2. Calculate QC metadata
  seu$sample <- sample_name
  seu$condition <- condition
  
  # Mitochondrial content 
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^mt-|^MT-")
  # Ribosomal content 
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern="^Rps|^Rpl|^RPS|^RPL")

  # Store pre-filtering metrics for comparison
  n_cells_before <- ncol(seu)
  total_umi_before <- sum(seu$nCount_RNA)
  median_genes_before <- median(seu$nFeature_RNA)
  
  # 3. Adaptive QC thresholds (based on data distribution)
  nf_low <- quantile(seu$nFeature_RNA, CFG_QC_QUANTILE_LOW)
  nf_high <- quantile(seu$nFeature_RNA, CFG_QC_QUANTILE_HIGH)
  mt_high <- quantile(seu$percent.mt, CFG_QC_QUANTILE_HIGH)
  
  # 4. Apply QC filters
  seu$qc_pass <- seu$nFeature_RNA > nf_low & 
                 seu$nFeature_RNA < nf_high & 
                 seu$percent.mt < mt_high
  
  seu_filtered <- subset(seu, subset = qc_pass == TRUE)
  
  if(ncol(seu_filtered) == 0) return(NULL)
  n_cells_after <- ncol(seu_filtered)
  
  # 5. Return results with comprehensive statistics
  return(list(
    obj = seu_filtered, 
    orig = seu,
    thresholds = data.frame(sample=sample_name, nf_low=nf_low, nf_high=nf_high, mt_high=mt_high),
    stats = list(
      sample = sample_name, 
      condition = condition,
      cells_before = n_cells_before, 
      cells_after = n_cells_after,
      umi_before = total_umi_before, 
      umi_after = sum(seu_filtered$nCount_RNA),
      median_genes_before = median_genes_before, 
      median_genes_after = median(seu_filtered$nFeature_RNA)
    )
  ))
}

#' Remove ambient RNA contamination using SoupX
#' 
#' Estimates and removes ambient RNA contamination by comparing raw (unfiltered)
#' and filtered count matrices. Uses quick clustering to estimate contamination
#' fraction and adjusts counts accordingly.
#' 
#' @param seu Filtered Seurat object
#' @param sample_name Sample identifier for logging
#' @param matrix_path Path to filtered matrix directory (parent contains raw/)
#' @param plots_dir Directory for plots (currently unused)
#' @return List with:
#'   - obj: Seurat object with corrected counts
#'   - stats: Contamination statistics (contamination %, UMI reduction %)
remove_ambient_rna <- function(seu, sample_name, matrix_path, plots_dir) {
  # Check if SoupX package is available
  if(!requireNamespace("SoupX", quietly=TRUE)) {
    log("SoupX not available, skipping ambient RNA removal")
    return(list(obj = seu, stats = NULL))
  }
  
  res <- list(obj = seu, stats = NULL)
  
  # SoupX expects parent directory containing raw/ and filtered/ subdirectories
  # Get parent directory from matrix_path (should be like .../filtered)
  parent_dir <- dirname(matrix_path)
  
  # Check if raw directory exists
  raw_dir <- file.path(parent_dir, "raw")
  if(!dir.exists(raw_dir)) {
    log("SoupX: No raw matrix directory found, skipping ambient RNA removal for", sample_name)
    return(list(obj = seu, stats = NULL))
  }
  
    # Load matrices manually
    log("  Loading filtered matrix...")
    toc <- Seurat::Read10X(matrix_path)
    log("  Loading raw matrix...")
    tod <- Seurat::Read10X(raw_dir)
    
    # Create SoupChannel with soup profile
    log("  Creating SoupChannel...")
    sc <- SoupX::SoupChannel(tod = tod, toc = toc, calcSoupProfile = TRUE)
    
    # Quick clustering directly on toc for SoupX
    log("  Performing quick clustering for contamination estimation...")
    seu_toc <- CreateSeuratObject(counts = toc)
    seu_toc <- NormalizeData(seu_toc, verbose = FALSE)
    seu_toc <- FindVariableFeatures(seu_toc, nfeatures = 3000, verbose = FALSE)
    seu_toc <- ScaleData(seu_toc, verbose = FALSE)
    seu_toc <- RunPCA(seu_toc, npcs = 50, verbose = FALSE)
    seu_toc <- FindNeighbors(seu_toc, dims = 1:50, verbose = FALSE)
    seu_toc <- FindClusters(seu_toc, resolution = 0.8, verbose = FALSE)
    
    # Extract clusters for toc cells
    clusters <- setNames(as.character(seu_toc$seurat_clusters), colnames(seu_toc))
    
    # Set clusters in SoupChannel
    sc <- SoupX::setClusters(sc, clusters)
    
    # Automatic contamination estimation with clusters
    log("  Estimating contamination with autoEstCont...")
    sc <- SoupX::autoEstCont(sc, doPlot = FALSE)
    
    # Extract contamination fraction
    contam <- sc$metaData$rho[1]
    log(sprintf("  Estimated contamination: %.1f%%", contam * 100))
    
    # Adjust counts
    log("  Adjusting counts...")
    corrected <- SoupX::adjustCounts(sc, roundToInt = TRUE)

    # Create NEW Seurat object from corrected counts (don't replace in filtered seu)
    # corrected has all toc cells, seu has filtered cells - need to filter corrected to match seu
    seu_cells <- colnames(seu)
    corrected_filtered <- corrected[, seu_cells]
    
    # Replace counts in Seurat object
    pre_umi <- sum(seu$nCount_RNA)
    seu[["RNA"]] <- CreateAssayObject(counts = corrected_filtered)
    post_umi <- sum(seu$nCount_RNA)

    # Calculate contamination metrics
    res$stats <- list(
      contamination = contam, 
      umi_reduction = 100 * (pre_umi - post_umi) / pre_umi
    )
    res$obj <- seu
    
    log(sprintf("SoupX: Removed %.2f%% ambient RNA from %s (UMI reduction: %.2f%%)", 
                contam*100, sample_name, res$stats$umi_reduction))

  return(res)
}



#' Detect and remove doublets using scDblFinder
#' 
#' Identifies potential doublets (two cells captured as one) using scDblFinder.
#' Simulates artificial doublets and trains a classifier to predict real doublets.
#' Returns singlet-only Seurat object.
#' 
#' @param seu Filtered Seurat object
#' @param sample_name Sample identifier for logging
#' @param plots_dir Directory for plots (currently unused)
#' @return List with:
#'   - obj: Filtered Seurat object (singlets only)
#'   - stats: Doublet detection statistics (n doublets, % removed)
detect_doublets <- function(seu, sample_name, plots_dir) {
  # Check if scDblFinder package is available
  if(!requireNamespace("scDblFinder", quietly=TRUE)) {
    log("scDblFinder not available, skipping doublet detection")
    return(list(obj = seu, stats = NULL))
  }
  
  res <- list(obj = seu, stats = NULL)
  library(SingleCellExperiment)
  
  # Convert to SingleCellExperiment format required by scDblFinder
  n_before <- ncol(seu)
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  # Run doublet detection (simulates doublets and trains classifier)
  sce <- scDblFinder::scDblFinder(sce)
  
  # Transfer doublet scores and classifications back to Seurat object
  seu$doublet_score <- colData(sce)$scDblFinder.score
  seu$scDblFinder.class <- colData(sce)$scDblFinder.class
  
  # Filter out predicted doublets
  n_doublets <- sum(seu$scDblFinder.class == "doublet")
  seu_filtered <- subset(seu, subset = scDblFinder.class == "singlet")
  
  # Compile statistics
  res$obj <- seu_filtered
  res$stats <- list(
    sample_name = sample_name,
    cells_before = n_before, 
    cells_after = ncol(seu_filtered),
    doublets_detected = n_doublets,
    percent_removed = 100 * n_doublets / n_before,
    raw_df = data.frame(
      doublet_score = seu$doublet_score,
      doublet_class = seu$scDblFinder.class,
      nFeature_RNA = seu$nFeature_RNA,
      nCount_RNA = seu$nCount_RNA,
      sample_name = sample_name
    )
  )
  
  log(sprintf("scDblFinder: Removed %d doublets (%.2f%%) from %s", 
              n_doublets, res$stats$percent_removed, sample_name))
  
  return(res)
}

#' Normalize data and select highly variable features
#' 
#' Performs log-normalization (scale factor 10,000) to correct for sequencing
#' depth differences, identifies highly variable genes (HVGs) that capture
#' biological variation, and applies z-score scaling for downstream PCA.
#' 
#' @param seu Filtered Seurat object
#' @param sample_name Sample identifier for logging
#' @param plots_dir Directory for plots (currently unused)
#' @return List with:
#'   - obj: Normalized and scaled Seurat object
#'   - stats: Normalization statistics (n HVGs)
normalize_and_select <- function(seu, sample_name, plots_dir) {
  # Log-normalization (default: scale factor 10,000)
  seu <- NormalizeData(seu, verbose = FALSE)
  
  # Identify highly variable genes (3000 HVGs)
  seu <- FindVariableFeatures(seu, nfeatures = 3000, verbose = FALSE)
  
  # Z-score scaling (mean=0, variance=1) for PCA
  seu <- ScaleData(seu, verbose = FALSE)
  
  return(list(
    obj = seu,
    stats = list(
      sample = sample_name,
      features_hvg = length(VariableFeatures(seu))
    )
  ))
}

#===============================================================================
# 2. INTEGRATION FUNCTIONS
#===============================================================================

#' Optimized Seurat CCA integration
#' 
#' Integrates multiple samples using Canonical Correlation Analysis (CCA).
#' Identifies integration anchors, performs integration, and generates
#' dimensionality reductions (PCA, UMAP, t-SNE) and clusters.
#' 
#' @param obj.list List of normalized Seurat objects to integrate
#' @param hvg_count Number of highly variable genes for integration (default: 3000)
#' @param n_dims Number of dimensions for CCA integration (default: 30)
#' @return List with:
#'   - obj: Integrated Seurat object with reductions and clusters
#'   - time: Integration time in minutes
#'   - success: Integration success status (TRUE/FALSE)
#'   - method: Method name ("Seurat_CCA")
run_optimized_seurat_integration <- function(obj.list, 
                                             hvg_count = 3000,
                                             n_dims = 30) {
  
  log(sprintf("Running optimized Seurat CCA integration with %d HVGs and %d dims", hvg_count, n_dims))
  start_time <- Sys.time()
  
  # 1. Select integration features
  log("Selecting integration features...")
  features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = hvg_count)
  log(sprintf("Selected %d integration features", length(features)))
  
  # 2. Scale data and run PCA on each object
  log("Scaling data and running PCA on each sample...")
  obj.list <- lapply(obj.list, function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, npcs = 50, verbose = FALSE)
    return(x)
  })
  
  # 3. Find integration anchors using CCA
  log("Finding integration anchors using CCA...")
  anchors <- FindIntegrationAnchors(
    object.list = obj.list,
    anchor.features = features,
    reduction = "cca",
    dims = 1:n_dims,
    verbose = FALSE
  )
  
  # 4. Integrate data
  log("Integrating data...")
  seurat_integrated <- IntegrateData(
    anchorset = anchors,
    dims = 1:n_dims,
    verbose = FALSE
  )
  
  # Set integrated assay as default
  DefaultAssay(seurat_integrated) <- "integrated"
  
  # Find variable features on integrated data
  log("Finding variable features on integrated data...")
  seurat_integrated <- FindVariableFeatures(seurat_integrated, 
                                           nfeatures = hvg_count, 
                                           verbose = FALSE)
  log(sprintf("Identified %d variable features", length(VariableFeatures(seurat_integrated))))
  
  # Scale for PCA
  log("Scaling for PCA computation...")
  seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)
  
  # Run PCA
  log("Running PCA on integrated data...")
  seurat_integrated <- RunPCA(seurat_integrated, npcs = 50, verbose = FALSE)
  
  # UMAP and clustering
  log("Running UMAP...")
  seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:n_dims,
                              seed.use = CFG_RANDOM_SEED, verbose = FALSE)
  
  log("Running t-SNE...")
  seurat_integrated <- RunTSNE(seurat_integrated, reduction = "pca", dims = 1:n_dims,
                              seed.use = CFG_RANDOM_SEED, verbose = FALSE, check_duplicates = FALSE)
  
  log("Finding neighbors and clustering...")
  seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", dims = 1:n_dims, verbose = FALSE)
  seurat_integrated <- FindClusters(seurat_integrated, resolution = CFG_CLUSTERING_RESOLUTION,
                                   random.seed = CFG_RANDOM_SEED, verbose = FALSE)
  log(sprintf("Identified %d clusters", length(unique(seurat_integrated$seurat_clusters))))
  
  integration_time <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
  log(sprintf("Seurat CCA integration completed in %.2f minutes", integration_time))
  
  return(list(
    obj = seurat_integrated,
    time = integration_time,
    success = TRUE,
    method = "Seurat_CCA"
  ))
}

#' Wrapper for Seurat CCA integration with visualization
#' 
#' Runs optimized Seurat integration and generates diagnostic plots comparing
#' dimensionality reductions.
#' 
#' @param sample_objects List of normalized Seurat objects
#' @param plots_dir Directory for static plots
#' @param dash_dir Directory for interactive plots
#' @return List with integrated object, timing, and success status
integrate_seurat <- function(sample_objects, plots_dir, dash_dir) {
  result <- run_optimized_seurat_integration(obj.list = sample_objects, n_dims = 30)
  
  # Create plots after integration
  if (result$success) {
    log("Creating Seurat integration plots...")
    plot_dimred_comparison(result$obj, "seurat_integrated", dash_dir, plots_dir)
  }
  
  return(result)
}

#' Optimized STACAS integration
#' 
#' STACAS (Semi-supervised Integration with Reference) uses semi-supervised
#' rPCA to preserve biological variation better than standard CCA. Can leverage
#' cell type annotations to guide integration.
#' 
#' @param obj.list List of normalized Seurat objects to integrate
#' @param hvg_count Number of highly variable genes (default: 3000)
#' @param annotation_col Metadata column with cell type labels (NULL = unsupervised)
#' @param n_dims Number of dimensions for integration (default: 30)
#' @return List with:
#'   - obj: Integrated Seurat object with reductions and clusters
#'   - time: Integration time in minutes
#'   - success: Integration success status (TRUE/FALSE)
#'   - method: Method name ("STACAS_optimized")
run_best_integration <- function(obj.list, 
                                 hvg_count = 3000,
                                 annotation_col = NULL, 
                                 n_dims = 30) {
  
  log(sprintf("Running optimized STACAS integration with %d HVGs and %d dims", hvg_count, n_dims))
  start_time <- Sys.time()
  
  # 1. Select integration features HVG
  log("Selecting integration features...")
  features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = hvg_count)
  log(sprintf("Selected %d integration features", length(features)))
  
  # 2. Run STACAS integration
  log("Running STACAS integration...")
  log(sprintf("  Parameters: alpha=1.0, k.anchor=5, dims=1:%d", n_dims))
  
  if (!is.null(annotation_col)) {
    log(sprintf("  Using semi-supervised mode with annotation column: %s", annotation_col))
  } else {
    log("  Using unsupervised mode (similar to Seurat rPCA)")
  }
  
  stacas_integrated <- STACAS::Run.STACAS(
    object.list = obj.list,
    anchor.features = features,     
    cell.labels = annotation_col,                         
    dims = 1:n_dims
  )
  log("STACAS integration completed - keeping unscaled data to preserve variance")
  
  # Set integrated assay as default
  DefaultAssay(stacas_integrated) <- "integrated"
  
  # Find variable features on integrated data
  log("Finding variable features on integrated data...")
  stacas_integrated <- FindVariableFeatures(stacas_integrated, 
                                           nfeatures = hvg_count, 
                                           verbose = FALSE)
  log(sprintf("Identified %d variable features", length(VariableFeatures(stacas_integrated))))
  
  # Scale only for PCA (but preserve original integrated slot unscaled)
  log("Scaling for PCA computation...")
  stacas_integrated <- ScaleData(stacas_integrated, verbose = FALSE)
  
  # Run PCA
  log("Running PCA on integrated data...")
  stacas_integrated <- RunPCA(stacas_integrated, npcs = 50, verbose = FALSE)
  
  # UMAP and clustering
  log("Running UMAP...")
  stacas_integrated <- RunUMAP(stacas_integrated, reduction = "pca", dims = 1:n_dims,
                              seed.use = CFG_RANDOM_SEED, verbose = FALSE)
  
  log("Running t-SNE...")
  stacas_integrated <- RunTSNE(stacas_integrated, reduction = "pca", dims = 1:n_dims,
                              seed.use = CFG_RANDOM_SEED, verbose = FALSE, check_duplicates = FALSE)
  
  log("Finding neighbors and clustering...")
  stacas_integrated <- FindNeighbors(stacas_integrated, reduction = "pca", dims = 1:n_dims, verbose = FALSE)
  stacas_integrated <- FindClusters(stacas_integrated, resolution = CFG_CLUSTERING_RESOLUTION,
                                   random.seed = CFG_RANDOM_SEED, verbose = FALSE)
  log(sprintf("Identified %d clusters", length(unique(stacas_integrated$seurat_clusters))))
  
  integration_time <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
  log(sprintf("STACAS integration completed in %.2f minutes", integration_time))
  
  return(list(
    obj = stacas_integrated,
    time = integration_time,
    success = TRUE,
    method = "STACAS_optimized"
  ))
}

#' Wrapper for STACAS integration with SingleR annotation
#' 
#' Annotates cell types using SingleR before integration, then runs STACAS
#' in semi-supervised mode using these annotations. Falls back to unsupervised
#' mode if SingleR/celldex are unavailable.
#' 
#' @param sample_objects List of normalized Seurat objects
#' @param plots_dir Directory for static plots
#' @param dash_dir Directory for interactive plots
#' @return List with integrated object, timing, and success status
integrate_stacas <- function(sample_objects, plots_dir, dash_dir) {
  if(!requireNamespace("STACAS", quietly=TRUE)) {
    log("STACAS package not available, skipping")
    return(list(obj = NULL, time = NA, success = FALSE, method = "STACAS_not_available"))
  }
  
    log("Preparing STACAS integration with SingleR cell type annotation...")
    library(STACAS)
    
    # Annotate cell types using SingleR before integration
    if(!requireNamespace("SingleR", quietly=TRUE) || !requireNamespace("celldex", quietly=TRUE)) {
      log("SingleR/celldex not available, running STACAS in unsupervised mode")
      annotation_col <- NULL
      sample_objects_annotated <- sample_objects
    } else {
      log("Annotating cell types with SingleR (per sample)...")
      library(SingleR)
      library(celldex)
      
      # Load appropriate reference (adjust based on species)
      log("Loading MouseRNAseqData reference...")
      ref <- celldex::MouseRNAseqData()
      log("Reference loaded successfully")
      
      sample_objects_annotated <- lapply(names(sample_objects), function(sname) {
        log(paste("  Annotating", sname, "..."))
        x <- sample_objects[[sname]]
        
        # Convert to SingleCellExperiment for SingleR
        sce <- as.SingleCellExperiment(x)
        
        # Run SingleR annotation
        singler_results <- SingleR(test = sce, ref = ref, labels = ref$label.main)
        
        # Add cell type labels to Seurat object
        x$cell.type <- singler_results$labels
        log(sprintf("    Annotated %d cells with %d unique cell types", 
                   ncol(x), length(unique(singler_results$labels))))
        
        return(x)
      })
      names(sample_objects_annotated) <- names(sample_objects)
      log("SingleR annotation completed for all samples")
      annotation_col <- "cell.type"
    }
    
    # Run optimized STACAS integration
    integration_result <- run_best_integration(
      obj.list = sample_objects_annotated,
      annotation_col = annotation_col, # NULL if no SingleR, else "cell.type"
      n_dims = 30                     # Extended dims for robust integration
    )
    
    # Create plots after integration
    if (integration_result$success) {
      log("Creating STACAS integration plots...")
      plot_dimred_comparison(integration_result$obj, "stacas_integrated", dash_dir, plots_dir)
    }
    
    return(integration_result)
    
    return(list(obj = NULL, time = NA, success = FALSE, method = "STACAS_error"))
}

# Finalize integration - keep both methods for comparison
# Both Seurat CCA and STACAS are run and saved for comparison
# STACAS is used as primary output (better bio-conservation)
# Returns: List with both integrated objects and metadata
finalize_integration <- function(seurat_result, stacas_result) {
  
  log("=== Integration Comparison Summary ===")
  
  # Report Seurat CCA results
  if(seurat_result$success && !is.null(seurat_result$obj)) {
    log(sprintf("Seurat CCA Integration:"))
    log(sprintf("  Method: %s", seurat_result$method))
    log(sprintf("  Time: %.2f minutes", seurat_result$time))
    log(sprintf("  Cells: %d", ncol(seurat_result$obj)))
    log(sprintf("  Clusters: %d", length(unique(seurat_result$obj$seurat_clusters))))
  } else {
    log("Seurat CCA Integration: FAILED")
  }
  
  # Report STACAS results
  if(stacas_result$success && !is.null(stacas_result$obj)) {
    log(sprintf("STACAS Integration:"))
    log(sprintf("  Method: %s", stacas_result$method))
    log(sprintf("  Time: %.2f minutes", stacas_result$time))
    log(sprintf("  Cells: %d", ncol(stacas_result$obj)))
    log(sprintf("  Clusters: %d", length(unique(stacas_result$obj$seurat_clusters))))
    
    # Check if cell type annotations exist (from STACAS+SingleR)
    if("cell.type" %in% colnames(stacas_result$obj@meta.data)) {
      n_celltypes <- length(unique(stacas_result$obj$cell.type))
      log(sprintf("  Cell types (SingleR): %d", n_celltypes))
    }
  } else {
    log(sprintf("STACAS Integration: %s", 
                ifelse(!is.null(stacas_result$method), stacas_result$method, "FAILED")))
  }
  
  # Choose primary method for downstream analysis --> STACAS
  if(stacas_result$success && !is.null(stacas_result$obj)) {
    merged_seu <- stacas_result$obj
    primary_method <- stacas_result$method
    log(sprintf("\nUsing %s as primary method for downstream analysis", primary_method))
  } else {
    merged_seu <- seurat_result$obj
    primary_method <- seurat_result$method
    log(sprintf("\nUsing %s as primary method (STACAS unavailable)", primary_method))
  }
  
  log("=== End Integration Comparison ===\n")
  
  return(list(
    # Primary output (for downstream analysis)
    merged_seu = merged_seu,
    primary_method = primary_method,
    
    # Both integration results (for comparison)
    seurat_obj = seurat_result$obj,
    stacas_obj = stacas_result$obj,
    
    # Timing and success stats
    seurat_time = seurat_result$time,
    stacas_time = stacas_result$time,
    seurat_success = seurat_result$success,
    stacas_success = stacas_result$success
  ))
}

#' Main integration wrapper
#' 
#' Runs both Seurat CCA and STACAS integration methods for comparison.
#' Prefers STACAS for downstream analysis (better biological preservation)
#' but keeps both results for benchmarking.
#' 
#' @param sample_objects List of normalized Seurat objects
#' @param qc_data_after QC metadata after filtering (currently unused)
#' @param plots_dir Directory for static plots
#' @param dash_dir Directory for interactive plots
#' @return List with:
#'   - merged_seu: Primary integrated object (STACAS preferred)
#'   - primary_method: Method name for primary object
#'   - seurat_obj: Seurat CCA integrated object
#'   - stacas_obj: STACAS integrated object
#'   - seurat_time, stacas_time: Integration times
#'   - seurat_success, stacas_success: Success flags
perform_integration <- function(sample_objects, qc_data_after, plots_dir, dash_dir) {
  log("Running Seurat CCA integration method...")
  seurat_result <- integrate_seurat(sample_objects, plots_dir, dash_dir)

  log("Running STACAS integration method...")
  stacas_result <- integrate_stacas(sample_objects, plots_dir, dash_dir)

  # Prefer STACAS for downstream analysis if available
  if (!is.null(stacas_result$obj)) {
    merged_seu <- stacas_result$obj
    primary_method <- stacas_result$method
  } else {
    merged_seu <- seurat_result$obj
    primary_method <- seurat_result$method
  }

  # Return results in expected format (both objects for comparison)
  return(list(
    merged_seu = merged_seu,
    primary_method = primary_method,
    seurat_time = seurat_result$time,
    stacas_time = stacas_result$time,
    seurat_obj = seurat_result$obj,
    stacas_obj = stacas_result$obj,
    seurat_success = seurat_result$success,
    stacas_success = stacas_result$success
  ))
}

#===============================================================================
# 3. SAMPLE PROCESSING LOOP
#===============================================================================
# Per-sample workflow:
#   1. Load count matrix from STARsolo output
#   2. QC and filtering (adaptive thresholds)
#   3. Ambient RNA removal (SoupX)
#   4. Doublet detection (scDblFinder)
#   5. Normalization and HVG selection
#===============================================================================

log("=== READING SAMPLESHEET ===")

# Read samplesheet
ss <- read.table(opt$samplesheet, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
log(sprintf("Loaded samplesheet with %d samples", nrow(ss)))

# Initialize storage lists
sample_objects <- list()                # Normalized Seurat objects (all genes)
normalized_sample_objects <- list()     # Normalized Seurat objects (for integration)
qc_data_before <- list()                # QC metadata before filtering
qc_data_after <- list()                 # QC metadata after filtering
filtering_stats <- list()               # Cell filtering statistics
normalization_stats <- list()           # Normalization statistics
soupx_stats <- list()                   # SoupX contamination statistics
doublet_stats <- list()                 # Doublet detection statistics
successful_samples <- c()               # Successfully processed samples
skipped <- c()                          # Skipped samples

log("=== SAMPLE PROCESSING LOOP ===")

# Process each sample
for (i in 1:nrow(ss)) {
  sample_name <- ss$sample_name[i]
  condition <- ss$condition[i]
  
  log(sprintf("Processing sample %d/%d: %s (condition: %s)", i, nrow(ss), sample_name, condition))
  
  # 1. Load matrix
  # STARsolo directory structure: starsolo_dir/condition/original_sample/sample_name_Solo.out/Gene/filtered
  # Extract original sample name from sample_name (last part after underscore)
  sample_parts <- strsplit(sample_name, "_", fixed = TRUE)[[1]]
  original_sample <- sample_parts[length(sample_parts)]
  
  matrix_path <- file.path(opt$starsolo_dir, condition, original_sample, paste0(sample_name, "_Solo.out"), "Gene", "filtered")
  
  if (!dir.exists(matrix_path)) {
    log(sprintf("  WARNING: Matrix directory not found at: %s", matrix_path))
    skipped <- c(skipped, sample_name)
    next
  }
  
  counts <- load_sample_matrix(matrix_path, sample_name)
  
  if (is.null(counts)) {
    log(sprintf("  WARNING: Skipping %s - matrix could not be loaded", sample_name))
    skipped <- c(skipped, sample_name)
    next
  }
  
  log(sprintf("  Loaded %d genes x %d cells", nrow(counts), ncol(counts)))
  
  # 2. QC and Filtering
  log("  Performing QC and filtering...")
  qc_result <- perform_qc_and_filtering(counts, sample_name, condition, plots_dir, data_dir)
  
  if (is.null(qc_result) || is.null(qc_result$obj)) {
    log(sprintf("  WARNING: Skipping %s - no cells passed QC", sample_name))
    skipped <- c(skipped, sample_name)
    next
  }
  
  seu <- qc_result$obj
  filtering_stats[[sample_name]] <- qc_result$stats
  qc_data_before[[sample_name]] <- as.data.frame(qc_result$orig@meta.data)
  
  log(sprintf("  Retained %d/%d cells after QC (%.1f%%)", 
              qc_result$stats$cells_after, 
              qc_result$stats$cells_before,
              100 * qc_result$stats$cells_after / qc_result$stats$cells_before))
  
  # 3. Ambient RNA removal (SoupX)
  log("  Removing ambient RNA with SoupX...")
  soupx_result <- remove_ambient_rna(seu, sample_name, matrix_path, plots_dir)
  seu <- soupx_result$obj
  soupx_stats[[sample_name]] <- soupx_result$stats
  
  # 4. Doublet detection (scDblFinder)
  log("  Detecting doublets with scDblFinder...")
  doublet_result <- detect_doublets(seu, sample_name, plots_dir)
  seu <- doublet_result$obj
  doublet_stats[[sample_name]] <- doublet_result$stats
  
  if (ncol(seu) == 0) {
    log(sprintf("  WARNING: Skipping %s - no cells after doublet removal", sample_name))
    skipped <- c(skipped, sample_name)
    next
  }
  
  # 5. Normalization and feature selection
  log("  Normalizing and selecting features...")
  norm_result <- normalize_and_select(seu, sample_name, plots_dir)
  seu <- norm_result$obj
  normalization_stats[[sample_name]] <- norm_result$stats
  
  # Save normalized object to data directory
  normalized_rds_path <- file.path(data_dir, paste0(sample_name, "_normalized.rds"))
  saveRDS(seu, normalized_rds_path)
  log(sprintf("  Saved normalized object: %s", normalized_rds_path))
  
  # Save normalized expression matrix as data frame (HVG only)
  normalized_matrix <- GetAssayData(seu, assay = "RNA", slot = "data")
  hvg_genes <- VariableFeatures(seu)
  normalized_matrix_hvg <- as.data.frame(as.matrix(normalized_matrix[hvg_genes, ]))
  normalized_matrix_hvg$gene <- rownames(normalized_matrix_hvg)
  
  normalized_csv_path <- file.path(data_dir, paste0(sample_name, "_normalized_hvg.csv"))
  write.csv(normalized_matrix_hvg, normalized_csv_path, row.names = FALSE)
  log(sprintf("  Saved normalized HVG matrix: %s (%d genes x %d cells)", 
              normalized_csv_path, nrow(normalized_matrix_hvg), ncol(normalized_matrix_hvg)-1))
  
  # Store processed object
  sample_objects[[sample_name]] <- seu
  normalized_sample_objects[[sample_name]] <- seu  # Also store for integration
  qc_data_after[[sample_name]] <- as.data.frame(seu@meta.data)
  successful_samples <- c(successful_samples, sample_name)
  
  log(sprintf("  ✓ Sample %s processed successfully: %d cells, %d HVGs", 
              sample_name, ncol(seu), length(VariableFeatures(seu))))
}

log(sprintf("=== SAMPLE PROCESSING COMPLETE: %d/%d samples successful ===", 
            length(successful_samples), nrow(ss)))

if (length(skipped) > 0) {
  log(sprintf("Skipped samples: %s", paste(skipped, collapse = ", ")))
}

if (length(sample_objects) == 0) {
  stop("ERROR: No samples were successfully processed. Cannot proceed with integration.")
}


# Create combined normalized Seurat object (with only 3000 HVGs, no annotation yet)
# This is for cellcomm and geneprogramm discovery - normalized data only, no integration
log("Creating combined normalized Seurat object (3000 HVGs only, no integration)...")
normalized_combined <- NULL
if(length(normalized_sample_objects) > 0) {
  # Merge all normalized objects
  normalized_combined <- merge(
    x = normalized_sample_objects[[1]],
    y = if(length(normalized_sample_objects) > 1) normalized_sample_objects[2:length(normalized_sample_objects)] else NULL,
    add.cell.ids = names(normalized_sample_objects),
    project = "normalized_combined"
  )
  
  # Subset to only the HVG genes (should be 3000 or fewer across all samples)
  hvg_union <- unique(unlist(lapply(normalized_sample_objects, function(obj) VariableFeatures(obj))))
  normalized_combined <- normalized_combined[hvg_union, ]
  
  log(sprintf("  Combined normalized object: %d cells, %d HVG genes", 
              ncol(normalized_combined), nrow(normalized_combined)))
  log(sprintf("  Samples: %s", paste(unique(normalized_combined$sample), collapse = ", ")))
  
  seurat_normalized_hvg_file <- file.path(data_dir, "seurat_normalized_hvg.rds")
  saveRDS(normalized_combined, seurat_normalized_hvg_file)
  log(sprintf("Saved seurat_normalized_hvg (HVGs only): %s", seurat_normalized_hvg_file))
}

#===============================================================================
# 4. COMPREHENSIVE QC VISUALIZATION
#===============================================================================
# Generate comprehensive visualization of QC results:
#   - Before/after QC metrics (boxplots, scatter, violin, heatmap)
#   - Filtering summary statistics per sample
#   - Doublet detection results and distributions
#   - Ambient RNA contamination levels (SoupX)
#   - Highly variable gene (HVG) selection per sample and condition
#   - Pre-integration dimensionality reduction (PCA, UMAP, t-SNE)
#===============================================================================

log("=== GENERATING COMPREHENSIVE QC PLOTS ===")

# Generate QC plots (boxplots and scatter - existing function)
if (length(qc_data_before) > 0) {
  log("Generating before-QC plots...")
  generate_integrated_qc_plots(qc_data_before, "before_QC", plots_dir, dash_dir)
  
  # Additional violin plots
  generate_violin_plots(qc_data_before, "before", plots_dir, dash_dir)
  
  # QC heatmap
  if(requireNamespace("pheatmap", quietly = TRUE)) {
    generate_qc_heatmap(qc_data_before, "before", plots_dir)
  }
}

if (length(qc_data_after) > 0) {
  log("Generating after-QC plots...")
  generate_integrated_qc_plots(qc_data_after, "after_QC", plots_dir, dash_dir)
  
  # Additional violin plots
  generate_violin_plots(qc_data_after, "after", plots_dir, dash_dir)
  
  # QC heatmap
  if(requireNamespace("pheatmap", quietly = TRUE)) {
    generate_qc_heatmap(qc_data_after, "after", plots_dir)
  }
}

# Filtering summary
if (length(filtering_stats) > 0) {
  log("Generating filtering summary plot...")
  generate_filtering_summary(filtering_stats, plots_dir)
}

# Doublet removal summary
if (length(doublet_stats) > 0) {
  log("Generating doublet detection summary...")
  generate_doublet_summary_from_stats(doublet_stats, plots_dir)
}

# Ambient RNA contamination summary
if (length(soupx_stats) > 0) {
  log("Generating ambient RNA contamination plot...")
  generate_ambient_rna_plot(soupx_stats, plots_dir)
}

# HVG plots
log("Generating HVG selection plots...")
generate_hvg_plots(normalized_sample_objects, plots_dir)
generate_hvg_per_condition_plot(normalized_sample_objects, plots_dir)
generate_hvg_dotplot(normalized_sample_objects, plots_dir, top_n = CFG_TOP_HVG_DOTPLOT)

# Pre-integration dimensionality reduction comparison
log("Generating pre-integration comparison plot...")
pre_integrated_obj <- generate_dimred_comparison_pre(normalized_sample_objects, plots_dir)

#===============================================================================
# 5. INTEGRATION ANALYSIS
#===============================================================================
# Integrate normalized samples using Seurat CCA and STACAS methods
# Primary method (STACAS) is used for downstream analysis
#===============================================================================

log("=== INTEGRATION PHASE (using normalized, 3000 HVG data) ===")

# Perform integration on NORMALIZED data
integration_results <- perform_integration(normalized_sample_objects, qc_data_after, plots_dir, dash_dir)
merged_seu <- integration_results$merged_seu
primary_method <- integration_results$primary_method
seurat_time <- integration_results$seurat_time
stacas_time <- integration_results$stacas_time
seurat_obj <- integration_results$seurat_obj
stacas_obj <- integration_results$stacas_obj

log(sprintf("Integration completed using method: %s (on normalized 3000 HVG data)", primary_method))

# Post-integration PCA comparison: Seurat vs STACAS side by side
log("Generating post-integration PCA comparison plot (Seurat vs STACAS)...")
if(!is.null(seurat_obj) && !is.null(stacas_obj)) {
  # Seurat PCA
  p_seurat_pca <- DimPlot(seurat_obj, reduction = "pca", group.by = "condition") +
    labs(title = "Seurat CCA - PCA") +
    theme_minimal(base_size = 16)
  
  # STACAS PCA
  p_stacas_pca <- DimPlot(stacas_obj, reduction = "pca", group.by = "condition") +
    labs(title = "STACAS - PCA") +
    theme_minimal(base_size = 16)
  
  # Save side by side comparison
  p_combined <- p_seurat_pca + p_stacas_pca
  ggsave(file.path(plots_dir, "pca_comparison_seurat_vs_stacas.png"), 
         p_combined, width = 16, height = 6, dpi = 150)
  log("Saved PCA comparison plot: Seurat vs STACAS")
} else if(!is.null(stacas_obj)) {
  # Only STACAS available
  log("Only STACAS integration available - STACAS plots already generated")
} else if(!is.null(seurat_obj)) {
  # Only Seurat available
  generate_dimred_comparison_seurat(seurat_obj, plots_dir)
}

# Generate comprehensive 3x3 integration comparison grid
log("Generating 3x3 integration comparison grid (Pre/Seurat/StacAS x PCA/UMAP/tSNE)...")
generate_integration_3x3_grid(
  pre_obj = pre_integrated_obj,
  seurat_obj = seurat_obj,
  stacas_obj = stacas_obj,
  plots_dir = plots_dir
)

#===============================================================================
# 6. QUALITY METRICS CALCULATION
#===============================================================================
# Calculate comprehensive quality metrics:
#   1. Processing stats (samples, cells, genes)
#   2. Sample-level QC metrics
#   3. Integration quality (silhouette scores)
#   4. Final dataset quality assessment
#===============================================================================

log("=== QUALITY METRICS CALCULATION ===")

# Calculate comprehensive quality metrics
quality_metrics <- list()

# 1. Overall processing metrics
quality_metrics$processing <- list(
  total_samples_attempted = nrow(ss),
  successful_samples = length(successful_samples), 
  success_rate = round(100 * length(successful_samples) / nrow(ss), 2),
  skipped_samples = length(skipped),
  total_cells_final = ncol(merged_seu),
  total_genes_detected = nrow(merged_seu)
)

# 2. Sample-level metrics (from stored stats)
if(exists("filtering_stats") && length(filtering_stats) > 0) {
  quality_metrics$filtering <- filtering_stats
}
if(exists("normalization_stats") && length(normalization_stats) > 0) {
  quality_metrics$normalization <- normalization_stats  
}

# 3. Integration quality metrics (if applicable)
if(length(successful_samples) > 1) {
  # Calculate condition mixing in UMAP space (subsample for performance)
  condition_silhouette <- try({
    n_cells <- ncol(merged_seu)
    if (n_cells > CFG_SILHOUETTE_SUBSAMPLE) {
      log(paste("Subsampling", CFG_SILHOUETTE_SUBSAMPLE, "cells for silhouette calculation (performance)"))
      cell_idx <- sample(n_cells, CFG_SILHOUETTE_SUBSAMPLE)
      umap_coords <- Embeddings(merged_seu, "umap")[cell_idx, ]
      conditions <- merged_seu$condition[cell_idx]
    } else {
      umap_coords <- Embeddings(merged_seu, "umap")
      conditions <- merged_seu$condition
    }
    
    # Calculate silhouette scores
    sil <- cluster::silhouette(as.numeric(as.factor(conditions)), dist(umap_coords))
    sil
  })
  
  if(!inherits(condition_silhouette, "try-error")) {
    sil_mean <- mean(condition_silhouette[,3])
    sil_median <- median(condition_silhouette[,3])
    
    quality_metrics$integration <- list(
      method = primary_method,  # Correct method name (Seurat RPCA or STACAS)
      conditions = length(unique(merged_seu$condition)),
      condition_silhouette_mean = sil_mean,
      condition_silhouette_median = sil_median,
      integration_quality = ifelse(sil_mean < CFG_SILHOUETTE_THRESHOLD, "Good", "Poor"),
      integration_success = sil_mean < CFG_SILHOUETTE_THRESHOLD  # Good mixing if < threshold
    )
    
    log(paste("Integration quality (Silhouette):", quality_metrics$integration$integration_quality))
    log(paste("  Mean silhouette score:", round(sil_mean, 3)))
    log(paste("  Median silhouette score:", round(sil_median, 3)))
  } else {
    log("Warning: Silhouette score calculation failed")
  }
}

# 4. Final data quality assessment
quality_metrics$final_quality <- list(
  median_genes_per_cell = median(merged_seu$nFeature_RNA),
  median_umi_per_cell = median(merged_seu$nCount_RNA),
  median_mt_percent = median(merged_seu$percent.mt),
  n_clusters_detected = length(unique(merged_seu$seurat_clusters)),
  highly_variable_genes = length(VariableFeatures(merged_seu))
)

log("Quality metrics calculated successfully")

#-------------------------------------------------------------------------------
# 6.1. Quality Metrics Visualization
#-------------------------------------------------------------------------------
# Generate final summary plots:
#   - Overall quality metrics (cells, genes, conditions, samples, HVGs)
#   - Per-condition QC statistics (mean features, counts, MT%)
#-------------------------------------------------------------------------------
generate_final_quality_plots(merged_seu, primary_method, plots_dir, dash_dir)

#===============================================================================
# 7. FINAL OUTPUT & EXPORT
#===============================================================================
# Save outputs:
#   - seurat_qc.rds: Unnormalized merged object after QC filtering
#   - seurat_normalized.rds: Normalized merged object (all genes, no integration)
#   - seurat_normalized_hvg.rds: Normalized merged object (HVGs only, no integration)
#   - integrated_seurat.rds: Seurat CCA integrated object
#   - integrated_stacas.rds: STACAS integrated object
#   - integrated_primary.rds: Primary integrated object for downstream analysis
#   - integration_comparison.rds: Benchmarking data 
#   - integration_comparison_table.csv: Summary comparison table
#   - quality_metrics.rds: Comprehensive QC and integration metrics
#   - integrated_metadata.csv: Cell-level metadata with all QC metrics
#===============================================================================

log("=== FINAL OUTPUT GENERATION ===")

# 1. Save merged QC object (unnormalized raw counts after QC)
log("Creating merged QC object (unnormalized raw counts)...")
# Extract raw counts from normalized objects and create unnormalized merge
unnormalized_objects <- lapply(names(sample_objects), function(sn) {
  obj <- sample_objects[[sn]]
  # Create new object with raw counts only
  raw_obj <- CreateSeuratObject(
    counts = GetAssayData(obj, assay = "RNA", slot = "counts"),
    meta.data = obj@meta.data,
    project = sn
  )
  return(raw_obj)
})
names(unnormalized_objects) <- names(sample_objects)

seurat_qc <- merge(x = unnormalized_objects[[1]], 
                   y = unnormalized_objects[-1], 
                   merge.data = TRUE,
                   project = "seurat_qc")
seurat_qc_file <- file.path(data_dir, "seurat_qc.rds")
saveRDS(seurat_qc, seurat_qc_file)
log(sprintf("Saved seurat_qc (unnormalized): %d cells x %d genes", 
            ncol(seurat_qc), nrow(seurat_qc)))

# 2. Save merged normalized object (all genes, no integration)
log("Creating merged normalized object (all genes)...")
seurat_normalized <- merge(x = sample_objects[[1]], 
                           y = sample_objects[-1], 
                           merge.data = TRUE,
                           project = "seurat_normalized")
DefaultAssay(seurat_normalized) <- "RNA"
seurat_normalized_file <- file.path(data_dir, "seurat_normalized.rds")
saveRDS(seurat_normalized, seurat_normalized_file)
log(sprintf("Saved seurat_normalized (all genes): %d cells x %d genes", 
            ncol(seurat_normalized), nrow(seurat_normalized)))

# 3. Save Seurat CCA integration
if(exists("seurat_obj") && !is.null(seurat_obj)) {
  integrated_seurat_file <- file.path(data_dir, "integrated_seurat.rds")
  saveRDS(seurat_obj, integrated_seurat_file)
  log("Saved integrated_seurat (Seurat CCA):", integrated_seurat_file)
}

# 4. Save STACAS integration
if(exists("stacas_obj") && !is.null(stacas_obj)) {
  integrated_stacas_file <- file.path(data_dir, "integrated_stacas.rds")
  saveRDS(stacas_obj, integrated_stacas_file)
  log("Saved integrated_stacas:", integrated_stacas_file)
}

# 5. Save PRIMARY integrated object (used for downstream - whichever was successful)
integrated_primary_file <- file.path(data_dir, "integrated_primary.rds") 
saveRDS(merged_seu, integrated_primary_file)
log(sprintf("Saved PRIMARY integrated object (method: %s):", primary_method), integrated_primary_file)

# Create comprehensive integration comparison for benchmarking
if(exists("seurat_obj") && exists("stacas_obj") && !is.null(seurat_obj) && !is.null(stacas_obj)) {
  integration_comparison <- list(
    # Timing data
    seurat_time = seurat_time,
    stacas_time = stacas_time,
    seurat_success = !is.null(seurat_obj),
    stacas_success = !is.null(stacas_obj),
    primary_method = primary_method,
    
    # Detailed metrics
    seurat_cca = list(
      n_cells = ncol(seurat_obj),
      n_genes = nrow(seurat_obj),
      n_clusters = length(unique(seurat_obj$seurat_clusters)),
      time_minutes = seurat_time,
      method = "Seurat_CCA_optimized",
      object_size_mb = as.numeric(object.size(seurat_obj)) / 1024^2
    ),
    stacas = list(
      n_cells = ncol(stacas_obj),
      n_genes = nrow(stacas_obj),
      n_clusters = length(unique(stacas_obj$seurat_clusters)),
      n_celltypes = if("cell.type" %in% colnames(stacas_obj@meta.data)) length(unique(stacas_obj$cell.type)) else NA,
      time_minutes = stacas_time,
      method = "STACAS_optimized",
      object_size_mb = as.numeric(object.size(stacas_obj)) / 1024^2
    )
  )
  
  # Add quality metrics if available
  if(exists("quality_metrics") && !is.null(quality_metrics$integration)) {
    integration_comparison$silhouette_score <- quality_metrics$integration$condition_silhouette_mean
    integration_comparison$integration_quality <- quality_metrics$integration$integration_quality
  }
  
  # Save as RDS for benchmarking
  comparison_file <- file.path(data_dir, "integration_comparison.rds")
  saveRDS(integration_comparison, comparison_file)
  log("Saved integration comparison report to:", comparison_file)
  
  # Also create CSV summary table
  comparison_table <- data.frame(
    Method = c("Seurat_CCA", "STACAS"),
    Time_Minutes = c(seurat_time, stacas_time),
    Success = c(TRUE, TRUE),
    N_Cells = c(ncol(seurat_obj), ncol(stacas_obj)),
    N_Clusters = c(
      length(unique(seurat_obj$seurat_clusters)),
      length(unique(stacas_obj$seurat_clusters))
    ),
    Object_Size_MB = c(
      integration_comparison$seurat_cca$object_size_mb,
      integration_comparison$stacas$object_size_mb
    ),
    stringsAsFactors = FALSE
  )
  
  write.csv(comparison_table, file.path(data_dir, "integration_comparison_table.csv"), row.names=FALSE)
  log("Saved integration comparison table to CSV")
}

# Save quality metrics
quality_file <- file.path(data_dir, "quality_metrics.rds")
saveRDS(quality_metrics, quality_file) 
log("Saved quality metrics to:", quality_file)

# Export metadata with all QC metrics
metadata_export <- merged_seu@meta.data
metadata_file <- file.path(data_dir, "integrated_metadata.csv")
write.csv(metadata_export, metadata_file, row.names = TRUE)
log("Exported metadata to:", metadata_file)

log("=== QC & INTEGRATION PIPELINE COMPLETED SUCCESSFULLY ===")
log(sprintf("Final dataset: %d cells x %d genes", ncol(merged_seu), nrow(merged_seu)))
log(sprintf("Conditions: %d (%s)", 
            length(unique(merged_seu$condition)), 
            paste(unique(merged_seu$condition), collapse=", ")))
log(sprintf("Integration method: %s", primary_method))
log(sprintf("Output directory: %s", data_dir))
log("Check plots_dir and dash_dir for visualizations")

# End performance tracking
end_performance_tracking(
  perf_tracker,
  success = TRUE,
  additional_metrics = list(
    n_cells = ncol(merged_seu),
    n_genes = nrow(merged_seu),
    n_conditions = length(unique(merged_seu$condition)),
    integration_method = primary_method
  )
)

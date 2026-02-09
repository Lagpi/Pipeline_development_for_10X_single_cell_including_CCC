#!/usr/bin/env Rscript
################################################################################
# 04b_nichenet.R
#
# NicheNet Analysis - Ligand-receptor prioritization
# 
# This script performs NicheNet analysis following the official vignette:
# - Predicts active ligands based on DE genes
# - Identifies ligand-target gene regulatory networks
# - Uses sender-agnostic approach (all potential ligands)
# - Analyzes each receiver cell type across all conditions
#
# Required inputs:
# - Clustered Seurat object with cell type annotations
# - DE results from step 03 (comparing conditions)
# - NicheNet reference files (lr_network, ligand_target_matrix)
#
# Outputs:
# - Ligand activity predictions per receiver
# - Target gene networks
# - Comprehensive visualizations
################################################################################

################################################################################
# SETUP
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(nichenetr)
  library(patchwork)
  library(plotly)
  library(htmlwidgets)
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
source(file.path(script_dir, "utils", "04_intercellular_plots.R"))

# Font configuration
CFG_BASE_SIZE <- 18
CFG_TITLE_SIZE <- 24
CFG_AXIS_TITLE_SIZE <- 20
CFG_LEGEND_SIZE <- 14

set.seed(CFG_RANDOM_SEED)

################################################################################
# ARGUMENT PARSING
################################################################################
option_list <- list(
  make_option("--seurat", type="character", default=NULL, help="Seurat object with annotations (.rds)"),
  make_option("--normalized", type="character", default=NULL, help="Normalized Seurat with all genes (.rds)"),
  make_option("--out", type="character", default="04_cellcomm", help="Output directory"),
  make_option("--species", type="character", default="mouse", help="Species (mouse/human)"),
  make_option("--nichenet_dir", type="character", default=NULL, help="Directory with NicheNet reference files")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$seurat)) {
  print_help(opt_parser)
  stop("Input Seurat object is required (--seurat).", call.=FALSE)
}

if (is.null(opt$nichenet_dir)) {
  stop("NicheNet reference directory is required (--nichenet_dir).", call.=FALSE)
}

# Setup directory structure
base_dir <- opt$out
dirs <- setup_analysis_directories(base_dir, "04_cellcomm")
data_dir <- dirs$data_dir
plots_dir <- dirs$plots_dir
dash_dir <- dirs$dash_dir
rds_dir <- dirs$rds_dir

# Parse nextflow config
config_path <- file.path(dirname(script_dir), "nextflow.config")
nf_params <- parse_nextflow_config(config_path)
pipeline_output_dir <- nf_params$output_dir %||% dirname(dirname(base_dir))

# Setup logging
log_file <- file.path(base_dir, "04b_nichenet.log")
setup_logging(log_file)

# Start performance tracking
perf_tracker <- start_performance_tracking("04b_nichenet", base_dir)

log_section("NicheNet Analysis (Official Vignette Approach)")
log("Input Seurat (annotations):", opt$seurat)
log("Input Seurat (normalized, all genes):", opt$normalized)
log("Output Dir:", base_dir)
log("Species:", opt$species)
log("NicheNet Dir:", opt$nichenet_dir)

################################################################################
# LOAD DATA
################################################################################

# Load annotated Seurat for cell type annotations
seu <- readRDS(opt$seurat)
log(sprintf("Loaded annotated Seurat object: %d cells, %d genes", ncol(seu), nrow(seu)))

# Ensure cell_type annotation exists
if (!"cell_type" %in% colnames(seu@meta.data)) {
  if ("seurat_clusters" %in% colnames(seu@meta.data)) {
    seu$cell_type <- paste0("Cluster_", seu$seurat_clusters)
  } else {
    stop("No cell_type or seurat_clusters found in Seurat object.")
  }
}

# Detect condition field
cond_field <- NULL
for (f in c("condition", "treatment", "group", "orig.ident")) {
  if (f %in% colnames(seu@meta.data)) {
    if (length(unique(seu@meta.data[[f]])) > 1) {
      cond_field <- f
      break
    }
  }
}
if (is.null(cond_field)) cond_field <- "orig.ident"

# Load normalized Seurat with ALL genes for expression analysis
if (!is.null(opt$normalized) && file.exists(opt$normalized)) {
  seu_full <- readRDS(opt$normalized)
  log(sprintf("Loaded normalized Seurat with all genes: %d cells, %d genes", ncol(seu_full), nrow(seu_full)))
  
  # Transfer annotations and reductions from clustered to normalized
  if (all(colnames(seu) %in% colnames(seu_full))) {
    seu_full$cell_type <- seu$cell_type[match(colnames(seu_full), colnames(seu))]
    if (cond_field %in% colnames(seu@meta.data)) {
      seu_full@meta.data[[cond_field]] <- seu@meta.data[[cond_field]][match(colnames(seu_full), colnames(seu))]
    }
    
    # Transfer dimensionality reductions (UMAP, PCA) for visualization
    # Create new DimReduc objects to avoid assay mismatch warnings
    if ("umap" %in% names(seu@reductions)) {
      umap_coords <- seu[["umap"]]@cell.embeddings
      seu_full[["umap"]] <- CreateDimReducObject(
        embeddings = umap_coords,
        key = "UMAP_",
        assay = DefaultAssay(seu_full)
      )
      log("  Transferred UMAP reduction")
    }
    if ("pca" %in% names(seu@reductions)) {
      pca_coords <- seu[["pca"]]@cell.embeddings
      seu_full[["pca"]] <- CreateDimReducObject(
        embeddings = pca_coords,
        key = "PC_",
        assay = DefaultAssay(seu_full)
      )
      log("  Transferred PCA reduction")
    }
    
    # Use the full normalized object for analysis
    seu <- seu_full
    log("Using normalized Seurat with all genes for analysis")
  } else {
    log("WARNING: Cell names don't match between objects, using annotated object")
  }
} else {
  log("WARNING: No normalized Seurat provided, using clustered object with limited genes")
}

log(sprintf("Final Seurat for analysis: %d cells, %d genes", ncol(seu), nrow(seu)))

# Log gene ID format in Seurat
seurat_genes <- rownames(seu)
log("Gene ID format in Seurat:")
log("  First 10 genes:", paste(head(seurat_genes, 10), collapse=", "))
log("  Last 10 genes:", paste(tail(seurat_genes, 10), collapse=", "))
log("  Example gene check - are these gene symbols?")
log("    Contains 'ENSMUSG':", any(grepl("^ENSMUSG", seurat_genes)))
log("    Contains typical symbols (Actb, Gapdh, etc.):", any(seurat_genes %in% c("Actb", "Gapdh", "Cd3e", "Cd4", "Cd8a")))

Idents(seu) <- "cell_type"
log("Using condition field:", cond_field)

conditions <- unique(seu@meta.data[[cond_field]])
log("Found", length(conditions), "conditions:", paste(conditions, collapse = ", "))

################################################################################
# LOAD DE GENES FROM STEP 03
################################################################################

log("Loading differential expression results from Step 03...")

ds_results <- NULL
da_dir_candidates <- c(
  file.path(dirname(dirname(opt$seurat)), "03_differential"),
  file.path(dirname(opt$seurat), "..", "03_differential"),
  file.path(dirname(base_dir), "03_differential"),
  file.path(dirname(dirname(base_dir)), "03_differential"),
  file.path(pipeline_output_dir, "03_differential")
)

da_dir <- NULL
for (candidate in da_dir_candidates) {
  if (dir.exists(candidate)) {
    da_dir <- candidate
    log("Found differential analysis results at:", da_dir)
    break
  }
}

if (!is.null(da_dir)) {
  ds_file <- file.path(da_dir, "data", "DS_seurat_results.csv")
  if (!file.exists(ds_file)) {
    alt_ds <- list.files(da_dir, pattern = "^DS_seurat_results\\.csv$", recursive = TRUE, full.names = TRUE)
    if (length(alt_ds) > 0) {
      ds_file <- alt_ds[1]
      log("Found nested DE results at:", ds_file)
    }
  }
  if (file.exists(ds_file)) {
    ds_results <- read.csv(ds_file, stringsAsFactors = FALSE)
    
    if ("cluster" %in% colnames(ds_results)) {
      de_genes_by_cluster <- list()
      for (clust in unique(ds_results$cluster)) {
        cluster_genes <- ds_results %>%
          filter(cluster == clust, p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
          pull(gene) %>%
          unique()
        de_genes_by_cluster[[as.character(clust)]] <- cluster_genes
      }
      log("Loaded DE genes for", length(de_genes_by_cluster), "clusters from Step 03")
      
      total_de <- ds_results %>%
        filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
        pull(gene) %>%
        unique() %>%
        length()
      log("Total unique DE genes:", total_de)
    }
  }
}

################################################################################
# LOAD NICHENET REFERENCE DATA
################################################################################

log("Loading NicheNet reference files...")

lr_network_file <- NULL
ligand_target_matrix_file <- NULL

if (file.exists(file.path(opt$nichenet_dir, "lr_network.rds"))) {
  lr_network_file <- file.path(opt$nichenet_dir, "lr_network.rds")
} else if (file.exists(file.path(opt$nichenet_dir, "lr_network_mouse.rds"))) {
  lr_network_file <- file.path(opt$nichenet_dir, "lr_network_mouse.rds")
} else if (file.exists(file.path(opt$nichenet_dir, "lr_network_human.rds"))) {
  lr_network_file <- file.path(opt$nichenet_dir, "lr_network_human.rds")
}

ligand_target_matrix_file <- file.path(opt$nichenet_dir, "ligand_target_matrix.rds")

if (is.null(lr_network_file) || !file.exists(lr_network_file)) {
  stop("NicheNet lr_network file not found in: ", opt$nichenet_dir)
}

if (!file.exists(ligand_target_matrix_file)) {
  stop("NicheNet ligand_target_matrix file not found in: ", opt$nichenet_dir)
}

lr_network <- readRDS(lr_network_file)
ligand_target_matrix <- readRDS(ligand_target_matrix_file)

# Convert gene symbols to mouse format if needed (First letter uppercase, rest lowercase)
# Check if genes are in human format (all uppercase)
sample_genes_lr <- head(unique(c(lr_network$from, lr_network$to)), 100)
if (mean(grepl("^[A-Z0-9]+$", sample_genes_lr)) > 0.5) {
  log("  Detected human gene format in LR network (e.g., UPPERCASE). Converting to mouse format...")
  
  # Convert function: First letter uppercase, rest lowercase (e.g., CXCR2 -> Cxcr2)
  convert_to_mouse <- function(gene) {
    paste0(toupper(substr(gene, 1, 1)), tolower(substr(gene, 2, nchar(gene))))
  }
  
  lr_network$from <- sapply(lr_network$from, convert_to_mouse, USE.NAMES = FALSE)
  lr_network$to <- sapply(lr_network$to, convert_to_mouse, USE.NAMES = FALSE)
  
  log("  Converted LR network gene symbols to mouse format")
  log("  Example conversion: CXCR2 -> Cxcr2, CCR6 -> Ccr6")
  
  # CRITICAL: Also convert ligand-target matrix row/column names
  log("  Checking ligand-target matrix gene format...")
  sample_matrix_genes <- head(c(rownames(ligand_target_matrix), colnames(ligand_target_matrix)), 100)
  if (mean(grepl("^[A-Z0-9]+$", sample_matrix_genes)) > 0.5) {
    log("  Converting ligand-target matrix gene names to mouse format...")
    rownames(ligand_target_matrix) <- sapply(rownames(ligand_target_matrix), convert_to_mouse, USE.NAMES = FALSE)
    colnames(ligand_target_matrix) <- sapply(colnames(ligand_target_matrix), convert_to_mouse, USE.NAMES = FALSE)
    log("  Converted matrix: ", nrow(ligand_target_matrix), " targets x ", ncol(ligand_target_matrix), " ligands")
  } else {
    log("  Matrix gene names already in correct format")
  }
}

log("Loaded NicheNet reference files:")
log("  - LR network: ", nrow(lr_network), " interactions")
log("  - Ligand-target matrix: ", nrow(ligand_target_matrix), " targets x ", 
    ncol(ligand_target_matrix), " ligands")

log("\nGene ID format in NicheNet ligand-target matrix:")
matrix_genes <- rownames(ligand_target_matrix)
log("  First 10 genes:", paste(head(matrix_genes, 10), collapse=", "))
log("  Last 10 genes:", paste(tail(matrix_genes, 10), collapse=", "))
log("  Contains 'ENSMUSG':", any(grepl("^ENSMUSG", matrix_genes)))
log("  Contains typical symbols:", any(matrix_genes %in% c("Actb", "Gapdh", "Cd3e", "Cd4", "Cd8a")))

log("\nGene ID overlap between Seurat and NicheNet matrix:")
overlap_genes <- intersect(seurat_genes, matrix_genes)
log("  Seurat genes:", length(seurat_genes))
log("  Matrix genes:", length(matrix_genes))
log("  Overlap:", length(overlap_genes), sprintf("(%.1f%% of Seurat genes)", 100 * length(overlap_genes) / length(seurat_genes)))
if (length(overlap_genes) > 0) {
  log("  Example overlapping genes:", paste(head(overlap_genes, 15), collapse=", "))
} else {
  log("  WARNING: NO OVERLAP! Gene ID mismatch detected!")
}

# Show non-matching examples
seurat_only <- setdiff(seurat_genes[1:20], matrix_genes)
matrix_only <- setdiff(matrix_genes[1:20], seurat_genes)
if (length(seurat_only) > 0) {
  log("  Example genes in Seurat but NOT in matrix:", paste(head(seurat_only, 10), collapse=", "))
}
if (length(matrix_only) > 0) {
  log("  Example genes in matrix but NOT in Seurat:", paste(head(matrix_only, 10), collapse=", "))
}

# Optional: Load weighted networks
weighted_networks_file <- file.path(opt$nichenet_dir, "weighted_networks.rds")
if (file.exists(weighted_networks_file)) {
  weighted_networks <- readRDS(weighted_networks_file)
  log("  - Weighted networks loaded")
} else {
  weighted_networks <- NULL
  log("  - Weighted networks not found (optional)")
}

################################################################################
# NICHENET ANALYSIS
################################################################################

log("Creating NicheNet overview plots...")
create_nichenet_overview_celltypes(seu, plots_dir, CFG_BASE_SIZE, CFG_TITLE_SIZE, CFG_LEGEND_SIZE)
create_nichenet_overview_conditions(seu, plots_dir, CFG_BASE_SIZE, CFG_TITLE_SIZE, CFG_LEGEND_SIZE)

cell_types <- levels(Idents(seu))
all_ligand_activities <- list()
all_ligand_activities_by_condition <- list()  # Store condition-specific activities

skip_nichenet <- FALSE
if (is.null(ds_results) || !"cluster" %in% colnames(ds_results)) {
  log("WARNING: No differential expression results from Step 03 found.")
  log("         NicheNet requires DE genes comparing conditions within cell types.")
  log("         Skipping NicheNet analysis.")
  skip_nichenet <- TRUE
}

# Run NicheNet per receiver cell type (following official vignette approach)
# Note: Uses ALL conditions for each receiver (NicheNet compares across conditions)
if (!skip_nichenet) {
  log("Starting NicheNet ligand activity analysis...")
  log("Approach: Sender-agnostic (as per official vignette)")
  log("")
  
  for (receiver in cell_types) {
    log("=== Analyzing receiver:", receiver, "===")
    
    # Use all cells from all conditions for this receiver (as per vignette)
    seu_receiver <- subset(seu, idents = receiver)
    
    if (ncol(seu_receiver) < 20) {
      log("  Too few cells for", receiver, "across all conditions (n=", ncol(seu_receiver), ") - skipping.")
      next
    }
    
    log("  Receiver cells:", ncol(seu_receiver), "across", length(conditions), "conditions")
    
    # Get DE genes for this receiver from Step 03
    # Take all significant genes from any condition comparison for this cluster
    cluster_de <- ds_results %>% 
      filter(cluster == receiver)
    
    if (nrow(cluster_de) == 0) {
      log("  No DE genes found in Step 03 for", receiver, "- skipping.")
      log("  (DE genes should compare conditions within this cell type)")
      next
    }
    
    # Get genes that are significant in at least one comparison
    # Use the 'significant' column if available, otherwise filter by thresholds
    if ("significant" %in% colnames(cluster_de)) {
      geneset_oi <- cluster_de %>% 
        filter(significant == TRUE) %>%
        pull(gene) %>% 
        unique()
      log("  Significant DE genes (any comparison):", length(geneset_oi))
      log("  Example DE genes:", paste(head(geneset_oi, 10), collapse=", "))
    } else {
      geneset_oi <- cluster_de %>% 
        filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
        pull(gene) %>% 
        unique()
      log("  Significant DE genes (p_adj<0.05, |LFC|>0.25):", length(geneset_oi))
      log("  Example DE genes:", paste(head(geneset_oi, 10), collapse=", "))
    }
    
    # Filter to genes in ligand-target matrix
    geneset_oi_original <- geneset_oi
    geneset_oi <- geneset_oi[geneset_oi %in% rownames(ligand_target_matrix)]
    log("  Gene set after filtering to ligand-target matrix:", length(geneset_oi))
    log("  Genes lost in filtering:", length(geneset_oi_original) - length(geneset_oi))
    
    if (length(geneset_oi) < length(geneset_oi_original)) {
      genes_lost <- setdiff(geneset_oi_original, geneset_oi)
      log("  Examples of lost genes:", paste(head(genes_lost, 10), collapse=", "))
    }
    if (length(geneset_oi) > 0) {
      log("  Examples of retained genes:", paste(head(geneset_oi, 10), collapse=", "))
    }
    
    # If too few genes remain, try using top genes by p-value regardless of significance flag
    if (length(geneset_oi) < 10) {
      log("  Too few genes after matrix filtering. Trying top genes by p-value...")
      geneset_oi <- cluster_de %>%
        arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
        head(200) %>%
        pull(gene) %>%
        unique()
      geneset_oi <- geneset_oi[geneset_oi %in% rownames(ligand_target_matrix)]
      log("  Gene set after taking top 200 by p-value:", length(geneset_oi))
    }
    
    log("  Final gene set of interest (DE genes):", length(geneset_oi))
    
    if (length(geneset_oi) < 5) {
      log("  Too few DE genes (", length(geneset_oi), ") - skipping. Minimum: 5")
      log("  This may indicate gene ID mismatch between Seurat and NicheNet matrix")
      next
    }
    
    # Define background: all expressed genes in receiver across all conditions
    # Use normalized data layer for expression check
    expr_data <- GetAssayData(seu_receiver, slot = "data", assay = "RNA")
    expressed_genes_receiver <- rownames(expr_data)[Matrix::rowSums(expr_data > 0) > ncol(seu_receiver) * 0.05]
    
    log("  Total genes in receiver Seurat:", nrow(expr_data))
    log("  Expressed genes (>5% cells):", length(expressed_genes_receiver))
    log("  Example expressed genes:", paste(head(expressed_genes_receiver, 10), collapse=", "))
    
    # Background genes (as per vignette: expressed genes that are in ligand_target_matrix)
    background_expressed_genes <- expressed_genes_receiver[expressed_genes_receiver %in% rownames(ligand_target_matrix)]
    log("  Background genes:", length(background_expressed_genes))
    log("  Ratio (background:geneset):", round(length(background_expressed_genes) / length(geneset_oi), 1), "x")
    
    # Check reasonable size
    if (length(background_expressed_genes) < 1000) {
      log("  Warning: Small background (", length(background_expressed_genes), "). Recommended: 5000-10000")
    }
    
    if (length(background_expressed_genes) < length(geneset_oi) * 3) {
      log("  Warning: Background too small relative to geneset (recommended: >5x)")
    }
    
    # Define expressed receptors
    all_receptors <- unique(lr_network$to)
    expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
    log("  Total receptors in LR network:", length(all_receptors))
    log("  Expressed receptors:", length(expressed_receptors))
    
    if (length(expressed_receptors) > 0) {
      log("  Example expressed receptors:", paste(head(expressed_receptors, min(10, length(expressed_receptors))), collapse=", "))
    } else {
      log("  WARNING: No receptors found! Checking receptor-gene overlap...")
      log("    Receptors in LR network (first 10):", paste(head(all_receptors, 10), collapse=", "))
      log("    Expressed genes in receiver (first 10):", paste(head(expressed_genes_receiver, 10), collapse=", "))
      log("    Check for gene ID format mismatch!")
    }
    
    # Define potential ligands (sender-agnostic: all ligands whose receptors are expressed)
    potential_ligands <- lr_network %>% 
      filter(to %in% expressed_receptors) %>% 
      pull(from) %>% 
      unique()
    log("  Potential ligands (sender-agnostic):", length(potential_ligands))
    
    # CRITICAL: Check overlaps before prediction
    log("  Checking input overlaps for predict_ligand_activities...")
    log("    Ligand-target matrix dimensions:", nrow(ligand_target_matrix), "targets x", ncol(ligand_target_matrix), "ligands")
    log("    Geneset size:", length(geneset_oi))
    log("    Background size:", length(background_expressed_genes))
    log("    Potential ligands:", length(potential_ligands))
    
    # Check geneset overlap with matrix rows (targets)
    geneset_in_matrix <- intersect(geneset_oi, rownames(ligand_target_matrix))
    log("    Geneset overlap with matrix targets:", length(geneset_in_matrix), "/", length(geneset_oi))
    
    # Check background overlap with matrix rows
    background_in_matrix <- intersect(background_expressed_genes, rownames(ligand_target_matrix))
    log("    Background overlap with matrix targets:", length(background_in_matrix), "/", length(background_expressed_genes))
    
    # Check potential_ligands overlap with matrix columns (ligands)
    ligands_in_matrix <- intersect(potential_ligands, colnames(ligand_target_matrix))
    log("    Potential ligands overlap with matrix ligands:", length(ligands_in_matrix), "/", length(potential_ligands))
    
    if (length(geneset_in_matrix) < 5) {
      log("  ERROR: Geneset has insufficient overlap with ligand-target matrix (<5 genes)")
      next
    }
    if (length(ligands_in_matrix) < 1) {
      log("  ERROR: No potential ligands found in ligand-target matrix columns")
      log("    Example potential ligands:", paste(head(potential_ligands, 10), collapse=", "))
      log("    Example matrix ligands (first 10):", paste(head(colnames(ligand_target_matrix), 10), collapse=", "))
      next
    }
    
    # CRITICAL: Use only ligands that are in the matrix
    potential_ligands_filtered <- ligands_in_matrix
    log("    Using", length(potential_ligands_filtered), "ligands present in matrix for prediction")
    
    # Ligand activity prediction
    log("  Predicting ligand activities...")
    log("    Calling predict_ligand_activities with:")
    log("      - geneset: ", length(geneset_oi), " genes")
    log("      - background: ", length(background_expressed_genes), " genes")
    log("      - potential_ligands: ", length(potential_ligands_filtered), " ligands")
    log("      - matrix: ", nrow(ligand_target_matrix), " x ", ncol(ligand_target_matrix))
    
    ligand_activities <- tryCatch({
      result <- predict_ligand_activities(
        geneset = geneset_oi, 
        background_expressed_genes = background_expressed_genes, 
        ligand_target_matrix = ligand_target_matrix, 
        potential_ligands = potential_ligands_filtered
      )
      log("    Result class:", paste(class(result), collapse=", "))
      if (!is.null(result) && nrow(result) > 0) {
        log("    Result dimensions:", nrow(result), "rows x", ncol(result), "cols")
        log("    Column names:", paste(colnames(result), collapse=", "))
        log("    SUCCESS: Returning data.frame with", nrow(result), "ligands")
      } else {
        log("    WARNING: Result is NULL or empty!")
      }
      result  # Return as last expression, not with return()
    }, error = function(e) {
      log("  ERROR in predict_ligand_activities:", e$message)
      NULL
    })
    
    if (is.null(ligand_activities)) {
      log("  ERROR: ligand_activities is NULL after tryCatch - skipping.")
      next
    }
    if (nrow(ligand_activities) == 0) {
      log("  ERROR: ligand_activities has 0 rows - skipping.")
      next
    }
    
    ligand_activities <- ligand_activities %>% arrange(-aupr_corrected)
    log("  Found", nrow(ligand_activities), "ligands. Top AUPR:", 
        round(ligand_activities$aupr_corrected[1], 3))
    
    write.csv(ligand_activities, 
              file.path(data_dir, paste0("nichenet_ligands_", receiver, ".csv")), 
              row.names = FALSE)
    
    # Store for summary
    all_ligand_activities[[receiver]] <- ligand_activities %>% 
      mutate(receiver = receiver) %>%
      select(test_ligand, aupr_corrected, receiver)
    
    log("  Generating visualizations...")
    
    # Top 10 ligands for plotting
    top_ligands_10 <- ligand_activities %>% 
      top_n(10, aupr_corrected) %>% 
      arrange(-aupr_corrected) %>% 
      pull(test_ligand)
    
    ## === NICHENET PLOTS ===
    log("  Creating NicheNet plots...")
    
    # 1. Ligand activity plots per condition + combined figure
    log("  Creating condition-specific ligand activity plots...")
    ligand_activity_plots <- list()
    all_ligands_by_condition <- list()
    
    for (cond in conditions) {
      tryCatch({
        # Get condition-specific DE genes (if available in cluster_de)
        # Look for comparison columns like "control_vs_aPD-1"
        cond_comparisons <- grep(cond, colnames(cluster_de), value = TRUE)
        
        # Use all DE genes for now (condition comparisons)
        geneset_cond <- geneset_oi  # Same geneset across conditions
        
        # Get condition-specific expressed receptors
        receiver_cells_cond <- rownames(seu@meta.data)[
          seu@meta.data[[cond_field]] == cond & 
          seu@meta.data$cell_type == receiver
        ]
        
        if (length(receiver_cells_cond) >= 10 && length(geneset_cond) >= 10) {
          seu_receiver_cond <- subset(seu, cells = receiver_cells_cond)
          expr_data_cond <- GetAssayData(seu_receiver_cond, layer = "data")
          
          # Expressed genes in this condition
          expressed_genes_cond <- rownames(expr_data_cond)[
            Matrix::rowSums(expr_data_cond > 0) > ncol(expr_data_cond) * 0.05
          ]
          background_genes_cond <- expressed_genes_cond[expressed_genes_cond %in% rownames(ligand_target_matrix)]
          
          # Expressed receptors in this condition
          expressed_receptors_cond <- intersect(expressed_genes_cond, unique(lr_network$to))
          
          # Potential ligands in this condition
          potential_ligands_cond <- lr_network %>%
            filter(to %in% expressed_receptors_cond) %>%
            pull(from) %>%
            unique()
          potential_ligands_cond <- intersect(potential_ligands_cond, colnames(ligand_target_matrix))
          
          if (length(potential_ligands_cond) >= 5) {
            # Predict ligand activities for this condition
            ligand_activities_cond <- predict_ligand_activities(
              geneset = geneset_cond,
              background_expressed_genes = background_genes_cond,
              ligand_target_matrix = ligand_target_matrix,
              potential_ligands = potential_ligands_cond
            )
            
            if (!is.null(ligand_activities_cond) && nrow(ligand_activities_cond) > 0) {
              ligand_activities_cond <- ligand_activities_cond %>% arrange(-aupr_corrected)
              all_ligands_by_condition[[cond]] <- ligand_activities_cond
              
              # Store in global list with receiver info
              ligand_activities_cond_with_receiver <- ligand_activities_cond %>%
                mutate(receiver = receiver)
              
              if (is.null(all_ligand_activities_by_condition[[cond]])) {
                all_ligand_activities_by_condition[[cond]] <- list()
              }
              all_ligand_activities_by_condition[[cond]][[receiver]] <- ligand_activities_cond_with_receiver
              
              # Create plot for this condition
              top_10_cond <- head(ligand_activities_cond, 10)
              
              p <- ggplot(top_10_cond, aes(x = reorder(test_ligand, aupr_corrected), y = aupr_corrected)) +
                geom_col(fill = "#0072B2") +
                coord_flip() +
                labs(title = paste0(receiver, "\n", cond),
                     x = "Ligand", y = "AUPR") +
                theme_minimal(base_size = 12) +
                theme(
                  axis.text = element_text(size = 10),
                  axis.title = element_text(size = 12),
                  plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
                )
              
              ligand_activity_plots[[cond]] <- p
              log("    - Ligand activity for", cond, ":", nrow(ligand_activities_cond), "ligands predicted")
            }
          }
        }
      }, error = function(e) {
        log("    WARNING: Ligand activity for", cond, "failed:", e$message)
      })
    }
    
    # Combine all conditions into one figure
    if (length(ligand_activity_plots) > 0) {
      tryCatch({
        n_conds <- length(ligand_activity_plots)
        ncols <- min(2, n_conds)
        nrows <- ceiling(n_conds / ncols)
        
        combined_activity <- patchwork::wrap_plots(ligand_activity_plots, ncol = ncols)
        activity_file <- file.path(plots_dir, paste0("nichenet_ligand_activity_", 
                                                      gsub(" ", "_", receiver), ".png"))
        ggsave(activity_file, combined_activity,
               width = 8 * ncols,
               height = 6 * nrows,
               dpi = 300, bg = "white")
        log("    - Combined ligand activity saved:", basename(activity_file))
        log("      Contains", n_conds, "conditions in", nrows, "×", ncols, "layout")
      }, error = function(e) {
        log("    WARNING: Combined ligand activity plot failed:", e$message)
      })
    }
    
    # Use global top 10 ligands for downstream analyses
    top_ligands_10 <- ligand_activities %>% 
      top_n(10, aupr_corrected) %>% 
      arrange(-aupr_corrected) %>% 
      pull(test_ligand)
    
    # 2. LR network per condition + combined figure
    log("  Creating condition-specific LR network plots...")
    
    # First pass: Collect all expressed receptors across all conditions
    all_expressed_receptors_union <- c()
    for (cond in conditions) {
      receiver_cells_cond <- rownames(seu@meta.data)[
        seu@meta.data[[cond_field]] == cond & 
        seu@meta.data$cell_type == receiver
      ]
      if (length(receiver_cells_cond) >= 10) {
        seu_receiver_cond <- subset(seu, cells = receiver_cells_cond)
        expr_data_cond <- GetAssayData(seu_receiver_cond, layer = "data")
        expressed_receptors_cond <- rownames(expr_data_cond)[
          Matrix::rowSums(expr_data_cond > 0) > ncol(expr_data_cond) * 0.05
        ]
        expressed_receptors_cond <- intersect(expressed_receptors_cond, unique(lr_network$to))
        all_expressed_receptors_union <- union(all_expressed_receptors_union, expressed_receptors_cond)
      }
    }
    
    # Filter to receptors that interact with top ligands
    lr_network_all <- lr_network %>%
      filter(from %in% top_ligands_10 & to %in% all_expressed_receptors_union)
    all_receptors_in_network <- unique(lr_network_all$to)
    
    log("  Using", length(all_receptors_in_network), "receptors across all conditions for consistent x-axis")
    
    # Second pass: Create plots with consistent x-axis
    lr_plots_by_condition <- list()
    
    for (cond in conditions) {
      tryCatch({
        # Subset receiver cells for this condition
        receiver_cells_cond <- rownames(seu@meta.data)[
          seu@meta.data[[cond_field]] == cond & 
          seu@meta.data$cell_type == receiver
        ]
        
        if (length(receiver_cells_cond) >= 10) {
          # Get expressed receptors in this condition
          seu_receiver_cond <- subset(seu, cells = receiver_cells_cond)
          expr_data_cond <- GetAssayData(seu_receiver_cond, layer = "data")
          expressed_receptors_cond <- rownames(expr_data_cond)[
            Matrix::rowSums(expr_data_cond > 0) > ncol(expr_data_cond) * 0.05
          ]
          expressed_receptors_cond <- intersect(expressed_receptors_cond, unique(lr_network$to))
          
          # Create LR network for this condition using consistent receptor set
          lr_network_subset <- lr_network %>%
            filter(from %in% top_ligands_10 & to %in% expressed_receptors_cond & to %in% all_receptors_in_network)
          
          if (nrow(lr_network_subset) > 0 || TRUE) {  # Always create plot for consistent layout
            # Create plot with all receptors on x-axis, but only show tiles for expressed ones
            vis_lr_df <- lr_network_subset %>%
              group_by(from, to) %>%
              summarize(val = 1, .groups = "drop")
            
            # Create complete grid with all combinations, then mark which are present
            all_combinations <- expand.grid(
              from = top_ligands_10,
              to = all_receptors_in_network,
              stringsAsFactors = FALSE
            )
            all_combinations <- all_combinations %>%
              left_join(vis_lr_df, by = c("from", "to")) %>%
              mutate(val = ifelse(is.na(val), 0, val))
            
            p <- ggplot(all_combinations, aes(x = to, y = from, fill = factor(val))) +
              geom_tile(color = "white", linewidth = 0.5) +
              scale_fill_manual(values = c("0" = "grey95", "1" = "#D55E00"), guide = "none") +
              labs(title = paste0(receiver, "\n", cond), 
                   x = "Receptor", y = "Ligand") +
              theme_minimal(base_size = 12) +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
                axis.text.y = element_text(size = 10),
                axis.title = element_text(size = 12),
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                panel.grid = element_blank()
              )
            
            lr_plots_by_condition[[cond]] <- p
            log("    - LR network for", cond, ":", sum(all_combinations$val == 1), "interactions")
          }
        }
      }, error = function(e) {
        log("    WARNING: LR network for", cond, "failed:", e$message)
      })
    }
    
    # Combine all conditions into one figure
    if (length(lr_plots_by_condition) > 0) {
      tryCatch({
        n_conds <- length(lr_plots_by_condition)
        ncols <- min(2, n_conds)
        nrows <- ceiling(n_conds / ncols)
        
        combined_lr <- patchwork::wrap_plots(lr_plots_by_condition, ncol = ncols)
        combined_file <- file.path(plots_dir, paste0("nichenet_lr_network_", 
                                                      gsub(" ", "_", receiver), ".png"))
        ggsave(combined_file, combined_lr, 
               width = 8 * ncols, 
               height = 6 * nrows, 
               dpi = 300, bg = "white")
        log("    - Combined LR network saved:", basename(combined_file))
        log("      Contains", n_conds, "conditions in", nrows, "×", ncols, "layout")
      }, error = function(e) {
        log("    WARNING: Combined LR network failed:", e$message)
      })
    }
    
    # 3. Target gene network (top 50 targets from top 10 ligands)
    tryCatch({
      ligand_target_df <- as.data.frame(ligand_target_matrix[geneset_oi, top_ligands_10])
      active_targets <- geneset_oi[apply(ligand_target_df, 1, max) > 0.2]
      if (length(active_targets) > 0) {
        create_nichenet_target_network(active_targets[1:min(50, length(active_targets))], 
                                        receiver, plots_dir, base_size = CFG_BASE_SIZE)
        log("    - Target network plot saved (top 50 targets)")
      }
    }, error = function(e) {
      log("    WARNING: Target network plot failed:", e$message)
    })
    
    log("  Completed", receiver)
    log("")
  }
  
  # Generate summary visualizations
  if (length(all_ligand_activities) > 0) {
    log("Generating summary visualizations...")
    
    # Create condition-specific summary heatmaps
    log("  Creating condition-specific summary heatmaps...")
    log("    Available conditions:", paste(conditions, collapse=", "))
    log("    all_ligand_activities_by_condition structure:")
    for (c in names(all_ligand_activities_by_condition)) {
      log("      -", c, ":", length(all_ligand_activities_by_condition[[c]]), "receivers")
    }
    summary_heatmap_plots <- list()
    
    for (cond in conditions) {
      log("    Processing condition:", cond)
      if (!is.null(all_ligand_activities_by_condition[[cond]]) && 
          length(all_ligand_activities_by_condition[[cond]]) > 0) {
        log("      Condition", cond, "has", length(all_ligand_activities_by_condition[[cond]]), "receivers")
        tryCatch({
          # Combine all receivers for this condition
          log("        Attempting to bind_rows for condition", cond)
          combined_activities_cond <- bind_rows(all_ligand_activities_by_condition[[cond]])
          log("        Combined activities has", nrow(combined_activities_cond), "rows")
          
          if (nrow(combined_activities_cond) > 0) {
            log("        Getting top 10 ligands...")
            # Get top 10 ligands across all receivers in this condition
            top_ligands_cond <- combined_activities_cond %>%
              group_by(test_ligand) %>%
              summarize(mean_aupr = mean(aupr_corrected), .groups = "drop") %>%
              arrange(-mean_aupr) %>%
              head(10) %>%
              pull(test_ligand)
            
            log("        Top 10 ligands:", paste(head(top_ligands_cond, 3), collapse=", "), "...")
            
            # Create matrix for heatmap
            # Handle potential duplicates by taking max AUPR per ligand-receiver pair
            heatmap_data_cond <- combined_activities_cond %>%
              filter(test_ligand %in% top_ligands_cond) %>%
              select(test_ligand, receiver, aupr_corrected) %>%
              group_by(test_ligand, receiver) %>%
              summarize(aupr_corrected = max(aupr_corrected, na.rm = TRUE), .groups = "drop") %>%
              pivot_wider(names_from = receiver, values_from = aupr_corrected, values_fill = 0)
            
            log("        Heatmap data dimensions:", nrow(heatmap_data_cond), "×", ncol(heatmap_data_cond) - 1)
            
            # Convert to matrix
            heatmap_matrix_cond <- as.matrix(heatmap_data_cond[, -1])
            rownames(heatmap_matrix_cond) <- heatmap_data_cond$test_ligand
            
            # Create heatmap as ggplot for combining
            heatmap_long <- heatmap_data_cond %>%
              pivot_longer(-test_ligand, names_to = "receiver", values_to = "aupr")
            
            p <- ggplot(heatmap_long, aes(x = receiver, y = test_ligand, fill = aupr)) +
              geom_tile(color = "white") +
              scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                                   midpoint = 0.5, limits = c(0, 1),
                                   name = "AUPR") +
              labs(title = cond, x = "Receiver", y = "Ligand") +
              theme_minimal(base_size = 12) +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                axis.text.y = element_text(size = 10),
                axis.title = element_text(size = 12),
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
              )
            
            summary_heatmap_plots[[cond]] <- p
            
            # Save individual condition-specific heatmap
            safe_cond_name <- gsub("[^A-Za-z0-9_-]", "_", cond)
            cond_summary_file <- file.path(plots_dir, sprintf("nichenet_summary_heatmap_%s.png", safe_cond_name))
            ggsave(cond_summary_file, p,
                   width = 8,
                   height = 6,
                   dpi = 300, bg = "white")
            
            log("    - Summary heatmap for", cond, ":", nrow(heatmap_matrix_cond), "ligands ×", 
                ncol(heatmap_matrix_cond), "receivers")
            log("      Saved:", basename(cond_summary_file))
          } else {
            log("        Skipping", cond, "- combined_activities_cond has 0 rows")
          }
        }, error = function(e) {
          log("    WARNING: Summary heatmap for", cond, "failed:", e$message)
          log("      Traceback:", paste(as.character(sys.calls()), collapse=" -> "))
        })
      } else {
        log("      Condition", cond, "skipped - no data available")
        if (is.null(all_ligand_activities_by_condition[[cond]])) {
          log("        Reason: all_ligand_activities_by_condition[[", cond, "]] is NULL")
        } else if (length(all_ligand_activities_by_condition[[cond]]) == 0) {
          log("        Reason: all_ligand_activities_by_condition[[", cond, "]] is empty")
        }
      }
    }
    
    # Combine all condition-specific summary heatmaps
    if (length(summary_heatmap_plots) > 0) {
      tryCatch({
        n_conds <- length(summary_heatmap_plots)
        ncols <- min(2, n_conds)
        nrows <- ceiling(n_conds / ncols)
        
        combined_summary <- patchwork::wrap_plots(summary_heatmap_plots, ncol = ncols)
        summary_file <- file.path(plots_dir, "nichenet_summary_heatmap.png")
        ggsave(summary_file, combined_summary,
               width = 8 * ncols,
               height = 6 * nrows,
               dpi = 300, bg = "white")
        log("  - Combined summary heatmap saved:", basename(summary_file))
        log("    Contains", n_conds, "conditions in", nrows, "×", ncols, "layout")
      }, error = function(e) {
        log("  WARNING: Combined summary heatmap failed:", e$message)
      })
    }
    
    # Original summary heatmap (all conditions combined into single heatmap)
    combined_activities <- bind_rows(all_ligand_activities)
    create_nichenet_summary_heatmap(
      combined_activities = combined_activities,
      plots_dir = plots_dir,
      n_top = 10,
      base_size = CFG_BASE_SIZE,
      title_size = CFG_TITLE_SIZE,
      axis_title_size = CFG_AXIS_TITLE_SIZE,
      legend_size = CFG_LEGEND_SIZE,
      condition = "all_combined"
    )
    log("  - Overall summary heatmap (all conditions combined) saved as: nichenet_summary_heatmap_all_combined.png")
    
    # LR interaction heatmap
    create_nichenet_lr_interaction_heatmap(
      seu = seu,
      lr_network = lr_network,
      all_ligand_activities = all_ligand_activities,
      plots_dir = plots_dir,
      n_top = 10,
      base_size = CFG_BASE_SIZE
    )
    log("  - LR interaction heatmap")
    
    log("NicheNet analysis completed successfully!")
    log("Analyzed", length(all_ligand_activities), "receiver cell types")
  } else {
    log("No receivers analyzed (insufficient DE genes or cells)")
  }
  
} else {
  log("NicheNet analysis skipped (no DE results from Step 03)")
}

################################################################################
# COMPLETION
################################################################################

saveRDS(list(
  status = "completed", 
  time = Sys.time(), 
  step = "nichenet",
  n_receivers_analyzed = length(all_ligand_activities)
), file.path(data_dir, "nichenet_summary.rds"))

log("NicheNet analysis script completed.")

end_performance_tracking(
  perf_tracker,
  success = TRUE,
  additional_metrics = list(
    n_cells = ncol(seu),
    n_cell_types = length(unique(seu$cell_type)),
    n_conditions = length(conditions),
    n_receivers_analyzed = length(all_ligand_activities)
  )
)

#!/usr/bin/env Rscript
################################################################################
# 04a_cellchat.R
#
# CellChat Analysis - Cell-cell communication inference
# 
# This script performs:
# - Comparison mode: Differential signaling between conditions
# - Single mode: Communication network for single condition
#
# Required inputs:
# - Clustered Seurat object with cell type annotations
#
# Outputs:
# - CellChat objects and interaction networks
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
  library(ggplot2)
  library(CellChat)
  library(patchwork)
  library(ComplexHeatmap)
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
source(file.path(script_dir, "utils", "pvalue_helpers.R"))

# Font configuration
CFG_BASE_SIZE <- 18
CFG_TITLE_SIZE <- 24
CFG_AXIS_TITLE_SIZE <- 20

set.seed(CFG_RANDOM_SEED)

################################################################################
# ARGUMENT PARSING
################################################################################
option_list <- list(
  make_option("--seurat", type="character", default=NULL, help="Seurat object (.rds)"),
  make_option("--out", type="character", default="04_cellcomm", help="Output directory"),
  make_option("--species", type="character", default="mouse", help="Species (mouse/human)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$seurat)) {
  print_help(opt_parser)
  stop("Input Seurat object is required (--seurat).", call.=FALSE)
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
log_file <- file.path(base_dir, "04a_cellchat.log")
setup_logging(log_file)

# Start performance tracking
perf_tracker <- start_performance_tracking("04a_cellchat", base_dir)

log_section("CellChat Analysis")
log("Input Seurat:", opt$seurat)
log("Output Dir:", base_dir)
log("Species:", opt$species)

################################################################################
# LOAD DATA
################################################################################

seu <- readRDS(opt$seurat)
log(sprintf("Loaded Seurat object: %d cells, %d genes", ncol(seu), nrow(seu)))

# Ensure cell_type annotation exists
if (!"cell_type" %in% colnames(seu@meta.data)) {
  if ("seurat_clusters" %in% colnames(seu@meta.data)) {
    seu$cell_type <- paste0("Cluster_", seu$seurat_clusters)
  } else {
    stop("No cell_type or seurat_clusters found in Seurat object.")
  }
}
Idents(seu) <- "cell_type"

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
log("Using condition field:", cond_field)

conditions <- unique(seu@meta.data[[cond_field]])
log("Found", length(conditions), "conditions:", paste(conditions, collapse = ", "))

################################################################################
# LOAD DE GENES FROM STEP 03 (for metadata in plots)
################################################################################

de_genes_all <- c()
de_source <- "None"

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
    if (length(alt_ds) > 0) ds_file <- alt_ds[1]
  }
  if (file.exists(ds_file)) {
    ds_results <- read.csv(ds_file, stringsAsFactors = FALSE)
    de_genes_all <- ds_results %>%
      filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
      pull(gene) %>%
      unique()
    log("Total", length(de_genes_all), "unique DE genes from Step 03")
    de_source <- "Seurat-pseudobulk"
  }
}

################################################################################
# CELLCHAT COMPARISON ANALYSIS (MULTI-CONDITION)
################################################################################

if (length(conditions) >= 2) {
  log("Running CellChat comparison analysis between conditions...")
  
  cellchat_list <- list()
  
  for (cond in conditions) {
    log("Processing condition:", cond)
    
    seu_cond <- subset(seu, subset = !!sym(cond_field) == cond)
    cellchat <- createCellChat(object = seu_cond, group.by = "cell_type", assay = "RNA")
    
    if (tolower(opt$species) == "human") {
      CellChatDB <- CellChatDB.human
    } else {
      CellChatDB <- CellChatDB.mouse
    }
    cellchat@DB <- CellChatDB
    
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat, thresh.pc = 0.1, thresh.fc = 0.05)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    cellchat_list[[cond]] <- cellchat
    log("  -", cond, ":", nrow(cellchat@net$count), "cell groups,", 
            sum(cellchat@net$count > 0), "interactions")
  }

  create_cellchat_network_circles_by_condition(
    cellchat_list = cellchat_list,
    plots_dir = plots_dir,
    log = log
  )
  
  if (length(cellchat_list) >= 2) {
    log("Computing network centrality for comparison...")
    for (i in seq_along(cellchat_list)) {
      cellchat_list[[i]] <- netAnalysis_computeCentrality(cellchat_list[[i]], slot.name = "netP")
    }
    
    log("Merging CellChat objects for comparison...")
    cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
    
    # Update metadata to use actual condition names instead of "Dataset 1/2"
    condition_names <- names(cellchat_list)
    if ("datasets" %in% colnames(cellchat_merged@meta)) {
      # Map generic Dataset labels to actual condition names
      current_levels <- levels(factor(cellchat_merged@meta$datasets))
      if (any(grepl("^Dataset", current_levels))) {
        log("  Replacing generic 'Dataset' labels with actual condition names...")
        for (i in seq_along(condition_names)) {
          cellchat_merged@meta$datasets[cellchat_merged@meta$datasets == paste0("Dataset ", i)] <- condition_names[i]
        }
      }
      cellchat_merged@meta$datasets <- factor(cellchat_merged@meta$datasets, 
                                               levels = condition_names)
      log(sprintf("  Dataset labels: %s", paste(levels(cellchat_merged@meta$datasets), collapse = ", ")))
    }
    
    saveRDS(cellchat_merged, file.path(rds_dir, "cellchat_merged.rds"))
    saveRDS(cellchat_list, file.path(rds_dir, "cellchat_object_list.rds"))
    
    log("Creating CellChat comparison visualizations...")
    
    log("Inspecting merged CellChat object structure...")
    log(sprintf("  Number of datasets: %d", length(cellchat_merged@net)))
    for (i in seq_along(cellchat_merged@net)) {
      net_name <- names(cellchat_merged@net)[i]
      log(sprintf("  Dataset %d (%s):", i, net_name))
      log(sprintf("    - Count matrix dim: %s", paste(dim(cellchat_merged@net[[i]]$count), collapse=" x ")))
      log(sprintf("    - Weight matrix dim: %s", paste(dim(cellchat_merged@net[[i]]$weight), collapse=" x ")))
      log(sprintf("    - Total interactions (count): %d", sum(cellchat_merged@net[[i]]$count > 0)))
      log(sprintf("    - Total interactions (weight): %.2f", sum(cellchat_merged@net[[i]]$weight)))
    }
    
    create_cellchat_interaction_comparison(
      cellchat_merged = cellchat_merged,
      de_source = de_source,
      de_count = length(de_genes_all),
      plots_dir = plots_dir
    )
    log("  Saved total interactions comparison")
    
    create_cellchat_diff_circle(
      cellchat_merged = cellchat_merged,
      plots_dir = plots_dir,
      log = log
    )
    
    create_cellchat_diff_heatmap(
      cellchat_merged = cellchat_merged,
      plots_dir = plots_dir,
      log = log
    )

    create_cellchat_diff_heatmaps_all(
      cellchat_list = cellchat_list,
      plots_dir = plots_dir,
      log = log
    )
    
    create_cellchat_signaling_roles(
      cellchat_list = cellchat_list,
      de_source = de_source,
      de_count = length(de_genes_all),
      plots_dir = plots_dir
    )
    log("  Saved signaling role scatter plots")
    
    create_cellchat_functional_similarity(
      cellchat_merged = cellchat_merged,
      plots_dir = plots_dir,
      log = log
    )
    
    create_cellchat_pathway_ranking(
      cellchat_merged = cellchat_merged,
      de_source = de_source,
      de_count = length(de_genes_all),
      plots_dir = plots_dir,
      log = log
    )
    
    create_cellchat_outgoing_heatmap(
      cellchat_list = cellchat_list,
      plots_dir = plots_dir,
      log = log
    )
    
    log("CellChat comparison completed successfully.")
    
  } else {
    log("Less than 2 conditions processed successfully - skipping comparison plots")
  }
  
} else {
  log("Only 1 condition found - running single CellChat analysis")
}

################################################################################
# CELLCHAT SINGLE ANALYSIS (SINGLE CONDITION)
################################################################################

if (length(conditions) == 1) {
  log("Running single CellChat analysis...")
  
  cellchat <- createCellChat(object = seu, group.by = "cell_type", assay = "RNA")
  
  if (tolower(opt$species) == "human") {
    CellChatDB <- CellChatDB.human
  } else {
    CellChatDB <- CellChatDB.mouse
  }
  cellchat@DB <- CellChatDB
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  saveRDS(cellchat, file.path(rds_dir, "cellchat_object.rds"))
  
  df.net <- subsetCommunication(cellchat)
  write.csv(df.net, file.path(data_dir, "cellchat_interactions.csv"), row.names = FALSE)
  
  log("Generating CellChat visualizations...")
  
  create_cellchat_network_circles(cellchat, plots_dir)
  log("  Saved circle plots (count + weight).")

  create_cellchat_bubble_plot(cellchat, plots_dir)
  log("  Saved bubble plot.")

  create_cellchat_heatmap_count(cellchat, plots_dir)
  log("  Saved count heatmap.")
  
  create_cellchat_heatmap_weight(cellchat, plots_dir)
  log("  Saved weight heatmap.")

  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  create_cellchat_centrality_heatmap(cellchat, plots_dir, n_pathways = 10, log = log)
  log("  Saved centrality heatmaps.")
  
  create_cellchat_signaling_aggregate(cellchat, plots_dir, n_pathways = 4)
  log("  Saved signaling aggregate plots.")
  
  log("CellChat single analysis completed successfully.")
}

################################################################################
# COMPLETION
################################################################################

saveRDS(list(status="completed", time=Sys.time(), step="cellchat"), 
        file.path(data_dir, "cellchat_summary.rds"))
log("CellChat analysis completed.")

end_performance_tracking(
  perf_tracker,
  success = TRUE,
  additional_metrics = list(
    n_cells = ncol(seu),
    n_cell_types = length(unique(seu$cell_type)),
    n_conditions = length(conditions)
  )
)

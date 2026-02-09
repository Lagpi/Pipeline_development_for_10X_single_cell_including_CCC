#!/usr/bin/env Rscript
#===============================================================================
# 06_gene_program_discovery_geneNMF.R
#
# This script uses GeneNMF to discover robust gene programs across samples
# and compares meta-program activities between conditions
#===============================================================================

################################################################################
# SETUP
################################################################################

if (!requireNamespace("GeneNMF", quietly = TRUE)) {
  message("Installing GeneNMF package...")
  tryCatch({
    install.packages("GeneNMF", repos = "https://cloud.r-project.org/")
    message("GeneNMF installed successfully")
  }, error = function(e) {
    message("Failed to install GeneNMF: ", e$message)
    stop("GeneNMF package is required. Install manually with: install.packages('GeneNMF')")
  })
}

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(RColorBrewer)
  library(tidyr)
  library(GeneNMF)
  library(plotly)
  library(htmlwidgets)
})

safe_lib <- function(pkgs) {
  for (p in pkgs) {
    suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE)))
  }
}
safe_lib(c("msigdbr", "fgsea"))

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
source(file.path(script_dir, "utils", "05_gene_program_plots.R"))
source(file.path(script_dir, "utils", "pvalue_helpers.R"))

# Font configuration (override defaults)
CFG_BASE_SIZE <- 18       # Base font size (was 16)
CFG_TITLE_SIZE <- 24      # Title font size (was 20)
CFG_AXIS_TITLE_SIZE <- 20 # Axis title size (was 18)

set.seed(CFG_RANDOM_SEED)

################################################################################
# ARGUMENT PARSING
################################################################################

option_list <- list(
  make_option("--seurat", type="character", default=NULL, 
              help="Seurat object (.rds) - can be integrated or normalized"),
  make_option("--normalized", type="character", default=NULL,
              help="Normalized Seurat object (.rds) - alternative to --seurat"),
  make_option("--annotated", type="character", default=NULL,
              help="Annotated Seurat object (.rds) - alternative to --seurat (with cell_type)"),
  make_option("--annotations", type="character", default=NULL,
              help="Annotations file (.rds) to add cell_type to normalized data"),
  make_option("--out", type="character", default="06_gene_programs", 
              help="Output directory"),
  make_option("--min_rank", type="integer", default=4, 
              help="Minimum number of programs (k) for NMF"),
  make_option("--max_rank", type="integer", default=9, 
              help="Maximum number of programs (k) for NMF"),
  make_option("--min_k", type="integer", default=NULL, 
              help="Alias for min_rank (deprecated)"),
  make_option("--max_k", type="integer", default=NULL, 
              help="Alias for max_rank (deprecated)"),
  make_option("--n_mp", type="integer", default=13, 
              help="Target number of meta-programs"),
  make_option("--min_exp", type="numeric", default=0.05, 
              help="Minimum expression fraction for genes"),
  make_option("--species", type="character", default="mouse", 
              help="Species (mouse/human) for GSEA"),
  make_option("--condition_field", type="character", default=NULL, 
              help="Metadata field for condition comparison")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Determine which Seurat object to use
# Priority: (normalized + annotations) > annotated > normalized > seurat
if (!is.null(opt$normalized) && !is.null(opt$annotations)) {
  seurat_path <- opt$normalized
  annotations_path <- opt$annotations
  use_combined <- TRUE
} else {
  seurat_path <- opt$annotated %||% opt$normalized %||% opt$seurat
  annotations_path <- NULL
  use_combined <- FALSE
}

if (is.null(seurat_path)) {
  print_help(opt_parser)
  stop("Input Seurat object is required (--seurat, --annotated, --normalized, or --normalized + --annotations).", call.=FALSE)
}

if (!is.null(opt$min_k) && is.null(opt$min_rank)) opt$min_rank <- opt$min_k
if (!is.null(opt$max_k) && is.null(opt$max_rank)) opt$max_rank <- opt$max_k

base_dir <- opt$out
dirs <- setup_analysis_directories(base_dir, "06_gene_programs")
data_dir <- dirs$data_dir
plots_dir <- dirs$plots_dir
dash_dir <- dirs$dash_dir
rds_dir <- dirs$rds_dir

log_file <- file.path(base_dir, "06_gene_program_discovery.log")
setup_logging(log_file)

# Start performance tracking
perf_tracker <- start_performance_tracking("06_gene_program_discovery", base_dir)

################################################################################
# FUNCTIONS
################################################################################

# Load Seurat object and prepare metadata fields for GeneNMF analysis
# Also integrates DE genes from Step 03 and ligand genes from Step 04
prepare_seurat_object <- function(seurat_path, condition_field = NULL, annotations_path = NULL) {
  log("Loading Seurat object from:", seurat_path)
  
  if (!file.exists(seurat_path)) {
    stop("Seurat file not found: ", seurat_path)
  }
  
  seu <- tryCatch({
    readRDS(seurat_path)
  }, error = function(e) {
    stop("Failed to load Seurat object: ", e$message)
  })
  
  if (!inherits(seu, "Seurat")) {
    stop("Loaded object is not a Seurat object")
  }
  
  # If annotations provided, merge cell_type metadata
  if (!is.null(annotations_path)) {
    log("Loading annotations from:", annotations_path)
    seu_anno <- readRDS(annotations_path)
    
    # Transfer cell_type metadata
    if ("cell_type" %in% colnames(seu_anno@meta.data)) {
      # Match cells by name and transfer cell_type
      shared_cells <- intersect(colnames(seu), colnames(seu_anno))
      
      if (length(shared_cells) == 0) {
        # Try to map via barcode extraction if direct matching fails
        log("No direct cell name match - attempting barcode-based mapping...")
        log("  Normalized cells:", ncol(seu))
        log("  Annotated cells:", ncol(seu_anno))
        log("  Example normalized:", head(colnames(seu), 2))
        log("  Example annotated:", head(colnames(seu_anno), 2))
        
        # Extract barcodes from both objects
        # HVG format: sample_BARCODE -> extract BARCODE
        # Clustered format: BARCODE_1 -> extract BARCODE
        seu_barcodes <- sub(".*_([ACGT]{16})$", "\\1", colnames(seu))
        anno_barcodes <- sub("^([ACGT]{16})_[0-9]+$", "\\1", colnames(seu_anno))
        
        # Create mapping table
        names(anno_barcodes) <- colnames(seu_anno)
        names(seu_barcodes) <- colnames(seu)
        
        # Match barcodes
        barcode_match <- match(seu_barcodes, anno_barcodes)
        matched_cells <- !is.na(barcode_match)
        
        if (sum(matched_cells) > 0) {
          log(sprintf("Successfully mapped %d/%d cells via barcode extraction", 
                      sum(matched_cells), ncol(seu)))
          
          # Transfer annotations using barcode mapping
          seu$cell_type <- NA
          matched_anno_names <- colnames(seu_anno)[barcode_match[matched_cells]]
          seu@meta.data[matched_cells, "cell_type"] <- 
            seu_anno@meta.data[matched_anno_names, "cell_type"]
          
          # Copy other important metadata if available
          for (meta_field in c("seurat_clusters", "data_driven_celltype", "condition", "sample")) {
            if (meta_field %in% colnames(seu_anno@meta.data) && 
                !meta_field %in% colnames(seu@meta.data)) {
              seu@meta.data[matched_cells, meta_field] <- 
                seu_anno@meta.data[matched_anno_names, meta_field]
              log(sprintf("  Transferred %s metadata", meta_field))
            }
          }
        } else {
          log("  First 3 cell names in normalized:", paste(head(colnames(seu), 3), collapse=", "))
          log("  First 3 cell names in annotated:", paste(head(colnames(seu_anno), 3), collapse=", "))
          stop("Cannot merge annotations - barcode extraction failed")
        }
      } else {
        # Direct matching worked
        seu$cell_type <- NA
        seu@meta.data[shared_cells, "cell_type"] <- seu_anno@meta.data[shared_cells, "cell_type"]
        log(sprintf("Transferred cell_type to %d cells (direct match)", length(shared_cells)))
      }
    } else {
      stop("No cell_type annotation found in annotations file")
    }
    
    log("Merged normalized data with annotations")
  }
  
  log("Successfully loaded Seurat object")
  log("  - Cells:", ncol(seu))
  log("  - Genes:", nrow(seu))
  
  # Ensure cell_type field exists in metadata (required for downstream analysis)
  if (!"cell_type" %in% colnames(seu@meta.data)) {
    if ("singler_main" %in% colnames(seu@meta.data)) {
      seu$cell_type <- seu$singler_main
    } else if ("seurat_clusters" %in% colnames(seu@meta.data)) {
      seu$cell_type <- paste0("Cluster_", seu$seurat_clusters)
    } else {
      stop("No cell_type, singler_main, or seurat_clusters found in Seurat object")
    }
  }
  
  # Auto-detect condition field if not provided (e.g., condition, treatment, group)
  if (is.null(condition_field)) {
    for (f in c("condition", "treatment", "group", "orig.ident")) {
      if (f %in% colnames(seu@meta.data) && length(unique(seu@meta.data[[f]])) > 1) {
        condition_field <- f
        break
      }
    }
  }
  
  if (is.null(condition_field)) {
    condition_field <- "orig.ident"
    log("Warning: No condition field detected, using orig.ident")
  }
  
  log("Condition field:", condition_field)
  
  # Detect sample field for splitting object into individual samples (required by multiNMF)
  sample_field <- NULL
  for (f in c("sample", "Sample", "patient", "orig.ident")) {
    if (f %in% colnames(seu@meta.data)) {
      sample_field <- f
      break
    }
  }
  
  if (is.null(sample_field)) {
    # Create artificial sample IDs if no sample field exists (assigns unique ID per condition)
    seu$geneNMF_sample <- paste0(seu@meta.data[[condition_field]], "_", 
                                  ave(seq_len(ncol(seu)), seu@meta.data[[condition_field]], 
                                      FUN = seq_along))
    sample_field <- "geneNMF_sample"
    log("Created artificial sample field:", sample_field)
  }
  
  log("Sample field:", sample_field)
  log("Number of samples:", length(unique(seu@meta.data[[sample_field]])))
  
  # ===========================================================================
  # LOAD DE GENES FROM STEP 03 (DIFFERENTIAL ANALYSIS)
  # Integrate differentially expressed genes to prioritize biologically relevant programs
  # ===========================================================================
  de_genes_all <- c()
  da_dir_candidates <- c(
    file.path(dirname(dirname(seurat_path)), "04_differential"),
    file.path(dirname(seurat_path), "..", "04_differential")
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
    if (file.exists(ds_file)) {
      ds_results <- read.csv(ds_file, stringsAsFactors = FALSE)
      de_genes_all <- ds_results %>%
        filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
        pull(gene) %>%
        unique()
      log("Loaded", length(de_genes_all), "DE genes from Step 03")
    }
  } else {
    log("Warning: No differential analysis results found")
  }
  
  # ===========================================================================
  # LOAD LIGAND GENES FROM STEP 04 (INTERCELLULAR COMMUNICATION)
  # Integrate top ligands from NicheNet to capture cell communication programs
  # ===========================================================================
  lr_genes_all <- c()
  comm_dir_candidates <- c(
    file.path(dirname(dirname(seurat_path)), "05_cellcomm"),
    file.path(dirname(seurat_path), "..", "05_cellcomm")
  )
  
  comm_dir <- NULL
  for (candidate in comm_dir_candidates) {
    if (dir.exists(candidate)) {
      comm_dir <- candidate
      log("Found cell communication results at:", comm_dir)
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
          log("Warning: Failed to load", basename(nf))
        })
      }
      lr_genes_all <- unique(unlist(lr_genes_list))
      log("Loaded", length(lr_genes_all), "ligand genes from Step 04 (NichNet)")
    }
  } else {
    log("Warning: No cell communication results found")
  }
  
  log("Integration summary:")
  log("  - DE genes from Step 03:", length(de_genes_all))
  log("  - Ligand genes from Step 04:", length(lr_genes_all))
  
  list(
    seu = seu,
    condition_field = condition_field,
    sample_field = sample_field,
    de_genes = de_genes_all,
    lr_genes = lr_genes_all
  )
}

#===============================================================================
# 2. RUN MULTINMF
# Perform NMF decomposition across multiple samples to identify recurrent gene programs
# Tests k=4 to k=9 ranks to find optimal number of programs per sample
# IMPORTANT: Uses log-normalized data from RNA assay with HVG genes only
# GeneNMF cannot handle negative values (e.g., from integrated assays)
#===============================================================================
run_multinmf <- function(seu, sample_field, min_k, max_k, min_exp, rds_dir) {
  log("Splitting object by sample...")
  # Use RNA assay (log-normalized) - GeneNMF cannot handle negative values from integrated data
  DefaultAssay(seu) <- "RNA"
  seu.list <- SplitObject(seu, split.by = sample_field)
  
  log("Running multiNMF across", length(seu.list), "samples...")
  log("Sample names:", paste(names(seu.list), collapse = ", "))
  log("k range:", min_k, "to", max_k)
  log("Min expression:", min_exp)
  log("This may take some time depending on dataset size...")
  
  # Check if samples have enough cells and genes
  for (s in names(seu.list)) {
    n_cells <- ncol(seu.list[[s]])
    n_genes <- nrow(seu.list[[s]])
    log(sprintf("  Sample %s: %d cells, %d genes", s, n_cells, n_genes))
    if (n_cells < 50) {
      log(sprintf("  WARNING: Sample %s has only %d cells - may be too few for NMF", s, n_cells))
    }
  }
  
  # Determine number of features to use for NMF (max 3000 highly variable genes)
  # Using variable features reduces computational cost and noise
  var_features <- VariableFeatures(seu)
  if (length(var_features) > 0) {
    nfeatures_to_use <- min(3000, length(var_features))
    log("Using", nfeatures_to_use, "variable features for NMF")
  } else {
    nfeatures_to_use <- min(3000, nrow(seu))
    log("No variable features set, using top", nfeatures_to_use, "genes")
  }
  
  geneNMF.programs <- tryCatch({
    withCallingHandlers({
      multiNMF(seu.list, 
               assay = "RNA", 
               k = min_k:max_k, 
               min.exp = min_exp,
               nfeatures = nfeatures_to_use)
    }, warning = function(w) {
      log("WARNING in multiNMF:", w$message)
      invokeRestart("muffleWarning")
    })
  }, error = function(e) {
    log("ERROR in multiNMF:")
    log("  Message:", e$message)
    log("  Call:", deparse(e$call))
    return(NULL)
  })
  
  if (is.null(geneNMF.programs)) {
    log("multiNMF failed - possible causes:")
    log("  1. Too few samples (need at least 2-3)")
    log("  2. Samples have too few cells")
    log("  3. Insufficient variable features")
    log("  4. k range too large for dataset size")
    stop("multiNMF failed - check logs above for details")
  }
  

  saveRDS(geneNMF.programs, file.path(rds_dir, "geneNMF_programs.rds"))
  log("Saved individual NMF programs")
  
  return(geneNMF.programs)
}

#===============================================================================
# 3. EXTRACT META-PROGRAMS
# Aggregate sample-specific programs into robust meta-programs using cosine similarity
# Meta-programs represent recurrent gene expression patterns across multiple samples
#===============================================================================
extract_metaprograms <- function(geneNMF.programs, n_mp, data_dir, rds_dir) {
  log("Extracting meta-programs...")
  log("Parameters:")
  log("  - metric: cosine")
  log("  - weight.explained: 0.8")
  log("  - nMP:", n_mp)
  log("  - min.confidence: 0.7")
  
  geneNMF.metaprograms <- tryCatch({
    getMetaPrograms(geneNMF.programs,
                    metric = "cosine",
                    weight.explained = 0.8,
                    nMP = n_mp,
                    min.confidence = 0.7)
  }, error = function(e) {
    log("Error in getMetaPrograms:", e$message)
    NULL
  })
  
  if (is.null(geneNMF.metaprograms)) {
    stop("getMetaPrograms failed")
  }
  
  # Save metaprograms
  saveRDS(geneNMF.metaprograms, file.path(rds_dir, "geneNMF_metaprograms.rds"))
  log("Saved meta-programs")
  
  # Print metrics
  log("Meta-program metrics:")
  print(geneNMF.metaprograms$metaprograms.metrics)
  write.csv(geneNMF.metaprograms$metaprograms.metrics, 
            file.path(data_dir, "metaprogram_metrics.csv"), 
            row.names = TRUE)
  
  # Save gene lists
  for (mp_name in names(geneNMF.metaprograms$metaprograms.genes)) {
    genes <- geneNMF.metaprograms$metaprograms.genes[[mp_name]]
    weights <- geneNMF.metaprograms$metaprograms.genes.weights[[mp_name]]
    
    df <- data.frame(
      gene = genes,
      weight = weights[genes]
    )
    
    write.csv(df, file.path(data_dir, paste0(mp_name, "_genes.csv")), row.names = FALSE)
  }
  
  log("Saved", length(geneNMF.metaprograms$metaprograms.genes), "meta-program gene lists")
  
  return(geneNMF.metaprograms)
}


#===============================================================================
# 5. CALCULATE MP SCORES USING SEURAT'S ADDMODULESCORE
# Score each cell for meta-program activity using Seurat's module scoring approach
# Scores represent relative expression of program genes vs. control gene sets
#===============================================================================
calculate_mp_scores <- function(seu, geneNMF.metaprograms, rds_dir) {
  log("Calculating meta-program scores with Seurat's AddModuleScore...")
  
  mp.genes <- geneNMF.metaprograms$metaprograms.genes
  
  # Determine number of control genes based on dataset size (min 10, max 50)
  # Control genes are randomly selected for background comparison
  n_genes <- nrow(seu)
  ctrl_genes <- min(50, max(10, floor(n_genes / 50)))
  log("Using", ctrl_genes, "control genes for module scoring")
  
  # Use Seurat's built-in AddModuleScore instead of UCell for better compatibility
  for (mp_name in names(mp.genes)) {
    tryCatch({
      seu <- AddModuleScore(
        object = seu,
        features = list(mp.genes[[mp_name]]),
        name = mp_name,
        ctrl = ctrl_genes,
        seed = 42
      )
      # Rename column to remove the "1" suffix that AddModuleScore adds
      col_name <- paste0(mp_name, "1")
      if (col_name %in% colnames(seu@meta.data)) {
        colnames(seu@meta.data)[colnames(seu@meta.data) == col_name] <- mp_name
      }
    }, error = function(e) {
      log("Warning: Failed to score", mp_name, ":", e$message)
    })
  }
  
  log("Added module scores for", length(mp.genes), "meta-programs")
  
  # Save Seurat object with MP scores
  saveRDS(seu, file.path(rds_dir, "seurat_with_MP_scores.rds"))
  
  list(seu = seu, mp.genes = mp.genes)
}

#===============================================================================
# 5b. PRINT META-PROGRAM INFORMATION
# Display top genes for each meta-program
#===============================================================================
print_mp_info <- function(mp.genes, log) {
  log("=== Meta-Program Gene Information ===")
  for (mp_name in names(mp.genes)) {
    top_genes <- head(mp.genes[[mp_name]], 10)
    log(sprintf("%s: %s", mp_name, paste(top_genes, collapse = ", ")))
  }
}

#===============================================================================
# 5c. RUN GSEA FOR META-PROGRAMS
# Helper function to run GSEA analysis using msigdbr and fgsea
#===============================================================================
runGSEA <- function(program, universe, category = "C8", subcategory = NULL, species = "mouse") {
  if (!requireNamespace("msigdbr", quietly = TRUE) || 
      !requireNamespace("fgsea", quietly = TRUE)) {
    return(NULL)
  }
  
  library(msigdbr)
  library(fgsea)
  
  # Get gene sets from msigdbr
  species_name <- if (tolower(species) == "mouse") "Mus musculus" else "Homo sapiens"
  
  if (!is.null(subcategory)) {
    gene_sets <- msigdbr(species = species_name, category = category, subcategory = subcategory)
  } else {
    gene_sets <- msigdbr(species = species_name, category = category)
  }
  
  # Convert to list format for fgsea
  pathways <- split(gene_sets$gene_symbol, gene_sets$gs_name)
  
  # Create ranked gene list (genes in program get score 1, others get 0)
  # This is a simple enrichment approach
  gene_ranks <- rep(0, length(universe))
  names(gene_ranks) <- universe
  gene_ranks[intersect(program, universe)] <- 1
  
  # Run fgsea
  tryCatch({
    fgsea_res <- fgsea(
      pathways = pathways,
      stats = gene_ranks,
      minSize = 5,
      maxSize = 500,
      nperm = 1000
    )
    
    # Sort by p-value and return
    fgsea_res <- fgsea_res[order(fgsea_res$pval), ]
    return(as.data.frame(fgsea_res))
  }, error = function(e) {
    message("GSEA failed: ", e$message)
    return(NULL)
  })
}

#===============================================================================
# 6. CONDITION COMPARISON
# Compare meta-program activities between experimental conditions
# Identifies condition-specific activation of gene programs
#===============================================================================
compare_conditions <- function(seu, mp.genes, condition_field, data_dir, plots_dir) {
  log("Comparing meta-program activities between conditions...")
  
  conditions <- unique(seu@meta.data[[condition_field]])
  log("Conditions:", paste(conditions, collapse = ", "))
  
  # Statistical comparison between conditions using Wilcoxon test
  # Only performed for pairwise comparisons (2 conditions)
  if (length(conditions) == 2) {
    log("Performing statistical comparison between 2 conditions...")
    
    comp_results <- list()
    
    for (mp_name in names(mp.genes)) {
      scores_cond1 <- seu@meta.data[seu@meta.data[[condition_field]] == conditions[1], mp_name]
      scores_cond2 <- seu@meta.data[seu@meta.data[[condition_field]] == conditions[2], mp_name]
      
      wtest <- wilcox.test(scores_cond1, scores_cond2)
      
      comp_results[[mp_name]] <- data.frame(
        MP = mp_name,
        condition1 = conditions[1],
        condition2 = conditions[2],
        mean_cond1 = mean(scores_cond1, na.rm = TRUE),
        mean_cond2 = mean(scores_cond2, na.rm = TRUE),
        fold_change = mean(scores_cond2, na.rm = TRUE) / mean(scores_cond1, na.rm = TRUE),
        p_value = wtest$p.value
      )
    }
    
    comp_df <- bind_rows(comp_results)
    comp_df$p_adj <- p.adjust(comp_df$p_value, method = "BH")
    comp_df$significant <- comp_df$p_adj < 0.05
    
    write.csv(comp_df, file.path(data_dir, "MP_condition_comparison.csv"), row.names = FALSE)
    
    log("Significant MPs:", sum(comp_df$significant))
    log("MP comparison results saved to:", file.path(data_dir, "MP_condition_comparison.csv"))
  }
}

#===============================================================================
# 7. CREATE MP-BASED UMAP
# Generate UMAP visualizations colored by meta-program activities
# Also creates integrated UMAP using MP scores as features
#===============================================================================
create_mp_umap <- function(seu, mp.genes, condition_field, plots_dir, rds_dir) {
  log("Saving Seurat object with MP scores...")
  
  saveRDS(seu, file.path(rds_dir, "seurat_with_MP_umaps.rds"))
  
  return(seu)
}

#===============================================================================
# 8. PERFORM GSEA (OPTIONAL)
# Functional annotation of meta-programs using Gene Ontology enrichment
# Helps interpret biological meaning of discovered programs
#===============================================================================
perform_gsea <- function(seu, geneNMF.metaprograms, data_dir, species = "mouse") {
  if (!("msigdbr" %in% rownames(installed.packages()) && 
        "fgsea" %in% rownames(installed.packages()))) {
    log("Skipping GSEA (msigdbr or fgsea not installed)")
    return(NULL)
  }
  
  log("Running GSEA for meta-program interpretation...")
  log("Species:", species)
  
  library(msigdbr)
  library(fgsea)
  
  gsea_results <- list()
  
  for (mp_name in names(geneNMF.metaprograms$metaprograms.genes)) {
    log("  GSEA for", mp_name)
    
    tryCatch({
      program <- geneNMF.metaprograms$metaprograms.genes[[mp_name]]
      
      # Run GSEA with C8 (cell type signatures) for cell type identification
      top_pathways_c8 <- runGSEA(program, 
                             universe = rownames(seu), 
                             category = "C8",
                             species = species)
      
      if (!is.null(top_pathways_c8) && nrow(top_pathways_c8) > 0) {
        top_pathways_c8$MP <- mp_name
        top_pathways_c8$category <- "C8_CellType"
        gsea_results[[paste0(mp_name, "_C8")]] <- head(top_pathways_c8, 10)
        
        write.csv(head(top_pathways_c8, 20), 
                 file.path(data_dir, paste0(mp_name, "_GSEA_C8.csv")), 
                 row.names = FALSE)
      }
      
      # Also run with GO Biological Process for functional annotation
      top_pathways_bp <- runGSEA(program, 
                             universe = rownames(seu), 
                             category = "C5", 
                             subcategory = "GO:BP",
                             species = species)
      
      if (!is.null(top_pathways_bp) && nrow(top_pathways_bp) > 0) {
        top_pathways_bp$MP <- mp_name
        top_pathways_bp$category <- "C5_GO_BP"
        gsea_results[[paste0(mp_name, "_GO_BP")]] <- head(top_pathways_bp, 10)
        
        write.csv(head(top_pathways_bp, 20), 
                 file.path(data_dir, paste0(mp_name, "_GSEA_GO_BP.csv")), 
                 row.names = FALSE)
      }
    }, error = function(e) {
      log("    GSEA failed for", mp_name, ":", e$message)
    })
  }
  
  if (length(gsea_results) > 0) {
    all_gsea <- bind_rows(gsea_results)
    write.csv(all_gsea, file.path(data_dir, "all_MP_GSEA_top10.csv"), row.names = FALSE)
    log("Saved GSEA results for", length(unique(all_gsea$MP)), "meta-programs")
    
    # Print summary
    log("=== GSEA Summary (Top Result per MP) ===")
    for (mp_name in unique(all_gsea$MP)) {
      mp_data <- all_gsea[all_gsea$MP == mp_name & all_gsea$category == "C8_CellType", ]
      if (nrow(mp_data) > 0) {
        top_hit <- mp_data[1, ]
        log(sprintf("%s: %s (p=%.3e)", mp_name, top_hit$pathway, top_hit$pval))
      }
    }
  }
  
  return(gsea_results)
}

#===============================================================================
# 8b. COMPARE WITH KNOWN CELL-TYPE SPECIFIC GENE SIGNATURES
# Check if discovered programs match known cell-type specific gene signatures
# from literature (NK cells, T cells, B cells, Macrophages, etc.)
#===============================================================================
compare_with_known_signatures <- function(geneNMF.metaprograms, seu, data_dir, species = "mouse") {
  log("Comparing meta-programs with known cell-type signatures...")
  
  # Define known gene signatures from literature
  known_signatures <- list()
  
  if (tolower(species) == "mouse") {
    known_signatures <- list(
      # NK cell signatures
      "NK_Cytotoxic" = c("Gzma", "Gzmb", "Prf1", "Nkg7", "Klrk1", "Klrd1", "Ncr1", "Klrb1c"),
      "NK_Maturation" = c("Cd27", "Cd11b", "Itgam", "Cd49b", "Klra4", "Klra8", "Klra9"),
      "NK_Activation" = c("Ifng", "Tnf", "Ccl3", "Ccl4", "Xcl1", "Cd69", "Fasl"),
      "NK_Cytokine_Response" = c("Stat1", "Stat4", "Irf1", "Irf8", "Il15ra", "Il18r1", "Il12rb2"),
      "NK_Memory" = c("Ly6c2", "Cxcr6", "Cd49a", "Itga1", "Cxcr3"),
      
      # T cell signatures
      "T_Cytotoxic" = c("Cd8a", "Cd8b1", "Gzmk", "Cxcr3", "Ccl5", "Prf1"),
      "T_Exhaustion" = c("Pdcd1", "Havcr2", "Lag3", "Tigit", "Ctla4", "Tox"),
      "T_Proliferation" = c("Mki67", "Top2a", "Stmn1", "Pcna", "Cdk1"),
      "Treg" = c("Foxp3", "Il2ra", "Ctla4", "Ikzf2", "Nt5e"),
      
      # B cell signatures
      "B_Activated" = c("Cd86", "Cd80", "H2-Ab1", "H2-Aa", "Cd83"),
      "B_Plasma" = c("Sdc1", "Xbp1", "Prdm1", "Jchain", "Mzb1"),
      
      # Myeloid signatures
      "Macrophage_M1" = c("Nos2", "Il1b", "Tnf", "Il6", "Cxcl10", "Cd86"),
      "Macrophage_M2" = c("Arg1", "Mrc1", "Cd163", "Il10", "Tgfb1", "Ccl17"),
      "Monocyte_Classical" = c("Ly6c2", "Ccr2", "Cd14", "S100a8", "S100a9"),
      "Monocyte_NonClassical" = c("Cx3cr1", "Nr4a1", "Ace", "Fcgr4"),
      "DC_Conventional" = c("Itgax", "Flt3", "Xcr1", "Clec9a", "Batf3"),
      
      # Granulocyte signatures
      "Neutrophil" = c("S100a8", "S100a9", "Cxcr2", "Csf3r", "Ly6g"),
      
      # Stromal signatures
      "Fibroblast" = c("Col1a1", "Col1a2", "Dcn", "Lum", "Pdgfra"),
      "Endothelial" = c("Pecam1", "Cdh5", "Vwf", "Kdr", "Eng")
    )
  } else {  # human
    known_signatures <- list(
      # NK cell signatures (human)
      "NK_Cytotoxic" = c("GZMA", "GZMB", "PRF1", "NKG7", "KLRK1", "KLRD1", "NCR1", "KLRB1"),
      "NK_Maturation" = c("CD27", "CD11B", "ITGAM", "CD49B", "KLRC1", "KLRC2"),
      "NK_Activation" = c("IFNG", "TNF", "CCL3", "CCL4", "XCL1", "CD69", "FASLG"),
      "NK_Memory" = c("LY6C", "CXCR6", "CD49A", "ITGA1", "CXCR3"),
      
      # T cell signatures (human)
      "T_Cytotoxic" = c("CD8A", "CD8B", "GZMK", "CXCR3", "CCL5", "PRF1"),
      "T_Exhaustion" = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX"),
      "T_Proliferation" = c("MKI67", "TOP2A", "STMN1", "PCNA", "CDK1"),
      "Treg" = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "NT5E"),
      
      # B cell signatures (human)
      "B_Activated" = c("CD86", "CD80", "HLA-DRA", "HLA-DRB1", "CD83"),
      "B_Plasma" = c("SDC1", "XBP1", "PRDM1", "JCHAIN", "MZB1"),
      
      # Myeloid signatures (human)
      "Macrophage_M1" = c("NOS2", "IL1B", "TNF", "IL6", "CXCL10", "CD86"),
      "Macrophage_M2" = c("ARG1", "MRC1", "CD163", "IL10", "TGFB1", "CCL17"),
      "Monocyte_Classical" = c("CD14", "CCR2", "S100A8", "S100A9", "FCGR1A"),
      "Monocyte_NonClassical" = c("FCGR3A", "CX3CR1", "NR4A1", "ACE"),
      "DC_Conventional" = c("ITGAX", "FLT3", "XCR1", "CLEC9A", "BATF3")
    )
  }
  
  log("Loaded", length(known_signatures), "known cell-type gene signatures")
  
  # Calculate overlap between meta-programs and known signatures
  comparison_results <- list()
  
  for (mp_name in names(geneNMF.metaprograms$metaprograms.genes)) {
    mp_genes <- geneNMF.metaprograms$metaprograms.genes[[mp_name]]
    
    for (sig_name in names(known_signatures)) {
      sig_genes <- known_signatures[[sig_name]]
      
      # Calculate overlap
      overlap <- intersect(mp_genes, sig_genes)
      n_overlap <- length(overlap)
      n_mp <- length(mp_genes)
      n_sig <- length(sig_genes)
      
      # Calculate Fisher's exact test p-value
      # Contingency table: MP vs signature overlap
      universe_size <- nrow(seu)
      
      # Fisher test
      cont_table <- matrix(c(
        n_overlap,                          # in both
        n_sig - n_overlap,                  # in signature only
        n_mp - n_overlap,                   # in MP only
        universe_size - n_mp - n_sig + n_overlap  # in neither
      ), nrow = 2)
      
      fisher_p <- tryCatch({
        fisher.test(cont_table, alternative = "greater")$p.value
      }, error = function(e) { 1 })
      
      # Jaccard index
      jaccard <- n_overlap / (n_mp + n_sig - n_overlap)
      
      comparison_results[[paste0(mp_name, "_", sig_name)]] <- data.frame(
        MetaProgram = mp_name,
        KnownSignature = sig_name,
        n_MP_genes = n_mp,
        n_Signature_genes = n_sig,
        n_Overlap = n_overlap,
        Overlap_Genes = paste(overlap, collapse = ", "),
        Jaccard_Index = jaccard,
        Fisher_pvalue = fisher_p,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine results
  comparison_df <- bind_rows(comparison_results)
  comparison_df$Fisher_padj <- p.adjust(comparison_df$Fisher_pvalue, method = "BH")
  comparison_df <- comparison_df %>% arrange(Fisher_pvalue)
  
  # Save results
  write.csv(comparison_df, file.path(data_dir, "MP_vs_KnownSignatures.csv"), row.names = FALSE)
  
  # Print significant matches
  log("=== Significant Matches with Known Signatures (FDR < 0.05) ===")
  sig_matches <- comparison_df %>% filter(Fisher_padj < 0.05, n_Overlap >= 3)
  
  if (nrow(sig_matches) > 0) {
    for (i in 1:min(20, nrow(sig_matches))) {
      row <- sig_matches[i, ]
      log(sprintf("%s <-> %s: %d/%d genes overlap (Jaccard=%.3f, p=%.2e)",
                 row$MetaProgram, row$KnownSignature, 
                 row$n_Overlap, row$n_Signature_genes,
                 row$Jaccard_Index, row$Fisher_pvalue))
    }
  } else {
    log("No significant matches found (FDR < 0.05)")
  }
  
  # Create heatmap of Jaccard similarities
  plot_mp_signatures_jaccard(comparison_df, file.path(dirname(data_dir), "plots"))
  
  return(comparison_df)
}

################################################################################
# MAIN ANALYSIS
################################################################################

run_gene_program_discovery <- function(seurat_path, output_dir, 
                                       min_k = 4, max_k = 9, n_mp = 10,
                                       min_exp = 0.05, condition_field = NULL,
                                       species = "mouse",
                                       data_dir = NULL, plots_dir = NULL, rds_dir = NULL,
                                       annotations_path = NULL) {
  
  log("=== Starting GeneNMF Gene Program Discovery ===")
  log("Input Seurat:", seurat_path)
  if (!is.null(annotations_path)) {
    log("Input Annotations:", annotations_path)
  }
  log("Output Dir:", output_dir)
  log("k range:", min_k, "to", max_k)
  log("Target MPs:", n_mp)
  log("Species:", species)
  
  # Use provided directories or create defaults
  if (is.null(data_dir)) data_dir <- file.path(output_dir, "data")
  if (is.null(plots_dir)) plots_dir <- file.path(output_dir, "plots")
  if (is.null(rds_dir)) rds_dir <- file.path(output_dir, "rds_objects")
  
  # 1. Prepare Seurat object and load upstream results
  prep <- prepare_seurat_object(seurat_path, condition_field, annotations_path)
  seu <- prep$seu
  condition_field <- prep$condition_field
  sample_field <- prep$sample_field
  de_genes <- prep$de_genes
  lr_genes <- prep$lr_genes
  
  # Combine upstream gene sets with HVGs for feature selection
  # This integrates results from differential analysis and cell communication
  feature_genes <- unique(c(de_genes, lr_genes))
  if (length(feature_genes) > 0) {
    log("Using", length(feature_genes), "genes from upstream analyses for NMF input")
    # Subset to genes present in Seurat
    feature_genes <- intersect(feature_genes, rownames(seu))
    log("  -", length(feature_genes), "genes available in dataset")
    
    write.csv(data.frame(gene = feature_genes, 
                        source = ifelse(feature_genes %in% de_genes, "DE", "Ligand")),
             file.path(data_dir, "feature_genes_used.csv"), 
             row.names = FALSE)
  } else {
    log("No upstream genes available, will use variable features from Seurat")
  }
  
  # 2. Run multiNMF
  geneNMF.programs <- run_multinmf(seu, sample_field, min_k, max_k, min_exp, rds_dir)
  
  # 3. Extract meta-programs
  geneNMF.metaprograms <- extract_metaprograms(geneNMF.programs, n_mp, data_dir, rds_dir)
  
  # 4. Visualize meta-programs
  visualize_metaprograms(geneNMF.metaprograms, plots_dir)
  
  # 5. Calculate MP scores
  mp_scores <- calculate_mp_scores(seu, geneNMF.metaprograms, rds_dir)
  seu <- mp_scores$seu
  mp.genes <- mp_scores$mp.genes
  
  # 5b. Print MP information
  print_mp_info(mp.genes, log)
  
  # 5c. Create combined metrics and genes table visualization
  create_mp_combined_metrics_genes_table(geneNMF.metaprograms, mp.genes, plots_dir, n_genes = 10, log)
  
  # 6. Compare conditions - create line plots
  compare_conditions(seu, mp.genes, condition_field, data_dir, plots_dir)
  
  # 6b. Create MP line plots by condition (all in one figure)
  create_mp_lines_by_condition(seu, mp.genes, condition_field, plots_dir, log)
  
  # 7. Create MP-based UMAP visualizations
  create_mp_umap_plots(seu, mp.genes, plots_dir, log)
  seu <- create_mp_umap(seu, mp.genes, condition_field, plots_dir, rds_dir)
  
  # 8. Perform GSEA
  gsea_results <- perform_gsea(seu, geneNMF.metaprograms, data_dir, species = species)
  
  # 8a. Compare with known cell-type specific signatures
  known_sig_comparison <- compare_with_known_signatures(geneNMF.metaprograms, seu, data_dir, species = species)
  
  # 8a2. Create table visualizations for known signatures
  if (!is.null(known_sig_comparison) && nrow(known_sig_comparison) > 0) {
    create_known_signatures_table(known_sig_comparison, plots_dir, n_top = 20, log)
    create_signature_overlap_genes_table(known_sig_comparison, plots_dir, n_top = 10, log)
  }
  
  # 8b. Create GSEA visualizations
  if (!is.null(gsea_results)) {
    create_gsea_dotplot(gsea_results, plots_dir, n_top = 5, log)
  }

  # 8c. MP scores: cell type comparison between two conditions (A vs B)
  create_mp_celltype_condition_comparison(
    seu, mp.genes, condition_field,
    plots_dir, data_dir,
    gsea_results,
    log
  )
  
  # 9. Save summary and final outputs
  # Generate comprehensive summary of analysis parameters and results
  summary_info <- list(
    timestamp = Sys.time(),
    n_cells = ncol(seu),
    n_genes = nrow(seu),
    n_samples = length(unique(seu@meta.data[[sample_field]])),
    n_conditions = length(unique(seu@meta.data[[condition_field]])),
    k_range = paste(min_k, "to", max_k),
    n_metaprograms = length(mp.genes),
    mp_names = names(mp.genes),
    n_de_genes_used = length(de_genes),
    n_lr_genes_used = length(lr_genes)
  )
  
  saveRDS(summary_info, file.path(data_dir, "analysis_summary.rds"))
  
  saveRDS(seu, file.path(data_dir, "seurat_with_gene_programs.rds"))
  log("Saved final Seurat object with gene programs")
  
  saveRDS(geneNMF.programs, file.path(data_dir, "nmf_result.rds"))
  log("Saved NMF results")
  
  log("=== Gene Program Discovery Completed ===")
  log("Summary:")
  log("  - Samples:", summary_info$n_samples)
  log("  - Conditions:", summary_info$n_conditions)
  log("  - Meta-programs:", summary_info$n_metaprograms)
  log("  - DE genes integrated:", summary_info$n_de_genes_used)
  log("  - L-R genes integrated:", summary_info$n_lr_genes_used)
  log("  - Output directory:", output_dir)
  
  return(list(
    seu = seu,
    metaprograms = geneNMF.metaprograms,
    mp.genes = mp.genes,
    summary = summary_info
  ))
}

################################################################################
# EXECUTE
################################################################################

log("=== Starting GeneNMF Gene Program Discovery ===")

results <- run_gene_program_discovery(
  seurat_path = seurat_path,
  annotations_path = if (use_combined) annotations_path else NULL,
  output_dir = base_dir,
  min_k = opt$min_rank,
  max_k = opt$max_rank,
  n_mp = opt$n_mp,
  min_exp = opt$min_exp,
  condition_field = opt$condition_field,
  species = opt$species,
  data_dir = data_dir,
  plots_dir = plots_dir,
  rds_dir = rds_dir
)

cat("Gene program discovery completed successfully!\n")

# End performance tracking
end_performance_tracking(
  perf_tracker,
  success = TRUE,
  additional_metrics = list(
    n_cells = ncol(results$seu),
    n_samples = length(unique(results$seu$orig.ident))
  )
)

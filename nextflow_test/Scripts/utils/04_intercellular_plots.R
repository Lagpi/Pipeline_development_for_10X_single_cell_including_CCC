################################################################################
# PLOTTING FUNCTIONS FOR INTERCELLULAR COMMUNICATION ANALYSIS (MODULE 04)
#
# This file contains all plotting functions for:
# - CellChat comparison analysis (multi-condition)
# - CellChat single analysis
# - NicheNet ligand-receptor prediction
#
# All functions follow consistent API:
# - Save to disk automatically
# - Return plot object for flexibility
# - Use CFG_* parameters for styling
################################################################################

# Configuration parameters (read from parent environment or use defaults)
CFG_RANDOM_SEED <- 42
CFG_BASE_SIZE <- if (exists("CFG_BASE_SIZE", envir = parent.frame(), inherits = TRUE)) {
  get("CFG_BASE_SIZE", envir = parent.frame(), inherits = TRUE)
} else { 18 }
CFG_TITLE_SIZE <- if (exists("CFG_TITLE_SIZE", envir = parent.frame(), inherits = TRUE)) {
  get("CFG_TITLE_SIZE", envir = parent.frame(), inherits = TRUE)
} else { 24 }
CFG_AXIS_TITLE_SIZE <- if (exists("CFG_AXIS_TITLE_SIZE", envir = parent.frame(), inherits = TRUE)) {
  get("CFG_AXIS_TITLE_SIZE", envir = parent.frame(), inherits = TRUE)
} else { 18 }
CFG_LEGEND_SIZE <- if (exists("CFG_LEGEND_SIZE", envir = parent.frame(), inherits = TRUE)) {
  get("CFG_LEGEND_SIZE", envir = parent.frame(), inherits = TRUE)
} else { 18 }
CFG_FDR_THRESHOLD <- 0.05
CFG_LOG2FC_THRESHOLD <- 0.25
CFG_MIN_CELLS_PER_TYPE <- 10
CFG_PLOT_WIDTH <- 12
CFG_PLOT_HEIGHT <- 8
CFG_PLOT_DPI <- 150

################################################################################
# CELLCHAT COMPARISON PLOTS (MULTI-CONDITION)
################################################################################

#' Plot 1: CellChat Total Interactions Comparison
#' 
#' Creates combined plot showing number and strength of interactions
#' across conditions
#' 
#' @param cellchat_merged Merged CellChat object
#' @param de_source String describing DE gene source
#' @param de_count Number of DE genes used
#' @param plots_dir Output directory
#' @return Combined patchwork plot object
create_cellchat_interaction_comparison <- function(cellchat_merged, de_source, de_count, plots_dir) {
  # Ensure condition names are used in metadata
  condition_names <- names(cellchat_merged@net)
  if ("datasets" %in% colnames(cellchat_merged@meta)) {
    cellchat_merged@meta$datasets <- factor(cellchat_merged@meta$datasets, 
                                             levels = condition_names)
  }
  
  p1 <- compareInteractions(cellchat_merged, show.legend = FALSE)
  p2 <- compareInteractions(cellchat_merged, show.legend = FALSE, measure = "weight")
  
  # Standardize axis labels and titles
  p1 <- p1 + 
    labs(title = "Number of Interactions", x = "Condition", y = "Number") +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14)
    )
  
  p2 <- p2 + 
    labs(title = "Interaction Strength", x = "Condition", y = "Strength") +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14)
    )
  
  de_note <- paste0("DE genes source: ", de_source, " (n=", de_count, ")")
  p_combined <- p1 + p2 + 
    patchwork::plot_annotation(
      title = "CellChat: Total Interactions Comparison",
      caption = de_note,
      theme = theme(
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        plot.caption = element_text(hjust = 0, size = 12, color = "gray30")
      )
    )
  
  save_dual_plot(p_combined, "cellchat_total_interactions_comparison", 
                plots_dir = plots_dir, width = 12, height = 6)
  
  return(p_combined)
}

#' Plot 2: CellChat Differential Interaction Circle Plots
#' 
#' Creates circle plots showing differential interactions between conditions
#' 
#' @param cellchat_merged Merged CellChat object
#' @param plots_dir Output directory
#' @param log Logging function
#' @return TRUE if successful, FALSE otherwise
create_cellchat_diff_circle <- function(cellchat_merged, plots_dir, log = cat) {
  # Check if there are differential edges
  net_diff <- cellchat_merged@net[[2]]$count - cellchat_merged@net[[1]]$count
  n_diff_edges <- sum(net_diff != 0)
  log(sprintf("  Number of differential edges: %d", n_diff_edges))
  
  if (n_diff_edges == 0) {
    log("  Skipping differential circle plots - no differential edges found")
    return(FALSE)
  }
  
  # Get condition names
  condition_names <- names(cellchat_merged@net)
  
  tryCatch({
    png(file.path(plots_dir, "cellchat_diff_circle.png"), 
        width = 14, height = 7, units = "in", res = 150)
    par(mfrow = c(1,2), xpd = TRUE)
    netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, 
                              comparison = c(1,2), show.legend = TRUE)
    netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight", 
                              comparison = c(1,2), show.legend = TRUE)
    dev.off()
    log("  Saved differential interaction circle plots")
    log(sprintf("  Comparing: %s vs %s", condition_names[2], condition_names[1]))
    return(TRUE)
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    log(sprintf("  Warning: Could not create differential circle plots: %s", e$message))
    log("  This is a known CellChat issue with certain data configurations")
    return(FALSE)
  })
}

#' Plot 3: CellChat Differential Interaction Heatmaps
#' 
#' Creates heatmaps showing differential interactions (count and weight)
#' 
#' @param cellchat_merged Merged CellChat object
#' @param plots_dir Output directory
#' @param log Logging function
#' @return TRUE if successful, FALSE otherwise
create_cellchat_diff_heatmap <- function(cellchat_merged, plots_dir, log = cat) {
  # Get condition names for informative logging
  condition_names <- names(cellchat_merged@net)
  
  tryCatch({
    p1 <- netVisual_heatmap(cellchat_merged, comparison = c(1,2))
    p2 <- netVisual_heatmap(cellchat_merged, measure = "weight", comparison = c(1,2))
    
    png(file.path(plots_dir, "cellchat_diff_heatmap.png"), 
        width = 12, height = 6, units = "in", res = 150)
    print(p1 + p2)
    dev.off()
    log("  Saved differential interaction heatmaps")
    log(sprintf("  Comparing: %s vs %s", condition_names[2], condition_names[1]))
    return(TRUE)
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    log(sprintf("  Warning: Could not create differential heatmaps: %s", e$message))
    return(FALSE)
  })
}

#' Plot 3b: CellChat Differential Interaction Heatmaps (All Pairs)
#'
#' Creates differential heatmaps for every pair of conditions in a combined figure
#'
#' @param cellchat_list List of CellChat objects per condition
#' @param plots_dir Output directory
#' @param log Logging function
#' @return TRUE if any plots were generated, FALSE otherwise
create_cellchat_diff_heatmaps_all <- function(cellchat_list, plots_dir, log = cat) {
  if (length(cellchat_list) < 2) return(FALSE)

  cond_names <- names(cellchat_list)
  if (is.null(cond_names) || any(cond_names == "")) {
    cond_names <- paste0("Condition", seq_along(cellchat_list))
    names(cellchat_list) <- cond_names
  }

  sanitize_name <- function(x) gsub("[^A-Za-z0-9_-]+", "_", x)
  
  log("  Creating differential heatmaps for all pairwise comparisons...")
  
  # First pass: calculate global color scale limits
  pair_idx <- combn(seq_along(cellchat_list), 2, simplify = FALSE)
  all_count_diffs <- c()
  all_weight_diffs <- c()
  
  for (idx in pair_idx) {
    pair_list <- cellchat_list[idx]
    merged_pair <- tryCatch({
      mergeCellChat(pair_list, add.names = names(pair_list))
    }, error = function(e) NULL)
    
    if (!is.null(merged_pair)) {
      count_mat1 <- merged_pair@net[[1]]$count
      count_mat2 <- merged_pair@net[[2]]$count
      all_count_diffs <- c(all_count_diffs, as.vector(count_mat2 - count_mat1))
      
      weight_mat1 <- merged_pair@net[[1]]$weight
      weight_mat2 <- merged_pair@net[[2]]$weight
      all_weight_diffs <- c(all_weight_diffs, as.vector(weight_mat2 - weight_mat1))
    }
  }
  
  # Calculate symmetric limits
  count_limit <- max(abs(all_count_diffs), na.rm = TRUE)
  weight_limit <- max(abs(all_weight_diffs), na.rm = TRUE)
  
  log(sprintf("  Global scale limits - Count: [%.1f, %.1f], Weight: [%.2f, %.2f]", 
              -count_limit, count_limit, -weight_limit, weight_limit))
  
  # Second pass: create all plots with consistent scales
  all_plots <- list()
  
  for (idx in pair_idx) {
    pair_list <- cellchat_list[idx]
    pair_names <- names(pair_list)
    comparison_title <- sprintf("%s vs %s", pair_names[1], pair_names[2])

    tryCatch({
      merged_pair <- mergeCellChat(pair_list, add.names = pair_names)

      # Create heatmaps with fixed color scales
      p1 <- netVisual_heatmap(merged_pair, color.use = c("blue", "white", "red"),
                              color.heatmap = c(-count_limit, count_limit))
      p2 <- netVisual_heatmap(merged_pair, measure = "weight", 
                              color.use = c("blue", "white", "red"),
                              color.heatmap = c(-weight_limit, weight_limit))
      
      combined_pair <- (p1 | p2) + 
        patchwork::plot_annotation(
          title = comparison_title,
          theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
        )
      
      all_plots[[length(all_plots) + 1]] <- combined_pair

      # Save individual comparison
      file_name <- paste0(
        "cellchat_diff_heatmap_",
        sanitize_name(pair_names[1]),
        "_vs_",
        sanitize_name(pair_names[2]),
        ".png"
      )

      png(file.path(plots_dir, file_name),
          width = 12, height = 6, units = "in", res = 150)
      print(combined_pair)
      dev.off()
      
      log(sprintf("  - Saved %s", comparison_title))

    }, error = function(e) {
      if (dev.cur() != 1) dev.off()
      log(sprintf("  Warning: Could not create heatmap for %s: %s",
                  comparison_title, e$message))
    })
  }

  # Create combined figure with all comparisons (2 per row)
  if (length(all_plots) > 0) {
    tryCatch({
      n_plots <- length(all_plots)
      n_cols <- 2  # Always 2 columns as requested
      n_rows <- ceiling(n_plots / n_cols)
      
      combined_all <- patchwork::wrap_plots(all_plots, ncol = n_cols, nrow = n_rows) +
        patchwork::plot_annotation(
          title = "CellChat: Differential Interactions Between All Conditions",
          theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
        )
      
      # Calculate appropriate figure size
      fig_width <- 24  # 12 inches per column × 2 columns
      fig_height <- 7 * n_rows  # 7 inches per row
      
      png(file.path(plots_dir, "cellchat_diff_heatmap_all_comparisons.png"),
          width = fig_width, height = fig_height, units = "in", res = 150)
      print(combined_all)
      dev.off()
      
      log(sprintf("  Saved combined figure with %d comparisons (%d rows × 2 columns)", 
                  n_plots, n_rows))
      return(TRUE)
    }, error = function(e) {
      if (dev.cur() != 1) dev.off()
      log(sprintf("  Error creating combined heatmap: %s", e$message))
      return(length(all_plots) > 0)
    })
  } else {
    log("  No plots were generated")
  }

  return(length(all_plots) > 0)
}

  #' Plot 3c: CellChat Network Circle Plots (Per Condition)
  #'
  #' Creates circle plots showing interaction count and weight for each condition
  #'
  #' @param cellchat_list List of CellChat objects per condition
  #' @param plots_dir Output directory
  #' @param log Logging function
  #' @return TRUE if any plots were generated, FALSE otherwise
  create_cellchat_network_circles_by_condition <- function(cellchat_list, plots_dir, log = cat) {
    if (length(cellchat_list) < 1) return(FALSE)

    cond_names <- names(cellchat_list)
    if (is.null(cond_names) || any(cond_names == "")) {
      cond_names <- paste0("Condition", seq_along(cellchat_list))
      names(cellchat_list) <- cond_names
    }

    sanitize_name <- function(x) gsub("[^A-Za-z0-9_-]+", "_", x)
    plotted_any <- FALSE

    for (cond in names(cellchat_list)) {
      cellchat <- cellchat_list[[cond]]
      safe_cond <- sanitize_name(cond)

      tryCatch({
        pdf(file.path(plots_dir, paste0("cellchat_network_circles_", safe_cond, ".pdf")),
            width = 20, height = 10)
        par(mfrow = c(1, 2), xpd = TRUE)
        netVisual_circle(cellchat@net$count,
                         vertex.weight = as.numeric(table(cellchat@idents)),
                         weight.scale = TRUE, label.edge = FALSE,
                         title.name = paste0(cond, ": Number of interactions"))
        netVisual_circle(cellchat@net$weight,
                         vertex.weight = as.numeric(table(cellchat@idents)),
                         weight.scale = TRUE, label.edge = FALSE,
                         title.name = paste0(cond, ": Interaction weights/strength"))
        dev.off()

        log("  Saved circle plots for condition:", cond)
        plotted_any <- TRUE
      }, error = function(e) {
        if (dev.cur() != 1) dev.off()
        log(sprintf("  Warning: Could not create circle plots for %s: %s", cond, e$message))
      })
    }

    return(plotted_any)
  }

#' Plot 4: CellChat Signaling Role Scatter Plots
#' 
#' Creates scatter plots showing signaling roles for each condition
#' 
#' @param cellchat_list List of CellChat objects per condition
#' @param de_source String describing DE gene source
#' @param de_count Number of DE genes used
#' @param plots_dir Output directory
#' @return Combined patchwork plot object
create_cellchat_signaling_roles <- function(cellchat_list, de_source, de_count, plots_dir) {
  # Calculate weight range for consistent scaling
  num_link <- sapply(cellchat_list, function(x) {
    rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
  })
  weight_minmax <- c(min(num_link), max(num_link))
  
  # Create scatter plot for each condition
  gg_list <- list()
  for (i in 1:length(cellchat_list)) {
    gg_list[[i]] <- netAnalysis_signalingRole_scatter(
      cellchat_list[[i]], 
      title = names(cellchat_list)[i],
      weight.MinMax = weight_minmax
    )
  }
  
  de_note <- paste0("DE genes source: ", de_source, " (n=", de_count, ")")
  p_combined <- patchwork::wrap_plots(plots = gg_list) +
    patchwork::plot_annotation(
      caption = de_note,
      theme = theme(plot.caption = element_text(hjust = 0, size = 9, color = "gray30"))
    )
  
  save_dual_plot(p_combined, "cellchat_signaling_roles_comparison", 
                plots_dir = plots_dir, width = 5 * length(cellchat_list), height = 5.3)
  
  return(p_combined)
}

#' Plot 5: CellChat Functional Similarity Analysis
#' 
#' Creates embedding and pathway distance plots for functional similarity
#' 
#' @param cellchat_merged Merged CellChat object
#' @param plots_dir Output directory
#' @param log Logging function
#' @return TRUE if successful, FALSE otherwise
create_cellchat_functional_similarity <- function(cellchat_merged, plots_dir, log = cat) {
  # Get condition names
  condition_names <- names(cellchat_merged@net)
  
  # Ensure condition names are used in metadata
  if ("datasets" %in% colnames(cellchat_merged@meta)) {
    cellchat_merged@meta$datasets <- factor(cellchat_merged@meta$datasets, 
                                             levels = condition_names)
  }
  
  tryCatch({
    cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "functional")
    
    # Use R-based UMAP if available
    if (requireNamespace("umap", quietly = TRUE)) {
      log("  Using R umap package for embedding")
      cellchat_merged <- netEmbedding(cellchat_merged, type = "functional", umap.method = "umap-learn")
    } else {
      log("  Attempting Python UMAP for embedding")
      cellchat_merged <- netEmbedding(cellchat_merged, type = "functional")
    }
    
    cellchat_merged <- netClustering(cellchat_merged, type = "functional")
    
    # Embedding plot with condition names
    p_embed <- netVisual_embeddingPairwise(cellchat_merged, type = "functional", 
                                            label.size = 3.5, group.names = condition_names)
    save_dual_plot(p_embed, "cellchat_functional_similarity", plots_dir = plots_dir, width = 10, height = 8)
    
    # Pathway distance with condition names
    p_dist <- rankSimilarity(cellchat_merged, type = "functional", group.names = condition_names)
    save_dual_plot(p_dist, "cellchat_pathway_distance", plots_dir = plots_dir, width = 5, height = 5)
    
    log("  Saved functional similarity plots")
    log(sprintf("  Comparing: %s", paste(condition_names, collapse = " vs ")))
    return(TRUE)
  }, error = function(e) {
    log(sprintf("  Warning: Could not create functional similarity plots: %s", e$message))
    if (grepl("UMAP", e$message, ignore.case = TRUE)) {
      log("  Install R umap package: install.packages('umap')")
      log("  Or Python umap-learn: pip install umap-learn")
    }
    return(FALSE)
  })
}

#' Plot 6: CellChat Pathway Ranking
#' 
#' Creates ranked pathway comparison plots (stacked and unstacked)
#' 
#' @param cellchat_merged Merged CellChat object
#' @param de_source String describing DE gene source
#' @param de_count Number of DE genes used
#' @param plots_dir Output directory
#' @param log Logging function
#' @return TRUE if successful, FALSE otherwise
create_cellchat_pathway_ranking <- function(cellchat_merged, de_source, de_count, plots_dir, log = cat) {
  # Ensure condition names are used in metadata
  condition_names <- names(cellchat_merged@net)
  if ("datasets" %in% colnames(cellchat_merged@meta)) {
    cellchat_merged@meta$datasets <- factor(cellchat_merged@meta$datasets, 
                                             levels = condition_names)
  }
  
  tryCatch({
    p1 <- rankNet(cellchat_merged, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE)
    p2 <- rankNet(cellchat_merged, mode = "comparison", measure = "weight", stacked = FALSE, do.stat = TRUE)
    
    de_note <- paste0("DE genes source: ", de_source, " (n=", de_count, ")")
    p_combined <- p1 + p2 +
      patchwork::plot_annotation(
        caption = de_note,
        theme = theme(plot.caption = element_text(hjust = 0, size = 9, color = "gray30"))
      )
    
    save_dual_plot(p_combined, "cellchat_pathway_ranking", plots_dir = plots_dir, width = 12, height = 6.3)
    log("  Saved pathway ranking plots")
    return(TRUE)
  }, error = function(e) {
    log(sprintf("  Warning: Could not create pathway ranking plots: %s", e$message))
    return(FALSE)
  })
}

#' Plot 7: CellChat Outgoing Signaling Heatmaps
#' 
#' Creates heatmaps showing outgoing signaling patterns for both conditions
#' 
#' @param cellchat_list List of CellChat objects per condition
#' @param plots_dir Output directory
#' @param log Logging function
#' @return TRUE if successful, FALSE otherwise
create_cellchat_outgoing_heatmap <- function(cellchat_list, plots_dir, log = cat) {
  tryCatch({
    pathway_union <- union(cellchat_list[[1]]@netP$pathways, cellchat_list[[2]]@netP$pathways)
    
    ht1 <- netAnalysis_signalingRole_heatmap(
      cellchat_list[[1]], 
      pattern = "outgoing",
      signaling = pathway_union,
      title = names(cellchat_list)[1],
      width = 5, 
      height = 8
    )
    
    ht2 <- netAnalysis_signalingRole_heatmap(
      cellchat_list[[2]],
      pattern = "outgoing",
      signaling = pathway_union,
      title = names(cellchat_list)[2],
      width = 5,
      height = 8
    )
    
    png(file.path(plots_dir, "cellchat_outgoing_heatmap.png"), 
        width = 12, height = 10, units = "in", res = 150)
    ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
    dev.off()
    
    log("  Saved outgoing signaling heatmaps")
    return(TRUE)
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    log(sprintf("  Warning: Could not create outgoing signaling heatmaps: %s", e$message))
    return(FALSE)
  })
}

################################################################################
# CELLCHAT SINGLE ANALYSIS PLOTS
################################################################################

#' Plot 8: CellChat Network Circle Plots (Single Analysis)
#' 
#' Creates circle plots showing interaction count and weight
#' 
#' @param cellchat CellChat object
#' @param plots_dir Output directory
create_cellchat_network_circles <- function(cellchat, plots_dir) {
  pdf(file.path(plots_dir, "cellchat_network_circles.pdf"), width = 20, height = 10)
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)), 
                   weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = as.numeric(table(cellchat@idents)), 
                   weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
  dev.off()
}

#' Plot 9: CellChat Bubble Plot (Top Interactions)
#' 
#' Creates bubble plot showing significant interactions
#' 
#' @param cellchat CellChat object
#' @param plots_dir Output directory
create_cellchat_bubble_plot <- function(cellchat, plots_dir) {
  pdf(file.path(plots_dir, "cellchat_bubble_top_interactions.pdf"), width = 12, height = 12)
  print(netVisual_bubble(cellchat, remove.isolate = FALSE, 
                        title.name = "Significant interactions (Overview)"))
  dev.off()
}

#' Plot 10: CellChat Heatmap (Count)
#' 
#' Creates heatmap showing interaction counts
#' 
#' @param cellchat CellChat object
#' @param plots_dir Output directory
create_cellchat_heatmap_count <- function(cellchat, plots_dir) {
  pdf(file.path(plots_dir, "cellchat_heatmap_count.pdf"), width = 10, height = 8)
  print(netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues"))
  dev.off()
}

#' Plot 11: CellChat Heatmap (Weight)
#' 
#' Creates heatmap showing interaction weights
#' 
#' @param cellchat CellChat object
#' @param plots_dir Output directory
create_cellchat_heatmap_weight <- function(cellchat, plots_dir) {
  pdf(file.path(plots_dir, "cellchat_heatmap_weight.pdf"), width = 10, height = 8)
  print(netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds"))
  dev.off()
}

#' Plot 12: CellChat Centrality Heatmaps
#' 
#' Creates network centrality heatmaps for top pathways
#' 
#' @param cellchat CellChat object with computed centrality
#' @param plots_dir Output directory
#' @param n_pathways Number of top pathways to plot
#' @param log Logging function
create_cellchat_centrality_heatmap <- function(cellchat, plots_dir, n_pathways = 10, log = cat) {
  pathways.show <- head(cellchat@netP$pathways, min(n_pathways, length(cellchat@netP$pathways)))
  
  if (length(pathways.show) == 0) {
    log("  No pathways available for centrality heatmap")
    return(FALSE)
  }
  
  pdf(file.path(plots_dir, "cellchat_centrality_heatmap.pdf"), width = 12, height = 8)
  par(mfrow = c(1, 1))
  for (pathway in pathways.show) {
    tryCatch({
      netAnalysis_signalingRole_network(cellchat, signaling = pathway, 
                                       width = 8, height = 2.5, font.size = 10)
    }, error = function(e) {
      log(paste("  Warning: Could not create centrality heatmap for pathway", pathway))
    })
  }
  dev.off()
  
  return(TRUE)
}

#' Plot 13: CellChat Signaling Aggregate Plots
#' 
#' Creates aggregate signaling plots for top pathways
#' 
#' @param cellchat CellChat object
#' @param plots_dir Output directory
#' @param n_pathways Number of top pathways to plot
create_cellchat_signaling_aggregate <- function(cellchat, plots_dir, n_pathways = 4) {
  pathways.show.agg <- head(cellchat@netP$pathways, n_pathways)
  
  if (length(pathways.show.agg) == 0) {
    return(FALSE)
  }
  
  pdf(file.path(plots_dir, "cellchat_signaling_roles_aggregate.pdf"), width = 12, height = 10)
  par(mfrow = c(min(2, length(pathways.show.agg)), min(2, length(pathways.show.agg))))
  for (i in 1:length(pathways.show.agg)) {
    netVisual_aggregate(cellchat, signaling = pathways.show.agg[i], layout = "circle")
  }
  dev.off()
  
  return(TRUE)
}

################################################################################
# NICHENET PLOTS
################################################################################

#' Plot 14: NicheNet Overview - Cell Types
#' 
#' Creates UMAP showing cell type distribution
#' 
#' @param seu Seurat object
#' @param plots_dir Output directory
#' @param base_size Base font size
#' @param title_size Title font size
#' @param legend_size Legend font size
#' @return ggplot object
create_nichenet_overview_celltypes <- function(seu, plots_dir, base_size = 16, 
                                               title_size = 20, legend_size = 14) {
  p <- DimPlot(seu, reduction = "umap", group.by = "ident", label = TRUE, repel = TRUE) +
    ggtitle("Cell Types in Dataset") +
    theme_minimal(base_size = base_size) +
    theme(plot.title = element_text(size = title_size, face = "bold"),
          legend.text = element_text(size = legend_size))
  
  save_dual_plot(p, "nichenet_overview_celltypes", plots_dir = plots_dir, width = 10, height = 8)
  
  return(p)
}

#' Plot 15: NicheNet Overview - Conditions
#' 
#' Creates UMAP showing experimental conditions
#' 
#' @param seu Seurat object
#' @param plots_dir Output directory
#' @param base_size Base font size
#' @param title_size Title font size
#' @param legend_size Legend font size
#' @return ggplot object or NULL if no condition field
create_nichenet_overview_conditions <- function(seu, plots_dir, base_size = 16, 
                                               title_size = 20, legend_size = 14) {
  if (!"condition" %in% colnames(seu@meta.data)) {
    return(NULL)
  }
  
  p <- DimPlot(seu, reduction = "umap", group.by = "condition", label = FALSE) +
    ggtitle("Experimental Conditions") +
    theme_minimal(base_size = base_size) +
    theme(plot.title = element_text(size = title_size, face = "bold"),
          legend.text = element_text(size = legend_size))
  
  save_dual_plot(p, "nichenet_overview_conditions", plots_dir = plots_dir, width = 10, height = 8)
  
  return(p)
}

#' Plot 16: NicheNet Ligand Activity Barplot
#' 
#' Creates barplot showing top predicted ligands for a receiver cell type
#' 
#' @param ligand_activities Data frame with ligand activities
#' @param receiver Receiver cell type name
#' @param plots_dir Output directory
#' @param n_top Number of top ligands to show
#' @param base_size Base font size
#' @param title_size Title font size
#' @param axis_title_size Axis title font size
#' @param legend_size Legend font size
#' @return ggplot object or NULL
create_nichenet_ligand_activity <- function(ligand_activities, receiver, plots_dir, 
                                            n_top = 20, base_size = 18, 
                                            title_size = 22, axis_title_size = 20, 
                                            legend_size = 16) {
  top_ligands <- ligand_activities %>% top_n(n_top, aupr_corrected) %>% pull(test_ligand)
  top_ligands_df <- ligand_activities %>% 
    filter(test_ligand %in% top_ligands) %>% 
    arrange(aupr_corrected)
  
  if (nrow(top_ligands_df) == 0) {
    return(NULL)
  }
  
  p <- ggplot(top_ligands_df, aes(x = reorder(test_ligand, aupr_corrected), 
                                  y = aupr_corrected, fill = aupr_corrected)) +
    geom_bar(stat = "identity", color = "white", size = 0.3) +
    coord_flip() +
    scale_fill_gradient(low = "#FFF5E6", high = "#FF6B35", name = "AUPR\n(Corrected)") +
    labs(title = paste("Top Predicted Ligands for", receiver), 
         x = "Ligand", y = "AUPR (Corrected)") +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.y = element_text(size = legend_size, face = "bold"),
      axis.text.x = element_text(size = legend_size),
      axis.title = element_text(size = axis_title_size, face = "bold"),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      legend.title = element_text(size = legend_size, face = "bold"),
      legend.text = element_text(size = legend_size),
      panel.grid.major.x = element_line(color = "gray90", size = 0.3),
      panel.grid.minor = element_blank()
    )
  
  save_dual_plot(p, paste0("nichenet_ligand_activity_", receiver), 
                plots_dir = plots_dir, width = 12, height = 10)
  
  return(p)
}

#' Plot 17: NicheNet Ligand-Receptor Network
#' 
#' Creates heatmap showing L-R pairs for top ligands
#' 
#' @param lr_network Ligand-receptor network data frame
#' @param top_ligands Vector of top ligand names
#' @param expressed_receptors Vector of expressed receptor names
#' @param receiver Receiver cell type name
#' @param plots_dir Output directory
#' @param base_size Base font size
#' @param title_size Title font size
#' @param axis_title_size Axis title font size
#' @param legend_size Legend font size
#' @return ggplot object or NULL
create_nichenet_lr_network <- function(lr_network, top_ligands, expressed_receptors, 
                                       receiver, plots_dir, base_size = 16, 
                                       title_size = 20, axis_title_size = 18, 
                                       legend_size = 14) {
  lr_network_top <- lr_network %>% 
    filter(from %in% top_ligands, to %in% expressed_receptors)
  
  if (nrow(lr_network_top) == 0) {
    return(NULL)
  }
  
  vis_lr_df <- lr_network_top %>% 
    group_by(from, to) %>% 
    summarize(val = 1, .groups = "drop")
  
  p <- ggplot(vis_lr_df, aes(x = to, y = from, fill = factor(val))) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("1" = "red"), guide = "none") +
    labs(title = paste("Ligand-Receptor pairs for", receiver), 
         x = "Receptor", y = "Ligand") +
    theme_minimal(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = legend_size),
          axis.text.y = element_text(size = legend_size),
          axis.title = element_text(size = axis_title_size),
          plot.title = element_text(size = title_size, face = "bold"))
  
  save_dual_plot(p, paste0("nichenet_lr_network_", receiver), 
                plots_dir = plots_dir, width = 12, height = 10)
  
  return(p)
}

#' Plot 18: NicheNet Target Network Heatmap
#' 
#' Creates heatmap showing ligand-target gene relationships
#' 
#' @param active_targets Data frame with ligand-target links
#' @param receiver Receiver cell type name
#' @param plots_dir Output directory
#' @param n_top Number of top targets to show per ligand
#' @param base_size Base font size
#' @param title_size Title font size
#' @param axis_title_size Axis title font size
#' @param legend_size Legend font size
#' @return ggplot object or NULL
create_nichenet_target_network <- function(active_targets, receiver, plots_dir, 
                                           n_top = 15, base_size = 16, 
                                           title_size = 20, axis_title_size = 18, 
                                           legend_size = 14) {
  if (is.null(active_targets) || nrow(active_targets) == 0) {
    return(NULL)
  }
  
  top_targets <- active_targets %>% 
    group_by(ligand) %>% 
    top_n(n_top, weight) %>% 
    ungroup() %>% 
    pull(target) %>% 
    unique()
  
  p <- ggplot(active_targets %>% filter(target %in% top_targets), 
             aes(x = target, y = ligand, fill = weight)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(title = paste("Ligand-Target links in", receiver), 
         x = "Target Gene", y = "Ligand") +
    theme_minimal(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = legend_size),
          axis.text.y = element_text(size = legend_size),
          axis.title = element_text(size = axis_title_size),
          plot.title = element_text(size = title_size, face = "bold"),
          legend.title = element_text(size = legend_size))
  
  save_dual_plot(p, paste0("nichenet_target_network_", receiver), 
                plots_dir = plots_dir, width = 14, height = 12)
  
  return(p)
}

#' Plot 19: NicheNet Circos Plot
#' 
#' Creates circos plot showing ligand-target communication flow
#' 
#' @param active_targets Data frame with ligand-target links
#' @param receiver Receiver cell type name
#' @param plots_dir Output directory
#' @param n_ligands Number of top ligands to include
#' @param n_targets Number of top targets to include
#' @param log Logging function
#' @return TRUE if successful, FALSE otherwise
create_nichenet_circos <- function(active_targets, receiver, condition = NULL, plots_dir, 
                                   n_ligands = 10, n_targets = 30, log = cat) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    log("    circlize package not available - skipping circos plot")
    return(FALSE)
  }
  
  if (is.null(active_targets) || nrow(active_targets) == 0) {
    log("    No active targets found for circos plot")
    return(FALSE)
  }
  
  library(circlize)
  
  # Get top targets
  top_targets_circos <- active_targets %>% 
    top_n(min(n_targets, nrow(active_targets)), weight) %>% 
    pull(target) %>% 
    unique()
  
  circos_links <- active_targets %>% 
    filter(target %in% top_targets_circos) %>%
    select(ligand, target, weight) %>%
    distinct() %>%
    filter(!is.na(ligand), !is.na(target), !is.na(weight), weight > 0)
  
  if (nrow(circos_links) == 0 || 
      length(unique(circos_links$ligand)) == 0 || 
      length(unique(circos_links$target)) == 0) {
    log("    No ligand-target links for circos plot")
    return(FALSE)
  }
  
  # Build filename with condition if provided
  file_name <- if (!is.null(condition)) {
    paste0("nichenet_circos_", condition, "_", receiver, ".pdf")
  } else {
    paste0("nichenet_circos_", receiver, ".pdf")
  }
  
  pdf(file.path(plots_dir, file_name), 
      width = 10, height = 10)
  
  circos.clear()
  
  chordDiagram(circos_links, 
              transparency = 0.5,
              directional = 1,
              direction.type = "arrows",
              link.arr.type = "big.arrow",
              annotationTrack = "grid",
              preAllocateTracks = list(track.height = 0.1))
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], 
               CELL_META$sector.index, 
               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
  }, bg.border = NA)
  
  plot_title <- if (!is.null(condition)) {
    paste("NicheNet Communication Flow for", receiver, "in", condition)
  } else {
    paste("NicheNet Communication Flow for", receiver)
  }
  title(plot_title)
  
  circos.clear()
  dev.off()
  
  log("    Saved circos plot:", file_name)
  return(TRUE)
}

#' Plot 20: NicheNet Legacy Target Heatmap
#' 
#' Creates heatmap for best ligand's target genes (legacy version)
#' 
#' @param active_targets Data frame with ligand-target links
#' @param best_ligand Name of best ligand
#' @param receiver Receiver cell type name
#' @param plots_dir Output directory
#' @param n_top Number of top targets to show
#' @param base_size Base font size
#' @return ggplot object or NULL
create_nichenet_legacy_heatmap <- function(active_targets, best_ligand, receiver, 
                                           plots_dir, n_top = 30, base_size = 16) {
  if (nrow(active_targets) <= 5) {
    return(NULL)
  }
  
  top_targets <- active_targets %>% top_n(n_top, weight)
  
  p <- ggplot(top_targets, aes(x = ligand, y = target, fill = weight)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(title = paste(best_ligand, "Target Genes in", receiver), 
         x = "Ligand", y = "Target Gene") +
    theme_minimal(base_size = base_size)
  
  save_dual_plot(p, paste0("nichenet_targets_heatmap_", receiver, "_", best_ligand), 
                plots_dir = plots_dir, width = 6, height = 9)
  
  return(p)
}

#' Plot 21: NicheNet Summary Heatmap (All Receivers)
#' 
#' Creates summary heatmap showing top ligands across all receiver cell types
#' 
#' @param combined_activities Data frame with combined ligand activities
#' @param plots_dir Output directory
#' @param n_top Number of top ligands per receiver
#' @param base_size Base font size
#' @param title_size Title font size
#' @param axis_title_size Axis title font size
#' @param legend_size Legend font size
#' @return ggplot object
create_nichenet_summary_heatmap <- function(combined_activities, plots_dir, 
                                            n_top = 10, base_size = 16, 
                                            title_size = 22, axis_title_size = 20, 
                                            legend_size = 16, condition = NULL) {
  top_ligands_all <- combined_activities %>%
    group_by(receiver) %>%
    top_n(n_top, aupr_corrected) %>%
    pull(test_ligand) %>%
    unique()
  
  heatmap_data <- combined_activities %>%
    filter(test_ligand %in% top_ligands_all)
  
  plot_title <- if (!is.null(condition) && nzchar(condition)) {
    paste0("NicheNet: Top Predicted Ligands across Receivers\n", condition)
  } else {
    "NicheNet: Top Predicted Ligands across Receivers"
  }

  p <- ggplot(heatmap_data, aes(x = receiver, y = test_ligand, fill = aupr_corrected)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient(low = "#FFF5E6", high = "#FF6B35", name = "AUPR\n(Corrected)",
                       limits = c(0, max(heatmap_data$aupr_corrected, na.rm = TRUE))) +
    labs(title = plot_title, 
         x = "Receiver Cell Type", y = "Ligand") +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = legend_size, face = "bold"),
      axis.text.y = element_text(size = legend_size),
      axis.title = element_text(size = axis_title_size, face = "bold"),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      legend.title = element_text(size = legend_size, face = "bold"),
      legend.text = element_text(size = legend_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  file_name <- if (!is.null(condition) && nzchar(condition)) {
    paste0("nichenet_summary_heatmap_", gsub("[^A-Za-z0-9_-]+", "_", condition))
  } else {
    "nichenet_summary_heatmap"
  }

  save_dual_plot(p, file_name, plots_dir = plots_dir, width = 12, height = 13)
  
  return(p)
}

#' Create comprehensive ligand-receptor interaction heatmap
#' 
#' Creates heatmap showing ligand-receptor pairs with expression levels and interaction strength
#' 
#' @param seu Seurat object
#' @param lr_network Ligand-receptor network database
#' @param all_ligand_activities List of ligand activities per cell type
#' @param plots_dir Output directory
#' @param n_top Number of top ligands per cell type
#' @param base_size Base font size
create_nichenet_lr_interaction_heatmap <- function(seu, lr_network, all_ligand_activities, 
                                                    plots_dir, n_top = 10, base_size = 12) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    message("pheatmap package not available - skipping LR interaction heatmap")
    return(invisible(NULL))
  }
  
  if (length(all_ligand_activities) == 0) {
    message("No ligand activities available for LR heatmap")
    return(invisible(NULL))
  }
  
  library(pheatmap)
  
  tryCatch({
    # Get top ligands across all receivers
    combined_activities <- bind_rows(all_ligand_activities)
    top_ligands <- combined_activities %>%
      group_by(receiver) %>%
      top_n(n_top, aupr_corrected) %>%
      pull(test_ligand) %>%
      unique()
    
    # Get expressed receptors
    expr_data <- GetAssayData(seu, layer = "data")
    all_receptors <- lr_network %>% pull(to) %>% unique()
    expressed_receptors <- intersect(all_receptors, rownames(seu))
    expressed_receptors <- expressed_receptors[Matrix::rowSums(expr_data[expressed_receptors, , drop=FALSE] > 0) > ncol(seu) * 0.05]
    
    # Filter LR network to top ligands and expressed receptors
    lr_subset <- lr_network %>%
      filter(from %in% top_ligands, to %in% expressed_receptors) %>%
      distinct(from, to)
    
    if (nrow(lr_subset) == 0) {
      message("No ligand-receptor pairs found for heatmap")
      return(invisible(NULL))
    }
    
    # Calculate average expression per cell type
    cell_types <- unique(Idents(seu))
    
    # Create matrix: rows = LR pairs, columns = cell types
    lr_pairs <- paste(lr_subset$from, lr_subset$to, sep = "_")
    expr_matrix <- matrix(0, nrow = length(lr_pairs), ncol = length(cell_types))
    rownames(expr_matrix) <- lr_pairs
    colnames(expr_matrix) <- cell_types
    
    for (i in 1:nrow(lr_subset)) {
      ligand <- lr_subset$from[i]
      receptor <- lr_subset$to[i]
      pair_name <- paste(ligand, receptor, sep = "_")
      
      for (ct in cell_types) {
        cells <- WhichCells(seu, idents = ct)
        if (length(cells) > 0) {
          # Average of ligand and receptor expression
          lig_expr <- if (ligand %in% rownames(expr_data)) {
            mean(expr_data[ligand, cells])
          } else { 0 }
          rec_expr <- if (receptor %in% rownames(expr_data)) {
            mean(expr_data[receptor, cells])
          } else { 0 }
          expr_matrix[pair_name, ct] <- (lig_expr + rec_expr) / 2
        }
      }
    }
    
    # Remove rows with all zeros
    expr_matrix <- expr_matrix[rowSums(expr_matrix) > 0, , drop = FALSE]
    
    if (nrow(expr_matrix) == 0) {
      message("No expressed LR pairs found")
      return(invisible(NULL))
    }
    
    # Limit to top 50 pairs by variance
    if (nrow(expr_matrix) > 50) {
      row_vars <- apply(expr_matrix, 1, var)
      top_pairs <- names(sort(row_vars, decreasing = TRUE)[1:50])
      expr_matrix <- expr_matrix[top_pairs, , drop = FALSE]
    }
    
    # Create heatmap
    png(file.path(plots_dir, "nichenet_lr_interaction_heatmap.png"),
        width = 12, height = max(10, nrow(expr_matrix) * 0.2), units = "in", res = 300)
    
    pheatmap(expr_matrix,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             scale = "row",
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
             main = "Ligand-Receptor Expression Across Cell Types",
             fontsize = base_size,
             fontsize_row = 8,
             fontsize_col = 10,
             angle_col = 45)
    
    dev.off()
    
    message("Saved ligand-receptor interaction heatmap")
    
  }, error = function(e) {
    message("ERROR creating LR interaction heatmap: ", e$message)
  })
  
  invisible(NULL)
}

#' Plot 35: NicheNet Circos Plot per Condition
#' 
#' Creates a combined circos plot for each condition showing all
#' ligand-receptor-target interactions across receiver cell types
#' 
#' @param seu Seurat object
#' @param condition_name Name of the condition
#' @param lr_network Ligand-receptor network data frame
#' @param all_ligand_activities List of ligand activities per receiver
#' @param plots_dir Output directory
#' @param n_top_ligands Number of top ligands per receiver to include
#' @param n_top_targets Number of top targets per ligand to include
#' @param log Logging function
#' @return TRUE if successful, FALSE otherwise
create_nichenet_condition_circos <- function(seu, condition_name, lr_network, 
                                             all_ligand_activities, plots_dir, 
                                             n_top_ligands = 8, n_top_targets = 20, 
                                             log = cat) {
  tryCatch({
    if (!requireNamespace("circlize", quietly = TRUE)) {
      log(sprintf("  circlize package not available - skipping condition circos for %s", condition_name))
      return(FALSE)
    }
    
    library(circlize)
    
    # Filter activities for this condition
    # Debug: check what activity names exist
    log(sprintf("  DEBUG: Looking for condition '%s' in %d activities", condition_name, length(all_ligand_activities)))
    log(sprintf("  DEBUG: Activity names: %s", paste(names(all_ligand_activities)[1:min(3, length(all_ligand_activities))], collapse=", ")))
    
    # More robust condition filtering
    cond_activities <- all_ligand_activities[grepl(paste0("^", condition_name, "_"), 
                                                     names(all_ligand_activities))]
    
    if (length(cond_activities) == 0) {
      # Try alternative naming patterns
      cond_activities <- all_ligand_activities[grepl(condition_name, names(all_ligand_activities))]
    }
    
    if (length(cond_activities) == 0) {
      log(sprintf("  WARNING: No ligand activities found for condition '%s'", condition_name))
      log(sprintf("  Available names: %s", paste(names(all_ligand_activities), collapse=", ")))
      return(FALSE)
    }
    
    log(sprintf("  Found %d activity entries for condition %s", length(cond_activities), condition_name))
    
    # Collect all L-R-T interactions for this condition
    all_interactions <- data.frame()
    
    for (activity_name in names(cond_activities)) {
      receiver <- sub(paste0("^", condition_name, "_"), "", activity_name)
      ligand_df <- cond_activities[[activity_name]]
      
      if (!is.data.frame(ligand_df) || nrow(ligand_df) == 0) {
        log(sprintf("    WARNING: Empty or invalid ligand data for receiver %s in condition %s", receiver, condition_name))
        next
      }
      
      top_ligands <- head(ligand_df$test_ligand, n_top_ligands)
      
      for (ligand in top_ligands) {
        # Get all receptors for this ligand, ranked by specificity/weight
        receptors_data <- lr_network %>%
          filter(from == ligand) %>%
          arrange(desc(weight)) %>%  # Sort by weight if available
          pull(to)
        
        if (length(receptors_data) > 0) {
          # Use multiple receptors for more diverse interactions
          n_receptors_to_use <- min(3, length(receptors_data))
          for (i in 1:n_receptors_to_use) {
            all_interactions <- rbind(all_interactions, 
                                     data.frame(
                                       ligand = ligand,
                                       receptor = receptors_data[i],  # Use top receptors
                                       receiver = receiver,
                                       stringsAsFactors = FALSE
                                     ))
          }
        }
      }
    }
    
    if (nrow(all_interactions) == 0) {
      log(sprintf("  No L-R-T interactions found for condition %s", condition_name))
      return(FALSE)
    }
    
    # Remove duplicates and limit for visualization
    all_interactions <- all_interactions %>% distinct()
    if (nrow(all_interactions) > 100) {
      all_interactions <- all_interactions %>% slice_sample(n = 100)
    }
    
    # Create interaction data for circos
    # Format: source -> target
    from_to <- data.frame(
      from = c(
        all_interactions$ligand,
        all_interactions$receptor
      ),
      to = c(
        all_interactions$receptor,
        all_interactions$receiver
      ),
      stringsAsFactors = FALSE
    )
    
    # Get unique elements for color assignment
    all_elements <- unique(c(from_to$from, from_to$to))
    
    # Assign colors: ligands = orange, receptors = green, receivers = blue
    colors <- structure(
      c(
        rep("#FF8C42", length(unique(all_interactions$ligand))),      # Orange for ligands
        rep("#52B788", length(unique(all_interactions$receptor))),    # Green for receptors
        rep("#4D94FF", length(unique(all_interactions$receiver)))     # Blue for receivers
      ),
      names = all_elements
    )
    
    # Create circos plot
    png(file.path(plots_dir, paste0("nichenet_circos_", condition_name, ".png")),
        width = 14, height = 14, units = "in", res = 150)
    
    circlize::circos.clear()
    circos.initialize(factors = all_elements, 
                     xlim = c(0, 1),
                     sector.width = rep(1/length(all_elements), length(all_elements)))
    
    circos.trackPlotRegion(factors = all_elements, ylim = c(0, 1),
                          track.height = 0.3,
                          panel.fun = function(x, y) {
                            circos.text(CELL_META$xcenter, CELL_META$ylim[2] + 0.1, 
                                      CELL_META$sector.index, cex = 0.8, adj = c(0.5, 0))
                          })
    
    # Draw links
    for (i in 1:nrow(from_to)) {
      circos.link(from_to$from[i], 0, from_to$to[i], 0,
                 col = adjustcolor(colors[from_to$from[i]], alpha.f = 0.4),
                 border = NA)
    }
    
    # Add title
    title(main = paste("NicheNet Ligand-Receptor-Target Network:", condition_name),
          cex.main = 1.8, font.main = 2)
    
    dev.off()
    circlize::circos.clear()
    
    log(sprintf("  Saved circos plot for condition: %s (%d interactions)", 
               condition_name, nrow(all_interactions)))
    
    return(TRUE)
    
  }, error = function(e) {
    log(sprintf("  Warning: Could not create circos plot for %s: %s", condition_name, e$message))
    if (dev.cur() != 1) dev.off()
    circlize::circos.clear()
    return(FALSE)
  })
}

#' Plot 36: NicheNet Bubble Plot per Condition
#' 
#' Creates bubble plot showing ligand-receptor interactions per condition
#' with sender (ligands) on top axis and receiver cell types represented
#' 
#' @param seu Seurat object
#' @param condition_name Name of the condition
#' @param all_ligand_activities List of ligand activities per receiver for this condition
#' @param lr_network Ligand-receptor network data frame
#' @param plots_dir Output directory
#' @param n_top_ligands Number of top ligands to show
#' @param base_size Base font size
#' @param title_size Title font size
#' @param axis_title_size Axis title font size
#' @param legend_size Legend font size
#' @param log Logging function
#' @return TRUE if successful, FALSE otherwise
create_nichenet_condition_bubble <- function(seu, condition_name, all_ligand_activities, 
                                            lr_network, plots_dir,
                                            n_top_ligands = 15, base_size = 14,
                                            title_size = 20, axis_title_size = 16,
                                            legend_size = 12, log = cat) {
  tryCatch({
    # Filter activities for this condition
    cond_activities <- all_ligand_activities[grepl(paste0("^", condition_name, "_"), 
                                                     names(all_ligand_activities))]
    
    if (length(cond_activities) == 0) {
      log(sprintf("  No ligand activities found for condition %s", condition_name))
      return(FALSE)
    }
    
    # Combine all receiver activities for this condition
    cond_data <- bind_rows(cond_activities)
    
    if (nrow(cond_data) == 0) {
      log(sprintf("  No data to plot for condition %s", condition_name))
      return(FALSE)
    }
    
    # Get top ligands overall in this condition
    top_ligands <- cond_data %>%
      group_by(test_ligand) %>%
      summarize(mean_aupr = mean(aupr_corrected, na.rm = TRUE), .groups = "drop") %>%
      top_n(n_top_ligands, mean_aupr) %>%
      pull(test_ligand)
    
    # Get receptor info for top ligands
    receptors_for_ligands <- lr_network %>%
      filter(from %in% top_ligands) %>%
      select(from, to) %>%
      distinct() %>%
      left_join(
        cond_data %>% select(test_ligand, receiver, aupr_corrected) %>% distinct(),
        by = c("from" = "test_ligand")
      )
    
    if (nrow(receptors_for_ligands) == 0) {
      log(sprintf("  No ligand-receptor pairs found for condition %s", condition_name))
      return(FALSE)
    }
    
    # Create bubble data
    bubble_data <- receptors_for_ligands %>%
      rename(ligand = from, receptor = to) %>%
      group_by(ligand, receptor, receiver) %>%
      summarize(aupr = mean(aupr_corrected, na.rm = TRUE), .groups = "drop") %>%
      filter(!is.na(aupr), aupr > 0)
    
    if (nrow(bubble_data) == 0) {
      log(sprintf("  No valid interactions for bubble plot in condition %s", condition_name))
      return(FALSE)
    }
    
    # Create LR pair labels
    bubble_data <- bubble_data %>%
      mutate(lr_pair = paste0(ligand, " - ", receptor))
    
    # Create bubble plot
    p <- ggplot(bubble_data, aes(x = receiver, y = lr_pair, size = aupr, fill = aupr)) +
      geom_point(shape = 21, color = "black", stroke = 0.5) +
      scale_size_continuous(name = "AUPR\n(Corrected)", range = c(2, 8)) +
      scale_fill_gradient(low = "#FFF5E6", high = "#FF6B35", name = "AUPR\n(Corrected)") +
      labs(
        title = paste("NicheNet: Ligand-Receptor-Receiver Interactions in", condition_name),
        x = "Receiver Cell Type",
        y = "Ligand - Receptor Interaction"
      ) +
      theme_minimal(base_size = base_size) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = legend_size, face = "bold"),
        axis.text.y = element_text(size = legend_size - 2),
        axis.title = element_text(size = axis_title_size, face = "bold"),
        plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
        legend.title = element_text(size = legend_size, face = "bold"),
        legend.text = element_text(size = legend_size),
        panel.grid.major.y = element_line(color = "gray90", size = 0.3),
        panel.grid.minor = element_blank()
      )
    
    # Save plot
    file_name <- paste0("nichenet_bubble_", gsub("[^A-Za-z0-9_-]+", "_", condition_name))
    save_dual_plot(p, file_name, plots_dir = plots_dir, width = 14, height = max(10, nrow(bubble_data) * 0.3))
    
    log(sprintf("  Saved bubble plot for condition: %s (%d interactions)", 
               condition_name, nrow(bubble_data)))
    
    return(TRUE)
    
  }, error = function(e) {
    log(sprintf("  Warning: Could not create bubble plot for %s: %s", condition_name, e$message))
    return(FALSE)
  })
}
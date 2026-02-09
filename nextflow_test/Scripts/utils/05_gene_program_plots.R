################################################################################
# plots.R
# Plotting functions for Module 05 - Gene Program Discovery
################################################################################

# Configuration parameters
CFG_RANDOM_SEED <- 42
CFG_BASE_SIZE <- 18
CFG_TITLE_SIZE <- 24
CFG_AXIS_TITLE_SIZE <- 20
CFG_AXIS_TEXT_SIZE <- 16
CFG_LEGEND_SIZE <- 18
CFG_PLOT_WIDTH <- 12
CFG_PLOT_HEIGHT <- 8
CFG_PLOT_DPI <- 150

################################################################################
# CONDITION COMPARISON PLOTS
################################################################################

#' Compare MP scores by cell type for two conditions and plot significance
#'
#' @param seu Seurat object with MP scores
#' @param mp.genes List of meta-program gene sets
#' @param condition_field Condition metadata field
#' @param plots_dir Output directory for plots
#' @param data_dir Output directory for tables
#' @param gsea_results Optional list of GSEA results per MP
#' @param log Logging function
create_mp_celltype_condition_comparison <- function(seu, mp.genes, condition_field, plots_dir, data_dir,
                                                    gsea_results = NULL, log = cat) {
  tryCatch({
    if (!"cell_type" %in% colnames(seu@meta.data)) {
      log("Warning: cell_type not found in metadata, skipping MP cell type comparison")
      return()
    }

    conditions <- seu@meta.data[[condition_field]]
    if (is.factor(conditions)) {
      condition_levels <- levels(conditions)
    } else {
      condition_levels <- unique(as.character(conditions))
    }

    if (length(condition_levels) < 2) {
      log("Warning: Need at least two conditions for MP comparison")
      return()
    }

    cond_a <- condition_levels[1]
    cond_b <- condition_levels[2]
    log("MP comparison: ", cond_a, " vs ", cond_b)

    mp_names <- names(mp.genes)
    mp_names <- mp_names[mp_names %in% colnames(seu@meta.data)]
    if (length(mp_names) == 0) {
      log("Warning: No MP score columns found in metadata")
      return()
    }

    cell_types <- unique(seu@meta.data$cell_type)

    results <- list()
    idx <- 1
    for (ct in cell_types) {
      for (mp in mp_names) {
        scores_a <- seu@meta.data[seu@meta.data$cell_type == ct &
                                   seu@meta.data[[condition_field]] == cond_a, mp]
        scores_b <- seu@meta.data[seu@meta.data$cell_type == ct &
                                   seu@meta.data[[condition_field]] == cond_b, mp]

        scores_a <- scores_a[!is.na(scores_a)]
        scores_b <- scores_b[!is.na(scores_b)]

        p_val <- NA_real_
        if (length(scores_a) >= 3 && length(scores_b) >= 3) {
          p_val <- tryCatch(wilcox.test(scores_a, scores_b)$p.value, error = function(e) NA_real_)
        }

        results[[idx]] <- data.frame(
          cell_type = ct,
          MP = mp,
          condition_a = cond_a,
          condition_b = cond_b,
          mean_a = mean(scores_a, na.rm = TRUE),
          mean_b = mean(scores_b, na.rm = TRUE),
          delta = mean(scores_b, na.rm = TRUE) - mean(scores_a, na.rm = TRUE),
          p_value = p_val,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1
      }
    }

    comp_df <- dplyr::bind_rows(results)
    comp_df$p_adj <- p.adjust(comp_df$p_value, method = "BH")
    comp_df$significant <- !is.na(comp_df$p_adj) & comp_df$p_adj < 0.05

    # Signature labels from GSEA if available
    signature_map <- setNames(rep("Unlabeled", length(mp_names)), mp_names)

    if (is.null(gsea_results)) {
      gsea_file <- file.path(data_dir, "all_MP_GSEA_top10.csv")
      if (file.exists(gsea_file)) {
        gsea_df <- read.csv(gsea_file, stringsAsFactors = FALSE)
        if ("MP" %in% colnames(gsea_df)) {
          for (mp in mp_names) {
            mp_rows <- gsea_df[gsea_df$MP == mp, , drop = FALSE]
            if (nrow(mp_rows) > 0) {
              if ("padj" %in% colnames(mp_rows)) {
                mp_rows <- mp_rows[order(mp_rows$padj), , drop = FALSE]
              }
              label_col <- intersect(c("pathway", "term", "Description"), colnames(mp_rows))
              if (length(label_col) > 0) {
                signature_map[mp] <- as.character(mp_rows[[label_col[1]]][1])
              }
            }
          }
        }
      }
    } else {
      for (mp in mp_names) {
        mp_res <- gsea_results[[mp]]
        if (!is.null(mp_res) && nrow(mp_res) > 0) {
          label_col <- intersect(c("pathway", "term", "Description"), colnames(mp_res))
          if (length(label_col) > 0) {
            mp_sorted <- mp_res
            if ("padj" %in% colnames(mp_res)) {
              mp_sorted <- mp_res[order(mp_res$padj), , drop = FALSE]
            }
            signature_map[mp] <- as.character(mp_sorted[[label_col[1]]][1])
          }
        }
      }
    }

    top3_map <- sapply(mp_names, function(mp) {
      genes <- mp.genes[[mp]]
      if (is.null(genes)) return("NA")
      paste(head(genes, 3), collapse = ", ")
    })

    label_map <- paste0(mp_names, "\n", signature_map[mp_names], "\n", top3_map[mp_names])
    comp_df$MP_label <- label_map[comp_df$MP]

    # Save label table
    label_df <- data.frame(
      MP = mp_names,
      signature = signature_map[mp_names],
      top3_genes = top3_map[mp_names],
      label = label_map,
      stringsAsFactors = FALSE
    )
    write.csv(label_df, file.path(data_dir, "MP_signature_top3_genes.csv"), row.names = FALSE)
    write.csv(comp_df, file.path(data_dir, "MP_celltype_condition_comparison.csv"), row.names = FALSE)

    # Plot heatmap of delta with significance
    p_cmp <- ggplot(comp_df, aes(x = MP_label, y = cell_type, fill = delta)) +
      geom_tile(color = "white", linewidth = 0.2) +
      geom_text(aes(label = ifelse(significant, "*", "")), size = 4) +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                           name = paste0("Δ score (", cond_b, " - ", cond_a, ")")) +
      labs(title = paste("MP Scores: ", cond_b, " vs ", cond_a, " (per Cell Type)", sep = ""),
           subtitle = "* FDR < 0.05",
           x = "Meta-Program (Signature + Top 3 Genes)",
           y = "Cell Type") +
      theme_bw(base_size = CFG_BASE_SIZE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = CFG_AXIS_TEXT_SIZE),
            axis.text.y = element_text(size = CFG_AXIS_TEXT_SIZE))

    ggsave(file.path(plots_dir, "MP_scores_celltype_condition_comparison.png"),
           p_cmp, width = 14, height = max(6, length(cell_types) * 0.4), dpi = CFG_PLOT_DPI)

    log("Saved MP cell type condition comparison plot")
  }, error = function(e) {
    log("Warning: MP cell type condition comparison failed:", e$message)
  })
}

#' Create violin plots of meta-program scores by condition
#'
#' @param seu Seurat object with MP scores
#' @param mp.genes List of meta-program gene sets
#' @param condition_field Condition metadata field
#' @param plots_dir Output directory
#' @param log Logging function
create_mp_violin_by_condition <- function(seu, mp.genes, condition_field, plots_dir, log = cat) {
  tryCatch({
    DefaultAssay(seu) <- "RNA"
    
    vln_plots <- list()
    for (mp in names(mp.genes)) {
      if (mp %in% colnames(seu@meta.data)) {
        vln_plots[[mp]] <- VlnPlot(seu, 
                                   features = mp, 
                                   group.by = condition_field,
                                   pt.size = 0) + 
          NoLegend() +
          theme(axis.title.x = element_blank())
      }
    }
    
    if (length(vln_plots) > 0) {
      p_violin <- wrap_plots(vln_plots, ncol = ceiling(sqrt(length(vln_plots))))
      
      ggsave(file.path(plots_dir, "MP_scores_by_condition.png"), 
             p_violin, width = 14, height = 10, dpi = CFG_PLOT_DPI)
      
      log("Saved violin plots per condition")
    }
  }, error = function(e) {
    log("Warning: Violin plots failed:", e$message)
  })
}

#' Create line plots showing MP activity across cell types
#'
#' @param seu Seurat object with MP scores
#' @param mp.genes List of meta-program gene sets
#' @param plots_dir Output directory
#' @param log Logging function
create_mp_celltype_lines <- function(seu, mp.genes, plots_dir, log = cat) {
  tryCatch({
    mp_celltype_data <- seu@meta.data %>%
      select(cell_type, all_of(names(mp.genes))) %>%
      pivot_longer(cols = all_of(names(mp.genes)), 
                  names_to = "MetaProgram", 
                  values_to = "Score") %>%
      group_by(cell_type, MetaProgram) %>%
      summarize(mean_score = mean(Score, na.rm = TRUE), .groups = "drop")
    
    p_line <- ggplot(mp_celltype_data, aes(x = cell_type, y = mean_score, 
                                            color = MetaProgram, group = MetaProgram)) +
      geom_line(size = 1.2) +
      geom_point(size = 3) +
      labs(title = "Meta-Program Activity Across Cell Types",
           x = "Cell Type", y = "Mean Score") +
      theme_bw(base_size = CFG_BASE_SIZE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = CFG_AXIS_TEXT_SIZE),
            axis.text.y = element_text(size = CFG_AXIS_TEXT_SIZE),
            legend.position = "right") +
      scale_color_brewer(palette = "Set3")
    
    ggsave(file.path(plots_dir, "MP_scores_by_celltype_lines.png"),
           p_line, width = 12, height = 6, dpi = CFG_PLOT_DPI)
    
    log("Saved line plot across cell types")
  }, error = function(e) {
    log("Warning: Line plot failed:", e$message)
  })
}

#' Create stacked barplot of cell type composition
#'
#' @param seu Seurat object
#' @param condition_field Condition metadata field
#' @param plots_dir Output directory
#' @param log Logging function
create_celltype_composition_stacked <- function(seu, condition_field, plots_dir, log = cat) {
  tryCatch({
    comp_data <- seu@meta.data %>%
      group_by(!!sym(condition_field), cell_type) %>%
      summarize(n = n(), .groups = "drop") %>%
      group_by(!!sym(condition_field)) %>%
      mutate(pct = n / sum(n) * 100)
    
    p_stacked <- ggplot(comp_data, aes(x = !!sym(condition_field), y = pct, fill = cell_type)) +
      geom_bar(stat = "identity") +
      labs(title = "Cell Type Composition by Condition",
           x = "Condition", y = "Percentage") +
      theme_bw(base_size = CFG_BASE_SIZE) +
      scale_fill_brewer(palette = "Set3") +
      theme(axis.text.x = element_text(size = CFG_AXIS_TEXT_SIZE),
            axis.text.y = element_text(size = CFG_AXIS_TEXT_SIZE),
            legend.position = "right")
    
    ggsave(file.path(plots_dir, "celltype_composition_stacked.png"),
           p_stacked, width = 8, height = 6, dpi = CFG_PLOT_DPI)
    
    log("Saved stacked composition barplot")
  }, error = function(e) {
    log("Warning: Stacked barplot failed:", e$message)
  })
}

#' Create fold-change barplot for MP activity comparison
#'
#' @param comp_df Data frame with comparison results
#' @param conditions Vector of condition names
#' @param plots_dir Output directory
#' @param log Logging function
create_mp_fold_change_plot <- function(comp_df, conditions, plots_dir, log = cat) {
  tryCatch({
    comp_df$MP <- factor(comp_df$MP, levels = comp_df$MP[order(comp_df$fold_change)])
    
    p_fc <- ggplot(comp_df, aes(x = MP, y = log2(fold_change), fill = significant)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      coord_flip() +
      scale_fill_manual(values = c("grey60", "red")) +
      labs(title = paste("Meta-Program Activity:", conditions[2], "vs", conditions[1]),
           x = "Meta-Program",
           y = "log2(Fold Change)") +
      theme_bw(base_size = CFG_BASE_SIZE) +
      theme(axis.text.x = element_text(size = CFG_AXIS_TEXT_SIZE),
            axis.text.y = element_text(size = CFG_AXIS_TEXT_SIZE))
    
    ggsave(file.path(plots_dir, "MP_fold_change_comparison.png"), 
           p_fc, width = 8, height = 6, dpi = CFG_PLOT_DPI)
    
    log("Saved fold change comparison plot")
  }, error = function(e) {
    log("Warning: Fold change plot failed:", e$message)
  })
}

################################################################################
# UMAP VISUALIZATION PLOTS
################################################################################

#' Create UMAP grid for all meta-programs
#'
#' @param seu Seurat object with MP scores
#' @param mp.genes List of meta-program gene sets
#' @param plots_dir Output directory
#' @param log Logging function
create_mp_umap_grid <- function(seu, mp.genes, plots_dir, log = cat) {
  tryCatch({
    mp_plots <- list()
    for (mp_name in names(mp.genes)) {
      if (mp_name %in% colnames(seu@meta.data)) {
        p <- FeaturePlot(seu, features = mp_name, reduction = "umap") +
          scale_color_viridis_c(option = "magma") +
          ggtitle(mp_name) +
          theme_void() +
          theme(plot.title = element_text(size = CFG_LEGEND_SIZE, face = "bold", hjust = 0.5),
                legend.position = "right",
                legend.key.size = unit(0.3, "cm"))
        
        mp_plots[[mp_name]] <- p
      }
    }
    
    if (length(mp_plots) > 0) {
      n_plots <- length(mp_plots)
      ncols <- min(4, ceiling(sqrt(n_plots)))
      
      p_grid <- wrap_plots(mp_plots, ncol = ncols)
      
      ggsave(file.path(plots_dir, "MP_UMAP_grid_all.png"),
             p_grid, width = ncols * 4, height = ceiling(n_plots / ncols) * 3.5, dpi = CFG_PLOT_DPI)
      
      log("Saved UMAP grid for all meta-programs")
    }
  }, error = function(e) {
    log("Warning: UMAP grid failed:", e$message)
  })
}

#' Create UMAP colored by cell types
#'
#' @param seu Seurat object
#' @param plots_dir Output directory
#' @param log Logging function
create_umap_celltypes <- function(seu, plots_dir, log = cat) {
  tryCatch({
    p_celltype <- DimPlot(seu, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
      ggtitle("Cell Types") +
      theme_void() +
      theme(plot.title = element_text(size = CFG_TITLE_SIZE, face = "bold"))
    
    ggsave(file.path(plots_dir, "UMAP_celltypes.png"),
           p_celltype, width = 10, height = 8, dpi = CFG_PLOT_DPI)
    
    log("Saved cell type UMAP")
  }, error = function(e) {
    log("Warning: Cell type UMAP failed:", e$message)
  })
}

#' Create UMAP colored by condition
#'
#' @param seu Seurat object
#' @param condition_field Condition metadata field
#' @param plots_dir Output directory
#' @param log Logging function
create_umap_conditions <- function(seu, condition_field, plots_dir, log = cat) {
  if (!is.null(condition_field) && condition_field %in% colnames(seu@meta.data)) {
    tryCatch({
      p_condition <- DimPlot(seu, reduction = "umap", group.by = condition_field) +
        ggtitle("Experimental Conditions") +
        theme_void() +
        theme(plot.title = element_text(size = CFG_TITLE_SIZE, face = "bold"))
      
      ggsave(file.path(plots_dir, "UMAP_conditions.png"),
             p_condition, width = 10, height = 8, dpi = CFG_PLOT_DPI)
      
      log("Saved condition UMAP")
    }, error = function(e) {
      log("Warning: Condition UMAP failed:", e$message)
    })
  }
}

#' Create combined MP activity UMAP
#'
#' @param seu Seurat object with MP scores
#' @param mp.genes List of meta-program gene sets
#' @param plots_dir Output directory
#' @param log Logging function
#' @return Updated Seurat object with MP_combined score
create_umap_mp_combined <- function(seu, mp.genes, plots_dir, log = cat) {
  tryCatch({
    mp_cols <- names(mp.genes)[names(mp.genes) %in% colnames(seu@meta.data)]
    if (length(mp_cols) > 0) {
      seu$MP_combined <- rowMeans(seu@meta.data[, mp_cols], na.rm = TRUE)
      
      p_combined <- FeaturePlot(seu, features = "MP_combined", reduction = "umap") +
        scale_color_viridis_c(option = "magma") +
        ggtitle("Combined Meta-Program Activity") +
        theme_void() +
        theme(plot.title = element_text(size = CFG_TITLE_SIZE, face = "bold"))
      
      ggsave(file.path(plots_dir, "UMAP_MP_combined.png"),
             p_combined, width = 10, height = 8, dpi = CFG_PLOT_DPI)
      
      log("Saved combined MP activity UMAP")
    }
  }, error = function(e) {
    log("Warning: Combined MP UMAP failed:", e$message)
  })
  
  return(seu)
}

################################################################################
# METAPROGRAM METRICS AND GSEA VISUALIZATION
################################################################################

#' Create combined metaprogram metrics and genes table visualization
#'
#' @param geneNMF.metaprograms GeneNMF metaprogram object
#' @param mp.genes List of gene sets per metaprogram
#' @param plots_dir Output directory
#' @param n_genes Number of top genes to show per MP
#' @param log Logging function
create_mp_combined_metrics_genes_table <- function(geneNMF.metaprograms, mp.genes, plots_dir, n_genes = 10, log = cat) {
  tryCatch({
    library(gridExtra)
    library(grid)
    library(dplyr)
    
    # Get metrics
    metrics <- geneNMF.metaprograms$metaprograms.metrics
    
    # Prepare combined data frame
    mp_names <- names(mp.genes)
    combined_data <- lapply(mp_names, function(mp_name) {
      genes <- head(mp.genes[[mp_name]], n_genes)
      metric_row <- metrics[mp_name, ]
      
      data.frame(
        MetaProgram = mp_name,
        N_Genes = length(mp.genes[[mp_name]]),
        Coverage = round(metric_row$sampleCoverage, 3),
        Silhouette = round(metric_row$silhouette, 3),
        MeanSim = round(metric_row$meanSimilarity, 3),
        Top_Genes = paste(genes, collapse = ", "),
        stringsAsFactors = FALSE
      )
    })
    
    combined_df <- do.call(rbind, combined_data)
    rownames(combined_df) <- NULL
    
    # Create table plot with larger dimensions
    png(file.path(plots_dir, "metaprogram_combined_info_table.png"),
        width = 18, height = max(10, nrow(combined_df) * 0.5), units = "in", res = 150)
    
    par(mar = c(2, 2, 4, 2))
    
    # Create table theme with CFG sizes
    tt <- ttheme_default(
      core = list(
        fg_params = list(fontsize = CFG_BASE_SIZE),
        bg_params = list(fill = "white", col = "gray80", lwd = 1.5)
      ),
      colhead = list(
        fg_params = list(fontsize = CFG_AXIS_TITLE_SIZE, fontface = "bold", col = "white"),
        bg_params = list(fill = "#2C3E50", col = "gray80", lwd = 1.5)
      )
    )
    
    g <- tableGrob(combined_df, rows = NULL, theme = tt)
    
    # Adjust column widths - more space for genes column
    g$widths <- unit(c(1.2, 0.8, 0.8, 0.8, 0.8, 6), "in")
    
    title <- textGrob("Meta-Programs: Metrics & Top Genes", 
                     gp = gpar(fontsize = CFG_TITLE_SIZE, fontface = "bold"))
    
    subtitle <- textGrob(
      sprintf("Coverage = Robust across samples | Silhouette = Internal consistency | MeanSim = Average similarity | Top %d genes per program", n_genes),
      gp = gpar(fontsize = max(10, CFG_BASE_SIZE - 2), col = "gray40")
    )
    
    grid.arrange(title, subtitle, g, 
                heights = unit(c(0.7, 0.5, nrow(combined_df) * 0.45), "in"),
                nrow = 3)
    
    dev.off()
    
    log("Saved combined metaprogram metrics and genes table visualization")
    
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    log("Warning: Combined MP table visualization failed:", e$message)
  })
}

#' Create metaprogram metrics table visualization (DEPRECATED - use create_mp_combined_metrics_genes_table)
#'
#' @param geneNMF.metaprograms GeneNMF metaprogram object
#' @param plots_dir Output directory
#' @param log Logging function
create_mp_metrics_table_plot <- function(geneNMF.metaprograms, plots_dir, log = cat) {
  log("Note: create_mp_metrics_table_plot is deprecated. Using combined metrics & genes table instead.")
  # Function deprecated - combined version is preferred
}

#' Create top genes table visualization per metaprogram (DEPRECATED - use create_mp_combined_metrics_genes_table)
#'
#' @param mp.genes List of gene sets per metaprogram
#' @param plots_dir Output directory
#' @param n_genes Number of top genes to show per MP
#' @param log Logging function
create_mp_genes_table_plot <- function(mp.genes, plots_dir, n_genes = 10, log = cat) {
  log("Note: create_mp_genes_table_plot is deprecated. Using combined metrics & genes table instead.")
  # Function deprecated - combined version is preferred
}

#' Create GSEA results visualization (dotplot style)
#'
#' @param gsea_results List of GSEA results per metaprogram
#' @param plots_dir Output directory
#' @param n_top Number of top pathways to show per MP
#' @param log Logging function
create_gsea_dotplot <- function(gsea_results, plots_dir, n_top = 5, log = cat) {
  if (is.null(gsea_results) || length(gsea_results) == 0) {
    log("No GSEA results available for visualization")
    return()
  }
  
  tryCatch({
    library(ggplot2)
    library(dplyr)
    
    # Combine all GSEA results
    all_gsea <- bind_rows(gsea_results)
    
    # Filter to C8 results only and get top pathways per MP
    if ("category" %in% colnames(all_gsea)) {
      gsea_c8 <- all_gsea %>% 
        filter(category == "C8_CellType") %>%
        group_by(MP) %>%
        slice_min(pval, n = n_top) %>%
        ungroup()
    } else {
      gsea_c8 <- all_gsea %>%
        group_by(MP) %>%
        slice_min(pval, n = n_top) %>%
        ungroup()
    }
    
    if (nrow(gsea_c8) == 0) {
      log("No C8 GSEA results to visualize")
      return()
    }
    
    # Simplify pathway names
    gsea_c8$pathway_short <- gsub("^.*?_", "", gsea_c8$pathway)
    gsea_c8$pathway_short <- gsub("_", " ", gsea_c8$pathway_short)
    gsea_c8$pathway_short <- substr(gsea_c8$pathway_short, 1, 50)
    
    # Add -log10(pval) for size
    gsea_c8$neglog10p <- -log10(gsea_c8$pval + 1e-300)
    
    p <- ggplot(gsea_c8, aes(x = MP, y = pathway_short)) +
      geom_point(aes(size = neglog10p, color = neglog10p)) +
      scale_color_gradient(low = "#FFF5E6", high = "#FF6B35", name = "-log10(p)") +
      scale_size_continuous(range = c(3, 10), name = "-log10(p)") +
      labs(
        title = "GSEA: Cell Type Signatures (C8) per Meta-Program",
        x = "Meta-Program",
        y = "Pathway"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "right",
        panel.grid.minor = element_blank()
      )
    
    ggsave(file.path(plots_dir, "GSEA_C8_dotplot.png"),
           p, width = max(12, length(unique(gsea_c8$MP)) * 1.2), 
           height = max(10, nrow(gsea_c8) * 0.3), 
           dpi = 150, limitsize = FALSE)
    
    log("Saved GSEA dotplot visualization")
    
    # Also create a table with top hit per MP
    top_hits <- gsea_c8 %>%
      group_by(MP) %>%
      slice_min(pval, n = 1) %>%
      ungroup() %>%
      select(MP, pathway_short, pval, NES = ES) %>%
      mutate(pval = sprintf("%.2e", pval))
    
    png(file.path(plots_dir, "GSEA_top_hits_table.png"),
        width = 14, height = max(6, nrow(top_hits) * 0.5), units = "in", res = 150)
    
    tt <- ttheme_default(
      core = list(
        fg_params = list(fontsize = 12),
        bg_params = list(fill = "white", col = "gray80", lwd = 1.5)
      ),
      colhead = list(
        fg_params = list(fontsize = 14, fontface = "bold"),
        bg_params = list(fill = "#4A90E2", col = "gray80", lwd = 1.5)
      )
    )
    
    g <- tableGrob(top_hits, rows = NULL, theme = tt)
    title <- textGrob("GSEA: Top Cell Type Signature per Meta-Program", 
                     gp = gpar(fontsize = 18, fontface = "bold"))
    
    grid.arrange(title, g, 
                heights = unit(c(0.8, nrow(top_hits) * 0.4), "in"),
                nrow = 2)
    
    dev.off()
    
    log("Saved GSEA top hits table")
    
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    log("Warning: GSEA visualization failed:", e$message)
  })
}
#' Create MP-based integrated UMAP
#'
#' @param seu Seurat object with MP scores
#' @param mp.genes List of meta-program gene sets
#' @param condition_field Condition metadata field
#' @param plots_dir Output directory
#' @param log Logging function
#' @return Updated Seurat object with MP-based UMAP
create_mp_integrated_umap <- function(seu, mp.genes, condition_field, plots_dir, log = cat) {
  tryCatch({
    matrix <- seu@meta.data[, names(mp.genes)]
    dimred <- as.matrix(matrix)
    colnames(dimred) <- paste0("MP_", seq(1, ncol(dimred)))
    
    seu@reductions[["MPsignatures"]] <- new("DimReduc",
                                            cell.embeddings = dimred,
                                            assay.used = "RNA",
                                            key = "MP_",
                                            global = FALSE)
    
    set.seed(123)
    seu <- RunUMAP(seu, 
                  reduction = "MPsignatures", 
                  dims = 1:ncol(dimred),
                  metric = "euclidean", 
                  reduction.name = "umap_MP")
    
    p_umap_cond <- DimPlot(seu, 
                          reduction = "umap_MP", 
                          group.by = condition_field) +
      theme(aspect.ratio = 1) +
      ggtitle("UMAP based on Meta-Program Scores")
    
    ggsave(file.path(plots_dir, "UMAP_MP_by_condition.png"), 
           p_umap_cond, width = 8, height = 7, dpi = CFG_PLOT_DPI)
    
    p_umap_celltype <- DimPlot(seu, 
                               reduction = "umap_MP", 
                               group.by = "cell_type") +
      theme(aspect.ratio = 1) +
      ggtitle("UMAP based on Meta-Program Scores")
    
    ggsave(file.path(plots_dir, "UMAP_MP_by_celltype.png"), 
           p_umap_celltype, width = 10, height = 7, dpi = CFG_PLOT_DPI)
    
    log("Saved MP-based UMAP plots")
  }, error = function(e) {
    log("Warning: MP-based UMAP failed:", e$message)
  })
  
  return(seu)
}
#' Create table visualization of known signature comparison
#'
#' @param comparison_df Data frame from compare_with_known_signatures
#' @param plots_dir Output directory
#' @param n_top Number of top matches to show
#' @param log Logging function
create_known_signatures_table <- function(comparison_df, plots_dir, n_top = 20, log = cat) {
  tryCatch({
    if (!requireNamespace("gridExtra", quietly = TRUE) || 
        !requireNamespace("grid", quietly = TRUE)) {
      log("gridExtra or grid package not available - skipping signature table plot")
      return(invisible(NULL))
    }
    
    library(gridExtra)
    library(grid)
    
    # Filter to significant and top matches
    sig_matches <- comparison_df %>%
      filter(Fisher_padj < 0.05, n_Overlap >= 2) %>%
      arrange(Fisher_pvalue) %>%
      head(n_top)
    
    if (nrow(sig_matches) == 0) {
      log("No significant matches to plot in table")
      return(invisible(NULL))
    }
    
    # Format table data
    table_data <- sig_matches %>%
      select(MetaProgram, KnownSignature, n_Overlap, n_Signature_genes, 
             Jaccard_Index, Fisher_pvalue) %>%
      mutate(
        Overlap = paste0(n_Overlap, "/", n_Signature_genes),
        Jaccard = sprintf("%.3f", Jaccard_Index),
        P_value = sprintf("%.2e", Fisher_pvalue)
      ) %>%
      select(MetaProgram, KnownSignature, Overlap, Jaccard, P_value)
    
    colnames(table_data) <- c("Meta-Program", "Known Signature", "Overlap", 
                              "Jaccard", "P-value")
    
    # Create plot
    png(file.path(plots_dir, "MP_KnownSignatures_table.png"),
        width = 14, height = max(6, nrow(table_data) * 0.4 + 1), units = "in", res = 150)
    
    # Color rows by cell type
    cell_types <- c("NK", "T_", "B_", "Macrophage", "Monocyte", "DC", "Neutrophil", 
                   "Fibroblast", "Endothelial", "Granulocyte")
    colors <- c("#E74C3C", "#3498DB", "#2ECC71", "#F39C12", "#9B59B6", 
               "#1ABC9C", "#E67E22", "#95A5A6", "#34495E", "#D35400")
    
    row_colors <- rep("white", nrow(table_data))
    for (i in 1:length(cell_types)) {
      matches <- grepl(cell_types[i], table_data$`Known Signature`)
      row_colors[matches] <- colors[i]
    }
    
    tt <- ttheme_default(
      core = list(
        fg_params = list(fontsize = CFG_BASE_SIZE),
        bg_params = list(fill = rep(c(row_colors, "gray95"), length.out = nrow(table_data)), 
                        col = "gray70", lwd = 1)
      ),
      colhead = list(
        fg_params = list(fontsize = CFG_AXIS_TITLE_SIZE, fontface = "bold", col = "white"),
        bg_params = list(fill = "#2C3E50", col = "gray70", lwd = 1.5)
      )
    )
    
    g <- tableGrob(table_data, rows = NULL, theme = tt)
    title <- textGrob("Top Meta-Program Matches with Known Cell-Type Signatures", 
                     gp = gpar(fontsize = CFG_TITLE_SIZE, fontface = "bold"))
    subtitle <- textGrob("(Significant matches: Fisher FDR < 0.05)", 
                        gp = gpar(fontsize = max(10, CFG_BASE_SIZE - 2), col = "gray30"))
    
    grid.arrange(title, subtitle, g, 
                heights = unit(c(0.6, 0.3, nrow(table_data) * 0.35), "in"),
                nrow = 3)
    
    dev.off()
    
    log("Saved known signatures comparison table")
    
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    log("Warning: Known signatures table visualization failed:", e$message)
  })
  
  invisible(NULL)
}

#' Create detailed gene overlap table for top signature matches
#'
#' @param comparison_df Data frame from compare_with_known_signatures
#' @param plots_dir Output directory
#' @param n_top Number of top matches to show
#' @param log Logging function
create_signature_overlap_genes_table <- function(comparison_df, plots_dir, n_top = 10, log = cat) {
  tryCatch({
    if (!requireNamespace("gridExtra", quietly = TRUE) || 
        !requireNamespace("grid", quietly = TRUE)) {
      log("gridExtra or grid package not available - skipping overlap genes table")
      return(invisible(NULL))
    }
    
    library(gridExtra)
    library(grid)
    
    # Get top matches with actual gene overlaps
    top_matches <- comparison_df %>%
      filter(Fisher_padj < 0.05, n_Overlap >= 3, nchar(Overlap_Genes) > 0) %>%
      arrange(Fisher_pvalue) %>%
      head(n_top)
    
    if (nrow(top_matches) == 0) {
      log("No matches with sufficient overlap to plot")
      return(invisible(NULL))
    }
    
    # Format for display - wrap long gene lists
    table_data <- top_matches %>%
      select(MetaProgram, KnownSignature, n_Overlap, Overlap_Genes) %>%
      mutate(
        Genes = sapply(Overlap_Genes, function(x) {
          genes <- strsplit(x, ", ")[[1]]
          if (length(genes) > 8) {
            paste0(paste(genes[1:8], collapse = ", "), "...")
          } else {
            x
          }
        })
      ) %>%
      select(MetaProgram, KnownSignature, n_Overlap, Genes)
    
    colnames(table_data) <- c("Meta-Program", "Known Signature", "N", "Overlapping Genes")
    
    # Create plot
    png(file.path(plots_dir, "MP_KnownSignatures_genes_table.png"),
        width = 16, height = max(6, nrow(table_data) * 0.6 + 1), units = "in", res = 150)
    
    tt <- ttheme_default(
      core = list(
        fg_params = list(fontsize = 10, hjust = 0, x = 0.05),
        bg_params = list(fill = c("white", "gray95"), col = "gray70", lwd = 1)
      ),
      colhead = list(
        fg_params = list(fontsize = 12, fontface = "bold", col = "white"),
        bg_params = list(fill = "#34495E", col = "gray70", lwd = 1.5)
      )
    )
    
    g <- tableGrob(table_data, rows = NULL, theme = tt)
    
    # Adjust column widths
    g$widths <- unit(c(0.15, 0.2, 0.08, 0.57), "npc")
    
    title <- textGrob("Overlapping Genes Between Meta-Programs and Known Signatures", 
                     gp = gpar(fontsize = 18, fontface = "bold"))
    
    grid.arrange(title, g, 
                heights = unit(c(0.6, nrow(table_data) * 0.5), "in"),
                nrow = 2)
    
    dev.off()
    
    log("Saved signature overlap genes table")
    
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    log("Warning: Overlap genes table visualization failed:", e$message)
  })
  
  invisible(NULL)
}

################################################################################
# META-PROGRAM VISUALIZATION
################################################################################

#' Generate heatmaps showing meta-program similarity and gene weights
#'
#' @param geneNMF.metaprograms GeneNMF metaprograms object
#' @param plots_dir Output directory for plots
visualize_metaprograms <- function(geneNMF.metaprograms, plots_dir) {
  log("Creating meta-program visualizations...")
  
  # 1. Similarity heatmap with larger fonts
  tryCatch({
    library(pheatmap)
    
    # Extract similarity matrix from GeneNMF object
    similarity_matrix <- geneNMF.metaprograms$metaprograms.similarity
    
    # Create annotation colors
    anno_colors <- brewer.pal(n = min(10, length(geneNMF.metaprograms$metaprograms.genes)), 
                              name = "Paired")
    names(anno_colors) <- names(geneNMF.metaprograms$metaprograms.genes)
    
    # Create annotation data frame
    annotation_df <- data.frame(
      MetaProgram = names(geneNMF.metaprograms$metaprograms.genes),
      row.names = names(geneNMF.metaprograms$metaprograms.genes)
    )
    
    annotation_colors_list <- list(MetaProgram = anno_colors)
    
    png(file.path(plots_dir, "metaprogram_similarity.png"), 
        width = 14, height = 12, units = "in", res = 150)
    
    pheatmap(similarity_matrix,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             color = viridis(100, option = "H", direction = 1),
             main = "Meta-Program Similarity",
             fontsize = 14,
             fontsize_row = 14,
             fontsize_col = 14,
             annotation_row = annotation_df,
             annotation_col = annotation_df,
             annotation_colors = annotation_colors_list,
             legend = TRUE,
             annotation_legend = TRUE,
             treeheight_row = 50,
             treeheight_col = 50,
             breaks = seq(0.1, 1, length.out = 100))
    
    dev.off()
    
    log("Saved meta-program similarity heatmap")
  }, error = function(e) {
    log("Warning: Meta-program similarity plot failed:", e$message)
  })
  
  # 2. Program-Gene Heatmap (Figure 2 style)
  tryCatch({
    mp_genes <- geneNMF.metaprograms$metaprograms.genes
    mp_weights <- geneNMF.metaprograms$metaprograms.genes.weights
    
    # Create matrix: rows = genes, columns = meta-programs
    all_genes <- unique(unlist(mp_genes))
    weight_matrix <- matrix(0, nrow = length(all_genes), ncol = length(mp_genes))
    rownames(weight_matrix) <- all_genes
    colnames(weight_matrix) <- names(mp_genes)
    
    for (mp_name in names(mp_genes)) {
      genes <- mp_genes[[mp_name]]
      weights <- mp_weights[[mp_name]][genes]
      weight_matrix[genes, mp_name] <- weights
    }
    
    # Plot heatmap
    png(file.path(plots_dir, "metaprogram_gene_weights.png"),
        width = 14, height = 12, units = "in", res = 150)
    
    library(pheatmap)
    pheatmap(weight_matrix,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             scale = "none",
             color = colorRampPalette(c("white", "orange", "red", "darkred"))(100),
             main = "Meta-Program Gene Weights",
             fontsize = 14,
             fontsize_row = 10,
             fontsize_col = 14,
             show_rownames = nrow(weight_matrix) < 200,
             legend = TRUE,
             legend_breaks = NA,
             legend_labels = NA,
             treeheight_row = 50,
             treeheight_col = 50)
    
    dev.off()
    log("Saved program-gene weight heatmap")
  }, error = function(e) {
    log("Warning: Gene weight heatmap failed:", e$message)
  })
}

# Plot Jaccard similarity heatmap between meta-programs and known signatures
plot_mp_signatures_jaccard <- function(comparison_df, plots_dir) {
  tryCatch({
    library(pheatmap)
    library(tidyr)
    library(dplyr)
    
    # Create matrix for heatmap
    jaccard_matrix <- comparison_df %>%
      select(MetaProgram, KnownSignature, Jaccard_Index) %>%
      pivot_wider(names_from = KnownSignature, values_from = Jaccard_Index) %>%
      column_to_rownames("MetaProgram") %>%
      as.matrix()
    
    png(file.path(plots_dir, "MP_vs_KnownSignatures_Jaccard.png"),
        width = 14, height = 10, units = "in", res = 150)
    
    pheatmap(jaccard_matrix,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
             main = "Meta-Programs vs. Known Cell-Type Signatures\n(Jaccard Similarity)",
             fontsize = 14,
             fontsize_row = 14,
             fontsize_col = 12,
             angle_col = 45,
             breaks = seq(0, max(jaccard_matrix), length.out = 101),
             legend = TRUE,
             treeheight_row = 50,
             treeheight_col = 50)
    
    dev.off()
    log("Saved Jaccard similarity heatmap")
  }, error = function(e) {
    log("Warning: Could not create Jaccard heatmap:", e$message)
  })
}
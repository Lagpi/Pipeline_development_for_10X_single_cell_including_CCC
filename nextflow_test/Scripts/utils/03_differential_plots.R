################################################################################
# Plotting Functions for Module 03 - Differential Analysis
# 
# All plotting functions for Differential Abundance (MiloR) and
# Differential State (Seurat + Muscat) analyses
################################################################################

# Global plotting constants (will be loaded from CFG_* in main script)
CFG_BASE_SIZE <- 18
CFG_TITLE_SIZE <- 24
CFG_AXIS_TITLE_SIZE <- 20
CFG_LEGEND_SIZE <- 16

#-------------------------------------------------------------------------------
# create_milo_purity_plot
# Neighborhood purity histogram with quality threshold line
#-------------------------------------------------------------------------------
create_milo_purity_plot <- function(purity_scores, mean_purity, plots_dir, 
                                   base_size = CFG_BASE_SIZE) {
  p <- ggplot(data.frame(purity = purity_scores), aes(x = purity)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = c(0.8, mean_purity), 
               linetype = c("dashed", "solid"), 
               color = c("red", "darkgreen"), 
               linewidth = 1) +
    theme_minimal(base_size = base_size) +
    labs(title = "MiloR: Neighborhood Purity Distribution",
         subtitle = paste0("Mean Purity: ", round(mean_purity, 3)),
         x = "Purity Score", y = "Count")
  
  ggsave(file.path(plots_dir, "milor_neighborhood_purity.png"), 
         p, width = 10, height = 7)
  
  return(p)
}

#-------------------------------------------------------------------------------
# create_milo_coverage_plot
# Cell type coverage barplot showing which cell types are covered
#-------------------------------------------------------------------------------
create_milo_coverage_plot <- function(coverage_df, plots_dir, 
                                     base_size = CFG_BASE_SIZE) {
  p_coverage <- ggplot(coverage_df, aes(x = reorder(cell_type, status), fill = status)) +
    geom_bar(stat = "count") +
    scale_fill_manual(values = c("Covered" = "steelblue", "Missing" = "coral")) +
    coord_flip() +
    theme_minimal(base_size = base_size) +
    labs(title = "MiloR: Cell Type Coverage in Neighborhoods",
         subtitle = paste0("Covered: ", sum(coverage_df$status == "Covered"), 
                          " | Missing: ", sum(coverage_df$status == "Missing")),
         x = "Cell Type", y = "Count", fill = "Status")
  
  ggsave(file.path(plots_dir, "milor_celltype_coverage.png"), 
         p_coverage, width = 12, height = 8)
  
  return(p_coverage)
}

#-------------------------------------------------------------------------------
# create_nhood_size_histograms
# Neighborhood size distributions across conditions
#-------------------------------------------------------------------------------
create_nhood_size_histograms <- function(nhood_sample_df, condition_colors, 
                                        plots_dir, base_size = CFG_BASE_SIZE) {
  plots <- list()
  
  # Overlay histogram
  p_hist_dodge <- ggplot(nhood_sample_df, aes(x = size, fill = condition)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    scale_fill_manual(values = condition_colors, name = "Condition") +
    theme_bw(base_size = base_size) +
    labs(title = "Neighborhood Size Distribution (Overlay)",
         x = "Neighborhood Size", y = "Count")
  
  ggsave(file.path(plots_dir, "DA_milo_nhood_size_hist_overlay.png"), 
         p_hist_dodge, width = 12, height = 8)
  plots$overlay <- p_hist_dodge
  
  # Density plot
  p_density <- ggplot(nhood_sample_df, aes(x = size, fill = condition, color = condition)) +
    geom_density(alpha = 0.4, size = 1) +
    scale_fill_manual(values = condition_colors, name = "Condition") +
    scale_color_manual(values = condition_colors, name = "Condition") +
    theme_bw(base_size = base_size) +
    labs(title = "Neighborhood Size Density",
         x = "Neighborhood Size", y = "Density")
  
  ggsave(file.path(plots_dir, "DA_milo_nhood_size_density.png"), 
         p_density, width = 12, height = 8)
  plots$density <- p_density
  
  # Faceted histogram
  p_hist_facet <- ggplot(nhood_sample_df, aes(x = size, fill = condition)) +
    geom_histogram(bins = 50, alpha = 0.8) +
    scale_fill_manual(values = condition_colors, name = "Condition") +
    facet_wrap(~ condition, ncol = 2, scales = "free_y") +
    theme_bw(base_size = base_size) +
    labs(title = "Neighborhood Size Distribution by Condition",
         x = "Neighborhood Size", y = "Count") +
    theme(strip.text = element_text(face = "bold", size = base_size))
  
  ggsave(file.path(plots_dir, "DA_milo_nhood_size_hist_facet.png"), 
         p_hist_facet, width = 12, height = 10)
  plots$facet <- p_hist_facet
  
  # Boxplot + violin
  p_box <- ggplot(nhood_sample_df, aes(x = condition, y = size, fill = condition)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1) +
    geom_violin(alpha = 0.3) +
    scale_fill_manual(values = condition_colors) +
    theme_bw(base_size = base_size) +
    labs(title = "Neighborhood Size Distribution",
         x = "Condition", y = "Neighborhood Size") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(plots_dir, "DA_milo_nhood_size_boxplot.png"), 
         p_box, width = 10, height = 8)
  plots$boxplot <- p_box
  
  return(plots)
}

#-------------------------------------------------------------------------------
# create_nhood_stacked_histogram
# Stacked histogram by dominant condition
#-------------------------------------------------------------------------------
create_nhood_stacked_histogram <- function(nhood_df, dominant_condition_colors, 
                                          plots_dir, base_size = CFG_BASE_SIZE) {
  p_hist <- ggplot(nhood_df, aes(x = size, fill = dominant_condition)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "stack") +
    scale_fill_manual(values = dominant_condition_colors, name = "Dominant\nCondition") +
    theme_bw(base_size = base_size) +
    labs(title = "Neighborhood Size by Dominant Condition",
         x = "Neighborhood Size", y = "Count") +
    theme(legend.position = "right")
  
  ggsave(file.path(plots_dir, "DA_milo_nhood_size_stacked.png"), 
         p_hist, width = 12, height = 8)
  
  return(p_hist)
}

#-------------------------------------------------------------------------------
# create_celltype_umap_overlay
# UMAP with cell type overlay for MiloR analysis
#-------------------------------------------------------------------------------
create_celltype_umap_overlay <- function(umap_coords, cell_types, plots_dir,
                                        base_size = CFG_BASE_SIZE) {
  cell_df <- data.frame(
    UMAP_1 = umap_coords[, 1],
    UMAP_2 = umap_coords[, 2],
    celltype = cell_types
  )
  
  p <- ggplot(cell_df, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
    geom_point(alpha = 0.3, size = 0.5) +
    theme_bw(base_size = base_size) +
    ggtitle("MiloR: Cell Types in UMAP Space") +
    labs(color = "Cell Type") +
    theme(legend.position = "right")
  
  ggsave(file.path(plots_dir, "DA_milo_celltype_umap.png"), 
         p, width = 12, height = 8)
  
  return(p)
}

#-------------------------------------------------------------------------------
# create_milo_volcano_plot
# Volcano plot for MiloR differential abundance results
#-------------------------------------------------------------------------------
create_milo_volcano_plot <- function(da_df, milo_fdr, plots_dir, dash_dir = NULL,
                                    base_size = CFG_BASE_SIZE,
                                    save_interactive = FALSE) {
  da_df$significant <- da_df$FDR < milo_fdr
  
  p <- ggplot(da_df, aes(x = logFC, y = -log10(FDR), color = significant)) +
    geom_point(alpha = 0.7, size = 2) +
    theme_bw(base_size = base_size) +
    ggtitle("MiloR: Differential Abundance Volcano Plot") +
    xlab("Log2 Fold Change") +
    ylab("-Log10(FDR)") +
    scale_color_manual(values = c("grey60", "red"), 
                      labels = c("Not Significant", "Significant")) +
    geom_hline(yintercept = -log10(milo_fdr), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30")
  
  ggsave(file.path(plots_dir, "DA_milo_volcano.png"), p, width = 10, height = 8)
  
  # Interactive version
  if (save_interactive && !is.null(dash_dir) && 
      "plotly" %in% rownames(installed.packages())) {
    p_plotly <- plotly::ggplotly(p)
    widget <- plotly::as_widget(p_plotly)
    htmlwidgets::saveWidget(widget, file.path(dash_dir, "DA_milo_volcano.html"), 
                           selfcontained = FALSE)
  }
  
  return(p)
}

#-------------------------------------------------------------------------------
# create_celltype_proportion_plots
# Stacked and side-by-side barplots of cell type proportions
#-------------------------------------------------------------------------------
create_celltype_proportion_plots <- function(props_df, condition_colors, 
                                            plots_dir, 
                                            base_size = CFG_BASE_SIZE) {
  plots <- list()
  
  # Stacked barplot
  p_stacked <- ggplot(props_df, aes(x = condition, y = proportion, fill = celltype)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(palette = "Set3") +
    theme_bw(base_size = base_size) +
    labs(title = "Cell Type Proportions by Condition (Stacked)",
         x = "Condition", y = "Proportion", fill = "Cell Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")
  
  ggsave(file.path(plots_dir, "celltype_proportions_stacked.png"), 
         p_stacked, width = 12, height = 8)
  plots$stacked <- p_stacked
  
  # Side-by-side barplot
  p_dodge <- ggplot(props_df, aes(x = celltype, y = proportion, fill = condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = condition_colors) +
    theme_bw(base_size = base_size) +
    labs(title = "Cell Type Proportions by Condition (Side-by-Side)",
         x = "Cell Type", y = "Proportion", fill = "Condition") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
  
  ggsave(file.path(plots_dir, "celltype_proportions_sidebyside.png"), 
         p_dodge, width = 14, height = 8)
  plots$sidebyside <- p_dodge
  
  return(plots)
}

#-------------------------------------------------------------------------------
# create_seurat_volcano_plot
# Volcano plot for Seurat differential state analysis
#-------------------------------------------------------------------------------
create_seurat_volcano_plot <- function(seurat_result, plot_fdr, plot_logfc,
                                      plots_dir, title = "Seurat DS: Volcano Plot",
                                      filename = "DS_seurat_volcano_all.png",
                                      base_size = CFG_BASE_SIZE,
                                      title_size = CFG_TITLE_SIZE,
                                      legend_size = CFG_LEGEND_SIZE) {
  
  # Add labels for top genes
  seurat_result$label <- ""
  seurat_result$regulation <- "Not Significant"
  sig_seurat <- seurat_result[seurat_result$significant, ]
  
  if (nrow(sig_seurat) > 0) {
    sig_seurat$score <- abs(sig_seurat$avg_log2FC) * -log10(sig_seurat$p_val_adj)
    sig_seurat <- sig_seurat[order(-sig_seurat$score), ]
    top_seurat <- head(sig_seurat, 12)
    seurat_result$label <- ifelse(seurat_result$gene %in% top_seurat$gene, 
                                  seurat_result$gene, "")
    seurat_result$regulation <- ifelse(!seurat_result$significant, "Not Significant",
                                      ifelse(seurat_result$avg_log2FC > 0, 
                                            "Up-regulated", "Down-regulated"))
  }
  
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("ggrepel package required for labeled volcano plots")
  }
  library(ggrepel)
  
  p <- ggplot(seurat_result, aes(x = avg_log2FC, y = -log10(p_val_adj), 
                                 color = regulation, label = label)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_text_repel(data = subset(seurat_result, label != ""), 
                   aes(label = label), 
                   size = 3.5, 
                   max.overlaps = 15,
                   box.padding = 0.4,
                   point.padding = 0.2,
                   segment.color = "grey50",
                   segment.size = 0.3,
                   min.segment.length = 0,
                   color = "black",
                   fontface = "bold") +
    theme_bw(base_size = base_size) +
    ggtitle(title) +
    xlab("Log2 Fold Change") +
    ylab("-Log10(Adjusted P-value)") +
    scale_color_manual(values = c("Up-regulated" = "red3", 
                                  "Down-regulated" = "blue3", 
                                  "Not Significant" = "grey60"), 
                      name = "Regulation") +
    geom_hline(yintercept = -log10(plot_fdr), linetype = "dashed", color = "darkgreen") +
    geom_vline(xintercept = c(-plot_logfc, plot_logfc), 
              linetype = "dashed", color = "grey30")
  
  ggsave(file.path(plots_dir, filename), p, width = 10, height = 8)
  
  return(p)
}

#-------------------------------------------------------------------------------
# create_muscat_volcano_plot
# Volcano plot for muscat differential state analysis
#-------------------------------------------------------------------------------
create_muscat_volcano_plot <- function(tbl, muscat_fdr, muscat_logfc, 
                                      cluster, comparison, plots_dir,
                                      base_size = CFG_BASE_SIZE,
                                      title_size = CFG_TITLE_SIZE,
                                      axis_title_size = CFG_AXIS_TITLE_SIZE,
                                      legend_size = CFG_LEGEND_SIZE) {
  
  tbl$significant <- tbl$FDR < muscat_fdr & abs(tbl$logFC) > muscat_logfc
  
  # Mark top genes
  tbl$label <- ""
  tbl$regulation <- "Not Significant"
  
  if (sum(tbl$significant, na.rm = TRUE) > 0) {
    sig_tbl <- tbl[tbl$significant, ]
    sig_tbl$score <- abs(sig_tbl$logFC) * -log10(sig_tbl$FDR)
    sig_tbl <- sig_tbl[order(-sig_tbl$score), ]
    top_n <- min(8, nrow(sig_tbl))
    
    if (top_n > 0) {
      top_genes_cluster <- head(sig_tbl, top_n)
      tbl$label <- ifelse(tbl$gene %in% top_genes_cluster$gene, tbl$gene, "")
    }
    
    tbl$regulation <- ifelse(!tbl$significant, "Not Significant",
                            ifelse(tbl$logFC > 0, "Up-regulated", "Down-regulated"))
  }
  
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("ggrepel package required for labeled volcano plots")
  }
  library(ggrepel)
  
  p <- ggplot(tbl, aes(x = logFC, y = -log10(FDR), color = regulation, label = label)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_text_repel(data = subset(tbl, label != ""), 
                   aes(label = label), 
                   size = 3.5, 
                   max.overlaps = 15,
                   box.padding = 0.4,
                   point.padding = 0.2,
                   segment.color = "grey50",
                   segment.size = 0.3,
                   min.segment.length = 0,
                   color = "black",
                   fontface = "bold") +
    theme_bw(base_size = base_size) +
    theme(axis.title = element_text(size = axis_title_size),
          plot.title = element_text(size = title_size, face = "bold"),
          legend.text = element_text(size = legend_size),
          legend.position = "bottom") +
    ggtitle(paste0(cluster, "\n(", comparison, ")")) +
    xlab("Log2 Fold Change") +
    ylab("-Log10(FDR)") +
    scale_color_manual(values = c("Up-regulated" = "red3", 
                                  "Down-regulated" = "blue3", 
                                  "Not Significant" = "grey60"), 
                      name = "Regulation") +
    geom_hline(yintercept = -log10(muscat_fdr), linetype = "dashed", color = "darkgreen") +
    geom_vline(xintercept = c(-muscat_logfc, muscat_logfc), 
              linetype = "dashed", color = "grey30")
  
  safe_cluster <- gsub("[^[:alnum:]_-]", "_", cluster)
  safe_comp <- gsub("[^[:alnum:]_-]", "_", comparison)
  ggsave(file.path(plots_dir, paste0("DS_muscat_volcano_", safe_cluster, "_", safe_comp, ".png")), 
        p, width = 10, height = 8)
  
  return(p)
}

#-------------------------------------------------------------------------------
# create_beeswarm_plot
# Beeswarm plot for gene expression across conditions
#-------------------------------------------------------------------------------
create_beeswarm_plot <- function(plot_data, cell_type, plots_dir, dash_dir = NULL,
                                base_size = CFG_BASE_SIZE,
                                title_size = CFG_TITLE_SIZE,
                                axis_title_size = CFG_AXIS_TITLE_SIZE,
                                legend_size = CFG_LEGEND_SIZE,
                                save_interactive = FALSE) {
  
  if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
    stop("ggbeeswarm package required for beeswarm plots")
  }
  library(ggbeeswarm)
  
  p <- ggplot(plot_data, aes(x = Gene, y = Expression, color = Condition)) +
    geom_quasirandom(alpha = 0.5, size = 1, dodge.width = 0.8) +
    stat_summary(fun = median, geom = "crossbar", 
                width = 0.5, color = "black", 
                position = position_dodge(width = 0.8)) +
    theme_bw(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = legend_size),
          axis.text.y = element_text(size = legend_size),
          axis.title = element_text(size = axis_title_size),
          plot.title = element_text(size = title_size, face = "bold"),
          legend.position = "top",
          legend.text = element_text(size = legend_size)) +
    ggtitle(paste("Top Gene Expression:", cell_type)) +
    ylab("Expression Level") +
    xlab("Gene") +
    scale_color_brewer(palette = "Set2")
  
  safe_ct_name <- gsub("[^[:alnum:]_-]", "_", cell_type)
  ggsave(file.path(plots_dir, paste0("DS_beeswarm_", safe_ct_name, ".png")),
         p, width = 12, height = 8)
  
  # Interactive version
  if (save_interactive && !is.null(dash_dir) && 
      "plotly" %in% rownames(installed.packages())) {
    p_interactive <- plotly::ggplotly(p)
    widget <- plotly::as_widget(p_interactive)
    htmlwidgets::saveWidget(widget, 
                           file.path(dash_dir, paste0("DS_beeswarm_", safe_ct_name, ".html")),
                           selfcontained = FALSE)
  }
  
  return(p)
}

#-------------------------------------------------------------------------------
# create_mean_ci_dotplot
# Mean ± 95% CI dot/point plot per condition and gene
# Optionally uses sample-level pseudo-bulk means (reduces pseudo-replication)
#-------------------------------------------------------------------------------
create_mean_ci_dotplot <- function(plot_data, cell_type, plots_dir, dash_dir = NULL,
                                  base_size = CFG_BASE_SIZE,
                                  title_size = CFG_TITLE_SIZE,
                                  axis_title_size = CFG_AXIS_TITLE_SIZE,
                                  legend_size = CFG_LEGEND_SIZE,
                                  save_interactive = FALSE,
                                  aggregate_by = c("auto", "sample", "cell")) {
  aggregate_by <- match.arg(aggregate_by)
  has_sample <- "Sample" %in% colnames(plot_data) && length(unique(plot_data$Sample)) > 1
  use_sample <- if (aggregate_by == "auto") has_sample else aggregate_by == "sample"
  if (use_sample && !has_sample) use_sample <- FALSE

  plot_data <- plot_data[complete.cases(plot_data[, c("Expression", "Condition", "Gene")]), ]
  if (nrow(plot_data) == 0) return(NULL)

  if (use_sample) {
    sample_means <- plot_data %>%
      dplyr::group_by(Gene, Condition, Sample) %>%
      dplyr::summarise(mean_expr = mean(Expression, na.rm = TRUE), .groups = "drop")

    summary_df <- sample_means %>%
      dplyr::group_by(Gene, Condition) %>%
      dplyr::summarise(
        mean_expr = mean(mean_expr, na.rm = TRUE),
        sd_expr = sd(mean_expr, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"
      )
  } else {
    sample_means <- NULL
    summary_df <- plot_data %>%
      dplyr::group_by(Gene, Condition) %>%
      dplyr::summarise(
        mean_expr = mean(Expression, na.rm = TRUE),
        sd_expr = sd(Expression, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"
      )
  }

  summary_df$se <- summary_df$sd_expr / sqrt(summary_df$n)
  summary_df$ci <- 1.96 * summary_df$se

  p <- ggplot(summary_df, aes(x = Condition, y = mean_expr, color = Condition)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = mean_expr - ci, ymax = mean_expr + ci),
                  width = 0.2, linewidth = 0.6) +
    theme_bw(base_size = base_size) +
    ggtitle(paste0("Mean ± 95% CI: ", cell_type,
                  if (use_sample) " (Sample-level)" else " (Cell-level)")) +
    xlab("Condition") +
    ylab("Mean Expression") +
    facet_wrap(~Gene, scales = "free_y") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = legend_size),
      axis.text.y = element_text(size = legend_size),
      axis.title = element_text(size = axis_title_size),
      plot.title = element_text(size = title_size, face = "bold"),
      legend.position = "none"
    )

  if (!is.null(sample_means)) {
    p <- p +
      geom_point(
        data = sample_means,
        aes(x = Condition, y = mean_expr),
        color = "grey40",
        size = 1.5,
        position = position_jitter(width = 0.15, height = 0)
      )
  }

  safe_ct_name <- gsub("[^[:alnum:]_-]", "_", cell_type)
  ggsave(file.path(plots_dir, paste0("DS_mean_ci_", safe_ct_name, ".png")),
         p, width = 12, height = 8)

  if (save_interactive && !is.null(dash_dir) &&
      "plotly" %in% rownames(installed.packages())) {
    p_interactive <- plotly::ggplotly(p)
    widget <- plotly::as_widget(p_interactive)
    htmlwidgets::saveWidget(widget,
                           file.path(dash_dir, paste0("DS_mean_ci_", safe_ct_name, ".html")),
                           selfcontained = FALSE)
  }

  return(p)
}

#-------------------------------------------------------------------------------
# create_combined_plot_grid
# Create combined grid of multiple plots using patchwork
#-------------------------------------------------------------------------------
create_combined_plot_grid <- function(plot_list, ncols, title, plots_dir, 
                                     filename, width_per_plot = 8, 
                                     height_per_plot = 7,
                                     title_size = CFG_TITLE_SIZE) {
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork package required for combined plots")
  }
  library(patchwork)
  
  n_plots <- length(plot_list)
  ncols <- min(ncols, n_plots)
  nrows <- ceiling(n_plots / ncols)
  
  p_combined <- wrap_plots(plot_list, ncol = ncols) +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(size = title_size + 4, 
                                              face = "bold", hjust = 0.5))
    )
  
  ggsave(file.path(plots_dir, filename),
         p_combined,
         width = ncols * width_per_plot,
         height = nrows * height_per_plot)
  
  return(p_combined)
}

#-------------------------------------------------------------------------------
# create_muscat_overall_volcano
# Overall volcano plot combining all muscat results
#-------------------------------------------------------------------------------
create_muscat_overall_volcano <- function(combined_data, muscat_fdr, muscat_logfc,
                                         plots_dir, 
                                         base_size = CFG_BASE_SIZE) {
  
  combined_data$significant <- combined_data$FDR < muscat_fdr & 
                               abs(combined_data$logFC) > muscat_logfc
  
  p <- ggplot(combined_data, aes(x = logFC, y = -log10(FDR), color = significant)) +
    geom_point(alpha = 0.4, size = 1.5) +
    theme_bw(base_size = base_size) +
    ggtitle(paste0("muscat DS: Volcano Plot (All Clusters, FDR = ", muscat_fdr, ")")) +
    xlab("Log2 Fold Change") +
    ylab("-Log10(FDR)") +
    scale_color_manual(values = c("grey60", "red"), 
                      labels = c("Not Significant", "Significant")) +
    geom_hline(yintercept = -log10(muscat_fdr), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-muscat_logfc, muscat_logfc), 
              linetype = "dashed", color = "grey30")
  
  ggsave(file.path(plots_dir, "DS_muscat_volcano_all.png"), 
         p, width = 12, height = 10)
  
  return(p)
}

#-------------------------------------------------------------------------------
# save_interactive_plot
# Helper function to save plotly interactive plots
#-------------------------------------------------------------------------------
save_interactive_plot <- function(p, filename) {
  widget <- if ("plotly" %in% class(p)) {
    plotly::as_widget(p)
  } else {
    plotly::as_widget(plotly::ggplotly(p))
  }
  htmlwidgets::saveWidget(widget, filename, selfcontained = FALSE)
}

#-------------------------------------------------------------------------------
# create_comparison_volcano_per_condition
# Create combined faceted volcano plot for all comparisons
#-------------------------------------------------------------------------------
create_comparison_volcano_per_condition <- function(seurat_result, plot_fdr, 
                                                   plot_logfc, plots_dir,
                                                   base_size = CFG_BASE_SIZE) {
  
  if (!"comparison" %in% colnames(seurat_result)) {
    return(NULL)
  }
  
  # Calculate dynamic dimensions based on number of comparisons
  n_comparisons <- length(unique(seurat_result$comparison))
  ncol_facet <- ceiling(sqrt(n_comparisons))
  nrow_facet <- ceiling(n_comparisons / ncol_facet)
  
  plot_width <- max(14, ncol_facet * 5)
  plot_height <- max(10, nrow_facet * 4)
  
  # Create combined faceted plot with fixed scales for direct comparison
  p_combined <- ggplot(seurat_result, aes(x = avg_log2FC, y = -log10(p_val_adj), 
                                          color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    facet_wrap(~ comparison, ncol = ncol_facet, scales = "fixed") +
    theme_bw(base_size = base_size) +
    ggtitle("Seurat DS (2.1): Differential Expression Volcano Plots") +
    xlab("Log2 Fold Change") +
    ylab("-Log10(Adjusted P-value)") +
    scale_color_manual(values = c("grey60", "red"), 
                      labels = c("Not Significant", "Significant"),
                      name = "Significance") +
    geom_hline(yintercept = -log10(plot_fdr), linetype = "dashed", 
               color = "blue", linewidth = 0.5) +
    geom_vline(xintercept = c(-plot_logfc, plot_logfc), 
              linetype = "dashed", color = "grey30", linewidth = 0.5) +
    theme(
      strip.text = element_text(size = base_size - 2, face = "bold"),
      legend.position = "bottom",
      panel.spacing = unit(0.8, "lines")
    )
  
  # Save combined plot
  ggsave(file.path(plots_dir, "DS_seurat_volcano_all_comparisons_combined.png"), 
        p_combined, width = plot_width, height = plot_height, dpi = 150)
  
  return(list(combined = p_combined))
}

#' Create MiloR Network Plot
#' 
#' Creates basic network visualization of differential abundance neighborhoods
#' 
#' @param milo Milo object with DA results
#' @param da Data frame with differential abundance results
#' @param plots_dir Directory to save plots
#' @param base_size Base font size for plots
#' @param title_size Title font size
#' @return ggplot object
create_milo_network_plot <- function(milo, da, plots_dir, 
                                    base_size = 16, title_size = 20) {
  
  p1 <- plotNhoodGraphDA(milo, da, layout="UMAP", alpha=0.1) +
    ggtitle("MiloR: Differential Abundance Neighborhoods") +
    theme_bw()
  
  ggsave(file.path(plots_dir, "DA_milo_network.png"), p1, width=10, height=8)
  
  return(p1)
}

#' Create MiloR Network Plot with Cell Type Annotation
#' 
#' Creates network visualization colored by cell type with combined UMAP overlay
#' 
#' @param milo Milo object
#' @param da Data frame with DA results
#' @param seu Seurat object
#' @param plots_dir Directory to save plots
#' @param base_size Base font size
#' @param log Logging function
#' @return List with celltype plot and combined plot
create_milo_network_with_celltype <- function(milo, da, seu, plots_dir,
                                             base_size = 16, log = cat) {
  
  # Determine which cell type column to use
  cell_type_col <- NULL
  if ("singler_main" %in% colnames(seu@meta.data) && 
      !all(grepl("^Cluster_", seu@meta.data$singler_main))) {
    cell_type_col <- "singler_main"
  } else if ("cell_type" %in% colnames(seu@meta.data) && 
            !all(grepl("^Cluster_", seu@meta.data$cell_type))) {
    cell_type_col <- "cell_type"
  } else if ("seurat_clusters" %in% colnames(seu@meta.data)) {
    cell_type_col <- "seurat_clusters"
  }
  
  if (is.null(cell_type_col)) {
    log("No valid cell type column found - skipping cell type annotation")
    return(NULL)
  }
  
  # Get cell type for each cell
  cell_types <- seu@meta.data[[cell_type_col]]
  
  # Assign cell types to neighborhoods based on majority voting
  nhood_index <- nhoodIndex(milo)
  
  # Check if nhood_index is valid
  if (is.null(nhood_index) || length(nhood_index) == 0 || 
      !is.matrix(nhood_index) || ncol(nhood_index) == 0 || 
      all(is.na(nhood_index))) {
    log("Warning: nhoodIndex is empty or invalid - skipping cell type annotation")
    return(NULL)
  }
  
  nhood_cell_types <- sapply(1:ncol(nhood_index), function(i) {
    cells_in_nhood <- which(nhood_index[, i] > 0)
    if (length(cells_in_nhood) == 0) return(NA)
    ct_table <- table(cell_types[cells_in_nhood])
    if (length(ct_table) == 0) return(NA)
    names(which.max(ct_table))
  })
  
  # Add to DA results
  da$celltype <- nhood_cell_types[match(da$Nhood, 1:length(nhood_cell_types))]
  
  # Plot with cell type colors
  p_celltype <- plotNhoodGraphDA(milo, da, layout="UMAP", alpha=0.1) +
    ggtitle("MiloR: Differential Abundance by Cell Type") +
    theme_bw()
  
  # Add cell type as additional layer
  p_celltype_overlay <- NULL
  p_combined_ct <- NULL
  
  if ("umap" %in% names(seu@reductions)) {
    umap_coords <- Embeddings(seu, reduction = "umap")
    
    p_celltype_overlay <- create_celltype_umap_overlay(
      umap_coords = umap_coords,
      cell_types = cell_types,
      plots_dir = plots_dir,
      base_size = base_size
    )
    
    # Create combined plot
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p_combined_ct <- p_celltype_overlay + p_celltype + 
        plot_layout(guides = "collect") +
        plot_annotation(
          title = "MiloR Differential Abundance with Cell Type Annotation",
          theme = theme(plot.title = element_text(size = 16, face = "bold"))
        )
      ggsave(file.path(plots_dir, "DA_milo_network_celltype.png"), 
            p_combined_ct, width=16, height=7)
      
      log("Created neighborhood plot with cell type annotation")
    }
  }
  
  ggsave(file.path(plots_dir, "DA_milo_network_celltype_simple.png"), 
        p_celltype, width=10, height=8)
  
  return(list(
    celltype_plot = p_celltype,
    combined_plot = p_combined_ct,
    da_updated = da
  ))
}

#' Create Combined UMAP + Neighborhood DA Plots per Condition
#' 
#' Creates side-by-side UMAP and neighborhood plots for each condition
#' 
#' @param milo Milo object
#' @param da Data frame with DA results
#' @param sample_metadata Sample metadata with condition information
#' @param plots_dir Directory to save plots
#' @param base_size Base font size
#' @return List of combined plots per condition
create_milo_combined_condition_plots <- function(milo, da, sample_metadata, plots_dir,
                                                base_size = 16) {
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("patchwork not available - skipping combined condition plots")
    return(NULL)
  }
  
  library(patchwork)
  
  all_conditions <- unique(sample_metadata$condition)
  combined_plots <- list()
  
  for (cond in all_conditions) {
    p_umap <- plotUMAP(milo) +
      ggtitle(paste0("UMAP: ", cond)) +
      theme_bw(base_size = base_size)
    
    p_nhood <- plotNhoodGraphDA(milo, da, layout="UMAP", alpha=0.05) +
      ggtitle(paste0("Neighborhood DA: ", cond)) +
      theme_bw(base_size = base_size)
    
    p_combined <- p_umap + p_nhood + 
      plot_layout(guides = "collect") +
      plot_annotation(
        title = paste0("MiloR Analysis: ", cond),
        theme = theme(plot.title = element_text(size = 16, face = "bold"))
      )
    
    safe_cond_name <- gsub("[^[:alnum:]_-]", "_", cond)
    ggsave(
      file.path(plots_dir, paste0("DA_milo_combined_", safe_cond_name, ".png")),
      p_combined,
      width = 16,
      height = 7
    )
    
    combined_plots[[cond]] <- p_combined
  }
  
  return(combined_plots)
}

#' Create Cell Type Stacked Barplot
#' 
#' Creates stacked barplot showing cell type proportions by condition
#' 
#' @param prop_df Data frame with celltype, condition, proportion columns
#' @param celltype_colors Named vector of colors for cell types
#' @param n_meta Total number of cells
#' @param plots_dir Directory to save plots
#' @param base_size Base font size
#' @param title_size Title font size
#' @param axis_title_size Axis title font size
#' @param legend_size Legend font size
#' @return ggplot object
create_celltype_stacked_barplot <- function(prop_df, celltype_colors, n_meta, plots_dir,
                                           base_size = 16, title_size = 20,
                                           axis_title_size = 18, legend_size = 14) {
  
  n_celltypes <- length(unique(prop_df$celltype))
  
  p_stacked <- ggplot(prop_df, aes(x = condition, y = proportion, fill = celltype)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = base_size),
      axis.title = element_text(size = axis_title_size, face = "bold"),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = legend_size, face = "bold"),
      legend.text = element_text(size = legend_size - 2)
    ) +
    scale_fill_manual(values = celltype_colors, name = "Cell Type") +
    labs(
      title = "Cell Type Proportions by Condition",
      subtitle = paste0("Total cells: ", n_meta),
      x = "Condition",
      y = "Proportion (%)"
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100))
  
  ggsave(
    file.path(plots_dir, "celltype_proportions_stacked.png"), 
    p_stacked, 
    width = 10 + (n_celltypes * 0.2), 
    height = 8
  )
  
  ggsave(
    file.path(plots_dir, "celltype_proportions_stacked.pdf"), 
    p_stacked, 
    width = 10 + (n_celltypes * 0.2), 
    height = 8
  )
  
  return(p_stacked)
}

#' Create Cell Type Side-by-Side Barplot
#' 
#' Creates side-by-side barplot showing cell type abundance by condition
#' 
#' @param count_df Data frame with celltype, condition, count columns
#' @param plots_dir Directory to save plots
#' @param base_size Base font size
#' @param title_size Title font size
#' @param axis_title_size Axis title font size
#' @param legend_size Legend font size
#' @return ggplot object
create_celltype_sidebyside_barplot <- function(count_df, plots_dir,
                                              base_size = 16, title_size = 20,
                                              axis_title_size = 18, legend_size = 14) {
  
  n_celltypes <- length(unique(count_df$celltype))
  
  p_sidebyside <- ggplot(count_df, aes(x = celltype, y = count, fill = condition)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = base_size - 2),
      axis.title = element_text(size = axis_title_size, face = "bold"),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      legend.position = "top",
      legend.title = element_text(size = legend_size, face = "bold"),
      legend.text = element_text(size = legend_size)
    ) +
    scale_fill_brewer(palette = "Set2", name = "Condition") +
    labs(
      title = "Cell Type Abundance by Condition",
      subtitle = "Absolute cell counts per cell type",
      x = "Cell Type",
      y = "Cell Count"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  ggsave(
    file.path(plots_dir, "celltype_abundance_sidebyside.png"), 
    p_sidebyside, 
    width = 12 + (n_celltypes * 0.3), 
    height = 8
  )
  
  ggsave(
    file.path(plots_dir, "celltype_abundance_sidebyside.pdf"), 
    p_sidebyside, 
    width = 12 + (n_celltypes * 0.3), 
    height = 8
  )
  
  return(p_sidebyside)
}

#' Create Cell Type Fold Change Plot
#' 
#' Creates barplot showing log2 fold changes in cell type proportions (2 conditions)
#' 
#' @param prop_wide Wide format data frame with proportions and fold changes
#' @param cond1 Name of condition 1
#' @param cond2 Name of condition 2
#' @param plots_dir Directory to save plots
#' @param base_size Base font size
#' @param title_size Title font size
#' @param axis_title_size Axis title font size
#' @param legend_size Legend font size
#' @return List with log2FC plot and absolute change plot
create_celltype_foldchange_plot <- function(prop_wide, cond1, cond2, plots_dir,
                                           base_size = 16, title_size = 20,
                                           axis_title_size = 18, legend_size = 14) {
  
  n_celltypes <- nrow(prop_wide)
  
  # Log2 fold change plot
  p_foldchange <- ggplot(prop_wide, aes(x = celltype, y = log2FC, fill = direction)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = base_size - 2),
      axis.title = element_text(size = axis_title_size, face = "bold"),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      legend.position = "top",
      legend.title = element_text(size = legend_size, face = "bold"),
      legend.text = element_text(size = legend_size)
    ) +
    scale_fill_manual(
      values = c("Increased" = "firebrick3", "Decreased" = "steelblue3"),
      name = paste0("Change in ", cond2, " vs ", cond1)
    ) +
    labs(
      title = "Cell Type Proportion Changes",
      subtitle = paste0(cond2, " vs ", cond1, " (Log2 Fold Change)"),
      x = "Cell Type",
      y = "Log2 Fold Change"
    )
  
  ggsave(
    file.path(plots_dir, "celltype_proportion_foldchange.png"), 
    p_foldchange, 
    width = 12 + (n_celltypes * 0.3), 
    height = 8
  )
  
  ggsave(
    file.path(plots_dir, "celltype_proportion_foldchange.pdf"), 
    p_foldchange, 
    width = 12 + (n_celltypes * 0.3), 
    height = 8
  )
  
  # Absolute change plot
  p_abschange <- ggplot(prop_wide, aes(x = celltype, y = abs_change, fill = direction)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = base_size - 2),
      axis.title = element_text(size = axis_title_size, face = "bold"),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      legend.position = "top",
      legend.title = element_text(size = legend_size, face = "bold"),
      legend.text = element_text(size = legend_size)
    ) +
    scale_fill_manual(
      values = c("Increased" = "firebrick3", "Decreased" = "steelblue3"),
      name = paste0("Change in ", cond2, " vs ", cond1)
    ) +
    labs(
      title = "Cell Type Proportion Changes (Absolute)",
      subtitle = paste0(cond2, " vs ", cond1, " (Percentage Point Difference)"),
      x = "Cell Type",
      y = "Change in Proportion (percentage points)"
    )
  
  ggsave(
    file.path(plots_dir, "celltype_proportion_absolute_change.png"), 
    p_abschange, 
    width = 12 + (n_celltypes * 0.3), 
    height = 8
  )
  
  return(list(
    log2fc_plot = p_foldchange,
    abschange_plot = p_abschange
  ))
}

#' Create Cell Type Proportion Heatmap
#' 
#' Creates heatmap showing cell type proportions across multiple conditions
#' 
#' @param prop_matrix Matrix with cell types as rows, conditions as columns
#' @param plots_dir Directory to save plots
#' @param log Logging function
#' @return TRUE if successful, FALSE otherwise
create_celltype_proportion_heatmap <- function(prop_matrix, plots_dir, log = cat) {
  
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    log("pheatmap not available - skipping heatmap")
    return(FALSE)
  }
  
  n_celltypes <- nrow(prop_matrix)
  
  png(
    file.path(plots_dir, "celltype_proportion_heatmap.png"),
    width = 10,
    height = 8 + (n_celltypes * 0.2),
    units = "in",
    res = 150
  )
  
  pheatmap::pheatmap(
    prop_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    treeheight_row = 0,
    treeheight_col = 0,
    scale = "none",
    color = colorRampPalette(c("white", "yellow", "red"))(100),
    main = "Cell Type Proportions Across Conditions (%)",
    fontsize = CFG_BASE_SIZE,
    fontsize_row = CFG_BASE_SIZE,
    fontsize_col = CFG_BASE_SIZE,
    fontsize_number = CFG_BASE_SIZE - 2,
    display_numbers = TRUE,
    number_format = "%.1f"
  )
  
  dev.off()
  
  log("Created proportion heatmap for multiple conditions")
  
  return(TRUE)
}

#' Create Muscat pbMDS Plot
#' 
#' Creates MDS plot for pseudobulk samples
#' 
#' @param pb Pseudobulk SingleCellExperiment object
#' @param plots_dir Directory to save plots
#' @return ggplot object
create_muscat_pbmds_plot <- function(pb, plots_dir) {
  
  p_pbmds <- pbMDS(pb)
  ggsave(file.path(plots_dir, "DS_muscat_pbMDS.png"), p_pbmds, width = 10, height = 8)
  
  return(p_pbmds)
}

#' Create Muscat pbHeatmap Plots
#' 
#' Creates multiple versions of pbHeatmap as per muscat vignette
#' 
#' @param pb Pseudobulk object
#' @param sce SingleCellExperiment object
#' @param res Muscat results object
#' @param plots_dir Directory to save plots
#' @param log Logging function
#' @return List with success status
create_muscat_pbheatmaps <- function(pb, sce, res, plots_dir, log = cat) {
  
  # Check if pb has the required metadata
  if (is.null(pb) || !"cluster_id" %in% names(colData(pb)) || !"sample_id" %in% names(colData(pb))) {
    log("Warning: pbHeatmap skipped - missing required metadata columns in pb object")
    return(list(success = FALSE))
  }
  
  # 1. Top-5 DS genes per cluster (default)
  png(file.path(plots_dir, "DS_muscat_pbHeatmap_top5.png"), width = 12, height = 10, units = "in", res = 150)
  tryCatch({
    pbHeatmap(sce, res, top_n = 5)
  }, error = function(e) {
    log(paste("Warning: pbHeatmap top-5 failed:", conditionMessage(e)))
  })
  dev.off()
  log("Created pbHeatmap (top-5 per cluster)")
  
  # 2. Cluster-specific heatmaps (top-20 for each cluster)
  if (!is.null(res$table)) {
    for (comparison in names(res$table)) {
      for (cluster in names(res$table[[comparison]])) {
        png(file.path(plots_dir, paste0("DS_muscat_pbHeatmap_", 
            gsub("[^[:alnum:]_-]", "_", cluster), ".png")), 
            width = 10, height = 8, units = "in", res = 150)
        tryCatch({
          pbHeatmap(sce, res, k = cluster, top_n = 20)
        }, error = function(e) {
          log(paste("Warning: pbHeatmap for cluster", cluster, "failed:", conditionMessage(e)))
        })
        dev.off()
      }
    }
    log("Created cluster-specific pbHeatmaps")
  }
  
  # 3. Gene-specific heatmap (single gene across all clusters)
  top_gene <- NULL
  if (!is.null(res$table)) {
    for (comparison in names(res$table)) {
      for (cluster in names(res$table[[comparison]])) {
        tbl <- res$table[[comparison]][[cluster]]
        if (!is.null(tbl) && is.data.frame(tbl) && nrow(tbl) > 0 && "FDR" %in% colnames(tbl)) {
          top_gene <- rownames(tbl)[which.min(tbl$FDR)]
          break
        }
      }
      if (!is.null(top_gene)) break
    }
  }
  
  if (!is.null(top_gene)) {
    png(file.path(plots_dir, paste0("DS_muscat_pbHeatmap_gene_", top_gene, ".png")), 
        width = 12, height = 8, units = "in", res = 150)
    tryCatch({
      pbHeatmap(sce, res, g = top_gene)
    }, error = function(e) {
      log(paste("Warning: pbHeatmap for gene", top_gene, "failed:", conditionMessage(e)))
    })
    dev.off()
    log(paste("Created gene-specific pbHeatmap for", top_gene))
  }
  
  return(list(success = TRUE, top_gene = top_gene))
}

#' Create Muscat UpSet Plot for Cluster Concordance
#' 
#' Creates UpSet plot showing DE gene overlap between clusters
#' 
#' @param res Muscat results object
#' @param muscat_fdr FDR threshold
#' @param muscat_logfc Log fold change threshold
#' @param plots_dir Directory to save plots
#' @param log Logging function
#' @return Number of clusters with DE genes
create_muscat_upset_plot <- function(res, muscat_fdr, muscat_logfc, plots_dir, log = cat) {
  
  if (!requireNamespace("UpSetR", quietly = TRUE) || is.null(res$table)) {
    return(0)
  }
  
  log("Creating UpSet plot for cluster concordance...")
  
  # Extract DE genes per cluster
  de_gs_by_k <- list()
  for (comparison in names(res$table)) {
    for (cluster in names(res$table[[comparison]])) {
      tbl <- res$table[[comparison]][[cluster]]
      if (!is.null(tbl) && is.data.frame(tbl) && nrow(tbl) > 0 && 
          "FDR" %in% colnames(tbl) && "logFC" %in% colnames(tbl)) {
        sig_genes <- rownames(tbl)[tbl$FDR < muscat_fdr & abs(tbl$logFC) > muscat_logfc]
        if (length(sig_genes) > 0) {
          de_gs_by_k[[cluster]] <- sig_genes
        }
      }
    }
  }
  
  if (length(de_gs_by_k) > 1) {
    png(file.path(plots_dir, "DS_muscat_upset_cluster_concordance.png"), 
        width = 14, height = 8, units = "in", res = 150)
    UpSetR::upset(UpSetR::fromList(de_gs_by_k), nsets = min(20, length(de_gs_by_k)))
    dev.off()
    log(paste("Created UpSet plot showing", length(de_gs_by_k), "clusters"))
  }
  
  return(length(de_gs_by_k))
}

#' Create Muscat UMAP/tSNE Feature Plots for Top Genes
#' 
#' Creates dimensionality reduction plots colored by expression of top DE genes
#' 
#' @param top_genes_viz Vector of gene names to plot
#' @param seu Seurat object
#' @param plots_dir Directory to save plots
#' @param base_size Base font size
#' @param title_size Title font size
#' @param log Logging function
#' @return Combined plot object
create_muscat_feature_plots <- function(top_genes_viz, seu, plots_dir, 
                                       base_size = 16, title_size = 20, log = cat) {
  
  if (length(top_genes_viz) < 4 || !("umap" %in% names(seu@reductions) || "tsne" %in% names(seu@reductions))) {
    return(NULL)
  }
  
  log("Creating dimensionality reduction plots colored by top DE gene expression...")
  
  reduction_to_use <- if ("umap" %in% names(seu@reductions)) "umap" else "tsne"
  
  # Get top 8 genes across all clusters
  top8_genes <- head(top_genes_viz, 8)
  top8_genes <- top8_genes[top8_genes %in% rownames(seu)]
  
  if (length(top8_genes) == 0) return(NULL)
  
  # Sample cells for faster plotting (max 100 per cluster)
  if (ncol(seu) > 800) {
    cells_per_cluster <- 100
    cluster_field <- if ("data_driven_celltype" %in% colnames(seu@meta.data)) 
                     "data_driven_celltype" else "seurat_clusters"
    
    # Sample from each cluster - use pmin to handle small clusters
    sampled_cells <- seu@meta.data %>%
      rownames_to_column("cell_id") %>%
      group_by(!!sym(cluster_field)) %>%
      {
        # For each group, sample the minimum of cells_per_cluster or group size
        do(., {
          sample_size <- min(cells_per_cluster, nrow(.))
          slice_sample(., n = sample_size)
        })
      } %>%
      pull(cell_id)
    
    seu_sampled <- seu[, sampled_cells]
  } else {
    seu_sampled <- seu
  }
  
  # Create feature plots
  plot_list <- lapply(top8_genes, function(g) {
    FeaturePlot(seu_sampled, features = g, reduction = reduction_to_use) +
      ggtitle(g) +
      theme_minimal(base_size = base_size) +
      theme(legend.position = "none",
            plot.title = element_text(size = title_size, face = "bold"))
  })
  
  # Combine plots
  combined_plot <- wrap_plots(plot_list, ncol = 4)
  ggsave(file.path(plots_dir, paste0("DS_muscat_", reduction_to_use, "_top_genes.png")), 
         combined_plot, width = 16, height = 8)
  log(paste("Created", reduction_to_use, "plots for top", length(top8_genes), "genes"))
  
  return(combined_plot)
}

#' Create Muscat Violin Plots per Cluster
#' 
#' Creates violin plots showing cell-level expression for top genes per cluster
#' 
#' @param res Muscat results object
#' @param sce SingleCellExperiment object
#' @param plots_dir Directory to save plots
#' @param legend_size Legend font size
#' @param log Logging function
#' @return Number of plots created
create_muscat_violin_plots <- function(res, sce, plots_dir, legend_size = 14, log = cat) {
  
  if (!requireNamespace("scater", quietly = TRUE) || is.null(sce) || is.null(res$table)) {
    return(0)
  }
  
  log("Creating violin plots for top DE genes per cluster...")
  
  n_plots <- 0
  
  for (comparison in names(res$table)) {
    for (cluster in names(res$table[[comparison]])) {
      tbl <- res$table[[comparison]][[cluster]]
      
      if (!is.null(tbl) && is.data.frame(tbl) && nrow(tbl) > 0 && 
          "FDR" %in% colnames(tbl)) {
        # Get top 6 genes for this cluster
        top6 <- rownames(tbl)[order(tbl$FDR)][1:min(6, nrow(tbl))]
        top6 <- top6[!is.na(top6)]
        
        if (length(top6) > 0) {
          # Subset SCE to cells in this cluster
          cluster_cells <- sce$cluster_id == cluster
          if (sum(cluster_cells) >= 10) {
            sce_cluster <- sce[, cluster_cells]
            
            # Create violin plot
            png(file.path(plots_dir, paste0("DS_muscat_violins_", 
                gsub("[^[:alnum:]_-]", "_", cluster), ".png")), 
                width = 12, height = 8, units = "in", res = 150)
            
            tryCatch({
              p <- scater::plotExpression(sce_cluster, 
                                         features = top6,
                                         x = "sample_id", 
                                         colour_by = "group_id",
                                         ncol = 3) +
                guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = legend_size),
                      axis.text.y = element_text(size = legend_size))
              print(p)
              n_plots <- n_plots + 1
            }, error = function(e) {
              log(paste("Warning: Violin plot for cluster", cluster, "failed:", conditionMessage(e)))
            })
            
            dev.off()
          }
        }
      }
    }
  }
  
  if (n_plots > 0) {
    log(paste("Created", n_plots, "violin plots per cluster"))
  }
  
  return(n_plots)
}

#-------------------------------------------------------------------------------
# create_muscat_beeswarm_per_celltype
# Create beeswarm plots for top genes per cell type
# Handles data preparation, plotting, and combined figure
#-------------------------------------------------------------------------------
create_muscat_beeswarm_per_celltype <- function(top_genes_viz, seu, plots_dir, dash_dir,
                                               base_size = CFG_BASE_SIZE,
                                               title_size = CFG_TITLE_SIZE,
                                               axis_title_size = CFG_AXIS_TITLE_SIZE,
                                               legend_size = CFG_LEGEND_SIZE,
                                               log = cat) {
  
  if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
    log("Warning: ggbeeswarm package not available - skipping beeswarm plots")
    return(NULL)
  }
  
  if (length(top_genes_viz) == 0) {
    log("Warning: No genes available for beeswarm plots")
    return(NULL)
  }
  
  library(ggbeeswarm)
  
  # Get cell type column
  cell_type_col <- if ("data_driven_celltype" %in% colnames(seu@meta.data)) {
    "data_driven_celltype"
  } else if ("cell.type" %in% colnames(seu@meta.data)) {
    "cell.type"
  } else {
    NULL
  }
  
  if (is.null(cell_type_col)) {
    log("Warning: No cell type column found - skipping beeswarm plots")
    return(NULL)
  }
  
  cell_types <- unique(seu@meta.data[[cell_type_col]])
  beeswarm_plots <- list()
  mean_ci_plots <- list()

  sample_col <- if ("sample" %in% colnames(seu@meta.data)) {
    "sample"
  } else if ("orig.ident" %in% colnames(seu@meta.data)) {
    "orig.ident"
  } else {
    NULL
  }
  
  # Create beeswarm plot for each cell type
  for (ct in cell_types) {
    # Subset to this cell type
    seu_subset <- subset(seu, subset = !!sym(cell_type_col) == ct)
    
    if (ncol(seu_subset) < 10) {
      log(paste("Skipping cell type", ct, "- too few cells:", ncol(seu_subset)))
      next
    }
    
    # Prepare data for plotting
    plot_data_list <- list()
    for (gene in top_genes_viz) {
      if (gene %in% rownames(seu_subset)) {
        vars_to_fetch <- c(gene, "condition")
        if (!is.null(sample_col)) {
          vars_to_fetch <- c(vars_to_fetch, sample_col)
        }
        expr <- FetchData(seu_subset, vars = vars_to_fetch)
        colnames(expr) <- c("Expression", "Condition", if (!is.null(sample_col)) "Sample")
        expr$Gene <- gene
        plot_data_list[[gene]] <- expr
      }
    }
    
    if (length(plot_data_list) > 0) {
      plot_data <- do.call(rbind, plot_data_list)
      
      # Create beeswarm plot using existing function
      p_beeswarm <- create_beeswarm_plot(
        plot_data = plot_data,
        cell_type = ct,
        plots_dir = plots_dir,
        dash_dir = dash_dir,
        base_size = base_size,
        title_size = title_size,
        axis_title_size = axis_title_size,
        legend_size = legend_size,
        save_interactive = TRUE
      )
      
      # Store for combined plot
      beeswarm_plots[[ct]] <- p_beeswarm
      
      log(paste("Created beeswarm plot for cell type:", ct))

      # Mean ± 95% CI dot plot (sample-level if available)
      p_mean_ci <- create_mean_ci_dotplot(
        plot_data = plot_data,
        cell_type = ct,
        plots_dir = plots_dir,
        dash_dir = dash_dir,
        base_size = base_size,
        title_size = title_size,
        axis_title_size = axis_title_size,
        legend_size = legend_size,
        save_interactive = TRUE,
        aggregate_by = "auto"
      )

      if (!is.null(p_mean_ci)) {
        mean_ci_plots[[ct]] <- p_mean_ci
        log(paste("Created mean ± CI dot plot for cell type:", ct))
      }
    }
  }
  
  # Create combined beeswarm figure
  if (length(beeswarm_plots) > 0) {
    p_combined_beeswarm <- create_combined_plot_grid(
      plot_list = beeswarm_plots,
      ncols = 3,
      title = "Gene Expression Across All Cell Types",
      plots_dir = plots_dir,
      filename = "DS_beeswarm_all_celltypes_combined.png",
      width_per_plot = 10,
      height_per_plot = 8,
      title_size = title_size
    )
    
    log(paste("Created combined beeswarm plot with", length(beeswarm_plots), "cell types"))
  }

  if (length(mean_ci_plots) > 0) {
    p_combined_mean_ci <- create_combined_plot_grid(
      plot_list = mean_ci_plots,
      ncols = 3,
      title = "Mean ± 95% CI by Condition (Top Genes)",
      plots_dir = plots_dir,
      filename = "DS_mean_ci_all_celltypes_combined.png",
      width_per_plot = 10,
      height_per_plot = 8,
      title_size = title_size
    )

    log(paste("Created combined mean ± CI plot with", length(mean_ci_plots), "cell types"))
  } else {
    p_combined_mean_ci <- NULL
  }

  if (length(beeswarm_plots) > 0 || length(mean_ci_plots) > 0) {
    return(list(
      individual_plots = beeswarm_plots,
      combined_plot = if (exists("p_combined_beeswarm")) p_combined_beeswarm else NULL,
      mean_ci_plots = mean_ci_plots,
      combined_mean_ci = p_combined_mean_ci
    ))
  }
  
  return(NULL)
}

#-------------------------------------------------------------------------------
# create_validation_funnel_plot
# Validation funnel showing Seurat discovery -> Muscat validation
#-------------------------------------------------------------------------------
create_validation_funnel_plot <- function(validation_summary, n_candidates, n_validated, 
                                         validation_rate, plots_dir,
                                         base_size = CFG_BASE_SIZE,
                                         title_size = CFG_TITLE_SIZE,
                                         axis_title_size = CFG_AXIS_TITLE_SIZE) {
  
  validation_summary_plot <- validation_summary[1:4, ]
  validation_summary_plot$Category <- factor(
    validation_summary_plot$Category,
    levels = c("Seurat Candidates (Discovery)", "Tested by Muscat", 
              "Validated (Robust)", "Seurat-only (Not Validated)")
  )
  
  p_funnel <- ggplot(validation_summary_plot, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
    geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), 
             vjust = -0.5, size = 5, fontface = "bold") +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = base_size - 2),
      axis.title = element_text(size = axis_title_size, face = "bold"),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      legend.position = "none"
    ) +
    ggtitle("Differential Gene Validation Funnel") +
    labs(
      subtitle = paste0("Discovery (Seurat): ", n_candidates, " genes | ",
                       "Validated (Muscat): ", n_validated, " genes (", validation_rate, "%)")
    ) +
    ylab("Number of Genes") +
    xlab("") +
    scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
  
  ggsave(file.path(plots_dir, "DS_validation_funnel.png"), p_funnel, width = 12, height = 9)
  
  return(p_funnel)
}

#-------------------------------------------------------------------------------
# create_validation_scatter_plot
# Scatter plot comparing Seurat vs Muscat log fold changes
#-------------------------------------------------------------------------------
create_validation_scatter_plot <- function(validated_genes, n_validated, n_tested, 
                                          plots_dir, dash_dir,
                                          base_size = CFG_BASE_SIZE,
                                          title_size = CFG_TITLE_SIZE,
                                          log = cat) {
  
  if (nrow(validated_genes) == 0) return(NULL)
  
  validated_genes$validation_label <- ifelse(
    validated_genes$final_validated,
    "Validated",
    ifelse(validated_genes$muscat_significant, "Muscat sig. but discordant", "Not validated")
  )
  
  cor_val <- cor(validated_genes$avg_log2FC, validated_genes$logFC_muscat, 
                use = "complete.obs", method = "pearson")
  
  p_scatter <- ggplot(validated_genes, 
                     aes(x = avg_log2FC, y = logFC_muscat, color = validation_label)) +
    geom_point(alpha = 0.6, size = 2.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
    theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = title_size, face = "bold"),
      plot.subtitle = element_text(size = base_size - 2),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    ) +
    ggtitle("Log2 Fold Change: Seurat Discovery vs Muscat Validation") +
    labs(
      subtitle = paste0("Pearson r = ", round(cor_val, 3), " | ",
                       n_validated, "/", n_tested, " genes validated")
    ) +
    xlab("Seurat Log2FC (Discovery)") +
    ylab("Muscat Log2FC (Validation)") +
    scale_color_manual(
      name = "Validation Status",
      values = c("Validated" = "#4DAF4A", 
                "Muscat sig. but discordant" = "#FF7F00",
                "Not validated" = "#999999")
    ) +
    annotate("text", x = Inf, y = -Inf, 
            label = paste("n =", nrow(validated_genes), "genes tested"), 
            hjust = 1.1, vjust = -0.5, size = 4.5)
  
  ggsave(file.path(plots_dir, "DS_validation_logFC_scatter.png"), p_scatter, width = 11, height = 9)
  
  # Interactive version
  if (requireNamespace("plotly", quietly = TRUE)) {
    validated_genes$hover_text <- paste0(
      "Gene: ", validated_genes$gene, "\n",
      "Cluster: ", validated_genes$cluster, "\n",
      "Seurat LogFC: ", round(validated_genes$avg_log2FC, 3), "\n",
      "Muscat LogFC: ", round(validated_genes$logFC_muscat, 3), "\n",
      "Status: ", validated_genes$validation_label
    )
    
    p_interactive <- plotly::ggplotly(p_scatter, tooltip = "hover_text")
    htmlwidgets::saveWidget(p_interactive, 
                          file.path(dash_dir, "DS_validation_logFC_scatter.html"),
                          selfcontained = FALSE)
  }
  
  log(paste("Created validation scatter plot (correlation:", round(cor_val, 3), ")"))
  
  return(p_scatter)
}

#-------------------------------------------------------------------------------
# create_validation_by_cluster_plot
# Validation rates per cluster/cell type
#-------------------------------------------------------------------------------
create_validation_by_cluster_plot <- function(validated_genes, plots_dir, data_dir,
                                             base_size = CFG_BASE_SIZE,
                                             title_size = CFG_TITLE_SIZE,
                                             axis_title_size = CFG_AXIS_TITLE_SIZE) {
  
  if (!"cluster" %in% colnames(validated_genes)) return(NULL)
  
  cluster_validation <- validated_genes %>%
    group_by(cluster) %>%
    summarise(
      n_candidates = n(),
      n_validated = sum(final_validated, na.rm = TRUE),
      validation_rate = round(100 * n_validated / n_candidates, 1),
      .groups = "drop"
    ) %>%
    arrange(desc(validation_rate))
  
  write.csv(cluster_validation, file.path(data_dir, "DS_validation_by_cluster.csv"), row.names = FALSE)
  
  # Plot validation rates by cluster
  p_cluster <- ggplot(cluster_validation, aes(x = reorder(cluster, validation_rate), 
                                              y = validation_rate)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_text(aes(label = paste0(n_validated, "/", n_candidates)), 
             hjust = -0.2, size = 3.5) +
    coord_flip() +
    theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = title_size, face = "bold"),
      axis.title = element_text(size = axis_title_size)
    ) +
    ggtitle("Validation Rate by Cell Type/Cluster") +
    xlab("Cluster") +
    ylab("Validation Rate (%)") +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.1)))
  
  ggsave(file.path(plots_dir, "DS_validation_by_cluster.png"), p_cluster, 
        width = 10, height = max(6, 0.4 * nrow(cluster_validation)))
  
  return(p_cluster)
}

#-------------------------------------------------------------------------------
# create_method_comparison_upset
# UpSet plot comparing Seurat and Muscat gene detection
#-------------------------------------------------------------------------------
create_method_comparison_upset <- function(seurat_genes, muscat_genes, plots_dir, log = cat) {
  
  if (!requireNamespace("UpSetR", quietly = TRUE)) return(NULL)
  
  library(UpSetR)
  
  all_genes <- unique(c(seurat_genes, muscat_genes))
  upset_data <- data.frame(
    Gene = all_genes,
    Seurat = as.integer(all_genes %in% seurat_genes),
    muscat = as.integer(all_genes %in% muscat_genes)
  )
  
  png(file.path(plots_dir, "method_comparison_upset.png"), width = 10, height = 6, units = "in", res = 150)
  print(upset(upset_data, sets = c("Seurat", "muscat"), 
             order.by = "freq",
             main.bar.color = "steelblue",
             sets.bar.color = c("red", "blue")))
  dev.off()
  
  log("Created UpSet plot for method comparison")
  
  return(invisible(NULL))
}

#-------------------------------------------------------------------------------
# create_method_comparison_barplot
# Barplot comparing Seurat and Muscat gene counts
#-------------------------------------------------------------------------------
create_method_comparison_barplot <- function(comparison_summary, seurat_genes, muscat_genes, 
                                            overlap_genes, seurat_only, muscat_only,
                                            plots_dir,
                                            muscat_fdr = 0.05,
                                            muscat_genes_at_fdr_threshold = NULL,
                                            base_size = CFG_BASE_SIZE,
                                            title_size = CFG_TITLE_SIZE,
                                            axis_title_size = CFG_AXIS_TITLE_SIZE) {
  
  comparison_summary_plot <- comparison_summary[1:5, ]
  
  # Build muscat label with FDR breakdown
  muscat_label <- paste0("muscat\n(", length(muscat_genes), " sig. genes")
  if (!is.null(muscat_genes_at_fdr_threshold) && 
      muscat_genes_at_fdr_threshold != length(muscat_genes)) {
    # Fallback occurred: show both counts
    genes_added <- length(muscat_genes) - muscat_genes_at_fdr_threshold
    muscat_label <- paste0(
      muscat_label, "\n",
      muscat_genes_at_fdr_threshold, " @ 0.05 + ",
      genes_added, " @ 0.1)"
    )
  } else {
    muscat_label <- paste0(muscat_label, ")")
  }
  
  comparison_summary_plot$Method <- factor(
    comparison_summary_plot$Method,
    levels = c("Seurat", "muscat", "Overlap", "Seurat Only", "muscat Only"),
    labels = c(
      paste0("Seurat\n(", length(seurat_genes), " sig. genes)"),
      muscat_label,
      paste0("Overlap\n(", length(overlap_genes), " genes)"),
      paste0("Seurat Only\n(", length(seurat_only), " genes)"),
      paste0("muscat Only\n(", length(muscat_only), " genes)")
    )
  )
  
  p_barplot <- ggplot(comparison_summary_plot, aes(x = Method, y = Count, fill = Method)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = Count), vjust = -0.5, size = 6, fontface = "bold") +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = base_size - 2),
      axis.title = element_text(size = axis_title_size, face = "bold"),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      legend.position = "none"
    ) +
    ggtitle("Differential Gene Detection: Method Comparison") +
    labs(
      subtitle = paste0(
        "Seurat: ", length(seurat_genes), " significant | ",
        "muscat: ", length(muscat_genes), " significant",
        if (!is.null(muscat_genes_at_fdr_threshold) && 
            muscat_genes_at_fdr_threshold != length(muscat_genes))
          paste0(" (", muscat_genes_at_fdr_threshold, " @ FDR 0.05 + ", 
                 length(muscat_genes) - muscat_genes_at_fdr_threshold, " @ FDR 0.1)") 
        else "",
        " | Overlap: ", length(overlap_genes), " genes"
      ),
      y = "Number of Significant Genes"
    ) +
    scale_fill_manual(values = c("red3", "blue3", "purple3", "pink3", "lightblue3")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  ggsave(file.path(plots_dir, "method_comparison_barplot.png"), p_barplot, width = 12, height = 9)
  
  return(p_barplot)
}

#-------------------------------------------------------------------------------
# create_method_comparison_logfc_scatter
# Scatter plot comparing log fold changes between methods
#-------------------------------------------------------------------------------
create_method_comparison_logfc_scatter <- function(seurat_sig, overlap_genes, muscat_data,
                                                   plots_dir, data_dir,
                                                   base_size = CFG_BASE_SIZE,
                                                   title_size = CFG_TITLE_SIZE,
                                                   log = cat) {
  
  if (length(overlap_genes) == 0 || length(muscat_data) == 0) {
    return(NULL)
  }
  
  # Combine muscat data
  muscat_combined <- dplyr::bind_rows(muscat_data)
  
  # Merge by gene
  seurat_fc <- seurat_sig[seurat_sig$gene %in% overlap_genes, c("gene", "avg_log2FC")]
  colnames(seurat_fc) <- c("gene", "Seurat_logFC")
  
  muscat_fc <- muscat_combined[muscat_combined$gene %in% overlap_genes, c("gene", "logFC")]
  colnames(muscat_fc) <- c("gene", "muscat_logFC")
  
  # Aggregate muscat by gene (take mean if multiple entries)
  muscat_fc <- aggregate(muscat_logFC ~ gene, data = muscat_fc, FUN = mean)
  
  # Merge
  fc_comparison <- merge(seurat_fc, muscat_fc, by = "gene")
  
  if (nrow(fc_comparison) == 0) return(NULL)
  
  # Calculate correlation
  cor_val <- cor(fc_comparison$Seurat_logFC, fc_comparison$muscat_logFC, use = "complete.obs")
  
  p_scatter <- ggplot(fc_comparison, aes(x = Seurat_logFC, y = muscat_logFC)) +
    geom_point(alpha = 0.6, size = 2, color = "purple") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = title_size, face = "bold")
    ) +
    ggtitle(paste0("Log Fold Change Comparison (Correlation: ", round(cor_val, 3), ")")) +
    xlab("Seurat Log2FC") +
    ylab("muscat Log2FC") +
    annotate("text", x = Inf, y = -Inf, 
            label = paste("n =", nrow(fc_comparison), "genes"), 
            hjust = 1.1, vjust = -0.5, size = 5)
  
  ggsave(file.path(plots_dir, "method_comparison_logFC_scatter.png"), p_scatter, width = 10, height = 8)
  write.csv(fc_comparison, file.path(data_dir, "method_comparison_logFC.csv"), row.names = FALSE)
  
  log(paste("Created scatter plot comparing log fold changes (correlation:", round(cor_val, 3), ")"))
  
  return(p_scatter)
}

#===============================================================================
# ADDITIONAL WRAPPER AND UTILITY FUNCTIONS
# (Moved from main script to ensure modularity)
#===============================================================================

#-------------------------------------------------------------------------------
# generate_colors
# Generate extended color palette for many groups using RColorBrewer
#-------------------------------------------------------------------------------
#' Generate extended color palette for many groups
#' @param n Number of colors needed
#' @return Character vector of colors
generate_colors <- function(n) {
  if (n <= 0) return(character(0))
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    return(grDevices::rainbow(n))
  }
  
  if (n <= 12) {
    return(RColorBrewer::brewer.pal(max(3, n), "Set3")[1:n])
  }
  
  if (n <= 24) {
    return(c(
      RColorBrewer::brewer.pal(12, "Set3"),
      RColorBrewer::brewer.pal(12, "Paired")
    )[1:n])
  }
  
  # For > 24 colors, interpolate
  base_colors <- c(
    RColorBrewer::brewer.pal(12, "Set3"),
    RColorBrewer::brewer.pal(12, "Paired"),
    RColorBrewer::brewer.pal(9, "Set1")
  )
  return(grDevices::colorRampPalette(base_colors)(n))
}

#-------------------------------------------------------------------------------
# create_neighborhood_validation_plots
# Create validation plots for MiloR neighborhoods
#-------------------------------------------------------------------------------
#' Create neighborhood validation plots
#' @param validation_results Validation results from validate_milor_neighborhoods
#' @param plots_dir Output directory for plots
#' @param dash_dir Output directory for interactive plots
#' @param log Logging function
create_neighborhood_validation_plots <- function(validation_results, plots_dir, dash_dir, log = cat) {
  if (!is.null(validation_results$reason)) return()
  
  # 1. Purity histogram
  p_purity_hist <- create_milo_purity_plot(
    purity_scores = validation_results$purity_scores,
    mean_purity = validation_results$mean_purity,
    plots_dir = plots_dir,
    base_size = CFG_BASE_SIZE
  )
  save_interactive_plot(p_purity_hist, 
                       file.path(dash_dir, "milor_neighborhood_purity.html"))
  
  # 2. Cell type coverage barplot
  if (length(validation_results$covered_cell_types) > 0) {
    coverage_df <- data.frame(
      cell_type = c(validation_results$covered_cell_types, 
                    validation_results$missing_cell_types),
      status = c(rep("Covered", length(validation_results$covered_cell_types)),
                rep("Missing", length(validation_results$missing_cell_types)))
    )
    
    p_coverage <- create_milo_coverage_plot(
      coverage_df = coverage_df,
      plots_dir = plots_dir,
      base_size = CFG_BASE_SIZE
    )
    
    ggsave(file.path(plots_dir, "milor_celltype_coverage.png"), 
           p_coverage, width = 10, height = max(6, length(unique(coverage_df$cell_type)) * 0.3), 
           dpi = 150)
    save_interactive_plot(p_coverage, 
                         file.path(dash_dir, "milor_celltype_coverage.html"))
  }
  
  # 3. Validation summary text plot
  summary_text <- paste0(
    "MiloR Neighborhood Validation Summary\n\n",
    "Purity Metrics:\n",
    "  • Mean Purity: ", round(validation_results$mean_purity * 100, 1), "%\n",
    "  • Median Purity: ", round(validation_results$median_purity * 100, 1), "%\n",
    "  • High Purity Neighborhoods (≥80%): ", round(validation_results$pct_high_purity, 1), "%\n\n",
    "Coverage Metrics:\n",
    "  • Cell Type Coverage: ", round(validation_results$coverage_pct, 1), "%\n",
    "  • Smallest Cluster: ", validation_results$smallest_cluster_name, 
    " (", validation_results$smallest_cluster_size, " cells)\n\n",
    "Validation Status: ", 
    ifelse(validation_results$validation_passed, "✓ PASSED", "✗ FAILED"), "\n\n"
  )
  
  if (length(validation_results$recommendations) > 0) {
    summary_text <- paste0(summary_text, "Recommendations:\n")
    for (i in seq_along(validation_results$recommendations)) {
      summary_text <- paste0(summary_text, "  ", i, ". ", 
                           validation_results$recommendations[i], "\n")
    }
  }
  
  # Save as text file
  writeLines(summary_text, file.path(plots_dir, "milor_validation_summary.txt"))
  
  log("Validation plots saved")
}

#-------------------------------------------------------------------------------
# create_milor_plots
# Create comprehensive MiloR visualization plots
#-------------------------------------------------------------------------------
#' Create MiloR visualization plots
create_milor_plots <- function(milo, da, seu, cond_field, sample_metadata, milo_fdr, 
                              plots_dir, dash_dir, validation_results = NULL, conditions = NULL, log = cat) {
  
  # Handle both old single-da and new list-of-da formats
  if (!is.list(da) || is.data.frame(da) || is.null(names(da))) {
    # Old format: single da object
    da_list <- list(da)
    names(da_list) <- "global_model"
    contrast_pairs <- list("global_model")
  } else if (length(da) > 0 && is.character(names(da)) && names(da)[1] != "") {
    # New format: list of pairwise contrasts
    da_list <- da
    contrast_pairs <- names(da)
  } else {
    log("No DA results available - skipping MiloR plots")
    return()
  }
  
  if (is.null(da) || (is.list(da) && length(da) == 0)) {
    log("No DA results available - skipping MiloR plots")
    return()
  }
  
  log("Creating MiloR visualization plots...")
  
  # Create validation plots if available
  if (!is.null(validation_results)) {
    create_neighborhood_validation_plots(validation_results, plots_dir, dash_dir, log)
  }
  
  # ========== INTRODUCTION: Basic MiloR Diagnostics ==========
  log("Creating Introduction: Basic MiloR diagnostics...")
  
  # Plot neighborhood size histogram (miloR built-in)
  p_nhood_size <- plotNhoodSizeHist(milo)
  ggsave(file.path(plots_dir, "DA_milo_intro_nhood_size_hist.png"), p_nhood_size, width = 10, height = 8)
  save_interactive_plot(ggplotly(p_nhood_size), 
                       file.path(dash_dir, "DA_milo_intro_nhood_size_hist.html"))
  log("Saved Introduction: Neighborhood size histogram")
  
  # Plot neighborhood distribution on UMAP (if available)
  # Check if UMAP reduction exists in reducedDims
  available_reductions <- SingleCellExperiment::reducedDimNames(milo)
  
  if (length(available_reductions) > 0 && "UMAP" %in% available_reductions) {
    # Create one UMAP with cell types at the top, then all neighborhood DA plots in a grid
    log("Creating combined figure: UMAP with cell types + all neighborhood DA plots...")
    
    # Plot 1: UMAP with annotated cell types (only once at the top)
    umap_coords <- SingleCellExperiment::reducedDim(milo, "UMAP")
    
    # Find cell type annotation column
    celltype_col <- if ("data_driven_celltype" %in% colnames(seu@meta.data)) {
      "data_driven_celltype"
    } else if ("cell.type" %in% colnames(seu@meta.data)) {
      "cell.type"
    } else if ("celltype" %in% colnames(seu@meta.data)) {
      "celltype"
    } else {
      "seurat_clusters"
    }
    
    plot_df <- data.frame(
      UMAP_1 = umap_coords[, 1],
      UMAP_2 = umap_coords[, 2],
      celltype = seu@meta.data[[celltype_col]][match(rownames(umap_coords), rownames(seu@meta.data))]
    )
    
    # Generate colors for cell types
    n_celltypes <- length(unique(plot_df$celltype))
    celltype_colors <- generate_colors(n_celltypes)
    
    # Calculate centroids for each cell type for label placement
    label_df <- plot_df %>%
      group_by(celltype) %>%
      summarize(
        UMAP_1 = median(UMAP_1),
        UMAP_2 = median(UMAP_2),
        .groups = "drop"
      )
    
    p_umap_celltype <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
      geom_point(size = 1.2, alpha = 0.7) +
      geom_text(data = label_df, aes(label = celltype), 
                color = "black", size = 5, fontface = "bold",
                hjust = 0.5, vjust = 0.5) +
      scale_color_manual(values = celltype_colors, name = "Cell Type") +
      theme_bw(base_size = 18) +
      ggtitle("UMAP: Annotated Cell Types") +
      labs(x = "UMAP 1", y = "UMAP 2") +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)
      )
    
    # Save standalone UMAP
    ggsave(file.path(plots_dir, "DA_milo_umap_celltypes.png"),
           p_umap_celltype, width = 12, height = 8)
    
    # Plot 2: Create all neighborhood DA plots in a list
    all_nhood_plots <- list()
    nhood_legend <- NULL
    
    for (i in seq_along(names(da_list))) {
      contrast_pair <- names(da_list)[i]
      da_for_plot <- as.data.frame(da_list[[contrast_pair]])
      contrast_safe <- make.names(contrast_pair)

      if (!all(c("SpatialFDR", "logFC") %in% colnames(da_for_plot))) {
        log(paste("Skipping neighborhood DA plot for", contrast_pair, "- missing SpatialFDR/logFC"))
        next
      }
      
      p_nhood <- tryCatch({
        plotNhoodGraphDA(milo, da_for_plot, layout="UMAP", alpha=0.05) +
          theme_bw(base_size = 11) +
          ggtitle(contrast_pair) +
          theme(
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            legend.position = "right"
          )
      }, error = function(e) {
        log(paste("Warning: Could not create neighborhood DA plot for", contrast_pair, "-", e$message))
        NULL
      })

      if (!is.null(p_nhood)) {
        # Extract legend from first plot
        if (is.null(nhood_legend)) {
          nhood_legend <- cowplot::get_legend(p_nhood)
        }
        
        # Remove legend from all plots
        p_nhood <- p_nhood + theme(legend.position = "none")
        all_nhood_plots[[contrast_pair]] <- p_nhood
        
        # Save individual neighborhood DA plot
        ggsave(file.path(plots_dir, paste0("DA_milo_nhood_", contrast_safe, ".png")),
               p_nhood, width = 8, height = 6)
      }
    }
    
    # Create combined figure: UMAP at top + grid of neighborhood DA plots below
    if (length(all_nhood_plots) > 0) {
      if (requireNamespace("patchwork", quietly = TRUE)) {
        # Calculate grid dimensions
        n_plots <- length(all_nhood_plots)
        n_cols <- min(3, n_plots)  # Max 3 columns
        n_rows <- ceiling(n_plots / n_cols)
        
        # Create grid of neighborhood plots
        nhood_grid <- wrap_plots(all_nhood_plots, ncol = n_cols)
        
        # Create legend panel with 3 columns if available
        if (!is.null(nhood_legend)) {
          legend_plot <- cowplot::plot_grid(nhood_legend, ncol = 3)
        } else {
          legend_plot <- plot_spacer()
        }
        
        # Combine: UMAP (2 cols) + Legend (1 col) on top row, then full grid below
        final_figure <- (p_umap_celltype | legend_plot) / nhood_grid + 
          plot_layout(heights = c(1, n_rows * 0.8), widths = c(2, 1))
        
        ggsave(file.path(plots_dir, "DA_milo_all_contrasts_combined.png"),
               final_figure, width = n_cols * 8, height = 8 + n_rows * 6, limitsize = FALSE)
        
        log(paste("Saved final combined figure: UMAP (2 cols) + Legend (1 col, 3 rows) + ", length(all_nhood_plots), "neighborhood DA plots"))
      } else {
        log("Warning: patchwork package not available, saving plots separately")
      }
    }
    
    # Plot 2: Neighborhoods colored by dominant condition
    nhood_counts <- nhoodCounts(milo)
    nhood_sizes <- colSums(nhood_counts > 0)
    
    # Determine dominant condition for each neighborhood
    dominant_conditions <- sapply(1:ncol(nhood_counts), function(i) {
      cells_in_nhood <- rownames(nhood_counts)[nhood_counts[, i] > 0]
      if (length(cells_in_nhood) == 0) return(NA)
      cond_table <- table(seu@meta.data[[cond_field]][match(cells_in_nhood, rownames(seu@meta.data))])
      if (length(cond_table) == 0) return(NA)
      names(which.max(cond_table))
    })
    
    # Create UMAP plot with neighborhood colors by dominant condition
    umap_coords <- SingleCellExperiment::reducedDim(milo, "UMAP")
    condition_colors <- generate_colors(length(unique(seu@meta.data[[cond_field]])))
    names(condition_colors) <- unique(seu@meta.data[[cond_field]])
    
    plot_df <- data.frame(
      UMAP_1 = umap_coords[, 1],
      UMAP_2 = umap_coords[, 2],
      condition = seu@meta.data[[cond_field]][match(rownames(umap_coords), rownames(seu@meta.data))]
    )
    
    p_cond_dist <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = condition)) +
      geom_point(size = 2, alpha = 0.6) +
      scale_color_manual(values = condition_colors, name = "Condition") +
      theme_bw() +
      ggtitle("Cell distribution by Condition") +
      labs(x = "UMAP 1", y = "UMAP 2")
    
    ggsave(file.path(plots_dir, "DA_milo_intro_conditions_distribution.png"), p_cond_dist, width = 10, height = 8)
    save_interactive_plot(ggplotly(p_cond_dist), 
                         file.path(dash_dir, "DA_milo_intro_conditions_distribution.html"))
    log("Saved Introduction: Cell distribution by 4 Conditions on UMAP")
  } else {
    log("Warning: UMAP reduction not available in Milo object, skipping UMAP-based plots")
  }
  
  # ========== Neighborhood size histogram
    nhood_counts <- nhoodCounts(milo)
    nhood_sizes <- colSums(nhood_counts > 0)
    
    # Debug: Check sample information
    log(paste("Sample column names in nhood_counts:", toString(head(colnames(nhood_counts), 10))))
    log(paste("Sample field in metadata:", if("sample" %in% colnames(seu@meta.data)) "sample" else "orig.ident"))
    
    temp_design_df <- data.frame(
      sample = seu@meta.data[[if("sample" %in% colnames(seu@meta.data)) "sample" else "orig.ident"]],
      condition = seu@meta.data[[cond_field]]
    )
    
    log(paste("Unique samples in metadata:", toString(unique(temp_design_df$sample))))
    log(paste("Unique conditions:", toString(unique(temp_design_df$condition))))
    
    # Create proper sample-to-condition mapping
    # Use row names from nhood_counts which correspond to cells
    cell_samples <- rownames(nhood_counts)
    sample_to_condition <- setNames(
      seu@meta.data[[cond_field]][match(cell_samples, rownames(seu@meta.data))],
      cell_samples
    )
    
    log(paste("Sample-to-condition mapping size:", length(sample_to_condition)))
    log(paste("NA values in mapping:", sum(is.na(sample_to_condition))))
    
    dominant_conditions <- apply(nhood_counts, 2, function(col) {
      if(all(col == 0)) return("Unknown")
      cells_in_nhood <- names(col)[col > 0]
      if (length(cells_in_nhood) == 0) return("Unknown")
      cell_conditions <- sample_to_condition[cells_in_nhood]
      cell_conditions <- cell_conditions[!is.na(cell_conditions)]
      if (length(cell_conditions) == 0) return("Unknown")
      cond_table <- table(cell_conditions)
      names(which.max(cond_table))
    })
    
    # Create expanded data for per-condition analysis
    # For each neighborhood, create entries for each condition it contains
    nhood_sample_data <- list()
    for (i in seq_along(nhood_sizes)) {
      nhood_col <- nhood_counts[, i]
      cells_in_nhood <- names(nhood_col)[nhood_col > 0]
      
      if (length(cells_in_nhood) > 0) {
        # Get conditions for cells in this neighborhood
        conditions_in_nhood <- sample_to_condition[cells_in_nhood]
        conditions_in_nhood <- conditions_in_nhood[!is.na(conditions_in_nhood)]
        
        if (length(conditions_in_nhood) > 0) {
          # Count cells per condition in this neighborhood
          cond_counts <- table(conditions_in_nhood)
          
          for (cond_name in names(cond_counts)) {
            nhood_sample_data[[paste(i, cond_name, sep = "_")]] <- data.frame(
              nhood_id = i,
              condition = cond_name,
              size = nhood_sizes[i],
              cell_count = as.numeric(cond_counts[cond_name]),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    
    # Check if we have data
    nhood_sample_df <- if (length(nhood_sample_data) == 0) {
      log("Warning: No valid neighborhood-sample data - skipping histogram plots")
      data.frame()
    } else {
      df <- do.call(rbind, nhood_sample_data)
      df <- df[complete.cases(df), ]
      log(paste("Created neighborhood-sample data:", nrow(df), "rows"))
      df
    }
    
    # Original data with dominant condition
    nhood_df <- data.frame(
      nhood_id = seq_along(nhood_sizes),
      size = nhood_sizes,
      dominant_condition = factor(dominant_conditions)
    )
    
    log(paste("Neighborhoods by condition:", toString(table(nhood_df$dominant_condition))))
    log(paste("Unique conditions:", toString(unique(nhood_sample_df$condition))))

    if (nrow(nhood_sample_df) > 0) {
        # Create all neighborhood size histograms using plot function
        nhood_plots <- create_nhood_size_histograms(
          nhood_sample_df = nhood_sample_df,
          condition_colors = condition_colors,
          plots_dir = plots_dir,
          base_size = CFG_BASE_SIZE
        )
        
        # Original stacked version (kept for compatibility) - works with nhood_df
        if ("dominant_condition" %in% colnames(nhood_df)) {
          dominant_condition_colors <- generate_colors(length(unique(nhood_df$dominant_condition)))
          names(dominant_condition_colors) <- unique(nhood_df$dominant_condition)
          
          p_hist_stacked <- create_nhood_stacked_histogram(
            nhood_df = nhood_df,
            dominant_condition_colors = dominant_condition_colors,
            plots_dir = plots_dir,
            base_size = CFG_BASE_SIZE
          )
        }
        
        if("plotly" %in% rownames(installed.packages())) {
          # Interactive overlay version for better exploration
          p_hist_interactive <- plotly::ggplotly(nhood_plots$overlay)
          save_interactive_plot(p_hist_interactive, file.path(dash_dir, "DA_milo_nhood_size_hist.html"))
        }
        
        log("Saved neighborhood size histogram")
    }
  
  # Network plot - skip for pairwise contrasts (use only first DA result if available)
  da_for_network <- if (is.list(da) && !is.data.frame(da) && length(da) > 0) {
    log("Note: Using first contrast for network/volcano plots (pairwise contrasts detected)")
    da[[1]]  # Use first pairwise contrast
  } else {
    da  # Use single DA result
  }
  
  if (!is.null(da_for_network) && is.data.frame(da_for_network)) {
    p1 <- NULL
    log("Skipping network plots for pairwise contrasts to avoid plotting issues")
  } else {
    log("Skipping network plot - no valid DA results")
    p1 <- NULL
  }
  
  # Network plot with cell type annotation - skip for pairwise contrasts
  celltype_result <- NULL
  log("Skipping celltype network plot for pairwise contrasts")
  
  # Update da with cell type if available
  if (!is.null(celltype_result) && !is.null(celltype_result$da_updated)) {
    da_for_network <- celltype_result$da_updated
  }
  
  # Skip combined condition and volcano plots for pairwise contrasts
  combined_cond_plots <- NULL
  p_volcano <- NULL
  log("Skipping combined condition and volcano plots for pairwise contrasts")
  
  log("Completed MiloR plots")
}

#-------------------------------------------------------------------------------
# create_celltype_proportion_plots (FULL WRAPPER)
# Generates stacked barplots, side-by-side barplots, fold change plots, and heatmaps
# showing cell type distribution across conditions
#-------------------------------------------------------------------------------
create_celltype_proportion_plots <- function(seu, cond_field, celltype_field = "data_driven_celltype", 
                                            plots_dir, data_dir, log = cat) {
  
  log("=== Creating Cell Type Proportion Plots ===")
  
  # Check if required columns exist
  if (!cond_field %in% colnames(seu@meta.data)) return()
  
  # Try different celltype field names
  celltype_candidates <- c(celltype_field, "individual_celltype", "cell_type", "celltype", "seurat_clusters")
  celltype_field_actual <- intersect(celltype_candidates, colnames(seu@meta.data))[1]
  if (is.na(celltype_field_actual)) return()
  
  log(paste("Using cell type field:", celltype_field_actual))
  log(paste("Using condition field:", cond_field))
  
  # Extract metadata
  meta <- seu@meta.data[, c(cond_field, celltype_field_actual)]
  colnames(meta) <- c("condition", "celltype")
  meta <- meta[!is.na(meta$condition) & !is.na(meta$celltype), ]
  if (nrow(meta) == 0) return()
  
  conditions <- unique(meta$condition)
  celltypes <- unique(meta$celltype)
  
  log(paste("Conditions:", length(conditions), "-", paste(conditions, collapse = ", ")))
  log(paste("Cell types:", length(celltypes)))
  
  # Calculate counts and proportions
  count_table <- table(meta$celltype, meta$condition)
  count_df <- as.data.frame(count_table)
  colnames(count_df) <- c("celltype", "condition", "count")
  
  # Calculate proportions per condition
  prop_df <- count_df %>%
    group_by(condition) %>%
    mutate(
      total = sum(count),
      proportion = count / total * 100
    ) %>%
    ungroup()
  
  # Save data
  write.csv(count_df, file.path(data_dir, "celltype_counts_by_condition.csv"), row.names = FALSE)
  write.csv(prop_df, file.path(data_dir, "celltype_proportions_by_condition.csv"), row.names = FALSE)
  
  log("Saved cell type count and proportion tables")
  
  # Generate colors
  n_celltypes <- length(celltypes)
  celltype_colors <- generate_colors(n_celltypes)
  names(celltype_colors) <- sort(celltypes)
  
  # ============================================================================
  # 1. STACKED BARPLOT - Proportions
  # ============================================================================
  p_stacked <- create_celltype_stacked_barplot(
    prop_df = prop_df,
    celltype_colors = celltype_colors,
    n_meta = nrow(meta),
    plots_dir = plots_dir,
    base_size = CFG_BASE_SIZE,
    title_size = CFG_TITLE_SIZE,
    axis_title_size = CFG_AXIS_TITLE_SIZE,
    legend_size = CFG_LEGEND_SIZE
  )
  
  log("Created stacked barplot of cell type proportions")
  
  # ============================================================================
  # 2. SIDE-BY-SIDE BARPLOT - Absolute counts
  # ============================================================================
  p_sidebyside <- create_celltype_sidebyside_barplot(
    count_df = count_df,
    plots_dir = plots_dir,
    base_size = CFG_BASE_SIZE,
    title_size = CFG_TITLE_SIZE,
    axis_title_size = CFG_AXIS_TITLE_SIZE,
    legend_size = CFG_LEGEND_SIZE
  )
  
  log("Created side-by-side barplot of cell type abundance")
  
  # ============================================================================
  # 3. PROPORTION CHANGE PLOT - Fold changes
  # ============================================================================
  if (length(conditions) == 2) {
      cond1 <- conditions[1]
      cond2 <- conditions[2]
      
      # Get proportions for each condition
      prop_wide <- prop_df %>%
        select(celltype, condition, proportion) %>%
        tidyr::pivot_wider(names_from = condition, values_from = proportion, values_fill = 0)
      
      # Calculate fold change (log2)
      prop_wide$log2FC <- log2((prop_wide[[cond2]] + 0.1) / (prop_wide[[cond1]] + 0.1))
      prop_wide$abs_change <- prop_wide[[cond2]] - prop_wide[[cond1]]
      
      # Sort by absolute change
      prop_wide <- prop_wide %>% arrange(desc(abs(abs_change)))
      prop_wide$celltype <- factor(prop_wide$celltype, levels = prop_wide$celltype)
      
      # Color by direction
      prop_wide$direction <- ifelse(prop_wide$log2FC > 0, "Increased", "Decreased")
      
      # Save fold change data
      write.csv(prop_wide, file.path(data_dir, "celltype_proportion_foldchanges.csv"), row.names = FALSE)
      
      # Create plots using function
      fc_plots <- create_celltype_foldchange_plot(
        prop_wide = prop_wide,
        cond1 = cond1,
        cond2 = cond2,
        plots_dir = plots_dir,
        base_size = CFG_BASE_SIZE,
        title_size = CFG_TITLE_SIZE,
        axis_title_size = CFG_AXIS_TITLE_SIZE,
        legend_size = CFG_LEGEND_SIZE
      )
      
      log(paste("Created fold change plot:", cond2, "vs", cond1))
      log("Created absolute change plot")
      
  } else if (length(conditions) > 2) {
    log("More than 2 conditions detected - creating heatmap of proportions")
    
      # Create proportion matrix for heatmap
      prop_matrix <- prop_df %>%
        select(celltype, condition, proportion) %>%
        tidyr::pivot_wider(names_from = condition, values_from = proportion, values_fill = 0) %>%
        as.data.frame()
      
      rownames(prop_matrix) <- prop_matrix$celltype
      prop_matrix$celltype <- NULL
      prop_matrix <- as.matrix(prop_matrix)
      
      # Save matrix
      write.csv(prop_matrix, file.path(data_dir, "celltype_proportion_matrix.csv"), row.names = TRUE)
      
      # Create heatmap using function
      heatmap_result <- create_celltype_proportion_heatmap(
        prop_matrix = prop_matrix,
        plots_dir = plots_dir,
        log = log
      )
  }
  
  log("=== Cell Type Proportion Plots Completed ===")
}

#-------------------------------------------------------------------------------
# create_muscat_plots
# Generates volcano plots, heatmaps, and expression plots for pseudobulk DE results
#-------------------------------------------------------------------------------
#' Create muscat visualization plots
#' 
#' Generates volcano plots, heatmaps, and expression plots for pseudobulk DE results
#' 
#' @param muscat_result Muscat results list
#' @param seu Seurat object
#' @param plots_dir Output directory for plots
#' @param dash_dir Output directory for interactive plots
#' @param log Logging function
create_muscat_plots <- function(muscat_result, seu, plots_dir, dash_dir, log = cat) {
  
  if (is.null(muscat_result$res)) {
    log("No muscat results available - skipping plots")
    return()
  }
  
  log("Creating muscat visualization plots...")
  
  res <- muscat_result$res
  pb <- muscat_result$pb
  sce <- muscat_result$sce
  muscat_fdr <- muscat_result$muscat_fdr
  muscat_logfc <- muscat_result$muscat_logfc
  
  # pbMDS plot
  p_pbmds <- create_muscat_pbmds_plot(
    pb = pb,
    plots_dir = plots_dir
  )
  
  # pbHeatmap - multiple versions as per muscat vignette
  heatmap_result <- create_muscat_pbheatmaps(
    pb = pb,
    sce = sce,
    res = res,
    plots_dir = plots_dir,
    log = log
  )
  
  # Extract top genes for additional plots
  top_genes_viz <- character(0)
  if (!is.null(res$table)) {
    for (comparison in names(res$table)) {
      for (cluster in names(res$table[[comparison]])) {
        tbl <- res$table[[comparison]][[cluster]]
        
        # Multiple safety checks
        if (is.null(tbl)) next
        if (!is.data.frame(tbl) && !is.matrix(tbl)) next
        if (nrow(tbl) == 0) next
        if (!"FDR" %in% colnames(tbl)) next
        
        # Try to get FDR column safely
        fdr_col <- tbl[["FDR"]]
        if (is.null(fdr_col)) next
        if (!is.numeric(fdr_col) && !is.vector(fdr_col)) next
        
        fdr_vals <- as.numeric(fdr_col)
        if (length(fdr_vals) == 0) next
        if (all(is.na(fdr_vals))) next
        if (all(is.infinite(fdr_vals))) next
        
        # Get row names safely
        gene_names <- rownames(tbl)
        if (is.null(gene_names) || length(gene_names) != length(fdr_vals)) next
        
        # Order and select top genes
        valid_idx <- !is.na(fdr_vals) & !is.infinite(fdr_vals)
        if (sum(valid_idx) == 0) next
        
        valid_fdr <- fdr_vals[valid_idx]
        valid_genes <- gene_names[valid_idx]
        
        ord <- order(valid_fdr)
        ordered_genes <- valid_genes[ord]
        top_genes_viz <- c(top_genes_viz, head(ordered_genes, 2))
      }
    }
  }
  
  # Process extracted genes
  top_genes_viz <- unique(top_genes_viz)
  if (length(top_genes_viz) > 0) {
    top_genes_viz <- top_genes_viz[1:min(6, length(top_genes_viz))]
    log(paste("Extracted", length(top_genes_viz), "top genes for visualization"))
  } else {
    log("No valid genes extracted from muscat results")
  }
  
  # If still no genes, use highly variable genes as fallback
  if (length(top_genes_viz) == 0) {
      hvg <- if ("RNA" %in% names(seu@assays)) VariableFeatures(seu) else character(0)
      if (length(hvg) > 0) top_genes_viz <- head(hvg, 6)
  }
  
  # ======================================================================
  # NEW: UpSet Plot - Between-cluster concordance
  # ======================================================================
  n_clusters_upset <- create_muscat_upset_plot(
    res = res,
    muscat_fdr = muscat_fdr,
    muscat_logfc = muscat_logfc,
    plots_dir = plots_dir,
    log = log
  )
  
  # ======================================================================
  # NEW: UMAP/t-SNE colored by expression - Top DS genes
  # ======================================================================
  feature_plots <- create_muscat_feature_plots(
    top_genes_viz = top_genes_viz,
    seu = seu,
    plots_dir = plots_dir,
    base_size = CFG_BASE_SIZE,
    title_size = CFG_TITLE_SIZE,
    log = log
  )
  
  # ======================================================================
  # NEW: Violin plots - Cell-level expression for top genes per cluster
  # ======================================================================
  n_violin_plots <- create_muscat_violin_plots(
    res = res,
    sce = sce,
    plots_dir = plots_dir,
    legend_size = CFG_LEGEND_SIZE,
    log = log
  )
  
  # ======================================================================
  # Beeswarm plots for top genes per cell type
  # ======================================================================
  beeswarm_result <- create_muscat_beeswarm_per_celltype(
    top_genes_viz = top_genes_viz,
    seu = seu,
    plots_dir = plots_dir,
    dash_dir = dash_dir,
    base_size = CFG_BASE_SIZE,
    title_size = CFG_TITLE_SIZE,
    axis_title_size = CFG_AXIS_TITLE_SIZE,
    legend_size = CFG_LEGEND_SIZE,
    log = log
  )
  
  log("Completed muscat plots")
}

#-------------------------------------------------------------------------------
# perform_combined_ds_validation
# Strategy: Use Seurat for high sensitivity candidate discovery, then validate with Muscat
# Returns: Data frame with validated genes and validation metrics
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# compare_seurat_muscat
# Compare Seurat and muscat differential expression results
#-------------------------------------------------------------------------------
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

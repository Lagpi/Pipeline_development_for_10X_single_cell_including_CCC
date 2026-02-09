#===============================================================================
# 02_clustering/plots.R
# Clustering and annotation specific plotting functions
#===============================================================================

#===============================================================================
# PLOTTING CONFIGURATION CONSTANTS
#===============================================================================

# Global plotting configuration - used across all plot functions
CFG_BASE_FONT_SIZE <- 18
CFG_TITLE_FONT_SIZE <- 22
CFG_AXIS_TITLE_SIZE <- 18
CFG_AXIS_TEXT_SIZE <- 16
CFG_LEGEND_TEXT_SIZE <- 16
CFG_PLOT_WIDTH <- 12
CFG_PLOT_HEIGHT <- 6

#===============================================================================
# DUAL PLOT SAVING (PNG + INTERACTIVE HTML)
#===============================================================================

# Save plot in both static (PNG) and interactive (HTML) formats
save_dual_plot <- function(plot, filename, plots_dir, dash_dir, 
                           width = CFG_PLOT_WIDTH, height = CFG_PLOT_HEIGHT, dpi = 150) {
  
  # Check if plot is a ggMarginal object (has ggExtra class attributes)
  is_ggmarginal <- inherits(plot, "ggExtraPlot") || 
                   (!inherits(plot, "gg") && inherits(plot, "gtable"))
  
  # Save static PNG
  png_path <- file.path(plots_dir, paste0(filename, ".png"))
  if (is_ggmarginal) {
    # ggMarginal objects need special handling
    png(png_path, width = width, height = height, units = "in", res = dpi)
    grid::grid.draw(plot)
    dev.off()
  } else {
    ggsave(png_path, plot, width = width, height = height, dpi = dpi)
  }
  
  # Save interactive HTML with larger font
  # For ggMarginal plots, extract the main plot without marginals for interactivity
  if (is_ggmarginal) {
    # Extract the main plot from ggMarginal object
    main_plot <- plot$layers[[1]]$ggplot  # Try to get main ggplot
    if (is.null(main_plot) || !inherits(main_plot, "gg")) {
      # Fallback: skip interactive for marginal plots
      return(invisible(NULL))
    }
    plot_for_html <- main_plot
  } else {
    plot_for_html <- plot
  }
  
  plot_interactive <- plotly::ggplotly(plot_for_html) %>%
    plotly::layout(
      font = list(size = CFG_BASE_FONT_SIZE),
      title = list(font = list(size = CFG_TITLE_FONT_SIZE)),
      xaxis = list(titlefont = list(size = CFG_AXIS_TITLE_SIZE)),
      yaxis = list(titlefont = list(size = CFG_AXIS_TITLE_SIZE))
    )
  html_path <- file.path(dash_dir, paste0(filename, ".html"))
  htmlwidgets::saveWidget(plot_interactive, html_path, selfcontained = TRUE)
}

#===============================================================================
# PCA VISUALIZATION
#===============================================================================

# Create comprehensive PCA visualization plots
# Returns: List of ggplot objects (variance, elbow, pca_condition, pca_sample, loadings)
create_pca_plots <- function(seu, pc_dims = NULL, jackstraw_used = FALSE) {
  log("Creating PCA visualization plots...")
  
  # Configuration constants
  CFG_PCA_SEARCH_SPACE <- 50
  CFG_TOP_LOADING_GENES <- 10
  base_font <- 16
  title_font <- 18
  axis_font <- 14
  axis_text <- 13
  
  plots <- list()
  
  # Calculate variance explained
  pca_obj <- seu[["pca"]]
  variance_explained <- (pca_obj@stdev)^2 / sum((pca_obj@stdev)^2) * 100
  cumulative_variance <- cumsum(variance_explained)
  
  # Variance plot
  var_data <- data.frame(
    PC = 1:length(variance_explained),
    Variance = variance_explained,
    Cumulative = cumulative_variance
  )
  
  # Build subtitle with JackStraw information
  var_subtitle <- paste("First", CFG_PCA_SEARCH_SPACE, "PCs explain", 
                       round(cumulative_variance[min(CFG_PCA_SEARCH_SPACE, length(cumulative_variance))], 1), 
                       "% of variance")
  if (!is.null(pc_dims)) {
    if (jackstraw_used) {
      var_subtitle <- paste0(var_subtitle, " | JackStraw selected ", length(pc_dims), " significant PCs")
    } else {
      var_subtitle <- paste0(var_subtitle, " | Using ", length(pc_dims), " PCs (selected)")
    }
  }
  
  # Main variance plot
  p_variance <- ggplot(var_data[1:min(CFG_PCA_SEARCH_SPACE, nrow(var_data)), ], aes(x = PC)) +
    geom_col(aes(y = Variance), fill = "steelblue", alpha = 0.7) +
    geom_line(aes(y = Cumulative / 5), color = "red", linewidth = 1) +
    scale_y_continuous(
      name = "Variance Explained (%)",
      sec.axis = sec_axis(~. * 5, name = "Cumulative Variance (%)", 
                         breaks = seq(0, 100, 20))
    ) +
    theme_minimal(base_size = base_font) +
    labs(title = "PCA Variance Analysis",
         subtitle = var_subtitle,
         x = "Principal Component") +
    theme(axis.title.y.right = element_text(color = "red", size = axis_font),
          axis.text.y.right = element_text(color = "red", size = axis_text),
          axis.text.x = element_text(size = axis_text),
          axis.text.y = element_text(size = axis_text),
          axis.title = element_text(size = axis_font))
  
  # Distribution plot on the side
  p_distribution <- ggplot(var_data[1:min(CFG_PCA_SEARCH_SPACE, nrow(var_data)), ], aes(x = Variance)) +
    geom_density(fill = "steelblue", alpha = 0.5) +
    geom_vline(xintercept = median(var_data$Variance[1:min(CFG_PCA_SEARCH_SPACE, nrow(var_data))]), 
               linetype = "dashed", color = "darkred", linewidth = 1) +
    coord_flip() +
    theme_minimal(base_size = base_font - 2) +
    labs(y = "Density", x = "Variance (%)") +
    theme(axis.text = element_text(size = axis_text - 2))
  
  # Combine plots side by side
  if (requireNamespace("patchwork", quietly = TRUE)) {
    plots$variance <- p_variance + p_distribution + patchwork::plot_layout(widths = c(4, 1))
  } else {
    plots$variance <- p_variance
  }
  
  # Elbow plot with JackStraw information
  if (!is.null(pc_dims)) {
    if (jackstraw_used) {
      subtitle_text <- paste("JackStraw selected", length(pc_dims), "significant PCs (p<0.05), capturing", 
                            round(cumulative_variance[max(pc_dims)], 1), "% variance")
    } else {
      subtitle_text <- paste("Using", length(pc_dims), "PCs (selected), capturing", 
                            round(cumulative_variance[max(pc_dims)], 1), "% variance")
    }
  } else {
    subtitle_text <- "Variance contribution per principal component"
  }
  
  plots$elbow <- Seurat::ElbowPlot(seu, ndims = min(CFG_PCA_SEARCH_SPACE, ncol(pca_obj))) +
    ggtitle("PCA Elbow Plot") +
    labs(subtitle = subtitle_text) +
    theme_minimal(base_size = base_font) +
    theme(axis.text.x = element_text(size = axis_text),
          axis.text.y = element_text(size = axis_text),
          axis.title = element_text(size = axis_font))
  
  # PCA scatter by condition
  if ("condition" %in% colnames(seu@meta.data)) {
    p_pca_cond <- Seurat::DimPlot(seu, reduction = "pca", group.by = "condition", pt.size = 0.5) +
      ggtitle("PCA by Condition") +
      labs(subtitle = paste0("PC1 (", round(variance_explained[1], 1), "%) vs PC2 (", 
                            round(variance_explained[2], 1), "%)")) +
      theme_minimal(base_size = base_font) +
      theme(axis.text.x = element_text(size = axis_text),
            axis.text.y = element_text(size = axis_text),
            axis.title = element_text(size = axis_font))
    
    # Add marginal distribution plots
    if (requireNamespace("ggExtra", quietly = TRUE)) {
      plots$pca_condition <- ggExtra::ggMarginal(p_pca_cond, type = "density", 
                                                  groupColour = TRUE, groupFill = TRUE,
                                                  alpha = 0.3, size = 4)
    } else {
      plots$pca_condition <- p_pca_cond
    }
  }
  
  # PCA scatter by sample
  if ("sample" %in% colnames(seu@meta.data)) {
    p_pca_samp <- Seurat::DimPlot(seu, reduction = "pca", group.by = "sample", pt.size = 0.5) +
      ggtitle("PCA by Sample") +
      labs(subtitle = paste0("PC1 (", round(variance_explained[1], 1), "%) vs PC2 (", 
                            round(variance_explained[2], 1), "%)")) +
      theme_minimal(base_size = base_font) +
      theme(axis.text.x = element_text(size = axis_text),
            axis.text.y = element_text(size = axis_text),
            axis.title = element_text(size = axis_font))
    
    # Add marginal distribution plots
    if (requireNamespace("ggExtra", quietly = TRUE)) {
      plots$pca_sample <- ggExtra::ggMarginal(p_pca_samp, type = "density", 
                                               groupColour = TRUE, groupFill = TRUE,
                                               alpha = 0.3, size = 4)
    } else {
      plots$pca_sample <- p_pca_samp
    }
  }
  
  # Top loading genes
  top_loadings <- pca_obj@feature.loadings[, 1:2, drop = FALSE]
  top_genes_pc1 <- rownames(top_loadings)[order(abs(top_loadings[, 1]), decreasing = TRUE)][1:CFG_TOP_LOADING_GENES]
  top_genes_pc2 <- rownames(top_loadings)[order(abs(top_loadings[, 2]), decreasing = TRUE)][1:CFG_TOP_LOADING_GENES]
  
  load_df <- data.frame(
    Gene = c(top_genes_pc1, top_genes_pc2),
    Loading = c(top_loadings[top_genes_pc1, 1], top_loadings[top_genes_pc2, 2]),
    PC = rep(c("PC1", "PC2"), each = CFG_TOP_LOADING_GENES)
  )
  
  plots$loadings <- ggplot(load_df, aes(x = reorder(Gene, abs(Loading)), y = Loading)) +
    geom_col(fill = "steelblue") +
    facet_wrap(~PC, scales = "free") +
    coord_flip() +
    theme_minimal(base_size = base_font) +
    labs(title = "Top Contributing Genes (PC1/PC2)", x = "Gene", y = "Loading")
  
  return(plots)
}

#===============================================================================
# UMAP VISUALIZATION
#===============================================================================

# Create standardized UMAP plot with consistent theming
# Returns: ggplot object
create_umap_plot <- function(seu, group_by, title, label = FALSE, repel = FALSE, 
                             colors = NULL, reduction = "umap", add_marginal = TRUE) {
  # Create base plot
  p <- Seurat::DimPlot(seu, group.by = group_by, reduction = reduction, 
                      label = label, repel = repel)
  
  # Add custom colors if provided
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }
  
  # Add consistent theming
  p <- p + 
    labs(title = title) + 
    theme_minimal(base_size = CFG_BASE_FONT_SIZE) +
    theme(plot.title = element_text(size = CFG_TITLE_FONT_SIZE, face = "bold"),
          axis.text.x = element_text(size = CFG_AXIS_TEXT_SIZE),
          axis.text.y = element_text(size = CFG_AXIS_TEXT_SIZE),
          axis.title = element_text(size = CFG_AXIS_TITLE_SIZE))
  
  # Add marginal distribution plots if requested
  if (add_marginal && requireNamespace("ggExtra", quietly = TRUE)) {
    p <- ggExtra::ggMarginal(p, type = "density", 
                             groupColour = TRUE, groupFill = TRUE,
                             alpha = 0.3, size = 4)
  }
  
  return(p)
}

# Create UMAP plot with custom labels from a different metadata column
# This allows showing one grouping (e.g., Phase) colored, but labeled by another (e.g., cell_type)
create_umap_plot_with_labels <- function(seu, group_by, label_by, title, 
                                         colors = NULL, reduction = "umap", 
                                         label_size = 6) {
  # Get embeddings
  embedding <- seu@reductions[[reduction]]@cell.embeddings
  df <- data.frame(
    UMAP_1 = embedding[, 1],
    UMAP_2 = embedding[, 2],
    group = seu@meta.data[[group_by]],
    label = seu@meta.data[[label_by]]
  )
  
  # Calculate label positions (centroids)
  label_pos <- df %>%
    group_by(label) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      .groups = "drop"
    )
  
  # Create plot
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = group)) +
    geom_point(size = 0.5, alpha = 0.6) +
    geom_text(data = label_pos, aes(label = label), 
              size = label_size, fontface = "bold", color = "black") +
    labs(title = title, color = group_by) +
    theme_minimal(base_size = CFG_BASE_FONT_SIZE) +
    theme(
      plot.title = element_text(size = CFG_TITLE_FONT_SIZE, face = "bold"),
      axis.text.x = element_text(size = CFG_AXIS_TEXT_SIZE),
      axis.text.y = element_text(size = CFG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = CFG_AXIS_TITLE_SIZE),
      legend.text = element_text(size = CFG_LEGEND_TEXT_SIZE),
      legend.title = element_text(size = CFG_AXIS_TITLE_SIZE)
    )
  
  # Add custom colors if provided
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }
  
  return(p)
}

#===============================================================================
# FEATURE PLOT VISUALIZATION
#===============================================================================

# Create standardized FeaturePlot with consistent theming
# Returns: ggplot object
create_feature_plot <- function(seu, features, title, label = FALSE, label.size = 4,
                                repel = FALSE, group.by = NULL, reduction = "umap",
                                color_option = "plasma") {
  # Create base plot (group.by not supported in this Seurat version)
  p <- Seurat::FeaturePlot(seu, features = features, reduction = reduction,
                          label = label, label.size = label.size, 
                          repel = repel)
  
  # Add color scale and theming
  p <- p +
    scale_color_viridis_c(option = color_option) +
    labs(title = title) + 
    theme_minimal(base_size = CFG_BASE_FONT_SIZE) +
    theme(plot.title = element_text(size = CFG_TITLE_FONT_SIZE, face = "bold"),
          axis.text.x = element_text(size = CFG_AXIS_TEXT_SIZE),
          axis.text.y = element_text(size = CFG_AXIS_TEXT_SIZE),
          axis.title = element_text(size = CFG_AXIS_TITLE_SIZE))
  
  return(p)
}

#===============================================================================
# GGPLOT THEMES
#===============================================================================

# Apply standard theme to ggplot objects
# Returns: theme object that can be added to ggplot
apply_standard_theme <- function(angle_x = 0) {
  theme_minimal(base_size = CFG_BASE_FONT_SIZE) +
    theme(plot.title = element_text(size = CFG_TITLE_FONT_SIZE, face = "bold"),
          axis.title = element_text(size = CFG_AXIS_TITLE_SIZE),
          axis.text.x = element_text(angle = angle_x, hjust = if(angle_x > 0) 1 else 0.5, 
                                     size = CFG_AXIS_TEXT_SIZE),
          axis.text.y = element_text(size = CFG_AXIS_TEXT_SIZE))
}

#===============================================================================
# CELL TYPE CENTERS VISUALIZATION
#===============================================================================

# Create cell type centers plot with labels
# Returns: ggplot object
create_celltype_centers_plot <- function(centers, title = "Cell Type Centers and Sizes") {
  p <- ggplot(centers, aes(x = UMAP1_center, y = UMAP2_center, label = cell_type)) +
    geom_point(aes(size = n_cells, color = cell_type), alpha = 0.85, show.legend = FALSE) +
    geom_text(size = 3, hjust = -0.2, vjust = -0.2) +
    scale_size_continuous(name = "Cells per type", range = c(3, 15)) +
    theme_minimal(base_size = CFG_BASE_FONT_SIZE) +
    labs(title = title, x = "UMAP_1", y = "UMAP_2") +
    theme(plot.title = element_text(size = CFG_TITLE_FONT_SIZE, face = "bold"),
          axis.title = element_text(size = CFG_AXIS_TITLE_SIZE))
  
  return(p)
}

#===============================================================================
# CLUSTERING QUALITY METRICS VISUALIZATION
#===============================================================================

# Create clustering quality metrics plots
# Returns: List with quality_metrics and quality_vs_clusters ggplot objects
create_quality_metrics_plots <- function(quality_metrics) {
  plots <- list()
  
  # Quality metrics comparison plot
  quality_long <- quality_metrics %>%
    dplyr::select(Resolution, AvgSilhouette, InertiaRatio) %>%
    tidyr::pivot_longer(cols = c(AvgSilhouette, InertiaRatio), 
                       names_to = "Metric", values_to = "Value") %>%
    dplyr::mutate(Metric = case_when(
      Metric == "AvgSilhouette" ~ "Average Silhouette Score",
      Metric == "InertiaRatio" ~ "Inertia Ratio (BCSS/TSS)",
      TRUE ~ Metric
    ))
  
  plots$quality_metrics <- ggplot(quality_long, aes(x = Resolution, y = Value, color = Metric)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~Metric, scales = "free_y") +
    labs(title = "Clustering Quality Metrics by Resolution",
         subtitle = "Higher values indicate better clustering quality",
         x = "Resolution", y = "Metric Value") +
    apply_standard_theme() +
    theme(legend.position = "none")
  
  # Cluster count vs quality
  plots$quality_vs_clusters <- ggplot(quality_metrics, aes(x = NumClusters, y = AvgSilhouette)) +
    geom_point(aes(size = InertiaRatio, color = Resolution), alpha = 0.7) +
    geom_smooth(method = "loess", se = FALSE, color = "gray50", linetype = "dashed") +
    scale_size_continuous(range = c(2, 8), guide = guide_legend(title = "Inertia Ratio")) +
    labs(title = "Clustering Quality: Silhouette Score vs Number of Clusters",
         subtitle = "Point size = Inertia Ratio | Ideal: high silhouette + optimal cluster count",
         x = "Number of Clusters", y = "Average Silhouette Score") +
    apply_standard_theme()
  
  return(plots)
}

#===============================================================================
# CLUSTERING STABILITY PLOT
#===============================================================================

# Create clustering stability plot
# Returns: ggplot object
create_stability_plot <- function(stab, selected_resolution, n_clusters) {
  p <- ggplot(stab, aes(x = Resolution, y = Clusters)) +
    geom_line(color = "steelblue", linewidth = 1.1) +
    geom_point(color = "darkblue", size = 3) +
    {if (!is.null(selected_resolution)) geom_vline(xintercept = selected_resolution, linetype = "dashed", color = "red", linewidth = 1)} +
    labs(
      title = "Clustering Stability Analysis",
      subtitle = paste0("Selected resolution: ", selected_resolution, " | ", n_clusters, " clusters"),
      x = "Resolution Parameter",
      y = "Number of Clusters Identified"
    ) +
    apply_standard_theme()
  
  return(p)
}

#===============================================================================
# CONFIDENCE DISTRIBUTION PLOT
#===============================================================================

# Create cell type confidence distribution plot
# Returns: ggplot object
create_confidence_plot <- function(metadata) {
  # Calculate mean and median confidence per cell type
  conf_summary <- metadata %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(
      mean_confidence = mean(celltype_confidence, na.rm = TRUE),
      median_confidence = median(celltype_confidence, na.rm = TRUE),
      sd_confidence = sd(celltype_confidence, na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(mean_confidence))
  
  # Reorder cell types by mean confidence
  conf_summary$cell_type <- factor(conf_summary$cell_type, 
                                   levels = conf_summary$cell_type)
  
  # Overall statistics
  overall_mean <- mean(metadata$celltype_confidence, na.rm = TRUE)
  
  p <- ggplot(conf_summary, aes(x = cell_type, y = mean_confidence)) +
    geom_col(aes(fill = mean_confidence), alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_confidence - sd_confidence, 
                      ymax = mean_confidence + sd_confidence),
                  width = 0.3, linewidth = 0.5) +
    geom_hline(yintercept = overall_mean, 
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_text(aes(label = n_cells), vjust = -0.5, size = 3.5) +
    scale_fill_gradient(low = "lightblue", high = "darkblue", 
                       name = "Mean\nConfidence") +
    labs(title = "Cell Type Assignment Confidence by Cell Type",
         subtitle = paste("Overall mean confidence:", round(overall_mean, 3), 
                         "| Numbers show cell counts"),
         x = "Cell Type", y = "Mean Confidence Score") +
    apply_standard_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
  
  return(p)
}

#===============================================================================
# ANNOTATION COMPARISON MATRIX PLOT
#===============================================================================

# Create annotation comparison matrix plot (confusion matrix with absolute counts + proportions)
# Returns: ggplot object
create_comparison_matrix_plot <- function(seu, percent_by = c("data_driven", "individual", "overall")) {
  percent_by <- match.arg(percent_by)

  # Create comparison table
  comparison_table <- table(
    Individual = seu$individual_celltype,
    DataDriven = seu$data_driven_celltype
  )

  # Convert to data frame for plotting
  comparison_df <- as.data.frame(comparison_table)

  # Convert factors to character to avoid level mismatch issues
  comparison_df$Individual <- as.character(comparison_df$Individual)
  comparison_df$DataDriven <- as.character(comparison_df$DataDriven)

  # Calculate proportions
  if (percent_by == "data_driven") {
    comparison_df <- comparison_df %>%
      dplyr::group_by(DataDriven) %>%
      dplyr::mutate(Proportion = ifelse(sum(Freq) > 0, Freq / sum(Freq), 0)) %>%
      dplyr::ungroup()
    percent_label <- "Proportion within data-driven type"
  } else if (percent_by == "individual") {
    comparison_df <- comparison_df %>%
      dplyr::group_by(Individual) %>%
      dplyr::mutate(Proportion = ifelse(sum(Freq) > 0, Freq / sum(Freq), 0)) %>%
      dplyr::ungroup()
    percent_label <- "Proportion within individual type"
  } else {
    total_count <- sum(comparison_df$Freq)
    comparison_df$Proportion <- ifelse(total_count > 0, comparison_df$Freq / total_count, 0)
    percent_label <- "Proportion of all cells"
  }

  # Create label with counts + percentages
  comparison_df$Label <- paste0(
    comparison_df$Freq, "\n(",
    round(comparison_df$Proportion * 100, 1), "%)"
  )

  # Calculate agreement (diagonal matches)
  # Use character comparison instead of factor comparison to avoid level mismatch
  diagonal_sum <- sum(comparison_df[comparison_df$Individual == comparison_df$DataDriven, "Freq"])
  total_cells <- sum(comparison_df$Freq)
  agreement_pct <- round(100 * diagonal_sum / total_cells, 1)

  p <- ggplot(comparison_df,
              aes(x = DataDriven, y = Individual, fill = Proportion)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = Label), size = 2.5, color = "black") +
    scale_fill_viridis_c(option = "plasma", name = "Proportion",
                         labels = scales::percent_format(accuracy = 1)) +
    labs(title = "Annotation Comparison Matrix (Cell Counts + Percentages)",
         subtitle = paste0(percent_label, " | Agreement Rate: ", diagonal_sum, " / ", 
                          total_cells, " cells (", agreement_pct, "%)"),
         x = "Data-Driven Cell Type (Cluster Consensus)",
         y = "Individual Cell Type (SingleR)") +
    apply_standard_theme(angle_x = 45) +
    theme(panel.grid = element_blank())

  return(p)
}

#===============================================================================
# CELL TYPE BAR PLOT COMPARISON
#===============================================================================

# Create side-by-side bar plot comparison of cell type distributions
# Returns: patchwork object combining two bar plots
create_celltype_barplot_comparison <- function(seu) {
  # Individual cell types (top 20)
  individual_counts <- as.data.frame(table(seu$individual_celltype)) %>%
    dplyr::arrange(desc(Freq)) %>%
    head(20) %>%
    dplyr::mutate(Type = "Individual (SingleR)")
  
  # Data-driven cell types (all)
  datadriven_counts <- as.data.frame(table(seu$data_driven_celltype)) %>%
    dplyr::mutate(Type = "Data-Driven (Consensus)")
  
  colnames(individual_counts)[1] <- "CellType"
  colnames(datadriven_counts)[1] <- "CellType"
  
  # Count total unique cell types
  n_individual_total <- length(unique(seu$individual_celltype))
  n_datadriven_total <- length(unique(seu$data_driven_celltype))
  
  p_bar_individual <- ggplot(individual_counts, aes(x = reorder(CellType, Freq), y = Freq)) +
    geom_bar(stat = "identity", fill = "coral", alpha = 0.8) +
    coord_flip() +
    labs(title = "Top 20 Individual Cell Types", 
         subtitle = paste("Total:", n_individual_total, "unique cell types"),
         x = "", y = "Number of Cells") +
    apply_standard_theme() +
    theme(axis.text = element_text(size = 11))
  
  p_bar_datadriven <- ggplot(datadriven_counts, aes(x = reorder(CellType, Freq), y = Freq)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    labs(title = "Data-Driven Cell Types",
         subtitle = paste("Total:", n_datadriven_total, "unique cell types"),
         x = "", y = "Number of Cells") +
    apply_standard_theme() +
    theme(axis.text = element_text(size = 11))
  
  # Combine with patchwork
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_combined <- p_bar_individual | p_bar_datadriven
    p_combined <- p_combined + 
      patchwork::plot_annotation(
        title = "Cell Type Distribution Comparison",
        theme = theme(plot.title = element_text(size = CFG_TITLE_FONT_SIZE + 2, 
                                                face = "bold", hjust = 0.5))
      )
    return(p_combined)
  } else {
    return(list(individual = p_bar_individual, datadriven = p_bar_datadriven))
  }
}

#===============================================================================
# INTERACTIVE CLUSTER PLOT
#===============================================================================

# Create interactive plotly cluster visualization
# Returns: plotly object
create_interactive_cluster_plot <- function(seu) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    warning("plotly package not available")
    return(NULL)
  }
  
  um <- as.data.frame(Seurat::Embeddings(seu, "umap"))
  colnames(um) <- c("UMAP_1", "UMAP_2") 
  plot_data <- cbind(seu@meta.data, um)
  
  p_interactive <- plotly::plot_ly(
    data = plot_data,
    x = ~UMAP_1, y = ~UMAP_2,
    color = ~seurat_clusters,
    colors = "Set3",
    type = "scatter", mode = "markers",
    marker = list(size = 3, opacity = 0.7),
    text = ~paste(
      "Cell:", rownames(plot_data),
      "<br>Cluster:", seurat_clusters, 
      "<br>Type:", cell_type,
      "<br>Genes:", nFeature_RNA,
      "<br>UMI:", nCount_RNA,
      "<br>Mito %:", round(percent.mt, 2)
    ),
    hovertemplate = "%{text}<extra></extra>"
  ) %>%
    plotly::layout(
      title = list(text = "Interactive Cluster Visualization", x = 0.5),
      xaxis = list(title = "UMAP 1"),
      yaxis = list(title = "UMAP 2")
    )
  
  return(p_interactive)
}

#===============================================================================
# CELL CYCLE DISTRIBUTION BARPLOT
#===============================================================================

# Create cell cycle distribution bar plot by cell type
# Returns: ggplot object
create_cell_cycle_barplot <- function(phase_by_cluster, phase_colors) {
  p <- ggplot(phase_by_cluster, aes(x = cell_type, y = percentage, fill = Phase)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = phase_colors) +
    labs(title = "Cell Cycle Distribution by Cell Type",
         x = "Cell Type", y = "Percentage (%)") +
    apply_standard_theme(angle_x = 45)
  
  return(p)
}

#===============================================================================
# INDIVIDUAL CELL CYCLE PLOTS
#===============================================================================

# Create combined cell cycle distribution plot for all cell types
# 
# This function generates a single faceted plot showing all cell types together
#
# @param seu Seurat object with Phase column
# @param phase_colors Named vector of colors for cell cycle phases
# @param plots_dir Directory to save static plots
# @param dash_dir Directory to save interactive plots
# @param group_by Column name to group by ("cell_type" or "seurat_clusters")
# @return Invisible NULL
create_cell_cycle_individual_plots <- function(seu, phase_colors, plots_dir, dash_dir, 
                                               group_by = "cell_type") {
  
  if (!group_by %in% colnames(seu@meta.data)) {
    log(paste("Warning: Column", group_by, "not found in metadata"))
    return(invisible(NULL))
  }
  
  if (!"Phase" %in% colnames(seu@meta.data)) {
    log("Warning: Phase column not found - cell cycle scoring may not have been performed")
    return(invisible(NULL))
  }
  
  # Prepare data
  phase_data <- seu@meta.data %>%
    dplyr::group_by(!!sym(group_by), Phase) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    dplyr::group_by(!!sym(group_by)) %>%
    dplyr::mutate(percentage = count / sum(count) * 100,
                  total_cells = sum(count))
  
  # Create overall barplot
  if (group_by == "seurat_clusters") {
    p_overall <- ggplot(phase_data, aes(x = !!sym(group_by), y = percentage, fill = Phase)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = phase_colors) +
      labs(title = "Cell Cycle Distribution by Seurat Cluster",
           x = "Cluster", y = "Percentage (%)") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    save_dual_plot(p_overall, "cell_cycle_by_seurat_cluster", plots_dir, dash_dir, width = 14, height = 7)
  }
  
  # Create combined faceted plot for all cell types in one view
  unique_groups <- unique(seu@meta.data[[group_by]])
  n_groups <- length(unique_groups)
  log(paste("Creating combined cell cycle plot for", n_groups, group_by, "groups"))
  
  # Calculate ncol and nrow for facet layout
  ncol_facet <- ceiling(sqrt(n_groups))
  nrow_facet <- ceiling(n_groups / ncol_facet)
  
  # Calculate dynamic plot dimensions
  plot_width <- max(14, ncol_facet * 4)
  plot_height <- max(10, nrow_facet * 3.5)
  
  # Create combined faceted plot
  if (group_by == "cell_type") {
    title_text <- "Cell Cycle Distribution Across All Cell Types"
    filename <- "cell_cycle_all_celltypes_combined"
  } else {
    title_text <- "Cell Cycle Distribution Across All Clusters"
    filename <- "cell_cycle_all_clusters_combined"
  }
  
  p_combined <- ggplot(phase_data, aes(x = Phase, y = percentage, fill = Phase)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.0f%%", percentage)), 
              vjust = -0.3, size = 2.5) +
    scale_fill_manual(values = phase_colors) +
    facet_wrap(as.formula(paste("~", group_by)), ncol = ncol_facet, scales = "fixed") +
    labs(title = title_text,
         x = "Cell Cycle Phase", y = "Percentage (%)") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      strip.text = element_text(size = 10, face = "bold"),
      panel.spacing = unit(0.5, "lines")
    ) +
    ylim(0, 100)
  
  # Add cell counts as facet labels
  p_combined <- p_combined +
    geom_text(data = phase_data %>% 
                dplyr::select(!!sym(group_by), total_cells) %>% 
                dplyr::distinct(),
              aes(x = 2, y = 95, label = paste0("n=", total_cells)),
              inherit.aes = FALSE, size = 3, color = "gray30")
  
  save_dual_plot(p_combined, filename, plots_dir, dash_dir, 
                width = plot_width, height = plot_height)
  
  log(paste("Created combined cell cycle plot for", n_groups, "groups"))
  return(invisible(NULL))
}

#===============================================================================
# ANNOTATION COMPARISON SIDE-BY-SIDE
#===============================================================================

# Create side-by-side comparison of individual vs data-driven cell type annotations
# Returns: patchwork object
create_annotation_comparison_sidebyside <- function(seu, individual_colors, consistent_colors) {
  p_individual <- create_umap_plot(seu, group_by = "individual_celltype", 
                                   title = "Individual Cell Type\n(SingleR Fine Labels)",
                                   colors = individual_colors, add_marginal = FALSE) +
    theme(legend.position = "none")
  
  p_datadriven <- create_umap_plot(seu, group_by = "data_driven_celltype", 
                                   title = "Data-Driven Cell Type\n(Cluster Consensus)",
                                   label = TRUE, repel = TRUE, colors = consistent_colors,
                                   add_marginal = FALSE)
  
  p_comparison <- p_individual | p_datadriven
  p_comparison <- p_comparison + 
    patchwork::plot_annotation(
      title = "Cell Type Annotation Comparison",
      theme = theme(plot.title = element_text(size = CFG_TITLE_FONT_SIZE + 2, 
                                               face = "bold", hjust = 0.5))
    )
  
  return(p_comparison)
}

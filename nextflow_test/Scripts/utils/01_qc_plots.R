#===============================================================================
# 01_qc/plots.R
# QC and Integration specific plotting functions
#===============================================================================

# Plot functions specific to QC and integration analysis
# - Quality control metrics visualization
# - Sample comparison plots
# - Integration diagnostics

#' Save QC plot with error handling
#' 
#' Specialized function for saving QC plots. Saves to both plots_dir and dash_dir.
#' Creates fallback empty plot if primary plotting fails.
#' 
#' @param p ggplot object
#' @param filename Output filename (will add .png extension if needed)
#' @param width Plot width in inches (default: 9)
#' @param height Plot height in inches (default: 4)
#' @param dpi Resolution for PNG (default: 300)
#' @param plots_dir Output directory for static plots (default: from parent environment)
#' @param dash_dir Output directory for dashboard (default: from parent environment)
#' @return Invisible path to PNG file
#' @examples
#' p <- ggplot(mtcars, aes(mpg, wt, color = factor(cyl))) + geom_point()
#' save_qc_plot(p, "qc_scatter.png", width = 10, height = 6)
save_qc_plot <- function(p, filename, width = 9, height = 4, dpi = 300,
                        plots_dir = NULL, dash_dir = NULL) {
  
  # Get directories from parent environment if not provided
  if (is.null(plots_dir)) {
    if (exists("plots_dir", envir = parent.frame(), inherits = FALSE)) {
      plots_dir <- get("plots_dir", envir = parent.frame())
    } else {
      plots_dir <- "."
    }
  }
  
  if (is.null(dash_dir)) {
    if (exists("dash_dir", envir = parent.frame(), inherits = FALSE)) {
      dash_dir <- get("dash_dir", envir = parent.frame())
    } else {
      dash_dir <- plots_dir  # Fallback to plots_dir if dash_dir not found
    }
  }
  
  # Ensure directories exist
  for (d in c(plots_dir, dash_dir)) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # Ensure .png extension
  if (!grepl("\\.(png|pdf)$", filename)) {
    filename <- paste0(filename, ".png")
  }
  
  tryCatch({
    # Save PNG to plots_dir
    png_path <- file.path(plots_dir, filename)
    suppressMessages(
      ggplot2::ggsave(png_path, plot = p, width = width, height = height, 
                     dpi = dpi, bg = "white")
    )
    
    # Copy PNG to dash_dir for dashboard (only if different from plots_dir)
    dash_path <- file.path(dash_dir, filename)
    if (normalizePath(png_path, mustWork = FALSE) != normalizePath(dash_path, mustWork = FALSE)) {
      file.copy(png_path, dash_path, overwrite = TRUE)
    }
    
    # Save PDF to plots_dir only
    pdf_filename <- sub("\\.png$", ".pdf", filename)
    pdf_path <- file.path(plots_dir, pdf_filename)
    suppressMessages(
      ggplot2::ggsave(pdf_path, plot = p, width = width, height = height, 
                     device = "pdf")
    )
    
    if (exists("log", mode = "function")) {
      log("Saved QC plot:", filename, "and", pdf_filename)
    }
    
    return(invisible(png_path))
    
  }, error = function(e) {
    if (exists("log", mode = "function")) {
      log("ERROR saving QC plot", filename, ":", e$message)
    } else {
      warning("Error saving QC plot ", filename, ": ", e$message)
    }
    
    # Create empty plot as fallback
    tryCatch({
      empty_plot <- ggplot2::ggplot() + 
        ggplot2::annotate("text", x = 0.5, y = 0.5, 
                         label = paste("Plot failed:", filename), hjust = 0.5) +
        ggplot2::theme_void()
      png_path <- file.path(plots_dir, filename)
      ggplot2::ggsave(png_path, plot = empty_plot, width = width, height = height, dpi = dpi)
    }, error = function(e2) {
      if (exists("log", mode = "function")) {
        log("ERROR: Could not even save fallback plot:", e2$message)
      }
    })
    
    return(invisible(NULL))
  })
}


#' Save plot in dual format (PNG + Interactive HTML)
#' 
#' Saves a ggplot as both a static PNG (in plots_dir) and an interactive 
#' HTML file (in dash_dir). This is the STANDARD plotting function for the pipeline.
#' 
#' @param p ggplot object
#' @param filename Base filename without extension
#' @param plots_dir Directory for static plots (default: from parent environment)
#' @param dash_dir Directory for interactive plots (default: from parent environment)
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 6)
#' @param dpi Resolution for PNG (default: 150)
#' @return Invisible list with paths to PNG and HTML files
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(mpg, wt)) + geom_point()
#' save_dual_plot(p, "scatter_plot")
save_dual_plot <- function(p, filename, plots_dir = NULL, dash_dir = NULL,
                           width = 12, height = 6, dpi = 150) {
  
  # Get directories from parent environment if not provided
  if (is.null(plots_dir)) {
    if (exists("plots_dir", envir = parent.frame(), inherits = FALSE)) {
      plots_dir <- get("plots_dir", envir = parent.frame())
    } else {
      plots_dir <- "."
    }
  }
  
  if (is.null(dash_dir)) {
    if (exists("dash_dir", envir = parent.frame(), inherits = FALSE)) {
      dash_dir <- get("dash_dir", envir = parent.frame())
    } else {
      dash_dir <- plots_dir
    }
  }
  
  # Ensure directories exist
  for (d in c(plots_dir, dash_dir)) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # Remove any existing extension from filename
  base_name <- tools::file_path_sans_ext(filename)
  
  # Save static PNG
  png_filename <- paste0(base_name, ".png")
  png_path <- file.path(plots_dir, png_filename)
  suppressMessages(
    ggplot2::ggsave(png_path, plot = p, width = width, height = height, 
                   dpi = dpi, bg = "white")
  )
  
  # Save interactive HTML
  html_filename <- paste0(base_name, ".html")
  html_path <- file.path(dash_dir, html_filename)
  
  if (requireNamespace("plotly", quietly = TRUE) && 
      requireNamespace("htmlwidgets", quietly = TRUE)) {
    
    plot_interactive <- plotly::ggplotly(p)
    htmlwidgets::saveWidget(plot_interactive, html_path, selfcontained = TRUE)
  }
  
  if (exists("log", mode = "function")) {
    log("Saved dual-format plot:", png_filename, "&", html_filename)
  }
  
  invisible(list(png = png_path, html = html_path))
}


#' Create dimensionality reduction comparison plots for integrated data
#'
#' Generates UMAP and tSNE plots showing different metadata variables
#' (condition, sample, clusters) for an integrated Seurat object.
#'
#' @param seurat_obj Integrated Seurat object with UMAP and tSNE reductions
#' @param method_name Integration method name (for filenames)
#' @param dash_dir Directory for interactive HTML plots
#' @param plots_dir Directory for static PNG plots
#' @return Invisible NULL
plot_dimred_comparison <- function(seurat_obj, method_name, dash_dir, plots_dir) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    if (exists("log", mode = "function")) {
      log("WARNING: Seurat package not available for plotting")
    }
    return(invisible(NULL))
  }
  
  # Check if UMAP and tSNE reductions exist
  has_umap <- "umap" %in% names(seurat_obj@reductions)
  has_tsne <- "tsne" %in% names(seurat_obj@reductions)
  
  if (!has_umap && !has_tsne) {
    if (exists("log", mode = "function")) {
      log("WARNING: No UMAP or tSNE reductions found in object")
    }
    return(invisible(NULL))
  }
  
  # Define metadata variables to plot
  plot_vars <- c("condition", "sample", "seurat_clusters")
  
  # Add cell.type if it exists (from SingleR)
  if ("cell.type" %in% colnames(seurat_obj@meta.data)) {
    plot_vars <- c(plot_vars, "cell.type")
  }
  
  # Create plots for each reduction type and metadata variable
  for (reduction in c("umap", "tsne")) {
    if ((reduction == "umap" && !has_umap) || (reduction == "tsne" && !has_tsne)) {
      next
    }
    
    reduction_name <- toupper(reduction)
    
    for (var in plot_vars) {
      if (!var %in% colnames(seurat_obj@meta.data)) {
        next
      }
      
      tryCatch({
        # Create DimPlot
        p <- Seurat::DimPlot(seurat_obj, 
                            reduction = reduction,
                            group.by = var,
                            label = (var == "seurat_clusters" || var == "cell.type"),
                            label.size = 3,
                            repel = TRUE,
                            pt.size = 0.5) +
          ggplot2::ggtitle(paste(method_name, "-", reduction_name, "by", var)) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        
        # Save plot
        filename <- paste0(method_name, "_", reduction, "_", var)
        save_dual_plot(p, filename, 
                      plots_dir = plots_dir, 
                      dash_dir = dash_dir,
                      width = 10, height = 8)
        
      }, error = function(e) {
        if (exists("log", mode = "function")) {
          log(paste("WARNING: Could not create", reduction, "plot for", var, "-", e$message))
        }
      })
    }
  }
  
  invisible(NULL)
}


#' Generate comprehensive QC plots for combined dataset
#' 
#' Creates boxplots and scatter plots for QC metrics across all samples.
#' Generates both static and interactive visualizations.
#' 
#' @param qc_list List of data frames containing QC metrics per sample
#' @param suffix Label suffix ("before_QC" or "after_QC")
#' @param plots_dir Directory for static plots
#' @param dash_dir Directory for interactive plots
#' @return Invisible NULL
generate_integrated_qc_plots <- function(qc_list, suffix, plots_dir, dash_dir) {
  if (length(qc_list) == 0) return(NULL)
  
  # Combine all sample data
  df <- do.call(rbind, qc_list)
  if (exists("log", mode = "function")) {
    log(sprintf("Generating %s-processing plots for %d cells across %d samples", 
                suffix, nrow(df), length(unique(df$sample))))
  }
  
  # Get font size configs from parent environment (or use defaults)
  CFG_BASE_FONT_SIZE <- if (exists("CFG_BASE_FONT_SIZE", envir = parent.frame(), inherits = TRUE)) {
    get("CFG_BASE_FONT_SIZE", envir = parent.frame(), inherits = TRUE)
  } else { 16 }
  
  CFG_AXIS_TITLE_SIZE <- if (exists("CFG_AXIS_TITLE_SIZE", envir = parent.frame(), inherits = TRUE)) {
    get("CFG_AXIS_TITLE_SIZE", envir = parent.frame(), inherits = TRUE)
  } else { 14 }
  
  CFG_TITLE_FONT_SIZE <- if (exists("CFG_TITLE_FONT_SIZE", envir = parent.frame(), inherits = TRUE)) {
    get("CFG_TITLE_FONT_SIZE", envir = parent.frame(), inherits = TRUE)
  } else { 18 }
  
  CFG_MAX_CELLS_INTERACTIVE <- if (exists("CFG_MAX_CELLS_INTERACTIVE", envir = parent.frame(), inherits = TRUE)) {
    get("CFG_MAX_CELLS_INTERACTIVE", envir = parent.frame(), inherits = TRUE)
  } else { 30000 }

  # 1. Feature and MT Boxplots (Colored by Condition)
  p_feat_box <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = nFeature_RNA, fill = condition)) +
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::theme_minimal(base_size = CFG_BASE_FONT_SIZE) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = CFG_BASE_FONT_SIZE),
          axis.text.y = ggplot2::element_text(size = CFG_BASE_FONT_SIZE),
          axis.title = ggplot2::element_text(size = CFG_AXIS_TITLE_SIZE),
          plot.title = ggplot2::element_text(size = CFG_TITLE_FONT_SIZE)) +
    ggplot2::labs(title = paste("nFeature_RNA Distribution (", suffix, ")", sep=""), y = "nFeature_RNA")
  save_dual_plot(p_feat_box, paste0("boxplot_nFeature_", suffix), plots_dir, dash_dir)
  
  p_mt_box <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = percent.mt, fill = condition)) +
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::theme_minimal(base_size = CFG_BASE_FONT_SIZE) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = CFG_BASE_FONT_SIZE),
          axis.text.y = ggplot2::element_text(size = CFG_BASE_FONT_SIZE),
          axis.title = ggplot2::element_text(size = CFG_AXIS_TITLE_SIZE),
          plot.title = ggplot2::element_text(size = CFG_TITLE_FONT_SIZE)) +
    ggplot2::labs(title = paste("Mitochondrial % Distribution (", suffix, ")", sep=""), y = "% MT")
  save_dual_plot(p_mt_box, paste0("boxplot_percent_mt_", suffix), plots_dir, dash_dir)

  # 2. Scatter: nFeature vs nCount (All Samples)
  p_scatter <- ggplot2::ggplot(df, ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, color = condition)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.4) +
    ggplot2::theme_minimal(base_size = CFG_BASE_FONT_SIZE) +
    ggplot2::theme(axis.title = ggplot2::element_text(size = CFG_AXIS_TITLE_SIZE),
          plot.title = ggplot2::element_text(size = CFG_TITLE_FONT_SIZE)) +
    ggplot2::labs(title = paste("nFeature vs nCount (", suffix, ")", sep=""))
  save_dual_plot(p_scatter, paste0("scatter_nFeature_vs_nCount_", suffix), plots_dir, dash_dir, width = 10, height = 8)

  # 3. Interactive Scatter (Color by Sample) using Plotly
  if(requireNamespace("plotly", quietly=TRUE)) {
      df_sub <- if(nrow(df) > CFG_MAX_CELLS_INTERACTIVE) df[sample(nrow(df), CFG_MAX_CELLS_INTERACTIVE), ] else df
      p_inter <- plotly::plot_ly(df_sub, x = ~nCount_RNA, y = ~nFeature_RNA, color = ~sample,
                                 text = ~paste("Sample:", sample, "<br>Cond:", condition),
                                 type = 'scatter', mode = 'markers', marker = list(size = 3, opacity = 0.6)) %>%
        plotly::layout(
          title = list(text = paste("Interactive QC (", suffix, "): nFeature vs nCount", sep=""),
                      font = list(size = CFG_TITLE_FONT_SIZE)),
          font = list(size = CFG_BASE_FONT_SIZE),
          xaxis = list(titlefont = list(size = CFG_AXIS_TITLE_SIZE)),
          yaxis = list(titlefont = list(size = CFG_AXIS_TITLE_SIZE))
        )
      htmlwidgets::saveWidget(p_inter, file.path(dash_dir, paste0("interactive_qc_scatter_", suffix, ".html")), selfcontained = TRUE)
    }
  
  invisible(NULL)
}


#' Generate final quality metrics visualization plots
#' 
#' Creates summary plots showing overall dataset quality metrics and 
#' per-condition QC statistics after integration.
#' 
#' @param merged_seu Integrated Seurat object
#' @param primary_method Integration method name (for subtitle)
#' @param plots_dir Directory for static plots
#' @param dash_dir Directory for interactive plots
#' @return Invisible NULL
generate_final_quality_plots <- function(merged_seu, primary_method, plots_dir, dash_dir) {
  
  if (exists("log", mode = "function")) {
    log("Creating quality metrics visualization...")
  }
  
  # Get font size configs from parent environment (or use defaults)
  CFG_BASE_FONT_SIZE <- if (exists("CFG_BASE_FONT_SIZE", envir = parent.frame(), inherits = TRUE)) {
    get("CFG_BASE_FONT_SIZE", envir = parent.frame(), inherits = TRUE)
  } else { 16 }
  
  CFG_AXIS_TITLE_SIZE <- if (exists("CFG_AXIS_TITLE_SIZE", envir = parent.frame(), inherits = TRUE)) {
    get("CFG_AXIS_TITLE_SIZE", envir = parent.frame(), inherits = TRUE)
  } else { 14 }
  
  CFG_TITLE_FONT_SIZE <- if (exists("CFG_TITLE_FONT_SIZE", envir = parent.frame(), inherits = TRUE)) {
    get("CFG_TITLE_FONT_SIZE", envir = parent.frame(), inherits = TRUE)
  } else { 18 }
  
  # 1. Overall quality metrics summary plot
  final_metrics_df <- data.frame(
    Metric = c("Total Cells", "Total Genes", "Conditions", "Samples", "Variable Features"),
    Value = c(
      ncol(merged_seu),
      nrow(merged_seu), 
      length(unique(merged_seu$condition)),
      length(unique(merged_seu$sample)),
      length(Seurat::VariableFeatures(merged_seu))
    )
  )
  
  p_final_metrics <- ggplot2::ggplot(final_metrics_df, ggplot2::aes(x = Metric, y = Value)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = format(Value, big.mark = ",")), 
                      vjust = -0.3, fontface = "bold", size = 5) +
    ggplot2::theme_minimal(base_size = CFG_BASE_FONT_SIZE) +
    ggplot2::labs(title = "Final Dataset Quality Metrics",
         subtitle = paste("Integration method:", ifelse(!is.null(primary_method), primary_method, "None")),
         x = "Quality Metric", y = "Count") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = CFG_BASE_FONT_SIZE),
          axis.text.y = ggplot2::element_text(size = CFG_BASE_FONT_SIZE),
          axis.title = ggplot2::element_text(size = CFG_AXIS_TITLE_SIZE),
          plot.title = ggplot2::element_text(size = CFG_TITLE_FONT_SIZE)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15)))
  
  ggplot2::ggsave(file.path(plots_dir, "final_quality_metrics.png"), 
                 p_final_metrics, width = 10, height = 6, dpi = 150)
  
  if (requireNamespace("plotly", quietly = TRUE) && 
      requireNamespace("htmlwidgets", quietly = TRUE)) {
    p_final_interactive <- plotly::ggplotly(p_final_metrics)
    htmlwidgets::saveWidget(p_final_interactive, 
                           file.path(dash_dir, "final_quality_metrics.html"), 
                           selfcontained = TRUE)
  }
  
  # 2. QC metrics by condition summary
  if (requireNamespace("dplyr", quietly = TRUE)) {
    qc_condition_summary <- merged_seu@meta.data %>%
      dplyr::group_by(condition) %>%
      dplyr::summarise(
        n_cells = dplyr::n(),
        mean_nFeature = mean(nFeature_RNA, na.rm = TRUE),
        mean_nCount = mean(nCount_RNA, na.rm = TRUE),
        mean_mt = mean(percent.mt, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Create condition summary plot
    p_qc_by_condition <- ggplot2::ggplot(qc_condition_summary, ggplot2::aes(x = condition)) +
      ggplot2::geom_col(ggplot2::aes(y = mean_nFeature), alpha = 0.8, fill = "steelblue") +
      ggplot2::labs(title = "Mean Features by Condition",
           x = "Condition", y = "Mean nFeature_RNA") +
      ggplot2::theme_minimal(base_size = CFG_BASE_FONT_SIZE) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = CFG_BASE_FONT_SIZE),
            axis.text.y = ggplot2::element_text(size = CFG_BASE_FONT_SIZE),
            axis.title = ggplot2::element_text(size = CFG_AXIS_TITLE_SIZE),
            plot.title = ggplot2::element_text(size = CFG_TITLE_FONT_SIZE))
    
    ggplot2::ggsave(file.path(plots_dir, "mean_features_by_condition.png"), 
                   p_qc_by_condition, width = 10, height = 6, dpi = 150)
    
    if (requireNamespace("plotly", quietly = TRUE) && 
        requireNamespace("htmlwidgets", quietly = TRUE)) {
      p_qc_interactive <- plotly::ggplotly(p_qc_by_condition)
      htmlwidgets::saveWidget(p_qc_interactive, 
                             file.path(dash_dir, "mean_features_by_condition.html"), 
                             selfcontained = TRUE)
    }
  }
  
  if (exists("log", mode = "function")) {
    log("Saved quality metrics visualization plots")
  }
  
  invisible(NULL)
}


#-------------------------------------------------------------------------------
# 4.1. Violin Plots for QC Metrics
#-------------------------------------------------------------------------------
# Create violin plots for QC metrics
# @param qc_list List of data frames with QC metrics per sample
# @param suffix String suffix for plot names ("before_QC" or "after_QC")
# @param plots_dir Directory to save static plots
# @param dash_dir Directory to save interactive plots

#===============================================================================
# additional_plots.R
# Additional plotting functions for Module 01 QC visualization
#===============================================================================

#-------------------------------------------------------------------------------
# 4.1. Violin Plots for QC Metrics
#-------------------------------------------------------------------------------
# Create violin plots for QC metrics
# @param qc_list List of data frames with QC metrics per sample
# @param suffix String suffix for plot names ("before_QC" or "after_QC")
# @param plots_dir Directory to save static plots
# @param dash_dir Directory to save interactive plots
generate_violin_plots <- function(qc_list, suffix, plots_dir, dash_dir) {
  df <- do.call(rbind, qc_list)
  
  # Violin plot: nCount_RNA
  p_count <- ggplot(df, aes(x = sample, y = nCount_RNA, fill = condition)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", alpha = 0.5) +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0("UMI Counts Distribution (", suffix, ")"), 
         y = "nCount_RNA", x = "Sample") +
    scale_y_log10()
  ggsave(file.path(plots_dir, paste0("violin_nCount_RNA_", suffix, ".png")), 
         p_count, width = 12, height = 6, dpi = 150)
  
  # Violin plot: nFeature_RNA
  p_feature <- ggplot(df, aes(x = sample, y = nFeature_RNA, fill = condition)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", alpha = 0.5) +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0("Gene Counts Distribution (", suffix, ")"), 
         y = "nFeature_RNA", x = "Sample") +
    scale_y_log10()
  ggsave(file.path(plots_dir, paste0("violin_nFeature_RNA_", suffix, ".png")), 
         p_feature, width = 12, height = 6, dpi = 150)
  
  # Violin plot: percent.mt
  p_mt <- ggplot(df, aes(x = sample, y = percent.mt, fill = condition)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", alpha = 0.5) +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0("Mitochondrial % Distribution (", suffix, ")"), 
         y = "% MT", x = "Sample")
  ggsave(file.path(plots_dir, paste0("violin_percent.mt_", suffix, ".png")), 
         p_mt, width = 12, height = 6, dpi = 150)
}

#-------------------------------------------------------------------------------
# 4.2. QC Heatmap
#-------------------------------------------------------------------------------
# Create heatmap matrix of QC metrics
# @param qc_list List of data frames with QC metrics
# @param suffix String suffix for plot name
# @param plots_dir Directory to save plots
generate_qc_heatmap <- function(qc_list, suffix, plots_dir) {
  # Calculate summary statistics per sample
  stats_df <- lapply(names(qc_list), function(sname) {
    df <- qc_list[[sname]]
    data.frame(
      sample = sname,
      condition = unique(df$condition)[1],
      n_cells = nrow(df),
      median_umi = median(df$nCount_RNA),
      median_genes = median(df$nFeature_RNA),
      median_mt = median(df$percent.mt),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  # Reshape for heatmap
  mat <- stats_df %>%
    select(sample, median_umi, median_genes, median_mt) %>%
    tibble::column_to_rownames("sample") %>%
    as.matrix() %>%
    t()
  
  # Scale each metric
  mat_scaled <- t(scale(t(mat)))
  
  # Plot heatmap
  png(file.path(plots_dir, paste0("heatmatrix_qc_", suffix, ".png")), 
      width = 10, height = 6, units = "in", res = 150)
  pheatmap::pheatmap(mat_scaled, 
                     cluster_rows = FALSE,
                     cluster_cols = TRUE,
                     color = colorRampPalette(c("blue", "white", "red"))(50),
                     main = paste0("QC Metrics Heatmap (", suffix, ")"),
                     fontsize = 12,
                     angle_col = 45)
  dev.off()
}

#-------------------------------------------------------------------------------
# 4.4. Filtering Summary
#-------------------------------------------------------------------------------
# Plot filtering summary statistics
# @param filtering_stats List of filtering stats per sample
# @param plots_dir Directory to save plots
generate_filtering_summary <- function(filtering_stats, plots_dir) {
  # Convert to data frame
  df <- lapply(names(filtering_stats), function(sname) {
    s <- filtering_stats[[sname]]
    data.frame(
      sample = s$sample,
      condition = s$condition,
      cells_before = s$cells_before,
      cells_after = s$cells_after,
      retention_rate = s$cells_after / s$cells_before * 100,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  # Bar plot of cell retention
  p <- ggplot(df, aes(x = sample, y = retention_rate, fill = condition)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.1f%%", retention_rate)), 
              vjust = -0.5, size = 3) +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Cell Retention After QC Filtering",
         y = "Retention Rate (%)", x = "Sample") +
    ylim(0, 100)
  
  ggsave(file.path(plots_dir, "filtering_summary.png"), 
         p, width = 12, height = 6, dpi = 150)
}

#-------------------------------------------------------------------------------
# 4.5. Doublet Detection Summary
#-------------------------------------------------------------------------------
# Plot doublet removal summary from saved statistics
# @param doublet_stats List of doublet statistics from detect_doublets
# @param plots_dir Directory to save plots
generate_doublet_summary_from_stats <- function(doublet_stats, plots_dir) {
  # Extract doublet info from stats
  df <- lapply(names(doublet_stats), function(sname) {
    stats <- doublet_stats[[sname]]
    if(!is.null(stats)) {
      data.frame(
        sample = stats$sample_name,
        cells_before = stats$cells_before,
        cells_after = stats$cells_after,
        doublets_detected = stats$doublets_detected,
        doublet_rate = stats$percent_removed,
        stringsAsFactors = FALSE
      )
    }
  }) %>% 
    bind_rows()
  
  if(nrow(df) == 0) {
    # No doublet data available
    return(invisible(NULL))
  }
  
  # Reshape for plotting
  df_plot <- df %>%
    select(sample, cells_after, doublets_detected) %>%
    pivot_longer(cols = c(cells_after, doublets_detected), 
                 names_to = "type", values_to = "count")
  
  df_plot$type <- factor(df_plot$type, 
                         levels = c("cells_after", "doublets_detected"),
                         labels = c("Singlets (retained)", "Doublets (removed)"))
  
  # Create stacked bar plot
  p <- ggplot(df_plot, aes(x = sample, y = count, fill = type)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(data = df, 
              aes(x = sample, y = cells_before, 
                  label = sprintf("%d doublets\n(%.1f%%)", doublets_detected, doublet_rate)),
              inherit.aes = FALSE, vjust = -0.3, size = 3, lineheight = 0.9) +
    scale_fill_manual(values = c("Singlets (retained)" = "steelblue", 
                                  "Doublets (removed)" = "coral"),
                      name = "Cell Type") +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
    labs(title = "Doublet Detection Summary (scDblFinder)",
         subtitle = "Showing cells retained vs doublets removed per sample",
         y = "Number of Cells", x = "Sample")
  
  ggsave(file.path(plots_dir, "doublet_removal_summary.png"), 
         p, width = 14, height = 7, dpi = 150)
}

#' Plot doublet removal summary
#' @param sample_objects List of Seurat objects after doublet removal
#' @param plots_dir Directory to save plots
generate_doublet_summary <- function(sample_objects, plots_dir) {
  # Extract doublet info from metadata
  # Note: After filtering, doublets are removed, but we need original counts
  # This function will show what was detected, even if doublets are already removed
  
  df <- lapply(names(sample_objects), function(sname) {
    seu <- sample_objects[[sname]]
    
    # Check if doublet metadata exists (it should if scDblFinder was run)
    if("scDblFinder.class" %in% colnames(seu@meta.data)) {
      # These are cells AFTER doublet removal (only singlets remain)
      n_singlets <- ncol(seu)
      n_doublets <- 0  # Already removed
      
      # Try to estimate doublets from score distribution
      # Typical doublet rate is 0.8-8% depending on cell count
      expected_doublet_rate <- min(0.08, 0.008 * (n_singlets / 1000))
      estimated_doublets <- round(n_singlets * expected_doublet_rate / (1 - expected_doublet_rate))
      
      data.frame(
        sample = sname,
        condition = unique(seu$condition)[1],
        cells_retained = n_singlets,
        doublets_removed = estimated_doublets,
        doublet_rate = estimated_doublet_rate * 100,
        note = "estimated",
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        sample = sname,
        condition = unique(seu$condition)[1],
        cells_retained = ncol(seu),
        doublets_removed = 0,
        doublet_rate = 0,
        note = "not_run",
        stringsAsFactors = FALSE
      )
    }
  }) %>% bind_rows()
  
  # Create plot showing before and after
  df_plot <- df %>%
    mutate(cells_before = cells_retained + doublets_removed) %>%
    select(sample, condition, cells_before, cells_retained, doublets_removed) %>%
    pivot_longer(cols = c(cells_retained, doublets_removed), 
                 names_to = "type", values_to = "count")
  
  df_plot$type <- factor(df_plot$type, levels = c("cells_retained", "doublets_removed"),
                         labels = c("Singlets (retained)", "Doublets (removed)"))
  
  p <- ggplot(df_plot, aes(x = sample, y = count, fill = type)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(data = df, 
              aes(x = sample, y = cells_retained + doublets_removed, 
                  label = sprintf("%.1f%%\nremoved", doublet_rate)),
              inherit.aes = FALSE, vjust = -0.3, size = 3, lineheight = 0.9) +
    scale_fill_manual(values = c("Singlets (retained)" = "steelblue", 
                                  "Doublets (removed)" = "coral"),
                      name = "Cell Type") +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
    labs(title = "Doublet Detection Summary",
         subtitle = "Estimated doublet rates based on typical scDblFinder performance",
         y = "Number of Cells", x = "Sample")
  
  ggsave(file.path(plots_dir, "doublet_removal_summary.png"), 
         p, width = 12, height = 7, dpi = 150)
}

#-------------------------------------------------------------------------------
# 4.6. Ambient RNA Contamination
#-------------------------------------------------------------------------------
# Plot ambient RNA contamination summary
# @param soupx_stats List of SoupX statistics per sample
# @param plots_dir Directory to save plots
generate_ambient_rna_plot <- function(soupx_stats, plots_dir) {
  # Filter out NULL entries and extract contamination info
  df <- lapply(names(soupx_stats), function(sname) {
    if(!is.null(soupx_stats[[sname]])) {
      data.frame(
        sample = sname,
        contamination = soupx_stats[[sname]]$contamination * 100,
        umi_reduction = soupx_stats[[sname]]$umi_reduction,
        stringsAsFactors = FALSE
      )
    }
  }) %>% 
    bind_rows()
  
  if(nrow(df) == 0) {
    # No SoupX data available
    return(invisible(NULL))
  }
  
  # Bar plot of contamination
  p <- ggplot(df, aes(x = sample, y = contamination)) +
    geom_bar(stat = "identity", fill = "darkseagreen") +
    geom_text(aes(label = sprintf("%.1f%%", contamination)), 
              vjust = -0.5, size = 3) +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Ambient RNA Contamination (SoupX)",
         y = "Contamination (%)", x = "Sample")
  
  ggsave(file.path(plots_dir, "ambient_rna_contamination.png"), 
         p, width = 12, height = 6, dpi = 150)
}

#-------------------------------------------------------------------------------
# 4.7. HVG Selection Plots
#-------------------------------------------------------------------------------
# Create HVG selection plots per condition
# @param sample_objects List of Seurat objects
# @param plots_dir Directory to save plots
generate_hvg_plots <- function(sample_objects, plots_dir) {
  # Get unique conditions
  conditions <- unique(sapply(sample_objects, function(x) unique(x$condition)[1]))
  
  # Plot HVGs for each condition
  for(cond in conditions) {
    # Get samples for this condition
    cond_samples <- names(sample_objects)[sapply(sample_objects, function(x) unique(x$condition)[1] == cond)]
    
    if(length(cond_samples) > 0) {
      seu <- sample_objects[[cond_samples[1]]]
      
      # Plot variable features
      top10 <- head(VariableFeatures(seu), 10)
      
      tryCatch({
        p <- VariableFeaturePlot(seu) +
          labs(title = paste("Variable Features -", cond)) +
          theme_minimal(base_size = 16)
        
        # Try to add labels if gene metadata exists
        if(requireNamespace("ggrepel", quietly = TRUE) && length(top10) > 0) {
          # Get HVF info from the correct slot
          hvf_info <- HVFInfo(seu)
          if(!is.null(hvf_info) && nrow(hvf_info) > 0 && all(top10 %in% rownames(hvf_info))) {
            label_data <- hvf_info[top10, , drop = FALSE]
            if("mean" %in% colnames(label_data) && "variance.standardized" %in% colnames(label_data)) {
              label_data$gene <- rownames(label_data)
              p <- p + ggrepel::geom_label_repel(
                data = label_data,
                aes(x = mean, y = variance.standardized, label = gene),
                size = 3, max.overlaps = 20
              )
            }
          }
        }
        
        ggsave(file.path(plots_dir, paste0("hvg_selection_", gsub("-", "_", cond), ".png")), 
               p, width = 10, height = 8, dpi = 150)
      }, error = function(e) {
        warning("Could not generate HVG plot for ", cond, ": ", e$message)
      })
    }
  }
  
  # Global HVG selection plot (from first sample as example)
  if(length(sample_objects) > 0) {
    seu <- sample_objects[[1]]
    top10 <- head(VariableFeatures(seu), 10)
    
    p_global <- VariableFeaturePlot(seu) +
      labs(title = "Variable Features - Global Selection") +
      theme_minimal(base_size = 16)
    
    ggsave(file.path(plots_dir, "hvg_selection_global.png"), 
           p_global, width = 10, height = 8, dpi = 150)
  }
}

#' Create HVG dotplot showing top variable genes per condition
#' @param sample_objects List of Seurat objects
#' @param plots_dir Directory to save plots
#' @param top_n Number of top HVGs to show
generate_hvg_dotplot <- function(sample_objects, plots_dir, top_n = 10) {
  # Get top HVGs per condition
  hvg_per_cond <- lapply(names(sample_objects), function(sname) {
    seu <- sample_objects[[sname]]
    cond <- unique(seu$condition)[1]
    top_hvgs <- head(VariableFeatures(seu), top_n)
    data.frame(
      gene = top_hvgs,
      condition = cond,
      rank = 1:length(top_hvgs),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows() %>%
    group_by(condition, gene) %>%
    summarize(rank = mean(rank), .groups = "drop")
  
  # Select overall top HVGs
  top_genes <- hvg_per_cond %>%
    group_by(gene) %>%
    summarize(mean_rank = mean(rank), .groups = "drop") %>%
    arrange(mean_rank) %>%
    head(top_n * 2) %>%
    pull(gene)
  
  # Create dotplot data
  df_plot <- hvg_per_cond %>%
    filter(gene %in% top_genes)
  
  p <- ggplot(df_plot, aes(x = condition, y = gene, size = -rank, color = -rank)) +
    geom_point(alpha = 0.7) +
    scale_size_continuous(range = c(2, 8), name = "Rank\n(higher = better)") +
    scale_color_gradient(low = "lightblue", high = "darkred", name = "Rank\n(higher = better)") +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10)) +
    labs(title = "Top Variable Genes Per Condition",
         x = "Condition", y = "Gene")
  
  ggsave(file.path(plots_dir, "hvg_dotplot.png"), 
         p, width = 10, height = 12, dpi = 150)
}

#' Create per-condition HVG overview plot
#' @param sample_objects List of Seurat objects
#' @param plots_dir Directory to save plots
generate_hvg_per_condition_plot <- function(sample_objects, plots_dir) {
  # Count HVGs per condition
  hvg_counts <- lapply(names(sample_objects), function(sname) {
    seu <- sample_objects[[sname]]
    data.frame(
      sample = sname,
      condition = unique(seu$condition)[1],
      n_hvgs = length(VariableFeatures(seu)),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  p <- ggplot(hvg_counts, aes(x = condition, y = n_hvgs, fill = condition)) +
    geom_boxplot(outlier.size = 1) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    theme_minimal(base_size = 16) +
    labs(title = "Number of Highly Variable Genes Per Condition",
         y = "Number of HVGs", x = "Condition")
  
  ggsave(file.path(plots_dir, "hvg_selection_per_condition.png"), 
         p, width = 10, height = 6, dpi = 150)
}

#-------------------------------------------------------------------------------
# 4.8. Pre-Integration Dimensionality Reduction
#-------------------------------------------------------------------------------
# Create dimensionality reduction comparison plots (pre-integration)
# @param sample_objects List of Seurat objects before integration
# @param plots_dir Directory to save plots
generate_dimred_comparison_pre <- function(sample_objects, plots_dir) {
  # Merge samples without integration for comparison
  seu_merged <- merge(
    sample_objects[[1]], 
    y = sample_objects[-1],
    add.cell.ids = names(sample_objects),
    project = "merged"
  )
  
  # Run standard workflow
  seu_merged <- NormalizeData(seu_merged, verbose = FALSE)
  seu_merged <- FindVariableFeatures(seu_merged, nfeatures = 2000, verbose = FALSE)
  seu_merged <- ScaleData(seu_merged, verbose = FALSE)
  seu_merged <- RunPCA(seu_merged, npcs = 30, verbose = FALSE)
  seu_merged <- RunUMAP(seu_merged, dims = 1:20, verbose = FALSE)
  seu_merged <- RunTSNE(seu_merged, dims = 1:20, verbose = FALSE, check_duplicates = FALSE)
  
  # Plot PCA and UMAP side by side
  p_pca <- DimPlot(seu_merged, reduction = "pca", group.by = "condition") +
    labs(title = "Pre-Integration PCA - by Condition") +
    theme_minimal(base_size = 16)
  
  p_umap <- DimPlot(seu_merged, reduction = "umap", group.by = "condition") +
    labs(title = "Pre-Integration UMAP - by Condition") +
    theme_minimal(base_size = 16)
  
  p_combined <- p_pca + p_umap
  ggsave(file.path(plots_dir, "dimred_comparison_pre_integration.png"), 
         p_combined, width = 16, height = 6, dpi = 150)
  
  # Return the merged object for use in 3x3 grid plot
  return(invisible(seu_merged))
}

#' Create 3x3 integration comparison grid (PCA, UMAP, tSNE x Pre/Seurat/StacAS)
#' 
#' Generates a comprehensive 3x3 grid showing dimensionality reductions across
#' integration methods. Single shared legend to minimize redundancy.
#' 
#' @param pre_obj Pre-integration merged Seurat object
#' @param seurat_obj Seurat-integrated object (can be NULL)
#' @param stacas_obj StacAS-integrated object (can be NULL)
#' @param plots_dir Directory to save plots
generate_integration_3x3_grid <- function(pre_obj, seurat_obj, stacas_obj, plots_dir) {
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    if (exists("log", mode = "function")) {
      log("patchwork not available, skipping 3x3 grid plot")
    }
    return(invisible(NULL))
  }
  
  library(patchwork)
  
  # Load cowplot if available, otherwise use patchwork only
  use_cowplot <- requireNamespace("cowplot", quietly = TRUE)
  
  # Helper function to create plots without individual legends
  create_plot <- function(obj, reduction, method_name) {
    if (is.null(obj)) {
      return(ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = paste(method_name, "not available"), 
                 size = 5, color = "gray50") +
        theme_void())
    }
    
    DimPlot(obj, reduction = reduction, group.by = "condition") +
      labs(title = paste(method_name, "-", toupper(reduction))) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none",  # Remove individual legends
            plot.title = element_text(size = 11, face = "bold"))
  }
  
  # Create all 9 plots
  # Row 1: Pre-integration
  p_pre_pca <- create_plot(pre_obj, "pca", "Pre")
  p_pre_umap <- create_plot(pre_obj, "umap", "Pre")
  p_pre_tsne <- create_plot(pre_obj, "tsne", "Pre")
  
  # Row 2: Seurat integration
  p_seurat_pca <- create_plot(seurat_obj, "pca", "Seurat")
  p_seurat_umap <- create_plot(seurat_obj, "umap", "Seurat")
  p_seurat_tsne <- create_plot(seurat_obj, "tsne", "Seurat")
  
  # Row 3: StacAS integration
  p_stacas_pca <- create_plot(stacas_obj, "pca", "StacAS")
  p_stacas_umap <- create_plot(stacas_obj, "umap", "StacAS")
  p_stacas_tsne <- create_plot(stacas_obj, "tsne", "StacAS")
  
  # Combine into 3x3 grid
  layout <- (p_pre_pca | p_pre_umap | p_pre_tsne) /
            (p_seurat_pca | p_seurat_umap | p_seurat_tsne) /
            (p_stacas_pca | p_stacas_umap | p_stacas_tsne)
  
  # Try to add shared legend using cowplot if available
  if (use_cowplot) {
    # Extract legend from one plot (use pre_obj which should always exist)
    p_legend <- DimPlot(pre_obj, reduction = "pca", group.by = "condition") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "right",
            legend.title = element_text(size = 11, face = "bold"),
            legend.text = element_text(size = 10))
    
    # Extract just the legend
    legend <- cowplot::get_legend(p_legend)
    
    # Add legend to the right
    final_plot <- layout + 
      plot_annotation(
        title = "Integration Comparison: Pre-Integration vs Seurat vs StacAS",
        subtitle = "Dimensionality Reductions: PCA, UMAP, tSNE",
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                      plot.subtitle = element_text(size = 12, hjust = 0.5))
      )
    
    # Combine with legend using cowplot
    final_with_legend <- cowplot::plot_grid(
      final_plot, legend, 
      ncol = 2, 
      rel_widths = c(1, 0.15)
    )
    
    # Save the plot
    ggsave(file.path(plots_dir, "integration_comparison_3x3_grid.png"), 
           final_with_legend, width = 18, height = 15, dpi = 150, bg = "white")
  } else {
    # Fallback: use patchwork guide area for legend
    final_plot <- layout + 
      plot_annotation(
        title = "Integration Comparison: Pre-Integration vs Seurat vs StacAS",
        subtitle = "Dimensionality Reductions: PCA, UMAP, tSNE",
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                      plot.subtitle = element_text(size = 12, hjust = 0.5))
      ) +
      plot_layout(guides = "collect")  # Collect all guides into one
    
    # Save the plot
    ggsave(file.path(plots_dir, "integration_comparison_3x3_grid.png"), 
           final_plot, width = 18, height = 15, dpi = 150, bg = "white")
  }
  
  if (exists("log", mode = "function")) {
    log("Saved 3x3 integration comparison grid plot")
  }
  
  return(invisible(NULL))
}

#-------------------------------------------------------------------------------
# 5.1. Post-Integration Dimensionality Reduction (Seurat)
#-------------------------------------------------------------------------------
# Create dimensionality reduction comparison plots (Seurat-integrated)
# @param integrated_seu Seurat object after integration
# @param plots_dir Directory to save plots  
generate_dimred_comparison_seurat <- function(integrated_seu, plots_dir) {
  # Plot PCA and UMAP side by side
  p_pca <- DimPlot(integrated_seu, reduction = "pca", group.by = "condition") +
    labs(title = "Seurat-Integrated PCA - by Condition") +
    theme_minimal(base_size = 16)
  
  p_umap <- DimPlot(integrated_seu, reduction = "umap", group.by = "condition") +
    labs(title = "Seurat-Integrated UMAP - by Condition") +
    theme_minimal(base_size = 16)
  
  p_combined <- p_pca + p_umap
  ggsave(file.path(plots_dir, "dimred_comparison_seurat_integrated.png"), 
         p_combined, width = 16, height = 6, dpi = 150)
}

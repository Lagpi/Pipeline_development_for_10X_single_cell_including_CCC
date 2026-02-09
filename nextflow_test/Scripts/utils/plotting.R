# ==============================================================================
# Plotting Utilities
# ==============================================================================
# Shared plotting functions for saving figures in multiple formats
# 
# Functions:
# - save_plot(): Save ggplot as PNG and PDF
# - save_interactive(): Save plotly/htmlwidget as HTML
# - save_plot_multi(): Save plot in multiple formats
# ==============================================================================

#' Save ggplot to PNG and PDF with error handling
#' 
#' Saves plot in both PNG (for quick viewing) and PDF (for publication).
#' Requires plots_dir variable in calling environment.
#' 
#' @param p ggplot object
#' @param filename Output filename (will add .png/.pdf extensions)
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @param dpi Resolution for PNG in dots per inch (default: 300)
#' @param plots_dir Output directory (default: from parent environment)
#' @return Invisible path to PNG file
#' @examples
#' p <- ggplot(mtcars, aes(mpg, wt)) + geom_point()
#' save_plot(p, "scatter_plot.png", width = 8, height = 6)
save_plot <- function(p, filename, width = 10, height = 8, dpi = 300, 
                     plots_dir = NULL) {
  
  # Get plots_dir from parent environment if not provided
  if (is.null(plots_dir)) {
    if (exists("plots_dir", envir = parent.frame(), inherits = FALSE)) {
      plots_dir <- get("plots_dir", envir = parent.frame())
    } else {
      plots_dir <- "."
    }
  }
  
  # Ensure directory exists
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Ensure .png extension for base filename
  if (!grepl("\\.(png|pdf)$", filename)) {
    filename <- paste0(filename, ".png")
  }
  
  tryCatch({
    # Save PNG
    png_path <- file.path(plots_dir, filename)
    suppressMessages(
      ggplot2::ggsave(png_path, plot = p, width = width, height = height, 
                     dpi = dpi, bg = "white")
    )
    
    # Save PDF
    pdf_filename <- sub("\\.png$", ".pdf", filename)
    pdf_path <- file.path(plots_dir, pdf_filename)
    suppressMessages(
      ggplot2::ggsave(pdf_path, plot = p, width = width, height = height, 
                     device = "pdf")
    )
    
    if (exists("log", mode = "function")) {
      log("Saved plot:", filename, "and", pdf_filename)
    }
    
    return(invisible(png_path))
    
  }, error = function(e) {
    if (exists("log", mode = "function")) {
      log("ERROR saving plot", filename, ":", e$message)
    } else {
      warning("Error saving plot ", filename, ": ", e$message)
    }
    return(invisible(NULL))
  })
}


#' Save interactive plotly widget as self-contained HTML
#' 
#' Saves interactive plot (plotly, htmlwidget) as standalone HTML file.
#' Requires dash_dir variable in calling environment.
#' 
#' @param widget Interactive widget object (plotly, leaflet, DT, etc.)
#' @param filename Output HTML filename
#' @param selfcontained Create self-contained HTML (default: TRUE)
#' @param dash_dir Output directory (default: from parent environment)
#' @return Invisible path to HTML file
#' @examples
#' library(plotly)
#' p <- plot_ly(mtcars, x = ~mpg, y = ~wt)
#' save_interactive(p, "scatter_interactive.html")
save_interactive <- function(widget, filename, selfcontained = TRUE, 
                            dash_dir = NULL) {
  
  # Check if htmlwidgets is available
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    if (exists("log", mode = "function")) {
      log("WARNING: htmlwidgets not available, cannot save:", filename)
    }
    return(invisible(NULL))
  }
  
  # Get dash_dir from parent environment if not provided
  if (is.null(dash_dir)) {
    if (exists("dash_dir", envir = parent.frame(), inherits = FALSE)) {
      dash_dir <- get("dash_dir", envir = parent.frame())
    } else {
      dash_dir <- "."
    }
  }
  
  # Ensure directory exists
  if (!dir.exists(dash_dir)) {
    dir.create(dash_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Ensure .html extension
  if (!grepl("\\.html$", filename)) {
    filename <- paste0(tools::file_path_sans_ext(filename), ".html")
  }
  
  tryCatch({
    html_path <- file.path(dash_dir, filename)
    
    htmlwidgets::saveWidget(
      widget = widget,
      file = html_path,
      selfcontained = selfcontained,
      title = tools::file_path_sans_ext(filename)
    )
    
    if (exists("log", mode = "function")) {
      log("Saved interactive HTML:", filename)
    }
    
    return(invisible(html_path))
    
  }, error = function(e) {
    if (exists("log", mode = "function")) {
      log("ERROR saving interactive plot", filename, ":", e$message)
    } else {
      warning("Error saving interactive plot ", filename, ": ", e$message)
    }
    return(invisible(NULL))
  })
}


#' Save plot in multiple formats simultaneously
#' 
#' Convenience wrapper to save in PNG, PDF, and optionally SVG
#' 
#' @param p ggplot object
#' @param basename Base filename (without extension)
#' @param formats Vector of formats to save (default: c("png", "pdf"))
#' @param ... Additional arguments passed to ggsave
#' @return Invisible list of saved file paths
#' @examples
#' p <- ggplot(mtcars, aes(mpg, wt)) + geom_point()
#' save_plot_multi(p, "my_plot", formats = c("png", "pdf", "svg"))
save_plot_multi <- function(p, basename, formats = c("png", "pdf"), ...) {
  
  saved_files <- list()
  plots_dir <- if (exists("plots_dir", envir = parent.frame())) {
    get("plots_dir", envir = parent.frame())
  } else {
    "."
  }
  
  for (fmt in formats) {
    filename <- paste0(basename, ".", fmt)
    filepath <- file.path(plots_dir, filename)
    
    tryCatch({
      suppressMessages(
        ggplot2::ggsave(filepath, plot = p, device = fmt, ...)
      )
      saved_files[[fmt]] <- filepath
    }, error = function(e) {
      if (exists("log", mode = "function")) {
        log("WARNING: Could not save", fmt, "format:", e$message)
      }
    })
  }
  
  invisible(saved_files)
}


#' Create consistent theme for all plots
#' 
#' Returns a ggplot2 theme with publication-ready styling
#' 
#' @param base_size Base font size (default: 12)
#' @param base_family Font family (default: "")
#' @return ggplot2 theme object
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) + 
#'   geom_point() + 
#'   theme_publication()
theme_publication <- function(base_size = 12, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # Grid
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      panel.grid.minor = ggplot2::element_blank(),
      
      # Background
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      
      # Titles
      plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.2),
      plot.subtitle = ggplot2::element_text(color = "grey40"),
      
      # Axes
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "grey20"),
      axis.ticks = ggplot2::element_line(color = "grey80"),
      
      # Legend
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      
      # Strip (facets)
      strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

#' Save QC plot to both plots and dashboard directories
#' 
#' Specialized plotting function for QC module that saves plots to both
#' the main plots directory and copies to dashboard/interactive directory
#' for easy viewing. Saves both PNG and PDF formats.
#' 
#' @param p ggplot object
#' @param filename Output filename (will add .png/.pdf extensions)
#' @param width Plot width in inches (default: 9)
#' @param height Plot height in inches (default: 4)
#' @param dpi Resolution for PNG in dots per inch (default: 300)
#' @param plots_dir Main plots directory (default: from parent environment)
#' @param dash_dir Dashboard directory (default: from parent environment)
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
    
    plot_interactive <- tryCatch(
      plotly::ggplotly(p),
      error = function(e) {
        if (exists("log", mode = "function")) {
          log("Warning: Interactive plot failed:", e$message)
        }
        NULL
      }
    )
    
    if (!is.null(plot_interactive)) {
      htmlwidgets::saveWidget(plot_interactive, html_path, selfcontained = TRUE)
    }
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


# ==============================================================================
# Meta-Program Plotting Functions
# ==============================================================================

#' Create line plots of MP scores by condition
#' 
#' Creates line plots showing metaprogram scores across cell types,
#' with separate lines for each condition, all in one figure.
#' 
#' @param seu Seurat object with MP scores in metadata
#' @param mp.genes List of metaprogram genes (to get MP names)
#' @param condition_field Column name for condition
#' @param plots_dir Output directory for plots
#' @param log_fn Logging function
create_mp_lines_by_condition <- function(seu, mp.genes, condition_field, plots_dir, log_fn = message) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  
  log_fn("Creating MP score line plots by condition...")
  
  mp_names <- names(mp.genes)
  conditions <- unique(seu@meta.data[[condition_field]])
  
  # Get cell type annotation
  celltype_col <- if ("cell_type" %in% colnames(seu@meta.data)) {
    "cell_type"
  } else if ("singler_labels" %in% colnames(seu@meta.data)) {
    "singler_labels"
  } else {
    "seurat_clusters"
  }
  
  tryCatch({
    # Create data for plotting
    plot_data_list <- list()
    
    for (mp_name in mp_names) {
      if (!mp_name %in% colnames(seu@meta.data)) {
        log_fn("Warning: MP", mp_name, "not found in metadata")
        next
      }
      
      # Calculate mean score per cell type per condition
      df <- seu@meta.data[, c(mp_name, celltype_col, condition_field)]
      colnames(df) <- c("score", "celltype", "condition")
      
      agg_data <- df %>%
        dplyr::group_by(celltype, condition) %>%
        dplyr::summarise(
          mean_score = mean(score, na.rm = TRUE),
          se = sd(score, na.rm = TRUE) / sqrt(n()),
          .groups = "drop"
        ) %>%
        dplyr::mutate(mp = mp_name)
      
      plot_data_list[[mp_name]] <- agg_data
    }
    
    # Combine all data
    plot_data <- dplyr::bind_rows(plot_data_list)
    
    # Create line plot with facets for each MP
    p <- ggplot(plot_data, aes(x = celltype, y = mean_score, color = condition, group = condition)) +
      geom_line(size = 1) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se), width = 0.2) +
      facet_wrap(~mp, scales = "free_y", ncol = 2) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold", size = 10)
      ) +
      labs(
        x = "Cell Type",
        y = "Mean MP Score",
        color = "Condition",
        title = "Meta-Program Scores Across Cell Types by Condition"
      ) +
      scale_color_brewer(palette = "Set1")
    
    # Save plot
    save_plot(p, "mp_scores_lines_by_condition.png", 
              width = 14, height = max(8, ceiling(length(mp_names)/2) * 3), 
              plots_dir = plots_dir)
    
    # Also create individual plots for each MP
    for (mp_name in mp_names) {
      mp_data <- plot_data %>% dplyr::filter(mp == mp_name)
      
      p_single <- ggplot(mp_data, aes(x = celltype, y = mean_score, color = condition, group = condition)) +
        geom_line(size = 1.5) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se), width = 0.2, size = 1) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 14)
        ) +
        labs(
          x = "Cell Type",
          y = "Mean MP Score",
          color = "Condition",
          title = paste("Meta-Program:", mp_name)
        ) +
        scale_color_brewer(palette = "Set1")
      
      safe_name <- gsub("[^[:alnum:]_]", "_", mp_name)
      save_plot(p_single, paste0("mp_scores_lines_", safe_name, ".png"),
                width = 10, height = 6, plots_dir = plots_dir)
    }
    
    log_fn("Saved MP line plots")
    
  }, error = function(e) {
    log_fn("ERROR creating MP line plots:", e$message)
  })
  
  invisible(NULL)
}


#' Create UMAP plots colored by MP scores
#' 
#' Creates UMAP visualizations showing metaprogram activity across cells.
#' Generates both a grid of all MPs and individual plots with MP gene lists.
#' 
#' @param seu Seurat object with UMAP reduction and MP scores
#' @param mp.genes List of metaprogram genes
#' @param plots_dir Output directory for plots
#' @param log_fn Logging function
create_mp_umap_plots <- function(seu, mp.genes, plots_dir, log_fn = message) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (!requireNamespace("Seurat", quietly = TRUE)) return(invisible(NULL))
  
  log_fn("Creating MP UMAP visualizations...")
  
  mp_names <- names(mp.genes)
  
  # Check if UMAP exists
  if (!"umap" %in% names(seu@reductions)) {
    log_fn("Warning: No UMAP reduction found, skipping MP UMAP plots")
    return(invisible(NULL))
  }
  
  tryCatch({
    # Create individual UMAP plots for each MP with gene names
    for (mp_name in mp_names) {
      if (!mp_name %in% colnames(seu@meta.data)) {
        log_fn("Warning: MP", mp_name, "not found in metadata")
        next
      }
      
      # Create UMAP plot
      p <- FeaturePlot(seu, features = mp_name, reduction = "umap") +
        viridis::scale_color_viridis(option = "magma") +
        theme_minimal() +
        labs(title = paste("Meta-Program:", mp_name))
      
      # Get top genes for this MP
      top_genes <- head(mp.genes[[mp_name]], 10)
      gene_text <- paste(top_genes, collapse = ", ")
      
      # Add gene names as subtitle
      p <- p + labs(subtitle = paste("Top genes:", gene_text))
      
      safe_name <- gsub("[^[:alnum:]_]", "_", mp_name)
      save_plot(p, paste0("mp_umap_", safe_name, ".png"),
                width = 8, height = 7, plots_dir = plots_dir)
    }
    
    # Create grid of all MPs
    plot_list <- list()
    for (mp_name in mp_names) {
      if (mp_name %in% colnames(seu@meta.data)) {
        p <- FeaturePlot(seu, features = mp_name, reduction = "umap") +
          viridis::scale_color_viridis(option = "magma") +
          theme_minimal() +
          theme(
            legend.position = "right",
            plot.title = element_text(size = 10, face = "bold")
          ) +
          labs(title = mp_name)
        
        plot_list[[mp_name]] <- p
      }
    }
    
    if (length(plot_list) > 0) {
      # Combine plots
      combined <- patchwork::wrap_plots(plot_list, ncol = 2)
      save_plot(combined, "mp_umap_grid.png",
                width = 12, height = max(8, ceiling(length(plot_list)/2) * 4),
                plots_dir = plots_dir)
    }
    
    log_fn("Saved MP UMAP plots")
    
  }, error = function(e) {
    log_fn("ERROR creating MP UMAP plots:", e$message)
  })
  
  invisible(NULL)
}


#' Print MP names and gene composition
#' 
#' Prints metaprogram names and their top genes to console and log file.
#' 
#' @param mp.genes List of metaprogram genes
#' @param log_fn Logging function
#' @param n_genes Number of top genes to display per MP (default: 15)
print_mp_info <- function(mp.genes, log_fn = message, n_genes = 15) {
  log_fn("\n=== META-PROGRAM INFORMATION ===\n")
  log_fn("Discovered", length(mp.genes), "meta-programs:\n")
  
  for (mp_name in names(mp.genes)) {
    genes <- mp.genes[[mp_name]]
    log_fn("\n", mp_name, " (", length(genes), " genes)")
    log_fn("  Top genes:", paste(head(genes, n_genes), collapse = ", "))
  }
  
  log_fn("\n================================\n")
  
  invisible(NULL)
}


# ==============================================================================
# Publication-Ready Theme
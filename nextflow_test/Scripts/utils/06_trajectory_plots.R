#===============================================================================
# Trajectory Analysis Plotting Functions
# Module 06: Trajectory Inference
#===============================================================================

# Configuration parameters for plots (read from parent environment or use defaults)
CFG_BASE_SIZE <- if (exists("CFG_BASE_SIZE", envir = parent.frame(), inherits = TRUE)) {
  get("CFG_BASE_SIZE", envir = parent.frame(), inherits = TRUE)
} else { 18 }
CFG_TITLE_SIZE <- if (exists("CFG_TITLE_SIZE", envir = parent.frame(), inherits = TRUE)) {
  get("CFG_TITLE_SIZE", envir = parent.frame(), inherits = TRUE)
} else { 24 }
CFG_AXIS_TITLE_SIZE <- if (exists("CFG_AXIS_TITLE_SIZE", envir = parent.frame(), inherits = TRUE)) {
  get("CFG_AXIS_TITLE_SIZE", envir = parent.frame(), inherits = TRUE)
} else { 18 }
CFG_PLOT_WIDTH <- 12
CFG_PLOT_HEIGHT <- 8
CFG_PLOT_DPI <- 150

#' Create projection demonstration plot
#' Shows original trajectory vs. projected cells
create_projection_demo_plot <- function(rd, rd_name, pto_original, newPCA, newPTO, pseudotime_curves, plots_dir, log) {
  tryCatch({
    colors <- colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral')[-6])(100)
    
    pdf(file.path(plots_dir, "slingshot_projection_demo.pdf"), width = 12, height = 6)
    par(mfrow = c(1, 2))
    
    # Original cells
    plotcol_orig <- colors[cut(pseudotime_curves[, 1], breaks = 100)]
    plot(rd, col = plotcol_orig, pch = 16, asp = 1, 
         main = "Original Trajectory", xlab = paste0(rd_name, "_1"), ylab = paste0(rd_name, "_2"))
    lines(pto_original, lwd = 2, col = 'black')
    
    # Projected cells
    plotcol_new <- colors[cut(slingshot::slingPseudotime(newPTO)[, 1], breaks = 100)]
    plot(rd, col = 'grey80', pch = 16, asp = 1,
         main = "Projected Cells onto Trajectory", 
         xlab = paste0(rd_name, "_1"), ylab = paste0(rd_name, "_2"),
         xlim = range(c(rd[, 1], newPCA[, 1])), ylim = range(c(rd[, 2], newPCA[, 2])))
    lines(pto_original, lwd = 2, col = 'black')
    points(newPCA, col = plotcol_new, pch = 16)
    legend("topright", legend = c("Original", "Projected"), pch = 16, col = c("grey80", "red"))
    
    dev.off()
    
    log("  - Cell projection demonstration completed")
  }, error = function(e) {
    log("  - Warning: Cell projection demo failed:", conditionMessage(e))
  })
}

#' Create lineages plot using getLineages
create_lineages_plot <- function(dimred, lineages, clustering, plots_dir, log, condition = NULL) {
  tryCatch({
    safe_cond <- if (!is.null(condition) && nzchar(condition)) {
      gsub("[^A-Za-z0-9_-]+", "_", condition)
    } else {
      ""
    }
    
    file_name <- if (nzchar(safe_cond)) {
      sprintf("lineages_getLineages_%s.png", safe_cond)
    } else {
      "lineages_getLineages.png"
    }
    
    title_str <- if (!is.null(condition) && nzchar(condition)) {
      sprintf("Lineages (getLineages) - %s", condition)
    } else {
      "Lineages (getLineages)"
    }
    
    png(file.path(plots_dir, file_name), width = 1600, height = 1200)
    par(mar = c(4,4,2,2))
    cols <- RColorBrewer::brewer.pal(max(3, min(9, length(unique(clustering)))), "Set1")
    plot(dimred, col = cols[as.factor(clustering)], asp = 1, pch = 16, main = title_str)
    
    # Add cell type/cluster labels
    tryCatch({
      unique_clusters <- unique(clustering)
      for (clust in unique_clusters) {
        idx <- which(clustering == clust)
        if (length(idx) > 0) {
          center_x <- median(dimred[idx, 1], na.rm = TRUE)
          center_y <- median(dimred[idx, 2], na.rm = TRUE)
          text(center_x, center_y, labels = as.character(clust), 
               cex = 1.5, font = 2, col = "black")
        }
      }
    }, error = function(e2) {
      log("  - Warning: Could not add cluster labels:", conditionMessage(e2))
    })
    
    # Plot lineages manually - extract lineage paths from S4 object
    tryCatch({
      if (inherits(lineages, "SlingshotDataSet") || inherits(lineages, "PseudotimeOrdering")) {
        # Use slingshot plotting method
        slingshot::lines(lineages, lwd = 3, col = 'black', type = 'lineages')
      } else if (is.list(lineages) && "lineages" %in% names(lineages)) {
        # Extract lineage information if available
        for (lin in lineages$lineages) {
          lines(dimred[lin, ], lwd = 3, col = 'black')
        }
      }
    }, error = function(e2) {
      # If plotting lineages fails, just save the scatter plot
      log("  - Warning: Could not plot lineages, saved scatter plot only:", conditionMessage(e2))
    })
    
    dev.off()
    
    log(sprintf("Saved getLineages object and plot for %s", if(nzchar(safe_cond)) condition else "all cells"))
  }, error = function(e) {
    log("getLineages plotting failed: ", conditionMessage(e))
  })
}

#' Create slingshot trajectory overview plot
create_slingshot_trajectory_overview <- function(viz_data, coord_cols, curves, plots_dir, log) {
  tryCatch({
    png(file.path(plots_dir, "slingshot_trajectory_overview.png"), width = 1800, height = 1400)
    par(mfrow = c(3, 3))
    
    # Determine what to use for coloring - prefer cell_type over cluster numbers
    color_column <- if ("cell_type" %in% colnames(viz_data)) "cell_type" else "cluster"
    unique_colors <- unique(viz_data[[color_column]])
    n_colors <- length(unique_colors)
    colors_map <- rainbow(n_colors)
    color_factor <- as.factor(viz_data[[color_column]])
    color_numeric <- as.numeric(color_factor)
    
    # Plot trajectory curves by cluster/cell_type
    plot(viz_data[[coord_cols[1]]], viz_data[[coord_cols[2]]],
         col = colors_map[color_numeric],
         pch = 16, cex = 0.5, main = "Trajectory Curves by Cell Type")
    
    # Add curves if available
    if (length(curves) > 0) {
      for (i in seq_along(curves)) {
        curve <- curves[[i]]
        lines(curve$s, lwd = 3, col = "black")
      }
    }
    legend("topright", legend = unique_colors, 
           col = colors_map, pch = 16, cex = 0.8)
    
    # Pseudotime plots for each curve
    pseudotime_cols <- colnames(viz_data)[grepl("^curve", colnames(viz_data)) & !grepl("weight", colnames(viz_data))]
    for (i in 1:min(6, length(pseudotime_cols))) {
      pt_col <- pseudotime_cols[i]
      plot(viz_data[[coord_cols[1]]], viz_data[[coord_cols[2]]],
           col = colorRampPalette(c("navy", "cyan", "yellow", "red"))(100)[cut(viz_data[[pt_col]], 100)],
           pch = 16, cex = 0.5, main = paste("Pseudotime -", pt_col))
      
      # Add curve
      if (i <= length(curves)) {
        curve <- curves[[i]]
        lines(curve$s, lwd = 3, col = "black")
      }
    }
    
    # Cell type if available
    if ("cell_type" %in% colnames(viz_data)) {
      plot(viz_data[[coord_cols[1]]], viz_data[[coord_cols[2]]],
           col = rainbow(length(unique(viz_data$cell_type)))[as.factor(viz_data$cell_type)],
           pch = 16, cex = 0.5, main = "Trajectory by Cell Type")
      legend("topright", legend = unique(viz_data$cell_type), 
             col = rainbow(length(unique(viz_data$cell_type))), pch = 16, cex = 0.7)
    }
    dev.off()
  }, error = function(e) {
    log("Slingshot overview plot failed:", conditionMessage(e))
  })
}

#' Create generic trajectory plot
create_generic_trajectory_plot <- function(viz_data, coord_cols, plots_dir, log) {
  tryCatch({
    png(file.path(plots_dir, "trajectory_generic_overview.png"), width = 1600, height = 1200)
    par(mfrow = c(2, 2))
    
    # By cluster
    plot(viz_data[[coord_cols[1]]], viz_data[[coord_cols[2]]],
         col = rainbow(length(unique(viz_data$cluster)))[as.factor(viz_data$cluster)],
         pch = 16, cex = 0.5, main = "Cells by Cluster")
    
    # By pseudotime
    if ("pseudotime" %in% colnames(viz_data)) {
      plot(viz_data[[coord_cols[1]]], viz_data[[coord_cols[2]]],
           col = colorRampPalette(c("blue", "red"))(100)[cut(viz_data$pseudotime, 100)],
           pch = 16, cex = 0.5, main = "Pseudotime")
    }
    
    # By cell type if available
    if ("cell_type" %in% colnames(viz_data)) {
      plot(viz_data[[coord_cols[1]]], viz_data[[coord_cols[2]]],
           col = rainbow(length(unique(viz_data$cell_type)))[as.factor(viz_data$cell_type)],
           pch = 16, cex = 0.5, main = "Cells by Type")
    }
    dev.off()
  }, error = function(e) {
    log("Generic trajectory plot failed:", conditionMessage(e))
  })
}

#' Create pseudotime distribution histograms
create_pseudotime_histograms <- function(viz_data, plots_dir, log) {
  tryCatch({
    pseudotime_cols <- colnames(viz_data)[grepl("^curve", colnames(viz_data)) & !grepl("weight", colnames(viz_data))]
    n_curves <- length(pseudotime_cols)
    png(file.path(plots_dir, "pseudotime_histograms.png"), width = 1600, height = max(800, 400 * ceiling(n_curves/2)))
    par(mfrow = c(ceiling(n_curves/2), 2))
    
    for (pt_col in pseudotime_cols) {
      hist(viz_data[[pt_col]], breaks = 50, main = paste("Distribution -", pt_col),
           xlab = "Pseudotime", col = "steelblue", border = "white")
      abline(v = median(viz_data[[pt_col]], na.rm = TRUE), col = "red", lwd = 2, lty = 2)
    }
    dev.off()
  }, error = function(e) {
    log("Pseudotime histogram failed:", conditionMessage(e))
  })
}

#' Create pseudotime by cell type boxplots
create_pseudotime_boxplots <- function(viz_data, plots_dir, log) {
  tryCatch({
    if (!"cell_type" %in% colnames(viz_data)) return()
    
    pseudotime_cols <- colnames(viz_data)[grepl("^curve", colnames(viz_data)) & !grepl("weight", colnames(viz_data))]
    n_curves <- length(pseudotime_cols)
    png(file.path(plots_dir, "pseudotime_boxplots_by_celltype.png"), width = 1800, height = max(900, 450 * ceiling(n_curves/2)))
    par(mfrow = c(ceiling(n_curves/2), 2))
    
    for (pt_col in pseudotime_cols) {
      boxplot(viz_data[[pt_col]] ~ viz_data$cell_type,
              main = pt_col, xlab = "", ylab = "Pseudotime",
              las = 2, col = rainbow(length(unique(viz_data$cell_type))),
              cex.axis = 0.8)
    }
    dev.off()
  }, error = function(e) {
    log("Pseudotime boxplot failed:", conditionMessage(e))
  })
}

#' Create curve weights heatmap
create_curve_weights_heatmap <- function(viz_data, plots_dir, log) {
  tryCatch({
    if (!any(grepl("weight", colnames(viz_data)))) return()
    
    weight_cols <- colnames(viz_data)[grepl("weight", colnames(viz_data))]
    if (length(weight_cols) <= 1) return()
    
    weight_matrix <- as.matrix(viz_data[, weight_cols])
    
    # Sample cells if too many
    if (nrow(weight_matrix) > 1000) {
      sample_idx <- sample(1:nrow(weight_matrix), 1000)
      weight_matrix <- weight_matrix[sample_idx, ]
    }
    
    pheatmap::pheatmap(
      t(weight_matrix),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      show_colnames = FALSE,
      color = viridis::viridis(100),
      main = "Trajectory Curve Weights per Cell",
      filename = file.path(plots_dir, "trajectory_curve_weights_heatmap.png")
    )
  }, error = function(e) {
    log("Curve weights heatmap failed:", conditionMessage(e))
  })
}

#' Create minimal slingshot lineages plot
create_minimal_slingshot_plot <- function(dimred, curves, clustering, plots_dir, log) {
  tryCatch({
    png(file.path(plots_dir, "lineages_slingshot.png"), width = 1400, height = 1000)
    par(mar = c(4,4,2,2))
    cols <- RColorBrewer::brewer.pal(max(3, min(9, length(unique(clustering)))), "Set1")
    plot(dimred, col = cols[as.factor(clustering)], asp = 1, pch = 16, main = "Slingshot Lineages")
    lines(curves, lwd = 3, col = 'black')
    dev.off()
    
    log("Minimal Slingshot: saved lineages, curves, and plot")
  }, error = function(e) {
    log("Minimal Slingshot plotting failed:", conditionMessage(e))
  })
}

#' Create interactive pseudotime by condition plot
create_interactive_pseudotime_by_condition <- function(viz_data, cond_field, pt_col, dash_dir, log) {
  tryCatch({
    if (!pt_col %in% colnames(viz_data)) return()
    
    pt_condition_data <- viz_data[!is.na(viz_data[[pt_col]]), c(cond_field, pt_col, "cell_type")]
    colnames(pt_condition_data) <- c("condition", "pseudotime", "cell_type")
    
    p_pt_cond <- plotly::plot_ly(
      data = pt_condition_data,
      x = ~condition,
      y = ~pseudotime,
      color = ~condition,
      type = "box",
      text = ~paste("Condition:", condition, "<br>Pseudotime:", round(pseudotime, 3),
                   "<br>Cell Type:", cell_type),
      hovertemplate = "%{text}<extra></extra>"
    ) %>%
    plotly::layout(
      title = paste("Pseudotime Distribution by Condition -", pt_col),
      xaxis = list(title = "Condition"),
      yaxis = list(title = "Pseudotime"),
      showlegend = TRUE
    )
    save_interactive(p_pt_cond, file.path(dash_dir, paste0("pseudotime_", pt_col, "_by_condition.html")))
  }, error = function(e) {
    log("Interactive pseudotime by condition failed:", conditionMessage(e))
  })
}

#' Create interactive trajectory by cell type plot
create_interactive_celltype_plot <- function(viz_data, coord_cols, dash_dir, log) {
  tryCatch({
    if (!"cell_type" %in% colnames(viz_data)) return()
    
    p_celltype <- plotly::plot_ly(
      data = viz_data,
      x = ~get(coord_cols[1]),
      y = ~get(coord_cols[2]),
      color = ~as.factor(cell_type),
      type = "scatter",
      mode = "markers",
      marker = list(size = 5),
      text = ~paste("Cell:", cell_id, "<br>Cluster:", cluster, "<br>Cell Type:", cell_type),
      hovertemplate = "%{text}<extra></extra>"
    ) %>%
    plotly::layout(
      title = "Trajectory Analysis - Cell Types",
      xaxis = list(title = coord_cols[1]),
      yaxis = list(title = coord_cols[2])
    )
    save_interactive(p_celltype, file.path(dash_dir, "trajectory_celltypes.html"))
  }, error = function(e) {
    log("Interactive celltype plot failed:", conditionMessage(e))
  })
}

#' Create interactive pseudotime visualization (slingshot)
create_interactive_pseudotime_slingshot <- function(viz_data, coord_cols, dash_dir, log) {
  tryCatch({
    pseudotime_cols <- colnames(viz_data)[grepl("^curve", colnames(viz_data)) & !grepl("weight", colnames(viz_data))]
    for (i in 1:min(5, length(pseudotime_cols))) {
      pt_col <- pseudotime_cols[i]
      
      p_pseudo <- plotly::plot_ly(
        data = viz_data,
        x = ~get(coord_cols[1]),
        y = ~get(coord_cols[2]),
        color = ~get(pt_col),
        type = "scatter",
        mode = "markers",
        marker = list(size = 5, colorscale = "Viridis"),
        text = ~paste("Cell:", cell_id, "<br>Pseudotime:", round(get(pt_col), 3), "<br>Cluster:", cluster, "<br>Cell Type:", cell_type),
        hovertemplate = "%{text}<extra></extra>"
      ) %>%
      plotly::layout(
        title = paste("Pseudotime -", pt_col),
        xaxis = list(title = coord_cols[1]),
        yaxis = list(title = coord_cols[2])
      )
      save_interactive(p_pseudo, file.path(dash_dir, paste0("pseudotime_", pt_col, ".html")))
    }
  }, error = function(e) {
    log("Interactive pseudotime slingshot failed:", conditionMessage(e))
  })
}

#' Create interactive pseudotime curve comparison
create_interactive_curve_comparison <- function(viz_data, dash_dir, log) {
  tryCatch({
    pseudotime_cols <- colnames(viz_data)[grepl("^curve", colnames(viz_data)) & !grepl("weight", colnames(viz_data))]
    if (length(pseudotime_cols) < 2) return()
    
    p_compare <- plotly::plot_ly(
      data = viz_data,
      x = ~get(pseudotime_cols[1]),
      y = ~get(pseudotime_cols[2]),
      color = ~cell_type,
      type = "scatter",
      mode = "markers",
      marker = list(size = 5),
      text = ~paste("Cell Type:", cell_type, 
                   "<br>", pseudotime_cols[1], ":", round(get(pseudotime_cols[1]), 3),
                   "<br>", pseudotime_cols[2], ":", round(get(pseudotime_cols[2]), 3)),
      hovertemplate = "%{text}<extra></extra>"
    ) %>%
    plotly::layout(
      title = "Pseudotime Comparison Between Curves",
      xaxis = list(title = pseudotime_cols[1]),
      yaxis = list(title = pseudotime_cols[2])
    )
    save_interactive(p_compare, file.path(dash_dir, "pseudotime_curve_comparison.html"))
  }, error = function(e) {
    log("Interactive curve comparison failed:", conditionMessage(e))
  })
}

#' Create interactive curve weights visualization
create_interactive_curve_weights <- function(viz_data, coord_cols, dash_dir, log) {
  tryCatch({
    weight_cols <- colnames(viz_data)[grepl("weight", colnames(viz_data))]
    if (length(weight_cols) == 0) return()
    
    for (wt_col in weight_cols[1:min(3, length(weight_cols))]) {
      p_weight <- plotly::plot_ly(
        data = viz_data,
        x = ~get(coord_cols[1]),
        y = ~get(coord_cols[2]),
        color = ~get(wt_col),
        type = "scatter",
        mode = "markers",
        marker = list(size = 5, colorscale = "Blues"),
        text = ~paste("Cell:", cell_id, "<br>Weight:", round(get(wt_col), 3), "<br>Cluster:", cluster),
        hovertemplate = "%{text}<extra></extra>"
      ) %>%
      plotly::layout(
        title = paste("Trajectory", wt_col),
        xaxis = list(title = coord_cols[1]),
        yaxis = list(title = coord_cols[2])
      )
      save_interactive(p_weight, file.path(dash_dir, paste0("trajectory_", wt_col, ".html")))
    }
  }, error = function(e) {
    log("Interactive curve weights failed:", conditionMessage(e))
  })
}

#' Create interactive pseudotime visualization (generic)
create_interactive_pseudotime_generic <- function(viz_data, coord_cols, dash_dir, log) {
  tryCatch({
    if (!"pseudotime" %in% colnames(viz_data)) return()
    
    p_pseudo <- plotly::plot_ly(
      data = viz_data,
      x = ~get(coord_cols[1]),
      y = ~get(coord_cols[2]),
      color = ~pseudotime,
      type = "scatter",
      mode = "markers",
      marker = list(size = 5, colorscale = "Plasma"),
      text = ~paste("Cell:", cell_id, "<br>Pseudotime:", round(pseudotime, 3), "<br>Cluster:", cluster),
      hovertemplate = "%{text}<extra></extra>"
    ) %>%
    plotly::layout(
      title = "Pseudotime Analysis",
      xaxis = list(title = coord_cols[1]),
      yaxis = list(title = coord_cols[2])
    )
    save_interactive(p_pseudo, file.path(dash_dir, "pseudotime_analysis.html"))
  }, error = function(e) {
    log("Interactive pseudotime generic failed:", conditionMessage(e))
  })
}

#' Create interactive violin plot by cell type
create_interactive_violin_by_celltype <- function(viz_data, trajectory_method_used, dash_dir, log) {
  tryCatch({
    if (!"cell_type" %in% colnames(viz_data)) return()
    
    if (trajectory_method_used == "slingshot") {
      pseudotime_cols <- colnames(viz_data)[grepl("^curve", colnames(viz_data)) & !grepl("weight", colnames(viz_data))]
      
      for (pt_col in pseudotime_cols[1:min(3, length(pseudotime_cols))]) {
        p_violin <- plotly::plot_ly(
          data = viz_data,
          x = ~cell_type,
          y = ~get(pt_col),
          color = ~cell_type,
          type = "violin",
          box = list(visible = TRUE),
          meanline = list(visible = TRUE)
        ) %>%
        plotly::layout(
          title = paste("Pseudotime Distribution by Cell Type -", pt_col),
          xaxis = list(title = "Cell Type", tickangle = 45),
          yaxis = list(title = "Pseudotime")
        )
        save_interactive(p_violin, file.path(dash_dir, paste0("violin_", pt_col, "_by_celltype.html")))
      }
    } else if ("pseudotime" %in% colnames(viz_data)) {
      p_violin <- plotly::plot_ly(
        data = viz_data,
        x = ~cell_type,
        y = ~pseudotime,
        color = ~cell_type,
        type = "violin",
        box = list(visible = TRUE),
        meanline = list(visible = TRUE)
      ) %>%
      plotly::layout(
        title = "Pseudotime Distribution by Cell Type",
        xaxis = list(title = "Cell Type", tickangle = 45),
        yaxis = list(title = "Pseudotime")
      )
      save_interactive(p_violin, file.path(dash_dir, "violin_pseudotime_by_celltype.html"))
    }
  }, error = function(e) {
    log("Interactive violin by celltype failed:", conditionMessage(e))
  })
}

#' Create interactive 3D trajectory visualization
create_interactive_3d_trajectory <- function(viz_data, coord_cols, dash_dir, log) {
  tryCatch({
    if (length(coord_cols) < 2) return()
    
    pseudotime_cols <- colnames(viz_data)[grepl("^curve", colnames(viz_data)) & !grepl("weight", colnames(viz_data))]
    if (length(pseudotime_cols) < 1) return()
    
    p_3d <- plotly::plot_ly(
      data = viz_data,
      x = ~get(coord_cols[1]),
      y = ~get(coord_cols[2]),
      z = ~get(pseudotime_cols[1]),
      color = ~cell_type,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 3),
      text = ~paste("Cell Type:", cell_type,
                   "<br>", coord_cols[1], ":", round(get(coord_cols[1]), 3),
                   "<br>", coord_cols[2], ":", round(get(coord_cols[2]), 3),
                   "<br>Pseudotime:", round(get(pseudotime_cols[1]), 3)),
      hovertemplate = "%{text}<extra></extra>"
    ) %>%
    plotly::layout(
      title = "3D Trajectory Visualization",
      scene = list(
        xaxis = list(title = coord_cols[1]),
        yaxis = list(title = coord_cols[2]),
        zaxis = list(title = "Pseudotime")
      )
    )
    save_interactive(p_3d, file.path(dash_dir, "trajectory_3d.html"))
  }, error = function(e) {
    log("Interactive 3D trajectory failed:", conditionMessage(e))
  })
}

#' Create pseudotime distribution plots for all conditions combined
#' Shows density distributions of pseudotime for all lineages and all conditions in one figure
create_pseudotime_by_condition_combined <- function(viz_data, condition_field, plots_dir, log = cat) {
  tryCatch({
    if (!condition_field %in% colnames(viz_data)) {
      log("Warning: condition field not found in visualization data")
      return()
    }
    
    library(ggplot2)
    library(patchwork)
    
    # Get pseudotime columns (curves)
    pt_cols <- colnames(viz_data)[grepl("^curve", colnames(viz_data)) & !grepl("weight", colnames(viz_data))]
    
    if (length(pt_cols) == 0) {
      log("Warning: No pseudotime curves found in data")
      return()
    }
    
    conditions <- unique(viz_data[[condition_field]])
    conditions <- conditions[!is.na(conditions)]
    
    # Color palette for conditions
    if (length(conditions) <= 4) {
      cond_colors <- c("#2ecc71", "#3498db", "#e74c3c", "#f39c12")[1:length(conditions)]
    } else {
      cond_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(conditions))
    }
    
    names(cond_colors) <- conditions
    
    # Create a plot for each lineage
    plots <- list()
    
    for (i in seq_along(pt_cols)) {
      pt_col <- pt_cols[i]
      lineage_name <- gsub("^curve", "Lineage ", pt_col)
      
      # Prepare data for this lineage
      plot_data <- viz_data %>%
        dplyr::filter(!is.na(!!sym(pt_col))) %>%
        dplyr::select(all_of(c(pt_col, condition_field))) %>%
        dplyr::rename(pseudotime = !!pt_col, condition = !!condition_field)
      
      if (nrow(plot_data) == 0) next
      
      # Create density plot with all conditions
      p <- ggplot(plot_data, aes(x = pseudotime, fill = condition, color = condition)) +
        geom_density(alpha = 0.5, size = 1) +
        scale_fill_manual(values = cond_colors, name = "Condition") +
        scale_color_manual(values = cond_colors, name = "Condition") +
        labs(
          title = lineage_name,
          x = "Pseudotime",
          y = "Density"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 14, face = "bold"),
          legend.position = "right",
          panel.grid.minor = element_blank()
        )
      
      plots[[i]] <- p
    }
    
    if (length(plots) == 0) {
      log("Warning: No lineages could be plotted")
      return()
    }
    
    # Combine all plots
    n_lineages <- length(plots)
    n_cols <- min(2, n_lineages)
    n_rows <- ceiling(n_lineages / n_cols)
    
    combined <- patchwork::wrap_plots(plots, ncol = n_cols) +
      patchwork::plot_annotation(
        title = "Pseudotime Distributions by Condition (All Lineages)",
        subtitle = paste("Comparing", length(conditions), "conditions across", n_lineages, "lineages"),
        theme = theme(
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 14, hjust = 0.5)
        )
      )
    
    # Save as PNG with appropriate size
    width <- max(10, n_cols * 6)
    height <- max(8, n_rows * 6)
    
    ggsave(
      file.path(plots_dir, "pseudotime_density_all_conditions.png"),
      combined,
      width = width,
      height = height,
      dpi = 150,
      limitsize = FALSE
    )
    
    log("Saved combined pseudotime distributions for all conditions")
    
  }, error = function(e) {
    log("Warning: Pseudotime condition combined plot failed:", conditionMessage(e))
  })
}

#' Create histogram plots for pseudotime by condition
#' Shows overlaid histograms for all conditions
create_pseudotime_histograms_by_condition <- function(viz_data, condition_field, plots_dir, log = cat) {
  tryCatch({
    if (!condition_field %in% colnames(viz_data)) {
      log("Warning: condition field not found in visualization data")
      return()
    }
    
    library(ggplot2)
    library(patchwork)
    
    # Get pseudotime columns (curves)
    pt_cols <- colnames(viz_data)[grepl("^curve", colnames(viz_data)) & !grepl("weight", colnames(viz_data))]
    
    if (length(pt_cols) == 0) {
      log("Warning: No pseudotime curves found in data")
      return()
    }
    
    conditions <- unique(viz_data[[condition_field]])
    conditions <- conditions[!is.na(conditions)]
    
    # Color palette for conditions
    if (length(conditions) <= 4) {
      cond_colors <- c("#2ecc71", "#3498db", "#e74c3c", "#f39c12")[1:length(conditions)]
    } else {
      cond_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(conditions))
    }
    
    names(cond_colors) <- conditions
    
    # Create a histogram plot for each lineage
    plots <- list()
    
    for (i in seq_along(pt_cols)) {
      pt_col <- pt_cols[i]
      lineage_name <- gsub("^curve", "Lineage ", pt_col)
      
      # Prepare data for this lineage
      plot_data <- viz_data %>%
        dplyr::filter(!is.na(!!sym(pt_col))) %>%
        dplyr::select(all_of(c(pt_col, condition_field))) %>%
        dplyr::rename(pseudotime = !!pt_col, condition = !!condition_field)
      
      if (nrow(plot_data) == 0) next
      
      # Create histogram with all conditions
      p <- ggplot(plot_data, aes(x = pseudotime, fill = condition)) +
        geom_histogram(alpha = 0.6, bins = 30, position = "identity") +
        scale_fill_manual(values = cond_colors, name = "Condition") +
        labs(
          title = lineage_name,
          x = "Pseudotime",
          y = "Count"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 14, face = "bold"),
          legend.position = "right",
          panel.grid.minor = element_blank()
        )
      
      plots[[i]] <- p
    }
    
    if (length(plots) == 0) {
      log("Warning: No lineages could be plotted")
      return()
    }
    
    # Combine all plots
    n_lineages <- length(plots)
    n_cols <- min(2, n_lineages)
    n_rows <- ceiling(n_lineages / n_cols)
    
    combined <- patchwork::wrap_plots(plots, ncol = n_cols) +
      patchwork::plot_annotation(
        title = "Pseudotime Histograms by Condition (All Lineages)",
        subtitle = paste("Comparing", length(conditions), "conditions across", n_lineages, "lineages"),
        theme = theme(
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 14, hjust = 0.5)
        )
      )
    
    # Save as PNG
    width <- max(10, n_cols * 6)
    height <- max(8, n_rows * 6)
    
    ggsave(
      file.path(plots_dir, "pseudotime_histograms_all_conditions.png"),
      combined,
      width = width,
      height = height,
      dpi = 150,
      limitsize = FALSE
    )
    
    log("Saved pseudotime histograms for all conditions")
    
  }, error = function(e) {
    log("Warning: Pseudotime histograms plot failed:", conditionMessage(e))
  })
}



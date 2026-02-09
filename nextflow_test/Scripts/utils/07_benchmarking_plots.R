#===============================================================================
# 07_benchmarking/plots.R
# Benchmarking and performance plotting functions
#===============================================================================

# Plot functions specific to benchmarking analysis
# - Runtime comparisons from Nextflow trace
# - Memory usage plots from Nextflow trace
# - Quality metric comparisons
# - Integration metrics (iLISI, celltype_ASW)
# - Method performance evaluation
# - Trajectory analysis integration

#===============================================================================
# INTEGRATION METRICS SCORING FUNCTION
#===============================================================================

#' Score Integration Quality with Two Key Metrics
#'
#' Computes iLISI (batch mixing) and celltype ASW (biological conservation)
#' using scIntegrationMetrics package
#'
#' @param obj Seurat object with integrated data
#' @param batch_col Name of batch metadata column
#' @param label_col Name of cell type metadata column
#' @param reduction Reduction to use (default: "pca")
#' @param dims Which dimensions to use (default: 1:30)
#' @param iLISI_perplexity Perplexity for iLISI calculation (default: 20)
#'
#' @return Data frame with iLISI and celltype_ASW columns
#'
score_integration_2metrics <- function(
  obj,
  batch_col = "batch",
  label_col = "celltype",
  reduction = "pca",
  dims = 1:30,
  iLISI_perplexity = 20
) {
  # Verify reduction exists
  if (!reduction %in% Reductions(obj)) {
    stop(sprintf("Reduction '%s' not found in object.", reduction))
  }
  
  emb <- Embeddings(obj, reduction = reduction)
  if (ncol(emb) < max(dims)) {
    stop("Not enough dimensions in embedding for selected dims.")
  }
  
  # Get integration metrics
  m <- scIntegrationMetrics::getIntegrationMetrics(
    obj,
    meta.batch = batch_col,
    meta.label = label_col,
    reduction = reduction,
    dims = dims,
    iLISI_perplexity = iLISI_perplexity
  )
  
  # Extract iLISI (handles version differences)
  iLISI <- unlist(m)[c("iLISI", "norm_iLISI")]
  iLISI <- iLISI[!is.na(iLISI)][1]
  
  # Extract celltype ASW (handles version differences)
  asw <- unlist(m)[c("celltype_ASW", "ASW_celltype", "celltype_silhouette")]
  asw <- asw[!is.na(asw)][1]
  
  out <- data.frame(
    iLISI = as.numeric(iLISI),
    celltype_ASW = as.numeric(asw),
    row.names = NULL
  )
  return(out)
}

#===============================================================================
# PLOTTING FUNCTIONS
#===============================================================================

#' Create runtime by step plot from Nextflow trace data
#'
#' Visualizes total and average runtime for each pipeline step
create_runtime_by_step_plot <- function(step_summary, plots_dir, log = cat) {
  tryCatch({
    library(ggplot2)
    library(dplyr)
    
    # Sort by total runtime
    step_summary <- step_summary %>% arrange(desc(total_duration_min))
    
    # Create bar plot
    p <- ggplot(step_summary, aes(x = reorder(step, -total_duration_min), y = total_duration_min)) +
      geom_col(fill = "#3498db", alpha = 0.8) +
      geom_text(aes(label = sprintf("%.1f min", total_duration_min)), 
                vjust = -0.5, size = 4, fontface = "bold") +
      labs(
        title = "Pipeline Runtime by Step",
        x = "Pipeline Step",
        y = "Total Runtime (minutes)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank()
      )
    
    ggsave(file.path(plots_dir, "runtime_by_step.png"), 
           p, width = 12, height = 7, dpi = 150)
    
    log("Saved runtime by step plot")
    
  }, error = function(e) {
    log("Warning: Runtime plot failed:", conditionMessage(e))
  })
}

#' Create memory usage by step plot from Nextflow trace data
#'
#' Visualizes total and average memory usage for each pipeline step
create_memory_by_step_plot <- function(step_summary, plots_dir, log = cat) {
  tryCatch({
    library(ggplot2)
    library(dplyr)
    
    # Sort by total memory
    step_summary <- step_summary %>% arrange(desc(total_memory_gb))
    
    # Create bar plot
    p <- ggplot(step_summary, aes(x = reorder(step, -total_memory_gb), y = total_memory_gb)) +
      geom_col(fill = "#e74c3c", alpha = 0.8) +
      geom_text(aes(label = sprintf("%.1f GB", total_memory_gb)), 
                vjust = -0.5, size = 4, fontface = "bold") +
      labs(
        title = "Pipeline Memory Usage by Step",
        x = "Pipeline Step",
        y = "Total Memory (GB)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank()
      )
    
    ggsave(file.path(plots_dir, "memory_by_step.png"), 
           p, width = 12, height = 7, dpi = 150)
    
    log("Saved memory by step plot")
    
  }, error = function(e) {
    log("Warning: Memory plot failed:", conditionMessage(e))
  })
}

#' Create combined runtime vs memory scatter plot
#'
#' Shows relationship between runtime and memory usage per step
create_runtime_vs_memory_plot <- function(step_summary, plots_dir, log = cat) {
  tryCatch({
    library(ggplot2)
    
    p <- ggplot(step_summary, aes(x = total_duration_min, y = total_memory_gb)) +
      geom_point(aes(size = n_tasks, color = step), alpha = 0.7) +
      geom_text(aes(label = step), vjust = 1.5, size = 3) +
      scale_size_continuous(name = "Number of Tasks", range = c(3, 10)) +
      labs(
        title = "Runtime vs Memory Usage by Pipeline Step",
        x = "Total Runtime (minutes)",
        y = "Total Memory (GB)",
        color = "Step"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "right"
      )
    
    ggsave(file.path(plots_dir, "runtime_vs_memory.png"), 
           p, width = 12, height = 8, dpi = 150)
    
    log("Saved runtime vs memory plot")
    
  }, error = function(e) {
    log("Warning: Runtime vs memory plot failed:", conditionMessage(e))
  })
}

#' Create analysis metrics comparison plot
#' 
#' Visual comparison of analysis complexity by pipeline step
create_analysis_metrics_comparison <- function(metrics_df, plots_dir) {
  
  tryCatch({
    
    if (nrow(metrics_df) == 0) return(invisible(NULL))
    
    p <- ggplot(metrics_df, aes(x = step, y = complexity_score, fill = step)) +
      geom_col(alpha = 0.8) +
      labs(
        title = "Pipeline Step Analysis Metrics",
        x = "Pipeline Step",
        y = "Complexity (cells × genes / 1M)",
        subtitle = "Includes trajectory analysis integration"
      ) +
      theme_minimal(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )
    
    ggsave(file.path(plots_dir, "analysis_metrics_comparison.png"), 
           plot = p, width = 12, height = 8, dpi = 300)
    
  }, error = function(e) {
    message("Error creating metrics comparison:", e$message)
  })
}

#' Create trajectory contribution visualization
#' 
#' Shows how trajectory analysis adds to overall benchmarking
create_trajectory_contribution_plot <- function(trajectory_metrics, integration_metrics, plots_dir) {
  
  tryCatch({
    
    contrib_data <- data.frame(
      Analysis = c(
        "Integration",
        "Trajectory",
        "Total"
      ),
      Runtime_Minutes = c(
        integration_metrics$seurat_time %||% 10,
        5,  # Estimated trajectory time
        (integration_metrics$seurat_time %||% 10) + 5
      ),
      Memory_MB = c(
        integration_metrics$memory_mb %||% 500,
        200,  # Estimated trajectory memory
        (integration_metrics$memory_mb %||% 500) + 200
      ),
      Visualizations = c(
        10,
        trajectory_metrics$total_trajectory_visualizations %||% 15,
        10 + (trajectory_metrics$total_trajectory_visualizations %||% 15)
      )
    )
    
    # Create stacked bar chart
    p1 <- ggplot(contrib_data, aes(x = Analysis, y = Runtime_Minutes, fill = Analysis)) +
      geom_col(alpha = 0.7) +
      geom_text(aes(label = paste(Runtime_Minutes, "min")), vjust = -0.5, size = 5) +
      labs(
        title = "Runtime Contribution by Analysis Component",
        y = "Runtime (minutes)"
      ) +
      theme_minimal(base_size = 16) +
      theme(legend.position = "bottom")
    
    # Visualizations contribution
    p2 <- ggplot(contrib_data, aes(x = Analysis, y = Visualizations, fill = Analysis)) +
      geom_col(alpha = 0.7) +
      geom_text(aes(label = paste(Visualizations, "plots")), vjust = -0.5, size = 5) +
      labs(
        title = "Visualization Output Count",
        y = "Number of Plots/Charts"
      ) +
      theme_minimal(base_size = 16) +
      theme(legend.position = "bottom")
    
    # Combine
    combined <- patchwork::wrap_plots(p1, p2, ncol = 2)
    ggsave(file.path(plots_dir, "trajectory_contribution.png"), 
           plot = combined, width = 14, height = 6, dpi = 300)
    
  }, error = function(e) {
    message("Error creating contribution plot:", e$message)
  })
}

#===============================================================================
# INTEGRATION METRICS VISUALIZATION FUNCTIONS
#===============================================================================

#' Create Integration Metrics Comparison Plot
#'
#' Visualizes iLISI (batch mixing) and celltype_ASW (biological conservation)
#' for different integration methods
#'
#' @param metrics_df Data frame with methods and iLISI/celltype_ASW scores
#' @param plots_dir Directory to save plots
#' @param log Logging function
#'
create_integration_metrics_comparison_plot <- function(metrics_df, plots_dir, log = cat) {
  tryCatch({
    library(ggplot2)
    library(gridExtra)
    
    # Plot 1: iLISI (batch mixing)
    p1 <- ggplot(metrics_df, aes(x = reorder(Method, iLISI), y = iLISI, fill = Method)) +
      geom_col(alpha = 0.8) +
      geom_text(aes(label = sprintf("%.2f", iLISI)), vjust = -0.5, size = 5, fontface = "bold") +
      labs(
        title = "Integration Quality: Batch Mixing (iLISI)",
        x = "Integration Method",
        y = "iLISI Score",
        subtitle = "Higher = better batch removal"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12, style = "italic"),
        legend.position = "none"
      ) +
      ylim(0, max(metrics_df$iLISI) * 1.15)
    
    # Plot 2: Celltype ASW (biological conservation)
    p2 <- ggplot(metrics_df, aes(x = reorder(Method, celltype_ASW), y = celltype_ASW, fill = Method)) +
      geom_col(alpha = 0.8) +
      geom_text(aes(label = sprintf("%.2f", celltype_ASW)), vjust = -0.5, size = 5, fontface = "bold") +
      labs(
        title = "Integration Quality: Biological Conservation (Celltype ASW)",
        x = "Integration Method",
        y = "Celltype ASW",
        subtitle = "Higher = better preservation of cell types"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12, style = "italic"),
        legend.position = "none"
      ) +
      ylim(0, max(metrics_df$celltype_ASW) * 1.15)
    
    # Combined plot
    p <- gridExtra::grid.arrange(p1, p2, ncol = 2)
    
    ggsave(file.path(plots_dir, "integration_metrics_comparison.png"),
           plot = p, width = 14, height = 6, dpi = 150)
    
    log("Saved integration metrics comparison plot")
    
  }, error = function(e) {
    log("Warning: Integration metrics plot failed:", conditionMessage(e))
  })
}

#' Create Integration Metrics Trade-off Scatter Plot
#'
#' Shows the trade-off between batch removal (iLISI) and biological conservation (ASW)
#'
create_integration_tradeoff_plot <- function(metrics_df, plots_dir, log = cat) {
  tryCatch({
    library(ggplot2)
    
    p <- ggplot(metrics_df, aes(x = iLISI, y = celltype_ASW, color = Method, size = 5)) +
      geom_point(alpha = 0.8) +
      geom_text(aes(label = Method), vjust = -1.5, size = 4, fontface = "bold") +
      labs(
        title = "Integration Trade-off: Batch Removal vs Biological Conservation",
        x = "iLISI (Batch Mixing) →",
        y = "Celltype ASW (Conservation) →",
        subtitle = "Ideal: High iLISI + High Celltype ASW (top-right corner)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, style = "italic", hjust = 0.5),
        legend.position = "none",
        panel.grid = element_line(color = "gray90")
      ) +
      xlim(min(metrics_df$iLISI) * 0.9, max(metrics_df$iLISI) * 1.1) +
      ylim(min(metrics_df$celltype_ASW) * 0.9, max(metrics_df$celltype_ASW) * 1.1)
    
    ggsave(file.path(plots_dir, "integration_tradeoff_plot.png"),
           plot = p, width = 10, height = 8, dpi = 150)
    
    log("Saved integration trade-off plot")
    
  }, error = function(e) {
    log("Warning: Integration trade-off plot failed:", conditionMessage(e))
  })
}

#' Create Integration Metrics Summary Table Visualization
#'
#' Creates a formatted table showing key integration metrics
#'
create_integration_metrics_table <- function(metrics_df, plots_dir, log = cat) {
  tryCatch({
    library(gridExtra)
    library(grid)
    
    # Prepare table data with formatting
    table_data <- metrics_df %>%
      mutate(
        Method = as.character(Method),
        iLISI = sprintf("%.3f", iLISI),
        celltype_ASW = sprintf("%.3f", celltype_ASW)
      ) %>%
      select(Method, iLISI, celltype_ASW)
    
    # Create table grob
    table_grob <- gridExtra::tableGrob(
      table_data,
      rows = NULL,
      theme = gridExtra::ttheme_default(
        base_size = 12,
        colhead = list(
          fg_params = list(fontface = "bold", cex = 1.2),
          bg_params = list(fill = "#3498db", alpha = 0.8)
        ),
        row = list(
          fg_params = list(fontface = "plain"),
          bg_params = list(fill = c("#ecf0f1", "white"))
        )
      )
    )
    
    # Add title
    title <- grid::textGrob(
      "Integration Method Quality Metrics",
      gp = grid::gpar(fontsize = 16, fontface = "bold"),
      hjust = 0.5
    )
    
    # Combine title and table
    p <- gridExtra::arrangeGrob(
      title,
      table_grob,
      heights = c(0.1, 0.9)
    )
    
    ggsave(file.path(plots_dir, "integration_metrics_table.png"),
           plot = p, width = 8, height = 4, dpi = 150)
    
    log("Saved integration metrics table")
    
  }, error = function(e) {
    log("Warning: Integration metrics table failed:", conditionMessage(e))
  })
}#===============================================================================
# Memory Analysis Plotting Functions
# Visualizations for memory profiling and efficiency metrics
#===============================================================================

#' Create memory progression plot
#' @param memory_data Memory tracking data frame
#' @param plots_dir Output directory
create_memory_progression_plot <- function(memory_data, plots_dir) {
  
  if (is.null(memory_data) || nrow(memory_data) == 0) return(invisible(NULL))
  
  tryCatch({
    # Parse timestamps
    memory_data$time_num <- as.numeric(as.POSIXct(memory_data$timestamp))
    memory_data$time_num <- memory_data$time_num - min(memory_data$time_num)
    
    p <- ggplot(memory_data, aes(x = time_num, y = system_memory_mb)) +
      geom_area(alpha = 0.3, fill = "steelblue") +
      geom_line(size = 1.2, color = "darkblue") +
      geom_point(size = 2, color = "darkblue") +
      geom_text(aes(label = checkpoint), vjust = -0.7, size = 3) +
      labs(
        title = "Memory Usage Progression",
        x = "Time (seconds)",
        y = "System Memory Used (MB)",
        subtitle = "Shows memory consumption across analysis steps"
      ) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(plots_dir, "memory_progression.png"), 
           plot = p, width = 12, height = 7, dpi = 300)
    
  }, error = function(e) {
    message("Error creating memory progression plot:", e$message)
  })
}

#' Create memory efficiency comparison
#' @param efficiency_data Efficiency metrics data frame
#' @param plots_dir Output directory
create_memory_efficiency_plot <- function(efficiency_data, plots_dir) {
  
  if (is.null(efficiency_data) || nrow(efficiency_data) == 0) return(invisible(NULL))
  
  tryCatch({
    p <- ggplot(efficiency_data, aes(x = reorder(step, efficiency_score), 
                                      y = efficiency_score, 
                                      fill = efficiency_score)) +
      geom_col(alpha = 0.8) +
      scale_fill_gradient(low = "red", high = "green") +
      geom_text(aes(label = sprintf("%.4f", efficiency_score)), 
                vjust = -0.5, size = 4, fontface = "bold") +
      labs(
        title = "Memory Efficiency by Pipeline Step",
        x = "Pipeline Step",
        y = "Efficiency Score (higher = better)",
        subtitle = "Score = (Cells × Genes) / Memory Used"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    ggsave(file.path(plots_dir, "memory_efficiency.png"), 
           plot = p, width = 12, height = 7, dpi = 300)
    
  }, error = function(e) {
    message("Error creating memory efficiency plot:", e$message)
  })
}

#' Create memory per cell/gene analysis
#' @param efficiency_data Efficiency metrics data frame
#' @param plots_dir Output directory
create_memory_per_unit_plot <- function(efficiency_data, plots_dir) {
  
  if (is.null(efficiency_data) || nrow(efficiency_data) == 0) return(invisible(NULL))
  
  tryCatch({
    # Create two subplots
    p1 <- ggplot(efficiency_data, aes(x = reorder(step, memory_per_cell_kb), 
                                       y = memory_per_cell_kb,
                                       fill = memory_per_cell_kb)) +
      geom_col(alpha = 0.8) +
      scale_fill_gradient(low = "green", high = "red") +
      geom_text(aes(label = sprintf("%.4f", memory_per_cell_kb)), 
                vjust = -0.5, size = 3) +
      labs(
        title = "Memory per Cell",
        x = "",
        y = "KB per Cell"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    p2 <- ggplot(efficiency_data, aes(x = reorder(step, memory_per_gene_kb), 
                                       y = memory_per_gene_kb,
                                       fill = memory_per_gene_kb)) +
      geom_col(alpha = 0.8) +
      scale_fill_gradient(low = "green", high = "red") +
      geom_text(aes(label = sprintf("%.4f", memory_per_gene_kb)), 
                vjust = -0.5, size = 3) +
      labs(
        title = "Memory per Gene",
        x = "",
        y = "KB per Gene"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    combined <- patchwork::wrap_plots(p1, p2, ncol = 2)
    
    ggsave(file.path(plots_dir, "memory_per_unit.png"), 
           plot = combined, width = 14, height = 7, dpi = 300)
    
  }, error = function(e) {
    message("Error creating memory per unit plot:", e$message)
  })
}

#' Create memory load analysis
#' @param memory_logs List of memory tracking logs
#' @param plots_dir Output directory
create_memory_load_analysis_plot <- function(memory_logs, plots_dir) {
  
  if (is.null(memory_logs) || length(memory_logs) == 0) return(invisible(NULL))
  
  tryCatch({
    load_data <- data.frame(
      file = sapply(memory_logs, function(x) x$file),
      load_time = sapply(memory_logs, function(x) x$load_time_sec),
      load_speed = sapply(memory_logs, function(x) x$load_speed_mb_per_sec),
      file_size = sapply(memory_logs, function(x) x$file_size_mb),
      stringsAsFactors = FALSE
    )
    
    p <- ggplot(load_data, aes(x = file_size, y = load_time, size = load_speed)) +
      geom_point(alpha = 0.6, color = "steelblue") +
      geom_text(aes(label = file), size = 3, vjust = -0.5) +
      scale_size(range = c(2, 8)) +
      labs(
        title = "Memory Loading Performance",
        x = "File Size (MB)",
        y = "Load Time (seconds)",
        size = "Load Speed\n(MB/sec)",
        subtitle = "Larger points = faster loading"
      ) +
      theme_minimal(base_size = 14)
    
    ggsave(file.path(plots_dir, "memory_load_analysis.png"), 
           plot = p, width = 12, height = 8, dpi = 300)
    
  }, error = function(e) {
    message("Error creating load analysis plot:", e$message)
  })
}

#' Create cumulative memory plot
#' @param memory_logs List of memory tracking logs
#' @param plots_dir Output directory
create_cumulative_memory_plot <- function(memory_logs, plots_dir) {
  
  if (is.null(memory_logs) || length(memory_logs) == 0) return(invisible(NULL))
  
  tryCatch({
    memory_data <- data.frame(
      file = sapply(memory_logs, function(x) x$file),
      memory = sapply(memory_logs, function(x) x$memory_increase_mb),
      stringsAsFactors = FALSE
    )
    
    memory_data <- memory_data[order(memory_data$memory, decreasing = TRUE), ]
    memory_data$cumsum <- cumsum(memory_data$memory)
    memory_data$file <- factor(memory_data$file, levels = memory_data$file)
    
    p <- ggplot(memory_data, aes(x = file, y = cumsum)) +
      geom_col(aes(fill = memory), alpha = 0.8) +
      geom_line(aes(group = 1), size = 1.2, color = "darkred") +
      geom_point(size = 3, color = "darkred") +
      geom_text(aes(label = sprintf("%.0f MB", cumsum)), 
                vjust = -0.5, size = 3.5) +
      scale_fill_gradient(low = "lightblue", high = "darkblue") +
      labs(
        title = "Cumulative Memory Usage",
        x = "Loading Order (by memory size)",
        y = "Cumulative Memory (MB)",
        fill = "Step Memory\n(MB)"
      ) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(plots_dir, "cumulative_memory.png"), 
           plot = p, width = 12, height = 7, dpi = 300)
    
  }, error = function(e) {
    message("Error creating cumulative memory plot:", e$message)
  })
}

#' Create memory utilization heatmap
#' @param efficiency_data Efficiency metrics data frame
#' @param plots_dir Output directory
create_memory_utilization_heatmap <- function(efficiency_data, plots_dir) {
  
  if (is.null(efficiency_data) || nrow(efficiency_data) == 0) return(invisible(NULL))
  
  tryCatch({
    # Prepare data for heatmap
    heatmap_data <- efficiency_data[, c("step", "memory_used_mb", 
                                         "memory_per_cell_kb", "memory_per_gene_kb",
                                         "efficiency_score")]
    rownames(heatmap_data) <- heatmap_data$step
    heatmap_data$step <- NULL
    
    # Normalize for visualization
    heatmap_norm <- sweep(heatmap_data, 2, colMeans(heatmap_data, na.rm = TRUE), "-")
    heatmap_norm <- sweep(heatmap_norm, 2, apply(heatmap_data, 2, sd, na.rm = TRUE), "/")
    
    png(file.path(plots_dir, "memory_utilization_heatmap.png"), 
        width = 10, height = 8, units = "in", res = 300)
    
    pheatmap::pheatmap(
      as.matrix(heatmap_norm),
      main = "Memory Utilization Heatmap\n(Standardized Values)",
      color = colorRampPalette(c("blue", "white", "red"))(50),
      breaks = seq(-3, 3, by = 0.12),
      display_numbers = round(as.matrix(heatmap_data), 2),
      fontsize = 10,
      cluster_rows = TRUE,
      cluster_cols = FALSE
    )
    
    dev.off()
    
  }, error = function(e) {
    message("Error creating heatmap:", e$message)
  })
}

#' Create memory snapshot comparison plot
#' @param comparisons List of memory snapshot comparisons
#' @param plots_dir Output directory
create_memory_snapshot_comparison <- function(comparisons, plots_dir) {
  
  if (is.null(comparisons) || length(comparisons) == 0) return(invisible(NULL))
  
  tryCatch({
    # Prepare comparison data
    comp_data <- data.frame(
      step = names(comparisons),
      memory_increase = sapply(comparisons, function(x) x$memory_delta_mb),
      stringsAsFactors = FALSE
    )
    
    comp_data <- comp_data[order(comp_data$memory_increase, decreasing = TRUE), ]
    comp_data$step <- factor(comp_data$step, levels = comp_data$step)
    
    p <- ggplot(comp_data, aes(x = step, y = memory_increase, 
                                fill = memory_increase > 0)) +
      geom_col(alpha = 0.8) +
      geom_text(aes(label = sprintf("%+.0f MB", memory_increase)), 
                vjust = ifelse(comp_data$memory_increase > 0, -0.5, 1.5), 
                size = 4, fontface = "bold") +
      scale_fill_manual(values = c("red", "green")) +
      labs(
        title = "Memory Change by Step",
        x = "Pipeline Step",
        y = "Memory Change (MB)",
        subtitle = "Green = increase | Red = decrease"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    ggsave(file.path(plots_dir, "memory_snapshot_comparison.png"), 
           plot = p, width = 12, height = 7, dpi = 300)
    
  }, error = function(e) {
    message("Error creating snapshot comparison plot:", e$message)
  })
}

# Plot top 10 longest running Nextflow processes
create_top_processes_plot <- function(process_summary, plots_dir, log) {
  if(nrow(process_summary) == 0) {
    log("No process data available for top processes plot")
    return()
  }
  
  top_processes <- process_summary %>% 
    head(10) %>%
    ggplot(aes(x = reorder(process_id, -duration_min), y = duration_min)) +
    geom_col(fill = "#2c3e50", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f m", duration_min)), 
              vjust = -0.3, size = 3.5, fontface = "bold") +
    labs(
      title = "Top 10 Longest Running Nextflow Processes",
      x = "Process", 
      y = "Runtime (minutes)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  
  suppressWarnings(ggsave(file.path(plots_dir, "top_10_processes_runtime.png"), 
                         top_processes, width = 12, height = 6, dpi = 150))
  log("Saved top 10 longest processes plot")
}

# Plot integration metrics comparison (bar plot with color coding)
create_integration_metrics_comparison <- function(metrics_comparison, plots_dir, log) {
  metrics_long <- metrics_comparison %>%
    dplyr::select(Method, norm_iLISI, celltype_ASW) %>%
    tidyr::pivot_longer(cols = -Method, names_to = "Metric", values_to = "Score") %>%
    dplyr::mutate(
      Metric = dplyr::case_when(
        Metric == "norm_iLISI" ~ "Batch Mixing",
        Metric == "celltype_ASW" ~ "Cell Type Sep.",
        TRUE ~ Metric
      )
    ) %>%
    dplyr::group_by(Metric) %>%
    dplyr::mutate(
      is_best = !is.na(Score) & Score == max(Score, na.rm = TRUE),
      color = dplyr::case_when(
        is_best & !is.na(Score) ~ "#27ae60",  # Green for best
        is.na(Score) ~ "#95a5a6",              # Gray for NA
        TRUE ~ "#e74c3c"                       # Red for others
      )
    ) %>%
    dplyr::ungroup()
  
  p_metrics <- ggplot(metrics_long, aes(x = Method, y = Score, fill = color)) +
    geom_col(alpha = 0.85, width = 0.6) +
    geom_text(aes(label = sprintf("%.4f", Score)), 
              vjust = -0.3, fontface = "bold", size = 4) +
    facet_wrap(~Metric, scales = "free_y", ncol = 2) +
    scale_fill_identity() +
    labs(
      title = "Integration Quality Metrics Comparison",
      subtitle = "Green = Best method | Red = Worse method",
      x = "Integration Method",
      y = "Score"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray50"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    )
  
  suppressWarnings(ggsave(file.path(plots_dir, "integration_metrics_comparison.png"), 
                         p_metrics, width = 12, height = 7, dpi = 150))
  log("Saved integration metrics comparison plot")
}

# Plot integration metrics heatmap
create_integration_metrics_heatmap <- function(metrics_comparison, plots_dir, log) {
  metrics_heatmap_df <- metrics_comparison %>%
    dplyr::select(Method, norm_iLISI, celltype_ASW) %>%
    tidyr::pivot_longer(cols = -Method, names_to = "Metric", values_to = "Score") %>%
    dplyr::mutate(
      Metric = dplyr::case_when(
        Metric == "norm_iLISI" ~ "Batch Mixing",
        Metric == "celltype_ASW" ~ "Cell Type Sep.",
        TRUE ~ Metric
      ),
      Method = factor(Method, levels = unique(Method))
    ) %>%
    dplyr::group_by(Metric) %>%
    dplyr::mutate(
      Best = Score == max(Score, na.rm = TRUE),
      Max_score = max(Score, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()
  
  p_heatmap <- ggplot(metrics_heatmap_df, aes(x = Metric, y = Method, fill = Score)) +
    geom_tile(color = "white", size = 1.5) +
    geom_text(aes(label = sprintf("%.4f", Score), 
                  color = ifelse(Best, "white", "black")),
              fontface = "bold", size = 5) +
    scale_fill_gradient(low = "#e74c3c", high = "#27ae60", 
                       limits = c(0, max(metrics_heatmap_df$Max_score, na.rm = TRUE) * 1.1),
                       na.value = "#95a5a6") +
    scale_color_identity() +
    labs(
      title = "Integration Metrics Heatmap",
      subtitle = "Darker Green = Better | Darker Red = Worse",
      fill = "Score"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray50"),
      axis.text = element_text(face = "bold"),
      axis.title = element_blank(),
      legend.position = "right"
    )
  
  suppressWarnings(ggsave(file.path(plots_dir, "integration_metrics_heatmap.png"), 
                         p_heatmap, width = 8, height = 6, dpi = 150))
  log("Saved integration metrics heatmap")
}

# Plot integration radar chart comparison
create_integration_radar_chart <- function(comp_df, plots_dir, log, CFG_PLOT_DPI) {
  library(fmsb)
  
  radar_data <- comp_df %>%
    select(Method, Overall_Score) %>%
    column_to_rownames("Method") %>%
    t()
  
  radar_data <- rbind(
    rep(1, ncol(radar_data)),  # Max
    rep(0, ncol(radar_data)),  # Min
    radar_data
  )
  
  png(file.path(plots_dir, "integration_radar_comparison.png"), 
      width=10, height=10, units="in", res=CFG_PLOT_DPI)
  par(mar=c(1,1,3,1))
  radarchart(radar_data, 
             axistype=1,
             pcol=c("#E41A1C", "#377EB8"),
             pfcol=c(rgb(0.9,0.1,0.1,0.3), rgb(0.2,0.5,0.8,0.3)),
             plwd=3,
             plty=1,
             cglcol="grey", cglty=1, axislabcol="grey20",
             caxislabels=seq(0,1,0.25), cglwd=0.8,
             vlcex=1.2,
             title="Integration Methods Performance Comparison")
  legend("topright", legend=comp_df$Method, 
         col=c("#E41A1C", "#377EB8"), lty=1, lwd=3, cex=1.2)
  dev.off()
  log("  Saved radar chart comparison")
}

# Plot overall integration score comparison
create_integration_overall_score <- function(comp_df, CFG_BASE_SIZE, CFG_PLOT_WIDTH, CFG_PLOT_HEIGHT, CFG_PLOT_DPI, save_dual_plot) {
  p_overall <- ggplot(comp_df, aes(x=Method, y=Overall_Score, fill=Method)) +
    geom_col(alpha=0.8, show.legend=FALSE) +
    geom_text(aes(label=sprintf("%.3f", Overall_Score)), 
              vjust=-0.5, size=5, fontface="bold") +
    labs(title="Integration Method Overall Performance Score", 
         x="Integration Method", 
         y="Overall Score (0-1, higher is better)",
         subtitle=sprintf("Best Method: %s | Score: %.3f",
                         comp_df$Method[which.max(comp_df$Overall_Score)],
                         max(comp_df$Overall_Score)),
         caption="Score = 40% Time + 30% Memory + 30% Quality") +
    scale_fill_brewer(palette="Set1") +
    ylim(0, 1) +
    theme_minimal(base_size=CFG_BASE_SIZE) +
    theme(plot.caption=element_text(hjust=0, size=12, color="gray40"))
  save_dual_plot(p_overall, "integration_overall_score", 
                width=CFG_PLOT_WIDTH, height=CFG_PLOT_HEIGHT, dpi=CFG_PLOT_DPI)
}

# Plot integration metrics table comparison
create_integration_metrics_table <- function(comp_df, CFG_BASE_SIZE, CFG_PLOT_DPI, save_dual_plot) {
  comp_long <- comp_df %>%
    select(Method, Integration_Time_Min, Memory_Used_MB, N_Clusters, Silhouette_Score) %>%
    tidyr::pivot_longer(cols=-Method, names_to="Metric", values_to="Value")
  
  p_metrics_table <- ggplot(comp_long, aes(x=Metric, y=Value, fill=Method)) +
    geom_col(position="dodge", alpha=0.8) +
    facet_wrap(~Metric, scales="free_y", ncol=2) +
    labs(title="Integration Methods: Key Metrics Comparison") +
    theme_minimal(base_size=CFG_BASE_SIZE) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.position="bottom")
  save_dual_plot(p_metrics_table, "integration_metrics_table", 
                width=14, height=10, dpi=CFG_PLOT_DPI)
}

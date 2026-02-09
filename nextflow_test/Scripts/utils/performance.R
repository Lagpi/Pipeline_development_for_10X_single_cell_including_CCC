#===============================================================================
# Performance Tracking Utilities
# Track runtime and memory usage during pipeline execution
#===============================================================================

#' Initialize performance tracking for a module
#' @param module_name Name of the pipeline module (e.g., "01_qc_integration")
#' @param output_dir Base output directory
#' @return List with start_time, start_memory, module_name, perf_file
start_performance_tracking <- function(module_name, output_dir) {
  perf_dir <- file.path(output_dir, "performance_logs")
  dir.create(perf_dir, showWarnings = FALSE, recursive = TRUE)
  
  perf_file <- file.path(perf_dir, paste0(module_name, "_performance.rds"))
  
  perf_data <- list(
    module_name = module_name,
    start_time = Sys.time(),
    start_memory_mb = as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo 2>/dev/null || echo 0", intern = TRUE)) / 1024,
    pid = Sys.getpid(),
    hostname = Sys.info()["nodename"],
    r_version = paste(R.version$major, R.version$minor, sep = ".")
  )
  
  # Save initial state
  saveRDS(perf_data, perf_file)
  
  return(list(
    start_time = perf_data$start_time,
    start_memory = perf_data$start_memory_mb,
    module_name = module_name,
    perf_file = perf_file,
    peak_memory = perf_data$start_memory_mb
  ))
}

#' Update peak memory during execution
#' @param tracker Performance tracker object from start_performance_tracking
update_peak_memory <- function(tracker) {
  current_mem <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo 2>/dev/null || echo 0", intern = TRUE)) / 1024
  if (current_mem < tracker$peak_memory) {
    tracker$peak_memory <- tracker$start_memory - current_mem
  }
  return(tracker)
}

#' Finalize performance tracking for a module
#' @param tracker Performance tracker object from start_performance_tracking
#' @param success Whether the module completed successfully
#' @param additional_metrics Optional list of additional metrics to save
end_performance_tracking <- function(tracker, success = TRUE, additional_metrics = NULL) {
  end_time <- Sys.time()
  end_memory_mb <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo 2>/dev/null || echo 0", intern = TRUE)) / 1024
  
  # Calculate metrics
  elapsed_time_sec <- as.numeric(difftime(end_time, tracker$start_time, units = "secs"))
  elapsed_time_min <- elapsed_time_sec / 60
  memory_used_mb <- tracker$start_memory - end_memory_mb
  peak_memory_mb <- tracker$peak_memory
  
  # Load existing data and update
  if (file.exists(tracker$perf_file)) {
    perf_data <- readRDS(tracker$perf_file)
  } else {
    perf_data <- list(module_name = tracker$module_name)
  }
  
  perf_data$end_time <- end_time
  perf_data$elapsed_time_sec <- elapsed_time_sec
  perf_data$elapsed_time_min <- elapsed_time_min
  perf_data$elapsed_time_formatted <- sprintf("%.2f min (%.1f sec)", elapsed_time_min, elapsed_time_sec)
  perf_data$memory_used_mb <- memory_used_mb
  perf_data$peak_memory_mb <- peak_memory_mb
  perf_data$end_memory_mb <- end_memory_mb
  perf_data$success <- success
  perf_data$completed_at <- as.character(end_time)
  
  # Add any additional metrics
  if (!is.null(additional_metrics)) {
    perf_data <- c(perf_data, additional_metrics)
  }
  
  # Save final state
  saveRDS(perf_data, tracker$perf_file)
  
  # Also write human-readable version
  csv_file <- gsub("\\.rds$", ".csv", tracker$perf_file)
  write.csv(as.data.frame(perf_data[sapply(perf_data, function(x) length(x) == 1)]), 
            csv_file, row.names = FALSE)
  
  # Print summary
  cat("\n=== Performance Summary ===\n")
  cat("Module:", tracker$module_name, "\n")
  cat("Runtime:", perf_data$elapsed_time_formatted, "\n")
  cat("Memory used:", round(memory_used_mb, 2), "MB\n")
  cat("Peak memory:", round(peak_memory_mb, 2), "MB\n")
  cat("Status:", ifelse(success, "SUCCESS", "FAILED"), "\n")
  cat("===========================\n\n")
  
  return(perf_data)
}

#' Get current memory usage in MB
get_current_memory_mb <- function() {
  mem_free <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo 2>/dev/null || echo 0", intern = TRUE)) / 1024
  mem_total <- as.numeric(system("awk '/MemTotal/ {print $2}' /proc/meminfo 2>/dev/null || echo 1", intern = TRUE)) / 1024
  return(mem_total - mem_free)
}

#' Get current memory usage in bytes
#'
#' Wrapper used by benchmarking scripts expecting bytes
get_memory_usage <- function() {
  get_current_memory_mb() * 1024^2
}

#' Detailed memory profiling during execution
#' Tracks memory usage at specific checkpoints
start_detailed_memory_tracking <- function(module_name, output_dir) {
  perf_dir <- file.path(output_dir, "performance_logs")
  dir.create(perf_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get initial system state
  initial_memory <- list(
    timestamp = Sys.time(),
    memory_used_mb = get_current_memory_mb(),
    gc_stats = gc()
  )
  
  # Also track R heap memory
  if (requireNamespace("pryr", quietly = TRUE)) {
    initial_memory$r_memory_mb <- as.numeric(pryr::mem_used()) / 1024^2
  }
  
  tracker <- list(
    module_name = module_name,
    memory_file = file.path(perf_dir, paste0(module_name, "_memory_detailed.csv")),
    checkpoints = list(),
    initial_memory = initial_memory,
    checkpoint_count = 0
  )
  
  # Initialize checkpoint file
  checkpoint_df <- data.frame(
    timestamp = format(initial_memory$timestamp, "%Y-%m-%d %H:%M:%S"),
    checkpoint = "START",
    system_memory_mb = initial_memory$memory_used_mb,
    r_memory_mb = if(!is.null(initial_memory$r_memory_mb)) initial_memory$r_memory_mb else NA,
    gc_collections = NA,
    description = "Initial state",
    stringsAsFactors = FALSE
  )
  
  write.csv(checkpoint_df, tracker$memory_file, row.names = FALSE)
  
  return(tracker)
}

#' Record memory checkpoint during execution
#' @param tracker Tracker from start_detailed_memory_tracking
#' @param checkpoint_name Name of the checkpoint
#' @param description Optional description
record_memory_checkpoint <- function(tracker, checkpoint_name, description = "") {
  
  current_time <- Sys.time()
  system_memory <- get_current_memory_mb()
  
  # Track R memory if available
  r_memory <- NA
  if (requireNamespace("pryr", quietly = TRUE)) {
    r_memory <- as.numeric(pryr::mem_used()) / 1024^2
  }
  
  # Force garbage collection and get stats
  gc_stats <- gc()
  
  # Create checkpoint record
  checkpoint_data <- data.frame(
    timestamp = format(current_time, "%Y-%m-%d %H:%M:%S"),
    checkpoint = checkpoint_name,
    system_memory_mb = round(system_memory, 2),
    r_memory_mb = if(!is.na(r_memory)) round(r_memory, 2) else NA,
    gc_collections = sum(gc_stats[, "used.cells"]),
    description = description,
    stringsAsFactors = FALSE
  )
  
  # Append to tracking file
  tryCatch({
    write.csv(checkpoint_data, tracker$memory_file, 
              row.names = FALSE, append = TRUE, col.names = FALSE)
  }, error = function(e) {
    warning("Failed to write memory checkpoint: ", e$message)
  })
  
  # Update tracker
  tracker$checkpoint_count <- tracker$checkpoint_count + 1
  tracker$checkpoints[[checkpoint_name]] <- checkpoint_data
  
  return(tracker)
}

#' Analyze memory usage across checkpoints
#' @param tracker Tracker from start_detailed_memory_tracking
#' @return Data frame with memory analysis
analyze_memory_profile <- function(tracker) {
  
  if (!file.exists(tracker$memory_file)) {
    warning("Memory tracking file not found: ", tracker$memory_file)
    return(NULL)
  }
  
  # Read memory checkpoints
  memory_data <- read.csv(tracker$memory_file, stringsAsFactors = FALSE)
  
  if (nrow(memory_data) == 0) {
    warning("No memory checkpoints recorded")
    return(NULL)
  }
  
  # Calculate memory deltas
  memory_data$memory_delta_mb <- c(0, diff(memory_data$system_memory_mb))
  memory_data$r_memory_delta_mb <- c(NA, diff(memory_data$r_memory_mb))
  
  # Calculate peak memory
  peak_memory <- max(memory_data$system_memory_mb, na.rm = TRUE)
  
  # Identify memory spikes (>10% increase from previous)
  memory_data$is_spike <- FALSE
  for (i in 2:nrow(memory_data)) {
    prev_mem <- memory_data$system_memory_mb[i-1]
    curr_mem <- memory_data$system_memory_mb[i]
    if (!is.na(prev_mem) && !is.na(curr_mem)) {
      pct_increase <- (curr_mem - prev_mem) / prev_mem * 100
      if (pct_increase > 10) {
        memory_data$is_spike[i] <- TRUE
      }
    }
  }
  
  return(memory_data)
}

#' Generate memory report
#' @param tracker Tracker from start_detailed_memory_tracking
#' @param output_file Path to save report
generate_memory_report <- function(tracker, output_file = NULL) {
  
  if (is.null(output_file)) {
    output_file <- gsub("_detailed.csv$", "_report.txt", tracker$memory_file)
  }
  
  analysis <- analyze_memory_profile(tracker)
  
  if (is.null(analysis)) {
    warning("Cannot generate memory report: analysis failed")
    return(NULL)
  }
  
  # Generate report
  report <- sprintf(
    "================================================================================\n
DETAILED MEMORY PROFILING REPORT: %s
================================================================================\n\n
SUMMARY STATISTICS:\n
  Peak Memory Usage: %.2f MB\n
  Final Memory Usage: %.2f MB\n
  Total Checkpoints: %d\n
  Memory Spikes Detected: %d\n\n
MEMORY BY CHECKPOINT:\n%s\n\n
MEMORY EFFICIENCY:\n
  Average Memory per Checkpoint: %.2f MB\n
  Std Dev Memory: %.2f MB\n
  Coefficient of Variation: %.2f%%\n\n
RECOMMENDATIONS:\n%s\n
================================================================================\n",
    tracker$module_name,
    max(analysis$system_memory_mb, na.rm = TRUE),
    tail(analysis$system_memory_mb, 1),
    nrow(analysis),
    sum(analysis$is_spike, na.rm = TRUE),
    
    # Memory by checkpoint
    paste(
      sprintf("  %s: %.2f MB (%+.2f MB from previous)",
              analysis$checkpoint,
              analysis$system_memory_mb,
              analysis$memory_delta_mb),
      collapse = "\n"
    ),
    
    # Efficiency metrics
    mean(analysis$system_memory_mb, na.rm = TRUE),
    sd(analysis$system_memory_mb, na.rm = TRUE),
    (sd(analysis$system_memory_mb, na.rm = TRUE) / mean(analysis$system_memory_mb, na.rm = TRUE)) * 100,
    
    # Recommendations
    generate_memory_recommendations(analysis)
  )
  
  # Write report
  writeLines(report, output_file)
  cat(report)
  
  return(report)
}

#' Generate memory optimization recommendations
#' @param analysis Memory analysis data frame
#' @return Character string with recommendations
generate_memory_recommendations <- function(analysis) {
  
  recommendations <- character()
  
  # Check for large spikes
  spikes <- which(analysis$is_spike)
  if (length(spikes) > 0) {
    spike_checkpoints <- analysis$checkpoint[spikes]
    recommendations <- c(recommendations,
      sprintf("  ⚠ Large memory spike(s) at: %s\n    Consider memory optimization for these steps",
              paste(spike_checkpoints, collapse = ", ")))
  }
  
  # Check for memory leak patterns
  memory_trend <- diff(analysis$system_memory_mb)
  if (sum(memory_trend > 0, na.rm = TRUE) > nrow(analysis) * 0.5) {
    recommendations <- c(recommendations,
      "  ⚠ Increasing memory trend detected - possible memory leak\n    Check for proper object cleanup and garbage collection")
  }
  
  # Check peak memory
  peak_mem <- max(analysis$system_memory_mb, na.rm = TRUE)
  if (peak_mem > 16000) {
    recommendations <- c(recommendations,
      sprintf("  ⚠ High peak memory usage: %.2f MB\n    Consider processing in chunks or using sparse matrices", peak_mem))
  }
  
  if (length(recommendations) == 0) {
    recommendations <- "  ✓ Memory usage appears normal - no major issues detected"
  }
  
  return(paste(recommendations, collapse = "\n"))
}

#' Collect all performance logs from completed modules
#' @param output_dir Base output directory containing performance_logs
#' @return Data frame with performance metrics for all modules
collect_performance_metrics <- function(output_dir) {
  perf_dir <- file.path(output_dir, "performance_logs")
  
  if (!dir.exists(perf_dir)) {
    warning("Performance logs directory not found: ", perf_dir)
    return(NULL)
  }
  
  perf_files <- list.files(perf_dir, pattern = "_performance\\.rds$", full.names = TRUE)
  
  if (length(perf_files) == 0) {
    warning("No performance log files found in: ", perf_dir)
    return(NULL)
  }
  
  perf_list <- lapply(perf_files, function(f) {
    tryCatch({
      data <- readRDS(f)
      # Convert to single-row data frame
      df_data <- data[sapply(data, function(x) length(x) == 1)]
      as.data.frame(df_data, stringsAsFactors = FALSE)
    }, error = function(e) {
      warning("Failed to read ", basename(f), ": ", e$message)
      NULL
    })
  })
  
  # Remove NULL entries
  perf_list <- perf_list[!sapply(perf_list, is.null)]
  
  if (length(perf_list) == 0) {
    return(NULL)
  }
  
  # Combine into single data frame
  perf_df <- do.call(rbind, perf_list)
  
  # Sort by module name
  if ("module_name" %in% colnames(perf_df)) {
    perf_df <- perf_df[order(perf_df$module_name), ]
  }
  
  rownames(perf_df) <- NULL
  
  return(perf_df)
}

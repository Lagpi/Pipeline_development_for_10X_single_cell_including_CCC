# ==============================================================================
# Common Utility Functions
# ==============================================================================
# Shared helper functions used across all analysis scripts
# 
# Functions:
# - setup_analysis_directories(): Create standard output directory structure
# - safe_lib(): Safely load optional packages
# - has_pkg(): Check if package is available
# - %||%: Null-coalescing operator
# ==============================================================================

#' Setup analysis output directories
#' 
#' Creates standardized directory structure for analysis outputs:
#' - data/: CSV files and tables
#' - plots/: Static visualizations
#' - interactive/: Interactive HTML reports
#' - rds_objects/: Serialized R objects
#' 
#' @param base_dir Base output directory path
#' @param step_name Analysis step name (e.g., "01_qc_integration")
#' @return Named list with paths: data_dir, plots_dir, dash_dir, rds_dir
#' @examples
#' dirs <- setup_analysis_directories("/path/to/outputs", "01_qc_integration")
#' write.csv(data, file.path(dirs$data_dir, "results.csv"))
setup_analysis_directories <- function(base_dir, step_name) {
  data_dir <- file.path(base_dir, "data")
  plots_dir <- file.path(base_dir, "plots")
  dash_dir <- file.path(base_dir, "interactive")
  rds_dir <- file.path(base_dir, "rds_objects")
  
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dash_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
  
  list(
    data_dir = data_dir,
    plots_dir = plots_dir,
    dash_dir = dash_dir,
    rds_dir = rds_dir
  )
}

#' Setup analysis environment (directories + logging)
#' 
#' Combines directory creation and logging setup in one call.
#' Automatically assigns directory paths to parent environment.
#' 
#' @param base_dir Base output directory path
#' @param step_name Analysis step name (e.g., "01_qc_integration")
#' @param assign_to_parent If TRUE, assigns dir variables to parent environment (default: TRUE)
#' @return Named list with paths: base_dir, data_dir, plots_dir, dash_dir, rds_dir, log_file
#' @examples
#' setup_analysis_env(opt$output_dir, "01_qc_integration")
#' # Now data_dir, plots_dir, dash_dir are available in current environment
setup_analysis_env <- function(base_dir, step_name, assign_to_parent = TRUE) {
  # Create directories
  dirs <- setup_analysis_directories(base_dir, step_name)
  
  # Setup logging
  log_file <- file.path(base_dir, paste0(step_name, ".log"))
  if (exists("setup_logging", mode = "function")) {
    setup_logging(log_file)
  }
  
  # Prepare result
  result <- c(list(base_dir = base_dir), dirs, list(log_file = log_file))
  
  # Optionally assign to parent environment
  if (assign_to_parent) {
    parent_env <- parent.frame()
    for (name in names(result)) {
      assign(name, result[[name]], envir = parent_env)
    }
  }
  
  invisible(result)
}


#' Null-coalescing operator
#' 
#' Returns the right-hand side if left-hand side is NULL
#' 
#' @param a Value to check
#' @param b Default value if a is NULL
#' @return a if not NULL, otherwise b
#' @examples
#' NULL %||% "default"  # Returns "default"
#' "value" %||% "default"  # Returns "value"
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}


#' Check if a package is available
#' 
#' Checks if a package is installed without loading it
#' 
#' @param pkg Package name as string
#' @return TRUE if package is installed, FALSE otherwise
#' @examples
#' has_pkg("ggplot2")
#' has_pkg("nonexistent_package")
has_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}


#' Safely load optional packages
#' 
#' Attempts to load one or more packages with suppressed messages.
#' If a package is not available, logs a warning and continues.
#' Useful for optional functionality that should not break the pipeline.
#' 
#' @param pkgs Character vector of package names (or single string)
#' @return Logical: TRUE if ALL packages successfully loaded, FALSE otherwise
#' @examples
#' safe_lib("ggplot2")
#' safe_lib(c("plotly", "htmlwidgets"))
#' if (safe_lib("ggrepel")) { # Use as conditional check
#'   # ggrepel-specific code
#' }
safe_lib <- function(pkgs) {
  all_loaded <- TRUE
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      warning(sprintf("Package '%s' not available - skipping optional functionality", p), 
              call. = FALSE)
      all_loaded <- FALSE
    } else {
      suppressWarnings(suppressMessages(
        require(p, character.only = TRUE, quietly = TRUE)
      ))
    }
  }
  return(all_loaded)
}


#' Get script directory
#' 
#' Returns the directory containing the currently executing R script.
#' Works in both interactive and batch mode (Rscript).
#' 
#' @return Character string with script directory path
#' @examples
#' script_dir <- get_script_dir()
#' source(file.path(script_dir, "utils", "common.R"))
get_script_dir <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
  } else {
    return(getwd())
  }
}


#' Create directory if it doesn't exist
#' 
#' Helper function to ensure directory exists
#' 
#' @param path Directory path to create
#' @param recursive Create parent directories as needed
#' @return Invisible TRUE if successful
ensure_dir <- function(path, recursive = TRUE) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = recursive, showWarnings = FALSE)
  }
  invisible(file.exists(path))
}


#' Format bytes to human-readable size
#' 
#' Converts byte count to KB, MB, GB, etc.
#' 
#' @param bytes Numeric value in bytes
#' @return Character string with formatted size
#' @examples
#' format_bytes(1024)  # "1.00 KB"
#' format_bytes(1048576)  # "1.00 MB"
format_bytes <- function(bytes) {
  # Handle invalid inputs
  if (length(bytes) == 0) return("N/A")
  if (length(bytes) > 1) bytes <- bytes[1]  # Take first element if vector
  if (!is.numeric(bytes)) return("N/A")
  if (is.na(bytes) || bytes < 0) return("N/A")
  if (bytes == 0) return("0 B")
  
  units <- c("B", "KB", "MB", "GB", "TB")
  # Use base::log to avoid collision with custom log() function
  power <- min(floor(base::log(bytes, 1024)), length(units) - 1)
  
  value <- bytes / (1024^power)
  sprintf("%.2f %s", value, units[power + 1])
}


#' Setup standard output directory structure
#' 
#' Creates standardized directory structure for pipeline outputs.
#' Common pattern: base/data, base/plots, base/interactive, base/rds_objects
#' 
#' @param base_dir Base output directory
#' @param subdirs Character vector of subdirectory names (default: standard set)
#' @return Named list of directory paths
#' @examples
#' dirs <- setup_output_dirs("output/02_qc")
#' # Returns: list(base="...", data="...", plots="...", dash="...", rds="...")
setup_output_dirs <- function(base_dir, 
                              subdirs = c("data", "plots", "interactive", "rds_objects")) {
  # Create base directory
  ensure_dir(base_dir)
  
  # Create subdirectories
  dir_paths <- list(base = base_dir)
  
  for (subdir in subdirs) {
    # Map common aliases
    subdir_name <- switch(subdir,
      "interactive" = "dash",
      "rds_objects" = "rds",
      subdir
    )
    
    full_path <- file.path(base_dir, subdir)
    ensure_dir(full_path)
    dir_paths[[subdir_name]] <- full_path
  }
  
  return(dir_paths)
}


#' Find matrix directory for STARsolo output
#' 
#' Searches for STARsolo output matrices in various directory structures.
#' Handles both new (condition/sample) and old (flat) directory layouts.
#' Automatically gzips uncompressed files if R.utils is available.
#' 
#' @param sample_name Full sample name
#' @param root Root directory for STARsolo outputs
#' @param condition Experimental condition (optional)
#' @param sample Original sample name (optional)
#' @return Path to matrix directory or NA if not found
#' @examples
#' matrix_dir <- find_matrix_dir("ctrl_sample1", "/data/starsolo")
find_matrix_dir <- function(sample_name, root, condition = NULL, sample = NULL) {
  # Parse sample_name if condition/sample not provided
  if (is.null(condition) || is.null(sample)) {
    parts <- strsplit(sample_name, "_", fixed = TRUE)[[1]]
    condition <- parts[1]
    sample <- paste(parts[-1], collapse = "_")
  }
  
  # Define search candidates (new structure first, then fallbacks)
  candidates <- c(
    # New structure: root/condition/sample/sample_Solo.out/Gene/filtered
    file.path(root, condition, sample, paste0(sample_name, "_Solo.out"), "Gene", "filtered"),
    file.path(root, condition, sample, paste0(sample_name, "_Solo.out"), "Gene", "raw"),
    # Old flat structure
    file.path(root, sample_name, "Solo.out", "Gene", "filtered"),
    file.path(root, sample_name, "Solo.out", "Gene", "raw"),
    file.path(root, sample_name, "Gene", "filtered_feature_bc_matrix"),
    file.path(root, sample_name, "Gene")
  )
  
  # Check each candidate
  required_files <- c("barcodes.tsv", "features.tsv", "matrix.mtx")
  
  for (candidate_dir in candidates) {
    # Auto-gzip uncompressed files if possible
    for (file in required_files) {
      fpath <- file.path(candidate_dir, file)
      gzpath <- paste0(fpath, ".gz")
      
      if (file.exists(fpath) && !file.exists(gzpath)) {
        if (has_pkg("R.utils")) {
          message(sprintf("Auto-gzipping %s", fpath))
          R.utils::gzip(fpath, destname = gzpath, overwrite = TRUE, remove = TRUE)
        }
      }
    }
    
    # Check if all required files exist (gzipped)
    all_files_present <- all(file.exists(file.path(candidate_dir, paste0(required_files, ".gz"))))
    
    if (all_files_present) {
      if (exists("log", mode = "function")) {
        log("Found matrix for", sample_name, "in:", candidate_dir)
      }
      return(candidate_dir)
    }
  }
  
  # Not found
  if (exists("log_warning", mode = "function")) {
    log_warning("No valid matrix found for", sample_name, "in", root)
  } else {
    warning("No valid matrix found for sample ", sample_name)
  }
  
  return(NA_character_)
}

#' Parse nextflow.config to extract params
#' 
#' Reads and parses nextflow.config file to extract parameter values
#' 
#' @param config_path Path to nextflow.config. If NULL, searches relative to script directory
#' @return Named list of parameter values
#' @examples
#' params <- parse_nextflow_config()
#' species <- params$species
parse_nextflow_config <- function(config_path = NULL) {
  if (is.null(config_path)) {
    # Try to find nextflow.config relative to script directory
    if (exists("get_script_dir", mode = "function")) {
      script_dir <- get_script_dir()
      config_path <- file.path(dirname(script_dir), "nextflow.config")
    } else {
      config_path <- "nextflow.config"
    }
  }
  
  if (!file.exists(config_path)) {
    warning("nextflow.config not found at: ", config_path)
    return(list())
  }
  
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Reading parameters from nextflow.config: ", config_path)
  
  # Read config file
  config_lines <- readLines(config_path, warn = FALSE)
  
  # Extract params block
  in_params <- FALSE
  params <- list()
  
  for (line in config_lines) {
    line <- trimws(line)
    
    # Start of params block
    if (grepl("^params\\s*\\{", line)) {
      in_params <- TRUE
      next
    }
    
    # End of params block
    if (in_params && grepl("^\\}", line)) {
      break
    }
    
    # Parse parameter lines
    if (in_params && grepl("=", line) && !grepl("^//", line)) {
      # Remove comments
      line <- sub("//.*$", "", line)
      line <- trimws(line)
      
      # Parse key = value
      parts <- strsplit(line, "=")[[1]]
      if (length(parts) >= 2) {
        key <- trimws(parts[1])
        value <- trimws(paste(parts[-1], collapse = "="))
        
        # Remove quotes and trailing comma
        value <- gsub('^"', '', value)
        value <- gsub('"$', '', value)
        value <- gsub(',$', '', value)
        value <- trimws(value)
        
        params[[key]] <- value
      }
    }
  }
  
  return(params)
}

#' Parse comma-separated resolution values into numeric vector
#' 
#' @param x String or numeric value(s) to parse
#' @return Numeric vector of resolutions
#' @examples
#' parse_resolutions("0.2,0.4,0.6")  # c(0.2, 0.4, 0.6)
#' parse_resolutions(c(0.2, 0.4))    # c(0.2, 0.4)
parse_resolutions <- function(x) {
  if (is.null(x)) return(numeric())
  if (is.numeric(x)) return(as.numeric(x))
  x <- as.character(x)
  # Remove spaces
  x <- gsub("\\s+", "", x)
  as.numeric(strsplit(x, ",")[[1]])
}

#' Parse PC dimension specification into integer vector
#' 
#' Handles multiple formats:
#' - "30" means first 30 PCs (1:30)
#' - "1:30" means range from 1 to 30
#' - "1,2,5" means specific PCs
#' - "10-50" means range from 10 to 50 (also parses hyphen)
#' 
#' @param x String or numeric value to parse
#' @return Integer vector of PC indices
#' @examples
#' parse_dims("30")       # 1:30
#' parse_dims("10:50")    # 10:50
#' parse_dims("10-50")    # 10:50
#' parse_dims("1,2,5")    # c(1, 2, 5)
parse_dims <- function(x) {
  if (is.null(x)) return(1:30)
  if (is.numeric(x)) return(as.integer(x))
  
  x <- gsub("\\s+", "", as.character(x))
  
  # Handle range notation with colon: "1:30" or "10:50"
  if (grepl(":", x)) {
    parts <- strsplit(x, ":", fixed = TRUE)[[1]]
    a <- as.integer(parts[1])
    b <- as.integer(parts[2])
    if (is.na(a) || is.na(b) || a < 1 || b < a) return(1:30)
    return(seq(a, b))
  }
  
  # Handle range notation with hyphen: "10-50" (as seen in config)
  if (grepl("-", x) && !grepl(",", x)) {
    parts <- strsplit(x, "-", fixed = TRUE)[[1]]
    a <- as.integer(parts[1])
    b <- as.integer(parts[2])
    if (is.na(a) || is.na(b) || a < 1 || b < a) return(1:30)
    return(seq(a, b))
  }
  
  # Handle comma-separated list: "1,2,5,10"
  if (grepl(",", x)) {
    v <- suppressWarnings(as.integer(strsplit(x, ",")[[1]]))
    v <- v[!is.na(v)]
    if (length(v) == 0) return(1:30)
    return(v)
  }
  
  # Handle single number: "30" means 1:30
  n <- suppressWarnings(as.integer(x))
  if (is.na(n) || n < 1) return(1:30)
  return(1:n)
}

# ==============================================================================
# Logging Utilities
# ==============================================================================
# Centralized logging functions for pipeline execution tracking
# 
# Functions:
# - setup_logging(): Initialize log file
# - log(): Write timestamped message to log
# - log_section(): Write formatted section header
# - log_error(): Write error message
# ==============================================================================

#' Setup logging for an analysis module
#' 
#' Creates log file and directory structure
#' 
#' @param log_path Full path to log file
#' @return Path to created log file
#' @examples
#' log_file <- setup_logging("output/logs/analysis.log")
setup_logging <- function(log_path) {
  # Create log directory if needed
  log_dir <- dirname(log_path)
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Initialize log file with header
  cat("==============================================\n", file = log_path, append = FALSE)
  cat(sprintf("Analysis started: %s\n", Sys.time()), file = log_path, append = TRUE)
  cat("==============================================\n\n", file = log_path, append = TRUE)
  
  return(log_path)
}


#' Write timestamped log message
#' 
#' Writes message to both console and log file with timestamp.
#' Requires log_file variable to be defined in calling environment.
#' 
#' @param ... Message components (will be pasted together)
#' @return Invisible NULL
#' @examples
#' log_file <- "test.log"
#' log("Processing started")
#' log("Found", 100, "cells")
log <- function(...) {
  message <- paste(..., collapse = " ")
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] %s\n", timestamp, message)
  
  # Print to console
  cat(log_entry)
  
  # Write to log file if it exists in parent environment
  if (exists("log_file", envir = parent.frame(), inherits = FALSE)) {
    log_file_path <- get("log_file", envir = parent.frame())
    cat(log_entry, file = log_file_path, append = TRUE)
  }
  
  invisible(NULL)
}


#' Write section header to log
#' 
#' Creates a visually distinct section break in logs
#' 
#' @param title Section title
#' @param char Character to use for border (default: "=")
#' @param width Width of border line
#' @return Invisible NULL
#' @examples
#' log_section("Data Loading")
#' log_section("QC Analysis", char = "-")
log_section <- function(title, char = "=", width = 80) {
  border <- paste(rep(char, width), collapse = "")
  log(border)
  log(toupper(title))
  log(border)
  invisible(NULL)
}


#' Write error message to log
#' 
#' Logs error with special formatting
#' 
#' @param ... Error message components
#' @return Invisible NULL
#' @examples
#' log_error("Failed to load file:", "data.rds")
log_error <- function(...) {
  log("ERROR:", ...)
  invisible(NULL)
}


#' Write warning message to log
#' 
#' Logs warning with special formatting
#' 
#' @param ... Warning message components
#' @return Invisible NULL
#' @examples
#' log_warning("Missing optional package:", "ggrepel")
log_warning <- function(...) {
  log("WARNING:", ...)
  invisible(NULL)
}


#' Write success message to log
#' 
#' Logs success with special formatting
#' 
#' @param ... Success message components
#' @return Invisible NULL
#' @examples
#' log_success("Analysis completed successfully")
log_success <- function(...) {
  log("SUCCESS:", ...)
  invisible(NULL)
}

# Bootstrap utility for finding script directory
# This must be sourced before any other utils files
# Usage: source("utils/bootstrap.R"); script_dir <- get_script_dir()

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

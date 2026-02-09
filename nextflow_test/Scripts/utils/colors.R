# ==============================================================================
# Color Palette Utilities
# ==============================================================================
# Functions for generating consistent color schemes across plots
# 
# Functions:
# - generate_colors(): Smart color palette for any number of categories
# - get_condition_colors(): Specific colors for experimental conditions
# - get_celltype_colors(): Large palette for cell type annotation
# ==============================================================================

#' Generate color palette for any number of categories
#' 
#' Intelligently selects color palette based on number of categories.
#' Uses RColorBrewer palettes when available, falls back to viridis/rainbow.
#' 
#' @param n Number of colors needed
#' @param palette Palette name (default: "auto" for automatic selection)
#' @return Character vector of hex colors
#' @examples
#' generate_colors(5)
#' generate_colors(20)
#' generate_colors(50)
generate_colors <- function(n, palette = "auto") {
  if (n <= 0) return(character(0))
  if (n == 1) return("#3B9AB2")  # Single nice blue
  
  # Try using RColorBrewer for small-medium numbers
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    if (n <= 12) {
      # Use Set3 for up to 12 categories
      colors <- RColorBrewer::brewer.pal(max(3, n), "Set3")
      return(colors[1:n])
      
    } else if (n <= 24) {
      # Combine Set3 and Paired for up to 24 categories
      c1 <- RColorBrewer::brewer.pal(12, "Set3")
      c2 <- RColorBrewer::brewer.pal(12, "Paired")
      return(c(c1, c2)[1:n])
      
    } else if (n <= 33) {
      # Add Set1 for up to 33 categories
      c1 <- RColorBrewer::brewer.pal(12, "Set3")
      c2 <- RColorBrewer::brewer.pal(12, "Paired")
      c3 <- RColorBrewer::brewer.pal(9, "Set1")
      return(c(c1, c2, c3)[1:n])
      
    } else {
      # For very large numbers, interpolate
      base_colors <- c(
        RColorBrewer::brewer.pal(12, "Set3"),
        RColorBrewer::brewer.pal(12, "Paired"),
        RColorBrewer::brewer.pal(9, "Set1")
      )
      color_func <- grDevices::colorRampPalette(base_colors)
      return(color_func(n))
    }
  }
  
  # Fallback: try viridis
  if (requireNamespace("viridis", quietly = TRUE)) {
    if (n <= 20) {
      return(viridis::turbo(n))
    } else {
      return(viridis::viridis(n))
    }
  }
  
  # Last resort: rainbow
  return(grDevices::rainbow(n))
}


#' Get colors for experimental conditions
#' 
#' Provides consistent colors for common experimental conditions
#' 
#' @param conditions Character vector of condition names
#' @return Named vector of hex colors
#' @examples
#' get_condition_colors(c("control", "treatment", "knockout"))
get_condition_colors <- function(conditions) {
  # Predefined colors for common condition names
  color_map <- c(
    "control" = "#808080",      # Gray
    "ctrl" = "#808080",
    "wt" = "#808080",
    "wildtype" = "#808080",
    
    "treatment" = "#E41A1C",    # Red
    "treated" = "#E41A1C",
    "drug" = "#E41A1C",
    
    "knockout" = "#377EB8",     # Blue
    "ko" = "#377EB8",
    "mutant" = "#377EB8",
    
    "overexpression" = "#4DAF4A",  # Green
    "oe" = "#4DAF4A",
    
    "knockdown" = "#984EA3",    # Purple
    "kd" = "#984EA3",
    "shrna" = "#984EA3"
  )
  
  # Match conditions to predefined colors
  colors <- sapply(tolower(conditions), function(cond) {
    if (cond %in% names(color_map)) {
      return(color_map[cond])
    } else {
      return(NA)
    }
  })
  
  # For unmatched conditions, assign from palette
  if (any(is.na(colors))) {
    n_unmatched <- sum(is.na(colors))
    extra_colors <- generate_colors(n_unmatched)
    colors[is.na(colors)] <- extra_colors
  }
  
  names(colors) <- conditions
  return(unname(colors))
}


#' Get large color palette for cell types
#' 
#' Specialized palette for cell type annotations (supports 50+ types)
#' Uses IGV palette from ggsci if available
#' 
#' @param n Number of cell types
#' @return Character vector of hex colors
#' @examples
#' get_celltype_colors(30)
get_celltype_colors <- function(n) {
  if (n <= 0) return(character(0))
  
  # Try using ggsci IGV palette (best for many categories)
  if (requireNamespace("ggsci", quietly = TRUE)) {
    igv_colors <- ggsci::pal_igv()(51)
    if (n <= 51) {
      return(igv_colors[1:n])
    } else {
      # Extend with viridis
      if (requireNamespace("viridis", quietly = TRUE)) {
        extra <- viridis::turbo(n - 51)
        return(c(igv_colors, extra))
      } else {
        # Interpolate IGV colors
        color_func <- grDevices::colorRampPalette(igv_colors)
        return(color_func(n))
      }
    }
  }
  
  # Fallback to standard palette
  return(generate_colors(n))
}


#' Create color scale for continuous values
#' 
#' Returns ggplot2 scale for continuous color mapping
#' 
#' @param palette Palette name: "viridis", "plasma", "magma", "inferno", "cividis"
#' @param direction Color direction: 1 or -1 (default: 1)
#' @return ggplot2 scale object
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt, color = hp)) + 
#'   geom_point() + 
#'   scale_color_continuous_custom("viridis")
scale_color_continuous_custom <- function(palette = "viridis", direction = 1) {
  if (requireNamespace("viridis", quietly = TRUE)) {
    return(viridis::scale_color_viridis_c(option = palette, direction = direction))
  } else {
    return(ggplot2::scale_color_gradient(low = "blue", high = "red"))
  }
}


#' Create named color vector for categories
#' 
#' Convenience function to create named color vectors
#' 
#' @param categories Character vector of category names
#' @param colors Color vector (if NULL, generates automatically)
#' @return Named character vector of colors
#' @examples
#' make_color_vector(c("A", "B", "C"))
make_color_vector <- function(categories, colors = NULL) {
  if (is.null(colors)) {
    colors <- generate_colors(length(categories))
  }
  
  if (length(colors) != length(categories)) {
    stop("Number of colors must match number of categories")
  }
  
  names(colors) <- categories
  return(colors)
}

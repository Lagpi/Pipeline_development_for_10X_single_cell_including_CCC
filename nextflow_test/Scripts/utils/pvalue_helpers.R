#===============================================================================
# pvalue_helpers.R
# Helper functions to add p-values and statistical annotations to plots
#===============================================================================

# Add p-value annotation to ggplot2 plots
add_pvalue_annotation <- function(p, pval, x_pos = NULL, y_pos = NULL, label_size = 4) {
  if (is.na(pval) || is.null(pval)) return(p)
  
  # Format p-value
  if (pval < 0.001) {
    pval_label <- "p < 0.001 ***"
  } else if (pval < 0.01) {
    pval_label <- sprintf("p = %.3f **", pval)
  } else if (pval < 0.05) {
    pval_label <- sprintf("p = %.3f *", pval)
  } else {
    pval_label <- sprintf("p = %.3f", pval)
  }
  
  # Get default positions if not provided
  if (is.null(x_pos) || is.null(y_pos)) {
    # Use top-right corner by default
    x_pos <- Inf
    y_pos <- Inf
  }
  
  p <- p + ggplot2::annotate("text", 
                            x = x_pos, y = y_pos,
                            label = pval_label,
                            hjust = 1.1, vjust = 1.1,
                            size = label_size,
                            fontface = "italic",
                            color = "darkred")
  
  return(p)
}

# Format p-value for display
format_pvalue <- function(pval) {
  if (is.na(pval) || is.null(pval)) return("p = NA")
  
  if (pval < 0.001) {
    return("p < 0.001 ***")
  } else if (pval < 0.01) {
    return(sprintf("p = %.3f **", pval))
  } else if (pval < 0.05) {
    return(sprintf("p = %.3f *", pval))
  } else {
    return(sprintf("p = %.3f", pval))
  }
}

# Add p-value and effect size to plot title or subtitle
add_stats_to_title <- function(title, pval, effect_size = NULL, method = NULL) {
  stats_parts <- c()
  
  if (!is.null(method)) {
    stats_parts <- c(stats_parts, method)
  }
  
  if (!is.null(pval)) {
    stats_parts <- c(stats_parts, format_pvalue(pval))
  }
  
  if (!is.null(effect_size)) {
    stats_parts <- c(stats_parts, sprintf("ES = %.2f", effect_size))
  }
  
  if (length(stats_parts) > 0) {
    stats_str <- paste(stats_parts, collapse = " | ")
    return(sprintf("%s\n%s", title, stats_str))
  }
  
  return(title)
}

# Add significance stars to data
add_significance_stars <- function(pval) {
  if (is.na(pval) || is.null(pval)) return("")
  if (pval < 0.001) return("***")
  if (pval < 0.01) return("**")
  if (pval < 0.05) return("*")
  return("ns")
}

# Create p-value annotation layer for plotly plots
create_pvalue_plotly_annotation <- function(pval, x_pos = NULL, y_pos = NULL) {
  if (is.na(pval) || is.null(pval)) return(NULL)
  
  pval_label <- format_pvalue(pval)
  
  # Default positions (top-right)
  if (is.null(x_pos)) x_pos <- "right"
  if (is.null(y_pos)) y_pos <- "top"
  
  list(
    text = pval_label,
    xref = "paper",
    yref = "paper",
    x = if (x_pos == "right") 0.99 else if (x_pos == "left") 0.01 else 0.5,
    y = if (y_pos == "top") 0.99 else if (y_pos == "bottom") 0.01 else 0.5,
    xanchor = if (x_pos == "right") "right" else if (x_pos == "left") "left" else "center",
    yanchor = if (y_pos == "top") "top" else if (y_pos == "bottom") "bottom" else "middle",
    showarrow = FALSE,
    font = list(size = 12, color = "darkred", family = "Arial Black")
  )
}

# Add p-values to comparison plots (e.g., boxplots, violins)
add_comparison_pvalues <- function(data, x_col, y_col, test_name = "Wilcoxon", 
                                   group1 = NULL, group2 = NULL) {
  result <- list()
  
  if (is.null(group1) || is.null(group2)) {
    # Auto-detect groups
    groups <- unique(data[[x_col]])
    if (length(groups) >= 2) {
      group1 <- groups[1]
      group2 <- groups[2]
    } else {
      return(result)
    }
  }
  
  # Get values for each group
  vals1 <- data[[y_col]][data[[x_col]] == group1]
  vals2 <- data[[y_col]][data[[x_col]] == group2]
  
  # Remove NAs
  vals1 <- vals1[!is.na(vals1)]
  vals2 <- vals2[!is.na(vals2)]
  
  if (length(vals1) == 0 || length(vals2) == 0) {
    return(result)
  }
  
  # Perform test
  test_result <- tryCatch({
    if (tolower(test_name) == "wilcoxon") {
      wilcox.test(vals1, vals2)
    } else if (tolower(test_name) == "ttest") {
      t.test(vals1, vals2)
    } else if (tolower(test_name) == "mann-whitney") {
      wilcox.test(vals1, vals2)
    } else {
      wilcox.test(vals1, vals2)
    }
  }, error = function(e) {
    return(list(p.value = NA, statistic = NA))
  })
  
  # Calculate effect size (Cohen's d)
  mean_diff <- mean(vals1, na.rm = TRUE) - mean(vals2, na.rm = TRUE)
  pooled_sd <- sqrt(((length(vals1) - 1) * sd(vals1, na.rm = TRUE)^2 + 
                     (length(vals2) - 1) * sd(vals2, na.rm = TRUE)^2) / 
                    (length(vals1) + length(vals2) - 2))
  cohens_d <- mean_diff / pooled_sd
  
  result$pvalue <- test_result$p.value
  result$statistic <- test_result$statistic
  result$effect_size <- cohens_d
  result$test_name <- test_name
  result$group1 <- group1
  result$group2 <- group2
  
  return(result)
}

# Generate summary statistics for hover text
generate_stats_hover <- function(pval, n, method = "test") {
  hover_text <- sprintf("%s\np-value: %s\nn: %d", 
                       method,
                       format_pvalue(pval),
                       n)
  return(hover_text)
}

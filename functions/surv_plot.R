create_gene_survival_plot <- function(expression_data,
                                      predictions,
                                      meta_data,
                                      gene_name,
                                      subtypes_include = NULL,
                                      expression_cutoff = 0.5,
                                      time_column = "os",
                                      event_column = "censOS",
                                      high_color = "#E31A1C",
                                      low_color = "#1F78B4",
                                      title = NULL,
                                      risk_table = TRUE,
                                      pval = TRUE,
                                      conf_int = FALSE,
                                      output_file = NULL,
                                      width = 10,
                                      height = 8) {
  
  # Check if gene exists in expression data
  if (!gene_name %in% rownames(expression_data)) {
    stop(paste("Gene", gene_name, "not found in expression data."))
  }
  
  # Check required columns in meta_data
  if (!time_column %in% colnames(meta_data)) {
    stop(paste("Time column", time_column, "not found in meta_data."))
  }
  if (!event_column %in% colnames(meta_data)) {
    stop(paste("Event column", event_column, "not found in meta_data."))
  }
  if (!"sample_id" %in% colnames(meta_data)) {
    stop("sample_id column not found in meta_data.")
  }
  
  # Get expression values for the gene - ensure it's a numeric vector
  gene_expression <- expression_data[gene_name, ]
  
  # Convert to numeric vector if needed
  if (is.list(gene_expression) || is.data.frame(gene_expression)) {
    gene_expression <- as.numeric(unlist(gene_expression))
    names(gene_expression) <- colnames(expression_data)
  } else if (is.matrix(gene_expression)) {
    gene_expression <- as.vector(gene_expression)
    names(gene_expression) <- colnames(expression_data)
  }
  
  # Ensure we have sample names
  if (is.null(names(gene_expression))) {
    names(gene_expression) <- colnames(expression_data)
  }
  
  cat("Gene expression converted to numeric vector of length:", length(gene_expression), "\n")
  
  # Get sample names
  expr_samples <- names(gene_expression)
  pred_samples <- names(predictions)
  
  # Find common samples
  common_samples <- intersect(expr_samples, pred_samples)
  cat("Common samples between expression and predictions:", length(common_samples), "\n")
  
  if (length(common_samples) == 0) {
    stop("No common samples found between expression data and predictions.")
  }
  
  # Create data frame with sample info - only use common samples
  sample_data <- data.frame(
    sample_id = common_samples,
    expression = gene_expression[common_samples],
    subtype = predictions[common_samples],
    stringsAsFactors = FALSE
  )
  
  cat("Sample data created with", nrow(sample_data), "samples\n")
  
  # Filter by subtypes if specified
  if (!is.null(subtypes_include)) {
    sample_data <- sample_data %>%
      filter(subtype %in% subtypes_include)
    
    cat("After subtype filtering:", nrow(sample_data), "samples\n")
    
    if (nrow(sample_data) == 0) {
      stop("No samples found for specified subtypes.")
    }
  }
  
  # Merge with survival data and clean thoroughly
  survival_data <- sample_data %>%
    inner_join(meta_data, by = "sample_id") %>%
    filter(!is.na(get(time_column)), 
           !is.na(get(event_column)), 
           !is.na(expression),
           get(time_column) > 0) %>%  # Remove zero or negative times
    mutate(
      time_var = as.numeric(get(time_column)),
      event_var = as.numeric(get(event_column))
    ) %>%
    filter(!is.na(time_var), !is.na(event_var)) %>%
    arrange(sample_id)  # Ensure consistent ordering
  
  cat("Final survival data:", nrow(survival_data), "samples\n")
  
  if (nrow(survival_data) == 0) {
    stop("No samples with complete survival and expression data found.")
  }
  
  # Create expression groups based on cutoff
  cutoff_value <- quantile(survival_data$expression, expression_cutoff, na.rm = TRUE)
  
  survival_data <- survival_data %>%
    mutate(
      expression_group = ifelse(expression >= cutoff_value, "High", "Low"),
      expression_group = factor(expression_group, levels = c("Low", "High"))
    )
  
  # Check if we have both groups
  group_counts <- table(survival_data$expression_group)
  cat("Expression groups - Low:", group_counts["Low"], "High:", group_counts["High"], "\n")
  
  if (length(unique(survival_data$expression_group)) < 2) {
    warning("Only one expression group found. All samples may have the same expression value.")
    return(NULL)
  }
  
  # Ensure minimum group sizes
  if (min(group_counts) < 3) {
    warning("One or both groups have very few samples (< 3). Results may be unreliable.")
  }
  
  # Create survival object
  fit <- survfit(Surv(time_var, event_var) ~ expression_group, data = survival_data)
  
  # Create title
  if (is.null(title)) {
    subtype_text <- if (is.null(subtypes_include)) {
      "All Subtypes"
    } else {
      paste(subtypes_include, collapse = ", ")
    }
    title <- paste0("Kaplan-Meier: ", gene_name, " Expression\n(", subtype_text, ")")
  }
  
  # Create survival plot with error handling
  tryCatch({
    p <- ggsurvplot(
      fit,
      data = survival_data,
      pval = pval,
      conf.int = conf_int,
      risk.table = risk_table,
      risk.table.col = "strata",
      linetype = "strata",
      surv.median.line = "hv",
      ggtheme = theme_bw(),
      palette = c(low_color, high_color),
      title = title,
      xlab = paste("Time (", time_column, ")"),
      ylab = "Overall Survival Probability",
      legend.title = paste(gene_name, "Expression"),
      legend.labs = c("Low", "High"),
      font.main = c(14, "bold"),
      font.x = c(12, "plain"),
      font.y = c(12, "plain"),
      font.tickslab = c(10, "plain"),
      font.legend = c(10, "plain")
    )
  }, error = function(e) {
    cat("Error in ggsurvplot:", e$message, "\n")
    cat("Trying simplified plot...\n")
    
    # Try a simpler version without some formatting options
    p <<- ggsurvplot(
      fit,
      data = survival_data,
      pval = pval,
      palette = c(low_color, high_color),
      title = title,
      legend.labs = c("Low", "High")
    )
  })
  
  # Perform log-rank test
  logrank_test <- survdiff(Surv(time_var, event_var) ~ expression_group, data = survival_data)
  logrank_pvalue <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
  
  # Get summary statistics
  summary_stats <- list(
    n_total = nrow(survival_data),
    n_low = sum(survival_data$expression_group == "Low"),
    n_high = sum(survival_data$expression_group == "High"),
    cutoff_value = cutoff_value,
    cutoff_percentile = expression_cutoff,
    logrank_pvalue = logrank_pvalue,
    subtypes_analyzed = if (is.null(subtypes_include)) unique(survival_data$subtype) else subtypes_include,
    median_survival = tryCatch(surv_median(fit), error = function(e) NULL)
  )
  
  # Print summary
  cat("Gene:", gene_name, "\n")
  cat("Subtypes analyzed:", paste(summary_stats$subtypes_analyzed, collapse = ", "), "\n")
  cat("Total samples:", summary_stats$n_total, "\n")
  cat("Low expression (n =", summary_stats$n_low, "):", "< ", round(cutoff_value, 3), "\n")
  cat("High expression (n =", summary_stats$n_high, "):", "≥ ", round(cutoff_value, 3), "\n")
  cat("Log-rank test p-value:", formatC(logrank_pvalue, format = "e", digits = 3), "\n")
  
  # Save plot if requested
  if (!is.null(output_file)) {
    tryCatch({
      ggsave(output_file, plot = p$plot, width = width, height = height)
      message(paste("Plot saved to:", output_file))
    }, error = function(e) {
      message("Could not save plot: ", e$message)
    })
  }
  
  # Return results
  return(list(
    survival_fit = fit,
    plot = p,
    data = survival_data,
    summary_stats = summary_stats,
    logrank_test = logrank_test
  ))
}
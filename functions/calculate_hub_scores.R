library(ComplexHeatmap)
library(circlize)
library(dplyr)

#' Calculate Hub Scores and Create Heatmap (Expression-Based Specificity)
#'
#' @param data Data frame containing gene network metrics and expression columns
#' @param system Character string: "5" or "7" for 5-class or 7-class system
#' @param top_n Integer: Number of top hub genes to display in heatmap
#' @param degree_weight Numeric: Weight for degree centrality (default 0.4)
#' @param betweenness_weight Numeric: Weight for betweenness (default 0.3)
#' @param eigenvector_weight Numeric: Weight for eigenvector centrality (default 0.3)
#' @param specificity_method Character: "expression_diff", "fold_change", or "zscore"
#' @param order_by Character: How to order columns - "alphabetical", "hub_score", or "none" (clustering)
#' @param color_low Character: Color for low hub scores
#' @param color_mid Character: Color for mid hub scores
#' @param color_high Character: Color for high hub scores
#' @param title Character: Custom title for the heatmap (NULL for default)
#' @param fontsize_row Numeric: Font size for row names
#' @param fontsize_col Numeric: Font size for column names
#' @param output_file Character: Path to save PDF (NULL for no save)
#' @param width Numeric: Width of saved PDF in inches
#' @param height Numeric: Height of saved PDF in inches
#' 
#' @return List containing hub_scores data frame and heatmap object

calculate_hub_scores <- function(data,
                                 system = c("7", "5"),
                                 top_n = 34,
                                 degree_weight = 0.4,
                                 betweenness_weight = 0.3,
                                 eigenvector_weight = 0.3,
                                 specificity_method = c("expression_diff", "fold_change", "zscore"),
                                 order_by = c("alphabetical", "hub_score", "none"),
                                 color_low = "#FFF5E6",
                                 color_mid = "#E85D04",
                                 color_high = "#370617",
                                 title = NULL,
                                 fontsize_row = 10,
                                 fontsize_col = 11,
                                 output_file = NULL,
                                 width = 10,
                                 height = 8) {
  
  # Match arguments
  system <- match.arg(system)
  order_by <- match.arg(order_by)
  specificity_method <- match.arg(specificity_method)
  
  # Validate weights sum to 1
  if (abs(degree_weight + betweenness_weight + eigenvector_weight - 1) > 0.01) {
    warning("Weights do not sum to 1. Normalizing weights.")
    total <- degree_weight + betweenness_weight + eigenvector_weight
    degree_weight <- degree_weight / total
    betweenness_weight <- betweenness_weight / total
    eigenvector_weight <- eigenvector_weight / total
  }
  
  # Determine expression columns and subtypes based on system
  if (system == "7") {
    expr_cols <- c("mean_expr_UroA", "mean_expr_UroB", "mean_expr_UroC", 
                   "mean_expr_GU", "mean_expr_BaSq")
    subtypes <- c("UroA", "UroB", "UroC", "GU", "BaSq")
    default_title <- "Subtype-Specific Hub Gene Scores (7-Class System)"
  } else if (system == "5") {
    expr_cols <- c("mean_expr_5_Uro", "mean_expr_5_GU", "mean_expr_5_BaSq")
    subtypes <- c("Uro", "GU", "BaSq")
    default_title <- "Subtype-Specific Hub Gene Scores (5-Class System)"
  }
  
  # Use default title if not provided
  if (is.null(title)) {
    title <- default_title
  }
  
  # Check required columns exist
  required_cols <- c("name", "degree", "betweenness", "eigenvector")
  missing_cols <- setdiff(c(required_cols, expr_cols), colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Calculate global statistics for each gene
  hub_scores <- data %>%
    rowwise() %>%
    mutate(
      global_mean_expr = mean(c_across(all_of(expr_cols)), na.rm = TRUE),
      expr_sd = sd(c_across(all_of(expr_cols)), na.rm = TRUE),
      expr_sd = ifelse(is.na(expr_sd) | expr_sd == 0, 0.01, expr_sd)  # Avoid division by zero
    ) %>%
    ungroup()
  
  # Normalize network metrics to 0-1 scale
  hub_scores <- hub_scores %>%
    mutate(
      norm_degree = (degree - min(degree)) / (max(degree) - min(degree)),
      norm_betweenness = (betweenness - min(betweenness)) / (max(betweenness) - min(betweenness)),
      norm_eigenvector = eigenvector
    )
  
  # Calculate specificity scores based on expression values
  for (i in seq_along(subtypes)) {
    subtype <- subtypes[i]
    target_col <- expr_cols[i]
    other_cols <- expr_cols[-i]
    
    if (specificity_method == "expression_diff") {
      # Method 1: Difference from mean of others, normalized by gene variability
      hub_scores <- hub_scores %>%
        rowwise() %>%
        mutate(
          !!paste0(subtype, "_diff") := .data[[target_col]] - 
            mean(c_across(all_of(other_cols)), na.rm = TRUE),
          # Normalize by gene's standard deviation (z-score like)
          !!paste0(subtype, "_spec_score_raw") := 
            .data[[paste0(subtype, "_diff")]] / expr_sd
        ) %>%
        ungroup() %>%
        mutate(
          # Scale to 0-1 range (clip negative values to 0)
          !!paste0(subtype, "_spec_score") := pmax(0, .data[[paste0(subtype, "_spec_score_raw")]]),
          # Normalize to 0-1 using all genes' scores
          !!paste0(subtype, "_spec_score") := 
            (.data[[paste0(subtype, "_spec_score")]] - 
               min(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE)) /
            (max(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE) - 
               min(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE) + 0.001)
        )
      
    } else if (specificity_method == "fold_change") {
      # Method 2: Fold-change relative to mean of others
      hub_scores <- hub_scores %>%
        rowwise() %>%
        mutate(
          other_mean = mean(c_across(all_of(other_cols)), na.rm = TRUE),
          # Add small constant to avoid log(0)
          !!paste0(subtype, "_fc") := (.data[[target_col]] + 1) / (other_mean + 1),
          !!paste0(subtype, "_spec_score_raw") := log2(.data[[paste0(subtype, "_fc")]])
        ) %>%
        ungroup() %>%
        select(-other_mean) %>%
        mutate(
          # Clip negative fold-changes to 0
          !!paste0(subtype, "_spec_score") := pmax(0, .data[[paste0(subtype, "_spec_score_raw")]]),
          # Normalize to 0-1
          !!paste0(subtype, "_spec_score") := 
            (.data[[paste0(subtype, "_spec_score")]] - 
               min(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE)) /
            (max(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE) - 
               min(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE) + 0.001)
        )
      
    } else if (specificity_method == "zscore") {
      # Method 3: Z-score relative to global mean
      hub_scores <- hub_scores %>%
        mutate(
          !!paste0(subtype, "_zscore") := 
            (.data[[target_col]] - global_mean_expr) / expr_sd,
          !!paste0(subtype, "_spec_score") := pmax(0, .data[[paste0(subtype, "_zscore")]]),
          # Normalize to 0-1
          !!paste0(subtype, "_spec_score") := 
            (.data[[paste0(subtype, "_spec_score")]] - 
               min(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE)) /
            (max(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE) - 
               min(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE) + 0.001)
        )
    }
    
    # Calculate hub score: network centrality × specificity
    hub_scores <- hub_scores %>%
      mutate(
        !!paste0(subtype, "_hub_score") := (
          degree_weight * norm_degree + 
            betweenness_weight * norm_betweenness + 
            eigenvector_weight * norm_eigenvector
        ) * .data[[paste0(subtype, "_spec_score")]]
      )
  }
  
  # Get top hubs overall
  hub_score_cols <- paste0(subtypes, "_hub_score")
  
  top_hubs <- hub_scores %>%
    mutate(max_hub_score = pmax(!!!syms(hub_score_cols), na.rm = TRUE)) %>%
    arrange(desc(max_hub_score)) %>%
    slice_head(n = top_n)
  
  # Create matrix of hub scores
  hub_matrix <- top_hubs %>%
    dplyr::select(name, all_of(hub_score_cols)) %>%
    column_to_rownames("name") %>%
    as.matrix()
  
  colnames(hub_matrix) <- subtypes
  
  # TRANSPOSE: Make subtypes as rows, genes as columns
  hub_matrix <- t(hub_matrix)
  
  # Create color vector for gene names based on is_seed status
  gene_colors <- top_hubs$is_seed
  names(gene_colors) <- top_hubs$name
  
  # Set colors: FALSE (not seed) = red, TRUE (seed) = black
  gene_name_colors <- ifelse(gene_colors, "black", "red")
  names(gene_name_colors) <- names(gene_colors)
  
  # Determine column order (genes)
  if (order_by == "alphabetical") {
    col_order_param <- order(colnames(hub_matrix))
    cluster_cols_param <- FALSE
    # Reorder gene colors to match column order
    gene_name_colors <- gene_name_colors[colnames(hub_matrix)[col_order_param]]
  } else if (order_by == "hub_score") {
    # Reverse order so highest hub score comes first
    col_order_param <- seq_len(ncol(hub_matrix))
    cluster_cols_param <- FALSE
    # Gene colors already in correct order since top_hubs is sorted by hub score
    gene_name_colors <- gene_name_colors[colnames(hub_matrix)]
  } else {
    col_order_param <- NULL
    cluster_cols_param <- TRUE
    # For clustering, we'll need to handle colors after heatmap creation
    gene_name_colors <- gene_name_colors[colnames(hub_matrix)]
  }
  
  # Color scheme
  col_fun <- colorRamp2(
    c(0, max(hub_matrix, na.rm = TRUE) / 2, max(hub_matrix, na.rm = TRUE)),
    c(color_low, color_mid, color_high)
  )
  
  # Create heatmap (note: rows are now subtypes, columns are genes)
  ht <- Heatmap(hub_matrix,
                name = "Hub\nScore",
                col = col_fun,
                cluster_rows = FALSE,  # Don't cluster subtypes
                cluster_columns = cluster_cols_param,
                column_order = col_order_param,
                show_row_names = TRUE,
                show_column_names = TRUE,
                rect_gp = gpar(col = "black", lwd = 2),
                row_names_side = "left",  # Row labels on left
                row_names_gp = gpar(fontsize = fontsize_row),
                column_names_gp = gpar(fontsize = fontsize_col, col = gene_name_colors),
                column_title = title,
                heatmap_legend_param = list(title = "Hub Score"))
  
  # Save to file if requested
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
    draw(ht)
    dev.off()
    message(paste("Heatmap saved to:", output_file))
  }
  
  # Return both hub scores and heatmap
  return(list(
    hub_scores = hub_scores,
    top_hubs = top_hubs,
    heatmap = ht,
    method = specificity_method
  ))
}
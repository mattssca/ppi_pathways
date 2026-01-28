library(ComplexHeatmap)
library(circlize)
library(dplyr)

#' Create Expression Rank Heatmap (Transposed: Subtypes as Rows)
#'
#' @param data Data frame containing gene names and rank columns
#' @param system Character string: "5" or "7" for 5-class or 7-class system
#' @param top_n Integer: Number of top genes to display (NULL for all genes)
#' @param order_by Character: How to order columns - "alphabetical", "degree", "hub_score", or "none" (clustering)
#' @param degree_weight Numeric: Weight for degree centrality (default 0.4, used if order_by = "hub_score")
#' @param betweenness_weight Numeric: Weight for betweenness (default 0.3, used if order_by = "hub_score")
#' @param eigenvector_weight Numeric: Weight for eigenvector centrality (default 0.3, used if order_by = "hub_score")
#' @param specificity_method Character: "expression_diff", "fold_change", or "zscore" (used if order_by = "hub_score")
#' @param show_values Logical: Whether to display rank values in cells
#' @param title Character: Custom title for the heatmap (NULL for default)
#' @param fontsize_row Numeric: Font size for row names
#' @param fontsize_col Numeric: Font size for column names
#' @param output_file Character: Path to save PDF (NULL for no save)
#' @param width Numeric: Width of saved PDF in inches
#' @param height Numeric: Height of saved PDF in inches
#' 
#' @return List containing heatmap object and hub_scores (if calculated)

plot_rank_heatmap <- function(data,
                              system = c("7", "5"),
                              top_n = NULL,
                              order_by = c("alphabetical", "degree", "hub_score", "none"),
                              degree_weight = 0.4,
                              betweenness_weight = 0.3,
                              eigenvector_weight = 0.3,
                              specificity_method = c("expression_diff", "fold_change", "zscore"),
                              show_values = FALSE,
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
  
  # Select appropriate rank columns and subtypes based on system
  if (system == "7") {
    rank_pattern <- "^rank_7_"
    subtypes_to_plot <- c("UroA", "UroB", "UroC", "GU", "BaSq")
    rank_cols <- paste0("rank_7_", subtypes_to_plot)
    expr_cols <- c("mean_expr_UroA", "mean_expr_UroB", "mean_expr_UroC", 
                   "mean_expr_GU", "mean_expr_BaSq")
    default_title <- "Gene Expression Ranks Across Subtypes (7-Class)"
    color_breaks <- c(1, 2, 3, 4, 5)
    color_palette <- c("#9E0142", "#F46D43", "#FFFFBF", "#3288BD", "#5E4FA2")
    legend_labels <- c("1 (Highest)", "2", "3", "4", "5 (Lowest)")
  } else if (system == "5") {
    rank_pattern <- "^rank_5_"
    subtypes_to_plot <- c("Uro", "GU", "BaSq")
    rank_cols <- paste0("rank_5_", subtypes_to_plot)
    expr_cols <- c("mean_expr_5_Uro", "mean_expr_5_GU", "mean_expr_5_BaSq")
    default_title <- "Gene Expression Ranks Across Subtypes (5-Class)"
    color_breaks <- c(1, 2, 3)
    color_palette <- c("#9E0142", "#F46D43", "#FFFFBF")
    legend_labels <- c("1 (Highest)", "2", "3 (Lowest)")
  }
  
  # Use default title if not provided
  if (is.null(title)) {
    title <- default_title
  }
  
  # Check if rank columns exist
  missing_cols <- setdiff(rank_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste0("Missing rank columns: ", paste(missing_cols, collapse = ", "),
                "\nAvailable columns: ", paste(colnames(data), collapse = ", ")))
  }
  
  # Calculate hub scores if needed for ordering or filtering
  hub_scores_data <- NULL
  if (order_by == "hub_score" || (!is.null(top_n) && order_by == "hub_score")) {
    
    # Check required columns for hub score calculation
    required_cols <- c("name", "degree", "betweenness", "eigenvector")
    missing_hub_cols <- setdiff(c(required_cols, expr_cols), colnames(data))
    if (length(missing_hub_cols) > 0) {
      stop(paste("Missing columns for hub score calculation:", 
                 paste(missing_hub_cols, collapse = ", ")))
    }
    
    # Validate weights sum to 1
    if (abs(degree_weight + betweenness_weight + eigenvector_weight - 1) > 0.01) {
      warning("Weights do not sum to 1. Normalizing weights.")
      total <- degree_weight + betweenness_weight + eigenvector_weight
      degree_weight <- degree_weight / total
      betweenness_weight <- betweenness_weight / total
      eigenvector_weight <- eigenvector_weight / total
    }
    
    # Calculate global statistics for each gene
    hub_scores_data <- data %>%
      rowwise() %>%
      mutate(
        global_mean_expr = mean(c_across(all_of(expr_cols)), na.rm = TRUE),
        expr_sd = sd(c_across(all_of(expr_cols)), na.rm = TRUE),
        expr_sd = ifelse(is.na(expr_sd) | expr_sd == 0, 0.01, expr_sd)
      ) %>%
      ungroup()
    
    # Normalize network metrics to 0-1 scale
    hub_scores_data <- hub_scores_data %>%
      mutate(
        norm_degree = (degree - min(degree)) / (max(degree) - min(degree)),
        norm_betweenness = (betweenness - min(betweenness)) / (max(betweenness) - min(betweenness)),
        norm_eigenvector = eigenvector
      )
    
    # Calculate specificity scores based on expression values
    for (i in seq_along(subtypes_to_plot)) {
      subtype <- subtypes_to_plot[i]
      target_col <- expr_cols[i]
      other_cols <- expr_cols[-i]
      
      if (specificity_method == "expression_diff") {
        hub_scores_data <- hub_scores_data %>%
          rowwise() %>%
          mutate(
            !!paste0(subtype, "_diff") := .data[[target_col]] - 
              mean(c_across(all_of(other_cols)), na.rm = TRUE),
            !!paste0(subtype, "_spec_score_raw") := 
              .data[[paste0(subtype, "_diff")]] / expr_sd
          ) %>%
          ungroup() %>%
          mutate(
            !!paste0(subtype, "_spec_score") := pmax(0, .data[[paste0(subtype, "_spec_score_raw")]]),
            !!paste0(subtype, "_spec_score") := 
              (.data[[paste0(subtype, "_spec_score")]] - 
                 min(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE)) /
              (max(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE) - 
                 min(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE) + 0.001)
          )
        
      } else if (specificity_method == "fold_change") {
        hub_scores_data <- hub_scores_data %>%
          rowwise() %>%
          mutate(
            other_mean = mean(c_across(all_of(other_cols)), na.rm = TRUE),
            !!paste0(subtype, "_fc") := (.data[[target_col]] + 1) / (other_mean + 1),
            !!paste0(subtype, "_spec_score_raw") := log2(.data[[paste0(subtype, "_fc")]])
          ) %>%
          ungroup() %>%
          select(-other_mean) %>%
          mutate(
            !!paste0(subtype, "_spec_score") := pmax(0, .data[[paste0(subtype, "_spec_score_raw")]]),
            !!paste0(subtype, "_spec_score") := 
              (.data[[paste0(subtype, "_spec_score")]] - 
                 min(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE)) /
              (max(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE) - 
                 min(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE) + 0.001)
          )
        
      } else if (specificity_method == "zscore") {
        hub_scores_data <- hub_scores_data %>%
          mutate(
            !!paste0(subtype, "_zscore") := 
              (.data[[target_col]] - global_mean_expr) / expr_sd,
            !!paste0(subtype, "_spec_score") := pmax(0, .data[[paste0(subtype, "_zscore")]]),
            !!paste0(subtype, "_spec_score") := 
              (.data[[paste0(subtype, "_spec_score")]] - 
                 min(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE)) /
              (max(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE) - 
                 min(.data[[paste0(subtype, "_spec_score")]], na.rm = TRUE) + 0.001)
          )
      }
      
      # Calculate hub score: network centrality × specificity
      hub_scores_data <- hub_scores_data %>%
        mutate(
          !!paste0(subtype, "_hub_score") := (
            degree_weight * norm_degree + 
              betweenness_weight * norm_betweenness + 
              eigenvector_weight * norm_eigenvector
          ) * .data[[paste0(subtype, "_spec_score")]]
        )
    }
    
    # Calculate composite hub score: maximum hub score across all subtypes
    hub_score_cols <- paste0(subtypes_to_plot, "_hub_score")
    hub_scores_data <- hub_scores_data %>%
      mutate(composite_hub_score = pmax(!!!syms(hub_score_cols), na.rm = TRUE))
    
    # Use hub scores data for filtering/ordering
    data <- hub_scores_data
  }
  
  # Filter top genes if requested
  if (!is.null(top_n)) {
    if (order_by == "degree" && "degree" %in% colnames(data)) {
      data <- data %>%
        arrange(desc(degree)) %>%
        slice_head(n = top_n)
    } else if (order_by == "hub_score" && "composite_hub_score" %in% colnames(data)) {
      data <- data %>%
        arrange(desc(composite_hub_score)) %>%
        slice_head(n = top_n)
    } else {
      data <- data %>%
        slice_head(n = top_n)
    }
  }
  
  # Select rank columns
  rank_columns <- data %>%
    dplyr::select(name, all_of(rank_cols))
  
  # Prepare data matrix (excluding the name column)
  rank_matrix <- as.matrix(rank_columns[, -1])
  rownames(rank_matrix) <- rank_columns$name
  colnames(rank_matrix) <- subtypes_to_plot
  
  # TRANSPOSE: Make subtypes as rows, genes as columns
  rank_matrix <- t(rank_matrix)
  
  # Color scheme
  col_fun <- colorRamp2(color_breaks, color_palette)
  
  # Determine column order (genes)
  if (order_by == "alphabetical") {
    col_order_param <- order(colnames(rank_matrix))
    cluster_cols_param <- FALSE
  } else if (order_by == "degree") {
    if ("degree" %in% colnames(data)) {
      # Reverse order so highest degree comes first
      col_order_param <- rev(seq_len(ncol(rank_matrix)))
      cluster_cols_param <- FALSE
    } else {
      warning("'degree' column not found. Using alphabetical order.")
      col_order_param <- order(colnames(rank_matrix))
      cluster_cols_param <- FALSE
    }
  } else if (order_by == "hub_score") {
    if ("composite_hub_score" %in% colnames(data)) {
      # Reverse order so highest composite hub score comes first
      col_order_param <- seq_len(ncol(rank_matrix))
      cluster_cols_param <- FALSE
    } else {
      warning("'composite_hub_score' column not found. Using alphabetical order.")
      col_order_param <- order(colnames(rank_matrix))
      cluster_cols_param <- FALSE
    }
  } else {
    col_order_param <- NULL
    cluster_cols_param <- TRUE
  }
  
  # Create cell function for displaying values
  if (show_values) {
    cell_fun_param <- function(j, i, x, y, width, height, fill) {
      grid.text(rank_matrix[i, j], x, y, gp = gpar(fontsize = 9))
    }
  } else {
    cell_fun_param <- NULL
  }
  
  # Create heatmap (note: rows are now subtypes, columns are genes)
  ht <- Heatmap(rank_matrix,
                name = "Rank",
                col = col_fun,
                cluster_rows = FALSE,  # Don't cluster subtypes
                cluster_columns = cluster_cols_param,
                column_order = col_order_param,
                show_row_names = TRUE,
                show_column_names = TRUE,
                row_names_side = "left",  # Row labels on left
                row_names_gp = gpar(fontsize = fontsize_row),
                column_names_gp = gpar(fontsize = fontsize_col),
                cell_fun = cell_fun_param,
                rect_gp = gpar(col = "black", lwd = 2),
                column_title = title,
                heatmap_legend_param = list(
                  title = "Rank\n(1=Highest\nExpression)",
                  at = color_breaks,
                  labels = legend_labels
                ))
  
  # Save to file if requested
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
    draw(ht)
    dev.off()
    message(paste("Heatmap saved to:", output_file))
  }
  
  # Return heatmap and hub scores (if calculated)
  if (!is.null(hub_scores_data)) {
    return(list(
      heatmap = ht,
      hub_scores = hub_scores_data,
      hub_score_cols = hub_score_cols
    ))
  } else {
    return(ht)
  }
}
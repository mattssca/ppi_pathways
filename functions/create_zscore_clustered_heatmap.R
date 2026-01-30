create_zscore_clustered_heatmap <- function() {
  
  # Prepare data
  network_genes <- network_subtype_analysis$gene
  expr_matrix <- uroscanseq_data$expr_df[network_genes, , drop = FALSE]
  expr_matrix <- apply(expr_matrix, 2, as.numeric)
  rownames(expr_matrix) <- network_genes
  
  # Z-score normalize by rows (genes)
  expr_matrix_zscore <- t(scale(t(expr_matrix)))
  
  # Handle any NAs from genes with zero variance
  expr_matrix_zscore[is.na(expr_matrix_zscore)] <- 0
  
  sample_subtypes <- uroscanseq_data$subtype_5_vector
  
  # Order samples by subtype
  subtype_order <- c("Uro", "GU", "BaSq")
  ordered_samples <- names(sample_subtypes)[order(match(sample_subtypes, subtype_order))]
  expr_matrix_ordered <- expr_matrix_zscore[, ordered_samples]
  
  # Colors
  subtype_colors <- c("Uro" = "#3cb44b", "GU" = "#4363d8", "BaSq" = "#CD2626")
  
  # Column annotation
  col_anno <- HeatmapAnnotation(
    Subtype = sample_subtypes[ordered_samples],
    col = list(Subtype = subtype_colors),
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    simple_anno_size = unit(0.4, "cm")
  )
  
  # Row annotations
  gene_stats <- network_subtype_analysis %>%
    arrange(match(gene, rownames(expr_matrix_ordered)))
  
  row_anno <- rowAnnotation(
    `Disc.\nPower` = gene_stats$max_diff,
    col = list(
      `Disc.\nPower` = colorRamp2(c(0, 0.6, 1.2), c("#FEE08B", "#FC8D59", "#D73027"))
    ),
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    simple_anno_size = unit(0.3, "cm")
  )
  
  # Main heatmap with z-score appropriate color scale
  ht <- Heatmap(
    expr_matrix_ordered,
    name = "Z-score",
    
    # Z-score appropriate color scheme (-3 to +3 standard deviations)
    col = circlize::colorRamp2(c(-2, 0, 2), c("#4DF76F", "black", "#F74D4D")),
    
    # Clustering
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D2",
    
    # Appearance
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 10),
    row_names_side = "left",
    
    # Annotations
    top_annotation = col_anno,
    right_annotation = row_anno,
    
    # Titles
    column_title = "Samples by Molecular Subtype",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_title = "Network Hub Genes\n(Clustered)",
    row_title_gp = gpar(fontsize = 11, fontface = "bold"),
    
    
    # Legend
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9),
      legend_height = unit(4, "cm"),
      border = "black",
      at = c(-3, -2, -1, 0, 1, 2, 3),
      labels = c("-3", "-2", "-1", "0", "1", "2", "3")
    )
  )
  
  return(ht)
}

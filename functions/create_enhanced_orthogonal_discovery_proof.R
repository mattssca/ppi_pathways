create_enhanced_orthogonal_discovery_proof <- function() {
  
  # Panel A: Zero overlap (make more compact)
  overlap_data <- data.frame(
    method = c("Network\nHubs", "Classifier\nDEGs"),
    n_genes = c(10, 725),
    unique_genes = c(10, 725)
  )
  
  panel_a <- ggplot(overlap_data, aes(x = method, y = unique_genes, fill = method)) +
    geom_col(alpha = 0.8, width = 0.6) +
    geom_text(aes(label = paste0("n = ", unique_genes)), vjust = -0.5, fontface = "bold", size = 3.5) +
    scale_fill_manual(values = c("#E31A1C", "#1F78B4")) +
    scale_y_continuous(trans = "log10", labels = scales::comma) +
    labs(
      title = "A. Orthogonal Discovery",
      subtitle = "0% overlap",
      x = "Method",
      y = "Genes (log scale)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "red", face = "bold"),
      axis.text.x = element_text(size = 9),
      axis.title = element_text(size = 9)
    ) +
    annotate("text", x = 1.5, y = 50, label = "0%\nOverlap", 
             size = 4, fontface = "bold", color = "black")
  
  # Panel B: Expression heatmap - ADJUSTED LIMITS
  heatmap_data <- network_subtype_analysis %>%
    select(gene, Uro_mean, GU_mean, BaSq_mean) %>%
    pivot_longer(-gene, names_to = "subtype", values_to = "expression") %>%
    mutate(
      subtype = gsub("_mean", "", subtype),
      subtype = factor(subtype, levels = c("Uro", "GU", "BaSq")),
      # SAME ORDER AS PANEL C: ordered by max_diff, highest at top
      gene = factor(gene, levels = network_subtype_analysis$gene[order(network_subtype_analysis$max_diff)])
    )
  
  # Calculate actual data range for proper limits
  data_range <- range(heatmap_data$expression, na.rm = TRUE)
  data_max <- max(abs(data_range))
  
  panel_b <- ggplot(heatmap_data, aes(x = subtype, y = gene, fill = expression)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient2(
      low = "#4DF76F", mid = "black", high = "#F74D4D",
      midpoint = 0, name = "Mean\nExpression",
      # Use symmetric limits based on actual data range
      limits = c(-0.8, 0.5)  # Adjusted to your data: min ≈ -0.75, max ≈ 0.46
    ) +
    labs(
      title = "B. Subtype Expression Patterns",
      subtitle = "Network genes show clear preferences",
      x = "Subtype",
      y = "Network Hub Genes"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.text.x = element_text(size = 9, face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 9),
      panel.grid = element_blank(),
      legend.key.size = unit(0.4, "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
  
  # Panel C: Discrimination strength - SAME ORDER
  panel_c <- network_subtype_analysis %>%
    mutate(
      # SAME ORDER AS PANEL B: ordered by max_diff for consistency
      gene = fct_reorder(gene, max_diff),
      discrimination_strength = case_when(
        max_diff > 1.0 ~ "Very Strong",
        max_diff > 0.5 ~ "Strong", 
        TRUE ~ "Moderate"
      )
    ) %>%
    ggplot(aes(x = gene, y = max_diff, fill = discrimination_strength)) +
    geom_col(alpha = 0.8) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40", alpha = 0.8) +
    scale_fill_manual(
      values = c("Very Strong" = "#D73027", "Strong" = "#FC8D59", "Moderate" = "#FEE08B"),
      name = "Strength"
    ) +
    coord_flip() +
    labs(
      title = "C. Discrimination Power",
      subtitle = "70% show strong differences",
      x = "Gene",
      y = "Max Expression Difference"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9),
      legend.key.size = unit(0.4, "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
  
  # Combine all three panels
  combined_enhanced <- panel_a + panel_b + panel_c +
    plot_layout(widths = c(1, 1.2, 1.3)) +
    plot_annotation(
      title = "Network Analysis Reveals Subtype Biology Beyond DEG Clustering",
      subtitle = "Pathway-informed hub analysis discovers functionally relevant subtype-defining genes",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5)
      )
    )
  
  return(combined_enhanced)
}

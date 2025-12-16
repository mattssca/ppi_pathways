plot_signature_correlations <- function(correlation_results, 
                                        network_name = "Network",
                                        show_combined = TRUE,
                                        show_summary = TRUE) {
  
  # Create significance and direction categories for stromal
  correlation_results$stromal_category <- case_when(
    correlation_results$stromal_pval >= 0.05 ~ "Not Significant",
    correlation_results$stromal_cor > 0 & correlation_results$stromal_pval < 0.05 ~ "Positive",
    correlation_results$stromal_cor < 0 & correlation_results$stromal_pval < 0.05 ~ "Negative"
  )
  
  # Create significance and direction categories for immune
  correlation_results$immune_category <- case_when(
    correlation_results$immune_pval >= 0.05 ~ "Not Significant",
    correlation_results$immune_cor > 0 & correlation_results$immune_pval < 0.05 ~ "Positive",
    correlation_results$immune_cor < 0 & correlation_results$immune_pval < 0.05 ~ "Negative"
  )
  
  # Stromal signature plot
  p1 <- ggplot(correlation_results, aes(x = stromal_cor, y = 1, color = stromal_category)) +
    geom_point(size = 2, alpha = 0.7, position = position_jitter(height = 0.1)) +
    scale_color_manual(values = c("Not Significant" = "gray70", 
                                  "Positive" = "red", 
                                  "Negative" = "blue")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = paste(network_name, "- Stromal Signature Correlations"),
         subtitle = "Red: Positive correlation, Blue: Negative correlation (p < 0.05)",
         x = "Stromal Correlation",
         y = "",
         color = "Correlation Type") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(0.5, 1.5)
  
  # Immune signature plot
  p2 <- ggplot(correlation_results, aes(x = immune_cor, y = 1, color = immune_category)) +
    geom_point(size = 2, alpha = 0.7, position = position_jitter(height = 0.1)) +
    scale_color_manual(values = c("Not Significant" = "gray70", 
                                  "Positive" = "red", 
                                  "Negative" = "blue")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = paste(network_name, "- Immune Signature Correlations"),
         subtitle = "Red: Positive correlation, Blue: Negative correlation (p < 0.05)",
         x = "Immune Correlation",
         y = "",
         color = "Correlation Type") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(0.5, 1.5)
  
  # Display individual plots
  print(p1)
  print(p2)
  
  # Optional: Combine both plots vertically
  if(show_combined) {
    library(patchwork)
    combined_plot <- p1 / p2
    print(combined_plot)
  }
  
  # Summary tables
  if(show_summary) {
    cat(paste("\n", network_name, "- Stromal Signature Summary:\n"))
    print(table(correlation_results$stromal_category))
    
    cat(paste("\n", network_name, "- Immune Signature Summary:\n"))
    print(table(correlation_results$immune_category))
  }
  
  # Return the plots and updated data
  return(list(
    stromal_plot = p1,
    immune_plot = p2,
    combined_plot = if(show_combined) p1 / p2 else NULL,
    data_with_categories = correlation_results
  ))
}

# Usage examples:
# For ERBB2
erbb2_plots <- plot_signature_correlations(
  correlation_results = erbb2_correlations,
  network_name = "ERBB2",
  show_combined = TRUE,
  show_summary = TRUE
)

egfr_plots <- plot_signature_correlations(
  correlation_results = egfr_correlations,
  network_name = "ERBB2",
  show_combined = TRUE,
  show_summary = TRUE
)

fgr3_plots <- plot_signature_correlations(
  correlation_results = fgfr3_correlations,
  network_name = "ERBB2",
  show_combined = TRUE,
  show_summary = TRUE
)

# For EGFR
# egfr_plots <- plot_signature_correlations(
#   correlation_results = egfr_correlations,
#   network_name = "EGFR"
# )

# For FGFR3
# fgfr3_plots <- plot_signature_correlations(
#   correlation_results = fgfr3_correlations,
#   network_name = "FGFR3"
# )
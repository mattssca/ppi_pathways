# Calculate gene means and create ranking
gene_means <- rowMeans(expr_data)

gene_rank_data <- data.frame(
  gene = names(gene_means),
  mean_expression = gene_means
) %>%
  arrange(desc(mean_expression)) %>%
  mutate(rank = row_number())

ggplot() +
  # Background layer: all non-network genes
  geom_point(data = gene_rank_data[gene_rank_data$in_network == FALSE, ],
             aes(x = rank, y = mean_expression),
             color = "grey", alpha = 1, size = 0.6) +

  # Foreground layer: network genes with better visibility
  geom_point(data = gene_rank_data[gene_rank_data$in_network == TRUE, ],
             aes(x = rank, y = mean_expression),
             color = "red", alpha = 1, size = 0.6) +

  # Add quantile lines
  geom_hline(yintercept = quantiles[1], color = "darkred", linetype = "dashed", alpha = 0.7) +
  geom_hline(yintercept = quantiles[2], color = "forestgreen", linetype = "dashed", alpha = 0.7) +
  geom_hline(yintercept = quantiles[3], color = "navyblue", linetype = "dashed", alpha = 0.7) +

  # Add annotations
  annotate("text", x = nrow(gene_rank_data) * 0.8, y = quantiles[1] + 0.02,
           label = paste("25th percentile:", round(quantiles[1], 3)),
           color = "black", size = 3) +
  annotate("text", x = nrow(gene_rank_data) * 0.8, y = quantiles[2] + 0.02,
           label = paste("Median:", round(quantiles[2], 3)),
           color = "black", size = 3) +
  annotate("text", x = nrow(gene_rank_data) * 0.8, y = quantiles[3] + 0.02,
           label = paste("75th percentile:", round(quantiles[3], 3)),
           color = "black", size = 3) +

  # Add manual legend
  annotate("point", x = nrow(gene_rank_data) * 0.05, y = max(gene_rank_data$mean_expression) * 0.9,
           color = "red", size = 1) +
  annotate("text", x = nrow(gene_rank_data) * 0.08, y = max(gene_rank_data$mean_expression) * 0.9,
           label = "Network genes", hjust = 0, size = 3) +
  annotate("point", x = nrow(gene_rank_data) * 0.05, y = max(gene_rank_data$mean_expression) * 0.85,
           color = "lightgray", size = 1) +
  annotate("text", x = nrow(gene_rank_data) * 0.08, y = max(gene_rank_data$mean_expression) * 0.85,
           label = "Other genes", hjust = 0, size = 3) +

  labs(
    title = "Gene Expression Ranking with Distribution Statistics",
    x = "Gene Rank (1 = Highest Expression)",
    y = "Mean Expression"
  ) +
  theme_minimal()

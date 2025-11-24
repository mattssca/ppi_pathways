#load data
load("ppi_signal/data/uroscanseq_data_Rdata")
load("ppi_signal/results_data/pathway_genes_df.Rdata")

#subset data object
expr_data = uroscanseq_data$expr_df

#examine expression distribution to choose appropriate thresholds
gene_means <- rowMeans(expr_data)

#summary statistics
summary(gene_means)

#visualize
hist(gene_means, breaks = 50, main = "Distribution of Gene Mean Expression")
abline(v = quantile(gene_means, c(0.25, 0.5, 0.75)), col = "red", lty = 2)

#check some specific pathway genes
pathway_gene_means <- gene_means[names(gene_means) %in% pathway_genes$gene]
summary(pathway_gene_means)

#remove bottom 25% of genes by expression
expression_cutoff <- quantile(gene_means, 0.25)  # -0.081
expressed_genes <- names(gene_means)[gene_means > expression_cutoff]

#filter pathway genes
filtered_pathway_genes <- pathway_genes[pathway_genes$gene %in% expressed_genes, ]

#filter expression data to keep only well-expressed genes
filtered_expr_data <- expr_data[rownames(expr_data) %in% expressed_genes, ]

cat("Expression cutoff (25th percentile):", round(expression_cutoff, 3), "\n")
cat("Genes kept:", length(expressed_genes), "out of", length(gene_means), "\n")

#visualize the cutoff
cutoff <- quantile(gene_means, 0.25)
pathway_gene_means <- gene_means[names(gene_means) %in% pathway_genes$gene]

ggplot(data.frame(x = gene_means), aes(x = x)) +
  geom_histogram(bins = 50, fill = "lightgray", color = "black", alpha = 0.7) +
  geom_vline(xintercept = cutoff, color = "blue", linewidth = 1.2) +
  geom_rug(data = data.frame(x = pathway_gene_means),
           aes(x = x), color = "red", linewidth = 1) +
  annotate("text", x = cutoff + 0.01, y = Inf,
           label = paste("Cutoff:", round(cutoff, 3)),
           color = "blue", vjust = 1.2, hjust = 0) +
  labs(
    title = "Gene Expression Distribution with Filter",
    x = "Mean Gene Expression",
    y = "Count"
  ) +
  theme_minimal()

save(filtered_pathway_genes, file = "ppi_signal/results_data/filtered_pathway_genes.Rdata")
save(filtered_expr_data, file = "ppi_signal/results_data/filtered_expr_data.Rdata")

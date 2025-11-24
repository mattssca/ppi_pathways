head(egfr_metrics_consoldiated$name)
head(fgfr3_metrics_consoldiated$name)
head(erbb2_metrics_consoldiated$name)

# =============================================================================
# PATHWAY GENE OVERLAP ANALYSIS
# =============================================================================

library(ggplot2)
library(dplyr)
library(VennDiagram)
library(UpSetR)
library(gridExtra)
library(RColorBrewer)

# =============================================================================
# 1. EXTRACT GENE SETS FROM EACH PATHWAY
# =============================================================================

# Extract gene lists
egfr_genes <- egfr_metrics_consoldiated$name
fgfr3_genes <- fgfr3_metrics_consoldiated$name
erbb2_genes <- erbb2_metrics_consoldiated$name

# Create named list for analysis
gene_sets <- list(
  EGFR = egfr_genes,
  FGFR3 = fgfr3_genes,
  ERBB2 = erbb2_genes
)

# Basic statistics
cat("=== PATHWAY GENE SET SIZES ===\n")
cat("EGFR network genes:", length(egfr_genes), "\n")
cat("FGFR3 network genes:", length(fgfr3_genes), "\n")
cat("ERBB2 network genes:", length(erbb2_genes), "\n")
cat("Total unique genes:", length(unique(c(egfr_genes, fgfr3_genes, erbb2_genes))), "\n\n")

# =============================================================================
# 2. PAIRWISE OVERLAPS
# =============================================================================

# Calculate pairwise intersections
egfr_fgfr3_overlap <- intersect(egfr_genes, fgfr3_genes)
egfr_erbb2_overlap <- intersect(egfr_genes, erbb2_genes)
fgfr3_erbb2_overlap <- intersect(fgfr3_genes, erbb2_genes)

# Three-way overlap
all_three_overlap <- intersect(intersect(egfr_genes, fgfr3_genes), erbb2_genes)

cat("=== PAIRWISE OVERLAPS ===\n")
cat("EGFR ∩ FGFR3:", length(egfr_fgfr3_overlap), "genes\n")
cat("EGFR ∩ ERBB2:", length(egfr_erbb2_overlap), "genes\n")
cat("FGFR3 ∩ ERBB2:", length(fgfr3_erbb2_overlap), "genes\n")
cat("All three pathways:", length(all_three_overlap), "genes\n\n")

# =============================================================================
# 3. CREATE OVERLAP MATRIX FOR ANALYSIS
# =============================================================================

# Create a comprehensive overlap matrix
all_genes <- unique(c(egfr_genes, fgfr3_genes, erbb2_genes))

overlap_matrix <- data.frame(
  gene = all_genes,
  in_EGFR = all_genes %in% egfr_genes,
  in_FGFR3 = all_genes %in% fgfr3_genes,
  in_ERBB2 = all_genes %in% erbb2_genes
) %>%
  mutate(
    num_pathways = as.numeric(in_EGFR) + as.numeric(in_FGFR3) + as.numeric(in_ERBB2),
    pathway_combination = case_when(
      in_EGFR & in_FGFR3 & in_ERBB2 ~ "All_Three",
      in_EGFR & in_FGFR3 & !in_ERBB2 ~ "EGFR_FGFR3",
      in_EGFR & !in_FGFR3 & in_ERBB2 ~ "EGFR_ERBB2",
      !in_EGFR & in_FGFR3 & in_ERBB2 ~ "FGFR3_ERBB2",
      in_EGFR & !in_FGFR3 & !in_ERBB2 ~ "EGFR_Only",
      !in_EGFR & in_FGFR3 & !in_ERBB2 ~ "FGFR3_Only",
      !in_EGFR & !in_FGFR3 & in_ERBB2 ~ "ERBB2_Only"
    )
  )

# Summary of overlap patterns
overlap_summary <- overlap_matrix %>%
  count(pathway_combination, num_pathways) %>%
  arrange(desc(n))

cat("=== OVERLAP PATTERNS ===\n")
print(overlap_summary)

# =============================================================================
# 4. VISUALIZATIONS
# =============================================================================

# A. Venn Diagram
venn_plot <- venn.diagram(
  x = gene_sets,
  category.names = c("EGFR", "FGFR3", "ERBB2"),
  filename = NULL,
  output = TRUE,

  # Appearance
  lwd = 2,
  lty = 'blank',
  fill = c("#FF6B6B", "#4ECDC4", "#45B7D1"),
  alpha = 0.6,

  # Numbers
  cex = 1.2,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),

  # Margins
  margin = 0.2
)

# Display Venn diagram
grid.newpage()
grid.draw(venn_plot)

# B. UpSet Plot for better visualization of complex overlaps
upset_data <- overlap_matrix %>%
  dplyr::select(in_EGFR, in_FGFR3, in_ERBB2) %>%
  as.data.frame()

# Convert logical to numeric for UpSetR
upset_data[] <- lapply(upset_data, as.numeric)
names(upset_data) <- c("EGFR", "FGFR3", "ERBB2")

upset_plot <- upset(upset_data,
                    sets = c("EGFR", "FGFR3", "ERBB2"),
                    sets.bar.color = c("#FF6B6B", "#4ECDC4", "#45B7D1"),
                    order.by = "freq",
                    empty.intersections = "on")

print(upset_plot)

# C. Bar plot of overlap patterns
p_overlap_bars <- overlap_summary %>%
  ggplot(aes(x = reorder(pathway_combination, n), y = n, fill = factor(num_pathways))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Gene Distribution Across Pathway Combinations",
       x = "Pathway Combination", y = "Number of Genes",
       fill = "Number of\nPathways") +
  theme_minimal() +
  scale_fill_brewer(type = "seq", palette = "Blues") +
  theme(legend.position = "bottom")

print(p_overlap_bars)

# =============================================================================
# 5. FUNCTIONAL ANALYSIS OF OVERLAPPING GENES
# =============================================================================

# Analyze functional characteristics of genes in different overlap categories
analyze_overlap_functions <- function(overlap_matrix, metrics_data) {

  # Combine all pathway metrics (you'll need to combine your three datasets)
  all_pathway_metrics <- bind_rows(
    egfr_metrics_consoldiated %>% mutate(source_pathway = "EGFR"),
    fgfr3_metrics_consoldiated %>% mutate(source_pathway = "FGFR3"),
    erbb2_metrics_consoldiated %>% mutate(source_pathway = "ERBB2")
  )

  # For genes in multiple pathways, take the maximum functional significance
  gene_func_summary <- all_pathway_metrics %>%
    group_by(name) %>%
    summarise(
      max_functional_significance = max(functional_significance, na.rm = TRUE),
      max_degree = max(degree, na.rm = TRUE),
      max_betweenness = max(betweenness, na.rm = TRUE),
      primary_functional_category = functional_category[which.max(functional_significance)],
      primary_pathway_role = pathway_role[which.max(functional_significance)],
      pathways_present = paste(unique(source_pathway), collapse = ","),
      .groups = 'drop'
    )

  # Merge with overlap information
  overlap_functional <- overlap_matrix %>%
    left_join(gene_func_summary, by = c("gene" = "name"))

  return(overlap_functional)
}

overlap_functional <- analyze_overlap_functions(overlap_matrix, egfr_metrics_consoldiated)

# Functional category distribution by overlap pattern
func_by_overlap <- overlap_functional %>%
  filter(!is.na(primary_functional_category)) %>%
  count(pathway_combination, primary_functional_category) %>%
  group_by(pathway_combination) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()

p_func_overlap <- func_by_overlap %>%
  ggplot(aes(x = pathway_combination, y = percentage, fill = primary_functional_category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Functional Composition of Gene Overlap Categories",
       x = "Pathway Combination", y = "Percentage",
       fill = "Functional Category") +
  theme_minimal() +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "bottom")

print(p_func_overlap)

# =============================================================================
# COMPREHENSIVE PPI NETWORK FUNCTIONAL ANALYSIS
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(corrplot)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(gridExtra)

# Load consolidated data (assuming you have this from all three pathways)
# load("ppi_signal/results_data/egfr_metrics_consolidated.Rdata")
# load("ppi_signal/results_data/fgfr3_metrics_consolidated.Rdata")
# load("ppi_signal/results_data/erbb2_metrics_consolidated.Rdata")

# Combine all pathway data for comparative analysis
# all_metrics <- rbind(egfr_metrics_consolidated, fgfr3_metrics_consolidated, erbb2_metrics_consolidated)

# For now, using egfr_metrics_consolidated as example
all_metrics <- egfr_metrics_consoldiated  # Note: fix typo in your variable name

# =============================================================================
# 1. FUNCTIONAL LANDSCAPE OVERVIEW
# =============================================================================

# A. Functional category distribution
p1_func_dist <- all_metrics %>%
  count(pathway, functional_category) %>%
  ggplot(aes(x = functional_category, y = n, fill = pathway)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Functional Category Distribution Across Pathways",
       x = "Functional Category", y = "Number of Genes",
       fill = "Pathway") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_d()

# B. Pathway role distribution
p1_role_dist <- all_metrics %>%
  count(pathway, pathway_role) %>%
  ggplot(aes(x = pathway_role, y = n, fill = pathway)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Pathway Role Distribution",
       x = "Pathway Role", y = "Number of Genes",
       fill = "Pathway") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_d()

print(p1_func_dist / p1_role_dist)

# =============================================================================
# 2. NETWORK CENTRALITY vs FUNCTIONAL IMPORTANCE
# =============================================================================

# Key insight: Are functionally important genes also topologically central?
p2_central_func <- all_metrics %>%
  ggplot(aes(x = degree, y = functional_significance,
             color = functional_category, size = betweenness)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(1, 8)) +
  labs(title = "Network Centrality vs Functional Importance",
       subtitle = "Size = Betweenness Centrality",
       x = "Degree (Number of Connections)",
       y = "Functional Significance (-log10 p-value)",
       color = "Functional Category") +
  theme_minimal() +
  scale_color_brewer(type = "qual", palette = "Set2") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed")

# Add correlation coefficient
cor_coef <- cor(all_metrics$degree, all_metrics$functional_significance, use = "complete.obs")
p2_central_func <- p2_central_func +
  annotate("text", x = Inf, y = Inf,
           label = paste("r =", round(cor_coef, 3)),
           hjust = 1.1, vjust = 1.5, size = 4)

print(p2_central_func)

# =============================================================================
# 3. EXPRESSION VARIABILITY vs NETWORK PROPERTIES
# =============================================================================

# Do network hubs have more variable expression across cancer subtypes?
p3_expr_var <- all_metrics %>%
  ggplot(aes(x = pathway_role, y = expression_range, fill = functional_category)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Expression Variability by Network Role",
       subtitle = "Across Cancer Subtypes",
       x = "Pathway Role", y = "Expression Range",
       fill = "Functional Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white")

print(p3_expr_var)

# Statistical test
library(broom)
expr_anova <- aov(expression_range ~ pathway_role, data = all_metrics)
cat("ANOVA: Expression Range ~ Pathway Role\n")
print(tidy(expr_anova))

# =============================================================================
# 4. FUNCTIONAL-STRUCTURAL IMPORTANCE ANALYSIS
# =============================================================================

# Identify top candidates for therapeutic targeting
top_candidates <- all_metrics %>%
  arrange(desc(functional_structural_importance)) %>%
  head(20) %>%
  mutate(gene_label = paste0(name, " (", functional_category, ")"))

p4_candidates <- top_candidates %>%
  ggplot(aes(x = reorder(gene_label, functional_structural_importance),
             y = functional_structural_importance,
             fill = pathway_role)) +
  geom_col() +
  coord_flip() +
  labs(title = "Top 20 Therapeutic Candidates",
       subtitle = "Based on Functional-Structural Importance",
       x = "Gene (Functional Category)",
       y = "Functional-Structural Importance Score",
       fill = "Pathway Role") +
  theme_minimal() +
  scale_fill_brewer(type = "qual", palette = "Set1")

print(p4_candidates)

# =============================================================================
# 5. SEED GENE ANALYSIS
# =============================================================================

# Are seed genes enriched in specific functional categories?
seed_analysis <- all_metrics %>%
  group_by(functional_category, is_seed) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(functional_category) %>%
  mutate(total = sum(count),
         percentage = count / total * 100) %>%
  filter(is_seed == TRUE)

p5_seed <- seed_analysis %>%
  ggplot(aes(x = reorder(functional_category, percentage), y = percentage)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  coord_flip() +
  labs(title = "Seed Gene Enrichment by Functional Category",
       x = "Functional Category",
       y = "Percentage of Genes that are Seeds") +
  theme_minimal() +
  geom_hline(yintercept = mean(all_metrics$is_seed) * 100,
             linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 1, y = mean(all_metrics$is_seed) * 100 + 2,
           label = "Overall Seed %", color = "red", size = 3)

print(p5_seed)

# Chi-square test for seed enrichment
seed_table <- table(all_metrics$functional_category, all_metrics$is_seed)
cat("Chi-square test: Functional Category ~ Seed Status\n")
print(chisq.test(seed_table))

# =============================================================================
# 6. COMMUNITY vs FUNCTIONAL ORGANIZATION
# =============================================================================

# Do network communities align with functional modules?
community_func <- all_metrics %>%
  count(community, functional_category) %>%
  group_by(community) %>%
  mutate(community_total = sum(n),
         percentage = n / community_total * 100)

p6_community <- community_func %>%
  ggplot(aes(x = factor(community), y = percentage, fill = functional_category)) +
  geom_bar(stat = "identity") +
  labs(title = "Functional Composition of Network Communities",
       x = "Network Community", y = "Percentage",
       fill = "Functional Category") +
  theme_minimal() +
  scale_fill_brewer(type = "qual", palette = "Set3")

print(p6_community)

# =============================================================================
# 7. EXPRESSION HEATMAP OF FUNCTIONAL HUBS
# =============================================================================

# Focus on top functional hubs and their expression patterns
functional_hubs <- all_metrics %>%
  filter(pathway_role == "Functional_Hub" | functional_structural_importance > 2) %>%
  arrange(desc(functional_structural_importance)) %>%
  head(20)

# Prepare expression matrix
expr_matrix <- functional_hubs %>%
  dplyr::select(name, mean_expr_UroB, mean_expr_GU, mean_expr_UroA,
         mean_expr_BaSq, mean_expr_UroC) %>%
  column_to_rownames("name") %>%
  as.matrix()

# Annotation for heatmap
annotation_row <- functional_hubs %>%
  dplyr::select(name, functional_category, pathway_role) %>%
  column_to_rownames("name")

# Create heatmap
pheatmap(expr_matrix,
         annotation_row = annotation_row,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Expression Patterns of Top Functional Hubs\nAcross Cancer Subtypes",
         fontsize = 10,
         angle_col = 45)

# =============================================================================
# 8. CORRELATION ANALYSIS OF KEY METRICS
# =============================================================================

# Select key numeric variables for correlation analysis
cor_vars <- all_metrics %>%
  dplyr::select(degree, betweenness, eigenvector, functional_significance,
         functional_diversity, expression_range, overall_sd,
         functional_structural_importance, global_mean) %>%
  na.omit()

cor_matrix <- cor(cor_vars)

# Create correlation plot
corrplot(cor_matrix, method = "color", type = "upper",
         order = "hclust", tl.col = "black", tl.srt = 45,
         title = "Correlation Matrix: Network & Functional Properties",
         mar = c(0,0,2,0))

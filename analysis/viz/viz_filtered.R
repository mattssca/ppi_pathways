library(UpSetR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggalluvial)
library(circlize)

# Define consistent color palette
pathway_colors <- c("EGFR" = "#E69F00",
                    "FGFR3" = "#56B4E9",
                    "ERBB2" = "#009E73")

# Load data
load("results_out/combined_metrics.Rdata")

# Get gene lists per pathway
egfr_genes <- combined_metrics %>% filter(pathway == "EGFR") %>% pull(name)
fgfr3_genes <- combined_metrics %>% filter(pathway == "FGFR3") %>% pull(name)
erbb2_genes <- combined_metrics %>% filter(pathway == "ERBB2") %>% pull(name)

# Create binary matrix
all_genes <- unique(combined_metrics$name)
upset_data <- data.frame(
  gene = all_genes,
  EGFR = as.integer(all_genes %in% egfr_genes),
  FGFR3 = as.integer(all_genes %in% fgfr3_genes),
  ERBB2 = as.integer(all_genes %in% erbb2_genes)
)

# Upset plot
upset(upset_data,
      sets = c("EGFR", "FGFR3", "ERBB2"),
      order.by = "freq",
      sets.bar.color = c(pathway_colors["EGFR"],
                         pathway_colors["FGFR3"],
                         pathway_colors["ERBB2"]),
      main.bar.color = "gray30",
      sets.x.label = "Pathway Gene Count",
      mainbar.y.label = "Gene Intersection Size",
      text.scale = c(1.5, 1.3, 1.3, 1.2, 1.5, 1.3))

# Add title using grid
library(grid)
grid.text("Gene Overlap Across RTK Signaling Pathways",
          x = 0.65, y = 0.95,
          gp = gpar(fontsize = 16, fontface = "bold"))

# Bubble/Scatter Plot
ggplot(combined_metrics, aes(x = degree, y = global_mean)) +
  geom_point(aes(size = betweenness,
                 color = pathway,
                 shape = hpa_expression_pattern),
             alpha = 0.7) +
  geom_text_repel(data = combined_metrics %>% filter(degree > 100),
                  aes(label = name)) +
  scale_color_manual(values = pathway_colors) +
  scale_size_continuous(range = c(2, 10)) +
  labs(title = "Network Hubs vs Expression Level",
       subtitle = "Node size represents betweenness centrality",
       x = "Degree (Network Connectivity)",
       y = "Mean Expression (All Subtypes)",
       color = "Pathway",
       shape = "HPA Expression",
       size = "Betweenness") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12))


# Alluvial plot
flow_data <- combined_metrics %>%
  count(pathway, hpa_expression_pattern, hpa_cancer_confidence)

ggplot(flow_data,
       aes(y = n, axis1 = pathway, axis2 = hpa_expression_pattern,
           axis3 = hpa_cancer_confidence)) +
  geom_alluvium(aes(fill = pathway), width = 1/12) +
  geom_stratum(width = 1/12, fill = "gray80", color = "white") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_manual(values = pathway_colors) +
  scale_x_discrete(limits = c("Pathway", "Expression Pattern", "Cancer Confidence")) +
  labs(title = "Gene Classification Flow Across Analysis Layers",
       subtitle = "Pathway → Protein Expression Pattern → Cancer Confidence",
       y = "Number of Genes",
       fill = "Pathway") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.text.x = element_text(size = 11, face = "bold"))


# Chord plot
# Create adjacency matrix
pathways <- c("EGFR", "FGFR3", "ERBB2")
overlap_matrix <- matrix(0, nrow = 3, ncol = 3,
                         dimnames = list(pathways, pathways))

for(i in 1:3) {
  for(j in 1:3) {
    genes_i <- combined_metrics %>% filter(pathway == pathways[i]) %>% pull(name)
    genes_j <- combined_metrics %>% filter(pathway == pathways[j]) %>% pull(name)
    overlap_matrix[i,j] <- length(intersect(genes_i, genes_j))
  }
}

# Set up plotting parameters
par(mar = c(1, 1, 3, 1))

chordDiagram(overlap_matrix,
             grid.col = pathway_colors,
             transparency = 0.5,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.3))

# Add title
title("Pathway Gene Overlap Network",
      cex.main = 1.5,
      font.main = 2)

# Add labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1],
              CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5), cex = 1.2)
}, bg.border = NA)

# Clear circlize plot
circos.clear()


# Waterfall plot
# Rank genes by expression range across subtypes
waterfall_data <- combined_metrics %>%
  arrange(desc(expression_range)) %>%
  mutate(rank = row_number())

ggplot(waterfall_data, aes(x = rank, y = expression_range)) +
  geom_segment(aes(xend = rank, yend = 0, color = pathway), size = 1) +
  geom_point(aes(color = pathway, shape = hpa_expression_pattern), size = 2) +
  scale_color_manual(values = pathway_colors) +
  labs(title = "Gene Expression Variability Across Bladder Cancer Subtypes",
       subtitle = "Genes ranked by expression range (max - min across subtypes)",
       x = "Gene Rank",
       y = "Expression Range",
       color = "Pathway",
       shape = "HPA Expression") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12))

# Calculate pathway overlap for each gene
pathway_overlap <- combined_metrics %>%
  group_by(name) %>%
  summarise(
    pathways = paste(sort(unique(pathway)), collapse = "+"),
    n_pathways = n_distinct(pathway)
  ) %>%
  ungroup()

# Add overlap info to main data
combined_metrics <- combined_metrics %>%
  left_join(pathway_overlap, by = "name")

# Check what pathway combinations actually exist
cat("\nActual pathway combinations in data:\n")
print(unique(combined_metrics$pathways))

# Define colors ONLY for combinations that exist
# First get the unique combinations
existing_combos <- unique(combined_metrics$pathways)

# Create color palette based on what exists
overlap_colors <- setNames(
  c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#0072B2", "#F0E442"),
  c("EGFR", "FGFR3", "ERBB2", "EGFR+ERBB2", "EGFR+FGFR3", "ERBB2+FGFR3", "EGFR+ERBB2+FGFR3")
)

# Keep only colors for combinations that exist in data
overlap_colors <- overlap_colors[names(overlap_colors) %in% existing_combos]

# Create label data AFTER adding n_pathways
label_data <- combined_metrics %>%
  filter(degree > 100 | n_pathways > 1)

# Plot with overlap-based colors
ggplot(combined_metrics, aes(x = degree, y = global_mean)) +
  geom_point(aes(size = betweenness,
                 color = pathways,
                 shape = hpa_expression_pattern),
             alpha = 0.8) +
  geom_text_repel(data = label_data,
                  aes(label = name, color = pathways),
                  size = 3,
                  fontface = ifelse(label_data$n_pathways == 3, "bold", "plain"),
                  show.legend = FALSE) +
  scale_color_manual(values = overlap_colors,
                     name = "Pathway Membership") +
  scale_size_continuous(range = c(2, 10)) +
  labs(title = "Network Hubs vs Expression Level",
       subtitle = "Color indicates pathway membership (single, paired, or core)",
       x = "Degree (Network Connectivity)",
       y = "Mean Expression (All Subtypes)",
       shape = "HPA Expression",
       size = "Betweenness") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.position = "right")

# Print summary of overlap scenarios
cat("\nPathway Overlap Summary:\n")
overlap_summary <- combined_metrics %>%
  distinct(name, pathways, n_pathways) %>%
  count(pathways, n_pathways) %>%
  arrange(n_pathways, pathways)

print(overlap_summary)

# List core genes if they exist
core_genes <- combined_metrics %>%
  filter(n_pathways == 3) %>%
  distinct(name) %>%
  pull(name)

if(length(core_genes) > 0) {
  cat("\nCore genes (all 3 pathways):\n")
  cat(paste(core_genes, collapse = ", "), "\n")
} else {
  cat("\nNo genes found in all 3 pathways.\n")
}

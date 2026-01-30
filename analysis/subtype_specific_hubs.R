library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(igraph)

load('ppi_pathways/data/uroscanseq_data_.Rdata')
load('ppi_pathways/results_out/erbb2_meta_consolidated.Rdata')
load('ppi_pathways/results_out/fgfr3_meta_consolidated.Rdata')
load("~/BIOINFORMATICS/basq_invest/phase_2/data/predicted_imvigor.Rdata")
load("~/BIOINFORMATICS/basq_invest/phase_1/data_out/meta_imvigor.Rdata")
source("ppi_pathways/functions/plot_rank_heatmap.R")
source("ppi_pathways/functions/calculate_hub_scores.R")
source("../../BIOINFORMATICS/GIT_HUB/ppi_pathways/functions/surv_plot.R")
source("../../BIOINFORMATICS/GIT_HUB/ppi_pathways/functions/create_enhanced_orthogonal_discovery_proof.R")
source("../../BIOINFORMATICS/GIT_HUB/ppi_pathways/functions/create_zscore_clustered_heatmap.R")

#prepare data
ppi_data <- erbb2_meta_consoldiated
node_data <- read.csv(file = "ppi_pathways/results_data/erbb_node_data.csv")
top_genes <- node_data %>% dplyr::select(shared.name)
ppi_data <- ppi_data %>% filter(name %in% top_genes$shared.name)

#### re-calculating network centralities
#get the filtered gene set
filtered_genes <- ppi_data$name
cat("Number of genes in filtered network:", length(filtered_genes), "\n")

#filter the edge list to only include connections between these genes
filtered_edges <- erbb2_edges %>%
  filter(from %in% filtered_genes & to %in% filtered_genes) %>%
  distinct(from, to, .keep_all = TRUE)

#create the filtered igraph network
filtered_network <- graph_from_data_frame(filtered_edges, directed = FALSE)

#remove self-loops if any
filtered_network <- simplify(filtered_network, remove.multiple = TRUE, remove.loops = TRUE)

#calculate new centralities on the filtered network
hits_scores <- hits_scores(filtered_network)

new_centralities <- data.frame(
  name = V(filtered_network)$name,
  filtered_degree = degree(filtered_network),
  filtered_betweenness = betweenness(filtered_network, normalized = TRUE),
  filtered_closeness = closeness(filtered_network, normalized = TRUE),
  filtered_eigenvector = eigen_centrality(filtered_network)$vector,
  filtered_hub_score = hits_scores$hub,
  filtered_authority = hits_scores$authority
)

#replace the old centralities
#select only the new centrality columns and rename them to match the old names
new_centralities_renamed <- new_centralities %>%
  dplyr::rename(
    degree = filtered_degree,
    betweenness = filtered_betweenness,
    closeness = filtered_closeness,
    eigenvector = filtered_eigenvector,
    hub_score = filtered_hub_score,
    authority_score = filtered_authority
  )

#remove the old centrality columns from ppi_data
ppi_data_updated <- ppi_data %>%
  dplyr::select(-degree, -betweenness, -closeness, -eigenvector, -hub_score, -authority_score)

#join the new centrality columns back in by 'name'
ppi_data <- ppi_data_updated %>%
  dplyr::left_join(new_centralities_renamed, by = "name")

#check that the update worked
dplyr::glimpse(ppi_data)

#add 5 class info also
ppi_data <- ppi_data %>%
  rename(
    rank_7_UroA = rank_UroA,
    rank_7_UroB = rank_UroB,
    rank_7_UroC = rank_UroC,
    rank_7_GU = rank_GU,
    rank_7_BaSq = rank_BaSq
  )

uro_samples <- names(uroscanseq_data$subtype_5_vector)[uroscanseq_data$subtype_5_vector == "Uro"]
gu_samples <- names(uroscanseq_data$subtype_5_vector)[uroscanseq_data$subtype_5_vector == "GU"]
basq_samples <- names(uroscanseq_data$subtype_5_vector)[uroscanseq_data$subtype_5_vector == "BaSq"]

genes_in_ppi <- ppi_data$name

mean_expr_uro <- rowMeans(uroscanseq_data$expr_df[genes_in_ppi, uro_samples], na.rm = TRUE)
mean_expr_gu_5 <- rowMeans(uroscanseq_data$expr_df[genes_in_ppi, gu_samples], na.rm = TRUE)
mean_expr_basq_5 <- rowMeans(uroscanseq_data$expr_df[genes_in_ppi, basq_samples], na.rm = TRUE)

ppi_data <- ppi_data %>%
  mutate(
    mean_expr_5_Uro = mean_expr_uro[name],
    mean_expr_5_GU = mean_expr_gu_5[name],
    mean_expr_5_BaSq = mean_expr_basq_5[name]
  )

ppi_data <- ppi_data %>%
  rowwise() %>%
  mutate(
    expr_5class = list(c(Uro = mean_expr_5_Uro, 
                         GU = mean_expr_5_GU, 
                         BaSq = mean_expr_5_BaSq)),
    ranks_5class = list(rank(-unlist(expr_5class), ties.method = "first")),
    rank_5_Uro = ranks_5class[[1]],
    rank_5_GU = ranks_5class[[2]],
    rank_5_BaSq = ranks_5class[[3]]
  ) %>%
  ungroup() %>%
  select(-expr_5class, -ranks_5class)

#plots!
####################################################################################################
plot_rank_heatmap(data = ppi_data, system = "5", order_by = "hub_score", top_n = 10)
plot_rank_heatmap(data = ppi_data, system = "7", order_by = "hub_score")

calculate_hub_scores(data = ppi_data, system = "5", order_by = "hub_score", top_n = 10)
calculate_hub_scores(data = ppi_data, system = "7", order_by = "hub_score")

####################################################################################################
#create survival plot with optimal cutoff (80th percentile)
create_gene_survival_plot(
  expression_data = sjodahl_pred$data,
  predictions = sjodahl_pred$predictions_5classes,
  meta_data = sjodahl_2017_meta,
  gene_name = "CTNNB1",
  subtypes_include = c("Uro", "GU"),
  expression_cutoff = 0.8,
  time_column = "surv_os_time",
  event_column = "surv_os_event",
  title = "CTNNB1 Expression in Uro & GU Subtypes"
)

create_gene_survival_plot(
  expression_data = sjodahl_pred$data,
  predictions = sjodahl_pred$predictions_5classes,
  meta_data = sjodahl_2017_meta,
  gene_name = "CTNNB1",
  subtypes_include = c("Uro"),
  expression_cutoff = 0.8,
  time_column = "surv_os_time",
  event_column = "surv_os_event",
  title = "CTNNB1 Expression Uro Subtype"
)

create_gene_survival_plot(
  expression_data = sjodahl_pred$data,
  predictions = sjodahl_pred$predictions_5classes,
  meta_data = sjodahl_2017_meta,
  gene_name = "CTNNB1",
  subtypes_include = c("GU"),
  expression_cutoff = 0.8,
  time_column = "surv_os_time",
  event_column = "surv_os_event",
  title = "CTNNB1 Expression GU Subtype"
)

####################################################################################################
# Create the enhanced figure
enhanced_proof_figure <- create_enhanced_orthogonal_discovery_proof()
print(enhanced_proof_figure)

####################################################################################################
# Create and draw z-score normalized heatmap
zscore_heatmap <- create_zscore_clustered_heatmap()

draw(zscore_heatmap,
     column_title = "Z-score Normalized Expression Reveals Gene Clustering Patterns", 
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     padding = unit(c(2, 2, 2, 12), "mm"))

####################################################################################################
# Example usage with custom genes
data = plot_rank_heatmap(data = ppi_data, system = "5", order_by = "hub_score", top_n = 10)
exact_gene_comparison <- create_hardcoded_gene_comparison(hub_scores_df = hub_scores_df, uro_genes = c("CTNNB1", "AKT1"), gu_genes = c("ERBB2", "ERBB3"), basq_genes = c("EGFR", "SHC1"))

# Draw the comparison
draw(exact_gene_comparison$heatmap,
     column_title = "Network Weighting Enhances Subtype Discrimination: Selected Hub Genes",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     padding = unit(c(2, 2, 2, 2), "mm"))

#load data
load(file = "results_out/node_meta_consolidated.Rdata")
load("data/uroscanseq_data_Rdata")

#source function
source("functions/expand_signature_network.R")

#load packages
library(dplyr)
library(igraph)
library(STRINGdb)
library(ggraph)
library(biomaRt)

#get unique genes in the network
unique_genes = unique(node_meta_consoldiated$name)

#get ppi network info
all_nodes = expand_signature_network(expr_data = uroscanseq_data$expr_df,
                                     seed_genes = unique_genes,
                                     return_data = TRUE,
                                     expr_summary = "mean",
                                     max_added_genes = 0,
                                     layout_method = "kk",
                                     min_degree = 1,
                                     string_score_threshold = 500,
                                     subtype_vector = uroscanseq_data$subtype_7_vector)

library(dplyr)

#get vectors of gene names for each pathway
egfr_genes   <- egfr_meta_consoldiated$name
fgfr3_genes  <- fgfr3_meta_consoldiated$name
erbb2_genes  <- erbb2_meta_consoldiated$name

#for each gene, check which pathways it appears in
all_nodes$node_metrics <- all_nodes$node_metrics %>%
  rowwise() %>%
  mutate(
    pathway_cond = paste(
      c(
        if (name %in% egfr_genes) "EGFR" else NULL,
        if (name %in% fgfr3_genes) "FGFR3" else NULL,
        if (name %in% erbb2_genes) "ERBB2" else NULL
      ),
      collapse = ","
    ),
    pathway_cond = ifelse(
      grepl(",", pathway_cond),
      pathway_cond,
      paste0("Unique:", pathway_cond)
    )
  ) %>%
  ungroup()

#susbet
collapsed_node_metrics = all_nodes$node_metrics

#get relevant columns
these_columns = c("degree", "is_seed", "betweenness", "closeness", "eigenvector", "hub_score",
                  "authority_score", "community", "pathway", "IHC tissue name", "is_amplified",
                  "is_homdel", "tcga_mutsig_q_value", "tcga_num_mutations", "tcga_num_samples_mutated",
                  "tcga_profiled_samples", "tcga_frequency", "tcga_mutation_frequency", "pathway_sharing",
                  "n_pathways", "pathways_list")

#remove redundant columns from node data
tmp_node_metrics = node_meta_consoldiated %>%
  select(-all_of(these_columns))

tmp_node_metrics_unique <- tmp_node_metrics %>%
  distinct(name, .keep_all = TRUE)

merged_metrics <- collapsed_node_metrics %>%
  left_join(tmp_node_metrics_unique, by = "name")

#export to cytoscape
write.table(
  merged_metrics,
  file = "cytoscape/collapsed_all_node_metrics.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

all_nodes_igraph <- igraph::as_data_frame(all_nodes$g, what = "edges")

#write edge data
write.table(
  all_nodes_igraph,
  file = "cytoscape/collapsed_all_edges.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

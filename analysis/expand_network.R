#load data
load("ppi_signal/data/uroscanseq_data_Rdata")
load("ppi_signal/results_data/pathway_results.Rdata")

#source scripts
source("ppi_signal/functions/expand_signature_network.R")

#load packages
library(dplyr)
library(igraph)
library(STRINGdb)
library(ggraph)
library(biomaRt)

#create nodes object
egfr_nodes = expand_signature_network(expr_data = uroscanseq_data$expr_df,
                                      seed_genes = pathway_results$egfr_genes,
                                      return_data = TRUE,
                                      expr_summary = "mean",
                                      max_added_genes = 100,
                                      layout_method = "kk",
                                      min_degree = 1,
                                      string_score_threshold = 500,
                                      subtype_vector = uroscanseq_data$subtype_7_vector)

fgfr3_nodes = expand_signature_network(expr_data = uroscanseq_data$expr_df,
                                       seed_genes = pathway_results$fgfr3_genes,
                                       return_data = TRUE,
                                       expr_summary = "mean",
                                       max_added_genes = 100,
                                       layout_method = "kk",
                                       min_degree = 1,
                                       string_score_threshold = 500,
                                       subtype_vector = uroscanseq_data$subtype_7_vector)

erbb2_nodes = expand_signature_network(expr_data = uroscanseq_data$expr_df,
                                       seed_genes = pathway_results$erbb2_genes,
                                       return_data = TRUE,
                                       expr_summary = "mean",
                                       max_added_genes = 100,
                                       layout_method = "kk",
                                       min_degree = 1,
                                       string_score_threshold = 500,
                                       subtype_vector = uroscanseq_data$subtype_7_vector)

#get edge data
egfr_edge <- igraph::as_data_frame(egfr_nodes$g, what = "edges")
fgfr3_edge <- igraph::as_data_frame(fgfr3_nodes$g, what = "edges")
erbb2_edge <- igraph::as_data_frame(erbb2_nodes$g, what = "edges")

#extract node metrics
egfr_node_metrics = egfr_nodes$node_metrics
fgfr3_node_metrics = fgfr3_nodes$node_metrics
erbb2_node_metrics = erbb2_nodes$node_metrics

#export results
#write as tab-delimited file
write.table(egfr_edge, file = "../GIT_REPOS/ppi_signal/results_data/egfr_edge.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fgfr3_edge, file = "../GIT_REPOS/ppi_signal/results_data/fgfr3_edge.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(erbb2_edge, file = "../GIT_REPOS/ppi_signal/results_data/erbb2_edge.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#save nodes object
save(egfr_nodes, file = "../GIT_REPOS/ppi_signal/results_data/egfr_nodes_object.Rdata")
save(fgfr3_nodes, file = "../GIT_REPOS/ppi_signal/results_data/fgfr3_nodes_object.Rdata")
save(erbb2_nodes, file = "../GIT_REPOS/ppi_signal/results_data/erbb2_nodes_object.Rdata")

#save node metrics
save(egfr_node_metrics, file = "../GIT_REPOS/ppi_signal/results_data/egfr_node_metrics.Rdata")
save(fgfr3_node_metrics, file = "../GIT_REPOS/ppi_signal/results_data/fgfr3_node_metrics.Rdata")
save(erbb2_node_metrics, file = "../GIT_REPOS/ppi_signal/results_data/erbb2_node_metrics.Rdata")

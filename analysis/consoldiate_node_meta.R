#load data
load("data/uroscanseq_data_Rdata")
load("results_data/egfr_node_metrics.Rdata")
load("results_data/fgfr3_node_metrics.Rdata")
load("results_data/erbb2_node_metrics.Rdata")

#add cohort info to node metrics
egfr_node_metrics$pathway = "EGFR"
fgfr3_node_metrics$pathway = "FGFR3"
erbb2_node_metrics$pathway = "ERBB2"

#load packges
library(purrr)

####################################################################################################
#MERGE NODE METRICS
#merge the node metric data
node_metrics_combined = rbind(egfr_node_metrics, fgfr3_node_metrics, erbb2_node_metrics)

####################################################################################################
#ABOVE THRESHOLD INFORMATION
#load data
load("results_data/filtered_pathway_genes.Rdata")

#add info on above threshold expression
node_metrics_combined$is_well_expressed <- node_metrics_combined$name %in% filtered_pathway_genes$gene

#check the results
table(node_metrics_combined$is_well_expressed)

#save
save(node_metrics_combined, file = "results_data/node_metrics_combined.Rdata")

####################################################################################################
#GLOBAL EXPRESSION MEAN
#load data
load("data/uroscanseq_data_Rdata")

#subset expression data
expr_data = uroscanseq_data$expr_df

#calculate global mean expression for each gene
global_means <- rowMeans(expr_data)

#add the global_mean column to node_metrics_combined
node_metrics_combined$global_mean <- global_means[node_metrics_combined$name]

####################################################################################################
#ANOVA AND SIZE EFFECTS
#load data
load("results_data/anova_results.Rdata")

#add anova results to node metrics
node_metrics_combined <- node_metrics_combined %>%
  left_join(anova_with_effects, by = c("name" = "gene"))

####################################################################################################
#SUBTYPE RANK EXPRESSION
#rank the subtypes based on highest mean expression (highest gets rank 1)
rank_matrix <- t(apply(node_metrics_combined[, c("mean_expr_UroA", "mean_expr_UroB",
                                                 "mean_expr_UroC", "mean_expr_GU",
                                                 "mean_expr_BaSq")], 1, function(x) {
                                                   rank(-x, ties.method = "min")  #negative for descending rank (highest = 1)
                                                 }))

#convert to data frame and add proper column names
rank_df <- as.data.frame(rank_matrix)

colnames(rank_df) <- c("rank_UroA", "rank_UroB", "rank_UroC", "rank_GU", "rank_BaSq")

#add to node metrics
node_metrics_combined <- cbind(node_metrics_combined, rank_df)

####################################################################################################
#SINGLE CELL CLASSIFIACTION
#load data
load("results_data/egfr_single_cell_class.Rdata")
load("results_data/fgfr3_single_cell_class.Rdata")
load("results_data/erbb2_single_cell_class.Rdata")

#break out pathway information
egfr_metrics = node_metrics_combined %>% filter(pathway == "EGFR")
fgfr3_metrics = node_metrics_combined %>% filter(pathway == "FGFR3")
erbb2_metrics = node_metrics_combined %>% filter(pathway == "ERBB2")

#add single cell info
egfr_metrics = egfr_metrics %>% left_join(egfr_single_cell_class, by = "name")
fgfr3_metrics = fgfr3_metrics %>% left_join(fgfr3_single_cell_class, by = "name")
erbb2_metrics = erbb2_metrics %>% left_join(erbb2_single_cell_class, by = "name")

#combine back
node_metrics_combined = rbind(egfr_metrics, fgfr3_metrics, erbb2_metrics)

####################################################################################################
#FUNCTIONAL ANNOTATIONS
#load data
load("results_data/egfr_functional_annotations.Rdata")
load("results_data/fgfr3_functional_annotations.Rdata")
load("results_data/erbb2_functional_annotations.Rdata")

egfr_metrics = egfr_metrics %>% left_join(egfr_functional_annotations, by = "name")
fgfr3_metrics = fgfr3_metrics %>% left_join(fgfr3_functional_annotations, by = "name")
erbb2_metrics = erbb2_metrics %>% left_join(erbb2_functional_annotations, by = "name")

node_metrics_combined = rbind(egfr_metrics, fgfr3_metrics, erbb2_metrics)

####################################################################################################
#IHC META
#load data
load("results_data/egfr_hpa.Rdata")
load("results_data/fgfr3_hpa.Rdata")
load("results_data/erbb2_hpa.Rdata")

egfr_metrics = egfr_metrics %>% left_join(egfr_hpa, by = "name")
fgfr3_metrics = fgfr3_metrics %>% left_join(fgfr3_hpa, by = "name")
erbb2_metrics = erbb2_metrics %>% left_join(erbb2_hpa, by = "name")

####################################################################################################
#TCGA MUTATION DATA
load("C:/Users/matts/Desktop/GIT_REPOS/ppi_pathways/results_data/tcga_mut_data.Rdata")

#join with node metrics
egfr_metrics = egfr_metrics %>% left_join(tcga_mut_data, by = "name")
fgfr3_metrics = fgfr3_metrics %>% left_join(tcga_mut_data, by = "name")
erbb2_metrics = erbb2_metrics %>% left_join(tcga_mut_data, by = "name")

#combine return
node_metrics_combined = rbind(egfr_metrics, fgfr3_metrics, erbb2_metrics)

####################################################################################################
#EXPORTS
#rename for clarity
egfr_meta_consoldiated = egfr_metrics
fgfr3_meta_consoldiated = fgfr3_metrics
erbb2_meta_consoldiated = erbb2_metrics
node_meta_consoldiated = rbind(egfr_meta_consoldiated, fgfr3_meta_consoldiated, erbb2_meta_consoldiated)

#save
save(egfr_meta_consoldiated, file = "results_out/egfr_meta_consolidated.Rdata")
save(fgfr3_meta_consoldiated, file = "results_out/fgfr3_meta_consolidated.Rdata")
save(erbb2_meta_consoldiated, file = "results_out/erbb2_meta_consolidated.Rdata")
save(node_meta_consoldiated, file = "results_out/node_meta_consolidated.Rdata")

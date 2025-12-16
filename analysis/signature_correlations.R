#load packages
library(LundTaxR)
library(tibble)

#load data
load("ppi_pathways/data/uroscanseq_data_.Rdata")

load("ppi_pathways/results_out/erbb2_meta_consolidated.Rdata")
load("ppi_pathways/results_out/egfr_meta_consolidated.Rdata")
load("ppi_pathways/results_out/fgfr3_meta_consolidated.Rdata")

#subset expressiondata
expr_df = uroscanseq_data$expr_df

#get stromal and immuen signature genes
stromal_genes = as.data.frame(signatures$signatures_plot %>% filter(signature == "Stromal141_UP")) %>% pull(hgnc_symbol)
immune_genes = as.data.frame(signatures$signatures_plot %>% filter(signature == "Immune141_UP")) %>% pull(hgnc_symbol)

#subset expression to signature genes
stroma_expr = expr_df %>% 
  rownames_to_column("name") %>% 
  filter(name %in% stromal_genes) %>% 
  column_to_rownames("name")

immune_expr = expr_df %>% 
  rownames_to_column("name") %>% 
  filter(name %in% immune_genes) %>% 
  column_to_rownames("name")

head(stroma_expr)[1:5]
head(immune_expr)[1:5]

#calcualte scores
stroma_sample_means <- colMeans(stroma_expr, na.rm = TRUE)
immune_sample_means <- colMeans(immune_expr, na.rm = TRUE)

#run helper
#ERBB2
erbb2_correlations <- calculate_signature_correlations(node_meta = erbb2_meta_consoldiated,
                                                       expr_df = expr_df,
                                                       stroma_sample_means = stroma_sample_means,
                                                       immune_sample_means = immune_sample_means)

#ERBB2
egfr_correlations <- calculate_signature_correlations(node_meta = egfr_meta_consoldiated,
                                                       expr_df = expr_df,
                                                       stroma_sample_means = stroma_sample_means,
                                                       immune_sample_means = immune_sample_means)

#ERBB2
fgfr3_correlations <- calculate_signature_correlations(node_meta = fgfr3_meta_consoldiated,
                                                       expr_df = expr_df,
                                                       stroma_sample_means = stroma_sample_means,
                                                       immune_sample_means = immune_sample_means)

save(erbb2_correlations, file = "ppi_pathways/results_data/erbb2_correlations.Rdata")
save(egfr_correlations, file = "ppi_pathways/results_data/egfr_correlations.Rdata")
save(fgfr3_correlations, file = "ppi_pathways/results_data/fgfr3_correlations.Rdata")

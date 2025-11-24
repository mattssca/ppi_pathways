#subset the expression_data to pathway genes
#load expression data
load("../GIT_REPOS/ppi_signal/data/uroscanseq_data_Rdata")

#load pathway genes
load("../GIT_REPOS/ppi_signal/results_data/pathway_results.Rdata")

#unpack pathway genes
egfr_genes = data.frame(pathway_results$egfr_genes) %>% 
  rename(gene = `pathway_results.egfr_genes`) %>% 
  mutate(pathway = "EGFR")

fgfr3_genes = data.frame(pathway_results$fgfr3_genes) %>%
  rename(gene = `pathway_results.fgfr3_genes`) %>%
  mutate(pathway = "FGFR3")

erbb2_genes = data.frame(pathway_results$erbb2_genes) %>% 
  rename(gene = `pathway_results.erbb2_genes`) %>%
  mutate(pathway = "ERBB2")

pathway_genes = rbind(egfr_genes, fgfr3_genes, erbb2_genes)

#what genes are in uroscaseq expression data
uroscan_genes = row.names(uroscanseq_data$expr_df)

#remove signature genes not available in the expression data
not_in_exp = setdiff(pathway_genes$gene, uroscan_genes)

#remove rows with gene names in not_in_exp from the expression data
uroscanseq_data$expr_df <- uroscanseq_data$expr_df[!rownames(uroscanseq_data$expr_df) %in% not_in_exp, ]

#add genes to signature object
sig_genes$FGFR_signal <- genes_sig_pathway

save(uroscanseq_data, file = "../GIT_REPOS/ppi_signal/data/uroscanseq_data_Rdata")
save(pathway_genes, file = "../GIT_REPOS/ppi_signal/results_data/pathway_genes_df.Rdata")

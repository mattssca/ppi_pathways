calculate_signature_correlations <- function(node_meta, 
                                             expr_df, 
                                             stroma_sample_means, 
                                             immune_sample_means,
                                             meta_name_col = "name") {
  
  #extract genes from the node metadata
  node_genes <- node_meta[[meta_name_col]]
  
  #subset expression data to only these genes
  node_expr <- expr_df %>% 
    rownames_to_column("name") %>% 
    filter(name %in% node_genes) %>% 
    column_to_rownames("name")
  
  #initialize results dataframe
  correlation_results <- data.frame(
    gene = character(),
    stromal_cor = numeric(),
    stromal_pval = numeric(),
    immune_cor = numeric(),
    immune_pval = numeric(),
    stringsAsFactors = FALSE
  )
  
  #loop through each gene and calculate correlations
  for(gene in rownames(node_expr)) {
    #get gene expression across samples
    gene_expr <- as.numeric(node_expr[gene, ])
    
    #calculate correlation with stromal scores
    stromal_cor_test <- cor.test(gene_expr, stroma_sample_means)
    
    #calculate correlation with immune scores
    immune_cor_test <- cor.test(gene_expr, immune_sample_means)
    
    #add results to dataframe
    correlation_results <- rbind(correlation_results, 
                                 data.frame(
                                   gene = gene,
                                   stromal_cor = stromal_cor_test$estimate,
                                   stromal_pval = stromal_cor_test$p.value,
                                   immune_cor = immune_cor_test$estimate,
                                   immune_pval = immune_cor_test$p.value
                                 ))
  }
  
  #add stromal correlation summary
  correlation_results$stromal_summary <- case_when(
    correlation_results$stromal_pval >= 0.05 ~ "non_sig",
    correlation_results$stromal_cor > 0 & correlation_results$stromal_pval < 0.05 ~ "pos",
    correlation_results$stromal_cor < 0 & correlation_results$stromal_pval < 0.05 ~ "neg"
  )
  
  #add immune correlation summary
  correlation_results$immune_summary <- case_when(
    correlation_results$immune_pval >= 0.05 ~ "non_sig",
    correlation_results$immune_cor > 0 & correlation_results$immune_pval < 0.05 ~ "pos",
    correlation_results$immune_cor < 0 & correlation_results$immune_pval < 0.05 ~ "neg"
  )
  
  #convert to factors
  correlation_results$stromal_summary <- factor(correlation_results$stromal_summary, 
                                                levels = c("pos", "neg", "non_sig"))
  correlation_results$immune_summary <- factor(correlation_results$immune_summary, 
                                               levels = c("pos", "neg", "non_sig"))
  
  return(correlation_results)
}

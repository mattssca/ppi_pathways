#get all genes expressed in normal urothelium and passing the quality control
normal_urothelium = node_meta_consoldiated %>% filter(hpa_normal_expressed) %>%  filter(hpa_normal_reliable %in% c("High", "Medium"))
cancer_urothelium = node_meta_consoldiated %>% filter(hpa_cancer_expressed)

#get genes that are expressed in normal or cancer tissue, or both
#subset to genes passing the expression threshold
#subset to subtype specific genes
meta_sub = node_meta_consoldiated %>% 
  filter(hpa_passes_filter) %>% 
  filter(is_well_expressed) #%>% 
  filter(significant_05)

library(reactable)
reactable(meta_sub, 
          filterable = TRUE, 
          searchable = TRUE, 
          defaultPageSize = 20,
          sortable = TRUE,
          resizable = TRUE)

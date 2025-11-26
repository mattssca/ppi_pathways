uroa_mut = mutations_tcga_combined$UroA_mut_data %>% select(name, UroA_mut_freq)
urob_mut = mutations_tcga_combined$UroB_mut_data %>% select(name, UroB_mut_freq)
uroc_mut = mutations_tcga_combined$UroC_mut_data %>% select(name, UroC_mut_freq)
gu_mut = mutations_tcga_combined$GU_mut_data %>% select(name, GU_mut_freq)
basq_mut = mutations_tcga_combined$BaSq_mut_data %>% select(name, BaSq_mut_freq)

uroa_cna = cna_tcga_combined$UroA_mut_cna %>% select(name, UroA_alt_freq)
urob_cna = cna_tcga_combined$UroB_mut_cna %>% select(name, UroB_alt_freq)
uroc_cna = cna_tcga_combined$UroC_mut_cna %>% select(name, UroC_alt_freq)
gu_cna = cna_tcga_combined$GU_mut_cna %>% select(name, GU_alt_freq)
basq_cna = cna_tcga_combined$BaSq_mut_cna %>% select(name, BaSq_alt_freq)

erbb2_meta_consoldiated = erbb2_meta_consoldiated %>%
  left_join(uroa_mut, by = "name") %>%
  left_join(uroa_cna , by = "name") %>%
  left_join(urob_mut, by = "name") %>%
  left_join(urob_cna , by = "name") %>%
  left_join(uroc_mut, by = "name") %>%
  left_join(uroc_cna, by = "name") %>%
  left_join(gu_mut, by = "name") %>%
  left_join(gu_cna, by = "name") %>%
  left_join(basq_mut, by = "name") %>%
  left_join(basq_cna, by = "name")

egfr_meta_consoldiated = egfr_meta_consoldiated %>%
  left_join(uroa_mut, by = "name") %>%
  left_join(uroa_cna , by = "name") %>%
  left_join(urob_mut, by = "name") %>%
  left_join(urob_cna , by = "name") %>%
  left_join(uroc_mut, by = "name") %>%
  left_join(uroc_cna, by = "name") %>%
  left_join(gu_mut, by = "name") %>%
  left_join(gu_cna, by = "name") %>%
  left_join(basq_mut, by = "name") %>%
  left_join(basq_cna, by = "name")

fgfr3_meta_consoldiated = fgfr3_meta_consoldiated %>%
  left_join(uroa_mut, by = "name") %>%
  left_join(uroa_cna , by = "name") %>%
  left_join(urob_mut, by = "name") %>%
  left_join(urob_cna , by = "name") %>%
  left_join(uroc_mut, by = "name") %>%
  left_join(uroc_cna, by = "name") %>%
  left_join(gu_mut, by = "name") %>%
  left_join(gu_cna, by = "name") %>%
  left_join(basq_mut, by = "name") %>%
  left_join(basq_cna, by = "name")

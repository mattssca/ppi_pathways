library(dplyr)
library(tidyr)

cna_tcga = read.table(file = "../../Downloads/blca_tcga_pub_2017/data_cna.txt", sep = "\t", header = TRUE)
pred_tcga = classify_samples(this_data = TCGA_object409$salmon39_geTMM)

subtypes_tcga = as.data.frame(pred_tcga$predictions_7classes) %>% 
  rename(subtype = `pred_tcga$predictions_7classes`)

head(cna_tcga)[1:5]
head(subtypes_tcga)


#fix sample names: replace dots with hyphens
colnames(cna_tcga) <- gsub("\\.", "-", colnames(cna_tcga))

#check the fix
head(colnames(cna_tcga))

# Convert to long format for easier analysis
cna_long <- cna_tcga %>%
  select(-Entrez_Gene_Id) %>%  # Remove if not needed
  pivot_longer(
    cols = -Hugo_Symbol, 
    names_to = "sample_id", 
    values_to = "cna_value"
  ) %>%
  filter(!is.na(cna_value))  # Remove missing values

# Preview
head(cna_long)

# Fix the sample ID mismatch by truncating subtypes to 15 characters
subtypes_df <- subtypes_df %>%
  mutate(sample_id = substr(sample_id, 1, 15))

# Now merge with CNA data
cna_with_subtypes_final <- cna_long %>%
  left_join(subtypes_df, by = "sample_id") %>%
  filter(!is.na(subtype))

# Verify the fix
cat("Successfully matched samples:", nrow(cna_with_subtypes_final))
cat("\nUnique samples with subtypes:", length(unique(cna_with_subtypes_final$sample_id)))

# Check subtype distribution
table(cna_with_subtypes_final$subtype)

# Create CNA classifications
cna_patterns <- cna_with_subtypes_final %>%
  mutate(
    cna_type = case_when(
      cna_value >= 2 ~ "Amplified",
      cna_value <= -2 ~ "Deep_Deletion", 
      cna_value == 1 ~ "Gain",
      cna_value == -1 ~ "Shallow_Deletion",
      cna_value == 0 ~ "Neutral",
      TRUE ~ "Other"
    ),
    is_altered = cna_value != 0,
    is_amplified = cna_value >= 2,
    is_deleted = cna_value <= -2,
    is_significantly_altered = abs(cna_value) >= 2
  )

# Get top altered genes by subtype
alterations_by_subtype <- cna_patterns %>%
  filter(is_significantly_altered) %>%
  group_by(subtype, Hugo_Symbol) %>%
  summarise(
    n_altered = n(),
    n_amplified = sum(is_amplified),
    n_deleted = sum(is_deleted),
    .groups = "drop"
  ) %>%
  # Add total samples per subtype
  left_join(
    cna_patterns %>% 
      group_by(subtype) %>% 
      summarise(total_samples = n_distinct(sample_id), .groups = "drop"),
    by = "subtype"
  ) %>%
  mutate(
    alteration_frequency = n_altered / total_samples,
    amplification_frequency = n_amplified / total_samples,
    deletion_frequency = n_deleted / total_samples
  )

#subset to subtype specific data frames
uroa_cna = alterations_by_subtype %>% filter(subtype == "UroA")
urob_cna = alterations_by_subtype %>% filter(subtype == "UroB")
uroc_cna = alterations_by_subtype %>% filter(subtype == "UroC")
gu_cna = alterations_by_subtype %>% filter(subtype == "GU")
basq_cna = alterations_by_subtype %>% filter(subtype == "BaSq")

#mutate columns
UroA_mut_cna = uroa_cna %>% 
  select(Hugo_Symbol, alteration_frequency, amplification_frequency, deletion_frequency) %>% 
  rename(name = Hugo_Symbol, UroA_alt_freq =  alteration_frequency, UroA_amp_freq = amplification_frequency, UroA_del_freq = deletion_frequency)

UroB_mut_cna = urob_cna %>% 
  select(Hugo_Symbol, alteration_frequency, amplification_frequency, deletion_frequency) %>% 
  rename(name = Hugo_Symbol, UroB_alt_freq =  alteration_frequency, UroB_amp_freq = amplification_frequency, UroB_del_freq = deletion_frequency)

UroC_mut_cna = uroc_cna %>% 
  select(Hugo_Symbol, alteration_frequency, amplification_frequency, deletion_frequency) %>% 
  rename(name = Hugo_Symbol, UroC_alt_freq =  alteration_frequency, UroC_amp_freq = amplification_frequency, UroC_del_freq = deletion_frequency)

GU_mut_cna = gu_cna %>% 
  select(Hugo_Symbol, alteration_frequency, amplification_frequency, deletion_frequency) %>% 
  rename(name = Hugo_Symbol, GU_alt_freq =  alteration_frequency, GU_amp_freq = amplification_frequency, GU_del_freq = deletion_frequency)

BaSq_mut_cna = basq_cna %>% 
  select(Hugo_Symbol, alteration_frequency, amplification_frequency, deletion_frequency) %>% 
  rename(name = Hugo_Symbol, BaSq_alt_freq =  alteration_frequency, BaSq_amp_freq = amplification_frequency, BaSq_del_freq = deletion_frequency)

cna_tcga_combined = list(UroA_mut_cna = UroA_mut_cna,
                         UroB_mut_cna = UroB_mut_cna,
                         UroC_mut_cna = UroC_mut_cna,
                         GU_mut_cna = GU_mut_cna,
                         BaSq_mut_cna = BaSq_mut_cna
                         )

save(cna_tcga_combined, file = "ppi_pathways/results_data/cna_tcga_subtypes.Rdata")

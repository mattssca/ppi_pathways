library(dplyr)
library(tidyr)
library(data.table)

#read data
mutations_tcga <- fread(
  "../../Downloads/blca_tcga_pub_2017/data_mutations.txt",
  sep = "\t",
  header = TRUE,
  fill = TRUE,
  quote = ""
)

#subset to columns of interest
mutation_tcga_sub = mutations_tcga %>% 
  select(Hugo_Symbol, Consequence, Variant_Classification, Variant_Type, Tumor_Sample_Barcode)

# Prepare sample IDs
mutation_tcga_sub <- mutation_tcga_sub %>%
  mutate(
    sample_id = gsub("\\.", "-", Tumor_Sample_Barcode),
    sample_id = substr(sample_id, 1, 15)
  )

#run classifier
pred_tcga = classify_samples(this_data = TCGA_object409$salmon39_geTMM)

#extract subtype info
subtypes_tcga = as.data.frame(pred_tcga$predictions_7classes) %>% 
  rename(subtype = `pred_tcga$predictions_7classes`) %>% 
  rownames_to_column("sample_id")

# Fix the sample ID mismatch by truncating subtypes to 15 characters
subtypes_df <- subtypes_tcga %>%
  mutate(sample_id = substr(sample_id, 1, 15))

# Add subtype info
mutations_with_subtypes <- mutation_tcga_sub %>%
  left_join(subtypes_df, by = "sample_id") %>%
  filter(!is.na(subtype))

# Calculate frequencies within each subtype
mutations_by_subtype <- mutations_with_subtypes %>%
  group_by(subtype) %>%
  mutate(n_subtype_samples = n_distinct(sample_id)) %>%
  group_by(subtype, Hugo_Symbol) %>%
  summarise(
    del_freq = sum(Variant_Type == "DEL") / first(n_subtype_samples),
    ins_freq = sum(Variant_Type == "INS") / first(n_subtype_samples),
    snp_freq = sum(Variant_Type == "SNP") / first(n_subtype_samples),
    .groups = "drop"
  )

# Split by subtype and rename columns
UroA_mut_data = mutations_by_subtype %>% filter(subtype == "UroA") %>%
  rename(
    name = Hugo_Symbol,
    UroA_del_freq = del_freq,
    UroA_ins_freq = ins_freq,
    UroA_snp_freq = snp_freq
  )

UroB_mut_data = mutations_by_subtype %>% filter(subtype == "UroB") %>%
  rename(
    name = Hugo_Symbol,
    UroB_del_freq = del_freq,
    UroB_ins_freq = ins_freq,
    UroB_snp_freq = snp_freq
  )

UroC_mut_data = mutations_by_subtype %>% filter(subtype == "UroC") %>%
  rename(
    name = Hugo_Symbol,
    UroC_del_freq = del_freq,
    UroC_ins_freq = ins_freq,
    UroC_snp_freq = snp_freq
  )

GU_mut_data = mutations_by_subtype %>% filter(subtype == "GU") %>%
  rename(
    name = Hugo_Symbol,
    GU_del_freq = del_freq,
    GU_ins_freq = ins_freq,
    GU_snp_freq = snp_freq
  )

BaSq_mut_data = mutations_by_subtype %>% filter(subtype == "BaSq") %>%
  rename(
    name = Hugo_Symbol,
    BaSq_del_freq = del_freq,
    BaSq_ins_freq = ins_freq,
    BaSq_snp_freq = snp_freq
  )

mutations_tcga_combined = list(
  UroA_mut_data = UroA_mut_data,
  UroB_mut_data = UroB_mut_data,
  UroC_mut_data = UroC_mut_data,
  GU_mut_data = GU_mut_data,
  BaSq_mut_data = BaSq_mut_data
)

# Save the results
save(mutations_tcga_combined, file = "ppi_pathways/results_data/mutations_tcga_subtypes_freqs.Rdata")

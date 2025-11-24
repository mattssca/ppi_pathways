#load packages
library(dplyr)

#read and clean
mutations <- read.table(
  "data/TCGA_mutated_genes.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

#clean column names
names(mutations) <- c(
  "name",
  "tcga_mutsig_q_value",
  "tcga_num_mutations",
  "tcga_num_samples_mutated",
  "tcga_profiled_samples",
  "tcga_frequency",
  "tcga_is_cancer_gene"
)

#convert frequency from "1.5%" to 0.015
tcga_mutations <- mutations %>%
  mutate(
    tcga_mutation_frequency = as.numeric(gsub("%", "", tcga_frequency)) / 100,
    tcga_is_cancer_gene = (tcga_is_cancer_gene == "Yes")
  )

#read CNA data
cna <- read.table(
  "data/TCGA_cna_genes.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

#clean column names
names(cna) <- c(
  "name",
  "tcga_gistic_q_value",
  "tcga_cytoband",
  "tcga_cna_type",
  "tcga_cna_samples",
  "tcga_cna_profiled_samples",
  "tcga_cna_frequency",
  "tcga_cna_is_cancer_gene"
)

#convert frequency and process
tcga_cna <- cna %>%
  group_by(name) %>%
  summarise(
    is_amplified = any(tcga_cna_type == "AMP"),
    is_homdel = any(tcga_cna_type == "HOMDEL"),
    .groups = "drop"
  )

#view cleaned data
head(tcga_cna)

#join the data
tcga_mut_data = tcga_cna %>% 
  left_join(tcga_mutations, by = "name")

#export data
save(tcga_mut_data, file = "results_data/tcga_mut_data.Rdata")

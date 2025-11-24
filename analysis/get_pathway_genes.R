#load pacakges
library(dplyr)
library(biomaRt)

#download Reactome pathway-gene mapping
download.file(
  url = "https://reactome.org/download/current/Ensembl2Reactome.txt",
  destfile = "../GIT_REPOS/ppi_signal/data/Ensembl2Reactome.txt"
)

#read the data
reactome_data <- read.delim("../GIT_REPOS/ppi_signal/data/Ensembl2Reactome.txt", 
                            header = FALSE, 
                            stringsAsFactors = FALSE)

colnames(reactome_data) <- c("gene_id", "pathway_id", "url", "pathway_name", "evidence", "species")

#filter for human pathways only
reactome_data <- reactome_data[reactome_data$species == "Homo sapiens", ]

#######################################################################################
#EGFR Pathway
egfr_all <- reactome_data %>%
  filter(grepl("EGFR", pathway_name, ignore.case = TRUE))

cat("EGFR: Found genes from these pathways:\n")
print(unique(egfr_all$pathway_name))

#######################################################################################
#FGFR3 pathway
fgfr3_all <- reactome_data %>%
  filter(grepl("FGFR3", pathway_name, ignore.case = TRUE))

cat("\n\nFGFR3: Found genes from these pathways:\n")
print(unique(fgfr3_all$pathway_name))

#######################################################################################
#ERBB2 GENES
erbb2_all <- reactome_data %>%
  filter(grepl("ERBB2", pathway_name, ignore.case = TRUE))

cat("\n\nRRBB2: Found genes from these pathways:\n")
print(unique(erbb2_all$pathway_name))

#######################################################################################
#convert Ensembl IDs to gene symbols using biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#EGFR genes
egfr_ensembl <- unique(egfr_all$gene_id)

egfr_symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      filters = "ensembl_gene_id",
                      values = egfr_ensembl,
                      mart = mart)

egfr_genes <- sort(unique(egfr_symbols$hgnc_symbol[egfr_symbols$hgnc_symbol != ""]))

#FGFR3 genes
fgfr3_ensembl <- unique(fgfr3_all$gene_id)

fgfr3_symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                       filters = "ensembl_gene_id",
                       values = fgfr3_ensembl,
                       mart = mart)

fgfr3_genes <- sort(unique(fgfr3_symbols$hgnc_symbol[fgfr3_symbols$hgnc_symbol != ""]))

#ERBB2 genes
erbb2_ensembl <- unique(erbb2_all$gene_id)

erbb2_symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      filters = "ensembl_gene_id",
                      values = erbb2_ensembl,
                      mart = mart)
erbb2_genes <- sort(unique(erbb2_symbols$hgnc_symbol[erbb2_symbols$hgnc_symbol != ""]))

#results
cat("\nEGFR pathway genes (", length(egfr_genes), "):\n", sep = "")
print(egfr_genes)

cat("\n\nFGFR3 pathway genes (", length(fgfr3_genes), "):\n", sep = "")
print(fgfr3_genes)

cat("\n\nERBB2 pathway genes (", length(erbb2_genes), "):\n", sep = "")
print(erbb2_genes)

#save results
pathway_results <- list(
  egfr_genes = egfr_genes,
  fgfr3_genes = fgfr3_genes,
  erbb2_genes = erbb2_genes
)

save(pathway_results, file =  "../GIT_REPOS/ppi_signal/results_data/pathway_results.Rdata")

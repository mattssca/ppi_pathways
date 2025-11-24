#load data
load("ppi_signal/results_data/egfr_nodes_object.Rdata")
load("ppi_signal/results_data/fgfr3_nodes_object.Rdata")
load("ppi_signal/results_data/erbb2_nodes_object.Rdata")

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)

# Create functional category mapping function
create_functional_categories <- function(go_term) {
  dplyr::case_when(
    grepl("cytoplasmic translation|ribosom|translation", go_term, ignore.case = TRUE) ~ "Translation",
    grepl("ERK|MAPK|cascade", go_term, ignore.case = TRUE) ~ "MAPK_Signaling",
    grepl("ERBB|peptide|hormone|growth", go_term, ignore.case = TRUE) ~ "Growth_Factor_Response",
    grepl("phosphorylation|tyrosine", go_term, ignore.case = TRUE) ~ "Phosphorylation",
    grepl("transferase|kinase activity", go_term, ignore.case = TRUE) ~ "Transferase_Regulation",
    TRUE ~ "Other"
  )
}

# Create pathway role mapping function
create_pathway_role <- function(functional_category, degree, betweenness) {
  dplyr::case_when(
    degree >= quantile(degree, 0.8, na.rm = TRUE) & functional_category != "Other" ~ "Functional_Hub",
    betweenness >= quantile(betweenness, 0.8, na.rm = TRUE) & functional_category != "Other" ~ "Functional_Bridge",
    functional_category %in% c("Translation", "MAPK_Signaling") ~ "Core_Process",
    functional_category == "Growth_Factor_Response" ~ "Signal_Receptor",
    functional_category == "Phosphorylation" ~ "Signal_Modifier",
    TRUE ~ "Peripheral"
  )
}

# Function to process functional annotations for a single pathway
process_pathway_functional_annotations <- function(metrics_data, pathway_name) {

  cat("Processing", pathway_name, "network...\n")

  # Gene symbol to Entrez ID conversion
  gene_entrez <- bitr(metrics_data$name, fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db, drop = FALSE)

  # GO enrichment analysis
  go_enrichment <- enrichGO(gene = gene_entrez$ENTREZID,
                            OrgDb = org.Hs.eg.db,
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)

  go_results <- as.data.frame(go_enrichment@result)
  significant_terms <- go_results[go_results$p.adjust < 0.05, ]

  cat("Found", nrow(significant_terms), "significant GO terms\n")

  # Create gene-GO term mappings
  gene_modules <- significant_terms %>%
    dplyr::mutate(
      network = pathway_name,
      genes_entrez = strsplit(geneID, "/")
    ) %>%
    tidyr::unnest(genes_entrez) %>%
    dplyr::left_join(gene_entrez, by = c("genes_entrez" = "ENTREZID")) %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::select(gene = SYMBOL, go_id = ID, go_term = Description, network,
                  p_value = pvalue, p_adjust = p.adjust)

  # Calculate functional scores per gene
  gene_functional_scores <- gene_modules %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      min_p_adjust = min(p_adjust, na.rm = TRUE),
      functional_diversity = dplyr::n_distinct(go_term),
      top_go_term = go_term[which.min(p_adjust)][1],
      .groups = 'drop'
    ) %>%
    dplyr::mutate(
      functional_significance = -log10(min_p_adjust),
      functional_category = create_functional_categories(top_go_term),
      functional_specificity_ratio = functional_significance / median(functional_significance, na.rm = TRUE)
    )

  # Create functional annotations for this pathway
  pathway_functional_annotations <- metrics_data %>%
    dplyr::left_join(gene_functional_scores %>% dplyr::select(-top_go_term), by = c("name" = "gene")) %>%
    dplyr::mutate(
      # Handle genes not in functional analysis
      functional_category = ifelse(is.na(functional_category), "Unannotated", functional_category),
      functional_significance = ifelse(is.na(functional_significance), 0, functional_significance),
      functional_diversity = ifelse(is.na(functional_diversity), 0, functional_diversity),
      functional_specificity_ratio = ifelse(is.na(functional_specificity_ratio), 0, functional_specificity_ratio),

      # Create pathway role based on network metrics + functional category
      pathway_role = create_pathway_role(functional_category, degree, betweenness),

      # Create combined functional-structural importance score
      functional_structural_importance = scale(functional_significance)[,1] + scale(degree)[,1] + scale(betweenness)[,1]
    ) %>%
    dplyr::select(name, min_p_adjust, functional_diversity, functional_significance,
                  functional_category, functional_specificity_ratio, pathway_role,
                  functional_structural_importance)

  # Summary
  cat("Functional categories for", pathway_name, ":\n")
  print(table(pathway_functional_annotations$functional_category))
  cat("\n")

  return(list(
    functional_annotations = pathway_functional_annotations,
    gene_modules = gene_modules,
    go_results = go_results
  ))
}

# Process each pathway individually
cat("=== Processing EGFR Network ===\n")
load("ppi_signal/results_data/egfr_nodes_object.Rdata")
egfr_results <- process_pathway_functional_annotations(egfr_metrics, "EGFR")
egfr_functional_annotations <- egfr_results$functional_annotations
save(egfr_functional_annotations, file = "ppi_signal/results_data/egfr_functional_annotations.Rdata")

cat("=== Processing FGFR3 Network ===\n")
load("ppi_signal/results_data/fgfr3_nodes_object.Rdata")
fgfr3_results <- process_pathway_functional_annotations(fgfr3_metrics, "FGFR3")
fgfr3_functional_annotations <- fgfr3_results$functional_annotations
save(fgfr3_functional_annotations, file = "ppi_signal/results_data/fgfr3_functional_annotations.Rdata")

cat("=== Processing ERBB2 Network ===\n")
load("ppi_signal/results_data/erbb2_nodes_object.Rdata")
erbb2_results <- process_pathway_functional_annotations(erbb2_metrics, "ERBB2")
erbb2_functional_annotations <- erbb2_results$functional_annotations
save(erbb2_functional_annotations, file = "ppi_signal/results_data/erbb2_functional_annotations.Rdata")

cat("=== All Pathway Processing Complete ===\n")
cat("Saved individual functional annotation files:\n")
cat("- egfr_functional_annotations.Rdata\n")
cat("- fgfr3_functional_annotations.Rdata\n")
cat("- erbb2_functional_annotations.Rdata\n")

# Optional: Create summary comparison
cat("\n=== Summary Across All Pathways ===\n")
all_results <- list(
  EGFR = egfr_results,
  FGFR3 = fgfr3_results,
  ERBB2 = erbb2_results
)

for(pathway in names(all_results)) {
  annotated_genes <- sum(all_results[[pathway]]$functional_annotations$functional_category != "Unannotated")
  total_genes <- nrow(all_results[[pathway]]$functional_annotations)
  cat(pathway, ": ", annotated_genes, "/", total_genes, " genes annotated (",
      round(annotated_genes/total_genes*100, 1), "%)\n", sep="")
}

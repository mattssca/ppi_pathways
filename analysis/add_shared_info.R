library(dplyr)

# Step 1: Get all unique genes from each pathway
egfr_genes <- unique(egfr_meta_consoldiated$name)
fgfr3_genes <- unique(fgfr3_meta_consoldiated$name)
erbb2_genes <- unique(erbb2_meta_consoldiated$name)

# Step 2: Create a function to classify gene sharing
classify_gene_sharing <- function(gene, pathway) {
  in_egfr <- gene %in% egfr_genes
  in_fgfr3 <- gene %in% fgfr3_genes
  in_erbb2 <- gene %in% erbb2_genes
  
  # Count how many pathways the gene appears in
  n_pathways <- sum(in_egfr, in_fgfr3, in_erbb2)
  
  # Determine sharing pattern
  if (n_pathways == 1) {
    return(paste0("Unique_to_", pathway))
  } else if (n_pathways == 3) {
    return("Shared_all_3")
  } else if (n_pathways == 2) {
    # Determine which two pathways
    if (in_egfr && in_fgfr3) {
      return("Shared_EGFR_FGFR3")
    } else if (in_egfr && in_erbb2) {
      return("Shared_EGFR_ERBB2")
    } else if (in_fgfr3 && in_erbb2) {
      return("Shared_FGFR3_ERBB2")
    }
  }
}

# Step 3: Add the new column to each dataframe
egfr_meta_consoldiated <- egfr_meta_consoldiated %>%
  mutate(
    pathway_sharing = sapply(name, classify_gene_sharing, pathway = "EGFR"),
    n_pathways = sapply(name, function(g) sum(g %in% egfr_genes, g %in% fgfr3_genes, g %in% erbb2_genes)),
    pathways_list = sapply(name, function(g) {
      pathways <- c()
      if (g %in% egfr_genes) pathways <- c(pathways, "EGFR")
      if (g %in% fgfr3_genes) pathways <- c(pathways, "FGFR3")
      if (g %in% erbb2_genes) pathways <- c(pathways, "ERBB2")
      paste(pathways, collapse = ",")
    })
  )

fgfr3_meta_consoldiated <- fgfr3_meta_consoldiated %>%
  mutate(
    pathway_sharing = sapply(name, classify_gene_sharing, pathway = "FGFR3"),
    n_pathways = sapply(name, function(g) sum(g %in% egfr_genes, g %in% fgfr3_genes, g %in% erbb2_genes)),
    pathways_list = sapply(name, function(g) {
      pathways <- c()
      if (g %in% egfr_genes) pathways <- c(pathways, "EGFR")
      if (g %in% fgfr3_genes) pathways <- c(pathways, "FGFR3")
      if (g %in% erbb2_genes) pathways <- c(pathways, "ERBB2")
      paste(pathways, collapse = ",")
    })
  )

erbb2_meta_consoldiated <- erbb2_meta_consoldiated %>%
  mutate(
    pathway_sharing = sapply(name, classify_gene_sharing, pathway = "ERBB2"),
    n_pathways = sapply(name, function(g) sum(g %in% egfr_genes, g %in% fgfr3_genes, g %in% erbb2_genes)),
    pathways_list = sapply(name, function(g) {
      pathways <- c()
      if (g %in% egfr_genes) pathways <- c(pathways, "EGFR")
      if (g %in% fgfr3_genes) pathways <- c(pathways, "FGFR3")
      if (g %in% erbb2_genes) pathways <- c(pathways, "ERBB2")
      paste(pathways, collapse = ",")
    })
  )

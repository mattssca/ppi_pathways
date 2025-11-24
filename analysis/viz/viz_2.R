# =============================================================================
# SUBTYPE-SPECIFIC FUNCTIONAL ENRICHMENT ANALYSIS
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(broom)

# =============================================================================
# 1. IDENTIFY SUBTYPE-SPECIFIC DIFFERENTIALLY EXPRESSED GENES
# =============================================================================

# Create subtype-specific gene rankings based on individual subtype expression
create_subtype_rankings <- function(data) {

  # Define subtypes
  subtypes <- c("UroB", "GU", "UroA", "BaSq", "UroC")

  # For each subtype, rank genes by their expression in that subtype
  subtype_rankings <- list()

  for(subtype in subtypes) {
    col_name <- paste0("mean_expr_", subtype)

    # Create ranking for this subtype
    ranking <- data %>%
      arrange(desc(!!sym(col_name))) %>%
      mutate(
        subtype = subtype,
        rank_in_subtype = row_number(),
        expression_value = !!sym(col_name),
        # Create categories based on expression levels
        expression_category = case_when(
          rank_in_subtype <= 50 ~ "Top_50_Expressed",
          rank_in_subtype <= 100 ~ "Top_100_Expressed",
          rank_in_subtype <= 200 ~ "Top_200_Expressed",
          expression_value > 0.5 ~ "Highly_Expressed",
          expression_value > 0 ~ "Moderately_Expressed",
          expression_value > -0.5 ~ "Low_Expressed",
          TRUE ~ "Very_Low_Expressed"
        )
      ) %>%
      dplyr::select(name, subtype, rank_in_subtype, expression_value, expression_category,
             functional_category, pathway_role, functional_significance,
             functional_structural_importance, sig_and_meaningful, effect_size_category)

    subtype_rankings[[subtype]] <- ranking
  }

  return(bind_rows(subtype_rankings))
}

# Create subtype rankings
subtype_data <- create_subtype_rankings(all_metrics)

# =============================================================================
# 2. SUBTYPE-SPECIFIC FUNCTIONAL ENRICHMENT TESTING
# =============================================================================

# Test for functional category enrichment in top expressed genes per subtype
test_functional_enrichment <- function(subtype_data, expression_cutoff = "Top_100_Expressed") {

  enrichment_results <- list()
  subtypes <- unique(subtype_data$subtype)

  for(st in subtypes) {

    subtype_subset <- subtype_data %>% filter(subtype == st)

    # Create contingency table for chi-square test
    # Top expressed vs. rest, by functional category
    contingency_data <- subtype_subset %>%
      mutate(
        is_top_expressed = expression_category == expression_cutoff,
        functional_cat_clean = ifelse(functional_category == "Unannotated", "Other", functional_category)
      ) %>%
      count(functional_cat_clean, is_top_expressed) %>%
      pivot_wider(names_from = is_top_expressed, values_from = n, values_fill = 0)

    # Calculate enrichment for each functional category
    func_enrichment <- contingency_data %>%
      mutate(
        total_in_category = `TRUE` + `FALSE`,
        total_top_expressed = sum(contingency_data$`TRUE`),
        total_genes = sum(contingency_data$`TRUE`) + sum(contingency_data$`FALSE`),
        expected_in_top = (total_in_category * total_top_expressed) / total_genes,
        fold_enrichment = `TRUE` / expected_in_top,
        percentage_in_top = (`TRUE` / total_top_expressed) * 100,
        percentage_of_category = (`TRUE` / total_in_category) * 100
      ) %>%
      filter(total_in_category >= 5) %>%  # Only test categories with ≥5 genes
      arrange(desc(fold_enrichment))

    func_enrichment$subtype <- st
    func_enrichment$expression_cutoff <- expression_cutoff

    enrichment_results[[st]] <- func_enrichment
  }

  return(bind_rows(enrichment_results))
}

# Test enrichment for different expression cutoffs
top50_enrichment <- test_functional_enrichment(subtype_data, "Top_50_Expressed")
top100_enrichment <- test_functional_enrichment(subtype_data, "Top_100_Expressed")
highly_enrichment <- test_functional_enrichment(subtype_data, "Highly_Expressed")

# =============================================================================
# 3. ANOVA-BASED SUBTYPE SPECIFICITY ANALYSIS
# =============================================================================

# Focus on genes with significant ANOVA results (differential across subtypes)
analyze_anova_functional_patterns <- function(data) {

  # Create categories based on ANOVA significance and effect sizes
  anova_categories <- data %>%
    mutate(
      anova_category = case_when(
        sig_and_meaningful & effect_size_category %in% c("Large", "Very Large") ~ "High_Subtype_Specificity",
        sig_and_meaningful & effect_size_category == "Medium" ~ "Medium_Subtype_Specificity",
        significant_05 ~ "Low_Subtype_Specificity",
        TRUE ~ "No_Subtype_Specificity"
      )
    )

  # Test functional enrichment in ANOVA categories
  anova_func_enrichment <- anova_categories %>%
    count(anova_category, functional_category) %>%
    group_by(anova_category) %>%
    mutate(
      total_in_anova_cat = sum(n),
      percentage = (n / total_in_anova_cat) * 100
    ) %>%
    ungroup()

  return(list(
    categorized_data = anova_categories,
    enrichment_results = anova_func_enrichment
  ))
}

anova_analysis <- analyze_anova_functional_patterns(all_metrics)

# =============================================================================
# 4. IDENTIFY SUBTYPE-SPECIFIC FUNCTIONAL SIGNATURES
# =============================================================================

# Alternative version using dplyr::first explicitly
find_subtype_specific_genes <- function(subtype_data, top_n = 50) {

  subtype_specific <- list()
  subtypes <- unique(subtype_data$subtype)

  for(st in subtypes) {

    # Get top N genes in this subtype
    top_in_subtype <- subtype_data %>%
      filter(subtype == st, rank_in_subtype <= top_n) %>%
      pull(name)

    # Check their ranks in other subtypes
    specificity_analysis <- subtype_data %>%
      filter(name %in% top_in_subtype) %>%
      group_by(name) %>%
      summarise(
        target_subtype = st,
        target_rank = rank_in_subtype[subtype == st],
        target_expression = expression_value[subtype == st],
        avg_rank_other_subtypes = mean(rank_in_subtype[subtype != st]),
        max_rank_other_subtypes = max(rank_in_subtype[subtype != st]),
        functional_category = dplyr::first(functional_category[subtype == st]),
        pathway_role = dplyr::first(pathway_role[subtype == st]),
        functional_significance = dplyr::first(functional_significance[subtype == st]),
        .groups = 'drop'
      ) %>%
      mutate(
        specificity_score = avg_rank_other_subtypes - target_rank,
        is_highly_specific = specificity_score > 100  # Much lower rank in other subtypes
      ) %>%
      arrange(desc(specificity_score))

    subtype_specific[[st]] <- specificity_analysis
  }

  return(bind_rows(subtype_specific))
}

# Try the corrected version
subtype_specific_genes <- find_subtype_specific_genes(subtype_data, top_n = 100)

# =============================================================================
# 5. VISUALIZATION: SUBTYPE-SPECIFIC FUNCTIONAL LANDSCAPES
# =============================================================================

# A. Heatmap of functional enrichment across subtypes
p1_enrichment_heatmap <- top100_enrichment %>%
  dplyr::select(functional_cat_clean, subtype, fold_enrichment) %>%
  pivot_wider(names_from = subtype, values_from = fold_enrichment, values_fill = 1) %>%
  column_to_rownames("functional_cat_clean") %>%
  as.matrix()

# Create heatmap
library(pheatmap)
pheatmap(p1_enrichment_heatmap,
         scale = "none",
         color = colorRampPalette(c("white", "red", "darkred"))(100),
         main = "Functional Category Enrichment in Top 100 Expressed Genes\nby Subtype",
         fontsize = 10,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# B. Bar plot of ANOVA-based functional enrichment
p2_anova_enrichment <- anova_analysis$enrichment_results %>%
  filter(anova_category %in% c("High_Subtype_Specificity", "Medium_Subtype_Specificity")) %>%
  ggplot(aes(x = reorder(functional_category, percentage), y = percentage,
             fill = anova_category)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Functional Enrichment in Subtype-Specific Genes",
       subtitle = "Based on ANOVA significance and effect sizes",
       x = "Functional Category", y = "Percentage of Genes",
       fill = "Subtype Specificity") +
  theme_minimal() +
  scale_fill_viridis_d()

print(p2_anova_enrichment)

# C. Subtype-specific gene signatures
p3_specific_signatures <- subtype_specific_genes %>%
  filter(is_highly_specific) %>%
  count(target_subtype, functional_category) %>%
  ggplot(aes(x = target_subtype, y = n, fill = functional_category)) +
  geom_bar(stat = "identity") +
  labs(title = "Subtype-Specific Gene Signatures",
       subtitle = "Functional composition of subtype-specific highly expressed genes",
       x = "Cancer Subtype", y = "Number of Specific Genes",
       fill = "Functional Category") +
  theme_minimal() +
  scale_fill_brewer(type = "qual", palette = "Set2")

print(p3_specific_signatures)

# =============================================================================
# 6. STATISTICAL ANALYSIS: SUBTYPE vs FUNCTION RELATIONSHIPS
# =============================================================================

# Test if functional categories show different expression patterns across subtypes
test_subtype_functional_differences <- function(data) {

  # Prepare data for statistical testing
  expr_long <- data %>%
    dplyr::select(name, functional_category, mean_expr_UroB, mean_expr_GU,
           mean_expr_UroA, mean_expr_BaSq, mean_expr_UroC) %>%
    pivot_longer(cols = starts_with("mean_expr_"),
                 names_to = "subtype", values_to = "expression") %>%
    mutate(subtype = gsub("mean_expr_", "", subtype))

  # ANOVA for each functional category
  func_categories <- unique(expr_long$functional_category)
  func_categories <- func_categories[func_categories != "Unannotated"]

  anova_results <- list()

  for(func_cat in func_categories) {

    func_data <- expr_long %>% filter(functional_category == func_cat)

    if(nrow(func_data) > 20) {  # Only test categories with sufficient genes

      anova_result <- aov(expression ~ subtype, data = func_data)
      anova_summary <- tidy(anova_result)

      # Post-hoc tests if significant
      if(anova_summary$p.value[1] < 0.05) {
        posthoc <- TukeyHSD(anova_result)
        posthoc_df <- tidy(posthoc) %>%
          filter(adj.p.value < 0.05) %>%
          arrange(adj.p.value)
      } else {
        posthoc_df <- NULL
      }

      anova_results[[func_cat]] <- list(
        functional_category = func_cat,
        anova_summary = anova_summary,
        significant = anova_summary$p.value[1] < 0.05,
        posthoc_results = posthoc_df
      )
    }
  }

  return(anova_results)
}

functional_anova_results <- test_subtype_functional_differences(all_metrics)


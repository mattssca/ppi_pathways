#load data
load("ppi_signal/results_data/node_metrics_combined.Rdata")
load("ppi_signal/data/uroscanseq_data_Rdata")

#source scripts
source("ppi_signal/functions/run_anova.R")

#get uniwue genes
unique_network_genes <- unique(node_metrics_combined$name)

#run analysis with effect sizes
anova_with_effects <- calculate_effect_sizes(expr_data = uroscanseq_data$expr_df,
                                             subtype_vector = uroscanseq_data$subtype_7_vector,
                                             genes = unique_network_genes)

#add corrections and thresholds
anova_with_effects <- anova_with_effects %>%
  filter(!is.na(p_value)) %>%
  mutate(
    #multiple testing correction
    p_adjusted_bh = p.adjust(p_value, method = "BH"),
    p_adjusted_bonf = p.adjust(p_value, method = "bonferroni"),

    #statistical significance
    significant_05 = p_adjusted_bh < 0.05,
    significant_01 = p_adjusted_bh < 0.01,

    #effect size categories (Cohen's conventions)
    effect_size_category = case_when(
      is.na(eta_squared) ~ "Unknown",
      eta_squared < 0.01 ~ "Negligible",      # < 1% variance explained
      eta_squared < 0.06 ~ "Small",           # 1-6% variance explained
      eta_squared < 0.14 ~ "Medium",          # 6-14% variance explained
      eta_squared >= 0.14 ~ "Large",          # >= 14% variance explained
      TRUE ~ "Unknown"
    ),

    #cohen's f categories
    cohens_f_category = case_when(
      is.na(cohens_f) ~ "Unknown",
      cohens_f < 0.1 ~ "Small",
      cohens_f < 0.25 ~ "Medium",
      cohens_f < 0.4 ~ "Large",
      cohens_f >= 0.4 ~ "Very Large",
      TRUE ~ "Unknown"
    ),

    #combined significance and effect size
    sig_and_meaningful = significant_05 & eta_squared >= 0.01,  #significant + at least 1% variance
    sig_and_large_effect = significant_05 & eta_squared >= 0.06  #significant + medium/large effect
  ) %>%
  arrange(desc(eta_squared))  #order by effect size


#print interpretation guide
cat("Effect Size Interpretation (eta-squared):\n")
cat("Negligible: < 1% of variance explained\n")
cat("Small: 1-6% of variance explained\n")
cat("Medium: 6-14% of variance explained\n")
cat("Large: >= 14% of variance explained\n\n")

#summary of results
summary_stats <- anova_with_effects %>%
  summarise(
    total_genes = n(),
    significant_05 = sum(significant_05, na.rm = TRUE),
    sig_with_small_effect = sum(sig_and_meaningful, na.rm = TRUE),
    sig_with_large_effect = sum(sig_and_large_effect, na.rm = TRUE),
    mean_eta_squared = mean(eta_squared, na.rm = TRUE),
    median_eta_squared = median(eta_squared, na.rm = TRUE)
  )

print(summary_stats)

#effect size distribution
cat("\nEffect size categories:\n")
print(table(anova_with_effects$effect_size_category, useNA = "ifany"))

#save results
save(anova_with_effects, file = "ppi_signal/results_data/anova_results.Rdata")

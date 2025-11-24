library(dplyr)
library(tidyr)

#function to calculate effect sizes for ANOVA
calculate_effect_sizes <- function(expr_data, subtype_vector, genes) {

  #filter to genes of interest
  expr_subset <- expr_data[rownames(expr_data) %in% genes, ]

  #convert to long format
  expr_long <- expr_subset %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather(key = "sample_id", value = "expression", -gene) %>%
    mutate(subtype = subtype_vector[sample_id]) %>%
    filter(!is.na(subtype))

  #calculate effect sizes for each gene
  results <- expr_long %>%
    group_by(gene) %>%
    do({
      tryCatch({
        #ANOVA
        anova_result <- aov(expression ~ subtype, data = .)
        anova_summary <- summary(anova_result)

        #Basic ANOVA results
        f_stat <- anova_summary[[1]]["subtype", "F value"]
        p_val <- anova_summary[[1]]["subtype", "Pr(>F)"]
        df_between <- anova_summary[[1]]["subtype", "Df"]
        df_within <- anova_summary[[1]]["Residuals", "Df"]

        #Effect size calculations
        ss_between <- anova_summary[[1]]["subtype", "Sum Sq"]
        ss_within <- anova_summary[[1]]["Residuals", "Sum Sq"]
        ss_total <- ss_between + ss_within

        #Eta-squared (proportion of variance explained)
        eta_squared <- ss_between / ss_total

        #Partial eta-squared (more commonly used)
        partial_eta_squared <- ss_between / (ss_between + ss_within)

        #Omega-squared (less biased than eta-squared)
        ms_between <- anova_summary[[1]]["subtype", "Mean Sq"]
        ms_within <- anova_summary[[1]]["Residuals", "Mean Sq"]
        omega_squared <- (ss_between - df_between * ms_within) / (ss_total + ms_within)
        omega_squared <- max(0, omega_squared)  # Can't be negative

        #Cohen's f (effect size for ANOVA)
        cohens_f <- sqrt(eta_squared / (1 - eta_squared))

        #Range of expression across subtypes (simple measure)
        subtype_means <- tapply(.$expression, .$subtype, mean, na.rm = TRUE)
        expression_range <- max(subtype_means, na.rm = TRUE) - min(subtype_means, na.rm = TRUE)

        #Standard deviation for context
        overall_sd <- sd(.$expression, na.rm = TRUE)

        data.frame(
          f_statistic = f_stat,
          p_value = p_val,
          df_between = df_between,
          df_within = df_within,
          eta_squared = eta_squared,
          partial_eta_squared = partial_eta_squared,
          omega_squared = omega_squared,
          cohens_f = cohens_f,
          expression_range = expression_range,
          overall_sd = overall_sd,
          effect_size_ratio = expression_range / overall_sd  #Range in SD units
        )
      }, error = function(e) {
        data.frame(
          f_statistic = NA, p_value = NA, df_between = NA, df_within = NA,
          eta_squared = NA, partial_eta_squared = NA, omega_squared = NA,
          cohens_f = NA, expression_range = NA, overall_sd = NA, effect_size_ratio = NA
        )
      })
    }) %>%
    ungroup()

  return(results)
}

#function to add HPA annotations and filter flags (NO FILTERING)
add_hpa_annotations <- function(pathway_metrics,
                                normal_ihc,
                                cancer_ihc) {

  hpa_data <- pathway_metrics %>%
    select(name) %>%
    distinct() %>%
    left_join(normal_ihc, by = "name") %>%
    left_join(cancer_ihc, by = "name") %>%

    #add all classification columns
    mutate(
      #classify normal expression
      hpa_normal_expressed = case_when(
        hpa_normal_level %in% c("High", "Medium") ~ TRUE,
        hpa_normal_level == "Low" ~ FALSE,
        hpa_normal_level == "Not detected" ~ FALSE,
        is.na(hpa_normal_level) ~ NA
      ),

      #classify normal reliability
      hpa_normal_reliable = case_when(
        hpa_normal_reliability %in% c("Enhanced", "Supported") ~ "High",
        hpa_normal_reliability == "Approved" ~ "Medium",
        hpa_normal_reliability == "Uncertain" ~ "Low",
        is.na(hpa_normal_reliability) ~ NA_character_
      ),

      #classify cancer expression
      hpa_cancer_expressed = case_when(
        (hpa_cancer_high + hpa_cancer_medium) > 0 ~ TRUE,
        hpa_cancer_low > 0 ~ FALSE,
        hpa_cancer_not_detected >= 10 ~ FALSE,
        is.na(hpa_cancer_high) ~ NA
      ),

      #cancer expression confidence
      hpa_cancer_confidence = case_when(
        (hpa_cancer_high + hpa_cancer_medium) >= 6 ~ "High",
        (hpa_cancer_high + hpa_cancer_medium) >= 3 ~ "Medium",
        (hpa_cancer_high + hpa_cancer_medium) >= 1 ~ "Low",
        TRUE ~ "None"
      ),

      #overall expression pattern
      hpa_expression_pattern = case_when(
        hpa_normal_expressed & hpa_cancer_expressed ~ "Constitutive",
        hpa_normal_expressed & !hpa_cancer_expressed ~ "Lost_in_cancer",
        !hpa_normal_expressed & hpa_cancer_expressed ~ "Cancer_upregulated",
        !hpa_normal_expressed & !hpa_cancer_expressed ~ "Not_expressed",
        is.na(hpa_normal_expressed) | is.na(hpa_cancer_expressed) ~ "Insufficient_data",
        TRUE ~ "Unknown"
      ),

      #overall data quality
      hpa_data_quality = case_when(
        !is.na(hpa_normal_level) & !is.na(hpa_cancer_high) &
          hpa_normal_reliable %in% c("High", "Medium") ~ "High",
        !is.na(hpa_normal_level) & !is.na(hpa_cancer_high) ~ "Medium",
        is.na(hpa_normal_level) | is.na(hpa_cancer_high) ~ "Low",
        TRUE ~ "None"
      ),

      #FILTER FLAGS
      #filter 1: Has both normal AND cancer data
      hpa_filter_has_data = hpa_data_quality %in% c("High", "Medium"),

      #filter 2a: Expressed in normal urothelium with good antibody
      hpa_filter_normal_quality = (
        hpa_normal_expressed == TRUE &
          hpa_normal_reliable %in% c("High", "Medium")
      ),

      #filter 2b: Expressed in at least 2 cancer patients (High/Medium)
      hpa_filter_cancer_quality = hpa_cancer_confidence %in% c("High", "Medium"),

      #combined filter: Passes overall quality criteria
      hpa_passes_filter = hpa_filter_has_data &
        (hpa_filter_normal_quality | hpa_filter_cancer_quality),

      #additional useful flags
      hpa_high_confidence = (
        hpa_data_quality == "High" &
          hpa_cancer_confidence == "High" &
          hpa_normal_reliable == "High"
      ),

      hpa_therapeutic_target = (
        hpa_expression_pattern == "Cancer_upregulated" &
          hpa_cancer_confidence %in% c("High", "Medium")
      ),

      hpa_potential_suppressor = (
        hpa_expression_pattern == "Lost_in_cancer" &
          hpa_normal_expressed == TRUE
      )
    )

  return(hpa_data)
}

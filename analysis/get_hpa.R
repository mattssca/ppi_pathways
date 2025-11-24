#load packages
library(dplyr)
library(readr)
library(tidyr)

#source functions
source("functions/download_hpa_data.R")
source("functions/add_hpa_annotations.R")

#NORMAL TISSUE IHC DATA
normal_ihc_url <- "https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip"
normal_ihc <- download_hpa_data(normal_ihc_url, "Normal Tissue IHC Data")

bladder_normal_ihc <- normal_ihc %>%
  filter(Tissue == "Urinary bladder") %>%
  rename(
    name = `Gene name`,
    hpa_normal_tissue = Tissue,
    hpa_normal_cell_type = `Cell type`,
    hpa_normal_level = Level,
    hpa_normal_reliability = Reliability
  ) %>%
  dplyr::select(-Gene, `IHC tissue name`)

#CANCER DATA
cancer_url <- "https://www.proteinatlas.org/download/tsv/cancer_data.tsv.zip"
cancer_data <- download_hpa_data(cancer_url, "Cancer Data")

bladder_cancer_ihc <- cancer_data %>%
  filter(Cancer == "urothelial cancer") %>%
  rename(
    name = `Gene name`,
    hpa_cancer_type = Cancer,
    hpa_cancer_high = High,
    hpa_cancer_medium = Medium,
    hpa_cancer_low = Low,
    hpa_cancer_not_detected = `Not detected`
  ) %>%
  dplyr::select(-Gene)

#export objects
save(bladder_normal_ihc, file = "results_data/bladder_normal_ihc.Rdata")
save(bladder_cancer_ihc, file = "results_data/bladder_cancer_ihc.Rdata")

#apply function to each pathway
egfr_hpa <- add_hpa_annotations(egfr_metrics, bladder_normal_ihc, bladder_cancer_ihc)
fgfr3_hpa <- add_hpa_annotations(fgfr3_metrics, bladder_normal_ihc, bladder_cancer_ihc)
erbb2_hpa <- add_hpa_annotations(erbb2_metrics, bladder_normal_ihc, bladder_cancer_ihc)

#save
save(egfr_hpa, file = "results_data/egfr_hpa.Rdata")
save(fgfr3_hpa, file = "results_data/fgfr3_hpa.Rdata")
save(erbb2_hpa, file = "results_data/erbb2_hpa.Rdata")

#load data
load("ppi_signal/results_data/egfr_node_metrics.Rdata")
load("ppi_signal/results_data/fgfr3_node_metrics.Rdata")
load("ppi_signal/results_data/erbb2_node_metrics.Rdata")

load("ppi_signal/results_data/bladder_normal_seurat_processed.Rdata")
load("ppi_signal/results_data/bladder_cancer_seurat_processed.Rdata")

#source functions
source("ppi_signal/functions/run_single_cell_classification.R")

#EGFR
egfr_sc_normal = run_single_cell_classification(seurat_object = bladder_normal_processed, my_genes = egfr_node_metrics$name, tissue_type = "normal")
egfr_sc_cancer = run_single_cell_classification(seurat_object = bladder_cancer_processed, my_genes = egfr_node_metrics$name, tissue_type = "cancer")
egfr_single_cell_class = egfr_sc_normal %>% left_join(egfr_sc_cancer, by = "name")

#FGFR
fgfr3_sc_normal = run_single_cell_classification(seurat_object = bladder_normal_processed, my_genes = fgfr3_node_metrics$name, tissue_type = "normal")
fgfr3_sc_cancer = run_single_cell_classification(seurat_object = bladder_cancer_processed, my_genes = fgfr3_node_metrics$name, tissue_type = "cancer")
fgfr3_single_cell_class = fgfr3_sc_normal %>% left_join(fgfr3_sc_cancer, by = "name")

#ERBB2
erbb2_sc_normal = run_single_cell_classification(seurat_object = bladder_normal_processed, my_genes = erbb2_node_metrics$name, tissue_type = "normal")
erbb2_sc_cancer = run_single_cell_classification(seurat_object = bladder_cancer_processed, my_genes = erbb2_node_metrics$name, tissue_type = "cancer")
erbb2_single_cell_class = erbb2_sc_normal %>% left_join(erbb2_sc_cancer, by = "name")

#save results
save(egfr_single_cell_class, file = "ppi_signal/results_data/egfr_single_cell_class.Rdata")
save(fgfr3_single_cell_class, file = "ppi_signal/results_data/fgfr3_single_cell_class.Rdata")
save(erbb2_single_cell_class, file = "ppi_signal/results_data/erbb2_single_cell_class.Rdata")

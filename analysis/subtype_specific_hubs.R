library(ComplexHeatmap)
library(circlize)
library(dplyr)

load('../../GIT_REPOS/ppi_pathways/data/uroscanseq_data_.Rdata')
load('../../GIT_REPOS/ppi_pathways/results_out/erbb2_meta_consolidated.Rdata')
load('../../GIT_REPOS/ppi_pathways/results_out/fgfr3_meta_consolidated.Rdata')
source("../../GIT_REPOS/ppi_pathways/functions/plot_rank_heatmap.R")

ppi_data <- fgfr3_meta_consoldiated

node_data <- read.csv(file = "../../analysis/ppi_presentation/erbb2/erbb_node_data.csv")
node_data <- read.csv(file = "../../analysis/ppi_presentation/erbb2/FGFR3_merged_cummunity default node.csv")

top_genes <- node_data %>% dplyr::select(shared.name)

ppi_data <- ppi_data %>% filter(name %in% top_genes$shared.name)

#add 5 class info also
ppi_data <- ppi_data %>%
  rename(
    rank_7_UroA = rank_UroA,
    rank_7_UroB = rank_UroB,
    rank_7_UroC = rank_UroC,
    rank_7_GU = rank_GU,
    rank_7_BaSq = rank_BaSq
  )

uro_samples <- names(uroscanseq_data$subtype_5_vector)[uroscanseq_data$subtype_5_vector == "Uro"]
gu_samples <- names(uroscanseq_data$subtype_5_vector)[uroscanseq_data$subtype_5_vector == "GU"]
basq_samples <- names(uroscanseq_data$subtype_5_vector)[uroscanseq_data$subtype_5_vector == "BaSq"]

genes_in_ppi <- ppi_data$name

mean_expr_uro <- rowMeans(uroscanseq_data$expr_df[genes_in_ppi, uro_samples], na.rm = TRUE)
mean_expr_gu_5 <- rowMeans(uroscanseq_data$expr_df[genes_in_ppi, gu_samples], na.rm = TRUE)
mean_expr_basq_5 <- rowMeans(uroscanseq_data$expr_df[genes_in_ppi, basq_samples], na.rm = TRUE)

ppi_data <- ppi_data %>%
  mutate(
    mean_expr_5_Uro = mean_expr_uro[name],
    mean_expr_5_GU = mean_expr_gu_5[name],
    mean_expr_5_BaSq = mean_expr_basq_5[name]
  )

ppi_data <- ppi_data %>%
  rowwise() %>%
  mutate(
    expr_5class = list(c(Uro = mean_expr_5_Uro, 
                         GU = mean_expr_5_GU, 
                         BaSq = mean_expr_5_BaSq)),
    ranks_5class = list(rank(-unlist(expr_5class), ties.method = "first")),
    rank_5_Uro = ranks_5class[[1]],
    rank_5_GU = ranks_5class[[2]],
    rank_5_BaSq = ranks_5class[[3]]
  ) %>%
  ungroup() %>%
  select(-expr_5class, -ranks_5class)

#plots!
plot_rank_heatmap(data = ppi_data, system = "5", order_by = "hub_score", top_n = 10)
plot_rank_heatmap(data = ppi_data, system = "7", order_by = "hub_score", top_n = 10)

calculate_hub_scores(data = ppi_data, system = "5", order_by = "hub_score", top_n = 10)
calculate_hub_scores(data = ppi_data, system = "7", order_by = "hub_score", top_n = 10)

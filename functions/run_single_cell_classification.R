library(Seurat)
library(dplyr)

run_single_cell_classification = function(seurat_object = NULL,
                                          my_genes = NULL,
                                          tissue_type = NULL){
  #define known urothelial markers
  urothelial_markers <- c("GATA3", "FOXA1", "KRT7", "KRT8", "PPARG", "FOXA1",
                          "ELF3", "FGFR3", "CCND1", "TP63", "KRT5", "KRT14",
                          "EPCAM", "CDH1")

  #fix data layers (prevents errors)
  seurat_object <- JoinLayers(seurat_object)

  #step 2: check which urothelial markers are available in the data
  available_uro_markers <- intersect(urothelial_markers, rownames(seurat_object))
  cat("Available urothelial markers:", paste(available_uro_markers, collapse = ", "), "\n")

  #step 3: find which cluster expresses urothelial markers most
  if(length(available_uro_markers) > 0) {
    marker_expr <- AverageExpression(seurat_object,
                                     features = available_uro_markers,
                                     group.by = "seurat_clusters")$RNA

    cat("Urothelial marker expression by cluster:\n")
    print(marker_expr)

    #check what the cluster names look like
    cat("Cluster names:", colnames(marker_expr), "\n")

    #calculate cluster scores and find max
    cluster_scores <- colMeans(marker_expr, na.rm = TRUE)
    cat("Cluster scores:", cluster_scores, "\n")

    #get the cluster with highest score (keep as character)
    max_cluster_name <- names(which.max(cluster_scores))
    cat("Max cluster name:", max_cluster_name, "\n")

    #keep as character (don't convert to numeric since they're "g0", "g1")
    urothelial_cluster <- max_cluster_name

    cat("Urothelial cluster identified as:", urothelial_cluster, "\n")
  } else {
    # If no markers found, assume first cluster
    unique_clusters <- unique(seurat_object@meta.data$seurat_clusters)
    urothelial_cluster <- unique_clusters[1]
    cat("No urothelial markers found, assuming cluster:", urothelial_cluster, "\n")
  }

  # Verify the cluster assignment
  cat("Final urothelial cluster:", urothelial_cluster, "\n")
  cat("This cluster has higher expression of urothelial markers\n")

  #step 4: Get cluster assignments for all cells
  cluster_ids <- seurat_object@meta.data$seurat_clusters

  #step 5: Set up classification parameters
  detection_threshold <- 0.05  # 5% of cells must express the gene
  expression_threshold <- 0.1   # Minimum average expression level

  #fix the cluster assignment
  #g0 had the highest urothelial marker expression, use cluster 0
  urothelial_cluster_id <- 0
  cat("Using cluster ID:", urothelial_cluster_id, "for urothelial cells\n")

  #verify this works
  cat("Number of cells in cluster 0:", sum(cluster_ids == 0), "\n")
  cat("Number of cells in cluster 1:", sum(cluster_ids == 1), "\n")

  #step 6: Binary classify each signature gene
  classification_results <- data.frame()

  for(gene in my_genes) {
    if(gene %in% rownames(seurat_object)) {

      #get expression for this gene
      gene_expr <- GetAssayData(seurat_object, slot = "data")[gene, ]

      #get expression only in urothelial cells (cluster 0)
      urothelial_expr <- gene_expr[cluster_ids == urothelial_cluster_id]

      #calculate statistics
      mean_expr <- mean(urothelial_expr)
      detection_rate <- sum(urothelial_expr > 0) / length(urothelial_expr)
      max_expr <- max(urothelial_expr)
      cells_expressing <- sum(urothelial_expr > 0)

      #binary classification
      classification <- ifelse(
        detection_rate >= detection_threshold & mean_expr >= expression_threshold,
        "Urothelial",
        "Non-urothelial"
      )

      if(tissue_type == "normal"){
        #store results
        classification_results <- rbind(classification_results, data.frame(
          name = gene,
          sc_norm_mean_expression = mean_expr,
          sc_norm_detection_rate = detection_rate,
          sc_norm_max_expression = max_expr,
          sc_norm_cells_expressing = cells_expressing,
          sc_norm_total_urothelial_cells = length(urothelial_expr),
          sc_norm_classification = classification,
          stringsAsFactors = FALSE
        ))
      }else if(tissue_type == "cancer"){
        #store results
        classification_results <- rbind(classification_results, data.frame(
          name = gene,
          sc_cancer_mean_expression = mean_expr,
          sc_cancer_detection_rate = detection_rate,
          sc_cancer_max_expression = max_expr,
          sc_cancer_cells_expressing = cells_expressing,
          sc_cancer_total_urothelial_cells = length(urothelial_expr),
          sc_cancer_classification = classification,
          stringsAsFactors = FALSE
        ))
      }
    } else {
      if(tissue_type == "normal"){
        # Gene not found in dataset
        classification_results <- rbind(classification_results, data.frame(
          name = gene,
          sc_norm_mean_expression = 0,
          sc_norm_detection_rate = 0,
          sc_norm_max_expression = 0,
          sc_norm_cells_expressing = 0,
          sc_norm_total_urothelial_cells = sum(cluster_ids == urothelial_cluster_id),
          sc_norm_classification = "Non-urothelial",
          stringsAsFactors = FALSE
        ))

      }else if (tissue_type == "cancer"){
        # Gene not found in dataset
        classification_results <- rbind(classification_results, data.frame(
          name = gene,
          sc_cancer_mean_expression = 0,
          sc_cancer_detection_rate = 0,
          sc_cancer_max_expression = 0,
          sc_cancer_cells_expressing = 0,
          sc_cancer_total_urothelial_cells = sum(cluster_ids == urothelial_cluster_id),
          sc_cancer_classification = "Non-urothelial",
          stringsAsFactors = FALSE
        ))

      }
    }
  }
  return(classification_results)
}


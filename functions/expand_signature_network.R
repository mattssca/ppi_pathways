expand_signature_network <- function(expr_data,
                                     subtype_vector,
                                     seed_genes = NULL,
                                     max_added_genes = 20,
                                     min_degree = 1,
                                     string_score_threshold = 400,
                                     out_dir = "viz/signature_networks", 
                                     verbose = TRUE,
                                     return_data = FALSE,
                                     layout_method = "kk",
                                     show_labels = TRUE,
                                     expr_summary = "mean",
                                     node_color = "lund",
                                     theme = "light",
                                     node_size = 10,
                                     node_degree = FALSE,
                                     plot_width = 15,
                                     plot_height = 10,
                                     color_scale_global = TRUE,
                                     genes_blacklist = NULL) {
  
  #suppress unused variable warnings
  all_ids <- NULL
  
  #create directory, if non-existing
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Get available genes in expression data
  available_genes <- rownames(expr_data)
  
  #get signature genes and filter to those in expression data
  signature_genes <- intersect(seed_genes, available_genes)
  
  if(length(signature_genes) == 0) {
    stop("None of the seed genes are present in the expression data")
  }
  
  if(verbose && length(signature_genes) < length(seed_genes)) {
    message(sprintf("Warning: %d seed genes not found in expression data", 
                    length(seed_genes) - length(signature_genes)))
  }
  
  ########### Step 1 ############
  if(verbose){message("1/8 Running StringDB...")}
  
  #run stringDB for retreiving PPI
  string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = string_score_threshold, input_directory = "")
  
  
  ########### Step 2 ############
  if(verbose){message("2/8 Mapping genes...")}
  
  #get STRING IDs for the genes in the selected signature
  mapped <- string_db$map(data.frame(gene = signature_genes), "gene", removeUnmappedRows = TRUE)
  signature_ids <- mapped$STRING_id
  
  
  ########### Step 3 ############
  if(verbose){message("3/8 Finding Neighbours...")}
  
  #get neighbors for the signature genes
  neighbors <- tryCatch(string_db$get_neighbors(signature_ids), error = function(e) character(0))
  if(length(neighbors) == 0 || all(is.na(neighbors))) stop("No neighbors found for signature genes in STRING.")
  
  #ensure neighbors are character and not NA/empty
  neighbors <- as.character(neighbors)
  neighbors <- neighbors[!is.na(neighbors) & neighbors != ""]
  if(length(neighbors) == 0) stop("No valid STRING neighbor IDs for alias mapping.")
  
  
  ########### Step 4 ############
  if(verbose){message("4/8 Mapping Protein IDs to Gene Symbols...")}
  
  #use biomaRt to map STRING protein IDs to gene symbols
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  neighbor_peptides <- sub("9606\\.", "", neighbors)
  
  mapping <- getBM(attributes = c("ensembl_peptide_id", "external_gene_name"),
                   filters = "ensembl_peptide_id",
                   values = neighbor_peptides,
                   mart = mart)
  
  #get unique neighbor genes that are not in the signature object
  neighbor_genes <- unique(mapping$external_gene_name)
  neighbor_genes <- setdiff(neighbor_genes, signature_genes)
  
  # FILTER: Keep only neighbors present in expression data
  neighbor_genes <- intersect(neighbor_genes, available_genes)
  
  if(verbose) {
    message(sprintf("Found %d neighbor genes present in expression data", length(neighbor_genes)))
  }
  
  
  ########### Step 5 ############
  if(verbose){message("5/8 Filtering Edges...")}
  
  #build initial network (signature genes)
  #filter by degree (number of connections to signature)
  edges <- suppressWarnings(string_db$get_interactions(signature_ids))
  
  if(nrow(edges) == 0) stop("No interactions found for signature genes in STRING.")
  
  #count the number of connections (degree) for each neighbor gene, 
  neighbor_counts <- table(c(edges$from, edges$to))
  
  #filter to include only valid neighbor genes, select the top neighbors 
  neighbor_counts <- neighbor_counts[names(neighbor_counts) %in% neighbor_genes]
  
  #exclude genes in the blacklist
  if(!is.null(genes_blacklist)){
    neighbor_counts <- neighbor_counts[!names(neighbor_counts) %in% genes_blacklist]
  }
  
  #based on degree (up to max_added_genes), and combine them with the 
  top_neighbors <- names(sort(neighbor_counts, decreasing=TRUE))[1:min(max_added_genes, length(neighbor_counts))]
  top_neighbors <- top_neighbors[!is.na(top_neighbors)]
  
  #original signature genes to create the expanded gene set.
  expanded_genes <- unique(c(signature_genes, top_neighbors))
  
  # FILTER: Ensure all expanded genes are in expression data
  expanded_genes <- intersect(expanded_genes, available_genes)
  
  #combine signature genes and neighbor genes, map them to STRING IDs, 
  all_genes <- unique(c(signature_genes, neighbor_genes))
  
  # FILTER: Keep only genes in expression data
  all_genes <- intersect(all_genes, available_genes)
  
  all_mapped <- suppressMessages(suppressWarnings(string_db$map(data.frame(gene = all_genes), "gene", removeUnmappedRows = TRUE)))
  all_ids <- all_mapped$STRING_id
  if(length(all_ids) == 0 || any(is.na(all_ids))) stop("No valid STRING IDs for expanded gene set.")
  
  
  ########### Step 6 ############
  if(verbose){message("6/8 Expanding Network...")}
  
  #get all interactions among expanded set
  network_edges <- string_db$get_interactions(all_ids)
  
  if(nrow(network_edges) == 0) stop("No network edges found for expanded gene set.")
  
  #count degrees for all nodes in the expanded network
  degree_table <- table(c(network_edges$from, network_edges$to))
  
  #map STRING IDs to gene symbols using biomaRt
  node_peptides <- sub("9606\\.", "", names(degree_table))
  
  node_mapping <- getBM(attributes = c("ensembl_peptide_id", "external_gene_name"),
                        filters = "ensembl_peptide_id",
                        values = node_peptides,
                        mart = mart)
  
  peptide_to_symbol <- setNames(node_mapping$external_gene_name, node_mapping$ensembl_peptide_id)
  
  degree_symbols <- peptide_to_symbol[node_peptides]
  
  #select top neighbors by degree, excluding signature genes
  top_neighbors <- setdiff(degree_symbols[order(degree_table, decreasing = TRUE)], signature_genes)
  top_neighbors <- unique(top_neighbors[!is.na(top_neighbors)])
  
  # FILTER: Keep only neighbors in expression data
  top_neighbors <- intersect(top_neighbors, available_genes)
  
  top_neighbors <- head(top_neighbors, max_added_genes)
  expanded_genes <- unique(c(signature_genes, top_neighbors))
  
  # FILTER: Final check - all expanded genes must be in expression data
  expanded_genes <- intersect(expanded_genes, available_genes)
  
  #build network for expanded genes
  expanded_mapped <- string_db$map(data.frame(gene = expanded_genes), "gene", removeUnmappedRows = TRUE)
  expanded_ids <- expanded_mapped$STRING_id
  network_edges <- string_db$get_interactions(expanded_ids)
  
  if(nrow(network_edges) == 0) stop("No network edges found for expanded gene set.")
  
  g <- graph_from_data_frame(network_edges, directed = FALSE)
  
  #map node names to gene symbols
  node_peptides <- sub("9606\\.", "", V(g)$name)
  
  node_mapping <- getBM(attributes = c("ensembl_peptide_id", "external_gene_name"),
                        filters = "ensembl_peptide_id",
                        values = node_peptides,
                        mart = mart)
  
  peptide_to_symbol <- setNames(node_mapping$external_gene_name, node_mapping$ensembl_peptide_id)
  
  V(g)$name <- peptide_to_symbol[node_peptides]
  
  # FILTER: Remove any nodes not in expression data
  nodes_to_keep <- V(g)$name %in% available_genes
  g <- induced_subgraph(g, which(nodes_to_keep))
  
  if(verbose){
    cat("Signature genes:", length(signature_genes), "\n")
    cat("Neighbors found:", length(neighbors), "\n")
    cat("Mapped neighbor genes:", length(neighbor_genes), "\n")
    cat("Expanded gene set:", length(expanded_genes), "\n")
    cat("Network nodes:", length(V(g)), "\n")
    cat("All nodes in expression data:", all(V(g)$name %in% available_genes), "\n")
    print(V(g)$name)
  }
  
  if (return_data) {
    message("Compiling data for return...")
    
    #compute additional node metrics
    n_nodes <- length(V(g))
    node_metrics <- data.frame(
      name = V(g)$name,
      degree = if(length(igraph::degree(g)) == n_nodes) igraph::degree(g) else rep(NA, n_nodes),
      is_seed = V(g)$name %in% signature_genes,
      betweenness = if(length(betweenness(g)) == n_nodes) betweenness(g) else rep(NA, n_nodes),
      closeness = if(length(closeness(g)) == n_nodes) closeness(g) else rep(NA, n_nodes),
      eigenvector = if(length(eigen_centrality(g)$vector) == n_nodes) eigen_centrality(g)$vector else rep(NA, n_nodes),
      hub_score = if(length(hub_score(g)$vector) == n_nodes) hub_score(g)$vector else rep(NA, n_nodes),
      authority_score = if(length(authority_score(g)$vector) == n_nodes) authority_score(g)$vector else rep(NA, n_nodes),
      community = if(length(membership(cluster_louvain(g))) == n_nodes) as.numeric(membership(cluster_louvain(g))) else rep(NA, n_nodes)
    )
    
    # Add mean expression for each subtype as new columns
    subtypes <- unique(subtype_vector)
    for (subtype in subtypes) {
      samples <- names(subtype_vector)[subtype_vector == subtype]
      valid_samples <- intersect(samples, colnames(expr_data))
      
      # Ensure we only extract genes that exist in expr_data
      genes_in_data <- intersect(node_metrics$name, rownames(expr_data))
      expr_sub <- expr_data[genes_in_data, valid_samples, drop = FALSE]
      
      # Create a named vector for safe matching
      mean_expr <- rowMeans(expr_sub, na.rm = TRUE)
      node_metrics[[paste0("mean_expr_", subtype)]] <- mean_expr[match(node_metrics$name, names(mean_expr))]
    }
    
    return(list(
      signature_genes = signature_genes,
      neighbors = neighbors,
      neighbor_genes = neighbor_genes,
      expanded_genes = expanded_genes,
      g = g,
      node_metrics = node_metrics
    ))
  } else {
    invisible(NULL)
  }
}

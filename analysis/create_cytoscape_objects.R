

#get edge data
egfr_edges <- igraph::as_data_frame(egfr_nodes$g, what = "edges")
fgfr3_edges <- igraph::as_data_frame(fgfr3_nodes$g, what = "edges")
erbb2_edges <- igraph::as_data_frame(erbb2_nodes$g, what = "edges")

#write edge data
write.table(
  egfr_edges,
  file = "cytoscape/egfr_edges.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  fgfr3_edges,
  file = "cytoscape/fgfr3_edges.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  erbb2_edges,
  file = "cytoscape/erbb2_edges.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#write node data
write.table(
  egfr_meta_consoldiated,
  file = "cytoscape/egfr_nodes.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  fgfr3_meta_consoldiated,
  file = "cytoscape/fgfr3_nodes.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  erbb2_meta_consoldiated,
  file = "cytoscape/erbb2_nodes.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

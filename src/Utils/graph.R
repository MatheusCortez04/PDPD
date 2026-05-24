library(igraph)
library(here)
source(here("src","Utils","utils.R"))


generate_graph_from_dataframe = function(dataframe){
cat("Creating graph from dataframe:\n")
graph = igraph::graph_from_data_frame(dataframe, directed = FALSE) %>%
  igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE)
print("Graph created successfully!")
return(graph)
}


graph_menu = function(){
  while(TRUE){
    clear_console()
    cat("--- Graph Menu ---\n\n")
    cat(" [1] Plot MDD Disease Module with LCC \n")
    cat(" [2] Plot BD Disease Module with LCC \n")
    cat(" [B] Back\n\n")

    input = readline(prompt = "Choice option: ")
    if (toupper(trimws(input)) == "B") break
     graph_function_mapper[[input]]()
  }
}

graph_function_mapper = list(
    '1' = function() {
        plot_disease_module("MDD")
    },
    '2' = function(){
         plot_disease_module("BD")
    }
)



extract_lcc_vertices= function(graph){
  total_vertices = igraph::vcount(graph)
  cat(sprintf("[INFO] Computing LCC. Initial graph contains %d vertices...\n", total_vertices))
  components = igraph::components(graph, mode="weak")
  biggest_component_id = which.max(components$csize)
  total_vertices_in_lcc = components$csize[biggest_component_id]

  percent_total_vertices_in_lcc = (total_vertices_in_lcc / total_vertices) * 100
  cat(sprintf("[INFO] LCC extraction done: N = %d (%.1f%% of input graph)\n", total_vertices_in_lcc, percent_total_vertices_in_lcc))
  vert_ids = igraph::V(graph)[components$membership == biggest_component_id]
  return(names(vert_ids))
}


get_disease_subgraph  = function(disease=c('MDD','BD')){
    disease = match.arg(disease)
    cat(sprintf("[INFO] Extracting disease subgraph for: %s\n", disease))
    disease_data <- list(
        MDD = list(name = "MDD", genes = get_mdd_disease_module()$entrez_id),
        BD = list(name = "BD", genes = get_bipolar_disease_module()$entrez_id)
    )
  ppi_df = import_ppi_interactions()
  global_graph = generate_graph_from_dataframe(ppi_df)
  lcc_global = extract_lcc_vertices(global_graph)
    
  disease_genes_df = disease_data[[disease]]$genes
  disease_genes_in_lcc = intersect(lcc_global, disease_genes_df)
  
  disease_subgraph = igraph::induced_subgraph(global_graph, vids = disease_genes_in_lcc)
  disease_subgraph = igraph::simplify(disease_subgraph, remove.multiple = TRUE, remove.loops = TRUE)
  return(disease_subgraph)
}

plot_disease_module = function(disease=c('MDD','BD')){
  set.seed(10)
  disease = match.arg(disease)
  graph_name= c("MDD" = "Disease Module Analysis - Major Depressive Disorder Subgraph",
                "BD" ="Disease Module Analysis - Bipolar Disorder Subgraph")
  
  edge_width = 0.7     
  color_edges = "#999999"
  color_common_genes  = "#0072B2"
  color_specific_genes = "#D55E00"
  vertex_size=4

  disease_subgraph = get_disease_subgraph(disease)
  disease_subgraph_vertex_names = igraph::V(disease_subgraph)$name
  total_vertices_in_disease_subgraph = igraph::vcount(disease_subgraph)
  disease_subgraph_lcc_vertices = extract_lcc_vertices(disease_subgraph)

  common_genes = get_common_gene()
  common_mask = disease_subgraph_vertex_names %in% common_genes
  total_common_genes_in_disease_subgraph = sum(common_mask)
  common_genes_in_disease_subgraph_percent = round(100 * total_common_genes_in_disease_subgraph / total_vertices_in_disease_subgraph, 1)
  common_genes_in_lcc_disease_subgraph = sum(disease_subgraph_lcc_vertices %in% common_genes)
  common_genes_in_lcc_disease_subgraph_percent = round(100 * common_genes_in_lcc_disease_subgraph / length(disease_subgraph_lcc_vertices), 1)

  total_disease_specific_genes_in_subgraph = (total_vertices_in_disease_subgraph-total_common_genes_in_disease_subgraph)
  specific_genes_in_subgraph_percent =round(100 * total_disease_specific_genes_in_subgraph / total_vertices_in_disease_subgraph, 1)

  igraph::V(disease_subgraph)$color = color_specific_genes
  igraph::V(disease_subgraph)$color[common_mask] = color_common_genes
  igraph::V(disease_subgraph)$size = vertex_size

  graph_layout = igraph::layout_with_fr(disease_subgraph)

  output_graph_dir= here('src','Graph',disease)
  create_dir(output_graph_dir)
  pdf_path = here(output_graph_dir,sprintf("Largest_Connected_Component_%s.pdf", disease))
  cat(sprintf("\n[INFO] Saving disease module plot to:\n%s\n", pdf_path))
  grDevices::pdf(file =pdf_path, width = 8, height = 8)

  plot(
    disease_subgraph,
    vertex.size = igraph::V(disease_subgraph)$size,
    vertex.label = NA,
    vertex.color = igraph::V(disease_subgraph)$color,
    edge.color = grDevices::adjustcolor(color_edges, alpha.f = 0.8),
    edge.width = edge_width,
    main =graph_name[disease],
    layout = graph_layout
  )
  graphics::legend(
  "topleft",
  inset = c(-0.07, -0.05),
  xpd = TRUE,
  legend = c(
      sprintf("Disease-Specific Genes (%d | %.1f%%)", total_disease_specific_genes_in_subgraph, specific_genes_in_subgraph_percent),
      sprintf("Common Genes (%d | %.1f%%)", total_common_genes_in_disease_subgraph, common_genes_in_disease_subgraph_percent)),
  col = c(color_specific_genes, color_common_genes),
  pch = 19,
  pt.cex = 1.5,
  bty = "n",
  cex = 0.9)
  score_filter = get_score_disease_gene_association()

  graphics::mtext(
    side = 1,
    line = 2, 
    text = sprintf("Total disease-associated genes (score filter %s): %d", score_filter, total_vertices_in_disease_subgraph),
    cex = 1.0,
    font = 1,
    col = "black")

  graphics::mtext(
      side = 1,
      line = 3.5,
      text = sprintf("Disease Module LCC size: %d | Common genes in LCC: %d (%.1f%%)", length(disease_subgraph_lcc_vertices),common_genes_in_lcc_disease_subgraph,common_genes_in_lcc_disease_subgraph_percent),
      cex = 1.0,
      font = 2, 
      col = "black"
    )
  grDevices::dev.off()
  invisible(disease_subgraph)
}


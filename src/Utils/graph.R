library(igraph)
library(here)
source(here("src","Utils","utils.R"))


generate_graph_from_dataframe = function(dataframe){
cat("Creating graph from dataframe:\n")
graph <- igraph::graph_from_data_frame(dataframe, directed = FALSE)
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
        plot_lcc("MDD")
    },
    '2' = function(){
         plot_lcc("BD")
    }
)



get_lcc_vertices= function(graph){
  total_nodes = igraph::vcount(graph)
  cat(sprintf("[INFO] Computing LCC. Initial graph contains %d nodes...\n", total_nodes))
  components = igraph::components(graph, mode="weak")
  biggest_cluster_id = which.max(components$csize)
  lcc_count = components$csize[biggest_cluster_id]

  percent = (lcc_count / total_nodes) * 100
  cat(sprintf("[INFO] LCC extraction done: N = %d (%.1f%% of input graph)\n", lcc_count, percent))
  vert_ids = igraph::V(graph)[components$membership == biggest_cluster_id]
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
  lcc_global = get_lcc_vertices(global_graph)
    
  disease_genes_df = disease_data[[disease]]$genes
  disease_genes_in_lcc = intersect(lcc_global, disease_genes_df)
  
  disease_subgraph = igraph::induced_subgraph(global_graph, vids = disease_genes_in_lcc)
  disease_subgraph = igraph::simplify(disease_subgraph, remove.multiple = TRUE, remove.loops = TRUE)
  return(disease_subgraph)
}

plot_lcc = function(disease=c('MDD','BD')){
  set.seed(10)
  disease = match.arg(disease)
  graph_name= c("MDD" = "Disease Module Analysis - Major Depressive Disorder Subgraph",
                "BD" ="Disease Module Analysis - Bipolar Disorder Subgraph")
  
  output_graph_dir= here('src','Graph',disease)
  disease_subgraph = get_disease_subgraph(disease)
  subgraph_lcc = get_lcc_vertices(disease_subgraph)

  total_nodes = igraph::vcount(disease_subgraph)
  common_genes = get_common_gene()

  common_in_network = sum(igraph::V(disease_subgraph)$name %in% common_genes)
  common_percent = round(100 * common_in_network / total_nodes, 1)
  common_in_lcc = sum(subgraph_lcc %in% common_genes)

  lcc_exclusive_count = length(subgraph_lcc) - common_in_lcc
  lcc_percent = round(100 * lcc_exclusive_count / total_nodes, 1)

  non_lcc_count = total_nodes - (common_in_network + lcc_exclusive_count)
  non_lcc_percent = round(100 - (common_percent + lcc_percent), 1)

  igraph::V(disease_subgraph)$color = "gray70"
  igraph::V(disease_subgraph)$size = 3

  lcc_mask = igraph::V(disease_subgraph)$name %in% subgraph_lcc
  igraph::V(disease_subgraph)$color[lcc_mask] = 'tomato'
  igraph::V(disease_subgraph)$size[lcc_mask] = 6

  common_mask = igraph::V(disease_subgraph)$name %in% common_genes
  igraph::V(disease_subgraph)$color[common_mask] = 'royalblue'
  igraph::V(disease_subgraph)$size[common_mask] = 5




  layout = igraph::layout_with_fr(disease_subgraph)
  create_dir(output_graph_dir)
  pdf_path = here(output_graph_dir,sprintf("Largest_Connected_Component_%s.pdf", disease))
  cat(sprintf("[INFO] Saving network plot directly to: %s\n", pdf_path))
  grDevices::pdf(file =pdf_path, width = 8, height = 8)
  
  plot(
    disease_subgraph,
    vertex.size = igraph::V(disease_subgraph)$size,
    vertex.label = NA,
    vertex.color = igraph::V(disease_subgraph)$color,
    edge.color = grDevices::adjustcolor("gray70", alpha.f = 0.4),
    edge.width = 0.5,
    main =graph_name[disease],
    layout = layout
  )
  graphics::legend(
  "topleft",
  legend = c(
      sprintf("Common Genes (%d | %.1f%%)", common_in_network, common_percent),
      sprintf("LCC Exclusive Genes (%d | %.1f%%)", lcc_exclusive_count, lcc_percent),
      sprintf("Other Genes (%d | %.1f%%)", non_lcc_count, non_lcc_percent)
  ),
  col = c('royalblue', 'tomato', "gray70"),
  pch = 19,
  pt.cex = 1.5,
  bty = "n",
  cex = 0.9)
  score_filter = get_score_disease_gene_association()

graphics::mtext(
    side = 1,
    line = 2, 
    text = sprintf("Total genes associated with the disease (score filter %s): %d", score_filter, total_nodes),
    cex = 1.0,
    font = 1,
    col = "black")

  graphics::mtext(
    side = 1,
    line = 3.5,
    text = sprintf("Total genes in the Largest Connected Component (LCC): %d - Common genes: %d", length(subgraph_lcc),common_in_lcc),
    cex = 1.0,
    font = 2, 
    col = "black"
  )

  grDevices::dev.off()
  invisible(disease_subgraph)
}


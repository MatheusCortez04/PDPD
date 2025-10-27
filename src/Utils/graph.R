library(igraph)

generate_graph_from_dataframe = function(dataframe){
cat("Creating graph from dataframe:\n")
graph <- graph_from_data_frame(dataframe, directed = FALSE)
print("Graph created successfully!")
return(graph)
}

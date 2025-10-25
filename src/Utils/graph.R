library(igraph)

generate_graph_from_dataframe = function(dataframe){
print("Criando Grafo com o dataframe:")
print(summary(dataframe))
graph <- graph_from_data_frame(dataframe, directed = FALSE)
print(" Grafo criado com sucesso!")
return(graph)
}

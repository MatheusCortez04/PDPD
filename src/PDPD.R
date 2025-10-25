library(utils)
library(here)
library(dplyr)

source(here("src","Utils","graph.R"))


main= function(){
    cat("Bibliotecas carregadas com sucesso.")
    cat("Realizando Leitura da PPI\n")
    ppi_df  = read.csv(here("src","Data","PPI_gysi.csv"), sep=",")
    cat("Leitura da PPI Finalizada\n")
    cat("Realizando criacao do dataframe somente com as proteinas\n")
    ppi_only_proteins_df =  ppi_df %>% select(proteinA_entrezid,proteinB_entrezid)
    graph = generate_graph_from_dataframe(ppi_only_proteins_df)
}


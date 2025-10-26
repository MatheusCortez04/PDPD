library(utils)
library(here)
library(dplyr)

source(here("src","Utils","graph.R"))

source(here("src","Utils","kernels.R"))


main= function(){
    cat("Libraries loaded successfully.")
    cat("Reading PPI data...\n")
    ppi_df  = read.csv(here("src","Data","PPI_gysi.csv"), sep=",")
    cat("PPI reading complete.\n")
    cat("Realizando criacao do dataframe somente com as proteinas\n")
    ppi_only_proteins_df =  ppi_df %>% select(proteinA_entrezid,proteinB_entrezid)
    graph = generate_graph_from_dataframe(ppi_only_proteins_df)
    difussion_kernel = generate_difussion_kernel(graph)
    pstep_kernel = generate_pstep_kernel(graph)
}


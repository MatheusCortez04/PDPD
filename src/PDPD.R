library(here)
library(igraph)
library(diffuStats)
library(dplyr)

print("Bibliotecas carregadas com sucesso.")
print("Realizando Leitura da PPI")

ppi_df  = read.csv(here("Data","PPI_gysi.csv"), sep=",")

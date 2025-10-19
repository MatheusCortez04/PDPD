library(here)
library(igraph)
library(diffuStats)
library(dplyr)

print("Bibliotecas carregadas com sucesso.")
print("Realizando Leitura da PPI")

ppi_df  = read.csv(here("Data","PPI_gysi.csv"), sep=",")
ppi_only_proteins_df = select(ppi_df,proteinA_entrezid,proteinB_entrezid)
write.csv(ppi_only_proteins_df,"Data/PPI_ONLY_PROTEINS.CSV")
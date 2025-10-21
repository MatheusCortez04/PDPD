library(here)
library(igraph)
library(diffuStats)
library(dplyr)
library(Matrix)

cat("Bibliotecas carregadas com sucesso.")
cat("Realizando Leitura da PPI\n")

ppi_df  = read.csv(here("Data","PPI_gysi.csv"), sep=",")
cat("Leitura da PPI Finalizada\n")
cat("Realizando criacao do dataframe somente com as proteinas\n")

ppi_only_proteins_df =  ppi_df %>% select(proteinA_entrezid,proteinB_entrezid)
print("Dataframe somente com as proteinas salvo com sucesso\n")

print("Criando Grafo com o dataframe contendo somente as proteinas\n")
graph <- graph_from_data_frame(ppi_only_proteins_df, directed = FALSE)

cat("Calculando o kernel diffusionKernel para o seu grafo...\n")
 K_diffusionKernel = diffusionKernel(graph,normalized = TRUE)

cat("Kernel calculado com sucesso!\n\n")
write.csv(K_diffusionKernel, file = "KERNEL_diffusionKernel.csv")
cat("--- Iniciando conversão para Matriz Esparsa  ---\n")
zeros_exatos <- sum(K_diffusionKernel == 0)
cat(paste("... Encontrados", zeros_exatos, "zeros exatos na matriz densa.\n"))

cat("... Removendo zeros exatos da estrutura (drop0)\n")
kernel_esparso <- drop0(K_diffusionKernel)
saveRDS(K_diffusionKernel, file = "KERNEL_diffusionKernel_ESPARSO.rds")
cat("Arquivo KERNEL_diffusionKernel_ESPARSO.rds salvo com sucesso.\n\n")

cat("--- Comparação de Tamanho em Memória ---\n")
print(paste("Tamanho do Kernel DENSO:", format(object.size(K_diffusionKernel), units = "auto")))
print(paste("Tamanho do Kernel ESPARSO:", format(object.size(kernel_esparso), units = "auto")))
cat("----------------------------------------\n\n")
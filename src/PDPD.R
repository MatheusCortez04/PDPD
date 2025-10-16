library(here)
library(igraph)
library(diffuStats)
library(Matrix)

print("Bibliotecas carregadas com sucesso.")
print("Realizando Leitura da PPI")

PPI_DF  = read.csv(here( "Data", "PPI_gysi.csv"), sep=",")
summary(PPI_DF)
print("Leitura da PPI finalizada")
grafo <- graph_from_data_frame(PPI_DF, directed = FALSE)

V(grafo)$name <- as.character(V(grafo)$name)
cat("Grafo criado com sucesso!\n")
cat("\n")

cat("Calculando o kernel de p-passos para o seu grafo...\n")
K = pStepKernel(graph = grafo)
cat("Kernel calculado com sucesso!\n\n")
write.csv(K, file = "kernel_pStep.csv")
K_esparso <- as(K, "sparseMatrix")

cat("--- Comparação de Tamanho em Memória ---\n")
print(paste("Tamanho do kernel denso:", format(object.size(K), units = "auto")))
print(paste("Tamanho do kernel esparso:", format(object.size(K_esparso), units = "auto")))
cat("\n")


cat("Arquivo 'kernel_pStep_esparso.rds' salvo com sucesso!\n")
cat("Este arquivo é muito menor e mais rápido de carregar do que o CSV.\n\n")

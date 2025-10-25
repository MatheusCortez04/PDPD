library(Matrix)


removed_0_in_matrix = function(kernel){
    cat("--- Iniciando conversão para Matriz Esparsa  ---\n")
    zeros_exatos <- sum(kernel == 0)
    cat(paste("... Encontrados", zeros_exatos, "zeros exatos na matriz densa.\n"))
    cat("... Removendo zeros exatos da estrutura (drop0)\n")
    kernel_esparso <- drop0(kernel)
    cat("--- Comparação de Tamanho em Memória ---\n")
    print(paste("Tamanho do Kernel DENSO:", format(object.size(kernel), units = "auto")))
    print(paste("Tamanho do Kernel ESPARSO:", format(object.size(kernel_esparso), units = "auto")))
    invisible(kernel_esparso)
}
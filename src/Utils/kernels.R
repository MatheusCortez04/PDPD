library(diffuStats)
source(here("src","Utils","Matrix.R"))
generate_difussion_kernel= function(graph,normalized=FALSE,output_path,save_rdata=FALSE){
    path_csv <- paste0(output_path, ".csv")
    path_rdata <- paste0(output_path, ".Rdata")
    cat("Calculando o kernel de difussao para o  grafo...\n")

    diffusionKernel = diffusionKernel(graph,normalized)
    cat("Kernel de difussão calculado com sucesso!\n\n")
    write.csv(diffusionKernel, file =path_csv)
    cat("Kernel de difussão salvo no arquivo:", path_csv, "\n")

    if(save_rdata){
       kernel_esparso = removed_0_in_matrix(diffusionKernel)
        save(kernel_esparso, file = path_rdata)
        cat("Objeto R (kernel_result) salvo em:", path_rdata, "\n")
    }
   invisible(diffusionKernel)
}
generate_p_kernel = function(graph,step=5,output_path,save_rdata=FALSE){
    path_csv <- paste0(output_path, ".csv")
    path_rdata <- paste0(output_path, ".Rdata")
    
    cat("Calculando o kernel de passeio aleatorio para o  grafo...\n")

    diffusionKernel = pStepKernel(graph,p=step)
    cat("Kernel  passeio aleatorio calculado com sucesso!\n\n")
    write.csv(diffusionKernel, file =path_csv)
    cat("Kernel  passeio aleatorio salvo no arquivo:", path_csv, "\n")

    if(save_rdata){
       kernel_esparso = removed_0_in_matrix(diffusionKernel)
        save(kernel_esparso, file = path_rdata)
        cat("Objeto R (kernel_result) salvo em:", path_rdata, "\n")
    }
    invisible(diffusionKernel)

}

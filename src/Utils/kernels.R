library(diffuStats)
generate_difussion_kernel= function(graph,normalized=TRUE,save_rdata=TRUE){
    output_path="diffusion_kernel"
    path_csv = here("src","Data","Kernels",paste0(output_path, ".csv"))
    path_rdata = here("src","Data","Kernels","RData",paste0(output_path, ".Rdata"))
    
    cat("Calculating diffusion kernel for the graph...\n")
    diffusion_kernel = diffusionKernel(graph,normalized,)

    cat("Diffusion kernel calculated successfully!\n\n")
    write.csv(diffusion_kernel, file =path_csv)
    cat("Diffusion kernel saved to file:", path_csv, "\n")

    if(save_rdata){
        save(diffusion_kernel, file = path_rdata)
        cat("R object 'diffusion_kernel' saved to:", path_rdata, "\n")
    }
   invisible(diffusion_kernel)
}
generate_pstep_kernel = function(graph,step=5,save_rdata=TRUE){
    output_path="pstep_kernel"
    path_csv = here("src","Data","Kernels",paste0(output_path, ".csv"))
    path_rdata = here("src","Data","Kernels","RData",paste0(output_path, ".Rdata"))
    
    cat("Calculating p-step (random walk) kernel for the graph...\n")

    pstep_kernel = pStepKernel(graph,p=step)
    cat("p-step kernel calculated successfully!\n\n")
    write.csv(pstep_kernel, file =path_csv)
    cat("p-step kernel saved to file:", path_csv, "\n")

    if(save_rdata){
        save(pstep_kernel, file = path_rdata)
        cat("R object 'pstep_kernel' saved to:", path_rdata, "\n")
    }
    invisible(pstep_kernel)

}
generate_regularised_laplacian_kernel = function(graph,normalized=TRUE,save_rdata=TRUE){
    output_path="regularised_laplacian_kernel"
    path_csv = here("src","Data","Kernels",paste0(output_path, ".csv"))
    path_rdata = here("src","Data","Kernels","RData",paste0(output_path, ".Rdata"))

    cat("Calculating regularised laplacian kernel for the graph...\n")
    regularised_laplacian_kernel = regularisedLaplacianKernel(graph,normalized)
    cat("Regularised laplacian kernel calculated successfully!\n\n")
    write.csv(regularised_laplacian_kernel, file =path_csv)
    if(save_rdata){
        save(regularised_laplacian_kernel, file = path_rdata)
        cat("R object 'regularised_laplacian_kernel' saved to:", path_rdata, "\n")
    }
     invisible(regularised_laplacian_kernel)
}
generate_commute_time_kernel = function(graph,normalized=TRUE,save_rdata=TRUE){
    output_path="commute_time_kernel"
    path_csv = here("src","Data","Kernels",paste0(output_path, ".csv"))
    path_rdata = here("src","Data","Kernels","RData",paste0(output_path, ".Rdata"))

    cat("Calculating commute time kernel for the graph...\n")
    commute_time_kernel = commuteTimeKernel(graph,normalized)
    cat("Commute time kernel calculated successfully!\n\n")
    write.csv(commute_time_kernel, file =path_csv)
    if(save_rdata){
        save(commute_time_kernel, file = path_rdata)
        cat("Commute time kernel saved to file:", path_csv, "\n")
    }
     invisible(commute_time_kernel)
}
generate_inverse_cosine_kernel = function(graph,save_rdata=TRUE){
    output_path="inverse_cosine_kernel"
    path_csv = here("src","Data","Kernels",paste0(output_path, ".csv"))
    path_rdata = here("src","Data","Kernels","RData",paste0(output_path, ".Rdata"))

     cat("Calculating inverse cosine kernel for the graph...\n")
    inverse_cosine_kernel = inverseCosineKernel(graph)
    cat("Inverse cosine kernel calculated successfully!\n\n")
    if(save_rdata){
        save(inverse_cosine_kernel, file = path_rdata)
        cat("R object 'inverse_cosine_kernel' saved to:", path_rdata, "\n")
    }
     invisible(inverse_cosine_kernel)
}
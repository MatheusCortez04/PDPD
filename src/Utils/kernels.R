library(diffuStats)
library(igraph)
source(here("src","Utils","utils.R"))
generate_difussion_kernel= function(graph,normalized=TRUE,save_rdata=TRUE){
    output_file_name="diffusion_kernel"
    output_dir=here("src","Data","Kernels")
    create_dir(output_dir)
    output_path_csv =here(output_dir,paste0(output_file_name, ".csv"))
    
    
    cat("Calculating diffusion kernel for the graph...\n")
    diffusion_kernel = diffusionKernel(graph,normalized)

    cat("Diffusion kernel calculated successfully!\n\n")
    write.csv(diffusion_kernel, file =output_path_csv)
    cat("Csv file  saved to:", output_path_csv, "\n")
    Sys.sleep(1.5)

    if(save_rdata){
        rdata_dir=here("src","Data","Kernels","RData")
        create_dir(rdata_dir)
        output_path_rdata = here(rdata_dir,paste0(output_file_name, ".Rdata"))
        save(diffusion_kernel, file = output_path_rdata)
        cat("R object 'diffusion_kernel' saved to:", output_path_rdata, "\n")
        Sys.sleep(1.5)
    }
   invisible(diffusion_kernel)
}
generate_pstep_kernel = function(graph,step=5,save_rdata=TRUE){
    output_file_name="pstep_kernel"    
    output_dir=here("src","Data","Kernels")
    create_dir(output_dir)
    output_path_csv =here(output_dir,paste0(output_file_name, ".csv"))
    cat("Calculating p-step (random walk) kernel for the graph...\n")

    pstep_kernel = pStepKernel(graph,p=step)
    cat("p-step kernel calculated successfully!\n\n")
    write.csv(pstep_kernel, file =output_path_csv)
    cat("Csv file  saved to:", output_path_csv, "\n")
    Sys.sleep(1.5)
    if(save_rdata){
        rdata_dir=here("src","Data","Kernels","RData")
        create_dir(rdata_dir)
        output_path_rdata = here(rdata_dir,paste0(output_file_name, ".Rdata"))
        save(pstep_kernel, file = output_path_rdata)
        cat("R object 'pstep_kernel' saved to:", output_path_rdata, "\n")
        Sys.sleep(1.5)
    }
    invisible(pstep_kernel)

}
generate_regularised_laplacian_kernel = function(graph,normalized=TRUE,save_rdata=TRUE){
    output_file_name="regularised_laplacian_kernel"
    output_dir=here("src","Data","Kernels")
    create_dir(output_dir)
    output_path_csv =here(output_dir,paste0(output_file_name, ".csv"))

    cat("Calculating regularised laplacian kernel for the graph...\n")
    regularised_laplacian_kernel = regularisedLaplacianKernel(graph,normalized)
    cat("Regularised laplacian kernel calculated successfully!\n\n")
    
    write.csv(regularised_laplacian_kernel, file =output_path_csv)
    cat("Csv file  saved to:", output_path_csv, "\n")
    Sys.sleep(1.5)
    if(save_rdata){
        rdata_dir=here("src","Data","Kernels","RData")
        create_dir(rdata_dir)
        output_path_rdata = here(rdata_dir,paste0(output_file_name, ".Rdata"))
        save(regularised_laplacian_kernel, file = output_path_rdata)
        cat("R object 'regularised_laplacian_kernel' saved to:", output_path_rdata, "\n")
        Sys.sleep(1.5)
    }
     invisible(regularised_laplacian_kernel)
}
generate_commute_time_kernel = function(graph,normalized=TRUE,save_rdata=TRUE){
    output_file_name="commute_time_kernel"
    output_dir=here("src","Data","Kernels")
    create_dir(output_dir)
    output_path_csv =here(output_dir,paste0(output_file_name, ".csv"))

    cat("Calculating commute time kernel for the graph...\n")
    commute_time_kernel = commuteTimeKernel(graph,normalized)
    cat("Commute time kernel calculated successfully!\n\n")
    min_value = min(commute_time_kernel)
    commute_time_kernel = commute_time_kernel - min_value
    write.csv(commute_time_kernel, file =output_path_csv)
    cat("Csv file  saved to:", output_path_csv, "\n")
    Sys.sleep(1.5)

    if(save_rdata){
        rdata_dir=here("src","Data","Kernels","RData")
        create_dir(rdata_dir)
        output_path_rdata = here(rdata_dir,paste0(output_file_name, ".Rdata"))
        save(commute_time_kernel, file = output_path_rdata)
        cat("R object 'commute_time_kernel' saved to:", output_path_rdata, "\n")
        Sys.sleep(1.5)
    }
     invisible(commute_time_kernel)
}
generate_inverse_cosine_kernel = function(graph,save_rdata=TRUE){
    output_file_name="inverse_cosine_kernel"
    output_dir=here("src","Data","Kernels")
    create_dir(output_dir)
    output_path_csv =here(output_dir,paste0(output_file_name, ".csv"))

    cat("Calculating inverse cosine kernel for the graph...\n")
    inverse_cosine_kernel = inverseCosineKernel(graph)
    write.csv(inverse_cosine_kernel, file =output_path_csv)
    cat("Inverse cosine kernel calculated successfully!\n\n")
    cat("Csv file  saved to:", output_path_csv, "\n")
    Sys.sleep(1.5)
    if(save_rdata){
        rdata_dir=here("src","Data","Kernels","RData")
        create_dir(rdata_dir)
        output_path_rdata = here(rdata_dir,paste0(output_file_name, ".Rdata"))
        save(inverse_cosine_kernel, file = output_path_rdata)
        cat("R object 'inverse_cosine_kernel' saved to:", output_path_rdata, "\n")
        Sys.sleep(1.5)
    }
     invisible(inverse_cosine_kernel)
}

kernel_function_mapper <- list(
    '1' = function(graph) {
        cat("\n--- Generating Diffusion Kernel ---\n")
        save_rdata = readline(prompt = "Enter save RData (default TRUE): ")
        is_valid_save_rdata=  is_valid_input_boolean(save_rdata)
        if(!is_valid_save_rdata){
             cat("\n--- Invalid input. Using default value ---\n")
           save_rdata = TRUE
        }
        generate_difussion_kernel(graph,save_rdata=save_rdata)
    },
    '2' = function(graph) {
        cat("\n--- Generating P-Step Kernel ---\n")
        step_input = readline(prompt = "Enter number of steps (default 5): ")
        save_rdata = readline(prompt = "Enter save RData (default TRUE): ")
        is_valid_save_rdata=  is_valid_input_boolean(save_rdata)
        if(!is_valid_save_rdata){
            cat("\n--- Invalid input. Using default value ---\n")
            save_rdata = TRUE
        }
        step <- as.integer(step_input)
        if (is.na(step) || step <= 0) { 
            step <- 5
            cat("(Using default: 5 steps)\n")
        }
        generate_pstep_kernel(graph, step = step,save_rdata)
    },
    '3' = function(graph) {
        cat("\n--- Generating Regularised Laplacian Kernel ---\n")
        save_rdata = readline(prompt = "Enter save RData (default TRUE): ")
        is_valid_save_rdata=  is_valid_input_boolean(save_rdata)
        if(!is_valid_save_rdata){
            cat("\n--- Invalid input. Using default value ---\n")
            save_rdata = TRUE
        }
        generate_regularised_laplacian_kernel(graph)
    },
    '4' = function(graph) {
        cat("\n--- Generating Commute Time Kernel ---\n")
        generate_commute_time_kernel(graph)
    },
    '5' = function(graph) {
        cat("\n--- Generating Inverse Cosine Kernel ---\n")
        save_rdata = readline(prompt = "Enter save RData (default TRUE): ")
        is_valid_save_rdata=  is_valid_input_boolean(save_rdata)
        if(!is_valid_save_rdata){
            cat("\n--- Invalid input. Using default value ---\n")
            save_rdata = TRUE
        }
        generate_inverse_cosine_kernel(graph,save_rdata)
    }
)

validate_kernel_input = function(input){
    valid_inputs = c("B","1","2","3","4","5")
    return (toupper(trimws(input)) %in% valid_inputs)
}

generate_kernel_menu <- function() {
    ppi_df  = read.csv(here("src","Data","PPI_gysi.csv"), sep=",")
    cat("PPI reading complete.\n")
    cat("Creating dataframe with only proteins\n")
    ppi_only_proteins_df =  ppi_df %>% select(proteinA_entrezid,proteinB_entrezid)
    ppi_graph = generate_graph_from_dataframe(ppi_only_proteins_df)
  while (TRUE) {
    clear_console()
    cat("--- Kernel Generation Menu ---\n\n")
    cat(" [1] Diffusion Kernel\n")
    cat(" [2] P-Step Kernel\n")
    cat(" [3] Regularised Laplacian Kernel\n")
    cat(" [4] Commute Time Kernel\n")
    cat(" [5] Inverse Cosine Kernel\n")
    cat(" [B] Back\n\n")
    input = readline(prompt = "Choice option: ")
    is_valid = validate_kernel_input(input)
    if (!is_valid) {
      cat("\n[Erro] Opção inválida. Tente novamente.\n")
      Sys.sleep(1.5)
      next
    }
    if (toupper(trimws(input)) == "B") break

     kernel_function_mapper[[input]](ppi_graph)
  }
}

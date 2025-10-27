library(utils)
library(here)
library(dplyr)

source(here("src","Utils","graph.R"))

source(here("src","Utils","kernels.R"))
source(here("src","Utils","utils.R"))
main= function(){
    cat("\nLibraries loaded successfully.\n")
    while(TRUE){
        clear_console()
        cat("--- R Kernel Generation Menu ---\n\n")
        cat("Select a kernel to generate:\n")
        cat(" [1] Diffusion Kernel\n")
        cat(" [2] P-Step Kernel\n")
        cat(" [3] Regularised Laplacian Kernel\n")
        cat(" [4] Commute Time Kernel\n")
        cat(" [5] Inverse Cosine Kernel\n")
        cat("\n [q] Quit Application\n\n")
    input =  readline(prompt = "Enter your choice: ")

    is_valid_input =validate_input(input)
    if(!is_valid_input){
        cat("\n[Error] Invalid option. Please try again.\n")
        Sys.sleep(1.5)
        next
    }
    if( toupper(trimws(input))=="Q"){
        break
    }
    ppi_df  = read.csv(here("src","Data","PPI_gysi.csv"), sep=",")
    cat("PPI reading complete.\n")
    cat("Creating dataframe with only proteins\n")
    ppi_only_proteins_df =  ppi_df %>% select(proteinA_entrezid,proteinB_entrezid)
    ppi_graph = generate_graph_from_dataframe(ppi_only_proteins_df)
    kernel_function_mapper[[input]](ppi_graph)
    
    }

    q(save="no")
    }









library(utils)
library(here)
library(tidyverse)
source(here("src","Utils","graph.R"))
source(here("src","Utils","kernels.R"))
source(here("src","Utils","drug.R"))
main= function(){

    cat("\nLibraries loaded successfully.\n")
        while(TRUE){
            cat("--- Main Menu ---\n\n")
            cat(" [1] Kernel Menu\n")
            cat(" [2] Drug Menu\n")
            cat(" [Q] Exit\n\n")
            input =  readline(prompt = "Enter your choice: ")
            if (toupper(trimws(input)) == "Q") break
            if (input == "1") generate_kernel_menu()
            if (input == "2") scoring_drug_disease_menu()
        }

    q(save="no") 
}








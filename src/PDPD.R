
library(utils)
library(here)
library(tidyverse)
source(here("src","Utils","graph.R"))
source(here("src","Utils","kernels.R"))
source(here("src","Utils","drug.R"))
source(here("src","Utils","evaluation.R"))
main= function(){

    # cat("\nLibraries loaded successfully.\n")
    #     while(TRUE){
    #         cat("--- Main Menu ---\n\n")
    #         cat(" [1] Kernel Menu\n")
    #         cat(" [2] Drug Menu\n")
    #         cat(" [3] Evaluation Menu\n")
    #         cat(" [Q] Exit\n\n")
    #         input =  readline(prompt = "Enter your choice: ")
    #         if (toupper(trimws(input)) == "Q") break
    #         if (input == "1") generate_kernel_menu()
    #         if (input == "2") scoring_drug_disease_menu()
    #         if (input == "3") evaluation_menu()
    #     }
    # mdd_genes_vector =get_mdd_genes()
#      calculate_drug_score()
#      generate_drug_rank()
# created_prediction_mdd_data()

#      generate_roc_curve_mdd()
calculate_drug_score()
generate_drug_rank("BD")
generate_drug_rank("MDD")
created_prediction_bipolar_data()
created_prediction_mdd_data()
generate_roc_curve_bipolar()
generate_roc_curve_mdd()
    # targets_validos <- intersect(alvos_teste, rownames(kernel))
    # disease_genes_validos <- intersect(mdd_genes_vector$entrez_id, rownames(kernel))
    # score_droga_teste = kernel[targets_validos,disease_genes_validos]
   

    q(save="no") 
}








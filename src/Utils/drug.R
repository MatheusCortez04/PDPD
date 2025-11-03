library(here)
source(here("src","Utils","utils.R"))

generate_ordered_drug_protein_matrix = function(drug_target_df,protein_nodes,drug_nodes){
  cat("\n--- Generating Drug Protein Matrix ---\n")
    output_path="drug_gene"
    path_csv = here("src","Data","Drug",paste0(output_path, ".csv"))
    path_rdata = here("src","Data","Drug","RData",paste0(output_path, ".Rdata"))

    matrix_drug_protein = drug_target_df %>% 
    filter( entrez_id %in% protein_nodes) %>%
    distinct(drugbank_id, entrez_id) %>%
    mutate(
      drugbank_id = factor(drugbank_id, levels = drug_nodes),
      entrez_id = factor(entrez_id, levels = protein_nodes),
      value=1
    ) %>%   
    pivot_wider(
        id_cols = drugbank_id, 
        names_from = entrez_id,
        values_from = value,
        values_fill = 0,
        names_expand = TRUE,
        id_expand = TRUE
    ) %>% arrange(drugbank_id)

    write.csv(matrix_drug_protein, file =path_csv,row.names=FALSE)
    cat("Csv file  saved to:", path_csv, "\n")
    Sys.sleep(1)
    matrix_drug_protein_numeric <- as.matrix(matrix_drug_protein[ , -1])
    rownames(matrix_drug_protein_numeric) <- matrix_drug_protein$drugbank_id
    storage.mode(matrix_drug_protein_numeric) <- "numeric"

    save(matrix_drug_protein_numeric, file = path_rdata)
    cat("R object 'drug_target_df' saved to:", path_rdata, "\n")
    Sys.sleep(1)
}



scoring_drug_disease_menu = function(){
  while(TRUE){
    clear_console()
    cat("--- Drug Score  Menu ---\n\n")
    cat(" [1] Generate Drug Protein Matrix\n")
    cat(" [2] Generate Kernel Score \n")
    cat(" [B] Back\n\n")

    input = readline(prompt = "Choice option: ")
    if (toupper(trimws(input)) == "B") break

     drug_function_mapper[[input]]()
  }
}

drug_function_mapper <- list(
    '1' = function() {
        cat("\n--- Generating and  Drug X Protein Matrix  ---\n")
        drug_target_df  = read.csv(here("src","Data","drug_targets_DrugBank_Gysi.csv"), sep=",")
        cat("Drug Target file reading complete.\n")
        ppi_df  = read.csv(here("src","Data","PPI_gysi.csv"), sep=",")
        cat("PPI  file reading complete.\n")
        protein_nodes = get_ppi_nodes(ppi_df)
        drug_nodes= get_drug_nodes(drug_target_df)
        generate_ordered_drug_protein_matrix(drug_target_df,protein_nodes,drug_nodes)
    },
    '2' = function(){
        
        calculate_kernel_score()
    }
)

calculate_kernel_score = function(){
drug_gene__matrix_file_path="drug_gene"
kernel_file_names = c('diffusion_kernel','pstep_kernel','regularised_laplacian_kernel','commute_time_kernel','inverse_cosine_kernel')
drug_protein_file_path = here("src","Data","Drug","RData",paste0(drug_gene__matrix_file_path, ".Rdata"))

drug_protein_file_already_exists = file.exists(drug_protein_file_path)
if(!drug_protein_file_already_exists){
    cat("\n[Error]: Required file Drug gene matrix not found. This file is necessary to calculate the score. Would you like to create it now\n")
     Sys.sleep(1.5)
    return()
}

drug_protein_matrix = load_rdata(drug_protein_file_path)
mdd_vector = build_protein_mdd_df()
bipolar_vector = mdd_vector = build_protein_bipolar_df()
for (kernel_name in kernel_file_names){
    cat("\n--- Generating Score Kernel ---\n")

    kernel_path_rdata = here("src","Data","Kernels","RData", paste0(kernel_name, ".Rdata"))
    kernel_file_exists = file.exists(kernel_path_rdata)
    if(!kernel_file_exists){
        cat("Kernel file",kernel_name,"no exists\n")
        Sys.sleep(1)
        next()
    }

    kernel = load_rdata(kernel_path_rdata)

    cat("Starting mdd matrix multiplication...\n")
    mdd_score = drug_protein_matrix %*% kernel %*% mdd_vector$is_disease 
    cat("Matrix multiplication completed!\n")
    
    output_mdd_file_path_Rdata = here("src","Data","Kernels","RData", paste0(kernel_name, "_mdd_score.Rdata"))
    output_mdd_file_path_csv = here("src","Data","Kernels", paste0(kernel_name, "_mdd_score.csv"))
    
    write.csv(mdd_score, file = output_mdd_file_path_csv)
    cat("CSV score saved to:", output_mdd_file_path_csv, "\n")
    save(mdd_score, file = output_mdd_file_path_Rdata)
    cat("R object score saved to:", output_mdd_file_path_Rdata, "\n")


    cat("Starting bipolar matrix multiplication...\n")
    bipolar_score = drug_protein_matrix %*% kernel %*% bipolar_vector$is_disease 
    cat("Matrix multiplication completed!\n")
    output_bipolar_file_path_Rdata = here("src","Data","Kernels","RData", paste0(kernel_name, "_bipolar_score.Rdata"))
    output_bipolar_file_path_csv = here("src","Data","Kernels", paste0(kernel_name, "_bipolar_score.csv"))
    
    write.csv(bipolar_score, file = output_bipolar_file_path_csv)
    cat("CSV score saved to:", output_bipolar_file_path_csv, "\n")
    save(bipolar_score, file = output_bipolar_file_path_Rdata)
    cat("R object score saved to:", output_bipolar_file_path_Rdata, "\n")
    
    }


} 

 



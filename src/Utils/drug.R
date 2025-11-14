library(here)
source(here("src","Utils","utils.R"))

generate_ordered_drug_protein_matrix = function(drug_target_df,protein_nodes,drug_nodes){
  cat("\n--- Generating Drug Target Matrix ---\n")
    output_file_name="drug_target_matrix"
    drug_dir =here("src","Data","Drug")
    drug_Rdata_dir =here("src","Data","Drug","RData")
    create_dir(drug_dir)
    create_dir(drug_Rdata_dir)
    
    output_path_csv = here(drug_dir,paste0(output_file_name, ".csv"))
    output_path_rdata = here(drug_Rdata_dir,paste0(output_file_name, ".Rdata"))


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

    write.csv(matrix_drug_protein, file =output_path_csv,row.names=FALSE)
    cat("Csv file  saved to:", output_path_csv, "\n")
    Sys.sleep(1)
    matrix_drug_protein_numeric <- as.matrix(matrix_drug_protein[ , -1])
    rownames(matrix_drug_protein_numeric) <- matrix_drug_protein$drugbank_id
    storage.mode(matrix_drug_protein_numeric) <- "numeric"

    save(matrix_drug_protein_numeric, file = output_path_rdata)
    cat("R object 'drug_target_df' saved to:", output_path_rdata, "\n")
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
        cat("\n--- Generating and  Drug  Target Matrix  ---\n")
        drug_target_df  = read.csv(here("src","Data","Drug","drug_targets_DrugBank_Gysi.csv"), sep=",")

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
    drug_gene__matrix_file_path="drug_target_matrix"
    kernel_file_names = c('diffusion_kernel','pstep_kernel','regularised_laplacian_kernel','commute_time_kernel','inverse_cosine_kernel')
    drug_target_rdata_file_path = here("src","Data","Drug","RData",paste0(drug_gene__matrix_file_path, ".Rdata"))

    drug_protein_file_already_exists = file.exists(drug_target_rdata_file_path)
    if(!drug_protein_file_already_exists){
        cat("\n[Error]: Required file Drug target matrix not found.\n This file is necessary to calculate the score. Would you like to create it now\n")
        Sys.sleep(1.5)
        return()
    }

    drug_protein_matrix = load_rdata(drug_target_rdata_file_path)
    mdd_vector = build_protein_mdd_df()
    bipolar_vector =build_protein_bipolar_df()

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
        mdd_score_df = data.frame(
            Drug_Id =  rownames(drug_protein_matrix),
            Kernel_Score = as.numeric(mdd_score)

        )
        mdd_score_df = mdd_score_df %>% arrange(desc(Kernel_Score))
        cat("Matrix multiplication completed!\n")


        kernel_score_dir = here("src","Data","Kernels","Score")
        mdd_score_dir =here(kernel_score_dir,"MDD")
        create_dir(mdd_score_dir)

        mdd_score_rdata_dir = here(mdd_score_dir,"RData")
        create_dir(mdd_score_rdata_dir)

        
        output_mdd_file_path_csv = here(mdd_score_dir, paste0(kernel_name, ".csv"))
        output_mdd_file_path_Rdata = here(mdd_score_rdata_dir,paste0(kernel_name, ".Rdata"))

        write.csv(mdd_score_df, file = output_mdd_file_path_csv,row.names=FALSE)
        cat("CSV score saved to:", output_mdd_file_path_csv, "\n")
        save(mdd_score_df, file = output_mdd_file_path_Rdata)
        cat("R object score saved to:", output_mdd_file_path_Rdata, "\n")


        cat("Starting bipolar matrix multiplication...\n")
        bipolar_score = drug_protein_matrix %*% kernel %*% bipolar_vector$is_disease 
        bipolar_score_df = data.frame(
            Drug_Id =  rownames(drug_protein_matrix),
            Kernel_Score = as.numeric(bipolar_score)

        )
        bipolar_score_df = bipolar_score_df %>% arrange(desc(Kernel_Score))
        cat("Matrix multiplication completed!\n")
        bipolar_score_dir = here(kernel_score_dir,"BD")
        create_dir(bipolar_score_dir)
        
        bipolar_score_rdata_dir = here(bipolar_score_dir,"RData")
        create_dir(bipolar_score_rdata_dir)

        output_bipolar_file_path_csv = here(bipolar_score_dir, paste0(kernel_name, ".csv"))
        output_bipolar_file_path_Rdata =here(bipolar_score_rdata_dir,paste0(kernel_name, ".Rdata")) 

        write.csv(bipolar_score_df, file = output_bipolar_file_path_csv,row.names=FALSE)
        cat("CSV score saved to:", output_bipolar_file_path_csv, "\n")
        save(bipolar_score_df, file = output_bipolar_file_path_Rdata)
        cat("R object score saved to:", output_bipolar_file_path_Rdata, "\n")

        }
} 


generate_drug_kernel_bipolar_score = function(){
    difussion_kernel_score = load_rdata(here("src","Data","Kernels","RData","Score", paste0("diffusion_kernel_bipolar_score.Rdata")))
    pstep_kernel_score = load_rdata(here("src","Data","Kernels","RData","Score", paste0("pstep_kernel_bipolar_score.Rdata")))
    regularised_laplacian_kernel_bipolar_score = load_rdata(here("src","Data","Kernels","RData","Score", paste0("regularised_laplacian_kernel_bipolar_score.Rdata")))
    inverse_cosine_kernel_bipolar_score = load_rdata(here("src","Data","Kernels","RData","Score", paste0("inverse_cosine_kernel_bipolar_score.Rdata")))
    commute_time_kernel_bipolar_score = load_rdata(here("src","Data","Kernels","RData","Score", paste0("commute_time_kernel_bipolar_score.Rdata")))
    merge_df = difussion_kernel_score %>%
        left_join(pstep_kernel_score,by = "Drug_Id") %>%
        left_join(regularised_laplacian_kernel_bipolar_score,by = "Drug_Id") %>%
        left_join(inverse_cosine_kernel_bipolar_score,by="Drug_Id") %>%
        left_join(commute_time_kernel_bipolar_score,by="Drug_Id") %>%
        rename(difussion_kernel_score=Kernel_Score.x,
                pstep_kernel_score=Kernel_Score.y,
                regularised_laplacian_kernel_score = Kernel_Score.x.x,
                inverse_cosine_kernel_score = Kernel_Score.y.y,
                commute_time_kernel_score=Kernel_Score
                )

    output_bipolar_Score_file_path_csv = here("src","Data","bipolar_kernel_score.csv")
    write.csv(merge_df,file=output_bipolar_Score_file_path_csv,row.names=FALSE)
    cat("CSV score saved to:", output_bipolar_Score_file_path_csv, "\n")
}

generate_drug_kernel_mdd_score = function(){
    difussion_kernel_score = load_rdata(here("src","Data","Kernels","RData","Score", paste0("diffusion_kernel_mdd_score.Rdata")))
    pstep_kernel_score = load_rdata(here("src","Data","Kernels","RData","Score", paste0("pstep_kernel_mdd_score.Rdata")))
    regularised_laplacian_kernel_bipolar_score = load_rdata(here("src","Data","Kernels","RData","Score", paste0("regularised_laplacian_kernel_mdd_score.Rdata")))
    inverse_cosine_kernel_bipolar_score = load_rdata(here("src","Data","Kernels","RData","Score", paste0("inverse_cosine_kernel_bipolar_score.Rdata")))
    commute_time_kernel_bipolar_score = load_rdata(here("src","Data","Kernels","RData","Score", paste0("commute_time_kernel_mdd_score.Rdata")))
    merge_df = difussion_kernel_score %>%
        left_join(pstep_kernel_score,by = "Drug_Id") %>%
        left_join(regularised_laplacian_kernel_bipolar_score,by = "Drug_Id") %>%
        left_join(inverse_cosine_kernel_bipolar_score,by="Drug_Id") %>%
        left_join(commute_time_kernel_bipolar_score,by="Drug_Id") %>%
        rename(difussion_kernel_score=Kernel_Score.x,
                pstep_kernel_score=Kernel_Score.y,
                regularised_laplacian_kernel_score = Kernel_Score.x.x,
                inverse_cosine_kernel_score = Kernel_Score.y.y,
                commute_time_kernel_score=Kernel_Score
                )

    output_mdd_score_file_path_csv = here("src","Data","mdd_kernel_score.csv")
    write.csv(merge_df,file=output_mdd_score_file_path_csv,row.names=FALSE)
    cat("CSV score saved to:", output_mdd_score_file_path_csv, "\n")
}

 



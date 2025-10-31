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
      cat("\n--- Generating and save Diseases Vector  ---\n")
       bipolar_gene= get_bipolar_disorder_genes()

            write.csv(result,"bipolar_gene.csv",)
    }
)



        #  drug_protein = load_rdata(here("src","Data","drug_gene.Rdata"))
        #  kernel  = load_rdata(here("src","Data","Kernels","RData","diffusion_kernel.Rdata"))

        # result <- drug_protein %*% kernel


        # write.csv(result,"teste.csv")
library(here)
source(here("src","Utils","utils.R"))

generate_ordered_drug_protein_matrix = function(drug_target_df,protein_nodes,drug_nodes){
    output_path="drug_gene"
    path_csv = here("src","Data",paste0(output_path, ".csv"))
    path_rdata = here("src","Data",paste0(output_path, ".Rdata"))

    matrix_final = drug_target_df %>% 
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

    write.csv(matrix_final, file =path_csv)
    cat("Csv file  saved to:", path_csv, "\n")

    save(matrix_final, file = path_rdata)
    cat("R object 'drug_target_df' saved to:", path_rdata, "\n")

}



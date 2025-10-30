library(here)

clear_console = function(){
    system("clear")
}

is_valid_input_boolean = function(input){
    return(toupper(input) == 'TRUE' || toupper(input) == 'FALSE')
}


get_ppi_nodes = function(ppi_df){
    proteinA_entrezid = ppi_df$proteinA_entrezid
    proteinB_entrezid = ppi_df$proteinB_entrezid
    all_protein_in_ppi_df = c(proteinA_entrezid,proteinB_entrezid)
    print(paste("All proteins in PPI Dataframe: ",length(all_protein_in_ppi_df)))
    unique_ordered_proteins_in_ppi_df = unique(all_protein_in_ppi_df)
    print(paste("Unique proteins in PPI Dataframe: ",length(unique_ordered_proteins_in_ppi_df)))
    invisible(unique_ordered_proteins_in_ppi_df)
}
get_drug_nodes = function(drug_target_df){
    drug_ids = drug_target_df$drugbank_id
    print(paste("All drugs in Drug to target Dataframe: ",length(drug_ids)))
    unique_drugs_in_df = unique(drug_ids)
    print(paste("Unique drugs in Drug to target Dataframe: ",length(unique_drugs_in_df)))
    invisible(unique_drugs_in_df)
}

main_menu_mapper = function(){
    cat("--- Main Menu ---\n\n")
    cat(" [1] Gerar e Salvar Kernel\n")
    cat(" [2] Scoring Droga-Doença\n")
    cat(" [q] Sair\n\n")
    input = readline(prompt = "Escolha a opção: ")
    is_valid = validate_input_main(input)
    if (!is_valid) {
        cat("\n[Erro] Opção inválida. Tente novamente.\n")
        Sys.sleep(1.5)
        next
    }
    if (toupper(trimws(input)) == "Q") break
    if (input == "1") gerar_kernel_menu()
    if (input == "2") scoring_drug_disease_menu()
}




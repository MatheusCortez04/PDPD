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

load_rdata <- function(path_file) {
  env <- new.env()
  load(path_file, envir = env)
  objs <- ls(env)
  if(length(objs) != 1) {
    stop("Mais de um objeto no arquivo RData")
  }
  return(env[[objs]])
}





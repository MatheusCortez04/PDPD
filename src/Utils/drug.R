library(here)
library(purrr)
library(dplyr)
scoring_drug_disease_menu = function(){
  while(TRUE){
    clear_console()
    cat("--- Drug Menu ---\n\n")
    cat(" [1] Generate Drug Score\n")
    cat(" [B] Back\n\n")

    input = readline(prompt = "Choice option: ")
    if (toupper(trimws(input)) == "B") break

     drug_function_mapper[[input]]()
  }
}

drug_function_mapper = list(
    '1' = function() {
    }
)

get_drugbank_ids = function(){
    drug_target_df  = read.csv(here("src","Data","Drug","drug_targets_DrugBank_Gysi.csv")) %>%
        dplyr::distinct()

    unique_drugbank_ids = drug_target_df %>%
        dplyr::pull(drugbank_id) %>%
        unique()

    print(paste("Unique drugs in Drug-Target Dataframe:", length(unique_drugbank_ids)))
    invisible(unique_drugbank_ids)
}
get_drug_targets_in_ppi <- function(drugbank_id,ppi_gene_nodes) {
  
    drug_target_df = read.csv(here("src","Data","Drug","drug_targets_DrugBank_Gysi.csv")) %>% 
        distinct() 
  

    drug_target_proteins = drug_target_df %>%
        filter(drugbank_id == !!drugbank_id) %>%
        distinct() %>%
        dplyr::mutate(entrez_id = as.character(entrez_id)) %>%
        pull(entrez_id)

    drug_targets_in_ppi <- intersect(drug_target_proteins, ppi_gene_nodes)

    invisible(drug_targets_in_ppi)
}

generate_drug_rank = function(){
    kernel_names <- c(
    "diffusion_kernel","pstep_kernel","regularised_laplacian_kernel",
    "commute_time_kernel","inverse_cosine_kernel"
  )


    drugbank_ids = get_drugbank_ids()
    ppi_gene_nodes <- get_ppi_nodes()
    drug_targets = purrr::map(drugbank_ids,get_drug_targets_in_ppi,ppi_gene_nodes)


    names(drug_targets) = drugbank_ids
    glimpse(drug_targets)

    # for(drugbank_id in drugbank_ids){

    # }

}
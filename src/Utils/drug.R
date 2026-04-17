library(here)
library(purrr)
library(dplyr)
source(here("src","Utils","utils.R"))
scoring_drug_disease_menu = function(){
  while(TRUE){
    clear_console()
    cat("--- Drug Menu ---\n\n")
    cat(" [1] Generate Drug Rank\n")
    cat(" [B] Back\n\n")

    input = readline(prompt = "Choice option: ")
    if (toupper(trimws(input)) == "B") break

     drug_function_mapper[[input]]()
  }
}

drug_function_mapper = list(
    '1' = function() {
        calculate_drug_score()
        generate_drug_rank("MDD")
        generate_drug_rank("BD")
    }
)

get_drugbank_ids = function(){
    drug_target_df  = read_drug_targets()

    unique_drugbank_ids = drug_target_df %>%
        dplyr::pull(drugbank_id) %>%
        unique()

    cat(sprintf("[INFO] Unique drugs in Drug-Target Dataframe: %d\n",length(unique_drugbank_ids)))
    invisible(unique_drugbank_ids)
}
get_drug_targets = function(drugbank_id,ppi_gene_nodes,drug_target_df) {

    drug_target_proteins = drug_target_df %>%
        filter(drugbank_id == !!drugbank_id) %>%
        distinct() %>%
        dplyr::mutate(
            drugbank_id = as.character(drugbank_id),
            entrez_id   = as.character(entrez_id)) %>%
        pull(entrez_id)

    drug_target_proteins = intersect(drug_target_proteins, ppi_gene_nodes)
    cat(sprintf("[INFO] For DrugBank ID %s, found %d unique targets.\n",drugbank_id,length(drug_target_proteins)))
    invisible(drug_target_proteins)
}

calculate_drug_score = function(){
    kernel_names = c(
        "diffusion_kernel","pstep_kernel","regularised_laplacian_kernel",
        "commute_time_kernel","inverse_cosine_kernel"
    )

    disease_data <- list(
        MDD = list(name = "MDD", genes = get_mdd_disease_module()$entrez_id),
        BD = list(name = "BD", genes = get_bipolar_disease_module()$entrez_id)
    )
    drugbank_ids = get_drugbank_ids()
    ppi_gene_nodes = get_ppi_nodes()
    drug_target_df=read_drug_targets()
    drug_targets = purrr::map(drugbank_ids, get_drug_targets, ppi_gene_nodes,drug_target_df)
    names(drug_targets) = drugbank_ids

}
read_drug_targets= function(){
    drug_target_df = read.csv(here("src","Data","Drug","drug_targets_DrugBank_Gysi.csv")) %>% 
        distinct() 
    invisible(drug_target_df)
}
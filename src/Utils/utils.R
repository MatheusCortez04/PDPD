library(here)
library(dplyr)
clear_console = function(){
    system("clear")
}

is_valid_input_boolean = function(input){
    return(toupper(input) == 'TRUE' || toupper(input) == 'FALSE')
}

get_ppi_nodes = function(){
    ppi_df  = read.csv(here("src","Data","PPI_gysi.csv"), sep=",")
    proteinA_entrezid = ppi_df$proteinA_entrezid
    proteinB_entrezid = ppi_df$proteinB_entrezid
    all_protein_in_ppi_df = c(proteinA_entrezid,proteinB_entrezid)
    cat(sprintf("[INFO] All proteins in PPI Dataframe: %d\n",length(all_protein_in_ppi_df)))
    unique_ordered_proteins_in_ppi_df = unique(all_protein_in_ppi_df)
    cat(sprintf("[INFO] Unique proteins in PPI Dataframe: %d\n",length(unique_ordered_proteins_in_ppi_df)))
    invisible(unique_ordered_proteins_in_ppi_df)
}

load_rdata <- function(path_file) {
  cat("Loading RData of file :", path_file, "\n")
  env <- new.env()
  load(path_file, envir = env)
  objs <- ls(env)
  return(env[[objs]])
}
get_mdd_disease_module = function() {
  
  mdd_disease_id = "C1269683"
  score_filter = 0.6
  
  disease_gene_df  = read.csv(here("src","Data","Disease","disease_genes.csv"), sep="\t")
  
  mdd_gene_module = disease_gene_df %>%
    filter(diseaseid == mdd_disease_id & score >= score_filter) %>%
    rename(entrez_id = geneid, disease_id = diseaseid) %>% 
    mutate(entrez_id = as.character(entrez_id)) %>% 
    select(entrez_id, disease_id, score)
  
  invisible(mdd_gene_module)
}
get_bipolar_disease_module = function() {
  
  bipolar_disease_id = "C0005586"
  score_filter = 0.6
  
  disease_gene_df  = read.csv(here("src","Data","Disease","disease_genes.csv"), sep="\t")
  
  bipolar_gene_module = disease_gene_df %>%
    filter(diseaseid == bipolar_disease_id & score >= score_filter) %>%
    rename(entrez_id = geneid, disease_id = diseaseid) %>% 
    mutate(entrez_id = as.character(entrez_id)) %>% 
    select(entrez_id, disease_id, score)
  
  invisible(bipolar_gene_module)
}

generate_drug_rank = function() {
    diseases = c("MDD", "BD")

    for(disease in diseases) {
        base::cat(base::paste("[WARN] Building Average Rank to Disease:", disease, "\n"))
        base_dir = here::here("src", "Data", "Drug", "Score", disease)
        
        pstep_df = read.csv(here::here(base_dir, "pstep_kernel.csv")) %>%
            dplyr::mutate(pstep_kernel = base::rank(dplyr::desc(max_gene_score), ties.method = "average")) %>%
            dplyr::select(drugbank_id, pstep_kernel)

        reg_lap_df = read.csv(here::here(base_dir, "regularised_laplacian_kernel.csv")) %>%
            dplyr::mutate(regularised_laplacian_kernel = base::rank(dplyr::desc(max_gene_score), ties.method = "average")) %>%
            dplyr::select(drugbank_id, regularised_laplacian_kernel)
        inv_cos_df = read.csv(here::here(base_dir, "inverse_cosine_kernel.csv")) %>%
           dplyr::mutate(inverse_cosine_kernel = base::rank(dplyr::desc(max_gene_score), ties.method = "average"))  %>%
            dplyr::select(drugbank_id, inverse_cosine_kernel)

        commute_df = read.csv(here::here(base_dir, "commute_time_kernel.csv")) %>%
             dplyr::mutate(commute_time_kernel = base::rank(dplyr::desc(max_gene_score), ties.method = "average"))  %>%
            dplyr::select(drugbank_id, commute_time_kernel)

        diffusion_df = read.csv(here::here(base_dir, "diffusion_kernel.csv")) %>%
             dplyr::mutate(diffusion_kernel = base::rank(dplyr::desc(max_gene_score), ties.method = "average"))  %>%
            dplyr::select(drugbank_id, diffusion_kernel)

        merge_df = pstep_df %>%
            dplyr::left_join(reg_lap_df, by = "drugbank_id") %>%
            dplyr::left_join(commute_df, by = "drugbank_id") %>%
            dplyr::left_join(inv_cos_df, by = "drugbank_id") %>%
            dplyr::left_join(diffusion_df, by = "drugbank_id") 
            

        final_df <- merge_df %>%
            dplyr::mutate(
                average_rank = rowMeans(
                    dplyr::select(., dplyr::ends_with("_kernel")), 
                    na.rm = TRUE
                )
            ) %>%
            dplyr::arrange(average_rank) %>%
            dplyr::select(drugbank_id,dplyr::ends_with("_kernel"),average_rank)

        write.csv(final_df, here(base_dir,paste0("average_rank_", disease, ".csv")), row.names = FALSE)
    }

}

create_dir <- function(path) {
  if (!dir.exists(path)) {
    message(paste("Creating Dir:", path))
    dir.create(path, recursive = TRUE, mode = "0777")
  }
}





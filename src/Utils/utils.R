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
    cat("PPI  file reading complete.\n")
    proteinA_entrezid = ppi_df$proteinA_entrezid
    proteinB_entrezid = ppi_df$proteinB_entrezid
    all_protein_in_ppi_df = c(proteinA_entrezid,proteinB_entrezid)
    print(paste("All proteins in PPI Dataframe: ",length(all_protein_in_ppi_df)))
    unique_ordered_proteins_in_ppi_df = unique(all_protein_in_ppi_df)
    print(paste("Unique proteins in PPI Dataframe: ",length(unique_ordered_proteins_in_ppi_df)))
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

summarise_drug_max_score = function(){
    diseases = c("MDD", "BD")
    kernel_names = c(
        "diffusion_kernel","pstep_kernel","regularised_laplacian_kernel",
        "commute_time_kernel","inverse_cosine_kernel"
    )
    for(disease in diseases){
        for(kernel in kernel_names){
            gene_score_path = here("src","Data","Kernels","Score",disease,paste0(kernel, ".csv"))
            if (!file.exists(gene_score_path)) next
             cat("\n✔Building Drug rank to Disease:", disease, "and kernel",kernel, "\n")
            gene_score_df = read.csv(gene_score_path)
            drug_target_df= load_drug_target_df()

            drug_gene_score_df = gene_score_df %>%
                dplyr::left_join(drug_target_df,by="entrez_id")%>%
                dplyr::filter(!is.na(drugbank_id) & drugbank_id != "") %>%
                dplyr::group_by(drugbank_id) %>%
                dplyr::summarise(
                    max_gene_score = max(gene_score, na.rm = TRUE),
                    top_entrez_id = entrez_id[which.max(gene_score)],
                    .groups = "drop" ) %>%  
                dplyr::arrange(desc(max_gene_score))
            out_dir  <- here("src","Data","Drug","Score", disease)
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
            out_path = here(out_dir,paste0(kernel, ".csv"))
            write.csv(drug_gene_score_df,out_path,row.names=FALSE)
            cat("\n✔ Drug Rank of Kernel:",kernel,"Saved in: ",out_path, "\n")
        }
        
    }
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


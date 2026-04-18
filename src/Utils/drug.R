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
        generate_drug_rank()
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
    invisible(drug_target_proteins)
}

calculate_drug_score = function(){
    
    cat(sprintf("\n[INFO] Starting drug scoring pipeline...\n"))
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

    cat(sprintf("[INFO] Total drugs: %d | PPI nodes: %d\n",
        length(drugbank_ids), length(ppi_gene_nodes)))

    drug_targets = purrr::map(drugbank_ids, get_drug_targets, ppi_gene_nodes,drug_target_df)
    names(drug_targets) = drugbank_ids

    for(kernel_name in kernel_names){
        
        cat(sprintf("\n[INFO] Processing kernel: %s\n", kernel_name))

        kernel_path = here("src","Data","Kernels","RData", paste0(kernel_name, ".Rdata"))

        if (!file.exists(kernel_path)) {
            cat("[WARN]Kernel:",kernel_name, "file not found. Skipping....\n")
            next()
        }

        kernel_matrix = load_rdata(kernel_path)
        for(disease in names(disease_data)){

            score_base_dir = here("src", "Data", "Drug", "Score", disease)
            create_dir(score_base_dir)

            disease_genes = disease_data[[disease]]$genes
            
            cat(sprintf("[INFO] Kernel: %s | Disease: %s | Genes: %d\n",kernel_name, disease, length(disease_genes)))

            scores = purrr::map_dbl(drug_targets,compute_drug_score,disease_genes,kernel_matrix)
            valid_scores = sum(!is.na(scores))
            na_scores = sum(is.na(scores))

            cat(sprintf("[INFO] Kernel: %s | Disease: %s | Valid scores: %d | NA scores: %d\n",kernel_name, disease, valid_scores, na_scores))


            result_df = data.frame(
            drugbank_id = names(scores),
            score = as.numeric(scores))
            
            output_file_name = here(score_base_dir, paste0(kernel_name, ".csv"))
            write.csv(result_df,file=output_file_name,row.names=FALSE)
            cat(sprintf("✔ [DONE] Kernel: %s | Disease: %s | Saved: %s\n",kernel_name, disease, output_file_name))
        }
    }

}

read_drug_targets = function(){
    drug_target_df = read.csv(here("src","Data","Drug","drug_targets_DrugBank_Gysi.csv")) %>% 
        distinct() 
    invisible(drug_target_df)
}

compute_drug_score = function(targets,disease_genes,kernel_matrix) {
  if (length(targets) <= 0) {
    return(NA_real_)
  }
  max(kernel_matrix[targets, disease_genes, drop = FALSE], na.rm = TRUE)
}

generate_drug_rank = function() {
    diseases = c("MDD", "BD")
    
    cat("\n[INFO] Starting drug ranking aggregation...\n")

    for(disease in diseases) {
        cat(paste("\n[INFO] Building Average Rank to Disease:", disease, "\n"))
        score_base_dir = here::here("src", "Data", "Drug", "Score", disease)
        
        pstep_df = read.csv(here::here(score_base_dir, "pstep_kernel.csv")) %>%
            dplyr::mutate(pstep_kernel = rank(dplyr::desc(score), ties.method = "average")) %>%
            dplyr::select(drugbank_id, pstep_kernel)

        reg_lap_df = read.csv(here::here(score_base_dir, "regularised_laplacian_kernel.csv")) %>%
            dplyr::mutate(regularised_laplacian_kernel = rank(dplyr::desc(score), ties.method = "average")) %>%
            dplyr::select(drugbank_id, regularised_laplacian_kernel)
        inv_cos_df = read.csv(here::here(score_base_dir, "inverse_cosine_kernel.csv")) %>%
           dplyr::mutate(inverse_cosine_kernel = rank(dplyr::desc(score), ties.method = "average"))  %>%
            dplyr::select(drugbank_id, inverse_cosine_kernel)

        commute_df = read.csv(here::here(score_base_dir, "commute_time_kernel.csv")) %>%
             dplyr::mutate(commute_time_kernel = rank(dplyr::desc(score), ties.method = "average"))  %>%
            dplyr::select(drugbank_id, commute_time_kernel)

        diffusion_df = read.csv(here::here(score_base_dir, "diffusion_kernel.csv")) %>%
             dplyr::mutate(diffusion_kernel = rank(dplyr::desc(score), ties.method = "average"))  %>%
            dplyr::select(drugbank_id, diffusion_kernel)

        merge_df = pstep_df %>%
            dplyr::left_join(reg_lap_df, by = "drugbank_id") %>%
            dplyr::left_join(commute_df, by = "drugbank_id") %>%
            dplyr::left_join(inv_cos_df, by = "drugbank_id") %>%
            dplyr::left_join(diffusion_df, by = "drugbank_id") 
        
        cat(sprintf(
            "[INFO] Merge complete | Rows: %d | Columns: %d\n",nrow(merge_df), ncol(merge_df)))

        final_df <- merge_df %>%
            dplyr::mutate(
                average_rank = rowMeans(
                    dplyr::select(., dplyr::ends_with("_kernel")), 
                    na.rm = TRUE
                )
            ) %>%
            dplyr::arrange(average_rank) %>%
            dplyr::select(drugbank_id,dplyr::ends_with("_kernel"),average_rank)
        output_file_name=here(score_base_dir,paste0("average_rank_", disease, ".csv"))        
        write.csv(final_df, file=output_file_name , row.names = FALSE)

        cat(sprintf("✔ [DONE] Disease: %s | Saved: %s\n", disease, output_file_name))
    }

}
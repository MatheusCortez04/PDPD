library(here)
library(purrr)
library(readr)
source(here("src","Utils","kernels.R"))
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
        run_drug_scoring()
        summarise_drug_max_score()
        generate_drug_rank()
    }
)

run_drug_scoring = function(){
    kernel_names = c(
        "diffusion_kernel","pstep_kernel","regularised_laplacian_kernel",
        "commute_time_kernel","inverse_cosine_kernel"
    )
    disease_data = list(
        MDD = list(name = "MDD", vector = build_protein_mdd_df()),
        BD = list(name = "BD", vector = build_protein_bipolar_df())
        
    )
    purrr::walk(kernel_names,function(kernel_name){
        kernel_path_rdata = here("src","Data","Kernels","RData", paste0(kernel_name, ".Rdata"))
        kernel_file_exists = file.exists(kernel_path_rdata)
        if(!kernel_file_exists){
            cat("Kernel file",kernel_name,"no exists\n")
            Sys.sleep(1)
            next()
        }
 
        kernel = load_rdata(kernel_path_rdata)
 
            purrr::walk(disease_data,score_genes_for_disease,
                 kernel = kernel,
                 kernel_name = kernel_name)
    })
}

generate_drug_rank = function() {
    diseases = c("MDD", "BD")

    for(disease in diseases) {
        base::cat(base::paste("[WARN] Building Average Rank to Disease:", disease, "\n"))
        base_dir = here::here("src", "Data", "Drug", "Score", disease)
        
        pstep_df = read.csv(here::here(base_dir, "pstep_kernel.csv")) %>%
            dplyr::mutate(rank_pstep = base::rank(dplyr::desc(max_gene_score), ties.method = "average")) %>%
            dplyr::select(drugbank_id, rank_pstep)

        reg_lap_df = read.csv(here::here(base_dir, "regularised_laplacian_kernel.csv")) %>%
            dplyr::mutate(rank_reg_lap = base::rank(dplyr::desc(max_gene_score), ties.method = "average")) %>%
            dplyr::select(drugbank_id, rank_reg_lap)
        inv_cos_df = read.csv(here::here(base_dir, "inverse_cosine_kernel.csv")) %>%
           dplyr::mutate(rank_inv_cos = base::rank(dplyr::desc(max_gene_score), ties.method = "average"))  %>%
            dplyr::select(drugbank_id, rank_inv_cos)

        commute_df = read.csv(here::here(base_dir, "commute_time_kernel.csv")) %>%
             dplyr::mutate(rank_commute = base::rank(dplyr::desc(max_gene_score), ties.method = "average"))  %>%
            dplyr::select(drugbank_id, rank_commute)

        diffusion_df = read.csv(here::here(base_dir, "diffusion_kernel.csv")) %>%
             dplyr::mutate(rank_diffusion = base::rank(dplyr::desc(max_gene_score), ties.method = "average"))  %>%
            dplyr::select(drugbank_id, rank_diffusion)

        merge_df = pstep_df %>%
            dplyr::left_join(reg_lap_df, by = "drugbank_id") %>%
            dplyr::left_join(commute_df, by = "drugbank_id") %>%
            dplyr::left_join(inv_cos_df, by = "drugbank_id") %>%
            dplyr::left_join(diffusion_df, by = "drugbank_id") 
            

        final_df <- merge_df %>%
            dplyr::mutate(
                average_rank = rowMeans(
                    dplyr::select(., dplyr::starts_with("rank_")), 
                    na.rm = TRUE
                )
            ) %>%
            dplyr::arrange(average_rank) %>%
            dplyr::select(drugbank_id,dplyr::starts_with("rank_"),average_rank)

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
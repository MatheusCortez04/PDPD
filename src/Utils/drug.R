library(here)
library(purrr)
library(readr)
source(here("src","Utils","kernels.R"))
scoring_drug_disease_menu = function(){
  while(TRUE){
    clear_console()
    cat("--- Drug Menu ---\n\n")
    cat(" [1] Generate Drug Score\n")
    cat(" [2] Generate MDD ROC \n")
    cat(" [3] Generate MDD Recall@K Graph \n")
    cat(" [4] Generate Bipolar ROC \n")
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
    },
    '2' = function(){
        generate_roc_curve_mdd()

    },
    '3'= function(){
        generate_recall_k_MDD()
    },
    '4'= function(){
        generate_roc_curve_bipolar()
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

        merge_df = pstep_df %>%
            dplyr::left_join(reg_lap_df, by = "drugbank_id") %>%
            dplyr::left_join(commute_df, by = "drugbank_id") %>%
            dplyr::left_join(inv_cos_df, by = "drugbank_id") 
            

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

generate_roc_curve_bipolar = function(){
    bipolar_repodb = read_tsv(here("src","Data","REPODB","BIPOLAR_REPODB.tsv"), show_col_types = FALSE)
    bipolar_repodb = bipolar_repodb %>% filter(status=="Approved")
    cat("Total de drogas valida pelo RepoDb: ", nrow(bipolar_repodb), "\n")
    bipolar_rank_file_path = here("src", "Data", "Drug", "Score", "BD","average_kernel_rank.csv")
    bipolar_average_rank = read.csv(bipolar_rank_file_path)

    pred_bipolar = bipolar_average_rank %>% 
        mutate(bipolar_repodb_validated = ifelse(drugbank_id %in% bipolar_repodb$drugbank_id, 1, 0))

    previstas_validas = pred_bipolar %>% filter(bipolar_repodb_validated==1)
    cat("Total de drogas previtas  validas: ", nrow(previstas_validas), "\n")
    pred_bipolar$bipolar_repodb_validated = factor(pred_bipolar$bipolar_repodb_validated,
                          levels = c(0,1),
                          labels = c("Not validated by repODB", "Validated by repODB"))

    if(nrow(previstas_validas)<=0){
        cat("Curva ROC indisponivel pois não foi previsto nenhuma droga:\n")
        Sys.sleep(1.5)
        return(NULL)
    }
    output_pdf_path = here("src", "Relatorios", "ROC", "bipolar_average_rank.pdf")
    table(pred_bipolar$bipolar_repodb_validated)
    dir.create(dirname(output_pdf_path), recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(output_pdf_path, width = 6, height = 6)
    roc_out = reportROC::reportROC(gold = pred_bipolar$bipolar_repodb_validated,
                predictor = -1*pred_bipolar$Mean_Rank,
                plot=FALSE)
    
    write.csv(pred_bipolar,file=here(output_dir,"prediction_bipolar.csv"),row.names=FALSE)
    plot(roc_out)
    grDevices::dev.off()
    message(sprintf("Salvando Curva ROC em: %s", output_pdf_path))
    Sys.sleep(1.5)
    return(roc_out)
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
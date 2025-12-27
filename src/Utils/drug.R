library(here)
library(purrr)
library(readr)

scoring_drug_disease_menu = function(){
  while(TRUE){
    clear_console()
    cat("--- Drug Menu ---\n\n")
    cat(" [1] Generate Drug Protein Matrix\n")
    cat(" [2] Generate Drug Score \n")
    cat(" [3] Generate MDD ROC \n")
    cat(" [4] Generate MDD Recall@K Graph \n")
    cat(" [5] Generate Bipolar ROC \n")
    cat(" [B] Back\n\n")

    input = readline(prompt = "Choice option: ")
    if (toupper(trimws(input)) == "B") break

     drug_function_mapper[[input]]()
  }
}

drug_function_mapper = list(
    '1' = function() {
    },
    '2' = function(){
        run_drug_scoring()
        summarise_drug_max_score()

    },
    '3'= function(){
        generate_roc_curve_mdd()
    },
    '4'= function(){
        generate_recall_k_MDD()
    },
    '5'= function(){
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
        dplyr::mutate(rank_pstep = base::rank(-max_gene_score, ties.method = "average")) %>%
        dplyr::select(drugbank_id, rank_pstep)

    reg_lap_df = read.csv(here::here(base_dir, "regularised_laplacian_kernel.csv")) %>%
        dplyr::mutate(rank_reg_lap = base::rank(-max_gene_score, ties.method = "average")) %>%
        dplyr::select(drugbank_id, rank_reg_lap)

    inv_cos_df = read.csv(here::here(base_dir, "inverse_cosine_kernel.csv")) %>%
        dplyr::mutate(rank_inv_cos = base::rank(-max_gene_score, ties.method = "average")) %>%
        dplyr::select(drugbank_id, rank_inv_cos)

    commute_df = read.csv(here::here(base_dir, "commute_time_kernel.csv")) %>%
        dplyr::mutate(rank_commute = base::rank(-max_gene_score, ties.method = "average")) %>%
        dplyr::select(drugbank_id, rank_commute)

    merge_df = pstep_df %>%
        dplyr::left_join(reg_lap_df, by = "drugbank_id") %>%
        dplyr::left_join(inv_cos_df, by = "drugbank_id") %>%
        dplyr::left_join(commute_df, by = "drugbank_id")

    final_df <- merge_df %>%
        dplyr::mutate(
            average_rank = rowMeans(
                dplyr::select(., dplyr::starts_with("rank_")), 
                na.rm = TRUE
            )
        ) %>%
        dplyr::arrange(average_rank) %>%
        dplyr::select(drugbank_id,rank_pstep,rank_reg_lap,rank_inv_cos,rank_commute,average_rank)

    write.csv(final_df, here(base_dir,paste0("average_rank_", disease, ".csv")), row.names = FALSE)
}

}


score_genes_for_disease = function(kernel,disease_info,kernel_name){
    compute_gene_scores_from_kernel(kernel,kernel_name,disease_info)
}
generate_roc_curve_mdd = function(){

    drug_target_df  = read.csv(here("src","Data","Drug","drug_targets_DrugBank_Gysi.csv"), sep=",")
    mdd_repodb = read_tsv(here("src","Data","REPODB","MDD_REPODB.tsv"), show_col_types = FALSE)
     mdd_repodb =  dplyr::semi_join(mdd_repodb,drug_target_df,by='drugbank_id')%>% 
        filter(status=="Approved")
        
    cat("Total de drogas validas pelo RepoDb: ", nrow(mdd_repodb), "\n")

    mdd_rank_file_path = here("src", "Data", "Drug", "Score", "MDD","average_kernel_rank.csv")
    mdd_average_rank = read.csv(mdd_rank_file_path)

    pred_mdd = mdd_average_rank %>% 
        mutate(mdd_repodb_validated = ifelse(drugbank_id %in% mdd_repodb$drugbank_id, 1, 0))
    previstas_validas = pred_mdd %>% filter(mdd_repodb_validated==1)
    cat("Total de drogas previtas  validas: ", nrow(previstas_validas), "\n")

    pred_mdd$mdd_repodb_validated = factor(pred_mdd$mdd_repodb_validated,
                          levels = c(0,1),
                          labels = c("Not validated by repODB", "Validated by repODB"))

    if(nrow(previstas_validas)<=0){
        cat("Curva ROC indisponivel pois não foi previsto nenhuma droga:\n")
        Sys.sleep(1.5)
        return(NULL)
    }
    output_dir = here("src", "Relatorios", "ROC", "MDD")
    output_graph_file_name = "mdd_roc.pdf" 
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(pred_mdd,file=here(output_dir,"prediction_mdd.csv"),row.names=FALSE)

    grDevices::pdf(here(output_dir,output_graph_file_name), width = 6, height = 6)
    roc_out = reportROC::reportROC(gold = pred_mdd$mdd_repodb_validated,
                                    predictor = -1*pred_mdd$Mean_Rank,
                                    plot=TRUE
                                    )
 
    grDevices::dev.off()
    message(sprintf("Salvando Curva ROC em: %s", here(output_dir,output_graph_file_name)))
    Sys.sleep(1.5)
    return(roc_out)
}
# COM A VALIDACAO ATUAL NÃO HÁ NENHUMA DROGA APROVADA PARA BIPOLARIDADE
# O QUE GERA UM PDF VAZIO, VERIFICAR 
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
generate_recall_k_MDD = function(){

    pred_mdd_rank_file_path = here("src", "Relatorios", "ROC", "MDD","prediction_mdd.csv")

        if (!file.exists(pred_mdd_rank_file_path)) {
            cat("[WARN] Drug Score :", pred_mdd_rank_file_path, "\n")
            next()
        }
    pred_mdd_rank = read.csv(pred_mdd_rank_file_path)
    valid_true <- pred_mdd_rank %>%
    filter(mdd_repodb_validated == "Validated by repODB") %>%
    pull(drugbank_id)

    pred_mdd_rank = pred_mdd_rank %>%
        select(drugbank_id,Mean_Rank)


    # Atribui valor 1 a toda droga prevista contida no vetor de drogas validas 
    valid_prediction = ifelse(pred_mdd_rank$drugbank_id %in% valid_true, 1, 0)

    #cumSum realiza a soma cumulariva e divide pelo total de verdadeiros positivos
    # recall = VP/VP+FN(neste caso nao tem FN a nao ser que seja inserido um valor de corte no rank)
    recall_values = cumsum(valid_prediction)/length(valid_true)

    recall_df= data.frame(K = 1:nrow(pred_mdd_rank),
                        Recall = recall_values)

        recall_df$highlight <- ifelse(recall_df$K %% 50 == 0, TRUE, FALSE)
        g = ggplot(recall_df, aes(x = K, y = Recall)) +
        geom_line(size = 1) +
        geom_point(
        data = subset(recall_df, highlight == TRUE),
        size = 3,
        color = "red"
        ) +
    geom_text(
    data = subset(recall_df, highlight == TRUE),
    aes(label = paste0("(", K, ", ", round(Recall, 3), ")")),
    vjust = -0.7,
    size = 3,
    check_overlap=TRUE
    )+
        labs(
            title = "Recall@K para MDD",
            x = "K (Top-K)",
            y = "Recall"
        ) +
        theme_minimal(base_size = 14)
    print(g)


    pdf_path <- here("src", "Relatorios", "ROC", "MDD", "recall_at_k_MDD.pdf")
    ggsave(pdf_path, plot = g, width = 8, height = 6)

}
compute_gene_scores_from_kernel= function(kernel,kernel_name,disease_info){
    disease_name = disease_info$name
    disease_vector = disease_info$vector

    score_vector =  kernel %*% disease_vector$is_disease
    protein_scores_df = data.frame(
        entrez_id=  rownames(kernel),
        gene_score = score_vector
    )
    output_dir = here("src","Data","Kernels","Score",disease_name)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_path_csv = here(output_dir, paste0(kernel_name, ".csv"))
    write.csv(protein_scores_df, file = output_path_csv, row.names = FALSE)
    cat("\n✔ Gene Score saved in:", output_path_csv, "\n")
}


load_drug_target_df = function(){
    drug_target_file_path = here("src","Data","Drug","drug_targets_DrugBank_Gysi.csv")
    drug_target_file_already_exists = file.exists(drug_target_file_path)
    if(!drug_target_file_already_exists){
        cat("\n[Error]: Required file Drug target file not found.\n This file is necessary to calculate the score. Would you like to create it now\n")
        Sys.sleep(1.5)
        return()
    }
    drug_target_df = read.csv(drug_target_file_path)
    drug_target_df %>% dplyr::distinct(drugbank_id, entrez_id) %>%
        dplyr::mutate(entrez_id = as.character(entrez_id))
    invisible(drug_target_df)
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
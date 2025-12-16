library(here)
library(purrr)
library(readr)
source(here("src","Utils","utils.R"))
generate_ordered_drug_protein_matrix = function(drug_target_df,protein_nodes,drug_nodes){
  cat("\n--- Generating Drug Target Matrix ---\n")
    output_file_name="drug_target_matrix"
    drug_dir =here("src","Data","Drug")
    drug_Rdata_dir =here("src","Data","Drug","RData")
    create_dir(drug_dir)
    create_dir(drug_Rdata_dir)
    
    output_path_csv = here(drug_dir,paste0(output_file_name, ".csv"))
    output_path_rdata = here(drug_Rdata_dir,paste0(output_file_name, ".Rdata"))

    matrix_drug_protein = build_drug_protein_matrix(drug_target_df,protein_nodes,drug_nodes)

    write.csv(matrix_drug_protein, file =output_path_csv,row.names=FALSE)
    cat("Csv file  saved to:", output_path_csv, "\n")
    Sys.sleep(1)
    matrix_drug_protein_numeric <- as.matrix(matrix_drug_protein[ , -1])
    rownames(matrix_drug_protein_numeric) <- matrix_drug_protein$drugbank_id
    storage.mode(matrix_drug_protein_numeric) <- "numeric"

    save(matrix_drug_protein_numeric, file = output_path_rdata)
    cat("R object 'drug_target_df' saved to:", output_path_rdata, "\n")
    Sys.sleep(1)
}

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
        cat("\n--- Generating and  Drug  Target Matrix  ---\n")
        drug_target_df  = read.csv(here("src","Data","Drug","drug_targets_DrugBank_Gysi.csv"), sep=",")

        cat("Drug Target file reading complete.\n")
        ppi_df  = read.csv(here("src","Data","PPI_gysi.csv"), sep=",")
        cat("PPI  file reading complete.\n")

        protein_nodes = get_ppi_nodes(ppi_df)
        drug_nodes= get_drug_nodes(drug_target_df)
        generate_ordered_drug_protein_matrix(drug_target_df,protein_nodes,drug_nodes)
    },
    '2' = function(){
        calculate_drug_score()
        generate_drug_rank("MDD")
        generate_drug_rank("BD")
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

calculate_drug_score = function(){
    kernel_names = c(
        "diffusion_kernel","pstep_kernel","regularised_laplacian_kernel",
        "commute_time_kernel","inverse_cosine_kernel"
    )
    drug_protein_matrix = load_drug_target_matrix()
    disease_data <- list(
        MDD = list(name = "MDD", vector = build_protein_mdd_df()),
        BD  = list(name = "BD",  vector = build_protein_bipolar_df())
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
        list_data = list(MDD = disease_data$MDD, BD = disease_data$BD)
        list_data %>%
            purrr::walk(process_scores_for_disease, 
                 kernel = kernel, 
                 drug_protein_matrix = drug_protein_matrix, 
                 kernel_name = kernel_name)
    })
} 

generate_drug_rank = function(disease = c("MDD", "BD")) {
    disease = match.arg(disease)

    cat("[WARN] Building Average Rank to Disease:", disease, "\n")

    kernel_names = c(
        "diffusion_kernel","pstep_kernel","regularised_laplacian_kernel",
        "commute_time_kernel","inverse_cosine_kernel"
    )

    score_base_dir = here("src", "Data", "Drug", "Score", disease)
    rdata_dir = here(score_base_dir, "RData")

    score_list = kernel_names %>%
        map(function(kernel){
            rdata_path = here(score_base_dir, "RData", paste0(kernel, ".Rdata"))
            if (!file.exists(rdata_path)) {
                cat("[WARN] Score RData not  found:", rdata_path, "\n")
                next()
        }
            df = load_rdata(rdata_path)
            df = df %>%
                mutate(!!kernel := rank(-drug_score,ties.method = "average")) %>%
                select(drugbank_id, all_of(kernel))
                return(df)
        }) %>% compact()



    final_rank_df = reduce(score_list, full_join, by = "drugbank_id")
    final_rank_df = final_rank_df %>%
        mutate(Mean_Rank = rowMeans(select(., all_of(kernel_names)), na.rm = TRUE)) %>%
        arrange(Mean_Rank)

    output_csv = here(score_base_dir, "average_kernel_rank.csv")
    write_csv(final_rank_df, output_csv)

    cat("\n✔ Average Kernel Rank saved in:", output_csv, "\n")
    Sys.sleep(1.5)
    return(final_rank_df)
}

build_drug_protein_matrix= function(data_frame,protein_nodes,drug_nodes){
  data_frame %>%
    filter(entrez_id %in% protein_nodes) %>%
    distinct(drugbank_id, entrez_id) %>%
    #O Value é utilizado para que cada iteração do dataframe original contenha o valor 1
    # e as demais lacunas sejam preenchidas com zero
    mutate(
      drugbank_id = factor(drugbank_id, levels = drug_nodes),
      entrez_id = factor(entrez_id, levels = protein_nodes),
      value = 1
    ) %>%
    pivot_wider(
      id_cols = drugbank_id,
      names_from = entrez_id,
      values_from = value,
      values_fill = 0,
      names_expand = TRUE,
      id_expand = TRUE
    ) %>%
    arrange(drugbank_id)
}

load_drug_target_matrix = function() {
    drug_target_rdata_file_path = here("src","Data","Drug","RData","drug_target_matrix.Rdata")
    drug_protein_file_already_exists = file.exists(drug_target_rdata_file_path)
    if(!drug_protein_file_already_exists){
        cat("\n[Error]: Required file Drug target matrix not found.\n This file is necessary to calculate the score. Would you like to create it now\n")
        Sys.sleep(1.5)
        return()
    }
  invisible(load_rdata(drug_target_rdata_file_path))
}

process_scores_for_disease = function(kernel,disease_info,drug_protein_matrix,kernel_name){
    
    disease_name <- disease_info$name
    disease_vector <- disease_info$vector

    cat(paste("\n--- Calculando", disease_name, "score para", kernel_name, "---\n"))
    score_vector = drug_protein_matrix %*% kernel %*% disease_vector$is_disease

    score_df = data.frame(
        drugbank_id = rownames(drug_protein_matrix),
        drug_score = as.numeric(score_vector)
    ) %>%
        arrange(desc(drug_score))

    output_dir <- here("src", "Data", "Drug", "Score", disease_name)
    rdata_sub_dir <- here(output_dir, "RData")

    dir.create(rdata_sub_dir, recursive = TRUE, showWarnings = FALSE)
    output_path_csv <- here(output_dir, paste0(kernel_name, ".csv"))
    output_path_rdata <- here(rdata_sub_dir, paste0(kernel_name, ".Rdata"))

    write.csv(score_df, file = output_path_csv, row.names = FALSE)

    save(score_df, file = output_path_rdata)
    cat(paste("   ✔ Score salvo em:", output_path_csv, "\n"))
    invisible(score_df)
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

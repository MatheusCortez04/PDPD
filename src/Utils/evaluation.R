source(here("src","Utils","drug.R"))
source(here("src","Utils","utils.R"))
source(here("src","Utils","open_targets_requests.R"))
library(tidyr)
library(purrr)
evaluation_menu = function(){
  while(TRUE){
    clear_console()
    cat("--- Evaluation Menu ---\n\n")
    cat(" [1] Generate MDD ROC \n")
    cat(" [2] Generate MDD Recall@K \n")
    cat(" [3] Generate BD ROC \n")
    cat(" [4] Generate BD Recall@K \n")
    cat(" [5] Generate Top n drug rank to MDD \n")
    cat(" [6] Generate Top n drug rank to BD \n")
    cat(" [B] Back\n\n")

    input = readline(prompt = "Choice option: ")
    if (toupper(trimws(input)) == "B") break
     evaluation_function_mapper[[input]]()
  }
}


evaluation_function_mapper = list(
    '1' = function() {
        generate_roc_curve("MDD")
    },
    '2' = function(){
         generate_recall_k("MDD")
    },
    '3'= function(){
       generate_roc_curve("BD")
    },
    '4'= function(){
        generate_recall_k("BD")
    },
    '5'= function(){
        n = readline(prompt = "Enter n value  to generate drug rank (Default 10): ")
        create_top_drugs_file("MDD", as.integer(n))
    },
    '6'= function(){
        n = readline(prompt = "Enter n value  to generate drug rank (Default 10): ")
        create_top_drugs_file("BD", as.integer(n))
    }
)

create_top_drugs_file= function(disease = c("MDD", "BD"),n=10){
    disease = match.arg(disease)

    output_dir = here("src", "Evaluation",disease)
    file_suffixes   = c("MDD" = "mdd", "BD" = "bipolar")
    file_suffix   = file_suffixes[disease]
    score_filter=get_score_disease_gene_association()
    prediction_disease_file_path = here(output_dir, paste0("prediction_", file_suffix,"_score_filter_",score_filter,"_.csv"))

    if(file.exists(prediction_disease_file_path)) {
        prediction_data = read.csv(prediction_disease_file_path) 
    }
    if(!file.exists(prediction_disease_file_path)) {
        cat("[WARN] Processed data not found. Running preparation first...\n")
        prediction_data = create_prediction_disease_info(disease)
    }
    drug_target_mapping  = import_drug_targets_df()
    ppi_gene_nodes= extract_ppi_genes()

    top_n_drugs= prediction_data %>% dplyr::slice_head(n=n) %>% 
        dplyr::select(drugbank_id,validation_status) %>%
        dplyr::rowwise()%>%
        dplyr::mutate(target_count = length(extract_targets_per_drug(drugbank_id,drug_target_mapping))) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(api_response = purrr::map(drugbank_id, get_drug_info_from_dbid)) %>%
        tidyr::unnest_wider(api_response) %>%
        select(drugbank_id,chembl_id,drug_name,target_count,validation_status)

        cat(sprintf("[INFO] Searching Drug Candidates by OpenTargets to  %s...\n", disease))
        drug_candidates <- get_drug_candidates(disease)
        top_n_joined = top_n_drugs %>% 
        dplyr::left_join(
            drug_candidates, 
            by = c("chembl_id" = "drug_id"),
        )

        approval_drugs = top_n_joined %>% filter(max_clinical_stage_open_targets=='APPROVAL')
        approval_drugs_final = approval_drugs %>%
            mutate(report_info = purrr::map(clinical_report_id, get_clinical_report_data)) %>%
                tidyr::unnest(report_info, names_sep = "_") %>% select(clinical_report_id,report_info_source,report_info_url,report_info_evidence_summary)
        glimpse(approval_drugs_final)
        top_n_joined = top_n_joined %>% dplyr::left_join(
            approval_drugs_final, 
            by = 'clinical_report_id',
        )
   

        
    output_file_name= here(output_dir,paste0("top_",n,"_drugs_",disease,"_score_filter_",score_filter,"_.csv"))

    write.csv(top_n_joined,output_file_name,row.names = FALSE)
    message(sprintf("[SUCCESS] Top Drug file  saved at: %s",output_file_name))

}

create_prediction_disease_info = function(disease = c("MDD", "BD")) {
    disease = match.arg(disease)

    repodb_prefixes = c("MDD" = "MDD", "BD" = "BIPOLAR")
    file_suffixes   = c("MDD" = "mdd", "BD" = "bipolar")
    
    repodb_prefix = repodb_prefixes[disease]
    file_suffix   = file_suffixes[disease]

    rank_file_path   = here("src", "Data", "Drug", "Score", disease, "average_kernel_rank.csv")
    repodb_file_path = here("src", "Data", "REPODB", paste0(repodb_prefix, "_REPODB.tsv"))
    output_dir       = here("src", "Evaluation", disease)
    score_filter = get_score_disease_gene_association()
    output_file_path = here(output_dir, paste0("prediction_", file_suffix,"_score_filter_",score_filter,"_.csv"))

    if (!file.exists(rank_file_path)) {
        cat("[WARN] Average Rank file not found. Building first...\n")
        drug_function_mapper[['1']]()
        drug_function_mapper[['2']]()
    }

    drug_target_mapping = import_drug_targets_df()
    gold_standard = read_tsv(repodb_file_path, show_col_types = FALSE)    
    gold_standard = gold_standard %>%  
        dplyr::semi_join(drug_target_mapping, by = 'drugbank_id')
 
    cat(sprintf("[INFO] Total valid drugs in RepoDB (Gold Standard): %d\n", nrow(gold_standard)))
    
    prediction_info = read.csv(rank_file_path)
    processed_predictions = prediction_info %>% 
        dplyr::mutate(
            validation_label = ifelse(drugbank_id %in% gold_standard$drugbank_id, 1, 0),
            validation_status = factor(validation_label,
                                       levels = c(0, 1),
                                       labels = c("Not validated by RepoDB", "Validated by RepoDB"))
        )

    validated_hits = processed_predictions %>% dplyr::filter(validation_label == 1)
    cat(sprintf("[INFO] Total predicted hits validated: %d\n", nrow(validated_hits)))
    
    if (nrow(validated_hits) <= 0) {
         cat("[WARN] ROC curve unavailable: No validated drugs were predicted in the ranking.\n")
         Sys.sleep(1.5)
         return(NULL)
    }

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(processed_predictions, file = output_file_path, row.names = FALSE)
    message(sprintf("[SUCCESS] prediction file saved at: %s", output_file_path))
    
    
    invisible(processed_predictions)
}

generate_roc_curve = function(disease = c("MDD", "BD")) {
    disease = match.arg(disease)

    file_suffixes = c("MDD" = "mdd", "BD" = "bipolar")
    file_suffix = file_suffixes[disease]

    output_dir = here("src", "Evaluation", disease)
    score_filter = get_score_disease_gene_association()
    
    prediction_file_path = here(output_dir, paste0("prediction_", file_suffix, "_score_filter_", score_filter, "_.csv"))

    if (file.exists(prediction_file_path)) {
        prediction_data = read.csv(prediction_file_path)
    }
    
    if (!file.exists(prediction_file_path)) {
        cat("[WARN] Processed data not found. Running preparation first...\n")
        prediction_data = create_prediction_disease_info(disease)
    }

    if (sum(prediction_data$validation_label) <= 0) {
        cat(sprintf("[WARN] Cannot plot ROC for %s: No hits found in the ranking.\n", disease))
        return(NULL)
    }
   
    output_roc_dir = here(output_dir, "ROC")
    dir.create(output_roc_dir, recursive = TRUE, showWarnings = FALSE)
    output_graph_file_name = paste0(file_suffix,"_score_filter_",score_filter,"_roc_curve.pdf")
    output_graph_path = here(output_roc_dir, output_graph_file_name)
    
    grDevices::pdf(output_graph_path, width = 6, height = 6)
    
    roc_results = reportROC::reportROC(
         gold = prediction_data$validation_label,
         predictor = -1 * prediction_data$average_rank,
         plot = TRUE
    )
    
    grDevices::dev.off()
    
    message(sprintf("[SUCCESS] ROC curve saved at: %s", output_graph_path))
    Sys.sleep(1.5)
    
    invisible(roc_results)
}

generate_recall_k = function(disease = c("MDD", "BD")) {
    disease = match.arg(disease)

    file_suffixes = c("MDD" = "mdd", "BD" = "bipolar")

    score_filter = get_score_disease_gene_association()
    plot_titles = c("MDD" = paste0("Recall@K - Transtorno Depressivo Maior GDA Score >= ",score_filter), 
                    "BD" = paste0("Recall@K - Transtorno Bipolar GDA Score >= ",score_filter))

    file_suffix = file_suffixes[disease]
    plot_title = plot_titles[disease]

    output_dir = here("src", "Evaluation", disease)
    prediction_file_path = here(output_dir, paste0("prediction_", file_suffix, ".csv"))

    if (file.exists(prediction_file_path)) {
        prediction_data = read.csv(prediction_file_path) 
    }
    if (!file.exists(prediction_file_path)) {
        cat("[WARN] Processed data not found. Running preparation first...\n")
        prediction_data = create_prediction_disease_info(disease)
    }

    total_hits_count = sum(prediction_data$validation_label)
    
    if (total_hits_count == 0) {
        cat(sprintf("[WARN] No hits found for recall calculation for %s.\n", disease))
        return(NULL)
    }
    
    recall_values = cumsum(prediction_data$validation_label) / total_hits_count

    recall_df = data.frame(
        K = 1:nrow(prediction_data),
        Recall = recall_values
    )

    recall_df$highlight = ifelse(recall_df$K %% 50 == 0, TRUE, FALSE)

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
            check_overlap = TRUE
        ) +
        labs(
            title = plot_title,
            x = "K (Top-K)",
            y = "Recall"
        ) +
        theme_minimal(base_size = 14)

    output_recall_dir = here(output_dir, "Recall")
    dir.create(output_recall_dir, recursive = TRUE, showWarnings = FALSE)
    
    pdf_filename = paste0("recall_at_k_", file_suffix, "_score_filter_",score_filter,"_.pdf")
    pdf_path = here(output_recall_dir, pdf_filename)
    
    ggsave(pdf_path, plot = g, width = 8, height = 6)
    message(sprintf("[SUCCESS] Recall@K graph saved at: %s", pdf_path))
    
    Sys.sleep(1.5)
    invisible(recall_df)
}


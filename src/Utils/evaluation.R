library(here)
library(dplyr)
library(reportROC)
source(here("src","Utils","utils.R"))
source(here("src","Utils","drug.R"))

evaluation_menu = function(){
  while(TRUE){
    clear_console()
    cat("--- Evaluation Menu ---\n\n")
    cat(" [1] Generate MDD ROC \n")
    cat(" [2] Generate MDD Recall@K \n")
    cat(" [3] Generate BD ROC \n")
    cat(" [4] Generate BD Recall@K \n")
    cat(" [5] Generate Kernel Roc \n")
    cat(" [B] Back\n\n")

    input = readline(prompt = "Choice option: ")
    if (toupper(trimws(input)) == "B") break

     evaluation_function_mapper[[input]]()
  }
}


evaluation_function_mapper = list(
    '1' = function() {
        generate_roc_curve_mdd()
        create_top_drugs_file("MDD",20)
    },
    '2' = function(){
         generate_recall_k_MDD()
         create_top_drugs_file("MDD",20)
    },
    '3'= function(){
       generate_roc_curve_bipolar()
    },
    '4'= function(){
       generate_recall_k_BD()
    },
    '5'= function(){
       generate_roc_to_kernel()
    }
)

created_prediction_bipolar_data = function(){

    drug_target_mapping = read_drug_targets()

    bipolar_gold_standard = read_tsv(
        here("src","Data","REPODB","BD_REPODB.tsv"),
        show_col_types = FALSE
    )    

    bipolar_gold_standard = bipolar_gold_standard %>%  
        dplyr::semi_join(drug_target_mapping, by = 'drugbank_id') 

    cat(sprintf("[INFO] Total valid drugs in RepoDB (Gold Standard): %d\n", nrow(bipolar_gold_standard)))

    rank_file_path = here("src", "Data", "Drug", "Score","BD", "average_rank_BD.csv")
    
    if (!file.exists(rank_file_path)) {
        cat("[WARN] Average Rank file not found. Building first...\n")
        drug_function_mapper[[1]]()
    }

    bipolar_prediction = read.csv(rank_file_path, stringsAsFactors = FALSE)

    processed_predictions = bipolar_prediction %>% 
        dplyr::mutate(
            validation_label = ifelse(drugbank_id %in% bipolar_gold_standard$drugbank_id, 1, 0),
            validation_status = factor(
                validation_label,
                levels = c(0, 1),
                labels = c("Not validated by RepoDB", "Validated by RepoDB")
            )
        ) %>%
        dplyr::arrange(average_rank)

    validated_hits = processed_predictions %>% filter(validation_label == 1)

    cat(sprintf("[INFO] Total predicted hits validated: %d\n", nrow(validated_hits)))

    if (nrow(validated_hits) <= 0) {
        cat("[WARN] ROC curve unavailable: No validated drugs were predicted.\n")
        return(NULL)
    }

    output_dir = here("src", "Evaluation", "BD")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    write.csv(
        processed_predictions,
        file = here(output_dir, "prediction_bipolar.csv"),
        row.names = FALSE
    )
    print(summary(processed_predictions$average_rank))
    print(head(processed_predictions[processed_predictions$validation_label == 1, ], 10))

    validated_hits <- processed_predictions %>% filter(validation_label == 1)

    print(head(validated_hits, 10))
    message(sprintf("[SUCCESS] prediction file saved at: %s",
        here(output_dir, "prediction_mdd.csv")
    ))

    invisible(processed_predictions)
}
generate_roc_curve_mdd = function(){

    output_dir = here("src", "Evaluation", "MDD")
    score_filter= get_score_disease_gene_association()
    prediction_mdd_file_path = here(output_dir, paste0("prediction_mdd_score_filter_",score_filter,"_.csv"))
    is_file_exists = file.exists(prediction_mdd_file_path)
    if (!is_file_exists) {
        cat("[WARN] Processed data not found. Running preparation first...\n")
        prediction_mdd = create_prediction_disease_info("MDD")
    } 
    if(is_file_exists) {
        prediction_mdd = read.csv(prediction_mdd_file_path)
    }

    if (sum(prediction_mdd$validation_label) <= 0) {
        cat("[WARN] Cannot plot ROC: No validated drugs found.\n")
        return(NULL)
    }

    output_roc_dir = here(output_dir, "ROC")
    dir.create(output_roc_dir, recursive = TRUE, showWarnings = FALSE)

    output_graph_file_name = paste0("mdd_score_filter_",score_filter,"_roc_curve.pdf")
    output_path = here(output_roc_dir, output_graph_file_name)

    grDevices::pdf(output_path, width = 6, height = 6)
    
    roc_results = reportROC::reportROC(
        gold = prediction_mdd$validation_label,
        predictor = -1 * prediction_mdd$average_rank,
        plot = TRUE
    )

    grDevices::dev.off()

    message(sprintf("[SUCCESS] ROC curve saved at: %s", output_path))
}

generate_roc_curve_bipolar = function(){

    output_dir = here("src", "Evaluation", "BD")
    prediction_bipolar_file_path = here(output_dir, "prediction_bipolar.csv")

    if (!file.exists(prediction_bipolar_file_path)) {
        cat("[WARN] Processed data not found. Running preparation first...\n")
        prediction_bipolar = created_prediction_bipolar_data()
    } else {
        prediction_bipolar = read.csv(prediction_bipolar_file_path)
    }

    if (sum(prediction_bipolar$validation_label) <= 0) {
        cat("[WARN] Cannot plot ROC: No validated drugs found.\n")
        return(NULL)
    }

    output_roc_dir = here(output_dir, "ROC")
    dir.create(output_roc_dir, recursive = TRUE, showWarnings = FALSE)

    output_graph_file_name = "bipolar_roc_curve.pdf"
    output_path = here(output_roc_dir, output_graph_file_name)

    grDevices::pdf(output_path, width = 6, height = 6)
    
    roc_results = reportROC::reportROC(
        gold = prediction_bipolar$validation_label,
        predictor = -1 * prediction_bipolar$average_rank,
        plot = TRUE
    )

    grDevices::dev.off()

    message(sprintf("[SUCCESS] ROC curve saved at: %s", output_path))
}
generate_recall_k_MDD = function(){
    output_dir = here("src", "Evaluation","MDD")
    score_filter=get_score_disease_gene_association()
    prediction_mdd_file_path =  here(output_dir, paste0("prediction_mdd_score_filter_",score_filter,"_.csv"))
    is_file_exists = file.exists(prediction_mdd_file_path)
    if(!is_file_exists){
        cat("[WARN] Processed data not found. Running preparation first...\n")
        prediction_data = create_prediction_disease_info("MDD")
    }
    if(is_file_exists){
        prediction_data = read.csv(prediction_mdd_file_path) 
    }

    total_hits_count = sum(prediction_data$validation_label)
    if (total_hits_count == 0) {
        cat("[WARN] No hits found for recall calculation.\n")
        return(NULL)
    }
    recall_values = cumsum(prediction_data$validation_label) / total_hits_count

    # #cumSum realiza a soma cumulariva e divide pelo total de verdadeiros positivos
    # # recall = VP/VP+FN(neste caso nao tem FN a nao ser que seja inserido um valor de corte no rank)

    recall_df = data.frame(
        K = 1:nrow(prediction_data),
        Recall = recall_values
    )

    recall_df$highlight <- ifelse(recall_df$K %% 50 == 0, TRUE, FALSE)
        g = ggplot(recall_df, aes(x = K, y = Recall)) +
        geom_line(size = 1) +
        geom_point(
            data = subset(recall_df, highlight == TRUE),
            size = 3,
            color = "red") +
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
        y = "Recall") +
    theme_minimal(base_size = 14)

    output_recall_dir= here(output_dir,"Recall")
    dir.create(output_recall_dir, recursive = TRUE, showWarnings = FALSE)
    pdf_path = here(output_recall_dir,"recall_at_k_MDD.pdf")
    ggsave(pdf_path, plot = g, width = 8, height = 6)

}
generate_recall_k_BD = function(){
    output_dir = here("src", "Evaluation","BD")
    prediction_bipolar_file_path = here(output_dir,"prediction_bipolar.csv")

    if(!file.exists(prediction_bipolar_file_path)){
        cat("[WARN] Processed data not found. Running preparation first...\n")
        prediction_data = created_prediction_bipolar_data()
    }
    else{
        prediction_data = read.csv(prediction_bipolar_file_path) 
    }

    total_hits_count = sum(prediction_data$validation_label)
    if (total_hits_count == 0) {
        cat("[WARN] No hits found for recall calculation.\n")
        return(NULL)
    }
    recall_values = cumsum(prediction_data$validation_label) / total_hits_count

    # #cumSum realiza a soma cumulariva e divide pelo total de verdadeiros positivos
    # # recall = VP/VP+FN(neste caso nao tem FN a nao ser que seja inserido um valor de corte no rank)

    recall_df = data.frame(
        K = 1:nrow(prediction_data),
        Recall = recall_values
    )

    recall_df$highlight <- ifelse(recall_df$K %% 50 == 0, TRUE, FALSE)
        g = ggplot(recall_df, aes(x = K, y = Recall)) +
        geom_line(size = 1) +
        geom_point(
            data = subset(recall_df, highlight == TRUE),
            size = 3,
            color = "red") +
    geom_text(
        data = subset(recall_df, highlight == TRUE),
        aes(label = paste0("(", K, ", ", round(Recall, 3), ")")),
        vjust = -0.7,
        size = 3,
        check_overlap=TRUE
    )+
    labs(
        title = "Recall@K Bipolar Disease",
        x = "K (Top-K)",
        y = "Recall") +
    theme_minimal(base_size = 14)

    output_recall_dir= here(output_dir,"Recall")
    dir.create(output_recall_dir, recursive = TRUE, showWarnings = FALSE)
    pdf_path = here(output_recall_dir,"recall_at_k_bd.pdf")
    ggsave(pdf_path, plot = g, width = 8, height = 6)

}


generate_roc_to_kernel = function(){
    diseases = c("MDD","BD")
    drug_score_dir =here("src","Data","Drug","Score")
     
    
    kernel_names = c(
        "diffusion_kernel","pstep_kernel","regularised_laplacian_kernel",
        "commute_time_kernel","inverse_cosine_kernel"
    )
    for(disease in diseases){
        gold_standard_file_path = here("src","Data","REPODB",paste0(disease,"_REPODB.tsv"))
        if (!file.exists(gold_standard_file_path)) {
            cat("[WARN] Drug Score to kernel:",kernel,"and disease:",disease,"\n")
            next()
            }
        gold_standard = read_tsv(gold_standard_file_path, show_col_types = FALSE)   
        for(kernel in kernel_names){
            input_file_path =here(drug_score_dir,disease,paste0(kernel,".csv")) 
            if (!file.exists(input_file_path)) {
                cat("[WARN] not found gold standard file to disease:",disease,"\n")
                next()
            }
            drug_score_kernel = read.csv(input_file_path)
            gold_standard = gold_standard %>%  dplyr::semi_join(drug_score_kernel,by='drugbank_id')

            kernel_predictions = drug_score_kernel %>% 
                dplyr::mutate(
                    validation_label =  ifelse(drugbank_id %in% gold_standard$drugbank_id, 1, 0),
                    validation_status = factor(validation_label,
                                             levels = c(0, 1),
                                             labels = c("Not validated by RepoDB", "Validated by RepoDB"))
            )

                validated_hits = kernel_predictions %>% filter(validation_label==1)
                cat(sprintf("[INFO] Total predicted hits validated: %d\n",nrow(validated_hits)))
                if (nrow(validated_hits) <= 0) {
                    cat("[WARN] ROC curve unavailable: No validated drugs were predicted in the ranking.\n")
                    Sys.sleep(1.5)
                    return(NULL)
                }
                output_roc_dir = here("src","Evaluation",disease,"ROC","KERNEL")
                output_graph_file_name = paste0(kernel,"_roc_curve.pdf")
                dir.create(output_roc_dir, recursive = TRUE, showWarnings = FALSE)

                grDevices::pdf(here(output_roc_dir,output_graph_file_name), width = 6, height = 6)
                roc_results = reportROC::reportROC(
                    gold = kernel_predictions$validation_label,
                    predictor =kernel_predictions$score,
                    plot = TRUE)

                
                grDevices::dev.off()
                message(sprintf("[SUCCESS] ROC curve saved at: %s",here(output_roc_dir, output_graph_file_name)))
                Sys.sleep(1.5)



        }
    }



}
create_top_drugs_file= function(disease = c("MDD", "BD"),n=10){
    disease = match.arg(disease)

    output_dir = here("src", "Evaluation",disease)
    file_suffixes = c("MDD" = "mdd", "BD" = "bipolar")
    file_suffix= file_suffixes[disease]
    score_filter=get_score_disease_gene_association()
    prediction_disease_file_path = here(output_dir, paste0("prediction_", file_suffix,"_score_filter_",score_filter,"_.csv"))
    if(file.exists(prediction_disease_file_path)) {
        prediction_data = read.csv(prediction_disease_file_path) 
    }
    if(!file.exists(prediction_disease_file_path)) {
        cat("[WARN] Processed data not found. Running preparation first...\n")
        prediction_data = create_prediction_disease_info(disease)
    }

    drug_target_mapping  = read_drug_targets()
    ppi_gene_nodes= get_ppi_nodes()
    top_n_drugs= prediction_data %>% dplyr::slice_head(n=n) %>% 
        dplyr::select(drugbank_id,validation_status) %>%
        dplyr::rowwise()%>%
        dplyr::mutate(target_count = length(get_drug_targets(drugbank_id,ppi_gene_nodes,drug_target_mapping))) %>%
        dplyr::ungroup()
    
    output_file_name= here(output_dir,paste0("top_",n,"_drugs_",disease,"_score_filter_",score_filter,"_.csv"))
    write.csv(top_n_drugs,output_file_name,row.names = FALSE)
    message(sprintf("[SUCCESS] Top Drug file  saved at: %s",output_file_name))
            Sys.sleep(1.5)
}


create_prediction_disease_info = function(disease = c("MDD", "BD")) {
    disease = match.arg(disease)

    repodb_files  = c("MDD" = "MDD_REPODB.tsv", "BD" = "BD_REPODB.tsv")
    rank_files    = c("MDD" = "average_rank_MDD.csv", "BD" = "average_rank_BD.csv")
    file_suffixes = c("MDD" = "mdd", "BD" = "bipolar")

    repodb_file    = repodb_files[disease]
    rank_file_name = rank_files[disease]
    file_suffix    = file_suffixes[disease]

    drug_target_mapping = read_drug_targets()


    gold_standard = readr::read_tsv(
        here("src", "Data", "REPODB", repodb_file),
        show_col_types = FALSE
    )    

    gold_standard = gold_standard %>%  
        dplyr::semi_join(drug_target_mapping, by = 'drugbank_id') 

    cat(sprintf("[INFO] Total valid drugs in RepoDB (Gold Standard) for %s: %d\n", disease, nrow(gold_standard)))

    rank_file_path = here("src", "Data", "Drug", "Score", disease, rank_file_name)
    print(file.exists(rank_file_path))
    if (!file.exists(rank_file_path)) {
        cat("[WARN] Average Rank file not found. Building first...\n")
        drug_function_mapper[[1]]()
    }

    prediction_data = read.csv(rank_file_path, stringsAsFactors = FALSE)

    processed_predictions = prediction_data %>% 
        dplyr::mutate(
            validation_label = ifelse(drugbank_id %in% gold_standard$drugbank_id, 1, 0),
            validation_status = factor(
                validation_label,
                levels = c(0, 1),
                labels = c("Not validated by RepoDB", "Validated by RepoDB")
            )
        ) %>%
        dplyr::arrange(average_rank)

    validated_hits = processed_predictions %>% dplyr::filter(validation_label == 1)

    cat(sprintf("[INFO] Total predicted hits validated: %d\n", nrow(validated_hits)))

    if (nrow(validated_hits) <= 0) {
        cat("[WARN] ROC curve unavailable: No validated drugs were predicted.\n")
        return(NULL)
    }


    output_dir = here("src", "Evaluation", disease)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    score_filter = get_score_disease_gene_association()
    

    output_file_name = paste0("prediction_", file_suffix, "_score_filter_", score_filter, "_.csv")
    output_file_path = here(output_dir, output_file_name)

    write.csv(
        processed_predictions,
        file = output_file_path,
        row.names = FALSE
    )

    message(sprintf("[SUCCESS] prediction file saved at: %s", output_file_path))
    
    invisible(processed_predictions)
}
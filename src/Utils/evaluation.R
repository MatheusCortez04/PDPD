
evaluation_menu = function(){
  while(TRUE){
    clear_console()
    cat("--- Evaluation Menu ---\n\n")
    cat(" [1] Generate MDD ROC \n")
    cat(" [2] Generate MDD Recall@K \n")
    cat(" [3] Generate Bipolar ROC \n")
    cat(" [4] Generate Bipolar Recall@K \n")
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
    },
    '2' = function(){
         generate_recall_k_MDD()
    },
    '3'= function(){
       generate_roc_curve_bipolar()
    },
    '4'= function(){
        generate_recall_k_bipolar()
    },
    '5' = function(){
        generate_roc_to_kernel()
    }
)


generate_recall_k_MDD = function(){
    prediction_mdd_dir = here("src", "Evaluation","MDD")
    prediction_mdd_file_name ="prediction_mdd.csv"
    pred_mdd_rank_file_path = here(prediction_mdd_dir,prediction_mdd_file_name)
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

    recall_path_dir = here(prediction_mdd_dir,"RECALL")
    pdf_path <- here(recall_path_dir,"recall_at_k_MDD.pdf")
    dir.create(recall_path_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(pdf_path, plot = g, width = 8, height = 6)

    message(sprintf("[SUCCESS] Recall graph saved at: %s",here(pdf_path)))
    Sys.sleep(1.5)
}

generate_roc_curve_mdd = function(){
    output_dir = here("src", "Evaluation","MDD")
    prediction_mdd_file_path = here(output_dir,"prediction_mdd.csv")

    if(!file.exists(prediction_mdd_file_path)){
        cat("[WARN] Processed data not found. Running preparation first...\n")
        prediction_mdd = created_prediction_mdd_data()
    }
    else{
        prediction_mdd = read.csv(prediction_mdd_file_path)
    }
 
    if (sum(prediction_mdd$validation_label) <= 0) {
        cat("[WARN] Cannot plot ROC: No hits found in the ranking.\n")
        return(NULL)
    }
   
    output_roc_dir = here(output_dir,"ROC")
    dir.create(output_roc_dir, recursive = TRUE, showWarnings = FALSE)
    output_graph_file_name = "mdd_roc_curve.pdf"
    grDevices::pdf(here(output_roc_dir,output_graph_file_name), width = 6, height = 6)
    
    roc_results = reportROC::reportROC(
         gold = prediction_mdd$validation_label,
         predictor = -1 * prediction_mdd$Mean_Rank,
         plot = TRUE
     )
     grDevices::dev.off()
     message(sprintf("[SUCCESS] ROC curve saved at: %s",here(output_dir, output_graph_file_name)))
     Sys.sleep(1.5)

}

created_prediction_mdd_data = function(){
    drug_target_mapping  = load_drug_target_df()
    mdd_gold_standard_repodb = read_tsv(here("src","Data","REPODB","MDD_REPODB.tsv"), show_col_types = FALSE)    
    mdd_gold_standard_repodb = mdd_gold_standard_repodb %>%  dplyr::semi_join(drug_target_mapping,by='drugbank_id') %>% 
      filter(status=="Approved")
 
    cat(sprintf("[INFO] Total valid drugs in RepoDB (Gold Standard): %d\n",nrow(mdd_gold_standard_repodb)))

    mdd_gold_standard_open_targets = read.csv(here("src","Data","REPODB","openTargets_final_mapeado.csv")) %>%
      dplyr::semi_join(drug_target_mapping, by='drugbank_id')

    drugbank_id_repodb = unique(mdd_gold_standard_repodb$drugbank_id)
    drugbank_id_op   = unique(mdd_gold_standard_open_targets$drugbank_id) # Aqui está o segredo
    
    merge = unique(c(drugbank_id_op, drugbank_id_repodb))

    cat(sprintf("[INFO] Unique drugs in RepoDB: %d\n", length(drugbank_id_repodb)))    
    cat(sprintf("[INFO] Unique drugs in OpenTargets: %d\n", length(drugbank_id_op)))  
    cat(sprintf("[INFO] Total unique drugs in MERGE: %d\n", length(merge)))
    
    rank_file_path = here("src", "Data", "Drug", "Score", "MDD", "average_kernel_rank.csv")
    
    if (!file.exists(rank_file_path)) {
      stop(base::sprintf("[ERROR] Average rank file not found at: %s", rank_file_path))
    }
     mdd_prediction = read.csv(rank_file_path)
    processed_predictions = mdd_prediction %>% 
        dplyr::mutate(
            validation_label =  ifelse(drugbank_id %in% merge, 1, 0),
            validation_status = factor(validation_label,
                                             levels = c(0, 1),
                                             labels = c("Not validated", "Validated by Gold Standard"))
            )

    validated_hits = processed_predictions %>% filter(validation_label==1)
    cat(sprintf("[INFO] Total predicted hits validated: %d\n",nrow(validated_hits)))
    if (nrow(validated_hits) <= 0) {
         cat("[WARN] ROC curve unavailable: No validated drugs were predicted in the ranking.\n")
         Sys.sleep(1.5)
         return(NULL)
    }
    output_dir =  here("src", "Evaluation","MDD")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(processed_predictions,file=here(output_dir,"prediction_mdd.csv"),row.names=FALSE)
    message(sprintf("[SUCCESS] prediction file  saved at: %s",here(output_dir,"prediction_mdd.csv")))

    invisible(processed_predictions)
}
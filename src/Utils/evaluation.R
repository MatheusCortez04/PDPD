library(here)
library(dplyr)
library(reportROC)
source(here("src","Utils","utils.R"))

evaluation_menu = function(){
  while(TRUE){
    clear_console()
    cat("--- Evaluation Menu ---\n\n")
    cat(" [1] Generate MDD ROC \n")
    cat(" [2] Generate MDD Recall@K \n")
    cat(" [4] Generate Bipolar ROC \n")
    cat(" [5] Generate Bipolar Recall@K \n")
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
    }
)


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
         predictor = -1 * prediction_mdd$average_rank,
         plot = TRUE
     )
     grDevices::dev.off()
     message(sprintf("[SUCCESS] ROC curve saved at: %s",here(output_dir, output_graph_file_name)))
     Sys.sleep(1.5)

}

created_prediction_mdd_data = function(){
    drug_target_mapping  = load_drug_target_df()
    mdd_gold_standard = read_tsv(here("src","Data","REPODB","MDD_REPODB.tsv"), show_col_types = FALSE)    
    mdd_gold_standard = mdd_gold_standard %>%  dplyr::semi_join(drug_target_mapping,by='drugbank_id') %>% 
      filter(status=="Approved")
 
    cat(sprintf("[INFO] Total valid drugs in RepoDB (Gold Standard): %d\n",nrow(mdd_gold_standard)))
    rank_file_path = here("src", "Data", "Drug", "Score", "MDD", "average_rank_MDD.csv")
    
    if (!file.exists(rank_file_path)) {
      stop(base::sprintf("[ERROR] Average rank file not found at: %s", rank_file_path))
    }
     mdd_prediction = read.csv(rank_file_path)
    processed_predictions = mdd_prediction %>% 
        dplyr::mutate(
            validation_label =  ifelse(drugbank_id %in% mdd_gold_standard$drugbank_id, 1, 0),
            validation_status = factor(validation_label,
                                             levels = c(0, 1),
                                             labels = c("Not validated by RepoDB", "Validated by RepoDB"))
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

generate_recall_k_MDD = function(){
    output_dir = here("src", "Evaluation","MDD")
    prediction_mdd_file_path = here(output_dir,"prediction_mdd.csv")

    if(!file.exists(prediction_mdd_file_path)){
        cat("[WARN] Processed data not found. Running preparation first...\n")
        prediction_data = created_prediction_mdd_data()
    }
    else{
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
    print(g)

    output_recall_dir= here(output_dir,"Recall")
    dir.create(output_recall_dir, recursive = TRUE, showWarnings = FALSE)
    pdf_path = here(output_recall_dir,"recall_at_k_MDD.pdf")
    ggsave(pdf_path, plot = g, width = 8, height = 6)

}
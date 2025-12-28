library(here)
library(dplyr)
library(reportROC)
source(here("src","Utils","utils.R"))

evaluation_menu = function(){
  while(TRUE){
    clear_console()
    cat("--- Evaluation Menu ---\n\n")
    cat(" [1] Generate MDD ROC \n")
    cat(" [2] Generate MDD Recall@K Graph \n")
    cat(" [4] Generate Bipolar ROC \n")
    cat(" [5] Generate Bipolar Recall@K Graph \n")
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
    },
    '3'= function(){
        generate_recall_k_MDD()
    },
    '4'= function(){
        generate_roc_curve_bipolar()
    }
)


generate_roc_curve_mdd = function(){
  drug_target_mapping  = load_drug_target_df()

  mdd_gold_standard = read_tsv(here("src","Data","REPODB","MDD_REPODB.tsv"), show_col_types = FALSE)    
    mdd_gold_standard = mdd_gold_standard %>%  dplyr::semi_join(mdd_gold_standard,drug_target_mapping,by='drugbank_id') %>% 
      filter(status=="Approved")
 
  cat(sprintf("[INFO] Total valid drugs in RepoDB (Gold Standard): %d\n",nrow(mdd_gold_standard)))

  rank_file_path = here::here("src", "Data", "Drug", "Score", "MDD", "average_rank_MDD.csv")

  if (!file.exists(rank_file_path)) {
      stop(base::sprintf("[ERROR] Average rank file not found at: %s", rank_file_path))
    }

    mdd_prediction = read.csv(rank_file_path)
    
    processed_predictions = mdd_prediction %>% 
        dplyr::mutate(validation_label = ifelse(drugbank_id %in% mdd_gold_standard$drugbank_id, 1, 0))
    glimpse(processed_predictions)
    validated_hits = processed_predictions %>% filter(validation_label==1)
    cat(sprintf("[INFO] Total predicted hits validated: %d\n",nrow(validated_hits)))
    processed_predictions$validation_label = factor(processed_predictions$validation_label,
                           levels = c(0,1),
                           labels = c("Not validated by repODB", "Validated by repODB"))

    if (nrow(validated_hits) <= 0) {
         cat("[WARN] ROC curve unavailable: No validated drugs were predicted in the ranking.\n")
         Sys.sleep(1.5)
         return(NULL)
    }
     output_dir <- here("src", "Evaluation", "ROC", "MDD")
     output_graph_file_name = "mdd_roc.pdf" 

     dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

     write.csv(processed_predictions,file=here(output_dir,"prediction_mdd.csv"),row.names=FALSE)

     grDevices::pdf(here(output_dir,output_graph_file_name), width = 6, height = 6)
    
    roc_results = reportROC::reportROC(
         gold = processed_predictions$validation_label,
         predictor = -1 * processed_predictions$average_rank,
         plot = TRUE
     )
 
     grDevices::dev.off()
     message(sprintf("[SUCCESS] ROC curve saved at: %s",here(output_dir, output_graph_file_name)))
     Sys.sleep(1.5)
    
     return(roc_results)
}
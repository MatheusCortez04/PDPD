
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
    output_dir =here("src","Evaluation","MDD")
    roc_out_dir= here("src","Evaluation","MDD","ROC")
    output_graph_file_name = "mdd_roc.pdf" 
    dir.create(roc_out_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(pred_mdd,file=here(output_dir,"prediction_mdd.csv"),row.names=FALSE)

    grDevices::pdf(here(roc_out_dir,output_graph_file_name), width = 6, height = 6)
    roc_out = reportROC::reportROC(gold = pred_mdd$mdd_repodb_validated,
                                    predictor = -1*pred_mdd$Mean_Rank,
                                    plot=TRUE
                                    )
 
    grDevices::dev.off()
    message(sprintf("[SUCCESS] ROC curve saved at: %s",here(roc_out_dir, output_graph_file_name)))
    Sys.sleep(1.5)
}

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

generate_roc_to_kernel = function(){
    diseases = c("MDD","BD")
    drug_score_dir =here("src","Data","Drug","Score")
     
    
    kernel_names = c(
        "diffusion_kernel","pstep_kernel","regularised_laplacian_kernel",
        "commute_time_kernel","inverse_cosine_kernel"
    )
    for(disease in diseases){
        gold_standard = read_tsv(here("src","Data","REPODB",paste0(disease,"_REPODB.tsv")), show_col_types = FALSE)   
        gold_standard = gold_standard %>% filter(status=="Approved")
        for(kernel in kernel_names){
            input_file_path =here(drug_score_dir,disease,paste0(kernel,".csv")) 
            if (!file.exists(input_file_path)) {
                cat("[WARN] Drug Score to kernel:",kernel,"and disease:",disease,"\n")
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
                    predictor =kernel_predictions$drug_score,
                    plot = TRUE)

                
                grDevices::dev.off()
                message(sprintf("[SUCCESS] ROC curve saved at: %s",here(output_roc_dir, output_graph_file_name)))
                Sys.sleep(1.5)



        }
    }



}
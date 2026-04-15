library(here)
library(purrr)
library(dplyr)
source(here("src","Utils","utils.R"))
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
    }
)

get_drugbank_ids = function(){
    drug_target_df  = read.csv(here("src","Data","Drug","drug_targets_DrugBank_Gysi.csv")) %>%
        dplyr::distinct()

    unique_drugbank_ids = drug_target_df %>%
        dplyr::pull(drugbank_id) %>%
        unique()

    print(paste("Unique drugs in Drug-Target Dataframe:", length(unique_drugbank_ids)))
    invisible(unique_drugbank_ids)
}
get_drug_targets_in_ppi <- function(drugbank_id,ppi_gene_nodes) {
  
    drug_target_df = read.csv(here("src","Data","Drug","drug_targets_DrugBank_Gysi.csv")) %>% 
        distinct() 
  

    drug_target_proteins = drug_target_df %>%
        filter(drugbank_id == !!drugbank_id) %>%
        distinct() %>%
        dplyr::mutate(entrez_id = as.character(entrez_id)) %>%
        pull(entrez_id)

    drug_targets_in_ppi <- intersect(drug_target_proteins, ppi_gene_nodes)

    invisible(drug_targets_in_ppi)
}

calculate_drug_score = function(){

        kernel_names = c(
        "diffusion_kernel","pstep_kernel","regularised_laplacian_kernel",
        "commute_time_kernel","inverse_cosine_kernel"
    )

    disease_data <- list(
        MDD = list(name = "MDD", genes = get_mdd_disease_module()$entrez_id),
        BD = list(name = "BD", genes = get_bipolar_disease_module()$entrez_id)
    )

    # 🔹 garantir que a pasta existe
    dir.create(here("src","Data","Drug","Score"), recursive = TRUE, showWarnings = FALSE)

    drugbank_ids = get_drugbank_ids()
    ppi_gene_nodes <- get_ppi_nodes()

    drug_targets = purrr::map(drugbank_ids, get_drug_targets_in_ppi, ppi_gene_nodes)
    names(drug_targets) = drugbank_ids

    for(kernel_name in kernel_names){

        kernel_path <- here("src","Data","Kernels","RData", paste0(kernel_name, ".Rdata"))

        if (!file.exists(kernel_path)) {
            cat("Kernel", kernel_name, "não encontrado\n")
            next
        }

        kernel_matrix <- load_rdata(kernel_path)

        for(disease in names(disease_data)){

            disease_name <- disease_data[[disease]]$name
            disease_genes <- disease_data[[disease]]$genes

            scores <- c()

            for(drug_id in names(drug_targets)){

                targets <- drug_targets[[drug_id]]

                if (length(targets) == 0) {
                    score <- NA
                } else {
                    valid_targets <- intersect(targets, rownames(kernel_matrix))
                    valid_disease <- intersect(disease_genes, colnames(kernel_matrix))
                    score <- mean(kernel_matrix[valid_targets, valid_disease, drop = FALSE])
                }

                scores[drug_id] <- score
            }

            # 🔹 criar dataframe FORA do loop
            result_df <- data.frame(
                drugbank_id = names(scores),
                score = as.numeric(scores),
                stringsAsFactors = FALSE
            )

            file_name <- paste0("drug_scores_", disease_name, "_", kernel_name, ".csv")
            file_path <- here("src", "Data", "Drug", "Score", file_name)

            write.csv(result_df, file_path, row.names = FALSE)

            cat("Salvo:", file_name, "\n")
        }
    }
}


generate_drug_rank = function(disease = c("MDD", "BD")) {

    disease = match.arg(disease)

    cat("[INFO] Building Average Rank to Disease:", disease, "\n")

    kernel_names = c(
        "diffusion_kernel","pstep_kernel","regularised_laplacian_kernel",
        "commute_time_kernel","inverse_cosine_kernel"
    )

    score_base_dir = here("src", "Data", "Drug", "Score")

    score_list = list()

    for (kernel in kernel_names) {

        file_name <- paste0("drug_scores_", disease, "_", kernel, ".csv")
        file_path <- here(score_base_dir, file_name)

        if (!file.exists(file_path)) {
            cat("[WARN] CSV not found:", file_path, "\n")
            next
        }

        df <- read.csv(file_path, stringsAsFactors = FALSE)

        # 🔬 calcular rank (quanto maior score, melhor)
        df[[kernel]] <- rank(-df$score, ties.method = "average")

        df <- df[, c("drugbank_id", kernel)]

        score_list[[kernel]] <- df
    }

    # 🔴 segurança
    if (length(score_list) == 0) {
        stop("Nenhum arquivo de score foi encontrado.")
    }

    # 🔬 merge de todos os kernels
    final_rank_df = Reduce(function(x, y) merge(x, y, by = "drugbank_id", all = TRUE), score_list)

    # 🔬 média dos ranks
    final_rank_df$average_rank = rowMeans(final_rank_df[, kernel_names], na.rm = TRUE)

    # 🔬 ordenar (menor rank = melhor)
    final_rank_df = final_rank_df[order(final_rank_df$average_rank), ]

    # 🔬 salvar
    output_csv = here(score_base_dir, paste0("average_rank_", disease, ".csv"))
    write.csv(final_rank_df, output_csv, row.names = FALSE)

    cat("\n✔ Average Rank salvo em:", output_csv, "\n")

    return(final_rank_df)
}


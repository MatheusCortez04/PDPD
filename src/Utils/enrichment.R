library(here)
library(dplyr)
library(tibble)
library(tidyr)
library(enrichR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
source(here("src","Utils","utils.R"))
source(here("src","Utils","graph_metrics.R"))

enrichment_menu = function(){
  while(TRUE){
    clear_console()
    cat("--- Evaluation Menu ---\n\n")
    cat(" [1] Drug Set Enrichment Analysis MDD \n")
    cat(" [2] Drug Set Enrichment Analysis BD \n")
    cat(" [3] Over Representation Analysis MDD \n")
    cat(" [4] Over Representation Analysis BD \n")
    cat(" [B] Back\n\n")

    input = readline(prompt = "Choice option: ")
    if (toupper(trimws(input)) == "B") break
     enrichment_function_mapper[[input]]()
  }
}

enrichment_function_mapper = list(
    '1' = function() {
        dsea_rank_by_disease("MDD")
    },
    '2' = function(){
        dsea_rank_by_disease("BD")
    },
    '3'= function(){
        over_representation_analysis_by_disease("MDD")
    },
    '4'= function(){
        over_representation_analysis_by_disease("BD")
    }
)

dsea_rank_by_disease = function(disease = c("MDD", "BD")) {
    disease = match.arg(disease)
    cat(sprintf("[INFO] Running GSEA enrichment pipeline for %s\n", disease))
    
    score_filter = get_score_disease_gene_association()
    ranking_dir = here("src", "Evaluation", disease)
    drug_ranking_filename = c("MDD" = sprintf("prediction_mdd_score_filter_%s_.csv", score_filter), 
                       "BD"  = sprintf("prediction_bipolar_score_filter_%s_.csv", score_filter))
    
    if(file.exists(here(ranking_dir, drug_ranking_filename[disease]))) {
        drug_ranking_df = read.csv(here(ranking_dir, drug_ranking_filename[disease])) %>%
            dplyr::select(drugbank_id, average_rank)
    }
    if(!file.exists(here(ranking_dir, drug_ranking_filename[disease]))) {
        cat("[WARN] Processed data not found. Running preparation first...\n")
        drug_ranking_df = create_prediction_disease_info(disease) %>%
            dplyr::select(drugbank_id, average_rank)
    }


    gsea_input =  prepare_drug_gsea_input(drug_ranking_df)
    drug_target_annotation_df =gsea_input$drug_target_annotation_df
    ranked_drug_vector =gsea_input$ranked_drug_vector

    enrichr_disease_response = run_enrichr_annotation(drug_target_annotation_df$symbol)
    enrichr_dbs = c("KEGG_2026","GO_Biological_Process_2026","Reactome_Pathways_2024")

    purrr::walk(enrichr_dbs, function(db_name) {
        cat(sprintf("\n[INFO] Processing database: %s\n",db_name))

        pathway_gene_annotation_df = expand_enrichr_results(enrichr_disease_response,db_name)

        drug_pathway_annotation =build_drug_pathway_annotation(drug_target_annotation_df,pathway_gene_annotation_df)

        pathway_summary_df = build_dug_pathway_summary(drug_pathway_annotation)

        output_dir <- here("src","Enrichment","Drugs",disease,db_name)
        create_dir(output_dir)
        write.csv(
            pathway_summary_df,
            file = here(output_dir,sprintf("%s_%s_pathway_summary.csv",disease,db_name)),
            row.names = FALSE)
        cat(sprintf("[INFO] Saved pathway summary: %s_%s_pathway_summary.csv\n",disease,db_name))

        cat("[INFO] Building TERM2Drug mapping...\n")
        term2Drug_df <- drug_pathway_annotation %>%
            dplyr::filter(Term %in% drug_pathway_annotation$Term) %>%
            dplyr::filter(!is.na(Term), 
                Term != "NA", 
                trimws(Term) != "", 
                !is.na(drugbank_id),
                trimws(drugbank_id) != ""
                ) %>%
            dplyr::select(term = Term, gene = drugbank_id) %>%
            dplyr::distinct() 
        
 
        cat("[INFO] Running GSEA...\n")


         gsea_result = clusterProfiler::GSEA(
             geneList      = ranked_drug_vector,
             TERM2GENE     = term2Drug_df,
             pvalueCutoff  = 0.05,
             eps           = 0,
             pAdjustMethod = "BH",
             nPermSimple   = 10000
         )
        
    
        if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
            cat(sprintf("[WARNING] No significant pathways identified for %s.\n", db_name))
            return() 
        }
        gsea_result_df = gsea_result@result %>% 
            dplyr::select(-log2err) %>% 
            dplyr::distinct()%>%
            dplyr::filter(
                p.adjust < 0.05,
                enrichmentScore > 0) %>% 
            dplyr::arrange(desc(NES))

        output_filename = sprintf("%s_gsea_%s_.csv", disease,db_name)
        write.csv(gsea_result_df, file=here(output_dir,output_filename), row.names = FALSE)
        cat(sprintf("[INFO] Saving results to: %s\n",output_filename))
        graph_output_filename = sprintf("%s_gsea_%s_.pdf", disease,db_name)
        grDevices::pdf(here(output_dir,graph_output_filename), width = 10, height = 8)
        print(
            enrichplot::gseaplot2(
            gsea_result,
            geneSetID = 1:5, 
            pvalue_table = FALSE, 
            subplots = 1:2
            )
        )
        grDevices::dev.off()
    })

    cat("[INFO] GSEA enrichment pipeline completed successfully.\n")

}


prepare_drug_gsea_input = function(drug_ranking_df){
    set.seed(42)

    cat("[INFO] Extracting drug targets...\n")
    drug_target_df = import_drug_targets_df()

    drug_target_list = purrr::map(
        drug_ranking_df$drugbank_id,
        extract_targets_per_drug,
        drug_target_df)

    names(drug_target_list) = drug_ranking_df$drugbank_id

    drug_target_annotation_df = tibble::tibble(
        drugbank_id = names(drug_target_list),
        entrez_id   = drug_target_list
    ) %>%tidyr::unnest(entrez_id)

    gene_symbols = AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys    = unique(as.character(drug_target_annotation_df$entrez_id)),
        column  = "SYMBOL",
        keytype = "ENTREZID"
    )

    drug_target_annotation_df = drug_target_annotation_df %>%
        dplyr::mutate(symbol = gene_symbols[as.character(entrez_id)]) %>%
        dplyr::filter(!is.na(symbol)) %>%
        dplyr::distinct(
            drugbank_id,
            entrez_id,
            symbol) %>%
        dplyr::left_join(
            drug_ranking_df,
            by = "drugbank_id"
        )
        
        cat("[INFO] Building ranked drug vector...\n")

        tie_breaker = runif(
            nrow(drug_ranking_df),
            min = -1e-6,
            max = 1e-6
        )

        drug_score_df = drug_ranking_df %>%
        dplyr::mutate(
            gsea_score = (-log10(average_rank + 1)) + tie_breaker
        ) %>%
        dplyr::arrange(desc(gsea_score))

    ranked_drug_vector = drug_score_df$gsea_score
    names(ranked_drug_vector) = drug_score_df$drugbank_id

    return(
        list(drug_target_annotation_df = drug_target_annotation_df,
            ranked_drug_vector = ranked_drug_vector))

}

run_enrichr_annotation = function(gene_symbols){

    gene_symbols = gene_symbols %>%
        unique() %>%
        stats::na.omit() %>%
        as.character()

    enrichr_databases = c(
        "KEGG_2026",
        "GO_Biological_Process_2026",
        "Reactome_Pathways_2024"
    )

    cat(sprintf(
        "[INFO] Running Enrichr annotation for %d unique genes...\n",
        length(gene_symbols)
    ))

    enrichr_response = enrichR::enrichr(
        gene_symbols,
        enrichr_databases
    )

    return(enrichr_response)
}
expand_enrichr_results = function(enrichr_response, db_name){
    enrichr_response[[db_name]] %>%
        tidyr::separate_rows(Genes, sep = ";") %>%
        dplyr::mutate(
            symbol = Genes
        ) %>%
        dplyr::select(-Genes) %>%
        dplyr::distinct()

}
build_drug_pathway_annotation = function(drug_target_annotation_df,pathway_gene_annotation_df){
    drug_target_annotation_df %>%
        dplyr::inner_join(
            pathway_gene_annotation_df,
            by = "symbol",
            relationship = "many-to-many"
        ) %>%
        dplyr::select(
            Term,
            P.value,
            Adjusted.P.value,
            Combined.Score,
            symbol,
            entrez_id,
            drugbank_id
        ) %>%
        dplyr::distinct()

}

build_dug_pathway_summary = function(drug_pathway_annotation){

    drug_pathway_annotation %>%
        dplyr::group_by(
            Term,
            Adjusted.P.value,
            P.value,
            Combined.Score
        ) %>%
        dplyr::summarise(
            target_count = dplyr::n_distinct(symbol),
            targets = paste(unique(symbol), collapse = ", "),
            drug_count = dplyr::n_distinct(drugbank_id),
            drugs = paste(unique(drugbank_id), collapse = ", "),
            .groups = "drop"
        ) %>%
        dplyr::arrange(desc(drug_count)) %>%
        dplyr::select(
            Term,
            Adjusted.P.value,
            P.value,
            Combined.Score,
            drug_count,
            target_count,
            drugs,
            targets
        ) %>%
        dplyr::distinct() %>%
        dplyr::filter(drug_count < 1000)

}
build_disease_pathway_summary = function(disease_gene_pathways){

    disease_gene_pathways %>%
        dplyr::group_by(
            Term,
            Adjusted.P.value,
            P.value,
            Combined.Score
        ) %>%
        dplyr::arrange(desc(drug_count)) %>%
        dplyr::select(
            Term,
            Adjusted.P.value,
            P.value,
            Combined.Score,
        ) %>%
        dplyr::distinct()

}

over_representation_analysis_by_disease = function(disease = c("MDD", "BD")){
    output_dir <- here("src","Enrichment",disease)

    disease <- match.arg(disease)
    cat(sprintf("[INFO] Running Over-Representation Analysis for %s\n", disease))
    if (!dir.exists(output_dir)) {
        create_dir(output_dir)
    }

    disease_gene_df = get_disease_genes(disease)

    gene_symbols = AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys    = unique(as.character(disease_gene_df$entrez_id)),
        column  = "SYMBOL",
        keytype = "ENTREZID")

    
    disease_gene_df = disease_gene_df %>% 
        dplyr::mutate(symbol = gene_symbols[as.character(entrez_id)]) %>%
        dplyr::filter(!is.na(symbol))%>%
        dplyr::distinct()

    
    enrichment_response = run_enrichr_annotation(disease_gene_df$symbol)
    disease_kegg_expand =expand_enrichr_results(enrichment_response,"KEGG_2026")
    
    disease_gene_pathways = disease_gene_df %>% 
        dplyr::inner_join(disease_kegg_expand,by='symbol') %>%
        dplyr::select(Term,
            disease_id,
            entrez_id,
            symbol,
            P.value,
            Adjusted.P.value,
            Combined.Score)
    graph_title_mapper   = c("MDD" = "Over Representation Analysis - Transtorno Depressivo", "BD" = "Over Representation Analysis - Transtorno Bipolar")
    graph_title   = graph_title_mapper[disease]
    enrichr_summary <- disease_gene_pathways %>%
        dplyr::group_by(Term) %>%
        dplyr::summarise(
            Count = n_distinct(symbol),
            Genes = paste(unique(symbol), collapse = "/"),
            Adjusted.P.value = min(Adjusted.P.value),
            Combined.Score = max(Combined.Score),
          .groups = "drop"
        ) %>% 
        dplyr::arrange(Adjusted.P.value)

    total_genes = n_distinct(disease_gene_df$symbol)
    plot_df <- enrichr_summary %>%
        dplyr::slice_head(n = 10) %>%
        dplyr::mutate(
            GeneRatio = Count / total_genes,
            logP = -log10(Adjusted.P.value)
        )
    p <- ggplot(plot_df, aes(
        x = GeneRatio,
        y = reorder(Term, GeneRatio)
    )) +
        geom_point(aes(
            size = Count,
            color = logP
        )) +
        scale_color_gradient(low = "red", high = "blue") +
        theme_minimal() +
        labs(
            title = graph_title,
            x = "Gene Ratio",
            y = "KEGG Pathway",
            color = "-log10(adj p-value)",
            size = "Gene Count"
        )

    plot_file <- file.path(output_dir,sprintf("%s_ORA_KEGG_dotplot.pdf", disease))
    
    grDevices::pdf(plot_file, width = 10, height = 8)
    print(p)
    grDevices::dev.off()
    
    cat(sprintf("[INFO] Plot saved at: %s\n", plot_file))
    table_file <- file.path(
        output_dir,
        sprintf("%s_ORA_KEGG_summary.csv", disease)
    )
    write.csv(enrichr_summary, table_file, row.names = FALSE)

    cat(sprintf("[INFO] Table saved at: %s\n", table_file))
}



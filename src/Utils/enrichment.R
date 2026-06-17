library(here)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
source(here("src","Utils","utils.R"))
source(here("src","Utils","graph_metrics.R"))
library(igraph)
library(tibble)
library(ggplot2)
library(tidyr)
library(enrichR)

calculate_go_enrichment = function(entrez_ids, disease, output_file_prefix) {
    output_dir = here("src", "Enrichment", disease)
    output_rdata_dir = here(output_dir, "RData")
    create_dir(output_dir)
    create_dir(output_rdata_dir)
    
    output_rdata_file = here(output_rdata_dir, sprintf("%s_GO_analysis.RData", output_file_prefix))
    

    cat(sprintf("[INFO] Running GO (Biological Process) enrichment analysis with %d genes...\n", length(entrez_ids)))

    enrich_go = clusterProfiler::enrichGO(
        gene = entrez_ids,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        ont = "BP",
        OrgDb = org.Hs.eg.db)

    if (is.null(enrich_go) || nrow(as.data.frame(enrich_go)) == 0) {
        cat("[WARNING] No significant GO terms identified.\n")
        return(invisible(NULL))
    }

    cat(sprintf("[INFO] Initial significant GO terms identified: %d\n", nrow(as.data.frame(enrich_go))))
    

    enrich_go_simplify = clusterProfiler::simplify(enrich_go, cutoff = 0.7, by = "p.adjust", select_fun = min)
    go_enrichment_simplify_df = as.data.frame(enrich_go_simplify)
    output_csv_file = here(output_dir, sprintf("%s_simplify_GO_analysis_%s_to_%s.csv", output_file_prefix,nrow(enrich_go),nrow(enrich_go_simplify)))
    cat(sprintf("[INFO] After simplification, unique GO terms remaining: %d\n", nrow(go_enrichment_simplify_df)))
    cat(sprintf("[INFO] Saving GO enrichment results to: %s\n", output_csv_file))
    
    save(enrich_go_simplify, file = output_rdata_file)
    write.csv(go_enrichment_simplify_df, file = output_csv_file, row.names = FALSE)
    
    return(enrich_go_simplify)
}

calculate_kegg_enrichment = function(entrez_ids, disease, output_file_prefix) {
    output_dir = here("src", "Enrichment", disease)
    output_rdata_dir = here(output_dir, "RData")
    create_dir(output_dir)
    create_dir(output_rdata_dir)
    
    output_rdata_file = here(output_rdata_dir, sprintf("%s_KEGG_analysis.RData", output_file_prefix))
    output_csv_file = here(output_dir, sprintf("%s_KEGG_analysis.csv", output_file_prefix))

    cat(sprintf("[INFO] Running KEGG pathway enrichment with %d genes...\n", length(entrez_ids)))
    
    enrich_kegg = clusterProfiler::enrichKEGG(
        gene = entrez_ids,
        organism = "hsa",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH"
    )
    
    if (is.null(enrich_kegg) || nrow(as.data.frame(enrich_kegg)) == 0) {
        cat("[WARNING] No significant KEGG pathways identified.\n")
        return(invisible(NULL))
    }
    
    enrich_kegg_df = as.data.frame(enrich_kegg)
    
    cat(sprintf("[INFO] Significant KEGG pathways identified: %d\n", nrow(enrich_kegg_df)))
    cat(sprintf("[INFO] KEGG enrichment results saved to: %s\n", output_csv_file))

    save(enrich_kegg, file = output_rdata_file)
    write.csv(enrich_kegg_df, file = output_csv_file, row.names = FALSE)
    

    return(enrich_kegg)

}

calculate_reactome_enrichment = function(entrez_ids, disease, output_file_prefix){
    output_dir = here("src", "Enrichment", disease)
    output_rdata_dir = here(output_dir, "RData")
    create_dir(output_dir)
    create_dir(output_rdata_dir)
    
    output_rdata_file = here(output_rdata_dir, sprintf("%s_Reactome_Analysis.RData", output_file_prefix))
    output_csv_file = here(output_dir, sprintf("%s_Reactome_Analysis.csv", output_file_prefix))

       cat(sprintf("[INFO] Running Reactome pathway enrichment with %d genes...\n", length(entrez_ids)))
    
    enrich_reactome = ReactomePA::enrichPathway(
        gene          = entrez_ids,
        pvalueCutoff  = 0.05,
        pAdjustMethod = "BH",
        maxGSSize     = 500
    )
    
    if (is.null(enrich_reactome) || nrow(as.data.frame(enrich_reactome)) == 0) {
        cat("[WARNING] No significant Reactome pathways identified.\n")
        return(invisible(NULL))
    }

    save(enrich_reactome, file = output_rdata_file)

    write.csv(as.data.frame(enrich_reactome), file = output_csv_file, row.names = FALSE)
    return(enrich_reactome)


}

calculate_enrichment_by_disease = function(disease = c("MDD", "BD")){

    disease = match.arg(disease)
    cat(sprintf("[INFO] Running enrichment pipeline for disease: %s\n",disease))

    disease_file_prefix = ifelse(disease == "MDD","mdd_all_genes","bipolar_all_genes")

    disease_genes = get_disease_genes(disease) %>%dplyr::pull(entrez_id)

    enrich_go= calculate_go_enrichment(entrez_ids= disease_genes,disease,output_file_prefix = disease_file_prefix)
    enrich_kegg= calculate_kegg_enrichment(entrez_ids= disease_genes,disease,output_file_prefix = disease_file_prefix)
    enrich_reactome= calculate_reactome_enrichment(entrez_ids= disease_genes,disease,output_file_prefix = disease_file_prefix)
    cat(sprintf("[INFO] Enrichment pipeline completed for disease: %s\n",disease))

    return(list(
        gene_ontology=enrich_go,
        kegg=enrich_kegg,
        reactome=enrich_reactome
    ))
}

calculate_enrichment_by_disease_specific_genes = function(disease = c("MDD", "BD")){

    disease = match.arg(disease)
    cat(sprintf("[INFO] Running enrichment pipeline for disease: %s\n",disease))

    specific_gene_file_prefix = ifelse(disease == "MDD","mdd_specific_genes","bipolar_specific_genes")

    disease_specific_genes = extract_specific_genes(disease) %>% pull(entrez_id)

    enrich_go = calculate_go_enrichment(entrez_ids= disease_specific_genes,disease,output_file_prefix = specific_gene_file_prefix)
    enrich_kegg =  calculate_kegg_enrichment(entrez_ids= disease_specific_genes,disease,output_file_prefix = specific_gene_file_prefix)
    enrich_reactome=calculate_reactome_enrichment(entrez_ids= disease_specific_genes,disease,output_file_prefix = specific_gene_file_prefix)
    cat(sprintf("[INFO] Enrichment pipeline completed for disease: %s\n",disease))

    return(list(
        gene_ontology=enrich_go,
        kegg=enrich_kegg,
        reactome=enrich_reactome
    ))
}

calculate_enrichment_common_genes = function(){
    cat(sprintf("[INFO] Running enrichment pipeline for Common Genes\n"))
    common_genes_file_prefix='common_genes'


    disease_common_genes = get_common_gene() 

    enrich_go = calculate_go_enrichment(entrez_ids= disease_common_genes,'Common_Genes',output_file_prefix = common_genes_file_prefix)
    enrich_kegg =  calculate_kegg_enrichment(entrez_ids= disease_common_genes,'Common_Genes',output_file_prefix = common_genes_file_prefix)
    enrich_reactome=calculate_reactome_enrichment(entrez_ids= disease_common_genes,'Common_Genes',output_file_prefix = common_genes_file_prefix)
    cat(sprintf("[INFO] Enrichment pipeline completed for Common Genes\n"))

    return(list(
        gene_ontology=enrich_go,
        kegg=enrich_kegg,
        reactome=enrich_reactome
    ))
}
comparate_cluster_go = function(entrez_ids,output_file_name){
    output_dir = here("src", "Enrichment",'Comparison','GeneOntology')
    output_rdata_dir = here(output_dir, "RData")
    output_csv_file= here(output_dir,sprintf("%s_comparison_go_simplify.csv",output_file_name))
    output_rdata_file = here(output_rdata_dir, sprintf("%s_comparison_go_simplify.RData", output_file_name))

    if (file.exists(output_rdata_file)) {
        cat("[INFO] Matrix already exists for this disease. Loading file...\n")
       load_cluster = load_rdata(output_rdata_file)
        return(invisible(load_cluster)) 
        
    }
    create_dir(output_dir)
    create_dir(output_rdata_dir)

    cat(sprintf("[INFO] Running Compare Cluster enrichment GO analysis with %d clusters...\n", length(entrez_ids)))

    compare_cluster = clusterProfiler::compareCluster(
        geneClusters = entrez_ids,
        fun= "enrichGO",
        OrgDb= org.Hs.eg.db,
        ont= "BP",
        pAdjustMethod= "BH",
        pvalueCutoff= 0.05,
        qvalueCutoff= 0.05)


        
    if (is.null(compare_cluster) || nrow(as.data.frame(compare_cluster)) == 0) {
        cat("[WARNING] No significant GeneOntology pathways identified.\n")
        return(invisible(NULL))
    }
    cat("[INFO] Simplifying redundant GO terms...\n")
    compare_cluster_simplify = simplify(compare_cluster, cutoff = 0.7, by = "p.adjust", select_fun = min)

    cat("[INFO] Saving RData and CSV files to disk...\n")
    save(compare_cluster_simplify, file = output_rdata_file)
    write.csv(as.data.frame(compare_cluster_simplify), file = output_csv_file, row.names = FALSE)
    return(invisible(compare_cluster_simplify))

}
comparate_cluster_kegg = function(entrez_ids,output_file_name){
    output_dir = here("src", "Enrichment",'Comparison','KEGG')
    output_rdata_dir = here(output_dir, "RData")
    output_csv_file= here(output_dir,sprintf("%s_comparison_kegg.csv",output_file_name))
    output_rdata_file = here(output_rdata_dir, sprintf("%s_comparison_kegg.RData", output_file_name))

    if (file.exists(output_rdata_file)) {
        cat("[INFO] Matrix already exists for this disease. Loading file...\n")
       load_cluster = load_rdata(output_rdata_file)
        return(invisible(load_cluster)) 
        
    }
    create_dir(output_dir)
    create_dir(output_rdata_dir)

    cat(sprintf("[INFO] Running Compare Cluster enrichment KEGG analysis with %d clusters...\n", length(entrez_ids)))

    compare_cluster_kegg = clusterProfiler::compareCluster(geneClusters = entrez_ids,
        fun= "enrichKEGG", 
        organism = "hsa", 
        keyType = "kegg", 
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH",
        maxGSSize = 500, 
        qvalueCutoff = 0.2)


    if (is.null(compare_cluster_kegg) || nrow(as.data.frame(compare_cluster_kegg)) == 0) {
        cat("[WARNING] No significant KEGG pathways identified.\n")
        return(invisible(NULL))
    }

    cat("[INFO] Saving RData and CSV files to disk...\n")
    save(compare_cluster_kegg, file = output_rdata_file)
    write.csv(as.data.frame(compare_cluster_kegg), file = output_csv_file, row.names = FALSE)
    return(invisible(compare_cluster_kegg))

}
comparate_cluster_reactome = function(entrez_ids,output_file_name){
    output_dir = here("src", "Enrichment",'Comparison','Reactome')
    output_rdata_dir = here(output_dir, "RData")
    output_csv_file= here(output_dir,sprintf("%s_comparison_reactome.csv",output_file_name))
    output_rdata_file = here(output_rdata_dir, sprintf("%s_comparison_reactome.RData", output_file_name))

    if (file.exists(output_rdata_file)) {
        cat("[INFO] Matrix already exists for this disease. Loading file...\n")
       load_cluster = load_rdata(output_rdata_file)
        return(invisible(load_cluster)) 
        
    }
    create_dir(output_dir)
    create_dir(output_rdata_dir)

    cat(sprintf("[INFO] Running Compare Cluster enrichment Reactome analysis with %d clusters...\n", length(entrez_ids)))

    compare_cluster_reactome = clusterProfiler::compareCluster(
        geneClusters = entrez_ids,
        fun= "enrichPathway", 
        organism = "human", 
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH",
        maxGSSize = 500, 
        qvalueCutoff = 0.2)


    if (is.null(compare_cluster_reactome) || nrow(as.data.frame(compare_cluster_reactome)) == 0) {
        cat("[WARNING] No significant Reactome pathways identified.\n")
        return(invisible(NULL))
    }

    cat("[INFO] Saving RData and CSV files to disk...\n")
    save(compare_cluster_reactome, file = output_rdata_file)
    write.csv(as.data.frame(compare_cluster_reactome), file = output_csv_file, row.names = FALSE)
    return(invisible(compare_cluster_reactome))

}


compare_disease_modules= function(){
    cat("[INFO] Preparing gene modules for cross-comparison...\n")

    mdd_specific = as.character(extract_specific_genes("MDD") %>% dplyr::pull(entrez_id))
    bd_specific = as.character(extract_specific_genes("BD") %>% dplyr::pull(entrez_id))
    common_genes = as.character(get_common_gene())

    module_list = list(
        "MDD (Specific)" = mdd_specific,
        "BD (Specific)"  = bd_specific,
        "Commons"        = common_genes
    )

    cat("[INFO] Running compareCluster (Gene Ontology). This may take a few minutes...\n")
    disease_comparison_go = comparate_cluster_go(module_list,"disease")
    disease_comparison_kegg = comparate_cluster_kegg(module_list,"disease")
    disease_comparison_reactome = comparate_cluster_reactome(module_list,"disease")
    return(invisible(
        list(
            go=disease_comparison_go,
            kegg=disease_comparison_kegg,
            reactome =disease_comparison_reactome
        )
    ))
}
sea_rank_drugs_by_disease = function(disease = c("MDD", "BD")) {
    disease = match.arg(disease)
    cat(sprintf("[INFO] Running GSEA enrichment pipeline for %s\n", disease))

    score_filter = get_score_disease_gene_association()
    rank_dir = here("src", "Evaluation", disease)
    rank_file_name = c("MDD" = sprintf("prediction_mdd_score_filter_%s_.csv", score_filter), 
                       "BD"  = sprintf("prediction_bipolar_score_filter_%s_.csv", score_filter))
    

    rank_drug_disease = read.csv(here(rank_dir, rank_file_name[disease])) %>%
        dplyr::select(drugbank_id, average_rank)

    cat("[INFO] Loading drug-target interactions...\n")
    drug_target_df = import_drug_targets_df()
    
    cat("[INFO] 3. Extraindo alvos específicos para os fármacos candidatos...\n")
    rank_drug_targets = purrr::map(rank_drug_disease$drugbank_id, extract_targets_per_drug, drug_target_df) 
    names(rank_drug_targets) = rank_drug_disease$drugbank_id

    rank_drug_targets_symbol = map_ranked_drug_targets_to_symbols(rank_drug_targets)

    cat("[INFO] Preparing ranked gene list for GSEA...\n")
    matriz_ranqueada <- rank_drug_disease %>%
        dplyr::mutate(
            escore_log = -log10(average_rank + 1),
            escore_gsea = as.numeric(scale(escore_log))
        ) %>%
        dplyr::arrange(desc(escore_gsea))

    fita_matematica <- matriz_ranqueada$escore_gsea
    names(fita_matematica) <- matriz_ranqueada$drugbank_id

    set.seed(42)
    ruido_estatistico <- runif(length(fita_matematica), min = -1e-6, max = 1e-6)
    fita_matematica_estabilizada <- sort(fita_matematica + ruido_estatistico, decreasing = TRUE)

    enrichr_dbs = c("KEGG_2026", "GO_Biological_Process_2026", "Reactome_Pathways_2024")
    all_unique_targets = unique(rank_drug_targets_symbol$symbol)
    
    cat(sprintf("[INFO] Running Enrichr annotation for %d unique target genes...\n",length(all_unique_targets)))
    enrichr_disease_response <- enrichR::enrichr(all_unique_targets, enrichr_dbs)


    purrr::walk(enrichr_dbs, function(db_name) {
        cat(sprintf("\n[INFO] Processing database: %s\n",db_name))
        enrichr_db_df <- enrichr_disease_response[[db_name]] %>%
            tidyr::separate_rows(Genes, sep = ";") %>%
            dplyr::mutate(symbol = Genes) %>%
            dplyr::select(-Genes) %>%
            dplyr::distinct() 

        cat("[INFO] Mapping drugs to enriched pathways...\n")
        drug_pathways <- rank_drug_targets_symbol %>%
            dplyr::inner_join(enrichr_db_df, by = "symbol", relationship = "many-to-many") %>%
            dplyr::select(Term, P.value, Adjusted.P.value, Combined.Score, symbol, drugbank_id) %>% 
            dplyr::distinct()

        matriz_agrupada <- drug_pathways %>%
            dplyr::group_by(Term, Adjusted.P.value, P.value, Combined.Score) %>%
            dplyr::summarise(
                target_count = n_distinct(symbol),
                targets = paste(unique(symbol), collapse = ", "),
                drug_count = n_distinct(drugbank_id),
                drugs = paste(unique(drugbank_id), collapse = ", "),
                .groups = "drop"
            ) %>%
            dplyr::arrange(desc(drug_count)) %>% 
            dplyr::select(Term, Adjusted.P.value, P.value, Combined.Score, drug_count, target_count, drugs, targets) %>%
            dplyr::distinct()
            
        matriz_agrupada_filtro = matriz_agrupada %>% dplyr::filter(drug_count < 1000)
        
        cat(sprintf("[INFO] Estrutura Inicial -> Termos brutos: %d | Termos filtrados (<1000): %d\n", 
                    nrow(matriz_agrupada), nrow(matriz_agrupada_filtro)))

        cat("[INFO] Building TERM2GENE mapping...\n")
        term2Gene_df <- drug_pathways %>%
            dplyr::filter(Term %in% matriz_agrupada_filtro$Term) %>%
            dplyr::filter(!is.na(Term), Term != "NA", trimws(Term) != "", !is.na(drugbank_id)) %>%
            dplyr::select(term = Term, gene = drugbank_id) %>%
            dplyr::distinct() 

        cat("[INFO] Running GSEA...\n")
        result_gsea_db = clusterProfiler::GSEA(
            geneList      = fita_matematica_estabilizada,
            TERM2GENE     = term2Gene_df,
            pvalueCutoff  = 0.05,
            eps           = 0,
            pAdjustMethod = "BH"
        )
        
    
        if (is.null(result_gsea_db) || nrow(result_gsea_db@result) == 0) {
            cat(sprintf("[WARNING] No significant pathways identified for %s.\n", db_name))
            return() 
        }

        result_gsea_db_raw_df = result_gsea_db@result %>% 
            dplyr::arrange(desc(NES)) %>%
            dplyr::select(-log2err) %>% 
            dplyr::distinct()
        
        
        output_dir= here("src","Enrichment","Drug",disease,db_name)
        create_dir(output_dir)

        raw_filename = sprintf("%s_gsea_%s_raw.csv",disease, db_name)
        write.csv(result_gsea_db_raw_df,file=here(output_dir,raw_filename), row.names = FALSE)
        cat(sprintf("[INFO] Matriz Primária Extraída: %d rotas mapeadas. (Arquivo: %s)\n", nrow(result_gsea_db_raw_df), raw_filename))


        vias_reprovadas <- result_gsea_db_raw_df %>%
            dplyr::filter(p.adjust >= 0.05 | enrichmentScore <= 0) %>%
            dplyr::select(ID, enrichmentScore, NES, p.adjust) %>%
            dplyr::arrange(enrichmentScore)
        
        cat(sprintf("[INFO] Filtered out %d pathways (FDR >= 0.05 or NES <= 0)\n",nrow(result_gsea_db_raw_df)))
        result_gsea_db_filter_df = result_gsea_db_raw_df %>% 
            dplyr::filter(p.adjust < 0.05, enrichmentScore > 0) %>% 
            dplyr::arrange(desc(NES))
            
        filter_filename = sprintf("%s_gsea_%s_filter.csv", disease,db_name)
        write.csv(result_gsea_db_filter_df, file=here(output_dir,filter_filename), row.names = FALSE)
        cat(sprintf("[INFO] Saving results to: %s\n",filter_filename))
    })

    
    cat("[INFO] GSEA enrichment pipeline completed successfully.\n")

}



plot_go_disease_comparison = function() {
    cat("\n[INFO] Loading the object with the 3 variables (MDD, BD, Commons)...\n")
    ora_go = load_rdata(here("src", "Enrichment", "Comparison", "GeneOntology", "RData", "disease_comparison_go_simplify.RData"))
    
    cat("[INFO] Generating the comparative dotplot...\n")
    
    graph = enrichplot::dotplot(ora_go, showCategory = 5) + 
        ggplot2::theme_bw() +
        ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = "black"),
            axis.text.y = ggplot2::element_text(size = 9, color = "black"),
            axis.title = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 15))
        ) +
        ggplot2::ggtitle("Comparative Analysis: MDD vs BD vs Commons")
    
    output_graph_dir = here("src", "Enrichment", "Comparison","GeneOntology", "Graph")
    create_dir(output_graph_dir)
    
    pdf_path = here(output_graph_dir, "Comparative_Pathway_Analysis_MDD_BD_Commons.pdf")
    cat(sprintf("[INFO] Saving comparative dotplot to:\n%s\n", pdf_path))
    

    grDevices::pdf(file = pdf_path, width = 10, height = 7)
    print(graph)
    grDevices::dev.off()
    cat("[INFO] PDF successfully exported.\n")
    return(invisible(graph))
}
plot_kegg_disease_comparison = function() {
    cat("\n[INFO] Loading the object with the 3 variables (MDD, BD, Commons)...\n")
    ora_kegg = load_rdata(here("src", "Enrichment", "Comparison", "KEGG", "RData", "disease_comparison_kegg.RData"))
    
    cat("[INFO] Generating the comparative dotplot...\n")
    
    graph = enrichplot::dotplot(ora_kegg, showCategory = 5) + 
        ggplot2::theme_bw() +
        ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = "black"),
            axis.text.y = ggplot2::element_text(size = 9, color = "black"),
            axis.title = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 15))
        ) +
        ggplot2::ggtitle("Comparative KEGG Analysis: MDD vs BD vs Commons")
    
    output_graph_dir = here("src", "Enrichment", "Comparison", "KEGG", "Graph")
    create_dir(output_graph_dir)
    
    pdf_path = here(output_graph_dir, "Comparative_KEGG_Analysis_MDD_BD_Commons.pdf")
    cat(sprintf("[INFO] Saving comparative dotplot to:\n%s\n", pdf_path))
    
    grDevices::pdf(file = pdf_path, width = 10, height = 7)
    print(graph)
    grDevices::dev.off()
    
    cat("[INFO] PDF successfully exported.\n")
    return(invisible(graph))
}

plot_reactome_disease_comparison = function() {
    cat("\n[INFO] Loading the object with the 3 variables (MDD, BD, Commons)...\n")
    ora_reactome = load_rdata(here("src", "Enrichment", "Comparison", "Reactome", "RData", "disease_comparison_reactome.RData"))
    
    cat("[INFO] Generating the comparative dotplot...\n")
    
    require(enrichplot, quietly = TRUE)
    
    graph = enrichplot::dotplot(ora_reactome, showCategory = 5) + 
        ggplot2::theme_bw() +
        ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = "black"),
            axis.text.y = ggplot2::element_text(size = 9, color = "black"),
            axis.title = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 15))
        ) +
        ggplot2::ggtitle("Comparative Reactome Analysis: MDD vs BD vs Commons")
    
    output_graph_dir = here("src", "Enrichment", "Comparison", "Reactome", "Graph")
    create_dir(output_graph_dir)
    
    pdf_path = here(output_graph_dir, "Comparative_Reactome_Analysis_MDD_BD_Commons.pdf")
    cat(sprintf("[INFO] Saving comparative dotplot to:\n%s\n", pdf_path))
    
    grDevices::pdf(file = pdf_path, width = 10, height = 7)
    print(graph)
    grDevices::dev.off()
    
    cat("[INFO] PDF successfully exported.\n")
    return(invisible(graph))
}

gsea_mdd_go = function() {
    set.seed(42) 
    output_dir = here("src", "Enrichment",'GSEA','GO','MDD')
    output_dir_rdata = here("src", "Enrichment",'GSEA','GO','MDD','RData')
    create_dir(output_dir_rdata)

    mdd_disease_graph = get_disease_subgraph('MDD')

    between_calc = igraph::betweenness(
        mdd_disease_graph,
        directed = FALSE
    )

    between_scores = between_calc
    names(between_scores) = igraph::V(mdd_disease_graph)$name

    geneList_mdd = sort(between_scores, decreasing = TRUE)
    
    cat("[INFO] Executando Enriquecimento Funcional Ranqueado...\n")

    resultado_gsea_go = clusterProfiler::gseGO(
        geneList      = geneList_mdd, 
        ont           = "BP", 
        OrgDb         = org.Hs.eg.db, 
        keytype       = "ENTREZID",
        exponent      = 1, 
        eps           = 0,            
        maxGSSize     = 500, 
        pvalueCutoff  = 0.05, 
        pAdjustMethod = "BH", 
        verbose       = FALSE, 
        seed          = TRUE
    )
    

    if (is.null(resultado_gsea_go) || nrow(as.data.frame(resultado_gsea_go)) == 0) {
        warning("Nenhuma via biológica significativa foi encontrada ou ocorreu falha de mapeamento de IDs.")
        return(NULL)
    }

    save(resultado_gsea_go, file = here(output_dir_rdata,'gsea_betweenness.RData'))
    
    tabela_mdd_limpa = as.data.frame(resultado_gsea_go) %>%
        dplyr::filter(!is.na(ID))
        
    write.table(
        tabela_mdd_limpa, 
        file = here(output_dir,'gsea_betweenness.csv'), 
        sep = ";",         
        dec = ",",         
        row.names = FALSE, 
        qmethod = "double"
    )
    cat("[INFO] Análise concluída com sucesso.\n")


    return(resultado_gsea_go)
}
gsea_mdd_kegg = function() {
    set.seed(42) 
    mdd_disease_graph = get_disease_subgraph('MDD')

    between_calc = igraph::betweenness(
        mdd_disease_graph,
        directed = FALSE
    )

    between_scores = between_calc
    names(between_scores) = igraph::V(mdd_disease_graph)$name



    geneList_mdd = sort(between_scores, decreasing = TRUE)
    
    cat("[INFO] Executando Enriquecimento Funcional Ranqueado...\n")

    resultado_gsea_kegg = clusterProfiler::gseKEGG(geneList=geneList_mdd, 
        organism = "hsa", 
        keyType = "kegg", 
        exponent = 1, 
        nPerm = 1000, 
        maxGSSize = 500, 
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH", 
        verbose = TRUE,
        se_internal_data = FALSE, 
        seed = TRUE)

    

    if (is.null(resultado_gsea_kegg) || nrow(as.data.frame(resultado_gsea_kegg)) == 0) {
        warning("[ALERTA] Nenhuma via biológica significativa foi encontrada ou ocorreu falha de mapeamento de IDs.")
        return(NULL)
    }
    output_dir = here("src", "Enrichment",'GSEA','KEGG','MDD')
    output_dir_rdata = here("src", "Enrichment",'GSEA','KEGG','MDD','RData')
    create_dir(output_dir_rdata)

    save(resultado_gsea_kegg, file = here(output_dir_rdata,'gsea_betweenness.RData'))
    
    tabela_mdd_limpa_kegg = as.data.frame(resultado_gsea_kegg) %>%
        dplyr::filter(!is.na(ID))
        
    write.table(
        tabela_mdd_limpa_kegg, 
        file = here(output_dir,'gsea_kegg_betweenness.csv'), 
        sep = ";",         
        dec = ",",         
        row.names = FALSE, 
        qmethod = "double"
    )
    cat("[INFO] Análise concluída com sucesso.\n")


    return(resultado_gsea_kegg)
}

gsea_mdd_reactome = function() {
    set.seed(42) 
    mdd_disease_graph = get_disease_subgraph('MDD')

    between_calc = igraph::betweenness(
        mdd_disease_graph,
        directed = FALSE
    )

    between_scores = between_calc
    names(between_scores) = igraph::V(mdd_disease_graph)$name

    geneList_mdd = sort(between_scores, decreasing = TRUE)
    
    cat("[INFO] Executando Enriquecimento Funcional Ranqueado...\n")

    resultado_reactome_kegg = ReactomePA::gsePathway(
        geneList=geneList_mdd, 
        organism = "human", 
        exponent = 1, 
        nPerm = 1000, 
        maxGSSize = 500, 
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH", 
        verbose = TRUE)


    

    if (is.null(resultado_reactome_kegg) || nrow(as.data.frame(resultado_reactome_kegg)) == 0) {
        warning("[ALERTA] Nenhuma via biológica significativa foi encontrada ou ocorreu falha de mapeamento de IDs.")
        return(NULL)
    }
    output_dir = here("src", "Enrichment",'GSEA','REACTOME')
    output_dir_rdata = here(output_dir,'RData')
    create_dir(output_dir_rdata)

    save(resultado_reactome_kegg, file = here(output_dir_rdata,'gsea_betweenness.RData'))
    
    tabela_mdd_limpa_kegg = as.data.frame(resultado_reactome_kegg) %>%
        dplyr::filter(!is.na(ID))
        
    write.table(
        tabela_mdd_limpa_kegg, 
        file = here(output_dir,'gsea_reactome_betweenness.csv'), 
        sep = ";",         
        dec = ",",         
        row.names = FALSE, 
        qmethod = "double"
    )
    cat("[INFO] Análise concluída com sucesso.\n")


    return(resultado_reactome_kegg)
}
gsea_bd_go = function() {
    set.seed(42) 

    bd_disease_graph = get_disease_subgraph('BD')

    between_calc = igraph::betweenness(
        bd_disease_graph,
        directed = FALSE
    )

    between_scores = between_calc
    names(between_scores) = igraph::V(bd_disease_graph)$name

    geneList_bd = sort(between_scores, decreasing = TRUE)
    
    cat("[INFO] Executando Enriquecimento Funcional Ranqueado...\n")

    resultado_gsea_go = clusterProfiler::gseGO(
        geneList      = geneList_bd, 
        ont           = "BP", 
        OrgDb         = org.Hs.eg.db, 
        keytype       = "ENTREZID",
        exponent      = 1, 
        eps           = 0,            
        maxGSSize     = 500, 
        pvalueCutoff  = 0.05, 
        pAdjustMethod = "BH", 
        verbose       = FALSE, 
        seed          = TRUE
    )
    

    if (is.null(resultado_gsea_go) || nrow(as.data.frame(resultado_gsea_go)) == 0) {
        warning("Nenhuma via biológica significativa foi encontrada ou ocorreu falha de mapeamento de IDs.")
        return(NULL)
    }
    output_dir = here("src", "Enrichment",'GSEA','GO','BD')
    output_dir_rdata = here(output_dir,'RData')
    create_dir(output_dir_rdata)


    save(resultado_gsea_go, file = here(output_dir_rdata,'gsea_betweenness.RData'))
    
    tabela_bd_limpa = as.data.frame(resultado_gsea_go) %>%
        dplyr::filter(!is.na(ID))
        
    write.table(
        tabela_bd_limpa, 
        file = here(output_dir,'gsea_betweenness.csv'), 
        sep = ";",         
        dec = ",",         
        row.names = FALSE, 
        qmethod = "double"
    )
    cat("[INFO] Análise concluída com sucesso.\n")


    return(resultado_gsea_go)
}
gsea_bd_kegg = function() {
    set.seed(42) 
    bd_disease_graph = get_disease_subgraph('BD')

    between_calc = igraph::betweenness(
        bd_disease_graph,
        directed = FALSE
    )

    between_scores = between_calc
    names(between_scores) = igraph::V(bd_disease_graph)$name



    geneList_bd = sort(between_scores, decreasing = TRUE)
    
    cat("[INFO] Executando Enriquecimento Funcional Ranqueado...\n")

    resultado_gsea_kegg = clusterProfiler::gseKEGG(geneList=geneList_bd, 
        organism = "hsa", 
        keyType = "kegg", 
        exponent = 1, 
        nPerm = 1000, 
        maxGSSize = 500, 
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH", 
        verbose = TRUE,
        se_internal_data = FALSE, 
        seed = TRUE)

    

    if (is.null(resultado_gsea_kegg) || nrow(as.data.frame(resultado_gsea_kegg)) == 0) {
        warning("[ALERTA] Nenhuma via biológica significativa foi encontrada ou ocorreu falha de mapeamento de IDs.")
        return(NULL)
    }
    output_dir = here("src", "Enrichment",'GSEA','KEGG','BD')
    output_dir_rdata = here(output_dir,'RData')
    create_dir(output_dir_rdata)

    save(resultado_gsea_kegg, file = here(output_dir_rdata,'gsea_betweenness.RData'))
    
    tabela_bd_limpa_kegg = as.data.frame(resultado_gsea_kegg) %>%
        dplyr::filter(!is.na(ID))
        
    write.table(
        tabela_bd_limpa_kegg, 
        file = here(output_dir,'gsea_kegg_betweenness.csv'), 
        sep = ";",         
        dec = ",",         
        row.names = FALSE, 
        qmethod = "double"
    )
    cat("[INFO] Análise concluída com sucesso.\n")


    return(resultado_gsea_kegg)
}

gsea_bd_reactome = function() {
    set.seed(42) 
    bd_disease_graph = get_disease_subgraph('BD')

    between_calc = igraph::betweenness(
        bd_disease_graph,
        directed = FALSE
    )

    between_scores = between_calc
    names(between_scores) = igraph::V(bd_disease_graph)$name

    geneList_bd = sort(between_scores, decreasing = TRUE)
    
    cat("[INFO] Executando Enriquecimento Funcional Ranqueado...\n")

    resultado_reactome = ReactomePA::gsePathway(
        geneList=geneList_bd, 
        organism = "human", 
        exponent = 1, 
        nPerm = 1000, 
        maxGSSize = 500, 
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH", 
        verbose = TRUE)


    

    if (is.null(resultado_reactome) || nrow(as.data.frame(resultado_reactome)) == 0) {
        warning("[ALERTA] Nenhuma via biológica significativa foi encontrada ou ocorreu falha de mapeamento de IDs.")
        return(NULL)
    }
    output_dir = here("src", "Enrichment",'GSEA','REACTOME')
    output_dir_rdata = here(output_dir,'RData')
    create_dir(output_dir_rdata)

    save(resultado_reactome, file = here(output_dir_rdata,'gsea_betweenness.RData'))
    
    tabela_bd_limpa_kegg = as.data.frame(resultado_reactome) %>%
        dplyr::filter(!is.na(ID))
        
    write.table(
        tabela_bd_limpa_kegg, 
        file = here(output_dir,'gsea_reactome_betweenness.csv'), 
        sep = ";",         
        dec = ",",         
        row.names = FALSE, 
        qmethod = "double"
    )
    cat("[INFO] Análise concluída com sucesso.\n")


    return(resultado_reactome)
}

map_ranked_drug_targets_to_symbols= function(ranked_drug_targets){
    cat(sprintf("[INFO] Ranked drugs: %d\n",length(ranked_drug_targets)))
   

    drug_target_df = tibble::tibble(
            drugbank_id = names(ranked_drug_targets),
            entrez_id   = ranked_drug_targets) %>%
        tidyr::unnest(entrez_id)

    cat(sprintf("[INFO] Drug-target interactions: %d\n",nrow(drug_target_df)))

    gene_symbols = AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys=unique(as.character(drug_target_df$entrez_id)),
        column  = "SYMBOL",
        keytype = "ENTREZID")

    drug_target_df = drug_target_df %>%
        dplyr::mutate(symbol = gene_symbols[as.character(entrez_id)]) %>%
        dplyr::filter(!is.na(symbol)) %>%
        dplyr::distinct(drugbank_id,entrez_id,symbol)
    
    cat(sprintf("[INFO] Expanded %d ranked drugs into %d drug-target unique interactions.\n",length(ranked_drug_targets),nrow(drug_target_df)))


    return(drug_target_df)

}

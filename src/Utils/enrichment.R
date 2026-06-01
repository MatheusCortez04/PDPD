library(here)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
source(here("src","Utils","utils.R"))
source(here("src","Utils","graph_metrics.R"))

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

enrichment_by_rank_drugs_by_disease = function(disease = c("MDD", "BD")) {
    disease = match.arg(disease)
    cat(sprintf("\n[INFO] Starting enrichment analysis for disease: %s\n", disease))
    
    score_filter = get_score_disease_gene_association()
    
    rank_dir = here("src", "Evaluation", disease)
    rank_file_name = c("MDD" = sprintf("prediction_mdd_score_filter_%s_.csv", score_filter), 
                       "BD" = sprintf("prediction_bipolar_score_filter_%s_.csv", score_filter))
    
    cat("[INFO] Loading the top 100 ranked drugs...\n")
    rank_drug_disease = read.csv(here(rank_dir, rank_file_name[disease])) %>% 
        dplyr::slice_head(n = 100)

    cat("[INFO] Importing drug targets mapping...\n")
    drug_target_df  = import_drug_targets_df()
    
    cat("[INFO] Extracting specific targets for the ranked drugs...\n")
    drug_targets = purrr::map(rank_drug_disease$drugbank_id, extract_targets_per_drug, drug_target_df)
    names(drug_targets) = rank_drug_disease$drugbank_id

    cat("[INFO] Running Gene Ontology (GO) pathway comparison...\n")
    comparate_cluster_go(drug_targets, sprintf("drug_pathways_rank_%s", disease))
    
    cat("[INFO] Running KEGG pathway comparison...\n")
    comparate_cluster_kegg(drug_targets, sprintf("drug_pathways_rank_%s", disease))
    
    cat("[INFO] Running Reactome pathway comparison...\n")
    comparate_cluster_reactome(drug_targets, sprintf("drug_pathways_rank_%s", disease))
    
    cat("[INFO] Enrichment analysis completed successfully!\n")
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
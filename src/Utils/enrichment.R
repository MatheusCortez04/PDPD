library(here)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
source(here("src","Utils","utils.R"))

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

    specific_gene_file_prefix = ifelse(disease == "MDD","mdd_specific_genes","bipolar_specific_genes")

    disease_genes = extract_genes_by_disease(disease = disease) %>%dplyr::pull(entrez_id)

    cat(sprintf("[INFO] Total disease-associated genes: %d\n",nrow(disease_genes)))

    calculate_go_enrichment(entrez_ids= disease_genes,output_file_prefix = disease_file_prefix)
    calculate_kegg_enrichment(entrez_ids= disease_genes,output_file_prefix = disease_file_prefix)
    calculate_go_enrichment(entrez_ids= disease_specific_genes,output_file_prefix = specific_gene_file_prefix)
    calculate_kegg_enrichment(entrez_ids= disease_specific_genes,output_file_prefix = specific_gene_file_prefix)

    cat(sprintf("[INFO] Enrichment pipeline completed for disease: %s\n",disease))

    return(invisible(NULL))
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

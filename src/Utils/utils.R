library(here)

clear_console = function(){
    system("clear")
}

is_valid_input_boolean = function(input){
    return(toupper(input) == 'TRUE' || toupper(input) == 'FALSE')
}


get_ppi_nodes = function(ppi_df){
    proteinA_entrezid = ppi_df$proteinA_entrezid
    proteinB_entrezid = ppi_df$proteinB_entrezid
    all_protein_in_ppi_df = c(proteinA_entrezid,proteinB_entrezid)
    print(paste("All proteins in PPI Dataframe: ",length(all_protein_in_ppi_df)))
    unique_ordered_proteins_in_ppi_df = unique(all_protein_in_ppi_df)
    print(paste("Unique proteins in PPI Dataframe: ",length(unique_ordered_proteins_in_ppi_df)))
    invisible(unique_ordered_proteins_in_ppi_df)
}
get_drug_nodes = function(drug_target_df){
    drug_ids = drug_target_df$drugbank_id
    print(paste("All drugs in Drug to target Dataframe: ",length(drug_ids)))
    unique_drugs_in_df = unique(drug_ids)
    print(paste("Unique drugs in Drug to target Dataframe: ",length(unique_drugs_in_df)))
    invisible(unique_drugs_in_df)
}

load_rdata <- function(path_file) {
  env <- new.env()
  load(path_file, envir = env)
  objs <- ls(env)
  if(length(objs) != 1) {
    stop("Mais de um objeto no arquivo RData")
  }
  return(env[[objs]])
}
get_mdd_genes = function(){
  major_depressive_disorder_id = "C1269683"
  score_filter = 0.6
  disease_gene_df  = read.csv(here("src","Data","disease_genes.csv"), sep="\t")
  disease_gene_df =disease_gene_df %>%
    filter(diseaseid ==major_depressive_disorder_id &
     score>=score_filter) %>%
     rename(gene_id=geneid,disease_id=diseaseid) %>% 
     select(gene_id,disease_id,score)
  invisible(disease_gene_df)

}
get_bipolar_disorder_genes = function(){
  bipolar_disorder_disorder_id = "C0005586"
  score_filter = 0.6
  disease_gene_df  = read.csv(here("src","Data","disease_genes.csv"), sep="\t")
  disease_gene_df =disease_gene_df %>%
    filter(diseaseid ==major_depressivve_disorder_id &
     score>=score_filter) %>%
     rename(gene_id=geneid,disease_id=diseaseid) %>% 
     select(gene_id,disease_id,score)
  invisible(disease_gene_df)

}



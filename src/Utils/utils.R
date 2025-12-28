library(here)
library(dplyr)
clear_console = function(){
    system("clear")
}

is_valid_input_boolean = function(input){
    return(toupper(input) == 'TRUE' || toupper(input) == 'FALSE')
}


get_ppi_nodes = function(){
    ppi_df  = read.csv(here("src","Data","PPI_gysi.csv"), sep=",")
    cat("PPI  file reading complete.\n")
    proteinA_entrezid = ppi_df$proteinA_entrezid
    proteinB_entrezid = ppi_df$proteinB_entrezid
    all_protein_in_ppi_df = c(proteinA_entrezid,proteinB_entrezid)
    print(paste("All proteins in PPI Dataframe: ",length(all_protein_in_ppi_df)))
    unique_ordered_proteins_in_ppi_df = unique(all_protein_in_ppi_df)
    print(paste("Unique proteins in PPI Dataframe: ",length(unique_ordered_proteins_in_ppi_df)))
    invisible(unique_ordered_proteins_in_ppi_df)
}
get_drug_nodes = function(){
    drug_target_df  = read.csv(here("src","Data","Drug","drug_targets_DrugBank_Gysi.csv"), sep=",")
    cat("Drug Target file reading complete.\n")
    drug_ids = drug_target_df$drugbank_id
    print(paste("All drugs in Drug to target Dataframe: ",length(drug_ids)))
    unique_drugs_in_df = unique(drug_ids)
    print(paste("Unique drugs in Drug to target Dataframe: ",length(unique_drugs_in_df)))
    invisible(unique_drugs_in_df)
}

load_rdata <- function(path_file) {
  cat("Loading RData of file :", path_file, "\n")
  env <- new.env()
  load(path_file, envir = env)
  objs <- ls(env)
  return(env[[objs]])
}
get_mdd_genes = function(){
  major_depressive_disorder_id = "C1269683"
  score_filter = 0.9
  disease_gene_df  = read.csv(here("src","Data","Disease","disease_genes.csv"), sep="\t")
  disease_gene_df =disease_gene_df %>%
    filter(diseaseid ==major_depressive_disorder_id &
     score>=score_filter) %>%
     rename(gene_id=geneid,disease_id=diseaseid) %>% 
     select(gene_id,disease_id,score)
  invisible(disease_gene_df)

}
get_bipolar_disorder_genes = function(){
  bipolar_disorder_id = "C0005586"
  score_filter =  0.9
  disease_gene_df  = read.csv(here("src","Data","Disease","disease_genes.csv"), sep="\t")
  disease_gene_df =disease_gene_df %>%
    filter(diseaseid ==bipolar_disorder_id &
     score>=score_filter) %>%
     rename(gene_id=geneid,disease_id=diseaseid) %>% 
     select(gene_id,disease_id,score)
  invisible(disease_gene_df)

}

build_protein_mdd_df= function(){
  all_proteins = get_ppi_nodes()
  mdd_genes = get_mdd_genes()
  cat("Creating MDD protein DataFrame....")
  mdd_vector = data.frame(gene_id=all_proteins)
  mdd_vector <- mdd_vector %>% mutate(
   is_disease = ifelse(gene_id %in% mdd_genes$gene_id, 1, 0))


  output_file_name ="mdd_genes_vector"
  output_path= here("src","Data","Disease")
  output_file_path_csv = here(output_path,paste0(output_file_name,".csv"))

  output_path_rdata =here(output_path,"Rdata")
  dir.create(output_path_rdata, recursive = TRUE, showWarnings = FALSE)
  output_file_path_rdata = here(output_path_rdata,paste0(output_file_name,".Rdata"))
  
  write.csv(mdd_vector,file=output_file_path_csv,row.names=FALSE )
  save(mdd_vector, file = output_file_path_rdata)
  invisible(mdd_vector)
}


build_protein_bipolar_df= function(){
  all_proteins = get_ppi_nodes()
  bipolar_genes = get_bipolar_disorder_genes()
  cat("Creating BD protein DataFrame....")
  bipolar_vector = data.frame(gene_id=all_proteins)
  bipolar_vector <- bipolar_vector %>% mutate(
   is_disease = ifelse(gene_id %in% bipolar_genes$gene_id, 1, 0))
  output_file_name ="bipolar_genes_vector"
  output_path= here("src","Data","Disease")
  output_file_path_csv = here(output_path,paste0(output_file_name,".csv"))

  output_path_rdata =here(output_path,"Rdata")
  dir.create(output_path_rdata, recursive = TRUE, showWarnings = FALSE)

  output_file_path_rdata = here(output_path_rdata,paste0(output_file_name,".Rdata"))
  
  write.csv(bipolar_vector,file=output_file_path_csv,row.names=FALSE)
  save(bipolar_vector, file = output_file_path_rdata)
  invisible(bipolar_vector)
}

load_drug_target_df = function(){
    drug_target_file_path = here("src","Data","Drug","drug_targets_DrugBank_Gysi.csv")
    drug_target_file_already_exists = file.exists(drug_target_file_path)
    if(!drug_target_file_already_exists){
        cat("\n[Error]: Required file Drug target file not found.\n This file is necessary to calculate the score. Would you like to create it now\n")
        Sys.sleep(1.5)
        return()
    }
    drug_target_df = read.csv(drug_target_file_path)
    drug_target_df %>% dplyr::distinct(drugbank_id, entrez_id) %>%
        dplyr::mutate(entrez_id = as.character(entrez_id))
    invisible(drug_target_df)
}
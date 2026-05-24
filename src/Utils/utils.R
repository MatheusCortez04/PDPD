library(here)
library(dplyr)
clear_console = function(){
    system("clear")
}

is_valid_input_boolean = function(input){
    return(toupper(input) == 'TRUE' || toupper(input) == 'FALSE')
}

extract_ppi_genes = function(){
    ppi_df  = import_ppi_interactions()
    proteinA_entrezid = ppi_df$proteinA_entrezid
    proteinB_entrezid = ppi_df$proteinB_entrezid
    all_protein_in_ppi_df = c(proteinA_entrezid,proteinB_entrezid)
    cat(sprintf("[INFO] All proteins in PPI Dataframe: %d\n",length(all_protein_in_ppi_df)))
    unique_ordered_proteins_in_ppi_df = unique(all_protein_in_ppi_df)
    cat(sprintf("[INFO] Unique proteins in PPI Dataframe: %d\n",length(unique_ordered_proteins_in_ppi_df)))
    invisible(unique_ordered_proteins_in_ppi_df)
}

load_rdata <- function(path_file) {
  cat("\n Loading RData of file :", path_file, "\n")
  env <- new.env()
  load(path_file, envir = env)
  objs <- ls(env)
  return(env[[objs]])
}
get_mdd_disease_module = function() {
  
  mdd_disease_id = "C1269683"
  score_filter = get_score_disease_gene_association()
  cat("\n[INFO] Building MDD disease module...\n")


  disease_gene_df  = read.csv(here("src","Data","Disease","disease_genes.csv"))
  ppi_nodes = extract_ppi_genes()
  mdd_gene_module = disease_gene_df %>%
    filter(diseaseid == mdd_disease_id & score >= score_filter) %>%
    rename(entrez_id = geneid, disease_id = diseaseid) %>% 
    mutate(entrez_id = as.character(entrez_id))

  cat(sprintf("[INFO] MDD genes after score filter: %d\n", nrow(mdd_gene_module)))

  mdd_gene_module=mdd_gene_module %>%
    filter(entrez_id %in% ppi_nodes)%>%
    select(entrez_id, disease_id, score)
  
  cat(sprintf("[INFO] MDD genes after PPI filter: %d\n", nrow(mdd_gene_module)))

  invisible(mdd_gene_module)
}
get_bipolar_disease_module = function() {
   cat("\n[INFO] Building Bipolar Disease disease module...\n")
  bipolar_disease_id = "C0005586"
  score_filter = get_score_disease_gene_association()

  
  ppi_nodes = extract_ppi_genes()
  disease_gene_df  = read.csv(here("src","Data","Disease","disease_genes.csv"))

 bipolar_gene_module= disease_gene_df %>%
    filter(diseaseid == bipolar_disease_id & score >= score_filter) %>%
    rename(entrez_id = geneid, disease_id = diseaseid) %>% 
    mutate(entrez_id = as.character(entrez_id))
  
  cat(sprintf("[INFO] Bipolar Disease genes after score filter: %d\n", nrow(bipolar_gene_module)))

   
  bipolar_gene_module=bipolar_gene_module %>%  
    filter(entrez_id %in% ppi_nodes)%>%
    select(entrez_id, disease_id, score)
  
  
  cat(sprintf("[INFO] Bipolar Disease genes after PPI filter: %d\n", nrow(bipolar_gene_module)))

  invisible(bipolar_gene_module)
}
create_dir <- function(path) {
  if (!dir.exists(path)) {
    message(paste("Creating Dir:", path))
    dir.create(path, recursive = TRUE, mode = "0777")
  }
}
get_score_disease_gene_association = function(){
  score_filter =  0.6
  invisible(score_filter)
}

get_common_gene = function(){
  mdd_disease_genes =get_mdd_disease_module()
  bipolar_disease_genes = get_bipolar_disease_module()
  common_genes = mdd_disease_genes %>% dplyr::semi_join(bipolar_disease_genes,by='entrez_id') %>% 
    distinct(entrez_id) %>%
    mutate(entrez_id= as.character(entrez_id)) %>% pull(entrez_id)
  return(common_genes)
}
import_ppi_interactions = function(){
    file_path=here("src","Data","PPI_gysi.csv")
    if(!file.exists(file_path)){
        stop(sprintf("[ERROR] PPI file not found:\n%s",file_path))
    }
    ppi_df = read.csv(file_path,stringsAsFactors = FALSE) %>%
      dplyr::distinct()

    return(ppi_df)
}

extract_ppi_genes = function(){
  ppi_df =import_ppi_interactions()
  ppi_genes = c(ppi_df$proteinA_entrezid,ppi_df$proteinB_entrezid) %>%unique()

  return(ppi_genes)
}
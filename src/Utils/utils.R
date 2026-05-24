library(here)

clear_console = function(){
    system("clear")
}

is_valid_input_boolean = function(input){
    return(toupper(input) == 'TRUE' || toupper(input) == 'FALSE')
}
load_rdata <- function(path_file) {
  cat("Loading RData of file :", path_file, "\n")
  env <- new.env()
  load(path_file, envir = env)
  objs <- ls(env)
  return(env[[objs]])
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

get_disease_genes = function(disease = c("MDD", "BD")) {
  disease = match.arg(disease)

  disease_ids = c(
    "MDD" = "C1269683", 
    "BD"  = "C0005586"
  )
  target_disease_id = disease_ids[disease]
  valid_genes = extract_ppi_genes()
  score_filter = get_score_disease_gene_association()
  disease_gene_df = read.csv(here("src", "Data", "Disease", "disease_genes.csv"))
  
  disease_gene_df = disease_gene_df %>%
    dplyr::filter(diseaseid == target_disease_id & score >= score_filter,geneid %in%valid_genes) %>%
    dplyr::rename(entrez_id = geneid, disease_id = diseaseid) %>% 
    dplyr::select(entrez_id, disease_id, score)
  
  invisible(disease_gene_df)
}

build_protein_disease_df = function(disease = c("MDD", "BD")) {
  disease = match.arg(disease)
  
  all_proteins = extract_ppi_genes()
  disease_genes = get_disease_genes(disease)
  

  cat(paste0("Creating ", disease, " protein DataFrame....\n"))
  
  disease_vector = data.frame(entrez_id = all_proteins)
  disease_vector = disease_vector %>% 
    dplyr::mutate(
      is_disease = ifelse(entrez_id %in% disease_genes$entrez_id, 1, 0)
    )

  file_prefix = ifelse(disease == "MDD", "mdd", "bipolar")
  score_filter = get_score_disease_gene_association()
  output_file_name = paste0(file_prefix, "_genes_vector_score_filter_",score_filter)
  
  output_path = here("src", "Data", "Disease")
  output_file_path_csv = here(output_path, paste0(output_file_name, "_.csv"))

  output_path_rdata = here(output_path, "Rdata")
  create_dir(output_path_rdata)
  output_file_path_rdata = here(output_path_rdata, paste0(output_file_name, "_.Rdata"))
  
  write.csv(disease_vector, file = output_file_path_csv, row.names = FALSE)
  save(disease_vector, file = output_file_path_rdata)
  
  invisible(disease_vector)
}

get_common_gene = function(){
  mdd_disease_genes =get_disease_genes("MDD")
  bipolar_disease_genes = get_disease_genes("BD")
  common_genes = mdd_disease_genes %>% dplyr::semi_join(bipolar_disease_genes,by='entrez_id') %>% 
    distinct(entrez_id) %>%
    mutate(entrez_id= as.character(entrez_id)) %>% pull(entrez_id)
  return(common_genes)
}

import_ppi_interactions = function(){
    file_path=here("src","Data","PPI_gysi.csv")
    if(!file.exists(file_path)){
      stop(sprintf("[ERROR] Required PPI file not found:\n%s
        \nPlease verify the file path or generate the interaction dataset before running the pipeline.",
        file_path))
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
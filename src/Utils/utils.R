library(here)

clear_console = function(){
    system("clear")
}

is_valid_input_boolean = function(input){
    return(toupper(input) == 'TRUE' || toupper(input) == 'FALSE')
}


get_ppi_nodes = function(){
    ppi_df  = get_ppi_dataframe()
    cat("PPI  file reading complete.\n")
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
  valid_genes = get_ppi_nodes()
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
  
  all_proteins = get_ppi_nodes()
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


get_ppi_dataframe = function(){
  ppi_df = read.csv(here("src","Data","PPI_gysi.csv"), sep=",") %>% distinct()

  return(ppi_df)
}

get_common_gene = function(){
  mdd_disease_genes =get_disease_genes("MDD")
  bipolar_disease_genes = get_disease_genes("BD")
  common_genes = mdd_disease_genes %>% dplyr::semi_join(bipolar_disease_genes,by='entrez_id') %>% 
    distinct(entrez_id) %>%
    mutate(entrez_id= as.character(entrez_id)) %>% pull(entrez_id)
  return(common_genes)
}
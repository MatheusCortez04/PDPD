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
  score_filter = 0.6
  cat("\n[INFO] Building MDD disease module...\n")


  disease_gene_df  = read.csv(here("src","Data","Disease","disease_genes.csv"), sep="\t")
  ppi_nodes = get_ppi_nodes()
  mdd_gene_module = disease_gene_df %>%
    filter(diseaseid == mdd_disease_id & score >= score_filter) %>%
    rename(entrez_id = geneid, disease_id = diseaseid) %>% 
    mutate(entrez_id = as.character(entrez_id)) %>% 

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
  score_filter = 0.6

  
  ppi_nodes = get_ppi_nodes()
  disease_gene_df  = read.csv(here("src","Data","Disease","disease_genes.csv"), sep="\t")

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

get_mdd_disease_module = function() {
  
  mdd_disease_id = "C1269683"
  score_filter = 0.6
  
  cat("\n[INFO] Building MDD disease module...\n")
  
  disease_gene_df  = read.csv(here("src","Data","Disease","disease_genes.csv"), sep="\t")
  cat(sprintf("[INFO] Total gene-disease associations loaded: %d\n", nrow(disease_gene_df)))
  
  ppi_nodes = get_ppi_nodes()
  cat(sprintf("[INFO] Total PPI nodes available: %d\n", length(ppi_nodes)))
  
  # ---- filtro inicial ----
  filtered_df = disease_gene_df %>%
    filter(diseaseid == mdd_disease_id)
  
  cat(sprintf("[INFO] MDD genes before score filter: %d\n", nrow(filtered_df)))
  
  # ---- filtro por score ----
  filtered_df = filtered_df %>%
    filter(score >= score_filter)
  
  cat(sprintf("[INFO] MDD genes after score >= %.2f: %d\n", score_filter, nrow(filtered_df)))
  
  # ---- transformação ----
  mdd_gene_module = filtered_df %>%
    rename(entrez_id = geneid, disease_id = diseaseid) %>% 
    mutate(entrez_id = as.character(entrez_id))
  
  # ---- projeção na PPI ----
  before_ppi <- nrow(mdd_gene_module)
  
  mdd_gene_module = mdd_gene_module %>%
    filter(entrez_id %in% ppi_nodes) %>%
    select(entrez_id, disease_id, score)
  
  after_ppi <- nrow(mdd_gene_module)
  
  cat(sprintf(
    "[INFO] MDD genes after PPI filter: %d (retention: %.2f%%)\n",
    after_ppi,
    100 * after_ppi / before_ppi
  ))
  
  # ---- sanity check ----
  if (after_ppi == 0) {
    cat("[WARN] MDD module is EMPTY after PPI filtering!\n")
  }
  
  invisible(mdd_gene_module)
}
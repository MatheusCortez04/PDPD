generate_ordered_drug_protein_matrix = function(drug_target_df,protein_nodes,drug_nodes){
    drug_target_df %>% 
    filter( entrez_id %in% protein_nodes) %>%
    distinct(drugbank_id, entrez_id) %>%
    mutate(
      drugbank_id = factor(drugbank_id, levels = drug_nodes),
      entrez_id = factor(entrez_id, levels = protein_nodes),
      value=1
    ) %>%   
    pivot_wider(
        id_cols = drugbank_id, 
        names_from = entrez_id,
        values_from = value,
        values_fill = 0,
        names_expand = TRUE,
        id_expand = TRUE
    ) %>% arrange(drugbank_id)

}
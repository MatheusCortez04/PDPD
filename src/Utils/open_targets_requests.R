library(httr2)
library(dplyr)

get_open_target_config = function(){
    return (
        list(
            url='https://api.platform.opentargets.org/api/v4/graphql',
            queries=list(
                drug_info_from_dbid='query search($queryString: String!) {
                        search(queryString: $queryString, entityNames: ["drug"]) {
                            hits {
                                id,
                                name
                            }
                        }
                    }',
                clinical_report='query RecordDetailQuery($clinicalReportId: String!) {
                    clinicalReport(clinicalReportId: $clinicalReportId) {
                        id
                        title
                        clinicalStage
                        source
                        url
                    }
                }'    
            
            )
            
    ))
}


get_drug_info_from_dbid = function (drugbank_id) {

  message(sprintf("[REQUEST] get_drug_info_from_dbid() | drugbank_id=%s", drugbank_id))
  opentarget_config = get_open_target_config()

  endpoint = opentarget_config$url
  graphql_query = opentarget_config$queries$drug_info_from_dbid
  message("[INFO] Sending search request to OpenTargets...")
  response = httr2::request(endpoint) %>%
    httr2::req_body_json(
      list(query = graphql_query,
        variables = list(queryString = drugbank_id))) %>%
    httr2::req_headers("Content-Type" = "application/json") %>%
    httr2::req_perform()

  status=httr2::resp_status(response)
  message(sprintf("[INFO] HTTP status = %s", status))
  data=response %>% httr2::resp_body_json()

  if (status >= 400) {
    message("[ERROR] Failed to retrieve ChEMBL ID from OpenTargets search.")
    return (list(chembl_id = NA, name = NA))
  }

  hits = purrr::pluck(data, "data", "search", "hits") %>% purrr::flatten()
  if (is.null(hits) || length(hits) == 0) {
    message("[WARN] No ChEMBL hits found for this DrugBank ID.")
    return (list(chembl_id = NA, name = NA))
  }
  chembl_id = hits$id
  message(sprintf("[INFO] Found ChEMBL ID(s) for %s: %s", drugbank_id, paste(chembl_id, collapse = ", ")))
  drug = list(
    chembl_id = hits$id,
    drug_name = hits$name
  )
  invisible(drug)
}
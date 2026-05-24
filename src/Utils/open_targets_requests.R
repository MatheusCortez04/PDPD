library(httr2)
library(dplyr)
library(purrr)

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
                drug_candidates='query drugsQuery($efoId: String!) {
                    disease(efoId: $efoId) {
                        name
                            drugAndClinicalCandidates {
                            count
                            rows {
                                id
                                maxClinicalStage
                                drug {
                                    id
                                    name
                                }
                                clinicalReports {
                                    id
                                }
                            }
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
            ),
            disease_ids = list("MDD" = "MONDO_0002009", "BD" = "MONDO_0004985")
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

get_drug_candidates = function (disease = c("MDD", "BD")) {
  disease = match.arg(disease)

  opentarget_config = get_open_target_config()
  clinical_stage_map = c(
    "PRECLINICAL" = 0,
    "PHASE_1" = 1,
    "PHASE_1_2" = 1.5,
    "PHASE_2" = 2,
    "PHASE_2_3" = 2.5,
    "PHASE_3" = 3,
    "PHASE_4" = 4,
    "APPROVED" = 5,
    "WITHDRAWN" = -1
  )
  
  disease_ids = opentarget_config$disease_ids
  endpoint = opentarget_config$url

  graphql_query = opentarget_config$queries$drug_candidates

  response = httr2::request(endpoint) %>%
    httr2::req_body_json(
      list(
        query = graphql_query,
        variables = list(efoId = disease_ids[[disease]])
      )
    ) %>%
    httr2::req_perform()

  data = response %>% resp_body_json()
  rows = purrr::pluck(
    data,
    "data",
    "disease",
    "drugAndClinicalCandidates",
    "rows"
  )
  if (is.null(rows) || length(rows) == 0) {
    return (NA)
  }

  parsed_data = purrr::map_dfr(rows, function (indication_drug) {
      first_report_id = ifelse(length(indication_drug$clinicalReports) > 0,
        indication_drug$clinicalReports[[1]] $id, NA_character_)
      tibble::tibble(
        chembl_id = indication_drug$drug$id %||% NA_character_,
        clinical_report_id = first_report_id,
        max_clinical_stage_open_targets = factor(
          indication_drug$maxClinicalStage,
          levels = c(
            "PRECLINICAL",
            "PHASE_1",
            "PHASE_1_2",
            "PHASE_2",
            "PHASE_2_3",
            "PHASE_3",
            "PHASE_4",
            "APPROVAL"
          )
        ),
      )
    }) %>%
    dplyr::distinct() %>%
    dplyr::arrange(desc(max_clinical_stage_open_targets))
  return (parsed_data)
}

get_clinical_report_data= function(clinicalReportId){
    opentarget_config = get_open_target_config()
    endpoint = opentarget_config$url
    graphql_query = opentarget_config$queries$clinical_report

    response = httr2::request(endpoint) %>%
        httr2::req_body_json(
            list(
                query = graphql_query,
                variables = list(clinicalReportId =clinicalReportId)
            )
        ) %>%
        httr2::req_perform()

    data = response %>% resp_body_json()
    report = data$data$clinicalReport
    if(is.null(report)) return(NULL)
    evidence_summary = sprintf("%s (%s)", report$source, report$url)
  return(tibble(
      clinical_report_id = clinicalReportId,
      source = report$source,
      url = report$url,
      evidence_summary
    ))
}


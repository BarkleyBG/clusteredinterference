

# # ##########################################
# # ###### Arg tests for data quality ########
# # ##########################################
# #
# #' Testing arguments for data quality
# #'
# #' @inheritParams policyFX
# #' @param formule 4-part formula
# #' @param ... dots
# #'
# #' @export
argTests <- function(
  data,
  formule,
  alphas,
  max_k_evals,
  ...
){
  ### ported from old program
  if( "tbl_df" %in% class(data) ) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  if ("group_prop_scores" %in% names(data)){
    stop("no columns in the data may be named 'group_prop_scores'")
  }

  stopifnot(0<= alphas & alphas <=1)
  stopifnot(0<= max_k_evals & is.numeric(max_k_evals))


  alphas <- sort(alphas)
  if (anyDuplicated(alphas)) {
    warning("Duplicate alphas supplied; duplicates will be removed")
    alphas <- alphas[!duplicated(alphas)]
  }

  formule_dims <- length(formule)
  len_lhs_formule <- formule_dims[1]
  len_rhs_formule <- formule_dims[2]
  modeling_formula <- stats::formula(
    stats::terms(formule, lhs = len_lhs_formule, rhs = -2))

  gp_terms <- stats::terms(formule, lhs = 0, rhs = len_rhs_formule)
  grouping_var_name <- attr(gp_terms, 'term.labels')
  grouping_var_vec <- data[[grouping_var_name]]

  if (any( grouping_var_vec != sort(grouping_var_vec))) {
    data <- data[order(grouping_var_vec), ]
  }

  out_list <- list(
    data = data,
    alphas = alphas,
    modeling_formula = modeling_formula,
    grouping_var_name = grouping_var_name
  )
}


# #' Splits Data by Clusters
# #'
# #' This function nicely splits the data up into important components by cluster
# #'
# #' @inheritParams policyFX
# #' @inheritParams argTests
# #' @param grouping_var_name the name of the cluster ID variable
# #' @param glmer_fit the fitted lme4::glmer() model object
# #'
# #' @export
splitData <- function(
  data,
  formule,
  glmer_fit,
  grouping_var_name
){
  len_lhs_formule <- length(formule)[1]
  len_rhs_formule <- length(formule)[2]
  model_matrix <- stats::model.matrix(glmer_fit)
  split_model_matrix <- by(model_matrix, data[[grouping_var_name]],
                           function(x){as.matrix(x)})
  split_otherdata <- by(data, data[[grouping_var_name]], function(group_data){
    outcome_vec <- Formula::model.part(
      formule, data = group_data, lhs = 1, drop = TRUE)
    treatment_vec <- Formula::model.part(
      formule, data = group_data, lhs = 2, drop = TRUE)
    if (len_lhs_formule > 2) {
      model_DV_vec <- Formula::model.part(formule, data = group_data,
                                          lhs = len_lhs_formule, drop = TRUE)
    } else {
      model_DV_vec <- treatment_vec
    }
    grouping_vec <- Formula::model.part(
      formule, data = group_data, rhs = 2, drop = TRUE)
    cluster_name <- grouping_vec[1]
    if(any(grouping_vec != cluster_name)) {
      stop("Error when determining unique cluster identifiers")
    }

    ##output from subfunction
    list(
      outcome_vec = outcome_vec,
      treatment_vec = treatment_vec,
      model_DV_vec = model_DV_vec,
      cluster_id_vec = grouping_vec,
      cluster_name = cluster_name
    )
  })

  split_data_list <-   mapply(
    FUN = function(model_mat, data_list){
      append(data_list, list(model_matrix=model_mat ))
    },
    model_mat = split_model_matrix,
    data_list = split_otherdata,
    SIMPLIFY = FALSE
  )
}


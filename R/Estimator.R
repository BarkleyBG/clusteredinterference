
#' Estimate causal effects (FX) of nonrandomized policies
#'
#' For the dots function, user may supply their own \code{target_grid} through the dots argument. The \code{target_grid} can be made through exported function \code{makeTargetGrid}
#'
#' @param data a dataframe (not a tibble)
#' @param formula a long formula call with outcome | treatment ~ covariates + (1|cluster_id) | cluster_id
#' @param alphas a vector of policies of interest. Each entry must be between 0 and 1
#' @param k_samps The maximum number of vectors to evaluate to estimate the omega term. Setting to 0 avoids approximation at the cost of increased computation time. Recommended to set <= 5.
#' @param ... The dots
#' @param nAGQ passed into lme4::glmer(). Recommended more than 1. Defaults to 2.
#' @param root_options passed to multiroot function.
#' @param return_matrices whether to return variance matrices
#' @param verbose if TRUE then will print (lots of) output. Defaults to FALSE.

#'
#' @return a list
#' @export
policyFX <- function(
  data,
  formula,
  alphas,
  k_samps=NULL,
  ...,
  verbose = FALSE,
  root_options = NULL,
  nAGQ=2,
  return_matrices = FALSE
){

  if (is.null(k_samps)) {
    stop("Please specify a numeric for `k_samps` argument. 0 is for full sample, but <=5 is recommended")
  }
  max_k_evals <- k_samps

  dots <- list(...)
  if ("target_grid" %in% names(dots)) {
    target_grid <- dots$target_grid
    alphas_in_grid <- sort(unique(target_grid$alpha1))
    if (missing("alphas")){
      alphas <- alphas_in_grid
    }
    if (!all(alphas == alphas_in_grid)){
      warning("target_grid has different alphas than those specified; the alphas from the target_grid argument will be used")
      alphas <- alphas_in_grid
    }
  } else {
    target_grid <- makeTargetGrid(alphas = alphas)
  }

  target_grid$k <- max_k_evals

  formule <- Formula::Formula(formula)

  arguments <- argTests(
    data = data,
    formule = formule,
    alphas = alphas,
    max_k_evals = max_k_evals
  )

  data <- arguments$data ##sorted
  alphas <- arguments$alphas ##sorted
  modeling_formula <- arguments$modeling_formula
  grouping_var_name<- arguments$grouping_var_name

  ## verbose
  if (verbose){print("I. Estimate GLMER model")}
  ## verbose

  glmer_fit <- lme4::glmer(
    data = data,
    formula = modeling_formula,
    family = stats::binomial,
    nAGQ = nAGQ
  )

  model_parms <- unlist(lme4::getME(glmer_fit, c("beta","theta")))


  split_data_list <- splitData(
    data = data,
    formule = formule,
    glmer_fit = glmer_fit,
    grouping_var_name = grouping_var_name
  )

  ## verbose
  if (verbose){print("IB. get MCFP vecs")}
  ## verbose

  vec_dfm_no_alphas <- sampleVecs4Trt(
    split_data_list = split_data_list,
    max_k_evals = max_k_evals
  )


  index_references <- makeIndexReferences(
    model_parms = model_parms,vec_dfm_no_alphas = vec_dfm_no_alphas)

  typical_args <- list(
    split_data_list = split_data_list,
    model_parms = model_parms,
    verbose = verbose
  )
  CFBI_args <- append(typical_args, list(
    alphas = alphas,
    glmer_fit = glmer_fit,
    root_options = root_options,
    vec_dfm_no_alphas = vec_dfm_no_alphas
  ))
  CFBI_MCFP_list <- do.call(what=estimateNuisance,args=CFBI_args)

  CFBI_estimates <-  CFBI_MCFP_list$CFBI_estimates
  MCFP_tidy_estimates <- CFBI_MCFP_list$MCFP_tidy_estimates


  target_args <- append(typical_args, list(
    CFBI_estimates = CFBI_estimates,
    MCFP_tidy_estimates = MCFP_tidy_estimates,
    alphas = alphas,
    target_grid = target_grid
  ))
  ## to output
  targets_and_cps <- do.call(what=estimateTargets,args=target_args)


  target_estimates <- targets_and_cps$target_grid

  var_args <- target_args
  var_args$target_grid <- NULL
  var_args$target_grid_estimates <- target_estimates


  ## Diagnostics
  MCFP_tidier_estimates <- tidyOutMCFP(MCFP_tidy_estimates)

  ## to output
  parameters <-  list(
    model_parms = model_parms,
    CFBI_estimates = CFBI_estimates,
    MCFP_estimates = MCFP_tidier_estimates,
    target_estimates = target_estimates
  )
  if (verbose){print('Beginning variance calculations. This may take some time.')}

  var_args$index_references <- index_references
  var_args$return_matrices <- return_matrices



  var_ests_list <- do.call(estimateVariance, args = var_args)

  ## to ouptut
  tidy_estimates <- tidyEstimates(var_ests_list$tidy_grid)
  ## to output
  if (return_matrices== TRUE){
    variance_matrices <-  var_ests_list$var_matrices
  } else{
    variance_matrices <- list()
  }



  out_list <- list(
    estimates = tidy_estimates,
    parameters = parameters,
    variance_matrices = variance_matrices,
    prop_scores = targets_and_cps$prop_scores
  )


}



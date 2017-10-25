
# #' Calculates Variance Estimates for all Target Estimators
# #'
# #' This calculates the variance in a 'modular' fashion,
# #' and then loops over the target estimators to calculate their variances
# #'
# #' @inheritParams estimateTargets
# #' @inheritParams calcTargetBreadMeat
# #' @inheritParams policyFX
# #' @param target_grid_estimates the point estimates for all targets
# #' @param pll_var_dir see \code{\link{policyFX}}
# #' @param index_references for the indices
# #' @param deriv_simple TRUE for simple, FALSE for Richardson
# #' @param quick_MCFP_deriv FALSE would be dumb, use TRUE default
# #'
# #' @export
estimateVariance <- function(
  split_data_list,
  alphas,
  model_parms,
  CFBI_estimates,
  MCFP_tidy_estimates ,
  target_grid_estimates,
  index_references,
  deriv_simple = TRUE,
  quick_MCFP_deriv = TRUE,
  return_matrices,
  verbose
){


  if (!any(target_grid_estimates$estVar==TRUE)){
    out_list <- list(
      tidy_grid = target_grid_estimates,
      var_matrices = list()
    )
    return(out_list)
  }

  deriv_opts_list <- getDerivOptions(deriv_simple)



  ## verbose
  ckpt_time <- Sys.time()
  if (verbose){print('V. Variance estimates')}
  ## verbose

  vec_dfm_no_alphas <- MCFP_tidy_estimates[
    MCFP_tidy_estimates$alpha == MCFP_tidy_estimates$alpha[[1]],
    c("n", "s", "all_probs", "trt_vecs_mats", "mult_factor")]


  AB_nuisance_list <- calcNuisanceBreadMeat(
    split_data_list = split_data_list,
    model_parms = model_parms,
    CFBI_estimates = CFBI_estimates,
    MCFP_tidy_estimates = MCFP_tidy_estimates,
    deriv_opts = deriv_opts_list$fast,
    index_references = index_references,
    quick_MCFP_deriv=quick_MCFP_deriv,
    alphas = alphas,
    verbose = verbose
  )

  ## verbose
  # ckpt_time <- Sys.time()
  ckpt_time <- checkpointTimer(time1=ckpt_time, verbose)
  if (verbose) {print(":target ests:")}
  ## verbose

  mats_out_list <-  list()


  for (row_num in 1:nrow(target_grid_estimates)){

    target_row <- target_grid_estimates[row_num,]

    if ((!is.na(target_row$alpha2_num)) &&
        (target_row$alpha1_num == target_row$alpha2_num)) {
      target_row$estVar <- FALSE
    }

    if (target_row$estVar==FALSE) {
      target_grid_estimates$var[row_num] <- 0

      if (return_matrices){
        mats_out_list[[row_num]] <- list(Meat_mat = NA, Bread_mat = NA)
      }

    } else {

      ## verbose
      # if (row_num >1){
      if (exists("ckpt_time")){
        ckpt_time <- checkpointTimer(time1=ckpt_time, verbose)
      } else{ckpt_time <- Sys.time()}
      if (verbose) {
        print(paste('Var_row_num', row_num, 'of', nrow(target_grid_estimates)))
      }
      ## verbose


      target_args <- list(
        split_data_list = split_data_list,
        model_parms = model_parms,
        CFBI_estimates = CFBI_estimates,
        MCFP_tidy_estimates = MCFP_tidy_estimates,
        deriv_opts = deriv_opts_list$fast,
        target_row = target_row,
        index_references = index_references,
        vec_dfm_no_alphas = vec_dfm_no_alphas
        # verbose = verbose
      )


      AB_target_list <- do.call(calcTargetBreadMeat, args = target_args)



      nuisanceABMats <- constructNuisanceABMats(
        AB_nuisance_list = AB_nuisance_list,
        target_row = target_row,
        index_references = index_references
        # verbose = verbose
      )


      ## verbose
      ckpt_time <- checkpointTimer(time1=ckpt_time, verbose)
      ## verbose

      mat_calcs <-  calcModVar(
        nuisanceABMats = nuisanceABMats,
        AB_target_list = AB_target_list,
        target_row = target_row
      )
      target_row <- mat_calcs$target_row

      ## verbose
      ckpt_time <- checkpointTimer(time1=ckpt_time, verbose)
      ## verbose

      target_grid_estimates$var[row_num] <- target_row$var

      if (return_matrices){
        mats_out_list[[row_num]] <- list(
          Meat_mat = mat_calcs$Meat_mat,
          Bread_mat = mat_calcs$Bread_mat
        )
      }
    }

    target_grid_estimates$se  <- sqrt(target_grid_estimates$var)
    target_grid_estimates$LCI <- target_grid_estimates$estimate -
      stats::qnorm(0.975)*target_grid_estimates$se
    target_grid_estimates$UCI <- target_grid_estimates$estimate +
      stats::qnorm(0.975)*target_grid_estimates$se

  }

  ## verbose
  if (verbose) {print(paste('Done with all variance calculations'))}
  ## verbose

  out_list <- list(
    tidy_grid = target_grid_estimates,
    var_matrices = mats_out_list
  )
}



# #' A Function To Get Derivatives
# #'
# #' @inheritParams modularVar
# #'
getDerivOptions <- function(deriv_simple){

  if (deriv_simple==TRUE){
    deriv_method_slow <- "simple"
    deriv_method_fast <- "simple"
  } else
    if (deriv_simple==FALSE){
      deriv_method <- "Richardson"
      deriv_method_slow <- "Richardson"
      deriv_method_fast <- "Richardson"
    } else
      if (deriv_simple=="parallel"){
        deriv_method_slow <- "Richardson"
        deriv_method_fast <- "simple"
      } else{
        stop("deriv_simple must be 'parallel' for mixed or  TRUE for 'simple' or FALSE for 'Richardson', which is slower")
      }
  deriv_opts_list <- list(
    slow = list(method = deriv_method_slow),
    fast = list(method = deriv_method_fast)
  )
}

# #' @title Calculates Sigma Matrix and Variance
# #'
# #' @description The last step in the variance calculations
# #'
# #' @param nuisanceABMats blah
# #' @param AB_target_list blahblah
# #' @param target_row blahblahblah
# #'
# #' @export
calcModVar <- function(
  nuisanceABMats,
  AB_target_list,
  target_row,
  verbose
){
  Bread_mat <- nuisanceABMats$A_mat
  Meat_list <- nuisanceABMats$B_list

  Bread_mat[nrow(Bread_mat),] <- AB_target_list$A_mat
  Bread_mat_inverse <- solve(Bread_mat)

  Meat_mat <- matrix(0, ncol = ncol(Bread_mat), nrow= nrow(Bread_mat))
  for (ii in 1:length(AB_target_list$B_list)) {
    Meat_list[[ii]][length( Meat_list[[ii]])] <- AB_target_list$B_list[[ii]]
    outer_psi <- Meat_list[[ii]] %*% t(Meat_list[[ii]])
    Meat_mat <- Meat_mat + outer_psi
  }

  Sigma_mat <- Bread_mat_inverse %*% Meat_mat %*% t(Bread_mat_inverse)
  target_row$var <- Sigma_mat[nrow(Sigma_mat), ncol(Sigma_mat)]


  list(
    target_row = target_row,
    Meat_mat = Meat_mat,
    Bread_mat = Bread_mat
  )
}


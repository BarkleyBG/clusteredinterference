
# #' Point estimates for mu and CE targets
# #'
# #' @inheritParams estimateNuisance
# #' @inheritParams argTests
# #' @inheritParams makeIndexReferences
# #' @param CFBI_estimates estimated CFBI nuisance params
# #' @param split_data_list all clusters' split data; see splitData
# #' @param MCFP_tidy_estimates estimated omegas
# #' @param target_grid grid of policies to estimate mu's and CE's
# #'
# #' @export
estimateTargets <- function(
  split_data_list,
  model_parms,
  CFBI_estimates,
  MCFP_tidy_estimates ,
  alphas,
  target_grid,
  verbose
){

  ## verbose
  if(verbose){print('IV. Target estimates')}
  ## verbose

  target_grid$estimate <- NA

  group_data_with_ps <-  ###IPW
    lapply(split_data_list, function(data_list){

      PS_onegroup <- targetEstFunLHS(
        return_prop_score_only = TRUE,
        model_parms = model_parms,
        outcome_vec =  data_list$outcome_vec,
        treatment_vec =  data_list$treatment_vec,
        model_mat =  data_list$model_matrix,
        MCFP = NULL
      )
      data.frame(
        CPS = PS_onegroup,
        cluster_name = data_list$cluster_name,
        stringsAsFactors = FALSE
      )
    })


  CPS_dfm <- do.call(rbind, group_data_with_ps)

  for (row_num in 1:nrow(target_grid)) {

    alpha1 <- target_grid$alpha1[row_num]
    alpha2 <- target_grid$alpha2[row_num]
    trt <- target_grid$trt[row_num]

    if ((is.na(alpha2)) || (alpha1 != alpha2)) {
      target_allgroups_list <-
        lapply(
          split_data_list,
          FUN = estimateGroupTarget,
          CPS_dfm = CPS_dfm,
          MCFP_tidy_estimates = MCFP_tidy_estimates,
          alpha1 = alpha1,
          alpha2 = alpha2,
          trt = trt
        )
      target_allgroups <- c(unlist(target_allgroups_list))

      target_estimate <- mean(target_allgroups, na.rm = TRUE)
    } else{
      target_estimate <- 0
    }

    target_grid$estimate[row_num] <- target_estimate
  }



  out_list <- list(target_grid = target_grid,
                   prop_scores = CPS_dfm)
}



# #' puts MCFP, Prop Score, and outcome vec together for group's contribution to target
# #'
# #' @inheritParams estimateTargets
# #' @inheritParams argTests
# #' @inheritParams makeTargetEstFun
# #' @param data_list one group's worth of split_data; see \code{\link{splitData}}
# #' @param alpha1 The first policy
# #' @param alpha2 the second policy (for CE) or NULL (for mu's)
# #' @param CPS_dfm the dataframe of cluster ID's and prop scores
# #'
# #' @export
estimateGroupTarget <- function(
  data_list,
  CPS_dfm,
  alpha1,
  alpha2,
  trt,
  MCFP_tidy_estimates
){

  CPS <- CPS_dfm$CPS[CPS_dfm$cluster_name==data_list$cluster_name]

  ### defense
  if (length(CPS)!=1){ stop("wrong CPS size")}
  ### defense

  MCFP_alpha1 <- getMCFPFromTidy(
    treatment_vec = data_list$treatment_vec,
    alpha = alpha1,
    MCFP_tidy_estimates = MCFP_tidy_estimates
  )

  if (!is.na(alpha2)){
    MCFP_alpha2 <- getMCFPFromTidy(
      treatment_vec = data_list$treatment_vec,
      alpha = alpha2,
      MCFP_tidy_estimates = MCFP_tidy_estimates
    )
    MCFP <- MCFP_alpha1-MCFP_alpha2
  } else {
    MCFP <- MCFP_alpha1
  }

  group_target_estimate <-
    targetEstFunLHS(
      MCFP = MCFP,
      outcome_vec = data_list$outcome_vec,
      treatment_vec = data_list$treatment_vec,
      trt = trt,
      CPS = CPS
    )
}

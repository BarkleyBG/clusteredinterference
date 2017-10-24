

# #' function for obtaining 'quasi-outcomes' from vector
# #'
# #' @param y the vector of outcomes
# #' @param a the vector of treatment
# #' @param t 0 for untreated, 1 for treated
# #'
# ##' @export
calcAvgOutcome <- function(y,a,t){
  if (t %in% 0:1){
    sub_y <- y[a==t]
    if (length(sub_y)==0) {
      # out <- 0 ##prior to August 6
      # out <- NA_real_ ## Starting August 6 in QE_NAs branch
      out <- 0 ##starting August 12th
    } else {
      out <- mean(sub_y)
    }
  } else {
    stop("t must equal 0 or 1")
  }
  return(out)
}

# # #' function for calculating outcome and quasi-outcomes from vector
# #'
# #' @param outcome_vec the vector of outcomes
# #' @param treatment_vec the vector of treatment
# #' @param trt 0 for untreated, 1 for treated, NA for overall
# #'
# #' @export
calcYBar <- function(outcome_vec, treatment_vec, trt){
 if (is.na(trt)){ avg_outcome <- mean(outcome_vec) } else{
    avg_outcome   <- calcAvgOutcome(y=outcome_vec, a = treatment_vec, t=trt)
 }
  avg_outcome
}










###########################
###### estimateMCFP ########
###########################

# #' Estimate the omegas, for point estimates
# #'
# #' @inheritParams getEstdCFBI_eeroot
# #' @param CFBI_estimates estimated fixed effect intercept in counterfactual worlds
# #' @param vec_dfm_no_alphas The vectors to evaluate for estimating omegas
# #'
# #' @export
estimateMCFP <- function(
  split_data_list,
  model_parms,
  CFBI_estimates,
  alphas,
  vec_dfm_no_alphas,
  verbose
){

  ## verbose
  ckpt_time <- Sys.time()
  if(verbose){print('III. Estimate MCFP')}
  ## verbose



  MCFP_onealpha_dfms <- lapply(1:length(alphas), function(alpha_num) {

    ## verbose
    if (alpha_num > 1) {  ckpt_time <- checkpointTimer(time1 = ckpt_time,verbose)}
    if(verbose){print(paste('III. alpha num', alpha_num, 'of', length(alphas)))}
    ## verbose

    MCFP_onealpha_dfm <- getMCFP1AlphaNum(
      alpha_num = alpha_num,
      alpha = alphas[alpha_num],
      CFBI_estimates = CFBI_estimates,
      split_data_list = split_data_list,
      model_parms = model_parms,
      vec_dfm_no_alphas = vec_dfm_no_alphas
    )
  })

  MCFP_estimates_dfm <- do.call(rbind,MCFP_onealpha_dfms)

  ### defense
  if (nrow(MCFP_estimates_dfm)!=(length(alphas)*nrow(vec_dfm_no_alphas))) {
    stop(
    'MCFP_estimates_dfm must have one row for each row in vec_dfm_no_alphas'
    )}
  ### defense

  ## verbose
  ckpt_time <- checkpointTimer(time1=ckpt_time, verbose)
  ## verbose

  row.names(MCFP_estimates_dfm) <- NULL
  MCFP_estimates_dfm
}


# #' innards of estimateMCFP
# #'
# #' @inheritParams estimateTargets
# #' @inheritParams makeIndexReferences
# #' @inheritParams argTests
# #' @inheritParams calcCFBIMCFP_onealpha
# #' @param alpha_num the number of the alpha parm
# #'
# #' @export
getMCFP1AlphaNum <- function(
  alpha_num,
  alpha,
  CFBI_estimates,
  split_data_list,
  model_parms,
  vec_dfm_no_alphas
){

  if(length(CFBI_estimates)==1){
    CFBI <- CFBI_estimates
  } else {
    CFBI <- CFBI_estimates[alpha_num]
  }
  MCFP_onealpha_allgroups_list <-
    lapply(split_data_list, function(data_list){

      force(data_list)
      MCFP_onealpha_onegroup_vec <-
        MCFP_EEqn_LHS( ##for one group
          model_parms = model_parms,
          CFBI = CFBI,
          vec_dfm_no_alphas = vec_dfm_no_alphas,
          model_DV_vec = data_list$model_DV_vec,
          model_mat = data_list$model_matrix
        )## returns a vector of each group's contributions to all MCFP's.

      ### defense
      if (nrow(vec_dfm_no_alphas)!= length(MCFP_onealpha_onegroup_vec)) {stop(
        'MCFP_onealpha_onegroup_vec must have one prob for each row in vec_dfm_no_alphas'
      )}
      ### defense


      MCFP_onealpha_onegroup_vec
    })

  MCFP_onealpha_allgroups <- do.call(cbind, MCFP_onealpha_allgroups_list)

  ### defense
  ##number of columns equals number of groups - very wide
  if (ncol(MCFP_onealpha_allgroups) !=
      length(split_data_list)) { stop(
         'MCFP_onealpha_allgroups must have one row for each group'
      )}
  ##number of rows equals number of probabilities
  if (nrow(MCFP_onealpha_allgroups)!=nrow(vec_dfm_no_alphas)) { stop(
    'MCFP_onealpha_allgroups must have one column for each alpha'
  )}
  ### defense

  MCFP_onealpha_avg <- apply(
    MCFP_onealpha_allgroups, 1, function(x){mean(stats::na.omit(x))})
  MCFP_onealpha_avg <- as.vector(MCFP_onealpha_avg)

  ### defense
  ## putting the MCFP_alpha_ests vertical for each alpha;
  ## alphas are columns
  if (length(MCFP_onealpha_avg)!=nrow(vec_dfm_no_alphas)){
    stop("MCFP_onealpha_avg should be a vector of nrow(vec_dfm_no_alphas) probs")
  }
  ### defense

  mcfp_onealpha_dfm <- vec_dfm_no_alphas
  mcfp_onealpha_dfm$alpha <- alpha
  mcfp_onealpha_dfm$alpha_num <- alpha_num
  mcfp_onealpha_dfm$CFBI <- CFBI
  mcfp_onealpha_dfm$MCFP <- MCFP_onealpha_avg

  mcfp_onealpha_dfm
}

# #' Makes output from \code{\link{estimateMCFP}} even more tidy for reader
# #'
# #' @param MCFP_tidy_estimates some un-tidied omegas
# #'
# #' @export
tidyOutMCFP <- function(MCFP_tidy_estimates){
  MCFP_list <- by(MCFP_tidy_estimates, MCFP_tidy_estimates$n, function(MCFP_n){
    if ( any(MCFP_n$all_probs == "removed_one")) {
      if (!( all(MCFP_n$all_probs=="removed_one"))) {
        stop("all should list 'removed one' if any does")
      }
      if (!( all(MCFP_n$n==MCFP_n$n[1]))) {
        stop("all sizes should be the same")
      }
      n_size <- MCFP_n$n[1]
      all_trts <- 0:n_size
      trts_estd <- MCFP_n$s[1:n_size] ## only for the first alpha
      trt_skipped <- all_trts[!(all_trts %in% trts_estd)]
      MCFP_n_out_list <- by(MCFP_n, MCFP_n$alpha, function(MCFP_n_a){
        MCFP_skipped <- 1-sum(MCFP_n_a$MCFP)
        MCFP_skipped_row <- data.frame(
          stringsAsFactors = FALSE,
          n = MCFP_n_a$n[1],
          s = trt_skipped,
          all_probs = "this_one_removed",
          trt_vecs_mats = NA,
          mult_factor = NA_real_,
          alpha = MCFP_n_a$alpha[1],
          alpha_num = MCFP_n_a$alpha_num[1],
          CFBI = MCFP_n_a$CFBI[1],
          MCFP = MCFP_skipped
        )

        MCFP_out <- rbind(MCFP_n_a, MCFP_skipped_row)
        return(MCFP_out)
      })
      MCFP_n_out <- do.call(rbind, MCFP_n_out_list)
      return(MCFP_n_out)
    } else {
      return(MCFP_n)
    }
  })

  MCFP_out <- do.call(rbind, MCFP_list)
  return(MCFP_out)
}


# #' Make Useful References for EEqn Var
# #'
# #' This should help to DRY with code
# #'
# #' @param model_parms the model parameter estimates
# #' @param vec_dfm_no_alphas provides the number of MCFP/omega estimates
# #'
# #' @export
makeIndexReferences <- function(
  model_parms,
  vec_dfm_no_alphas
){
  lmp <- length(model_parms)
  nrvna <- nrow(vec_dfm_no_alphas)

  mu_refs <- list(
    model = 1:lmp,
    CFBI = lmp+1,
    MCFP = lmp+1+1:nrvna,
    mu = lmp+1+nrvna+1,
    model_CFBI = 1:( lmp+1 ),
    model_CFBI_MCFP = 1:(lmp+1+nrvna)
  )

  ### defense
  if (mu_refs$mu!=(1+mu_refs$model_CFBI_MCFP[length(mu_refs$model_CFBI_MCFP)])){
    stop("refs are wrong")}
  ### defense

  CE_refs <- list(
    model = 1:lmp,
    CFBI1 = lmp+1,
    CFBI2 = lmp+2,
    MCFP1 = lmp+2+1:nrvna,
    MCFP2 = lmp+2+nrvna+1:nrvna,
    CE = lmp+2+2*nrvna+1,
    model_CFBI1 = 1:(lmp+1),
    model_CFBI2 = c(1:lmp,lmp+2),
    model_CFBI1_MCFP1 = c( 1:(lmp+1), lmp+2+1:nrvna),
    model_CFBI2_MCFP2 = c( 1:lmp,lmp+2, lmp+2+nrvna+1:nrvna)
  )


  ### defense
  if (length(CE_refs$CE)!=1) {stop("refs are wrong")}
  if (CE_refs$CE!=(1+CE_refs$model_CFBI2_MCFP2[length(CE_refs$model_CFBI2_MCFP2)])){
    stop("refs are wrong")}
  ### defense

  index_references <- list(
    mu_refs = mu_refs,
    CE_refs = CE_refs
  )
}


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

      PS_onegroup <- targetEEqnLHS(
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
    targetEEqnLHS(
      MCFP = MCFP,
      outcome_vec = data_list$outcome_vec,
      treatment_vec = data_list$treatment_vec,
      trt = trt,
      CPS = CPS
    )
}




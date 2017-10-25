
# #' Estimate Causal Nuisance Parameters
# #'
# #' This function can estimate the CFBI and MCFP parameters in several ways
# #'
# #' @inheritParams policyFX
# #' @inheritParams getCFBIMCFP
# #' @inheritParams estimateTargets
# #' @param vec_dfm_no_alphas bookkeeping for MCFP
# #'
# #' @export
estimateNuisance <- function(
  split_data_list,
  model_parms,
  alphas,
  glmer_fit,
  root_options,
  vec_dfm_no_alphas,
  verbose
){

  basic_arg_list <- list(
    split_data_list = split_data_list,
    model_parms = model_parms,
    alphas = alphas
  )
  CFBI_arg_list <- append(basic_arg_list, list(
    glmer_fit = glmer_fit,
    CFBI_numstart = 50,
    root_options = root_options,
    verbose = verbose
  ))
  CFBI_estimates <- do.call(getEstdCFBI_eeroot, CFBI_arg_list)

  MCFP_arg_list <- append(basic_arg_list, list(
    CFBI_estimates = CFBI_estimates,
    vec_dfm_no_alphas = vec_dfm_no_alphas,
    verbose = verbose
  ))
  MCFP_tidy_estimates <- do.call(estimateMCFP, MCFP_arg_list)


  out_list <- list(
    CFBI_estimates = CFBI_estimates,
    MCFP_tidy_estimates = MCFP_tidy_estimates
  )
}







# #' Estimate counterfactual fixed intercepts
# #'
# #' @inheritParams argTests
# #' @inheritParams estimateTargets
# #' @param model_parms estimated GLMER model parameters
# #' @param glmer_fit the fitted GLMER model object, for determining starting CFBI values
# #' @param CFBI_numstart number of CFBI's to test initially
# #' @param root_options additional options to pass
# #'
# #' @export
getEstdCFBI_eeroot <- function(
  split_data_list,
  model_parms,
  alphas,
  glmer_fit,
  CFBI_numstart,
  root_options = NULL,
  verbose
){
  if(verbose){print(paste("II. Estimate CFBI"))}
  ckpt_time <- Sys.time()

  beta_here <- model_parms[-length(model_parms)]
  theta_here <- model_parms[length(model_parms)]
  names(beta_here) <- names(theta_here) <- NULL

  new_param_list <- list(
    beta = beta_here,
    theta = theta_here
  )

  starting_dfm <- getStartingCFBI(
    new_param_list = new_param_list,
    glmer_fit = glmer_fit,
    alphas = alphas,
    CFBI_numstart = CFBI_numstart
  )


  CFBI_estimate <- rep(NA_real_, length(alphas))

  for (alpha_num in 1:length(alphas)) {

    alpha <- alphas[alpha_num]

    ## verbose
    if (alpha_num>1){ckpt_time <- checkpointTimer(time1=ckpt_time, verbose)}
    if (verbose) {print(paste("II. solving for alpha_num", alpha_num,
                              "of", length(alphas), "; alpha=", alpha))}
    ## verbose

    est_CFBI1Alpha_args <- list(
      split_data_list = split_data_list,
      model_parms = model_parms,
      alpha = alpha,
      starting_dfm = starting_dfm,
      root_options = root_options
    )


    CFBI_alpha_solved <- do.call(getCFBI1Alpha, est_CFBI1Alpha_args)

    CFBI_estimate[alpha_num] <- CFBI_alpha_solved$root
  }

  if (alpha_num == length(alphas)){
    ckpt_time <- checkpointTimer(time1=ckpt_time,verbose)
  }

  CFBI_estimate
}


# #' Determine starting points for CFBI estimation
# #'
# #' @inheritParams getEstdCFBI_eeroot
# #' @param new_param_list a list of estimated fixef's and raneff's from GLMER model
# #'
# #' @export
getStartingCFBI <- function(
  new_param_list,
  glmer_fit,
  alphas,
  CFBI_numstart
){
  starting_CFBIs <- seq(from=-3, to = 3,length.out=CFBI_numstart)
  starting_alphas <- guessStartingAlphas(
    starting_CFBIs = starting_CFBIs,
    new_param_list = new_param_list,
    glmer_fit = glmer_fit
  )
  starting_dfm <- data.frame(
    starting_alphas = unlist(starting_alphas),
    starting_CFBIs = starting_CFBIs
  )
  max_start_alpha <- starting_dfm$starting_alphas[nrow(starting_dfm)]
  min_start_alpha <- starting_dfm$starting_alphas[1]

  if ( sum(alphas>max_start_alpha)>=3) {
    print("expanding grid of starting values for CFBI")
    starting_CFBIs <- seq(from=2.5, to = 10, by = .75)

    starting_alphas <- guessStartingAlphas(
      starting_CFBIs = starting_CFBIs,
      new_param_list = new_param_list, glmer_fit = glmer_fit
    )
    starting_dfm2 <- data.frame(
      starting_alphas = unlist(starting_alphas),
      starting_CFBIs = starting_CFBIs
    )
    starting_dfm <- rbind(starting_dfm,starting_dfm2)
  }
  if ( sum(alphas<min_start_alpha)>=3) {
    print("expanding grid of starting values for CFBI")
    starting_CFBIs <- seq(from=-10, to = -2.5, by = .75)

    starting_alphas <- guessStartingAlphas(
      starting_CFBIs = starting_CFBIs,
      new_param_list = new_param_list, glmer_fit = glmer_fit
    )
    starting_dfm2 <- data.frame(
      starting_alphas = unlist(starting_alphas),
      starting_CFBIs = starting_CFBIs
    )
    starting_dfm <- rbind(starting_dfm,starting_dfm2)
  }
  max_start_alpha <- starting_dfm$starting_alphas[nrow(starting_dfm)]
  min_start_alpha <- starting_dfm$starting_alphas[1]
  if ( sum(alphas>max_start_alpha)>3) {
    warning("solver may be slow for CFBI's due to large values of alpha")
  }
  if ( sum(alphas<min_start_alpha)>3) {
    warning("solver may be slow for CFBI's due to small values of alpha")
  }
  return(starting_dfm)
}

# #' quickens CFBI, i think
# #'
# #' @inheritParams getStartingCFBI
# #' @param starting_CFBIs a vector of starting CFBIs
# #'
# #' @export
guessStartingAlphas <- function(
  starting_CFBIs,
  new_param_list,
  glmer_fit
){
  starting_alphas <- lapply(starting_CFBIs, function(start_val){
    new_param_list$beta[1] <- start_val
    indiv_probs <- suppressMessages(
      stats::predict(
        glmer_fit,
        newparams = new_param_list,
        type = "response"
      )
    )
    if ("try-error" %in% class(indiv_probs)){
      indiv_probs <- 99 ##no alpha will select this
    }

    mean(indiv_probs)
  })
}



# #' inner workings of getEstdCFBI_eeroot
# #'
# #' @inheritParams getEstdCFBI_eeroot
# #' @inheritParams argTests
# #' @inheritParams estimateTargets
# #' @param starting_dfm initial start
# #' @param alpha just the one alpha to use
# #'
# #' @export
getCFBI1Alpha <- function(
  alpha,
  split_data_list,
  model_parms,
  starting_dfm,
  root_options
){


  ordered_rows <- order(abs(starting_dfm$starting_alphas - alpha))
  starting_estimate <- mean(starting_dfm[ordered_rows,]$starting_CFBIs[1:2])

  CFBI_root_args <- list(
    geex_list = list(
      eeFUN = CFBI_only_EstFun,
      splitdt = split_data_list
    ),
    start = starting_estimate,
    root_options = root_options,
    alpha = alpha,
    model_parms = model_parms
  )
  do.call(calculateRoots,CFBI_root_args)
}



# #' Calculates Roots
# #'
# #' @description This estimates parameters. Thanks to package geex.
# #'
# #' @inheritParams calculateBreadMeat
# #' @param start the starting values for rootSolve::multiroot
# #' @param root_options passed into rootSolve::multiroot
# #' @param ... passed into rootsolver
# #'
# #' @export
calculateRoots <- function(geex_list, start, root_options, ...) {
  psi_i <- lapply(geex_list$splitdt, function(data_list_i) {
    force(data_list_i)
    geex_list$eeFUN(data_list = data_list_i, ...)
  })
  psi <- function(theta) {
    psii <- lapply(psi_i, function(f) {
      do.call(f, args = append(list(theta = theta), geex_list$ee_args))
    })
    apply(checkArray(simplify2array(psii)), 1, sum)
  }
  root_args <- append(root_options, list(f = psi, start = start))
  do.call(rootSolve::multiroot, args = root_args)
}



# # #' CFBI and MCFP together,  not by alphas
# # #'
# # #' @param split_data_list split data
# # #' @inheritParams argTests
# # #' @param model_parms estimated GLMER model parameters
# # #' @param glmer_fit the fitted GLMER model object, for determining starting CFBI values
# # #' @param CFBI_numstart number of CFBI's to test initially
# # #' @param root_options additional options to pass
# # #' @param vec_dfm_no_alphas The vectors to evaluate for estimating omegas
# # #'
# # #' @export
# getCFBIMCFP <- function(
#     split_data_list,
#     model_parms,
#     alphas,
#     glmer_fit,
#     CFBI_numstart,
#     root_options,
#     verbose,
#     vec_dfm_no_alphas
# ){
#
#   ## verbose
#   if(verbose) {print(paste("II. Estimate CFBI"))}
#   ckpt_time <- Sys.time()
#   ## verbose
#
#
#   beta_here <- model_parms[-length(model_parms)]
#   theta_here <- model_parms[length(model_parms)]
#   names(beta_here) <- names(theta_here) <- NULL
#
#   new_param_list <- list(
#     beta = beta_here,
#     theta = theta_here
#   )
#
#   starting_dfm <- getStartingCFBI(
#     new_param_list = new_param_list,
#     glmer_fit = glmer_fit,
#     alphas = alphas,
#     CFBI_numstart = CFBI_numstart
#   )
#
#
#   CFBI_estimates <- rep(NA_real_, length(alphas))
#
#   for (alpha_num in 1:length(alphas)) {
#
#     alpha <- alphas[alpha_num]
#
#     ## verbose
#     if (verbose) {
#       if (alpha_num>1){ ckpt_time <- checkpointTimer(time1=ckpt_time, verbose)}
#       print(paste( "II. solving for alpha_num",
#                    alpha_num, "of", length(alphas), "; alpha=", alpha))
#     }
#     ## verbose
#
#     CFBI_alpha_solved <- getCFBI1Alpha(
#       split_data_list = split_data_list,
#       model_parms = model_parms,
#       alpha = alpha,
#       starting_dfm = starting_dfm,
#       root_options = root_options
#     )
#
#     ## verbose
#     if (verbose) {
#       print(paste("II. solving for alpha_num", alpha_num,
#                   "of", length(alphas), "; alpha=", alpha))
#     }
#     ## verbose
#
#     CFBI_estimates[alpha_num] <- CFBI_alpha_solved$root
#   }
#
#   ## verbose
#   ckpt_time <- checkpointTimer(time1=ckpt_time,verbose)
#   ## verbose
#
#   print('III. Estimate MCFP')
#
#
#
#   MCFP_onealpha_dfms <- lapply(1:length(alphas), function(alpha_num) {
#     alpha <- alphas[alpha_num]
#
#     ## verbose
#     if (alpha_num > 1) {ckpt_time <- checkpointTimer(time1 = ckpt_time, verbose)}
#     if (verbose) {print(paste('III. alpha num', alpha_num, 'of', length(alphas)))}
#     ## verbose
#
#
#
#     MCFP_onealpha_dfm <- getMCFP1AlphaNum(
#       alpha_num = alpha_num,
#       alpha  = alpha ,
#       CFBI_estimates = CFBI_estimates,
#       split_data_list = split_data_list,
#       model_parms = model_parms,
#       vec_dfm_no_alphas = vec_dfm_no_alphas
#     )
#   })
#   MCFP_estimates_dfm <- do.call(rbind,MCFP_onealpha_dfms)
#
#   if (nrow(MCFP_estimates_dfm)!=(length(alphas)*nrow(vec_dfm_no_alphas))) {
#
#     stop(
#       'MCFP_estimates_dfm must have one row for each row in vec_dfm_no_alphas'
#     )}
#
#   ## verbose
#   ckpt_time <- checkpointTimer(time1=ckpt_time, verbose)
#   if (verbose) {print("done with CFBI and MCFP")}
#   ## verbose
#
#   row.names(MCFP_estimates_dfm) <- NULL
#
#
#   out_list <- list(
#     CFBI_estimates = CFBI_estimates,
#     MCFP_tidy_estimates = MCFP_estimates_dfm
#   )
# }






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
        MCFP_EstFun_LHS( ##for one group
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

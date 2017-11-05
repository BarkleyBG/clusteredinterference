



# #' Calculates the bread and meat matrices for the nuisance parameters
# #'
# #' @inheritParams modularVar
# #' @inheritParams calcTargetBreadMeat
# #'
# #' @export
calcNuisanceBreadMeat <- function(
  split_data_list,
  model_parms,
  CFBI_estimates,
  MCFP_tidy_estimates,
  deriv_opts,
  quick_MCFP_deriv,
  index_references,
  alphas,
  verbose
){

  ## verbose
  if (verbose) {print(":model stuff:")}
  ckpt_time <- Sys.time()
  ## verbose

  AB_model_list <- calculateBreadMeat(
    geex_list = list(
      eeFUN = modular_model_estfun,
      splitdt = split_data_list,
      ee_args = NULL),
    theta = model_parms,
    numDeriv_options = deriv_opts,
    silent = TRUE
  )



  ## verbose
  ckpt_time <- checkpointTimer(time1=ckpt_time, verbose)
  if (verbose) {print(":CFBI EstFun:")}
  ## verbose

  geex_list_CFBI <- list(
    eeFUN = modular_CFBI_estfun,
    splitdt = split_data_list,
    ee_args = NULL
  )



  AB_CFBI_list <-
    lapply(1:length(alphas), function(this_alpha_num) {

      ## verbose
      if (verbose){
        ckpt_time <- checkpointTimer(time1=ckpt_time, verbose)
        print(paste(':CFBI alpha_num', this_alpha_num, 'of', length(alphas)))
      }
      ## verbose

      one_CFBI_estimate <- CFBI_estimates[this_alpha_num]
      one_alpha <- alphas[this_alpha_num]
      model_and_one_CFBI_estimate <- c(model_parms, one_CFBI_estimate)

      AB_CFBI_stuff <- calculateBreadMeat(
        geex_list = geex_list_CFBI,
        theta = model_and_one_CFBI_estimate,
        numDeriv_options = deriv_opts,
        silent = TRUE,
        alpha = one_alpha
      )
    })



  if (verbose) {print(":MCFP EstFun:")}

  vec_dfm_no_alphas <-
    MCFP_tidy_estimates[MCFP_tidy_estimates$alpha == MCFP_tidy_estimates$alpha[[1]],
                        c("n", "s", "all_probs", "trt_vecs_mats", "mult_factor")]

  # geex_list_MCFP <- list(eeFUN = modular_MCFP_estfun,
  #                        splitdt = split_data_list,
  #                        ee_args = NULL)




  AB_MCFP_list <-
    lapply(1:length(alphas), function(this_alpha_num) {
      one_alpha <- alphas[this_alpha_num]
      one_CFBI_estimate <- CFBI_estimates[this_alpha_num]


      ## verbose
      if (verbose) {
        ckpt_time <- checkpointTimer(time1=ckpt_time, verbose)
        print(paste(':MCFP alpha_num', this_alpha_num, 'of', length(alphas)))
      }
      ## verbose


      one_dfm_MCFP_estimates <- MCFP_tidy_estimates[
        MCFP_tidy_estimates$alpha_num == this_alpha_num,]
      ## defense
      if (any(abs(one_dfm_MCFP_estimates$alpha - one_alpha) > 1e-3)) {
        stop("grabbing the wrong alphas in MCFP estimates")
      }
      if (any(abs(one_dfm_MCFP_estimates$CFBI - one_CFBI_estimate) > 1e-3)) {
        stop("grabbing the wrong alphas in MCFP estimates")
      }
      ## defense
      one_vec_MCFP_estimates <- one_dfm_MCFP_estimates$MCFP

      model_CFBI_MCFP_estimates <-
        c(model_parms, one_CFBI_estimate, one_vec_MCFP_estimates)
      model_CFBI_estimates_NO_MCFP <- c(model_parms, one_CFBI_estimate)


      # if (quick_MCFP_deriv) {

        A_SUBmat_MCFP_onealpha <- calculateBreadMeat(
          calc_which="Bread",
          geex_list = list(
            eeFUN = modular_MCFP_estfun_quick,
            splitdt = split_data_list,
            ee_args = NULL
          ),
          theta = model_CFBI_estimates_NO_MCFP,
          numDeriv_options = deriv_opts,
          silent = TRUE,
          alpha = one_alpha,
          num_model_parms = length(model_parms),
          MCFP_onealpha = one_vec_MCFP_estimates,
          vec_dfm_no_alphas = vec_dfm_no_alphas
        )
        A_SUBmat <- A_SUBmat_MCFP_onealpha$A_mat

        B_list_MCFP_onealpha <- calculateBreadMeat(
          calc_which="Meat",
          geex_list = list(
            eeFUN = modular_MCFP_estfun_quick,
            splitdt = split_data_list,
            ee_args = NULL
          ),
          theta = model_CFBI_estimates_NO_MCFP,
          numDeriv_options = deriv_opts,
          silent = TRUE,
          alpha = one_alpha,
          num_model_parms = length(model_parms),
          MCFP_onealpha = one_vec_MCFP_estimates,
          vec_dfm_no_alphas = vec_dfm_no_alphas
        )


        out <-  list(A_SUBmat = A_SUBmat,
                     B_list0 = B_list_MCFP_onealpha$B_list)


    #   } else {
    #     AB_MCFP_onealpha <- calculateBreadMeat(
    #       geex_list = list(
    #         eeFUN = modular_MCFP_estfun,
    #         splitdt = split_data_list,
    #         ee_args = NULL
    #       ),
    #       theta = model_CFBI_MCFP_estimates,
    #       numDeriv_options = deriv_opts,
    #       silent = TRUE,
    #       alpha = one_alpha,
    #       num_model_parms = length(model_parms),
    #       vec_dfm_no_alphas = vec_dfm_no_alphas
    #     )
    #     out <- AB_MCFP_onealpha
    #   }
      out
    })

  # if (quick_MCFP_deriv) {
    group_sizes <- sapply(split_data_list,
                          function(data_list) {
                            Nii <- nrow(data_list$model_matrix)
                          })
    group_size_tab <- table(sort(group_sizes))
    param_tab <- table(vec_dfm_no_alphas$n)

    ## defense
    stopifnot(all(unique(vec_dfm_no_alphas$n) == as.numeric(names(group_size_tab))))
    stopifnot(length(param_tab) == length(group_size_tab))
    ## defense
    diag_list <- list()
    for (size_num in 1:length(param_tab)) {
      diag_list[[size_num]] <- rep(group_size_tab[size_num],
                                   param_tab[size_num])
    }
    mat_diag <- c(unlist(diag_list))
    A_extramat <- diag(mat_diag)

    AB_MCFP_list$A_MCFP <- A_extramat
  # } ##end of MCFP_quick


  AB_nuisance_list <- list(model = AB_model_list,
                           CFBI = AB_CFBI_list,
                           MCFP = AB_MCFP_list,
                           vec_dfm_no_alphas = vec_dfm_no_alphas
  )
  # return(AB_nuisance_list)
}


# #' Puts together everything but the target param's Bread contribution.
# #'
# #' @param AB_nuisance_list ALl the nuisance bread
# #' @param target_row the estimate and info on alpha1, alpha2
# #' @param index_references references of indices.
# #'
# #' @export
constructNuisanceABMats <- function(
  AB_nuisance_list,
  target_row,
  index_references
){

  alpha1_num <- target_row$alpha1_num
  alpha2_num <- target_row$alpha2_num
  trt <- target_row$trt

  if (is.na(alpha2_num)){

    A_submat_model <- AB_nuisance_list$model$A_mat
    A_submat_CFBI <- AB_nuisance_list$CFBI[[alpha1_num]]$A_mat
    A_submat_MCFP_sub <- AB_nuisance_list$MCFP[[alpha1_num]]$A_SUBmat
    A_submat_MCFP_extra <- AB_nuisance_list$MCFP$A_MCFP

    model_indices <- index_references$mu_refs$model
    CFBI_index <- index_references$mu_refs$CFBI
    model_CFBI_indices <- index_references$mu_refs$model_CFBI
    MCFP_indices <- index_references$mu_refs$MCFP
    model_CFBI_MCFP_indices <- index_references$mu_refs$model_CFBI_MCFP

    mat_size <- index_references$mu_refs$mu
    A_mat <- matrix(NA, ncol = mat_size, nrow = mat_size)


    ### defense
    if (nrow(A_submat_model)!=ncol(A_submat_model)) {stop("error in matrix construction")}
    if (MCFP_indices[length(MCFP_indices)]+1!= mat_size) {stop("error in matrix construction")}
    ### defense

    A_mat[model_indices,model_indices] <- A_submat_model
    A_mat[model_indices,-model_indices] <- 0

    A_mat[CFBI_index,model_CFBI_indices] <- A_submat_CFBI
    A_mat[CFBI_index,-model_CFBI_indices] <- 0

    A_mat[MCFP_indices,model_CFBI_indices] <- A_submat_MCFP_sub
    A_mat[MCFP_indices,MCFP_indices] <- A_submat_MCFP_extra
    A_mat[MCFP_indices,-model_CFBI_MCFP_indices] <- 0

    A_mat[-nrow(A_mat), ncol(A_mat)] <- 0

    num_groups <- length(AB_nuisance_list$model$B_list)
    B_vecs_list <- list()

    B_vec <- rep(NA,length=mat_size)
    for (ii in 1:num_groups){
      B_vec[model_indices] <- AB_nuisance_list$model$B_list[[ii]]
      B_vec[CFBI_index] <- AB_nuisance_list$CFBI[[alpha1_num]]$B_list[[ii]]
      B_vec[MCFP_indices] <- AB_nuisance_list$MCFP[[alpha1_num]]$B_list[[ii]]

      ### defense
      if (any(!is.na(B_vec[length(B_vec)]))) {stop("all last entry needs to be NA")}
      if (any(is.na(B_vec[-length(B_vec)]))) {stop("only last entry can be NA")}
      ### defense

      B_vecs_list[[ii]] <- B_vec
    }




  } else {

    A_submat_model <- AB_nuisance_list$model$A_mat
    A_submat_CFBI1 <- AB_nuisance_list$CFBI[[alpha1_num]]$A_mat
    A_submat_CFBI2 <- AB_nuisance_list$CFBI[[alpha2_num]]$A_mat
    A_submat_MCFP1_sub <- AB_nuisance_list$MCFP[[alpha1_num]]$A_SUBmat
    A_submat_MCFP2_sub <- AB_nuisance_list$MCFP[[alpha2_num]]$A_SUBmat
    A_submat_MCFP_extra <- AB_nuisance_list$MCFP$A_MCFP



    model_indices <- index_references$CE_refs$model #1:nrow(A_submat_model)
    CFBI1_index <- index_references$CE_refs$CFBI1 #CFBI_index <- 1+nrow(A_submat_model)
    model_CFBI1_indices <- index_references$CE_refs$model_CFBI1
    MCFP1_indices <- index_references$CE_refs$MCFP1 #CFBI_index + 1:num_MCFP
    model_CFBI_MCFP1_indices <- index_references$CE_refs$model_CFBI1_MCFP1

    CFBI2_index <- index_references$CE_refs$CFBI2 #CFBI_index <- 1+nrow(A_submat_model)
    model_CFBI2_indices <- index_references$CE_refs$model_CFBI2
    MCFP2_indices <- index_references$CE_refs$MCFP2 #CFBI_index + 1:num_MCFP
    model_CFBI_MCFP2_indices <- index_references$CE_refs$model_CFBI2_MCFP2


    mat_size <- index_references$CE_refs$CE
    A_mat <- matrix(NA, ncol = mat_size, nrow = mat_size)

    ### defense
    if( nrow(A_submat_model)!=ncol(A_submat_model)) {stop("error in matrix construction")}
    if ((model_CFBI_MCFP2_indices[length(model_CFBI_MCFP2_indices)]+1)!= mat_size) {stop("error in matrix construction")}
    ### defense


    A_mat[model_indices,model_indices] <- A_submat_model
    A_mat[model_indices,-c(model_indices)] <- 0

    A_mat[CFBI1_index,model_CFBI1_indices] <- A_submat_CFBI1
    A_mat[CFBI1_index,-model_CFBI1_indices] <- 0


    A_mat[CFBI2_index,model_CFBI2_indices] <- A_submat_CFBI2
    A_mat[CFBI2_index,-model_CFBI2_indices] <- 0


    A_mat[MCFP1_indices,model_CFBI1_indices] <- A_submat_MCFP1_sub
    A_mat[MCFP1_indices,MCFP1_indices] <- A_submat_MCFP_extra
    A_mat[MCFP1_indices,-model_CFBI_MCFP1_indices] <- 0

    A_mat[MCFP2_indices,model_CFBI2_indices] <- A_submat_MCFP2_sub
    A_mat[MCFP2_indices,MCFP2_indices] <- A_submat_MCFP_extra
    A_mat[MCFP2_indices,-model_CFBI_MCFP2_indices] <- 0

    num_groups <- length(AB_nuisance_list$model$B_list)
    B_vecs_list <- list()
    B_vec <- rep(NA,length=mat_size)
    for (ii in 1:num_groups){

      B_vec[model_indices] <- AB_nuisance_list$model$B_list[[ii]]
      B_vec[CFBI1_index] <- AB_nuisance_list$CFBI[[alpha1_num]]$B_list[[ii]]
      B_vec[MCFP1_indices] <- AB_nuisance_list$MCFP[[alpha1_num]]$B_list[[ii]]
      B_vec[CFBI2_index] <- AB_nuisance_list$CFBI[[alpha2_num]]$B_list[[ii]]
      B_vec[MCFP2_indices] <- AB_nuisance_list$MCFP[[alpha2_num]]$B_list[[ii]]

      ### defense
      if (any(!is.na(B_vec[length(B_vec)]))) {stop("all last entry needs to be NA")}
      if (any(is.na(B_vec[-length(B_vec)]))) {stop("only last entry can be NA")}
      ### defense

      B_vecs_list[[ii]] <- B_vec
    }
  }

  ### defense
  if (any(!is.na(A_mat[nrow(A_mat),]))) {stop("all last row needs to be NA")}
  if (any(is.na(A_mat[-nrow(A_mat),]))) {stop("only last row can be NA")}
  ### defense

  AB_list <- list(
    A_mat = A_mat,
    B_list = B_vecs_list
  )
}


# #' calculates bread and Meat of mu or CE function
# #'
# #' @inheritParams modularVar
# #' @inheritParams estimateNuisance
# #' @param deriv_opts simple or richardson
# #' @param target_row row to estimate
# #' @param CFBI_estimates do i have to do this myself
# #' @param MCFP_tidy_estimates do i have to do this myself
# #'
# #' @export
calcTargetBreadMeat <- function(
  split_data_list,
  model_parms  ,
  CFBI_estimates ,
  MCFP_tidy_estimates ,
  deriv_opts ,
  index_references,
  target_row,
  vec_dfm_no_alphas
){

  force(split_data_list)

  alpha1 <- target_row$alpha1
  alpha2 <- target_row$alpha2
  trt <- target_row$trt
  alpha1_num <- target_row$alpha1_num
  alpha2_num <- target_row$alpha2_num


  CFBI_alpha1 <- CFBI_estimates[alpha1_num]
  MCFP_alpha1 <- MCFP_tidy_estimates$MCFP[
    (MCFP_tidy_estimates$alpha == alpha1)
    ]

  if (!is.na(alpha2)){
    CFBI_alpha2 <- CFBI_estimates[alpha2_num]
    MCFP_alpha2 <- MCFP_tidy_estimates$MCFP[
      (MCFP_tidy_estimates$alpha == alpha2)
      ]
  } else {
    CFBI_alpha2 <- NULL #NULL's are not concatenated
    MCFP_alpha2 <- NULL
  }

  CFBI_ests <- c(CFBI_alpha1,CFBI_alpha2)
  MCFP_ests <- c(MCFP_alpha1,MCFP_alpha2)

  all_point_ests <- c(
    model_parms,
    CFBI_ests,
    MCFP_ests,
    target_row$estimate
  )


  AB_target_list <-
    calculateBreadMeat(
      geex_list = list(
        eeFUN = makeTargetEstFun,
        splitdt = split_data_list,
        ee_args = NULL
      ),
      # is_target = TRUE,
      theta =  all_point_ests,
      numDeriv_options = deriv_opts,
      index_references = index_references,
      trt = trt,
      alpha1 = alpha1,
      alpha2 = alpha2,
      vec_dfm_no_alphas = vec_dfm_no_alphas
    )
}

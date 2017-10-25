



# #' Estimating equation function for the model parameters
# #'
# #' @inheritParams makeTargetEstFun
# #' @inheritParams modular_MCFP_estfun_quick
# #'
# #' @export
modular_model_estfun <- function( data_list ){

  force(data_list)

  f2 <- function(theta){

    I_model_vec <- derivLogLike(
      model_parms = theta,
      model_DV_vec =  data_list$model_DV_vec, ## this is participation
      model_mat = data_list$model_matrix
    )
  }
}

# #' Estimating equation function for the CFBI parameters
# #'
# #' @inheritParams makeTargetEstFun
# #' @inheritParams modular_MCFP_estfun_quick
# #'
# #' @export
modular_CFBI_estfun <- function( data_list, alpha ){

  force(data_list)

  f2 <- function(theta){

    II_alpha1_CFBI <- groupAvgProbTrt(
      model_parms = theta[-length(theta)],
      CFBI = theta[length(theta)],
      model_mat = data_list$model_matrix
    )
    CFBI_psi_alpha1 <- II_alpha1_CFBI - alpha

  }
}


# #' Estimating equation function for the MCFP parameters, quicker
# #'
# #' @inheritParams makeTargetEstFun
# #' @inheritParams estimateGroupTarget
# #' @param MCFP_onealpha estimates of MCFP for one alpha policy
# #' @inheritParams estimateNuisance
# #' @param alpha the policy
# #' @param num_model_parms this could be replaced with index_refernces
# #'
# #' @export
modular_MCFP_estfun_quick <- function(
  data_list,
  alpha,
  num_model_parms,
  vec_dfm_no_alphas,
  MCFP_onealpha
){

  force(data_list)


  f2 <- function(theta){

    theta_model <- theta[1:num_model_parms]
    theta_CFBI <- theta[num_model_parms+1]
    theta_MCFP <- MCFP_onealpha


    III_alpha1_MCFP_vec <- MCFP_EstFun_LHS(
      model_parms = theta_model,
      CFBI = theta_CFBI,
      vec_dfm_no_alphas = vec_dfm_no_alphas,
      model_DV_vec = data_list$model_DV_vec,
      model_mat = data_list$model_matrix
    )

    MCFP_psi_vec_alpha1 <- III_alpha1_MCFP_vec - theta_MCFP
    MCFP_psi_vec_alpha1[is.na(MCFP_psi_vec_alpha1)] <- 0
    MCFP_psi_vec_alpha1

  }
}





# #' Estimating equation function for just the CE parameter
# #'
# #' @inheritParams modularVar
# #' @inheritParams targetEstFunLHS
# #' @inheritParams calcYBar
# #' @inheritParams estimateGroupTarget
# #' @inheritParams makeIndexReferences
# #'
# #' @export
makeTargetEstFun <- function(
  data_list,
  alpha1,
  alpha2,
  trt,
  index_references,
  vec_dfm_no_alphas
){


  force(data_list)

  outcome_vec <- data_list$outcome_vec
  treatment_vec <- data_list$treatment_vec
  model_DV_vec <- data_list$model_DV_vec
  model_mat <- data_list$model_matrix


  if (!is.na(alpha2)){
    model_indices <- index_references$CE_refs$model
    CFBI_alpha1_index <- index_references$CE_refs$CFBI1
    CFBI_alpha2_index <- index_references$CE_refs$CFBI2
    CFBI_end_index <- CFBI_alpha2_index


    vec_dfm_alp1 <- vec_dfm_alp2 <- vec_dfm_no_alphas
    vec_dfm_alp1$alpha <- alpha1
    vec_dfm_alp2$alpha <- alpha2
    VDA <- rbind(vec_dfm_alp1,vec_dfm_alp2)
    rm(vec_dfm_alp1,vec_dfm_alp2)
  } else {
    model_indices <- index_references$mu_refs$model
    CFBI_alpha1_index <- index_references$mu_refs$CFBI
    CFBI_end_index <- CFBI_alpha1_index

    VDA <- vec_dfm_no_alphas
    VDA$alpha <- alpha1
  }



  f2 <- function(theta){

    MCFP_alpha1 <- getMCFPFromTheta(
      treatment_vec = treatment_vec,
      vec_dfm_alphas = VDA,
      alpha = alpha1,
      is_quick_deriv =TRUE,
      CFBI_end_position = CFBI_end_index,
      theta = theta
    )



    if (!is.na(alpha2)) {
      MCFP_alpha2 <- getMCFPFromTheta(
        treatment_vec = treatment_vec,
        vec_dfm_alphas = VDA,
        alpha = alpha2,

        is_quick_deriv =TRUE,

        CFBI_end_position = CFBI_alpha2_index,
        theta = theta
      )
      MCFP <- MCFP_alpha1-MCFP_alpha2
    } else {
      MCFP <- MCFP_alpha1
    }

    target_val <- targetEstFunLHS(
      model_parms = theta[model_indices],
      outcome_vec = outcome_vec,
      treatment_vec = treatment_vec,
      model_mat = model_mat,
      trt = trt,
      MCFP = MCFP,
      CPS = NULL
    )

    if (length(target_val)!=1) {browser(); stop("should be length 1")}

    QE_is_NA <- is.na(
      calcYBar(
        outcome_vec = outcome_vec,
        treatment_vec = treatment_vec,
        trt = trt
      )
    )
    if (QE_is_NA) {
      browser();
      stop("QE should not be NA")
      target_psi_vec <- 0
      names(target_psi_vec) <- "" ## for consistency
    } else {
      target_psi_vec <-  target_val - theta[length(theta)]
    }


    ### defense
    if (length(target_psi_vec)!=1) {
      browser(); stop("should be length 1")}
    ### defense

    target_psi_vec
  }

  f2
}







## ## ## ## ## ## ## ## ## ## ## ##
##### CFBI-only geexable EstFun #####




# #' Full estimating function for CFBI targets
# #'
# #' @inheritParams groupAvgProbTrt
# #' @inheritParams estimateGroupTarget
# #' @param alpha the policy
# #'
# #' @export
CFBI_only_EstFun <- function(
  data_list,
  model_parms,
  alpha
){
  force(data_list)
  force(model_parms)
  force(alpha)

  f2 <- function(theta){

    group_CFBI_contrib <- groupAvgProbTrt(
      model_parms = model_parms,
      CFBI = theta,
      model_mat = data_list$model_matrix
    )

    CFBI_psi <- group_CFBI_contrib - alpha
  }

}


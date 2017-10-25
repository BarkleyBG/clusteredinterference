

## ## ## ## ## ## ## ## ## ## ## ##
##### -2. get MCFP subset vec #####

# #' translates one group's treatment (via its MCFP indices) into the correct MCFP from the MCFP_tidy_estimates dfm.
# #'
# #' @param treatment_vec the vector of treatment to evaluate (info on n and s)
# #' @param alpha the policy of interest (info on alpha)
# #' @param MCFP_tidy_estimates the dataframe to grab this all from
# #'
# #' @export
getMCFPFromTidy <- function(
  treatment_vec,  alpha,  MCFP_tidy_estimates
){

  MCFP_indices <- getIndices4MCFP(
    treatment_vec = treatment_vec,
    vec_dfm_alphas = MCFP_tidy_estimates,
    alpha = alpha
  )

  if (all(MCFP_indices >=0 )) {
    MCFP <- as.vector(MCFP_tidy_estimates$MCFP %*% MCFP_indices)
  } else
    if (all(MCFP_indices <= 0 )) {
      MCFP_inv_factor <- MCFP_indices[MCFP_indices!=0][1]
      size_of_stratum <- -1/MCFP_inv_factor

      MCFP <- 1 + as.vector(MCFP_tidy_estimates$MCFP %*% (MCFP_indices*size_of_stratum))
      MCFP <- MCFP/size_of_stratum
    } else {
      stop("all MCFP_indices must be in (-1,0) or (0,1)")
    }

  MCFP
}

# #' returns which indices (rows from vec_dfm) are for one group's MCFP
# #'
# #' @inheritParams getMCFPFromTidy
# #' @param vec_dfm_alphas something about estimating omegas
# #'
# #' @export
getIndices4MCFP <- function(
  treatment_vec,
  vec_dfm_alphas,
  alpha
){

  group_size <- length(treatment_vec)
  sum_trt <- sum(treatment_vec)
  stratum_size <- choose(group_size, sum_trt)

  n_rows <- vec_dfm_alphas$n==group_size
  s_rows <- vec_dfm_alphas$s==sum_trt
  alpha_rows <- vec_dfm_alphas$alpha==alpha
  correct_row <- n_rows * s_rows * alpha_rows

   if (all(correct_row==0)) {
    possible_rows <- n_rows * alpha_rows
    if (all(vec_dfm_alphas$all_probs[
      which(possible_rows==1)] == "removed_one")) {
      MCFP_indices <- -1*possible_rows / stratum_size
    } else {  stop('no possible rows')  }   } else
    if (sum(correct_row)==1) {
      MCFP_indices <- correct_row / stratum_size
    } else {   stop('too many correct rows')    }

  ### defense
  if (!(all(MCFP_indices %in% ((0:1)/stratum_size) ))) {
    if (!(all(MCFP_indices %in% ((-1:0)/stratum_size )))) {
      stop("all MCFP indices must be (-1,0)/stratum_size or (0,1)/stratum_size")
    }}
  if (length(MCFP_indices)==0){ stop("MCFP indices must be positive length")}
  ### defense

  MCFP_indices
}

# #' translates one group's treatment (via its MCFP indices) into the correct theta value corresponding to that MCFP
# #'
# #' @inheritParams getMCFPFromTidy
# #' @inheritParams getIndices4MCFP
# #' @param CFBI_end_position bookkeeping variable
# #' @param theta bookkeeping i think?
# #' @param is_quick_deriv TRUE for quick deriv, default FALSE
# #'
# #' @export
getMCFPFromTheta <- function(
  treatment_vec,
  vec_dfm_alphas,
  alpha,
  is_quick_deriv = FALSE,
  CFBI_end_position,
  theta
){
  theta_length <- length(theta)

  MCFP_indices <- getIndices4MCFP(
    treatment_vec = treatment_vec,
    vec_dfm_alphas = vec_dfm_alphas,
    alpha = alpha
  )

  MCFP_indices_vec <- rep(0,theta_length)

  MCFP_start <- CFBI_end_position + 1
  MCFP_end <-  CFBI_end_position + length(MCFP_indices)
  if (theta_length != (MCFP_end +  1)) {
    stop('incorrect nrows of dfm or length of vec in alp1')
  }
  if (length(MCFP_indices) != (MCFP_end - MCFP_start + 1)) {
    stop("theta wont line up with indices")
  }

  MCFP_indices_vec[MCFP_start:MCFP_end] <- MCFP_indices

  if (all(MCFP_indices_vec >=0 )) {
    MCFP <- as.vector(theta %*% MCFP_indices_vec)
  } else
     if (all(MCFP_indices_vec <= 0 )) {
      MCFP_inv_factor <- MCFP_indices_vec[MCFP_indices_vec!=0][1]
      size_of_stratum <- -1/MCFP_inv_factor

      MCFP <- 1 + as.vector(theta %*% (MCFP_indices_vec*size_of_stratum))
      MCFP <- MCFP/size_of_stratum
    } else {
      stop("all MCFP_indices must be in (-1,0) or (0,1)")
    }
  return(MCFP)
}

# #' Get a joint probability
# #'
# #' @param fixefs the betas
# #' @param ranef_std_dev sigma
# #' @param model_mat L
# #' @param mult_factor keep this at 1 i think it's deprecated
# #' @param binary_vec what we're trying to find the joint probability of
# #' @inheritParams argTests
# #'
# #' @export
integrateIntegrand <- function(
  fixefs,
  ranef_std_dev,
  model_mat,
  binary_vec,
  mult_factor=1
){

  ### defense
  if (length(fixefs) != ncol(model_mat)) {stop(
    "number of columns in model_mat must equal number of fixed effects"
  )}
  if (length(ranef_std_dev) != 1) {stop(
    "ranef_std_dev must be of length 1"
  )}
  if (!is.matrix(model_mat)) {stop(
    "model_mat must be a matrix"
  )}
  if (!all(binary_vec %in% 0:1)) {
    stop('supplied binary_vec is not binary')
  }
  # if(mult_factor!=1){stop("mult_factor is deprecated")}
  ### defense

  integral <- try(cubature::adaptIntegrate(
    f = modelIntegrand,
    lowerLimit = -5*ranef_std_dev,
    upperLimit = 5*ranef_std_dev,
    fixefs = fixefs,
    ranef_std_dev = ranef_std_dev,
    model_mat = model_mat,
    binary_vec = binary_vec,
    mult_factor =mult_factor
  ))

  if("try-error" %in% class(integral)){
    return(0)
  } else if (is.na(integral$integral)){
    return(0)
  }

  integral$integral
}

# #' This function creates the integrand for the GLMER-type joint probabilities
# #'
# #' @inheritParams integrateIntegrand
# #' @param random_intercept this gets passed in via integration; the b's
# #'
# #' @export
modelIntegrand <- function(
  random_intercept,
  fixefs,
  ranef_std_dev,
  model_mat,
  binary_vec,
  mult_factor =1,
  RZC=1
){

  lin_pred <- model_mat%*%fixefs + random_intercept
  eta <- stats::plogis(lin_pred)

  if(mult_factor!=1){
    nthroot_mult_factor <- mult_factor^(1/length(binary_vec))
    eta_bern <- nthroot_mult_factor*(RZC*eta)^(binary_vec)*(1-(RZC*eta))^(1-binary_vec)
  } else {
    eta_bern <- (RZC*eta)^(binary_vec)*(1-(RZC*eta))^(1-binary_vec)
  }


  integrand <- prod(eta_bern)*
    stats::dnorm(random_intercept, mean=0, sd = ranef_std_dev)

}

# #' The scores? Derivative of log likelihood of treatment vectors
# #'
# #' @inheritParams MCFP_EstFun_LHS
# #' @param deriv_loglihood passed to numDeriv. 'simple' or 'Richardson'
# #'
# #' @export
derivLogLike <- function(
  model_parms,
  model_DV_vec,
  model_mat,
  deriv_loglihood = 'simple'
){
  deriv <- numDeriv::grad(
    func = modelLoglihood,
    x = model_parms,
    method = deriv_loglihood,
    model_mat = model_mat,
    model_DV_vec = model_DV_vec
  )
}

# #' computes the log_likelihood of a treatment vector
# #'
# #' @inheritParams MCFP_EstFun_LHS
# #' @param x model parameters (beta,sigma) passed in from derivLogLike
# #'
# #' @export
modelLoglihood <- function( x, model_mat, model_DV_vec ){
  ranef_std_dev <- x[length(x)]
  fixefs <- x[-length(x)]

  model_likelihood <- integrateIntegrand(
    fixefs = fixefs,
    ranef_std_dev = ranef_std_dev,
    model_mat = model_mat,
    binary_vec = model_DV_vec
  )

  log(model_likelihood)

}



## ## ## ## ## ## ## ## ## ## ## ## ##
## ##### II. CFBI EstFun builders #####

# #' One group's contribution to the objective function (and subtract alpha)
# #'
# #' @inheritParams MCFP_EstFun_LHS
# #'
# #' @export
groupAvgProbTrt <- function(model_parms,  CFBI, model_mat){

  ranef_std_dev <- model_parms[length(model_parms)]
  fixefs <- model_parms[-length(model_parms)]
  fixefs[1] <- CFBI
  CFBI_fixefs <- fixefs

  each_indiv_prob_trt <- lapply(
    1:nrow(model_mat), function(row_num){ ## per individual
      indiv_model_mat <- model_mat[row_num,,drop=FALSE]

      indiv_likelihood <- integrateIntegrand(
        fixefs = CFBI_fixefs,
        ranef_std_dev = ranef_std_dev,
        model_mat = indiv_model_mat,
        binary_vec = 1
      )
    })

  each_indiv_prob_trt_vec <- as.numeric(each_indiv_prob_trt)

  ### defense
  if (length(each_indiv_prob_trt_vec)!= nrow(model_mat)) {stop(
    "model_mat must have same number of rows as individual probabilites estimated"
  )}
  if (any(is.na(each_indiv_prob_trt_vec))) {stop(
    "NA values estimated for individual probabilities"
  )}
  ### defense

  group_avg_prob_trt <- mean(each_indiv_prob_trt_vec)
}

# ## ## ## ## ## ## ## ## ## ## ## ##
# ##### III. MCFP EstFun builders #####
#
# #' One group's contribution to the objective function for omegas
# #'
# #' @inheritParams estimateMCFP
# #' @param CFBI the estimated CFBI term
# #' @param model_parms theta, a vector of fixefs and raneff
# #' @param model_DV_vec the treatment vector
# #' @param model_mat the covariates (independent variables)
# #'
# #' @export
MCFP_EstFun_LHS <- function(
  model_parms,
  CFBI,
  vec_dfm_no_alphas,
  model_mat,
  model_DV_vec
){
  ranef_std_dev <- model_parms[length(model_parms)]
  fixefs <- model_parms[-length(model_parms)]
  fixefs[1] <- CFBI
  CFBI_fixefs <- fixefs

  group_size <- length(model_DV_vec)
  if (group_size != nrow(model_mat)) {stop("Cluster should be correct size")}
  sum_trt <- sum(model_DV_vec)
  if ((sum(1-model_DV_vec)+sum_trt) != group_size) {stop("data inconsistencies")}


  to_output <- lapply(
    1:nrow(vec_dfm_no_alphas),
    function(vec_dfm_row_num){
      vec_dfm_row <- vec_dfm_no_alphas[vec_dfm_row_num,]


      if (vec_dfm_row$n != group_size) {
        return(NA_real_) ## different group size.
      }

      vecs_to_eval_mat <- vec_dfm_row$trt_vecs_mats[[1]]
      mult_factor <- vec_dfm_row$mult_factor
      each_vec_contribution <- apply(
        X = vecs_to_eval_mat,
        MARGIN = 1,
        function(vec_to_eval){

          ### defense
          if (!is.vector(vec_to_eval)) {stop(
            'vec_to_eval must be a vector'
          )}
          if (!all(vec_to_eval%in% 0:1)) {stop(
            'vec_to_eval must be binary'
          )}
          ### defense

          group_likelihood <- integrateIntegrand(
            fixefs = CFBI_fixefs,
            ranef_std_dev = ranef_std_dev,
            model_mat = model_mat,
            binary_vec = vec_to_eval,
            mult_factor = mult_factor
          )
        })

      ### defense
      if (length(each_vec_contribution)!= nrow(vecs_to_eval_mat)) {
        stop('this could be right check the whole list-to-vec thing')
      }
      ### defense

      vec_contributions <- unlist(each_vec_contribution)

      ### defense
      if (!is.numeric(vec_contributions)) {stop(
        'vec_contributions must be numeric'
      )}
      if (any(1 <= vec_contributions)) {stop(
        'vec_contributions must be a vector of non-1 probabilities'
      )}
      if (length(vec_contributions)!= nrow(vecs_to_eval_mat)) {
        stop('this could be right; check the whole list-to-vec thing')
      }
      ### defense

      #this is sum_{a \in A(n,s)}Pr(A=a|L=l) = Pr(s(A) =s|L=l) / size(A(n,s))
      sum_all_MCFP_prob_ests <-  sum(vec_contributions)
    })

  group_MCFP_contributions_vec <- unlist(to_output) ## check if this works

  ### defense
  if (!is.numeric(group_MCFP_contributions_vec)) {stop(
    'group_MCFP_vec_contributions must be numeric'
  )}
  if (any(1 <= stats::na.omit(group_MCFP_contributions_vec))) {stop(
    'group_MCFP_vec_contributions must be a vector of non-1 probabilities'
  )}
  if (nrow(vec_dfm_no_alphas)!= length(group_MCFP_contributions_vec)) {stop(
    'group_MCFP_vec_contributions must have one prob for each row in vec_dfm_no_alphas'
  )}
  ### defense

 group_MCFP_contributions_vec
}



## #' One group's contribution to the Objective function for target estimands
## #'
## #' @inheritParams MCFP_EstFun_LHS
## #' @param MCFP the estimated omegas
## #' @param return_prop_score_only default FALSE
## #' @param treatment_vec the treatment vector
## #' @param outcome_vec the outcome vector
## #' @param trt equals to NA, 1 or 0
## #' @param CPS the cluster propensity score, for estimation not variance
## #'
## #' @export
targetEstFunLHS <- function(
  model_parms,
  outcome_vec,
  treatment_vec,
  model_mat,
  trt,
  CPS=NULL,
  return_prop_score_only = FALSE,
  MCFP
){
  if (!is.null(CPS)) {
    prop_score <- CPS
  } else{
    ranef_std_dev <- model_parms[length(model_parms)]
    fixefs <- model_parms[-length(model_parms)]

    prop_score <- integrateIntegrand(
      fixefs = fixefs,
      ranef_std_dev = ranef_std_dev,
      model_mat = model_mat,
      binary_vec = treatment_vec
    )

    if (return_prop_score_only) {
      return(prop_score)
    }
  }
  ### defense
  if (!is.vector(outcome_vec)) {
    stop("Outcome should be a vector")
  }
  if (!(all(outcome_vec == as.integer(outcome_vec)))) {
    stop("Outcome should be an integer")
  }
  if (is.na(prop_score)||prop_score==0){ prop_score <- 1}
  ### defense


  avg_outcome <- calcYBar(
    outcome_vec = outcome_vec,treatment_vec = treatment_vec, trt=trt)

  target_LHS_contribution <- avg_outcome * MCFP / prop_score ## may be NA

  ### defense
  if (length(target_LHS_contribution)!=1) { stop("Should return scalar")}
  ### defense

  target_LHS_contribution
}




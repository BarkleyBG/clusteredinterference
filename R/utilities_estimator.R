

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


# #' Make Useful References for EstFun Var
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






###########################################
###### makeTargetGrid of estimates ########
###########################################

## from Bradley Saul's inferference

#' Creates the grid of all estimands to estimate.
#'
#' @inheritParams policyFX
#' @param small_grid if TRUE then only estimates some policies. Default is FALSE.
#'
#' @export
makeTargetGrid <- function(alphas, small_grid = FALSE){

  trts <- c(NA_real_, 0, 1)

  mu_grid <- expand.grid(
    alpha1_num = 1:length(alphas),
    alpha2_num = NA_real_,
    trt = trts)

  for(row_num in 1:nrow(mu_grid)){
    mu_grid$alpha1[row_num] <- alphas[mu_grid$alpha1_num[row_num]]
    mu_grid$alpha2[row_num] <- NA
    mu_grid$estimand[row_num] <- ifelse(
      is.na(mu_grid$trt[row_num]), "mu", paste0("mu", mu_grid$trt[row_num]))
    mu_grid$effect_type[row_num] <- 'outcome'
  }


  CE_grid <- expand.grid(
    alpha1_num = 1:length(alphas),
    alpha2_num = 1:length(alphas),
    trt = trts)

  for(row_num in 1:nrow(CE_grid)){
    CE_grid$alpha1[row_num] <- alphas[CE_grid$alpha1_num[row_num]]
    CE_grid$alpha2[row_num] <- alphas[CE_grid$alpha2_num[row_num]]
    CE_grid$estimand[row_num] <- ifelse(
      is.na(CE_grid$trt[row_num]), "CE", paste0("QE", CE_grid$trt[row_num]))
    CE_grid$effect_type[row_num] <- 'contrast'
  }


  out_grid <- rbind(mu_grid,CE_grid)
  rownames(out_grid) <- NULL


  if (isTRUE(small_grid)) {out_grid <- out_grid[
    (is.na(out_grid$alpha2)|(out_grid$alpha1 > out_grid$alpha2)),
    ]}

  out_grid$estVar <- TRUE
  out_grid
}

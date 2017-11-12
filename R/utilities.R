
# #' Tidies the output estimates
tidyEstimates <- function(grid){
    grid$alpha1_num <- NULL
  grid$alpha2_num <- NULL
  grid$estVar <- NULL
  grid$estimand_type <- grid$estimand
  grid$k_samps <- grid$k
  grid$k <- NULL
  grid$estimand <- NA

  for (ii in 1:NROW(grid)){
    ## Quick fix for estimand names - issue #1
    grid$estimand_type[ii] <- renameEstimand(grid$estimand_type[ii])

    grid$estimand[ii] <- renameEstimandLong(
      estimand_type = grid$estimand_type[ii],
      alpha1 = grid$alpha1[ii],
      alpha2 = grid$alpha2[ii]
    )
  }
  if (!"var" %in% colnames(grid)) {grid$var <- NA}
  if (!"se"  %in% colnames(grid)) {grid$se  <- NA}
  if (!"LCI" %in% colnames(grid)) {grid$LCI <- NA}
  if (!"UCI" %in% colnames(grid)) {grid$UCI <- NA}
  colnames_ordered <- c(
    "estimand",
    "estimate", "var", "se", "LCI", "UCI",
    "alpha1", "alpha2", "trt",
    "estimand_type", "effect_type",
    "k_samps"
  )
  out <- grid[,colnames_ordered]
  row.names(out) <- NULL
  out
}

## Renaming estimands from old to new names. Quick fix for issue #1.
renameEstimand <- function(estimand_type){
  stopifnot(length(estimand_type)==1)
  num_chars <- nchar(estimand_type)
  new_sub <- changeSubEstimand(estimand_type)

  ## return new name
  paste0(new_sub,substr(estimand_type,3,max(num_chars,3)))
}

changeSubEstimand <- function(estimand_type){
  old_sub <- substr(estimand_type,0,2)
  if (old_sub %in% c("CE", "OE")) {new_sub <- "OE"} else
  if (old_sub %in% c("QE", "SE")) {new_sub <- "SE"} else
  if (old_sub == "mu") {new_sub <- "mu"} else {
    stop("no more possible substrings")
  }
  new_sub
}

renameEstimandLong <- function(estimand_type, alpha1, alpha2){
  stopifnot(length(estimand_type)==1)
  stopifnot(length(alpha1)==1)

  new_estimand <- renameEstimand(estimand_type)

  ifelse(
    substr(new_estimand,0,2)=="mu",
    paste0(new_estimand, "(", alpha1,  ")"),
    paste0(new_estimand, "(", alpha1, ",", alpha2, ")")
  )
}

# #' checkpointer
# #'
# #' @param time1  new time
# #'
# #' @export
checkpointTimer <- function(time1, verbose = FALSE){
  time2 <- Sys.time()
  dtime <- difftime(time2, time1, units="mins")
  if (verbose) {
    print(paste0("elapsed: ", round(dtime,2),
                 " minutes as of ",time2))
  }
  time2
}


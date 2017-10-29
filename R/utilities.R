
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
    grid$estimand[ii] <-
      ifelse(substr(grid$estimand_type[ii],0,2)=="mu",
             paste0(grid$estimand_type[ii], "(",
                    grid$alpha1[ii],  ")"
             ),
             paste0(grid$estimand_type[ii], "(",
                    grid$alpha1[ii], ",", grid$alpha2[ii], ")"
             )
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


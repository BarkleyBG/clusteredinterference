



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


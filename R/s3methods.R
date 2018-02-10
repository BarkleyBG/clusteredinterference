
#' Prints a summary of the estimates from a policyFX object
#'
#' @param x object of class "policyFX"
#' @param ... User may specify integer \code{nrows}.
#'
#' @method print policyFX
#'
#' @export
print.policyFX <- function(x,...){

  ests <- x$estimates
  keep_rows <- (ests$alpha1 != ests$alpha2)
  keep_rows[is.na(keep_rows)] <- 1
  ests <- ests[which(keep_rows==1), ]
  tot_rows <- NROW(ests)

  dots <- list(...)
  if (!("nrows" %in% names(dots))){
    nrows <- NROW(ests[is.na(ests$alpha2),]) + 2
  } else{
    nrows <- dots$nrows
  }
  nrows <- min(nrows, tot_rows)

  nrows_more <- tot_rows - nrows

  out_dfm <- ests[
    1:nrows,
    c("estimand", "estimate", "se", "LCI", "UCI")
    ]


  cat("------------- causal estimates --------------\n")
  print(out_dfm, row.names = FALSE, digits = 3)
  # cat('\n')
  if (nrows_more>0){
    cat('          ... and', nrows_more, 'more rows ...', "\n")
  }
  cat("---------------------------------------------\n")
}


#' Prints a summary of a policyFX object
#'
#' @param object object of class "policyFX"
#' @param ... User may specify integer \code{nrows}.
#'
#' @method summary policyFX
#'
#' @author Brian G. Barkley
#'
#' @export
summary.policyFX <- function(object, ...){


  ests <- object$estimates
  keep_rows <- (ests$alpha1 != ests$alpha2)
  keep_rows[is.na(keep_rows)] <- 1
  ests <- ests[which(keep_rows==1), ]
  tot_rows <- NROW(ests)

  dots <- list(...)
  if (!("nrows" %in% names(dots))){
    nrows <- NROW(ests[is.na(ests$alpha2),]) + 2
  } else{
    nrows <- dots$nrows
  }
  nrows <- min(nrows, tot_rows)

  nrows_more <- tot_rows - nrows

  out_dfm <- ests[
    1:nrows,
    c("estimand", "estimate", "se", "LCI", "UCI")
    ]

  model <- object$model
  # ps <- as.vector(object$propensity_scores$CPS)
  ps <- sprintf("%.3g",object$propensity_scores$CPS)
  names(ps) <- object$propensity_scores$cluster_name

  cat("------------- causal estimates --------------\n")
  print(out_dfm, row.names = FALSE, digits = 3)
  cat('\n')
  if (nrows_more>0){
    cat('          ... and', nrows_more, 'more rows ...', "\n")
  }
  cat('\n')
  cat("-------------- treatment model -------------\n")

  print(model)
  cat('\n')
  cat("------------- propensity scores -------------\n")

  print(ps, quote = FALSE)
  cat("---------------------------------------------\n")
}

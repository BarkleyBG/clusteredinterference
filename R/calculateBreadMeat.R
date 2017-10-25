

# #' Play Nicely with Arrays
# #'
# #' @param object array numeric matrix or o/w
# #'
# #' @export
checkArray <- function (object) {
  if (is.array(object)) {
    object
  }
  else if (is.numeric(object)) {
    array(object, dim = c(1, 1, length(object)))
  }
  else if (is.matrix(object)) {
    array(object, dim = c(1, 1, length(object)))
  }
  else {
    stop("Object is not an array, matrix, or numeric")
  }
}


# #' Calculations for the Bread Matrix (A) and Meat Matrix (B)
# #'
# #' @description this is designed to be modular
# #'
# #' @param geex_list estimating function and the split data
# #' @param theta Estimates to take derivs at
# #' @param numDeriv_options numDeriv options
# #' @param calc_which "Both" default, or "Bread" or "Meat" only
# #' @param ... passed into the Amat and B mat calcs
# #'
# #' @export
calculateBreadMeat <- function(
  geex_list,
  theta,
  calc_which="Both",
  numDeriv_options = list(method = "Richardson"),
  silent = TRUE, ...) {

  out_list <- list()
  if (is.null(geex_list$ee_args)) {
    ee_args <- NULL
  }
  ### defense
  if (!calc_which %in% c("Both", "Bread", "Meat")) {
    stop("calc_which must equal Both, Bread, or Meat ")
  }
  ### defense

  if (calc_which != "Meat") {
    A_mat <- with(geex_list, {
      psi_i <- lapply(splitdt, function(data_list_i) {
        force(data_list_i)

        eeFUN(data_list = data_list_i, ...)
      })
      A_i <- lapply(psi_i, function(ee) {
        force(ee)
        args <- append(list(fun = ee, x = theta), numDeriv_options)
        val <- do.call(numDeriv::jacobian, args = append(args,
                                                         ee_args))
        -val
      })
      A_i_array <- checkArray(simplify2array(A_i))
      A <- apply(A_i_array, 1:2, sum)

    })

     out_list$A_mat <- A_mat
  }
  if (calc_which != "Bread") {

    B_vecs_list <- with(geex_list, {
      psi_i <- lapply(splitdt, function(data_list_i) {
        force(data_list_i)

        eeFUN(data_list = data_list_i,...)
      })

      B_i <- lapply(psi_i, function(ee) {
        force(ee)
        ee_val <- do.call(ee, args = append(list(theta = theta), ee_args))
      })
    })
    out_list$B_list <-  B_vecs_list
  }

  out_list
}




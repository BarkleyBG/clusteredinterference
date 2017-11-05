
context('Estimating gamma_{0,alpha}, the CFBI')


load(file = quickLookup("test_CFBI.Rdata"))

est_CFBI1Alpha_args$split_data_list <-
  est_CFBI1Alpha_args$split_data_list[1:7]

set.seed(323232)
CFBI_alpha_solved <- do.call(getCFBI1Alpha,est_CFBI1Alpha_args)


testthat::test_that(
  desc = "getCFBI1Alpha() solves for the same gamma_{0,alpha} value",
  testthat::expect_equal(
    object =  CFBI_alpha_solved,
    expected = CFBI_alpha_baseline,
    tolerance = 1e-9,
    check.attributes = TRUE
  )
)

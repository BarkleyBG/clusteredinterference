
context('Estimate CFBI and MCFP values')

load(file = quickLookup("test_CFBI_MCFP.Rdata"))

suppressWarnings(RNGversion("3.5.0"))
set.seed(101010)

CFBI_MCFP_list <- do.call(what=estimateNuisance,args=CFBI_args)


testthat::test_that(
  desc = "Nuisance parameter estimation returns as expected",
  testthat::expect_equal(
    object =  CFBI_MCFP_list,
    expected = CFBI_MCFP_orig,
    tolerance = helper_tol,
    check.attributes = TRUE
  )
)

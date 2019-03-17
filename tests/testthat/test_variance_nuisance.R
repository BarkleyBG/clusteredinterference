
context('Estimating variance of nuisance parameters')

load(file=quickLookup("test_variance_nuisance.Rdata"))

suppressWarnings(RNGversion("3.5.0"))
set.seed(101010)

AB_nuisance_list <- do.call(calcNuisanceBreadMeat,nuisance_args)

testthat::test_that(
  desc = "Test that old code and new code produce same target estimates",
  testthat::expect_equal(
    object =  AB_nuisance_list,
    expected = AB_nuisance_list_orig,
    tolerance = helper_tol,
    check.attributes = TRUE
  )
)

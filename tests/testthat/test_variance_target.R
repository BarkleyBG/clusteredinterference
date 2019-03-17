
context('Estimating variance for target parameters')


load(file = quickLookup("test_variance_target_mu.Rdata"))

target_args$RZC <- NULL
suppressWarnings(RNGversion("3.5.0")) ## For backwards compatibility
set.seed(101010)
AB_target_list_mu <- do.call(calcTargetBreadMeat, args = target_args)

testthat::test_that(
  desc = "Bread and meat of mu return as expected",
  testthat::expect_equal(
    object =  AB_target_list_mu,
    expected = AB_target_list_orig,
    tolerance = 1e-9,
    check.attributes = FALSE
  )
)



load(file = quickLookup("test_variance_target_OE.Rdata"))

target_args$RZC <- NULL
set.seed(101010)
AB_target_list_OE <- do.call(calcTargetBreadMeat, args = target_args)

testthat::test_that(
  desc = "Bread and meat of overall effect return as expected",
  testthat::expect_equal(
    object =  AB_target_list_OE,
    expected = AB_target_list_orig,
    tolerance = 1e-9,
    check.attributes = FALSE
  )
)


context('Prepare CFBI solver')

load(file = quickLookup("test_CFBI_starting.Rdata"))

suppressWarnings(RNGversion("3.5.0"))
set.seed(101010)

starting_dfm <- getStartingCFBI(
  new_param_list = new_param_list,
  glmer_fit = glmer_fit,
  alphas = alphas,
  CFBI_numstart = CFBI_numstart
)


testthat::test_that(
  desc = "Test that starting CFBI dataframes are identical",
  testthat::expect_equal(
    object =  starting_dfm,
    expected = starting_dfm_orig,
    tolerance = helper_tol,
    check.attributes = TRUE
  )
)

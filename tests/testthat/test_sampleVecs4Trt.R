
context('Sampling vectors for valuation for MCFP/omegas')

load(file = quickLookup("test_sampleVecs4Trt.Rdata"))

suppressWarnings(RNGversion("3.5.0"))
set.seed(309)

vec_dfm_no_alphas <- sampleVecs4Trt(
  split_data_list = split_data_list,
  max_k_evals = max_k_evals
)


testthat::test_that(
  desc = "sampleVecs4Trt() returns same treatment vectors to evaluate",
  code = testthat::expect_identical(
    object =  vec_dfm_no_alphas[,1:4],
    expected = vec_dfm_no_alphas_orig
  )
)

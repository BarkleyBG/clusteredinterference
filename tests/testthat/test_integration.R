
context('Integration test')


filename <- quickLookup("test_integration.Rdata")
load(file=filename) ## loads est_args and baseline_effects

suppressWarnings(RNGversion("3.5.0"))
set.seed(101010)
pfx <- do.call( policyFX, args = est_args )

new_ests <- pfx$estimates
old_ests <- baseline_effects$estimates

for (ii in 1:nrow(old_ests)) {
  old_ests$estimand[[ii]] <- renameEstimand(old_ests$estimand[[ii]])
  # old_ests$estimand_type[[ii]] <- renameEstimand(old_ests$estimand_type[[ii]])
}

# OE2 <- old_ests[is.na(old_ests$trt),]
# NE2 <- new_ests[is.na(new_ests$trt), ]
#
# OEtrt <- old_ests[!is.na(old_ests$trt),]
# NEtrt <- new_ests[!is.na(new_ests$trt), ]
# NE2[, -(1:7)]

# testthat::test_that(
#   desc = "Whole estimation is identical",
#   testthat::expect_equal(
#     object =  OE2,
#     expected = NE2,
#     tolerance = helper_tol,
#     check.attributes = FALSE
#   )
# )
#
# testthat::test_that(
#   desc = "Whole estimation is identical for spillovers",
#   testthat::expect_equal(
#     object =  OEtrt,
#     expected = NEtrt,
#     tolerance = helper_tol,
#     check.attributes = FALSE
#   )
# )



testthat::test_that(
  desc = "Target ests are identical",
  testthat::expect_equal(
    object =  old_ests$estimate,
    expected = new_ests$estimate,
    tolerance = helper_tol,
    check.attributes = TRUE
  )
)

testthat::test_that(
  desc = "Target SEs are identical",
  testthat::expect_equal(
    object =  old_ests$se,
    expected = new_ests$se,
    tolerance = helper_tol,
    check.attributes = TRUE
  )
)
testthat::test_that(
  desc = "Model parms estimates are identical",
  testthat::expect_equal(
    object =  baseline_effects$parameters[[1]],
    expected = pfx$parameters[[1]],
    tolerance = helper_tol,
    check.attributes = TRUE
  )
)
testthat::test_that(
  desc = "CFBI  estimates are identical",
  testthat::expect_equal(
    object =  baseline_effects$parameters[[2]],
    expected = pfx$parameters[[2]],
    tolerance = helper_tol,
    check.attributes = TRUE
  )
)
testthat::test_that(
  desc = "MCFP estimates are identical",
  testthat::expect_equal(
    object =  baseline_effects$parameters[[3]],
    expected = pfx$parameters[[3]],
    tolerance = helper_tol,
    check.attributes = TRUE
  )
)
testthat::test_that(
  desc = "target estimates are identical",
  testthat::expect_equal(
    object =  baseline_effects$parameters[[4]],
    expected = pfx$parameters[[4]][,-(8:9)],
    tolerance = helper_tol,
    check.attributes = FALSE
  )
)
testthat::test_that(
  desc = "Var mats [1] are identical",
  testthat::expect_equal(
    object =  baseline_effects$variance_matrices[[1]],
    expected = pfx$variance_matrices[[1]],
    # tolerance = helper_tol,
    tolerance = 1e-5,
    check.attributes = TRUE
  )
)
testthat::test_that(
  desc = "Var mats [1]$Meat_mat are identical",
  testthat::expect_equal(
    object =  baseline_effects$variance_matrices[[1]]$Meat_mat,
    expected = pfx$variance_matrices[[1]]$Meat_mat,
    # tolerance = helper_tol,
    tolerance = 1e-5,
    check.attributes = TRUE
  )
)
testthat::test_that(
  desc = "Var mats [1]$Bread_mat are identical",
  testthat::expect_equal(
    object =  baseline_effects$variance_matrices[[1]]$Bread_mat,
    expected = pfx$variance_matrices[[1]]$Bread_mat,
    tolerance = helper_tol,
    check.attributes = TRUE
  )
)

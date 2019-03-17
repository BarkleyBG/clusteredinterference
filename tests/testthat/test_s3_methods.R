

context('Testing s3 print methods')


object <- readRDS(quickLookup("test_summary.Rds"))


testthat::test_that(
  "Print method outputs desired sections",
  {
    printed_output <- testthat::evaluate_promise(print(object), print = TRUE)

    testthat::expect_true(
      grepl("causal estimates", printed_output$output)
    )
  }
)

testthat::test_that(
  "Print method outputs desired sections",
  {

    printed_summary <- testthat::evaluate_promise(summary(object), print = TRUE)

    testthat::expect_true(
      grepl("causal estimates", printed_summary$output)
    )
    testthat::expect_true(
      grepl("treatment model", printed_summary$output)
    )
    testthat::expect_true(
      grepl("propensity scores", printed_summary$output)
    )
  }
)


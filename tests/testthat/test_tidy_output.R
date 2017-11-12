
context("tidying estimator output")


test_that(
  "Naming estimand_type is consistent", {
    expect_identical(renameEstimand("mu"),"mu")
    expect_identical(renameEstimand("mu1"),"mu1")
    expect_identical(renameEstimand("mu0"),"mu0")
    expect_identical(renameEstimand("CE"),"OE")
    expect_identical(renameEstimand("QE1"),"SE1")
    expect_identical(renameEstimand("QE0"),"SE0")
  }
)



test_that(
  "Naming estimand is consistent", {
    expect_identical(renameEstimand("mu(asdjklasd)"),"mu(asdjklasd)")
    expect_identical(renameEstimand("mu1(0.3333)"),"mu1(0.3333)")
    expect_identical(renameEstimand("mu0(foo,bar)"),"mu0(foo,bar)")
    expect_identical(renameEstimand("CE(0.22, 0.34)"),"OE(0.22, 0.34)")
    expect_identical(renameEstimand("QE1(0.5, 0.999999)"),"SE1(0.5, 0.999999)")
    expect_identical(renameEstimand("QE0(flarf)"),"SE0(flarf)")
  }
)



test_that(
  "Naming long estimand is consistent", {
    expect_identical(renameEstimandLong("mu", 1/2),"mu(0.5)")
    expect_identical(renameEstimandLong("mu1",0.2), "mu1(0.2)")
    expect_identical(renameEstimandLong("mu0", 0.33),"mu0(0.33)")
    expect_identical(renameEstimandLong("CE", 0.22, 0.34),"OE(0.22,0.34)")
    expect_identical(renameEstimandLong("QE1",0.5, 0.999999),"SE1(0.5,0.999999)")
    expect_identical(renameEstimandLong("QE0", 0.5, 1/2),"SE0(0.5,0.5)")
  }
)

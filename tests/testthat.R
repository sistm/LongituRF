library(testthat)

context("Tu vas me tester ça oui ?!")

test_that("Test about initialization", {
  expect_equal(length(runif(100)),100)
})

library(testthat)

# This test file lives at: modules/snagDecay/tests/testthat/
# The project root (containing modules/) is 4 levels up.
.projRoot <- normalizePath(file.path(getwd(), "..", "..", "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(.projRoot, "modules"))) {
  .projRoot <- getwd()
}
source(file.path(.projRoot, "modules", "snagDecay", "R", "transition.R"))

snagTransMat_test <- matrix(
  c(0.48, 0.38, 0.05, 0.00, 0.00,
    0.00, 0.52, 0.34, 0.05, 0.00,
    0.00, 0.00, 0.56, 0.30, 0.04,
    0.00, 0.00, 0.00, 0.62, 0.26,
    0.00, 0.00, 0.00, 0.00, 0.70),
  nrow = 5, byrow = TRUE
)

test_that("applyTransition returns integer vector of same length as input", {
  set.seed(42)
  result <- applyTransition(1:5, snagTransMat_test)
  expect_type(result, "integer")
  expect_length(result, 5L)
})

test_that("applyTransition output is always in 1:5", {
  set.seed(1)
  result <- applyTransition(rep(1:5, each = 100), snagTransMat_test)
  expect_true(all(result >= 1L & result <= 5L))
})

test_that("applyTransition never decreases DC (upper-triangular matrix)", {
  set.seed(99)
  dc_in <- rep(1:5, each = 500)
  dc_out <- applyTransition(dc_in, snagTransMat_test)
  expect_true(all(dc_out >= dc_in))
})

test_that("applyTransition DC5 input always returns DC5 (absorbing diagonal only)", {
  set.seed(7)
  result <- applyTransition(rep(5L, 200), snagTransMat_test)
  expect_true(all(result == 5L))
})

library(SpaDES.core)
library(data.table)
library(testthat)

snagTransMat_test <- matrix(
  c(0.48, 0.38, 0.05, 0.00, 0.00,
    0.00, 0.52, 0.34, 0.05, 0.00,
    0.00, 0.00, 0.56, 0.30, 0.04,
    0.00, 0.00, 0.00, 0.62, 0.26,
    0.00, 0.00, 0.00, 0.00, 0.70),
  nrow = 5, byrow = TRUE
)
snagFallProb_test <- c(DC1 = 0.09, DC2 = 0.09, DC3 = 0.10, DC4 = 0.12, DC5 = 0.30)

emptyCohorts <- data.table(
  pixelID = integer(), year = integer(),
  species = character(), B = numeric()
)

# testInit is not available in SpaDES.core >= 3.x; use simInit directly.
# This test file lives at: modules/DeadWood_snagDecay/tests/testthat/
# The project root (containing modules/) is 4 levels up.
.projRoot <- normalizePath(file.path(getwd(), "..", "..", "..", ".."), mustWork = FALSE)
# Validate; if the modules dir is not found 4 levels up, fall back to 0 levels (running from root)
if (!dir.exists(file.path(.projRoot, "modules"))) {
  .projRoot <- getwd()
}

testInit <- function(moduleName, params, objects, times = list(start = 0, end = 10)) {
  simInit(
    times   = times,
    modules = list(moduleName),
    params  = params,
    objects = objects,
    paths   = list(modulePath = file.path(.projRoot, "modules"))
  )
}

test_that("snagDecay init creates empty snagTable with correct schema", {
  sim <- testInit(
    "DeadWood_snagDecay",
    params = list(snagDecay = list(
      snagTransMat = snagTransMat_test,
      snagFallProb = snagFallProb_test,
      species      = "Pinus strobus"
    )),
    objects = list(cohortData = emptyCohorts)
  )
  sim <- spades(sim, events = "init")
  expect_s3_class(sim$snagTable, "data.table")
  expect_equal(nrow(sim$snagTable), 0L)
  expect_named(sim$snagTable, c("pixelID", "species", "DC", "ageInDC", "initBiomass"))
})

test_that("snagDecay init creates empty fallenSnags", {
  sim <- testInit(
    "DeadWood_snagDecay",
    params = list(snagDecay = list(
      snagTransMat = snagTransMat_test,
      snagFallProb = snagFallProb_test,
      species      = "Pinus strobus"
    )),
    objects = list(cohortData = emptyCohorts)
  )
  sim <- spades(sim, events = "init")
  expect_s3_class(sim$fallenSnags, "data.table")
  expect_equal(nrow(sim$fallenSnags), 0L)
  expect_named(sim$fallenSnags, c("pixelID", "species", "DC", "ageInDC", "initBiomass"))
})

test_that("snagDecay annual absorbs new mortality and populates snagTable", {
  cohorts <- data.table(
    pixelID = c(1L, 2L),
    year    = c(1L, 1L),
    species = "Pinus strobus",
    B       = c(10.0, 5.0)
  )
  sim <- testInit(
    "DeadWood_snagDecay",
    times  = list(start = 0, end = 1),
    params = list(snagDecay = list(
      snagTransMat = snagTransMat_test,
      snagFallProb = snagFallProb_test,
      species      = "Pinus strobus"
    )),
    objects = list(cohortData = cohorts)
  )
  set.seed(42)
  sim <- spades(sim, events = c("init", "annual"))
  expect_true(nrow(sim$snagTable) + nrow(sim$fallenSnags) == 2L)
})

test_that("snagDecay annual DC never decreases", {
  cohorts <- data.table(
    pixelID = 1:20,
    year    = rep(1L, 20),
    species = "Pinus strobus",
    B       = rep(5.0, 20)
  )
  sim <- testInit(
    "DeadWood_snagDecay",
    times  = list(start = 0, end = 10),
    params = list(snagDecay = list(
      snagTransMat = snagTransMat_test,
      snagFallProb = rep(0, 5),
      species      = "Pinus strobus"
    )),
    objects = list(cohortData = cohorts)
  )
  set.seed(1)
  sim <- spades(sim)
  expect_true(all(sim$snagTable$DC >= 1L & sim$snagTable$DC <= 5L))
})

test_that("snagDecay annual with 100% fall probability empties snagTable each year", {
  cohorts <- data.table(
    pixelID = 1L, year = 1L, species = "Pinus strobus", B = 5.0
  )
  sim <- testInit(
    "DeadWood_snagDecay",
    times  = list(start = 0, end = 1),
    params = list(snagDecay = list(
      snagTransMat = snagTransMat_test,
      snagFallProb = rep(1, 5),
      species      = "Pinus strobus"
    )),
    objects = list(cohortData = cohorts)
  )
  set.seed(5)
  sim <- spades(sim, events = c("init", "annual"))
  expect_equal(nrow(sim$snagTable), 0L)
  expect_equal(nrow(sim$fallenSnags), 1L)
})

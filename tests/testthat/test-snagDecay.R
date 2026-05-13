library(SpaDES.core)
library(data.table)
library(testthat)

snagTransMat_test <- list(
  "Pinus strobus" = matrix(
    c(0.48, 0.38, 0.05, 0.00, 0.00,
      0.00, 0.52, 0.34, 0.05, 0.00,
      0.00, 0.00, 0.56, 0.30, 0.04,
      0.00, 0.00, 0.00, 0.62, 0.26,
      0.00, 0.00, 0.00, 0.00, 0.70),
    nrow = 5, byrow = TRUE
  )
)
snagFallProb_test <- list(
  "Pinus strobus" = c(DC1 = 0.09, DC2 = 0.09, DC3 = 0.10, DC4 = 0.12, DC5 = 0.30)
)

snagTransMat_redpine <- matrix(c(
  0.109, 0.368, 0.471, 0.035, 0.017,
  0.000, 0.234, 0.434, 0.252, 0.081,
  0.000, 0.000, 0.704, 0.270, 0.027,
  0.000, 0.000, 0.000, 0.913, 0.087,
  0.000, 0.000, 0.000, 0.000, 1.000
), nrow = 5, byrow = TRUE)

snagTransMat_2sp <- list(
  "Pinus strobus"  = snagTransMat_test[["Pinus strobus"]],
  "Pinus resinosa" = snagTransMat_redpine
)
snagFallProb_2sp <- list(
  "Pinus strobus"  = c(DC1 = 0.09, DC2 = 0.09, DC3 = 0.10, DC4 = 0.12, DC5 = 0.30),
  "Pinus resinosa" = c(DC1 = 0.000, DC2 = 0.039, DC3 = 0.094, DC4 = 0.224, DC5 = 0.080)
)

emptyCohorts <- data.table(
  pixelID = integer(), year = integer(),
  species = character(), biomass = numeric()
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
    params = list(DeadWood_snagDecay = list(
      snagTransMat = snagTransMat_test,
      snagFallProb = snagFallProb_test,
      species      = c("Pinus strobus")
    )),
    objects = list(cohortData = emptyCohorts)
  )
  sim <- spades(sim, events = "init")
  expect_s3_class(sim$snagTable, "data.table")
  expect_equal(nrow(sim$snagTable), 0L)
  expect_named(sim$snagTable, c("pixelID", "species", "DC", "ageInDC", "initBiomass", "diameter_cm"))
})

test_that("snagDecay init creates empty fallenSnags", {
  sim <- testInit(
    "DeadWood_snagDecay",
    params = list(DeadWood_snagDecay = list(
      snagTransMat = snagTransMat_test,
      snagFallProb = snagFallProb_test,
      species      = c("Pinus strobus")
    )),
    objects = list(cohortData = emptyCohorts)
  )
  sim <- spades(sim, events = "init")
  expect_s3_class(sim$fallenSnags, "data.table")
  expect_equal(nrow(sim$fallenSnags), 0L)
  expect_named(sim$fallenSnags, c("pixelID", "species", "DC", "ageInDC", "initBiomass", "diameter_cm"))
})

test_that("snagDecay transition absorbs new mortality and populates snagTable", {
  cohorts <- data.table(
    pixelID = c(1L, 2L),
    year    = c(1L, 1L),
    species = "Pinus strobus",
    biomass = c(10.0, 5.0)
  )
  sim <- testInit(
    "DeadWood_snagDecay",
    times  = list(start = 0, end = 5),
    params = list(DeadWood_snagDecay = list(
      snagTransMat = snagTransMat_test,
      snagFallProb = snagFallProb_test,
      species      = c("Pinus strobus")
    )),
    objects = list(cohortData = cohorts)
  )
  set.seed(42)
  sim <- spades(sim, events = c("init", "transition"))
  expect_true(nrow(sim$snagTable) + nrow(sim$fallenSnags) == 2L)
})

test_that("snagDecay transition DC never decreases", {
  cohorts <- data.table(
    pixelID = 1:20,
    year    = rep(1L, 20),
    species = "Pinus strobus",
    biomass = rep(5.0, 20)
  )
  sim <- testInit(
    "DeadWood_snagDecay",
    times  = list(start = 0, end = 10),
    params = list(DeadWood_snagDecay = list(
      snagTransMat = snagTransMat_test,
      snagFallProb = list("Pinus strobus" = rep(0, 5)),
      species      = c("Pinus strobus")
    )),
    objects = list(cohortData = cohorts)
  )
  set.seed(1)
  sim <- spades(sim)
  expect_true(all(sim$snagTable$DC >= 1L & sim$snagTable$DC <= 5L))
})

test_that("snagDecay transition with 100% fall probability empties snagTable each step", {
  cohorts <- data.table(
    pixelID = 1L, year = 1L, species = "Pinus strobus", biomass = 5.0
  )
  sim <- testInit(
    "DeadWood_snagDecay",
    times  = list(start = 0, end = 5),
    params = list(DeadWood_snagDecay = list(
      snagTransMat = snagTransMat_test,
      snagFallProb = list("Pinus strobus" = rep(1, 5)),
      species      = c("Pinus strobus")
    )),
    objects = list(cohortData = cohorts)
  )
  set.seed(5)
  sim <- spades(sim, events = c("init", "transition"))
  expect_equal(nrow(sim$snagTable), 0L)
  expect_equal(nrow(sim$fallenSnags), 1L)
})

test_that("Init stops when snagTransMat is missing an entry for a listed species", {
  sim <- testInit(
    "DeadWood_snagDecay",
    params = list(DeadWood_snagDecay = list(
      snagTransMat = snagTransMat_test,           # only has Pinus strobus
      snagFallProb = snagFallProb_2sp,
      species      = c("Pinus strobus", "Pinus resinosa")
    )),
    objects = list(cohortData = emptyCohorts)
  )
  expect_error(spades(sim, events = "init"), regexp = "snagTransMat.*Pinus resinosa")
})

test_that("Init stops when snagFallProb is missing an entry for a listed species", {
  sim <- testInit(
    "DeadWood_snagDecay",
    params = list(DeadWood_snagDecay = list(
      snagTransMat = snagTransMat_2sp,
      snagFallProb = snagFallProb_test,           # only has Pinus strobus
      species      = c("Pinus strobus", "Pinus resinosa")
    )),
    objects = list(cohortData = emptyCohorts)
  )
  expect_error(spades(sim, events = "init"), regexp = "snagFallProb.*Pinus resinosa")
})

test_that("Transition absorbs mortality for both species from cohortData", {
  cohorts <- data.table(
    pixelID = c(1L, 2L, 3L),
    year    = c(1L, 1L, 1L),
    species = c("Pinus strobus", "Pinus resinosa", "Pinus strobus"),
    biomass = c(10.0, 8.0, 6.0)
  )
  sim <- testInit(
    "DeadWood_snagDecay",
    times  = list(start = 0, end = 5),
    params = list(DeadWood_snagDecay = list(
      snagTransMat = snagTransMat_2sp,
      snagFallProb = list(
        "Pinus strobus"  = rep(0, 5),
        "Pinus resinosa" = rep(0, 5)
      ),
      species = c("Pinus strobus", "Pinus resinosa")
    )),
    objects = list(cohortData = cohorts)
  )
  set.seed(1)
  sim <- spades(sim, events = c("init", "transition"))
  expect_equal(nrow(sim$snagTable[species == "Pinus strobus"]),  2L)
  expect_equal(nrow(sim$snagTable[species == "Pinus resinosa"]), 1L)
  expect_equal(nrow(sim$fallenSnags), 0L)
})

test_that("Species-specific fall probabilities are applied independently per species", {
  cohorts <- data.table(
    pixelID = 1:20,
    year    = rep(1L, 20),
    species = c(rep("Pinus strobus", 10), rep("Pinus resinosa", 10)),
    biomass = rep(5.0, 20)
  )
  sim <- testInit(
    "DeadWood_snagDecay",
    times  = list(start = 0, end = 5),
    params = list(DeadWood_snagDecay = list(
      snagTransMat = snagTransMat_2sp,
      snagFallProb = list(
        "Pinus strobus"  = rep(1, 5),
        "Pinus resinosa" = rep(0, 5)
      ),
      species = c("Pinus strobus", "Pinus resinosa")
    )),
    objects = list(cohortData = cohorts)
  )
  set.seed(42)
  sim <- spades(sim, events = c("init", "transition"))
  expect_equal(nrow(sim$fallenSnags[species == "Pinus strobus"]),  10L)
  expect_equal(nrow(sim$snagTable[species == "Pinus resinosa"]),   10L)
  expect_equal(nrow(sim$fallenSnags[species == "Pinus resinosa"]),  0L)
})

test_that("DC stays within [1,5] for both species across 50 years", {
  cohorts <- data.table(
    pixelID = 1:40,
    year    = rep(1L, 40),
    species = c(rep("Pinus strobus", 20), rep("Pinus resinosa", 20)),
    biomass = rep(5.0, 40)
  )
  sim <- testInit(
    "DeadWood_snagDecay",
    times  = list(start = 0, end = 50),
    params = list(DeadWood_snagDecay = list(
      snagTransMat = snagTransMat_2sp,
      snagFallProb = list(
        "Pinus strobus"  = rep(0, 5),
        "Pinus resinosa" = rep(0, 5)
      ),
      species = c("Pinus strobus", "Pinus resinosa")
    )),
    objects = list(cohortData = cohorts)
  )
  set.seed(7)
  sim <- spades(sim)
  expect_true(all(sim$snagTable$DC >= 1L & sim$snagTable$DC <= 5L))
})

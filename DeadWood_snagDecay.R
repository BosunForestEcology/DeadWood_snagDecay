# spatialExtent field omitted: removed from SpaDES.core API in version >= 3.0
defineModule(sim, list(
  name        = "DeadWood_snagDecay",
  description = "Advances standing dead White Pine trees through DC1-DC5 annually using
                 a Markov transition matrix, and stochastically transfers fallen snags
                 to sim$fallenSnags for consumption by DeadWood_DWDDecay.",
  keywords    = c("dead wood", "snag", "decay class", "Markov", "White Pine"),
  authors     = structure(list(list(given = "First", family = "Last",
                                    role = c("aut", "cre"),
                                    email = "email@example.com", comment = NULL)),
                           class = "person"),
  childModules = character(0),
  version     = list(DeadWood_snagDecay = "0.0.1"),
  timeframe   = as.POSIXlt(c(NA, NA)),
  timeunit    = "year",
  citation    = list(),
  documentation = list(),
  reqdPkgs    = list("data.table", "SpaDES.core (>= 3.0.0)"),
  parameters  = bindrows(
    defineParameter("snagTransMat", "matrix", matrix(0, 5, 5), NA, NA,
                    desc = "5x5 annual DC transition probability matrix for snags."),
    defineParameter("snagFallProb", "numeric", rep(0.1, 5), 0, 1,
                    desc = "Annual fall probability by DC (length 5)."),
    defineParameter("species", "character", "Pinus strobus", NA, NA,
                    desc = "Species to filter from cohortData.")
  ),
  inputObjects = bindrows(
    expectsInput("cohortData", "data.table",
                 desc = "Pixel-level cohort table with columns: pixelID, year, species, B (Mg/ha).")
  ),
  outputObjects = bindrows(
    createsOutput("snagTable", "data.table",
                  desc = "Current snag inventory: pixelID, species, DC, ageInDC, initBiomass."),
    createsOutput("fallenSnags", "data.table",
                  desc = "Snags that fell this timestep: same schema as snagTable.")
  )
))

doEvent.DeadWood_snagDecay <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, start(sim) + 1, "DeadWood_snagDecay", "annual", eventPriority = 1)
    },
    annual = {
      sim <- Annual(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "DeadWood_snagDecay", "annual", eventPriority = 1)
    },
    warning(paste("Undefined event type:", eventType, "in module snagDecay"))
  )
  return(invisible(sim))
}

Init <- function(sim) {
  if (all(P(sim)$snagTransMat == 0))
    stop("snagTransMat is the zero matrix — provide a real transition matrix in params.")
  if (length(P(sim)$snagFallProb) != 5L)
    stop("snagFallProb must have length 5 (one probability per decay class).")

  sim$snagTable <- data.table::data.table(
    pixelID     = integer(),
    species     = character(),
    DC          = integer(),
    ageInDC     = integer(),
    initBiomass = numeric()
  )
  sim$fallenSnags <- data.table::copy(sim$snagTable)
  return(invisible(sim))
}

Annual <- function(sim) {
  # Absorb new mortality for this year and this species
  newDead <- sim$cohortData[year == time(sim) & species == P(sim)$species]
  if (nrow(newDead) > 0) {
    sim$snagTable <- data.table::rbindlist(list(
      sim$snagTable,
      newDead[, .(pixelID, species, DC = 1L, ageInDC = 0L, initBiomass = B)]
    ))
  }

  if (nrow(sim$snagTable) == 0) {
    sim$fallenSnags <- data.table::copy(sim$snagTable)
    return(invisible(sim))
  }

  # Advance decay class via Markov transition
  oldDC <- sim$snagTable$DC
  sim$snagTable[, DC := applyTransition(DC, P(sim)$snagTransMat)]
  sim$snagTable[, ageInDC := data.table::fifelse(DC == oldDC, ageInDC + 1L, 0L)]

  # Stochastically simulate falls based on post-transition DC
  fallIdx <- sim$snagTable[, stats::rbinom(.N, 1L, P(sim)$snagFallProb[DC]) == 1L]
  sim$fallenSnags <- sim$snagTable[fallIdx]
  sim$snagTable   <- sim$snagTable[!fallIdx]

  return(invisible(sim))
}

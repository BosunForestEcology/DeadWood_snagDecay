# spatialExtent field omitted: removed from SpaDES.core API in version >= 3.0
defineModule(sim, list(
  name        = "DeadWood_snagDecay",
  description = "Advances standing dead trees (Pinus strobus, Pinus resinosa) through DC1-DC5
               every 5 years using species-specific Markov transition matrices, and stochastically
               transfers fallen snags to sim$fallenSnags for consumption by DeadWood_DWDDecay.",
  keywords    = c("dead wood", "snag", "decay class", "Markov", "White Pine", "Red Pine"),
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
    defineParameter("snagTransMat", "list",
                    list(
                      "Pinus strobus" = matrix(c(
                        0.109, 0.368, 0.471, 0.035, 0.017,  # from DC1
                        0.000, 0.348, 0.601, 0.051, 0.000,  # from DC2
                        0.000, 0.000, 0.704, 0.270, 0.027,  # from DC3
                        0.000, 0.000, 0.000, 0.913, 0.087,  # from DC4
                        0.000, 0.000, 0.000, 0.000, 1.000   # from DC5
                      ), nrow = 5, byrow = TRUE,
                      dimnames = list(paste0("from_DC", 1:5), paste0("to_DC", 1:5))),
                      "Pinus resinosa" = matrix(c(
                        0.109, 0.368, 0.471, 0.035, 0.017,  # from DC1
                        0.000, 0.234, 0.434, 0.252, 0.081,  # from DC2
                        0.000, 0.000, 0.704, 0.270, 0.027,  # from DC3
                        0.000, 0.000, 0.000, 0.913, 0.087,  # from DC4
                        0.000, 0.000, 0.000, 0.000, 1.000   # from DC5
                      ), nrow = 5, byrow = TRUE,
                      dimnames = list(paste0("from_DC", 1:5), paste0("to_DC", 1:5)))
                    ),
                    NA, NA,
                    desc = "Named list of 5x5 5-year conditional DC transition probability matrices,
                        one entry per species. Names must match entries in the species parameter.
                        Source: Vanderwel et al. 2006 Table 3."),
    defineParameter("snagFallProb", "list",
                    list(
                      "Pinus strobus"  = c(DC1 = 0.130, DC2 = 0.164, DC3 = 0.212, DC4 = 0.325, DC5 = 0.200),
                      "Pinus resinosa" = c(DC1 = 0.000, DC2 = 0.039, DC3 = 0.094, DC4 = 0.224, DC5 = 0.080)
                    ),
                    NA, NA,
                    desc = "Named list of 5-year snag fall probabilities by DC, one named numeric
                        vector (length 5) per species. Names must match entries in the species
                        parameter. Source: Vanderwel et al. 2006 Table 2."),
    defineParameter("species", "character",
                    c("Pinus strobus", "Pinus resinosa"),
                    NA, NA,
                    desc = "Character vector of species to process from cohortData. Each entry must
                        have a corresponding named entry in snagTransMat and snagFallProb."),
    defineParameter("defaultDiameter_cm", "numeric", 19.0, 7.5, NA,
                    desc = "Fallback diameter (cm) assigned to incoming snags when cohortData lacks a diameter_cm column. Minimum 7.5 cm.")
  ),
  inputObjects = bindrows(
    expectsInput("cohortData", "data.table",
                 desc = "Pixel-level cohort table with columns: pixelID, year, species, biomass (Mg/ha), diameter_cm (cm).")
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
      sim <- scheduleEvent(sim, start(sim) + 5, "DeadWood_snagDecay", "transition", eventPriority = 1)
    },
    transition = {
      sim <- Transition(sim)
      sim <- scheduleEvent(sim, time(sim) + 5, "DeadWood_snagDecay", "transition", eventPriority = 1)
    },
    warning(paste("Undefined event type:", eventType, "in module snagDecay"))
  )
  return(invisible(sim))
}

Init <- function(sim) {
  for (sp in P(sim)$species) {
    if (!sp %in% names(P(sim)$snagTransMat))
      stop(sprintf("snagTransMat has no entry for species '%s'. Add it to the parameter.", sp))
    if (!sp %in% names(P(sim)$snagFallProb))
      stop(sprintf("snagFallProb has no entry for species '%s'. Add it to the parameter.", sp))
    mat <- P(sim)$snagTransMat[[sp]]
    if (!is.matrix(mat) || !identical(dim(mat), c(5L, 5L)))
      stop(sprintf("snagTransMat[['%s']] must be a 5x5 matrix.", sp))
    if (any(rowSums(mat) < .Machine$double.eps))
      stop(sprintf("snagTransMat[['%s']] has an all-zero row — provide valid transition probabilities.", sp))
    fprob <- P(sim)$snagFallProb[[sp]]
    if (length(fprob) != 5L)
      stop(sprintf("snagFallProb[['%s']] must have length 5 (one probability per DC).", sp))
    if (any(fprob < 0) || any(fprob > 1))
      stop(sprintf("snagFallProb[['%s']] must have all values in [0, 1].", sp))
  }

  sim$snagTable <- data.table::data.table(
    pixelID     = integer(),
    species     = character(),
    DC          = integer(),
    ageInDC     = integer(),
    initBiomass = numeric(),
    diameter_cm = numeric()
  )
  sim$fallenSnags <- data.table::copy(sim$snagTable)
  return(invisible(sim))
}

Transition <- function(sim) {
  # Absorb mortality from the preceding 5-year interval for all tracked species
  newDead <- sim$cohortData[year > (time(sim) - 5) & year <= time(sim) &
                              species %in% P(sim)$species]
  if (nrow(newDead) > 0L) {
    if (!"diameter_cm" %in% names(newDead)) {
      warning("cohortData lacks diameter_cm — using defaultDiameter_cm (",
              P(sim)$defaultDiameter_cm, " cm) for all new snags.")
      newDead[, diameter_cm := P(sim)$defaultDiameter_cm]
    }
    if (anyNA(newDead$diameter_cm) || any(newDead$diameter_cm < 7.5, na.rm = TRUE))
      stop(sprintf(
        "cohortData contains %d snag(s) with diameter_cm < 7.5 cm. Minimum is 7.5 cm. Smallest: %.4g cm.",
        sum(newDead$diameter_cm < 7.5, na.rm = TRUE), min(newDead$diameter_cm, na.rm = TRUE)
      ))
    sim$snagTable <- data.table::rbindlist(list(
      sim$snagTable,
      newDead[, .(pixelID, species, DC = 1L, ageInDC = 0L, initBiomass = biomass, diameter_cm)]
    ))
  }

  if (nrow(sim$snagTable) == 0L) {
    sim$fallenSnags <- data.table::copy(sim$snagTable)
    return(invisible(sim))
  }

  if (anyNA(sim$snagTable$diameter_cm) || any(sim$snagTable$diameter_cm < 7.5))
    stop(sprintf(
      "snagTable contains piece(s) with diameter_cm < 7.5 cm or NA. Smallest: %.4g cm. Check cohortData$diameter_cm for 0 or missing values.",
      min(sim$snagTable$diameter_cm, na.rm = TRUE)
    ))

  oldDC   <- sim$snagTable$DC
  fallVec <- logical(nrow(sim$snagTable))

  for (sp in P(sim)$species) {
    idx <- which(sim$snagTable$species == sp)
    if (length(idx) == 0L) next

    spMat      <- P(sim)$snagTransMat[[sp]]
    spFallProb <- P(sim)$snagFallProb[[sp]]

    sim$snagTable[idx, DC := applyTransition(DC, spMat)]
    fallVec[idx] <- stats::rbinom(length(idx), 1L, spFallProb[oldDC[idx]]) == 1L
  }

  sim$snagTable[, ageInDC := data.table::fifelse(DC == oldDC, ageInDC + 5L, 0L)]

  sim$fallenSnags <- sim$snagTable[ fallVec, .(pixelID, species, DC, ageInDC, initBiomass, diameter_cm)]
  sim$snagTable   <- sim$snagTable[!fallVec, .(pixelID, species, DC, ageInDC, initBiomass, diameter_cm)]

  return(invisible(sim))
}

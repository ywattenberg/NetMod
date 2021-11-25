#' run a \code{siena07} algorithm until convergence
#' adapted from section 6.3.3 RSiena manual
#'
#' @param alg A control object, of class \code{sienaAlgorithm}
#' @param dat a siena data object as returned by \code{RSiena::sienaDataCreate}
#' @param eff a siena effects object as returned by \code{RSiena::getEffects}
#' @param cpu int number of cpus used for the estimation
#' @param useCluster logical, wheter to use a cluster of processes 
#'   (useful if multiple processors are available)
#' @param maxRuns int maximum number of times that \code{siena07} would run
#' @param max.value double, value used to check convergence of the parameters
#' @param batch desire interface.\code{FALSE} gives a gui, 
#'   \code{TRUE} print a small report to the console
#' @param returnDeps logical. Wheter to return the simulated networks in Phase 3
#' @param ... Arguments for the simulation function, see \code{RSiena::simstats0c}
#'
#' @return
#' @noRd
#'
#' @examples
siena07ToConvergence <- function(alg, dat, eff, cpu = 2,
                                 useCluster = TRUE, maxRuns = 10, max.value = 0.1,
                                 batch = FALSE, returnDeps = TRUE, ...) {
  numr <- 0
  # useCluster <- ifelse(cpu > 1,T,F)
  ans <- siena07(alg,
    data = dat, effects = eff,
    batch = batch, useCluster = useCluster, nbrNodes = cpu, # initC=T,
    returnDeps = returnDeps, ...
  ) # note that we return dependent variables from first run
  repeat {
    numr <- numr + 1 # count number of repeated runs
    maxt <- max(abs(ans$tconv[!ans$effects$fix[ans$effects$requested]]))
    tm <- ans$tconv.max
    # convergence indicator, excluding the fixed effects
    cat(numr, maxt, tm, "\n") # report how far we are
    if (maxt < max.value & tm < 0.2) {
      break
    } # success
    if (maxt > 3) {
      break
    } # divergence without much hope
    if (any(abs(ans$theta) > 20)) {
      break
    } # Block Elmer extra rule
    # of returning to good parameter values
    if (numr > maxRuns) {
      break
    } # now it has lasted too long
    ans <- siena07(alg,
      data = dat, effects = eff,
      batch = batch, useCluster = useCluster, nbrNodes = cpu, # initC=T,
      returnDeps = returnDeps, prevAns = ans, ...
    ) # later runs start
  }
  if (maxt > 0.10 & tm > 0.2) {
    cat("Warning: convergence inadequate.\n")
  }
  ans
}

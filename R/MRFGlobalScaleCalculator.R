#' Sets a global scale parameter for a GMRF or HSMRF model given a prior mean number of effective shifts.
#'
#' This function finds the global scale parameter value that produces the desired prior mean number of "effective" rate shifts.
#' Given a specified magnitude for an effective shift, shift.size, an effective shift occurs when two adjacent values are more than shift.size-fold apart from each other.
#' That is, an effective shift is the event that rate[i+1]/rate[i] > shift.size or rate[i+1]/rate[i] < 1/shift.size.
#'
#' Finding these values for a HSMRF model can take several seconds for large values of n.episodes because of the required numerical integration.
#'
#' @param n.episodes (numeric; no default) The number of episodes in the random field (the parameter vector will be this long).
#' @param model (character; no default) What model should the global scale parameter be set for? Options are "GMRF" and "HSMRF" for first-order models (also allowable: "GMRF1" and "HSMRF1") and "GMRF2" and HSMRF2" for second-order models.
#' @param prior.n.shifts (numeric; log(2)) The desired prior mean number of shifts.
#' @param shift.size (numeric; 2) The magnitude of change that defines an effective shift (measured as a fold-change).
#' @return The hyperprior.
#'
#' @references
#'
#' Magee et al. (2019) Locally adaptive Bayesian birth-death model successfully detects slow and rapid rate shifts.
#' doi: https://doi.org/10.1101/853960
#'
#' @examples
#' # Get global scale for a HSMRF model with 100 episodes.
#' gs <- setMRFGlobalScaleHyperpriorNShifts(100,"HSMRF")
#'
#' # Plot a draw from this HSMRF distribution

#' trajectory <- simulateMRF(n.episodes=100,model="HSMRF",global.scale.hyperprior=gs)
#' plot(1:100,rev(trajectory),type="l",xlab="time",ylab="speciation rate")
#'
#' @export

setMRFGlobalScaleHyperpriorNShifts <- function(n.episodes,model,prior.n.shifts=log(2),shift.size=2) {
  if ( model == "GMRF" ) {
    hyperprior <- setGMRFGlobalScaleExpectedNumberOfJumps(n.episodes,prior.n.shifts,shift.size)
  } else if ( model == "HSMRF" ) {
    hyperprior <- setHSMRFGlobalScaleExpectedNumberOfJumps(n.episodes,prior.n.shifts,shift.size)
  } else {
    stop("Unrecognized option for \"model\"")
  }
  return(hyperprior$hyperprior)
}



# Calculates global scale parameter for a Gaussian Markov random fielf from the prior mean number of "effective shifts" in the rate.
setHSMRFGlobalScaleExpectedNumberOfJumps <- function(n.episodes,prior.n.shifts=log(2),shift.size=2) {
  # We treat the change between each grid cell as a Bernoulli RV, so the collection of changes becomes binomial
  # From this we can calculate the expected number of cells where a shift occurs

  # recover()

  # Move to log-scale
  shift <- log(shift.size)

  # Probability of a shift for a value of zeta
  # We average the conditional p(shift | gamma) over p(gamma)
  quants <- seq(0.0001,0.9999,length.out=2000)

  # Transform so we can look up quantiles under regular cauchy distribution
  quants <- 1.0 - (1.0 - quants)/2.0
  probs <- 1/length(quants) # we're using quantiles, each gamma is equally likely

  # Function to optimize
  fn <- function(zeta) {
    # Grid of gammas
    gammas <- qcauchy(quants,0,zeta)
    # Number of expected shifts for each value of sigma
    num_expected_shifts <- sapply(gammas,function(x) {
      p_shift_one_cell_this_gamma <- pRightTailHorseshoeGrid(shift,x,grid.size=2000)/0.5
      return(p_shift_one_cell_this_gamma * (n.episodes-1))
    })
    # Average the per-sigma E(n_shifts) over p(sigma) to get overall expectation given zeta
    this_expected_num_shifts <- sum(probs * num_expected_shifts)
    return( (log(this_expected_num_shifts) - log(prior.n.shifts))^2 ) # Distance to target
  }

  # Find best value of zeta
  opts <- optimize(fn,c(0,1))
  zeta <- opts$minimum

  # Compute the prior on number of shifts for this zeta (to show user how well we approximated the target)
  gammas <- qcauchy(quants,0,zeta)
  num_expected_shifts <- sapply(gammas,function(x) {
    p_shift_one_cell_this_gamma <- pRightTailHorseshoeGrid(shift,x,grid.size=2000)/0.5
    return(p_shift_one_cell_this_gamma * (n.episodes-1))
  })

  # Estimate the error of our chosen global scale hyperprior
  computed_num_expected_shifts <- sum(probs * num_expected_shifts)
  return(list(hyperprior=zeta,E.n=computed_num_expected_shifts))
}

# Calculates global scale parameter for a Gaussian Markov random fielf from the prior mean number of "effective shifts" in the rate.
setGMRFGlobalScaleExpectedNumberOfJumps <- function(n.episodes,prior.n.shifts=log(2),shift.size=2) {
  # We treat the change between each grid cell as a Bernoulli RV, so the collection of changes becomes binomial
  # From this we can calculate the expected number of cells where a shift occurs

  # recover()

  # Move to log-scale
  shift <- log(shift.size)

  # Probability of a shift for a value of zeta
  # We average the conditional p(shift | sigma) over p(sigma)
  quants <- seq(0.0001,0.9999,length.out=2000)

  # Transform so we can look up quantiles under regular cauchy distribution
  quants <- 1.0 - (1.0 - quants)/2.0
  probs <- 1/length(quants) # we're using quantiles, each gamma is equally likely

  # Function to optimize
  fn <- function(zeta) {
    # Grid of sigmas
    sigmas <- qcauchy(quants,0,zeta)
    # Number of expected shifts for each value of sigma
    num_expected_shifts <- sapply(sigmas,function(x) {
      p_shift_one_cell_this_sigma <- pnorm(shift,0,x,lower.tail=FALSE)/0.5
      return(p_shift_one_cell_this_sigma * (n.episodes-1))
    })
    # Average the per-sigma E(n_shifts) over p(sigma) to get overall expectation given zeta
    this_expected_num_shifts <- sum(probs * num_expected_shifts)
    return( (log(this_expected_num_shifts) - log(prior.n.shifts))^2 ) # Distance to target
  }

  # Find best value of zeta
  opts <- optimize(fn,c(0,1))
  zeta <- opts$minimum

  # Compute the prior on number of shifts for this zeta (to show user how well we approximated the target)
  sigmas <- qcauchy(quants,0,zeta)
  num_expected_shifts <- sapply(sigmas,function(x) {
    p_shift_one_cell_this_sigma <- pnorm(shift,0,x,lower.tail=FALSE)/0.5
    return(p_shift_one_cell_this_sigma * (n.episodes-1))
  })

  # Estimate the error of our chosen global scale hyperprior
  computed_num_expected_shifts <- sum(probs * num_expected_shifts)
  return(list(hyperprior=zeta,E.n=computed_num_expected_shifts))
}

#' Simulates a single Markov random field trajectory.
#'
#' This function simulates a draw from a HSMRF or GMRF distribution given a user-specified global scale parameter.
#' The MRF can be taken to be on the log-scale (such as for a birth rate) or the real-scale.
#' The first value must be specified
#'
#' @param n.episodes (numeric; no default) The number of episodes in the random field (the parameter vector will be this long).
#' @param model (character; no default) What model should the global scale parameter be set for? Options are "GMRF" and "HSMRF".
#' @param global.scale.hyperprior (numeric; no default) The hyperprior on the global scale parameter.
#' @param initial.value (numeric; NULL) The first value in the MRF. If no value is specified, the field is assumed to start at 0 (if exponentiate=FALSE) or 1 (if exponentiate=TRUE).
#' @param exponentiate (logical; TRUE) If TRUE, the MRF model is taken to be on the log-scale and the values are returned on the real-scale (note this means that the specified initial value will be the log of the true initial value). If FALSE, the model is taken to be on the real scale.
#' @return A vector drawn from the specified MRF model on the specified (log- or real-) scale.
#'
#' @references
#'
#' Magee et al. (2019) Locally adaptive Bayesian birth-death model successfully detects slow and rapid rate shifts.
#' doi: https://doi.org/10.1101/853960
#'
#' Faulkner, James R., and Vladimir N. Minin. Locally adaptive smoothing with Markov random fields and shrinkage priors.
#' \emph{Bayesian analysis}, \bold{13 (1)}, 225.
#'
#' @examples
#' # Simulate a 100-episode HSMRF model for a speciation-rate through time
#' trajectory <- simulateMRF(n.episodes=100,model="HSMRF",global.scale.hyperprior=0.0021)
#' plot(1:100,rev(trajectory),type="l",xlab="time",ylab="speciation rate")
#'
#' @export

simulateMRF <- function(n.episodes,model,global.scale.hyperprior,initial.value=NULL,exponentiate=TRUE) {
  ndiffs <- n.episodes - 1
  if ( toupper(model) == "GMRF" ) {
    local_scales <- rep(1.0,ndiffs)
  } else if ( toupper(model) == "HSMRF" ) {
    local_scales <- abs(rcauchy(ndiffs,0.0,1.0))
  } else {
    stop("Unrecognized option for \"model\"")
  }

  global_scale <- abs(rcauchy(1,0.0,1.0))

  delta <- rnorm(ndiffs,0.0,local_scales * global_scale * global.scale.hyperprior)

  x <- numeric(n.episodes)
  if ( is.numeric(initial.value) ) {
    x[1] <- initial.value
  } else {
    x[1] <- 0.0
  }

  x[2:n.episodes] <- x[1] + cumsum(delta)

  if ( exponentiate == TRUE ) {
    x <- exp(x)
  }

  return(x)
}

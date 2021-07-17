#' Sets a global scale parameter for a GMRF or HSMRF model given a prior
#' mean number of effective shifts.
#'
#' This function finds the global scale parameter value that produces the
#' desired prior mean number of "effective" rate shifts. Given a specified
#' magnitude for an effective shift, shift_size, an effective shift occurs
#' when two adjacent values are more than shift_size-fold apart from each other.
#' That is, an effective shift is the event that rate[i+1]/rate[i] > shift_size
#' or rate[i+1]/rate[i] < 1/shift_size.
#'
#' Finding these values for a HSMRF model can take several seconds for large
#' values of n_episodes because of the required numerical integration.
#'
#' @param n_episodes (numeric; no default) The number of episodes in the random
#' field (the parameter vector will be this long).
#' @param model (character; no default) What model should the global scale
#' parameter be set for? Options are "GMRF" and "HSMRF" for first-order models
#' (also allowable: "GMRF1" and "HSMRF1") and "GMRF2" and HSMRF2" for
#' second-order models.
#' @param prior_n_shifts (numeric; log(2)) The desired prior mean number of
#' shifts.
#' @param shift_size (numeric; 2) The magnitude of change that defines an
#' effective shift (measured as a fold-change).
#' @return The hyperprior.
#'
#' @references
#'
#' Magee et al. (2019) Locally adaptive Bayesian birth-death model
#' successfully detects slow and rapid rate shifts.
#' doi: https://doi.org/10.1101/853960
#'
#' @examples
#' \donttest{
#' # Get global scale for a HSMRF model with 100 episodes.
#' gs <- setMRFGlobalScaleHyperpriorNShifts(100, "HSMRF")
#'
#' # Plot a draw from this HSMRF distribution
#'
#' trajectory <- simulateMRF(n_episodes = 100,
#'                           model = "HSMRF",
#'                           global_scale_hyperprior = gs)
#'
#' plot(1:100,
#'      rev(trajectory),
#'      type = "l",
#'      xlab = "time",
#'      ylab = "speciation rate")
#' }
#' @export


setMRFGlobalScaleHyperpriorNShifts <-
  function(n_episodes,
           model,
           prior_n_shifts = log(2),
           shift_size = 2) {
    if (model == "GMRF") {
      hyperprior <-
        .setGMRFGlobalScaleExpectedNumberOfJumps(n_episodes,
                                                 prior_n_shifts,
                                                 shift_size)
    } else if (model == "HSMRF") {
      hyperprior <-
        .setHSMRFGlobalScaleExpectedNumberOfJumps(n_episodes,
                                                  prior_n_shifts,
                                                  shift_size)
    } else {
      stop("Unrecognized option for \"model\"")
    }
    return(hyperprior$hyperprior)
  }

#' Simulates a single Markov random field trajectory.
#'
#' This function simulates a draw from a HSMRF or GMRF distribution given a
#' user-specified global scale parameter. The MRF can be taken to be on the
#' log-scale (such as for a birth rate) or the real-scale. The first value
#' must be specified
#'
#' @param n_episodes (numeric; no default) The number of episodes in the random
#' field (the parameter vector will be this long).
#' @param model (character; no default) What model should the global scale
#' parameter be set for? Options are "GMRF" and "HSMRF".
#' @param global_scale_hyperprior (numeric; no default) The hyperprior on the
#' global scale parameter.
#' @param initial_value (numeric; NULL) The first value in the MRF. If no value
#' is specified, the field is assumed to start at 0 (if exponentiate=FALSE) or
#' 1 (if exponentiate=TRUE).
#' @param exponentiate (logical; TRUE) If TRUE, the MRF model is taken to be on
#' the log-scale and the values are returned on the real-scale (note this means
#' that the specified initial value will be the log of the true initial value).
#' If FALSE, the model is taken to be on the real scale.
#' @return A vector drawn from the specified MRF model on the specified
#' (log- or real-) scale.
#'
#' @references
#'
#' Magee et al. (2020) Locally adaptive Bayesian birth-death model
#' successfully detects slow and rapid rate shifts.
#' \emph{PLoS Computational Biology}, \bold{16 (10)}: e1007999.
#'
#' Faulkner, James R., and Vladimir N. Minin. Locally adaptive smoothing with
#' Markov random fields and shrinkage priors.
#' \emph{Bayesian analysis}, \bold{13 (1)}, 225.
#'
#' @examples
#' \donttest{
#' # Simulate a 100-episode HSMRF model for a speciation-rate through time
#' trajectory <- simulateMRF(n_episodes = 100,
#'                           model = "HSMRF",
#'                           global_scale_hyperprior = 0.0021)
#' plot(1:100,
#'      rev(trajectory),
#'      type = "l",
#'      xlab = "time",
#'      ylab = "speciation rate")
#' }
#' @export

simulateMRF <-
  function(n_episodes,
           model,
           global_scale_hyperprior,
           initial_value = NULL,
           exponentiate = TRUE) {
    ndiffs <- n_episodes - 1
    if (toupper(model) == "GMRF") {
      local_scales <- rep(1.0, ndiffs)
    } else if (toupper(model) == "HSMRF") {
      local_scales <- abs(rcauchy(ndiffs, 0.0, 1.0))
    } else {
      stop("Unrecognized option for \"model\"")
    }

    global_scale <- abs(rcauchy(1, 0.0, 1.0))

    delta <-
      rnorm(ndiffs,
            0.0,
            local_scales * global_scale * global_scale_hyperprior)

    x <- numeric(n_episodes)
    if (is.numeric(initial_value)) {
      x[1] <- initial_value
    } else {
      x[1] <- 0.0
    }

    x[2:n_episodes] <- x[1] + cumsum(delta)

    if (exponentiate == TRUE) {
      x <- exp(x)
    }

    return(x)
  }

#' Bayes Factors in support of a shift in diversification rates over a given
#' time interval.
#'
#' This function computes the Bayes Factor in favor of a rate-shift between
#' time t1 and t2 (t1 < t2).
#' The default assumption (suitable to standard HSMRF and GMRF models) is
#' that the prior probability of a shift is 0.5.
#'
#' @param rate_trace (list; no default) The processed Rev output of the rate
#' of interest through time for computation (output of readTrace()).
#' @param time_trace (list; no default) The processed Rev output of the
#' change/interval times of the rate of interest through time for computation
#' (output of readTrace()).
#' @param rate_name (character; no default) The name of the parameter
#' (e.g. "speciation") for which Bayes Factor is to be calculated.
#' @param time_name (character; no default) The name of the interval times
#' (e.g. "interval_times) for the rate change times.
#' @param t1 (numeric; no default) Support will be assessed for a shift between
#' time t1 and time t2 (t1 < t2).
#' @param t2 (numeric; no default) Support will be assessed for a shift between
#' time t1 and time t2 (t1 < t2).
#' @param decrease (logical; default TRUE) Should support be assessed for a
#' decrease in the parameter (if TRUE) or an increase (if FALSE) between
#' t1 and t2?
#' @param prior_prob (numeric; 0.5) The prior probability of a shift over this
#' interval (default of 0.5 applies to standard HSMRF- and GMRF-based models).
#' @param return_2lnBF (logical; TRUE) Should the 2ln(BF) be returned
#' (if TRUE) or simply the BF (if FALSE)?
#'
#' @return The Bayes Factor.
#'
#' @references
#'
#' Kass and Raftery (1995) Bayes Factors.
#' \emph{JASA}, \bold{90 (430)}, 773-795.
#'
#' @examples
#'
#' \donttest{
#' #' # download the example datasets to working directory
#' url_times <-
#'    "https://revbayes.github.io/tutorials/intro/data/primates_EBD_speciation_times.log"
#' dest_path_times <- "primates_EBD_speciation_times.log"
#' download.file(url_times, dest_path_times)
#'
#' url_rates <-
#'    "https://revbayes.github.io/tutorials/intro/data/primates_EBD_speciation_rates.log"
#' dest_path_rates <- "primates_EBD_speciation_rates.log"
#' download.file(url_rates, dest_path_rates)
#'
#' # to run on your own data, change this to the path to your data file
#' speciation_time_file <- dest_path_times
#' speciation_rate_file <- dest_path_rates
#'
#' speciation_times <- readTrace(speciation_time_file, burnin = 0.25)
#' speciation_rate <- readTrace(speciation_rate_file, burnin = 0.25)
#'
#' calculateShiftBayesFactor(speciation_rate,
#'                           speciation_times,
#'                           "speciation",
#'                           "interval_times",
#'                           0.0,40.0,
#'                           decrease=FALSE)
#'
#' # remove file
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_times, dest_path_rates)
#' }
#'
#' @export


calculateShiftBayesFactor <-
  function(rate_trace,
           time_trace,
           rate_name,
           time_name,
           t1,
           t2,
           prior_prob = 0.5,
           decrease = TRUE,
           return_2lnBF = TRUE) {
    # Make sure times are in correct order
    times <- sort(unlist(c(t1, t2)))
    t1 <- times[1]
    t2 <- times[2]

    # Condense log lists into data frames
    rate_log <- do.call(rbind, rate_trace)
    time_log <- do.call(rbind, time_trace)

    # Find parameter and times, remove rest of trace
    is_rate <- grepl(paste0(rate_name, "["), names(rate_log), fixed = TRUE)
    is_time <- grepl(paste0(time_name, "["), names(time_log), fixed = TRUE)

    if (sum(is_rate) < 2) {
      stop(
        paste0(
          "Cannot find diversification rate parameter named \"",
          rate_name,
          "\" in rate_trace"
        )
      )
    }

    if (sum(is_time) < 1) {
      stop(
        paste0(
          "Cannot find rate change times parameter named \"",
          time_name,
          "\" in time_trace"
        )
      )
    }

    rate_log <- rate_log[, is_rate]
    time_log <- time_log[, is_time]

    # Add times[0] = 0 to time_log
    time_log$`interval_times[0]` <- rep(0.0, dim(time_log)[1])

    if (dim(time_log)[1] != dim(rate_log)[1]) {
      stop("Different number of samples for parameter and parameter times.")
    }

    if (dim(time_log)[2] != dim(rate_log)[2]) {
      stop("Number of interval times does not match dimension of parameter.")
    }

    # Check if the times are constant
    times_are_constant <- FALSE
    if (all(apply(time_log, 2, function(x) {
      max(x) - min(x) < .Machine$double.eps
    }))) {
      times_are_constant <- TRUE
    }

    # Vectors to store posterior distribution of rates at t1 and t2
    niter <- dim(time_log)[1]
    rate1 <- numeric(niter)
    rate2 <- numeric(niter)

    if (times_are_constant) {
      sorted_times <- sort(unlist(time_log[1, ]))

      # Avoid issues finding the maximum time
      if (t2 == max(sorted_times)) {
        t2 <- t2 * 0.9999
      }

      # Indices for rates at t1 and t2
      index1 <- findInterval(t1, sorted_times)
      index2 <- findInterval(t2, sorted_times)

      rate1 <- rate_log[, index1]
      rate2 <- rate_log[, index2]
    } else {
      for (i in 1:niter) {
        t1 <- times[1]
        t2 <- times[2]

        sorted_times <- sort(unlist(time_log[i, ]))

        # Avoid issues finding the maximum time
        if (t2 == max(sorted_times)) {
          t2 <- t2 * 0.9999
        }

        # Indices for rates at t1 and t2
        index1 <- findInterval(t1, sorted_times)
        index2 <- findInterval(t2, sorted_times)

        rate1[i] <- rate_log[i, index1]
        rate2[i] <- rate_log[i, index2]
      }
    }

    # Posterior probability of shift
    if (all(rate1 == rate2)) {
      stop("t1 and t2 are in the same time interval")
    }

    shift_prob <- sum(rate1 < rate2) / length(rate1)

    # Avoid infinite Bayes Factors
    if (shift_prob == 0.0) {
      shift_prob <- 1 / (length(rate1) + 1)
    } else if (shift_prob == 1.0) {
      shift_prob <- length(rate1) / (length(rate1) + 1)
    }

    if (decrease == FALSE) {
      shift_prob <- 1 - shift_prob
    }

    # Prior odds are 1 under a MRF model
    posterior_odds <- shift_prob / (1 - shift_prob)
    prior_odds <- prior_prob / (1 - prior_prob)
    BF <- posterior_odds / prior_odds

    if (return_2lnBF == TRUE) {
      BF <- 2 * log(BF)
    }

    return(BF)
  }

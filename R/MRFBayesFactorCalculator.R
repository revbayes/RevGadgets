#' Bayes Factors in support of a shift in diversification rates over a given time interval.
#'
#' This function computes the Bayes Factor in favor of a rate-shift between time t1 and t2 (t1 < t2).
#' The default assumption (suitable to standard HSMRF and GMRF models) is that the prior probability of a shift is 0.5.
#'
#' @param output (list; no default) The processed output for computation (output of processDivRates()).
#' @param parameter (character; no default) The name of the parameter (list element in output) for which Bayes Factor is to be calculated.
#' @param parameter.times (character; no default) The name of the interval times (list element in output) for the parameter.
#' @param t1 (numeric; no default) Support will be assesed for a shift between time t1 and time t2 (t1 < t2).
#' @param t2 (numeric; no default) Support will be assesed for a shift between time t1 and time t2 (t1 < t2).
#' @param decrease (logical; default TRUE) Should support be assessed for a decrease in the parameter (if TRUE) or an increase (if FALSE) between t1 and t2?
#' @param prior.prob (numeric; 0.5) The prior probability of a shift over this interval (default of 0.5 applies to standard HSMRF- and GMRF-based models).
#' @param lnBF (logical; TRUE) Should the 2ln(BF) be returned (if TRUE) or simply the BF (if FALSE)?
#'
#' @return The Bayes Factor.
#'
#' @references
#'
#' Kass and Raftery (1995) Bayes Factors.
#' \emph{JASA}, \bold{90 (430)}, 773-795.
#'
#' @examples
#'\dontrun{
#'
#' speciation_time_file <- system.file("extdata",
#'     "epi_bd/primates_EBD_speciation_times.p", package="RevGadgets")
#' speciation_rate_file <- system.file("extdata",
#'     "epi_bd/primates_EBD_speciation_rates.p", package="RevGadgets")
#' extinction_time_file <- system.file("extdata",
#'     "epi_bd/primates_EBD_extinction_times.p", package="RevGadgets")
#' extinction_rate_file <- system.file("extdata",
#'     "epi_bd/primates_EBD_extinction_rates.p", package="RevGadgets")
#'
#' primates <- processDivRates(speciation_time_log = speciation_time_file,
#'                             speciation_rate_log = speciation_rate_file,
#'                             extinction_time_log = extinction_time_file,
#'                             extinction_rate_log = extinction_rate_file,
#'                             burnin = 0.25)
#' calculateShiftBayesFactor(primates,"speciation rate","speciation time",0.0,40.0,decrease=FALSE)
#'}
#'
#' @export

calculateShiftBayesFactor <- function(output,parameter,parameter.times,t1,t2,prior.prob=0.5,decrease=TRUE,return.2lnBF=TRUE) {

  # Make sure times are in correct order
  times <- sort(c(t1,t2))
  t1 <- times[1]
  t2 <- times[2]

  # Find parameter and times
  parameter_index <- which(names(output) == parameter)
  if ( all(names(output) != parameter) ) {
    stop(paste0("Cannot find parameter named \"",parameter,"\" in output"))
  }

  times_index <- which(names(output) == parameter.times)
  if ( all(names(output) != parameter) ) {
    stop(paste0("Cannot find parameter named \"",parameter,"\" in output"))
  }

  # Get only what we need out of times and rates
  output[[parameter_index]] <- output[[parameter_index]][,!grepl("iteration",tolower(names(output[[parameter_index]])))]
  output[[parameter_index]] <- output[[parameter_index]][,!grepl("posterior",tolower(names(output[[parameter_index]])))]
  output[[parameter_index]] <- output[[parameter_index]][,!grepl("likelihood",tolower(names(output[[parameter_index]])))]
  output[[parameter_index]] <- output[[parameter_index]][,!grepl("prior",tolower(names(output[[parameter_index]])))]

  output[[times_index]] <- output[[times_index]][,!grepl("iteration",tolower(names(output[[times_index]])))]
  output[[times_index]] <- output[[times_index]][,!grepl("posterior",tolower(names(output[[times_index]])))]
  output[[times_index]] <- output[[times_index]][,!grepl("likelihood",tolower(names(output[[times_index]])))]
  output[[times_index]] <- output[[times_index]][,!grepl("prior",tolower(names(output[[times_index]])))]

  if ( dim(output[[times_index]])[1] != dim(output[[parameter_index]])[1] ) {
    stop("Different number of samples for parameter and parameter times.")
  }

  if ( dim(output[[times_index]])[2] != dim(output[[parameter_index]])[2] ) {
    stop("Number interval times does not match dimension of parameter.")
  }

  # Check if the times are constant
  times_are_constant <- FALSE
  if ( all(apply(output[[times_index]],2,function(x){max(x) - min(x) < .Machine$double.eps})) ) {
    times_are_constant <- TRUE
  }

  # Vectors to store posterior distribution of rates at t1 and t2
  niter <- dim(output[[times_index]])[1]
  rate1 <- numeric(niter)
  rate2 <- numeric(niter)

  if ( times_are_constant ) {
    sorted_times <- sort(output[[times_index]][1,])

    # Avoid issues finding the maximum time
    if ( t2 == max(sorted_times) ) {
      t2 <- t2 * 0.9999
    }

    # Indices for rates at t1 and t2
    index1 <- findInterval(t1,sorted_times)
    index2 <- findInterval(t2,sorted_times)

    rate1 <- output[[parameter_index]][,index1]
    rate2 <- output[[parameter_index]][,index2]
  } else {
    for (i in 1:niter) {
      t1 <- times[1]
      t2 <- times[2]

      sorted_times <- sort(output[[times_index]][i,])

      # Avoid issues finding the maximum time
      if ( t2 == max(sorted_times) ) {
        t2 <- t2 * 0.9999
      }

      # Indices for rates at t1 and t2
      index1 <- findInterval(t1,sorted_times)
      index2 <- findInterval(t2,sorted_times)

      rate1[i] <- output[[parameter_index]][i,index1]
      rate2[i] <- output[[parameter_index]][i,index2]
    }
  }

  # Posterior probability of shift
  if ( all(rate1 == rate2) ) {
    stop("t1 and t2 are in the same time interval")
  }

  shift_prob <- sum(rate1 < rate2)/length(rate1)

  # Avoid infinite Bayes Factors
  if ( shift_prob == 0.0 ) {
    shift_prob <- 1 / (length(rate1) + 1)
  } else if ( shift_prob == 1.0 ) {
    shift_prob <- length(rate1) / (length(rate1) + 1)
  }

  if ( decrease == FALSE ) {
    shift_prob <- 1 - shift_prob
  }

  # Prior odds are 1 under a MRF model
  posterior_odds <- shift_prob / (1 - shift_prob)
  prior_odds <- prior.prob / (1 - prior.prob)
  BF <- posterior_odds / prior_odds

  if ( return.2lnBF == TRUE ) {
    BF <- 2 * log(BF)
  }

  return(BF)
}


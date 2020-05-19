#' Process Diversification Rates
#'
#'
#' Processing the output of a episodic diversification rate analysis with mass-extinction events.
#'
#' For processing the output of an episodic diversification rate analysis. processDivRates()
#' assumes that the epochs are fixed rather than inferred. Additionally, it assumes
#' that times correspond to rates such that the first rate parameter (i.e. speciation[1])
#' corresponds to the present. Conversely, the first time parameter
#' (i.e. interval_times[1]) corresponds to the first time interval after the present,
#' moving backwards in time. processDivRates() relies on readTrace and produces a list
#' object that can be read by plotDivRates() to vizualize the results. For now,
#' only one log file per parameter type is accepted (i.e. log files from multiple runs
#' must be combined before reading into the function).
#'
#'@param speciation_time_log (vector of character strings or single character string; "") Path to speciation times log file(s)
#'@param speciation_rate_log (vector of character strings or single character string; "") Path to speciation rates log file(s)
#'@param extinction_time_log (vector of character strings or single character string; "") Path to extinction times log file(s)
#'@param extinction_rate_log (vector of character strings or single character string; "") Path to extinction rates log file(s)
#'@param fossilization_time_log (vector of character strings or single character string; "") Path to fossilization times log file(s)
#'@param fossilization_rate_log (vector of character strings or single character string; "") Path to fossilization rates log file(s)

#'@param burnin (single numeric value; default = 0) Fraction of generations to
#' discard (if value provided is between 0 and 1) or number of generations (if
#' value provided is greater than 1). Passed to readTrace().
#'
#'@return List object with processed rate and time parameters.
#'
#'@examples
#'
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
#'}
#'
#' @export

processDivRates <- function(speciation_time_log = "",
                            speciation_rate_log = "",
                            extinction_time_log = "",
                            extinction_rate_log = "",
                            fossilization_time_log = "",
                            fossilization_rate_log = "",
                            burnin = 0.25) {

  # enforce argument matching
  if (is.character(speciation_time_log) == FALSE) stop("speciation_time_log must be a character string or vector of strings")
  if (is.character(speciation_rate_log) == FALSE) stop("speciation_rate_log must be a character string or vector of strings")
  if (is.character(extinction_time_log) == FALSE) stop("extinction_time_log must be a character string or vector of strings")
  if (is.character(extinction_rate_log) == FALSE) stop("extinction_rate_log must be a character string or vector of strings")
  if (is.character(fossilization_time_log) == FALSE) stop("fossilization_time_log must be a character string vector of strings")
  if (is.character(fossilization_rate_log) == FALSE) stop("fossilization_rate_log must be a character string vector of strings")

    # check that both fozzilization files are either provided or not provided
  if (fossilization_time_log != "" &
      fossilization_rate_log == "") stop("both please provide both fossilization rates and times, or neither")
  if (fossilization_time_log == "" &
      fossilization_rate_log != "") stop("both please provide both fossilization rates and times, or neither")

  # check if speciation times log file(s) exist
  do_speciation_time_log_exist <- file.exists(speciation_time_log)
  if ( any(do_speciation_time_log_exist == FALSE) == TRUE ) {
    # print out paths to files that don't exist
    cat( "Some speciation_time_log files do not exist:",
         paste0("\t",speciation_time_log[do_speciation_time_log_exist == FALSE]), sep="\n")
    stop()
  }

  # check if speciation rates log file(s) exist
  do_speciation_rate_log_exist <- file.exists(speciation_rate_log)
  if ( any(do_speciation_rate_log_exist == FALSE) == TRUE ) {
    # print out paths to files that don't exist
    cat( "Some speciation_rate_log files do not exist:",
         paste0("\t",speciation_rate_log[do_speciation_rate_log_exist == FALSE]), sep="\n")
    stop()
  }

  # check if extinction times log file(s) exist
  do_extinction_time_log_exist <- file.exists(extinction_time_log)
  if ( any(do_extinction_time_log_exist == FALSE) == TRUE ) {
    # print out paths to files that don't exist
    cat( "Some extinction_time_log files do not exist:",
         paste0("\t",extinction_time_log[do_extinction_time_log_exist == FALSE]), sep="\n")
    stop()
  }

  # check if extinction rates log file(s) exist
  do_extinction_rate_log_exist <- file.exists(extinction_rate_log)
  if ( any(do_extinction_rate_log_exist == FALSE) == TRUE ) {
    # print out paths to files that don't exist
    cat( "Some extinction_rate_log files do not exist:",
         paste0("\t",extinction_rate_log[do_extinction_rate_log_exist == FALSE]), sep="\n")
    stop()
  }

  # check if fossilization times log file(s) exist if provided
  if (fossilization_time_log != "") {
   do_fossilization_time_log_exist <- file.exists(fossilization_time_log)
   if ( any(do_fossilization_time_log_exist == FALSE) == TRUE ) {
     # print out paths to files that don't exist
     cat( "Some fossilization_time_log files do not exist:",
          paste0("\t",fossilization_time_log[do_fossilization_time_log_exist == FALSE]), sep="\n")
     stop()
   }
  }

  # check if fossilization rates log file(s) exist if provided
  if (fossilization_rate_log != "") {
   do_fossilization_rate_log_exist <- file.exists(fossilization_rate_log)
   if ( any(do_fossilization_rate_log_exist == FALSE) == TRUE ) {
     # print out paths to files that don't exist
     cat( "Some fossilization_rate_log files do not exist:",
          paste0("\t",fossilization_rate_log[do_fossilization_rate_log_exist == FALSE]), sep="\n")
     stop()
   }
  }

  # read in log files as lists of data.frames with readTrace()
  speciation_time <- readTrace(path = speciation_time_log,
                                burnin = burnin)
  speciation_rate <- readTrace(path = speciation_rate_log,
                                burnin = burnin)
  extinction_time <- readTrace(path = extinction_time_log,
                                burnin = burnin)
  extinction_rate <- readTrace(path = extinction_rate_log,
                                burnin = burnin)
  if (fossilization_time_log != "") {
    fossilization_time <- readTrace(path = fossilization_time_log,
                                 burnin = burnin)
    fossilization_rate <- readTrace(path = fossilization_rate_log,
                                 burnin = burnin)
  } else {
    fossilization_time <- NULL
    fossilization_rate <- NULL}

  # check if all parameter types have the same number of log files
  if (is.null(fossilization_time)) {
  trace_lengths_same <- identical(length(speciation_time),
                                  length(speciation_rate),
                                  length(extinction_time),
                                  length(extinction_rate))
  } else {
    trace_lengths_same <- identical(length(speciation_time),
                                    length(speciation_rate),
                                    length(extinction_time),
                                    length(extinction_rate),
                                    length(fossilization_time),
                                    length(fossilization_rate))
  }

  if (trace_lengths_same == FALSE) {
    stop("You must provide the same number of log files for each parameter type.")}
  else if (trace_lengths_same == TRUE) {
    if (length(speciation_time) == 0) {
      stop("You must provide at least one log file per parameter type.")
    } else if (length(speciation_time) > 1) {
      stop("Currently, only one log file per parameter type is supported.")
    } else if (length(speciation_time) == 1) {

      # convert single item lists to data frames
      speciation_time <- speciation_time[[1]]
      speciation_rate <- speciation_rate[[1]]
      extinction_time <- extinction_time[[1]]
      extinction_rate <- extinction_rate[[1]]

      if (!is.null(fossilization_time)) {
        fossilization_time <- fossilization_time[[1]]
        fossilization_rate <- fossilization_rate[[1]]
      }

      # add in dummy distribution of 0 for time in the present
      speciation_time$`interval_times[0]` <- rep(0, nrow(speciation_time))
      extinction_time$`interval_times[0]` <- rep(0, nrow(extinction_time))
      if (!is.null(fossilization_time)) {fossilization_time$`interval_times[0]` <- rep(0, nrow(fossilization_time))}

      # Calculate the net-diversification and relative-extinction rates
      net_diversification_rate <- as.matrix(speciation_rate[,grepl("speciation", colnames(speciation_rate))]) -
                                   as.matrix(extinction_rate[,grepl("extinction", colnames(extinction_rate))])

      colnames(net_diversification_rate) <- paste(rep("net_div", times = ncol(net_diversification_rate)),
                                                   rep("[", times = ncol(net_diversification_rate)),
                                                   1:ncol(net_diversification_rate),
                                                   rep("]", times = ncol(net_diversification_rate)), sep = "")

      relative_extinction_rate <- as.matrix(extinction_rate[,grepl("extinction", colnames(extinction_rate))]) /
                                   as.matrix(speciation_rate[,grepl("speciation", colnames(speciation_rate))])

      colnames(relative_extinction_rate) <- paste(rep("rel_ext", times = ncol(relative_extinction_rate)),
                                                   rep("[", times = ncol(relative_extinction_rate)),
                                                   1:ncol(relative_extinction_rate),
                                                   rep("]", times = ncol(relative_extinction_rate)), sep = "")

        # return a list of processed data frames
        res <- list("speciation rate" = speciation_rate,
                    "extinction rate" = extinction_rate,
                    "net-diversification rate" = net_diversification_rate,
                    "relative-extinction rate" = relative_extinction_rate,
                    "fossilization rate" = fossilization_rate,
                    "speciation time" = speciation_time,
                    "extinction time" = extinction_time,
                    "fossilization time" = fossilization_time)
        return(res)
    }
  }
}

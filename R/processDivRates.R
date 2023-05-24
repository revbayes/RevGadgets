#' Process Diversification Rates
#'
#'
#' Processing the output of a episodic diversification rate analysis with
#' mass-extinction events.
#'
#' For processing the output of an episodic diversification rate analysis.
#' processDivRates() assumes that the epochs are fixed rather than inferred.
#' Additionally, it assumes that times correspond to rates such that the first
#' rate parameter (i.e. speciation[1]) corresponds to the present. Conversely,
#' the first time parameter (i.e. interval_times[1]) corresponds to the first
#' time interval after the present, moving backwards in time. processDivRates()
#' relies on readTrace and produces a list object that can be read by
#' plotDivRates() to visualize the results. For now, only one log file per
#' parameter type is accepted (i.e. log files from multiple runs must be
#' combined before reading into the function).
#'
#' @param speciation_time_log (vector of character strings or
#' single character string; "") Path to speciation times log file(s)
#' @param speciation_rate_log (vector of character strings or
#' single character string; "") Path to speciation rates log file(s)
#' @param extinction_time_log (vector of character strings or
#' single character string; "") Path to extinction times log file(s)
#' @param extinction_rate_log (vector of character strings or
#' single character string; "") Path to extinction rates log file(s)
#' @param fossilization_time_log (vector of character strings or
#' single character string; "") Path to fossilization times log file(s)
#' @param fossilization_rate_log (vector of character strings or
#' single character string; "") Path to fossilization rates log file(s)
#' @param burnin (single numeric value; default = 0) Fraction of generations to
#'  discard (if value provided is between 0 and 1) or number of generations (if
#'  value provided is greater than 1). Passed to readTrace().
#' @param probs (numeric vector; c(0.025, 0.975)) a vector of length two
#' containing the upper and lower bounds for the confidence intervals.
#' @param summary typically "mean" or "median"; the metric to summarize the
#' posterior distribution. Defaults to "Mean"
#' @return List object with processed rate and time parameters.
#'
#' @examples
#'
#' \donttest{
#' # download the example datasets to working directory
#'
#' url_ex_times <-
#'    "https://revbayes.github.io/tutorials/intro/data/primates_EBD_extinction_times.log"
#' dest_path_ex_times <- "primates_EBD_extinction_times.log"
#' download.file(url_ex_times, dest_path_ex_times)
#'
#' url_ex_rates <-
#'    "https://revbayes.github.io/tutorials/intro/data/primates_EBD_extinction_rates.log"
#' dest_path_ex_rates <- "primates_EBD_extinction_rates.log"
#' download.file(url_ex_rates, dest_path_ex_rates)
#'
#' url_sp_times <-
#'    "https://revbayes.github.io/tutorials/intro/data/primates_EBD_speciation_times.log"
#' dest_path_sp_times <- "primates_EBD_speciation_times.log"
#' download.file(url_sp_times, dest_path_sp_times)
#'
#' url_sp_rates <-
#'    "https://revbayes.github.io/tutorials/intro/data/primates_EBD_speciation_rates.log"
#' dest_path_sp_rates <- "primates_EBD_speciation_rates.log"
#' download.file(url_sp_rates, dest_path_sp_rates)
#'
#' # to run on your own data, change this to the path to your data file
#' speciation_time_file <- dest_path_sp_times
#' speciation_rate_file <- dest_path_sp_rates
#' extinction_time_file <- dest_path_ex_times
#' extinction_rate_file <- dest_path_ex_rates
#'
#' rates <- processDivRates(speciation_time_log = speciation_time_file,
#'                          speciation_rate_log = speciation_rate_file,
#'                          extinction_time_log = extinction_time_file,
#'                          extinction_rate_log = extinction_rate_file,
#'                          burnin = 0.25)
#'
#' # remove files
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_sp_times, dest_path_ex_times,
#'             dest_path_sp_rates, dest_path_ex_rates)
#' }
#'
#' @export

processDivRates <- function(speciation_time_log = "",
                            speciation_rate_log = "",
                            extinction_time_log = "",
                            extinction_rate_log = "",
                            fossilization_time_log = "",
                            fossilization_rate_log = "",
                            burnin = 0.25,
                            probs = c(0.025, 0.975),
                            summary = "mean") {
  # enforce argument matching
  if (is.character(speciation_time_log) == FALSE)
    stop("speciation_time_log must be a character string or vector of strings")
  if (is.character(speciation_rate_log) == FALSE)
    stop("speciation_rate_log must be a character string or vector of strings")
  if (is.character(extinction_time_log) == FALSE)
    stop("extinction_time_log must be a character string or vector of strings")
  if (is.character(extinction_rate_log) == FALSE)
    stop("extinction_rate_log must be a character string or vector of strings")
  if (is.character(fossilization_time_log) == FALSE)
    stop("fossilization_time_log must be a character string vector of strings")
  if (is.character(fossilization_rate_log) == FALSE)
    stop("fossilization_rate_log must be a character string vector of strings")

  # check that both fozzilization files are either provided or not provided
  if (fossilization_time_log != "" &
      fossilization_rate_log == "")
    stop("both please provide both fossilization rates and times, or neither")
  if (fossilization_time_log == "" &
      fossilization_rate_log != "")
    stop("both please provide both fossilization rates and times, or neither")

  # check if speciation times log file(s) exist
  do_speciation_time_log_exist <- file.exists(speciation_time_log)
  if (any(do_speciation_time_log_exist == FALSE) == TRUE) {
    # print out paths to files that don't exist
    stop(
      paste0(
      "Some speciation_time_log files do not exist:",
      paste0("\t", speciation_time_log[do_speciation_time_log_exist == FALSE]),
      sep = "\n"
    )
    )
  }

  # check if speciation rates log file(s) exist
  do_speciation_rate_log_exist <- file.exists(speciation_rate_log)
  if (any(do_speciation_rate_log_exist == FALSE) == TRUE) {
    # print out paths to files that don't exist
    stop(
      paste0(
      "Some speciation_rate_log files do not exist:",
      paste0("\t", speciation_rate_log[do_speciation_rate_log_exist == FALSE]),
      sep = "\n"
    )
    )
  }

  # check if extinction times log file(s) exist
  do_extinction_time_log_exist <- file.exists(extinction_time_log)
  if (any(do_extinction_time_log_exist == FALSE) == TRUE) {
    # print out paths to files that don't exist
    stop(
      paste0(
      "Some extinction_time_log files do not exist:",
      paste0("\t", extinction_time_log[do_extinction_time_log_exist == FALSE]),
      sep = "\n"
    )
    )
  }

  # check if extinction rates log file(s) exist
  do_extinction_rate_log_exist <- file.exists(extinction_rate_log)
  if (any(do_extinction_rate_log_exist == FALSE) == TRUE) {
    # print out paths to files that don't exist
    stop(
      paste0(
      "Some extinction_rate_log files do not exist:",
      paste0("\t", extinction_rate_log[do_extinction_rate_log_exist == FALSE]),
      sep = "\n"
    )
    )
  }

  # check if fossilization times log file(s) exist if provided
  if (fossilization_time_log != "") {
    do_fossilization_time_log_exist <-
      file.exists(fossilization_time_log)
    if (any(do_fossilization_time_log_exist == FALSE) == TRUE) {
      # print out paths to files that don't exist
      stop(
        paste0(
        "Some fossilization_time_log files do not exist:",
        paste0("\t",
               fossilization_time_log[do_fossilization_time_log_exist ==
                                        FALSE]),
        sep = "\n"
      )
      )
    }
  }

  # check if fossilization rates log file(s) exist if provided
  if (fossilization_rate_log != "") {
    do_fossilization_rate_log_exist <-
      file.exists(fossilization_rate_log)
    if (any(do_fossilization_rate_log_exist == FALSE) == TRUE) {
      # print out paths to files that don't exist
      stop(
        paste0(
        "Some fossilization_rate_log files do not exist:",
        paste0("\t",
               fossilization_rate_log[do_fossilization_rate_log_exist ==
                                        FALSE]),
        sep = "\n"
      )
      )
    }
  }

  # read in log files as lists of data.frames with readTrace()
  speciation_time <- readTrace(paths = speciation_time_log,
                               burnin = burnin)
  speciation_rate <- readTrace(paths = speciation_rate_log,
                               burnin = burnin)
  extinction_time <- readTrace(paths = extinction_time_log,
                               burnin = burnin)
  extinction_rate <- readTrace(paths = extinction_rate_log,
                               burnin = burnin)
  if (fossilization_time_log != "") {
    fossilization_time <- readTrace(paths = fossilization_time_log,
                                    burnin = burnin)
    fossilization_rate <- readTrace(paths = fossilization_rate_log,
                                    burnin = burnin)
  } else {
    fossilization_time <- NULL
    fossilization_rate <- NULL
  }

  # check if all parameter types have the same number of log files
  if (is.null(fossilization_time)) {
    trace_lengths_same <- identical(
      length(speciation_time),
      length(speciation_rate),
      length(extinction_time),
      length(extinction_rate)
    )
  } else {
    trace_lengths_same <- identical(
      length(speciation_time),
      length(speciation_rate),
      length(extinction_time),
      length(extinction_rate),
      length(fossilization_time),
      length(fossilization_rate)
    )
  }

  if (trace_lengths_same == FALSE) {
    stop("You must provide the same number of log files
         for each parameter type.")
  }
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
      speciation_time$`interval_times[0]` <-
        rep(0, nrow(speciation_time))
      extinction_time$`interval_times[0]` <-
        rep(0, nrow(extinction_time))
      if (!is.null(fossilization_time)) {
        fossilization_time$`interval_times[0]` <-
          rep(0, nrow(fossilization_time))
      }

      # Calculate the net-diversification and relative-extinction rates
      net_diversification_rate <-
        as.matrix(speciation_rate[, grepl("speciation",
                                          colnames(speciation_rate))]) -
        as.matrix(extinction_rate[, grepl("extinction",
                                          colnames(extinction_rate))])

      colnames(net_diversification_rate) <-
        paste(
          rep("net_div", times = ncol(net_diversification_rate)),
          rep("[", times = ncol(net_diversification_rate)),
          seq_len(ncol(net_diversification_rate)),
          rep("]", times = ncol(net_diversification_rate)),
          sep = ""
        )

      relative_extinction_rate <-
        as.matrix(extinction_rate[, grepl("extinction",
                                          colnames(extinction_rate))]) /
        as.matrix(speciation_rate[, grepl("speciation",
                                          colnames(speciation_rate))])

      colnames(relative_extinction_rate) <-
        paste(
          rep("rel_ext", times = ncol(relative_extinction_rate)),
          rep("[", times = ncol(relative_extinction_rate)),
          seq_len(ncol(relative_extinction_rate)),
          rep("]", times = ncol(relative_extinction_rate)),
          sep = ""
        )

      # return a list of processed data frames
      rates <- list(
        "speciation rate" = speciation_rate,
        "extinction rate" = extinction_rate,
        "net-diversification rate" = net_diversification_rate,
        "relative-extinction rate" = relative_extinction_rate,
        "fossilization rate" = fossilization_rate,
        "speciation time" = speciation_time,
        "extinction time" = extinction_time,
        "fossilization time" = fossilization_time
      )

      plotdata <-
        .makePlotData(rates = rates,
                      probs = probs,
                      summary = summary)
      return(plotdata)
    }
  }
}

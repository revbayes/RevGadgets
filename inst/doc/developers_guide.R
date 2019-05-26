## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ---- eval = FALSE-------------------------------------------------------
#  
#  #' Process Diversification Rates
#  #'
#  #'
#  #' Processing the output of a episodic diversification rate analysis with mass-extinction events.
#  #'
#  #' For processing the output of an episodic diversification rate analysis. processDivRates()
#  #' assumes that the epochs are fixed rather than inferred. Additionally, it assumes
#  #' that times correspond to rates such that the first rate parameter (i.e. speciation[1])
#  #' corresponds to the present. Conversely, the first time parameter
#  #' (i.e. interval_times[1]) corresponds to the first time interval after the present,
#  #' moving backwards in time. processDivRates() relies on readTrace and produces a list
#  #' object that can be read by plotDivRates() to vizualize the results. For now,
#  #' only one log file per parameter type is accepted (i.e. log files from multiple runs
#  #' must be combined before reading into the function).
#  #'
#  #'@param speciation_time_log (vector of character strings or single character string; "") Path to speciation times log file(s)
#  #'@param speciation_rate_log (vector of character strings or single character string; "") Path to speciation rates log file(s)
#  #'@param extinction_time_log (vector of character strings or single character string; "") Path to extinction times log file(s)
#  #'@param extinction_rate_log (vector of character strings or single character string; "") Path to extinction rates log file(s)
#  #'@param burnin (single numeric value; default = 0) Fraction of generations to
#  #' discard (if value provided is between 0 and 1) or number of generations (if
#  #' value provided is greater than 1). Passed to readTrace().
#  #'
#  #'@return List object with processed rate and time parameters.
#  #'
#  #'@examples
#  #'
#  #'\dontrun{
#  #'
#  #' speciation_time_file <- system.file("extdata",
#  #'     "epi_bd/primates_EBD_speciation_times.p", package="RevGadgets")
#  #' speciation_rate_file <- system.file("extdata",
#  #'     "epi_bd/primates_EBD_speciation_rates.p", package="RevGadgets")
#  #' extinction_time_file <- system.file("extdata",
#  #'     "epi_bd/primates_EBD_extinction_times.p", package="RevGadgets")
#  #' extinction_rate_file <- system.file("extdata",
#  #'     "epi_bd/primates_EBD_extinction_rates.p", package="RevGadgets")
#  #'
#  #' primates <- processDivRates(speciation_time_log = speciation_time_file,
#  #'                             speciation_rate_log = speciation_rate_file,
#  #'                             extinction_time_log = extinction_time_file,
#  #'                             extinction_rate_log = extinction_rate_file,
#  #'                             burnin = 0.25)
#  #'}
#  #'
#  #' @export
#  
#  processDivRates <- function(speciation_time_log = "",
#                              speciation_rate_log = "",
#                              extinction_time_log = "",
#                              extinction_rate_log = "",
#                              burnin = 0.25) {
#  
#    # enforce argument matching
#    if (is.character(speciation_time_log) == FALSE) stop("speciation_time_log must be a single character string")
#    if (is.character(speciation_rate_log) == FALSE) stop("speciation_rate_log must be a single character string")
#    if (is.character(extinction_time_log) == FALSE) stop("extinction_time_log must be a single character string")
#    if (is.character(extinction_rate_log) == FALSE) stop("extinction_rate_log must be a single character string")
#  
#    # check if speciation times log file(s) exist
#    do_speciation_time_log_exist <- file.exists(speciation_time_log)
#    if ( any(do_speciation_time_log_exist == FALSE) == TRUE ) {
#      # print out paths to files that don't exist
#      cat( "Some speciation_time_log files do not exist:",
#           paste0("\t",speciation_time_log[do_speciation_time_log_exist == FALSE]), sep="\n")
#      stop()
#    }
#  
#    # check if speciation rates log file(s) exist
#    do_speciation_rate_log_exist <- file.exists(speciation_rate_log)
#    if ( any(do_speciation_rate_log_exist == FALSE) == TRUE ) {
#      # print out paths to files that don't exist
#      cat( "Some speciation_rate_log files do not exist:",
#           paste0("\t",speciation_rate_log[do_speciation_rate_log_exist == FALSE]), sep="\n")
#      stop()
#    }
#  
#    # check if extinction times log file(s) exist
#    do_extinction_time_log_exist <- file.exists(extinction_time_log)
#    if ( any(do_extinction_time_log_exist == FALSE) == TRUE ) {
#      # print out paths to files that don't exist
#      cat( "Some extinction_time_log files do not exist:",
#           paste0("\t",extinction_time_log[do_extinction_time_log_exist == FALSE]), sep="\n")
#      stop()
#    }
#  
#    # check if extinction rates log file(s) exist
#    do_extinction_rate_log_exist <- file.exists(extinction_rate_log)
#    if ( any(do_extinction_rate_log_exist == FALSE) == TRUE ) {
#      # print out paths to files that don't exist
#      cat( "Some extinction_rate_log files do not exist:",
#           paste0("\t",extinction_rate_log[do_extinction_rate_log_exist == FALSE]), sep="\n")
#      stop()
#    }
#  
#    # read in log files as lists of data.frames with readTrace()
#    speciation_time <- readTrace(path = speciation_time_log,
#                                  burnin = burnin)
#    speciation_rate <- readTrace(path = speciation_rate_log,
#                                  burnin = burnin)
#    extinction_time <- readTrace(path = extinction_time_log,
#                                  burnin = burnin)
#    extinction_rate <- readTrace(path = extinction_rate_log,
#                                  burnin = burnin)
#  
#    # check if all parameter types have the same number of log files
#    trace_lengths_same <- identical(length(speciation_time),
#                                    length(speciation_rate),
#                                    length(extinction_time),
#                                    length(extinction_rate))
#  
#    if (trace_lengths_same == FALSE) {
#      stop("You must provide the same number of log files for each parameter type.")}
#    else if (trace_lengths_same == TRUE) {
#      if (length(speciation_time) == 0) {
#        stop("You must provide at least one log file per parameter type.")
#      } else if (length(speciation_time) > 1) {
#        stop("Currently, only one log file per parameter type is supported.")
#      } else if (length(speciation_time) == 1) {
#  
#        # convert single item lists to data frames
#        speciation_time <- speciation_time[[1]]
#        speciation_rate <- speciation_rate[[1]]
#        extinction_time <- extinction_time[[1]]
#        extinction_rate <- extinction_rate[[1]]
#  
#        # add in dummy distribution of 0 for time in the present
#        speciation_time$`interval_times[0]` <- rep(0, nrow(speciation_time))
#        extinction_time$`interval_times[0]` <- rep(0, nrow(extinction_time))
#  
#        # Calculate the net-diversification and relative-extinction rates
#        net_diversification_rate <- as.matrix(speciation_rate[,5:ncol(speciation_rate)]) -
#                                     as.matrix(extinction_rate[,5:ncol(extinction_rate)])
#  
#        colnames(net_diversification_rate) <- paste(rep("net_div", times = ncol(net_diversification_rate)),
#                                                     rep("[", times = ncol(net_diversification_rate)),
#                                                     1:ncol(net_diversification_rate),
#                                                     rep("]", times = ncol(net_diversification_rate)), sep = "")
#  
#        relative_extinction_rate <- as.matrix(extinction_rate[,5:ncol(extinction_rate)]) /
#                                     as.matrix(speciation_rate[,5:ncol(speciation_rate)])
#  
#        colnames(relative_extinction_rate) <- paste(rep("rel_ext", times = ncol(relative_extinction_rate)),
#                                                     rep("[", times = ncol(relative_extinction_rate)),
#                                                     1:ncol(relative_extinction_rate),
#                                                     rep("]", times = ncol(relative_extinction_rate)), sep = "")
#  
#          # return a list of processed data frames
#          res <- list("speciation rate" = speciation_rate,
#                      "extinction rate" = extinction_rate,
#                      "net-diversification rate" = net_diversification_rate,
#                      "relative-extinction rate" = relative_extinction_rate,
#                      "speciation time" = speciation_time,
#                      "extinction time" = extinction_time)
#          return(res)
#      }
#    }
#  }
#  
#  

## ---- eval = FALSE-------------------------------------------------------
#  #' Plot Diversification Rates
#  #'
#  #' Plots the output of a episodic diversification rate analysis
#  #'
#  #' Plots the output of diversification rate analyses. Takes as input the
#  #' output of processDivRates() and plotting parameters. Does not return an
#  #' object. For now, valid figure types include: "speciation rate",
#  #' "extinction rate","net-diversification rate", and "relative-extinction rate".
#  #' If colors are not provided (parameter col), the plot defaults to preset colors.
#  #'
#  #'
#  #' @param output (list; no default) The processed output for plotting (output of processDivRates()).
#  #' @param fig_types (character vector, c("speciation rate", "extinction rate", "net-diversification rate", "relative-extinction rate")) Which aspects of the model to visualize. See details for a complete description.
#  #' @param xlab (character; "million years ago") The label of the x-axis.
#  #' @param col (character; NULL) Colors used for printing in hex code. Must be of same length as fig_types.
#  #' @param col_alpha (character; "50") Alpha channel parameter for credible intervals plotting. May range from 00 (completely transparent) to FF (completely opaque).
#  #' @param offset (numeric; 0) Value to offset the x-axis by.
#  #' @param xaxt (character; "n") The type of x-axis to plot. By default, no x-axis is plotted (recommended).
#  #' @param yaxt (character; "s") The type of y-axis to plot.
#  #' @param ... Parameters delegated to various plotting functions.
#  #'
#  #' @return Plots diversification rates, does not return an object.
#  #'
#  #' @examples
#  #'
#  #' \dontrun{
#  #' # first run processDivRates()
#  #' speciation_time_file <- system.file("extdata",
#  #'     "epi_bd/primates_EBD_speciation_times.p", package="RevGadgets")
#  #' speciation_rate_file <- system.file("extdata",
#  #'     "epi_bd/primates_EBD_speciation_rates.p", package="RevGadgets")
#  #' extinction_time_file <- system.file("extdata",
#  #'     "epi_bd/primates_EBD_extinction_times.p", package="RevGadgets")
#  #' extinction_rate_file <- system.file("extdata",
#  #'     "epi_bd/primates_EBD_extinction_rates.p", package="RevGadgets")
#  #'
#  #' primates <- processDivRates(speciation_time_log = speciation_time_file,
#  #'                             speciation_rate_log = speciation_rate_file,
#  #'                             extinction_time_log = extinction_time_file,
#  #'                             extinction_rate_log = extinction_rate_file,
#  #'                             burnin = 0.25)
#  #' # then plot results:
#  #' plotDivRates(output = primates)
#  #' }
#  #'
#  #' @export
#  #' @importFrom graphics plot polygon
#  #' @importFrom stats quantile
#  
#  plotDivRates <- function(output,
#                           fig_types = c("speciation rate",
#                                         "extinction rate",
#                                         "net-diversification rate",
#                                         "relative-extinction rate"),
#                           xlab = "million years ago",
#                           col = NULL,
#                           col_alpha = "50",
#                           offset = 0,
#                           xaxt = "n",
#                           yaxt = "s",
#                           ...){
#    # Enforce argument matching
#    if (is.list(output) == FALSE) stop("output must be list of processed rates")
#    if (is.character(fig_types) == FALSE) stop("fig_types must be a character string or vector of strings")
#    if (is.character(xlab) == FALSE) stop("xlab must be a single character string")
#  
#    # check length of color vector and ensure colors are formatted as hex codes
#    if (is.null(col) == FALSE) {
#      if (length(col) != length(fig_types))
#        stop("col and fig_types must be of the same length")
#      bad_cols <- rep("FALSE", times = length(col))
#      for (i in 1:length(col)) {
#        col_split <- unlist(strsplit(col[i], split = ""))
#        if (length(col_split) != 7 | col_split[1] != "#") {
#          bad_cols[i] <- TRUE
#        }
#      }
#      if ( any(bad_cols == TRUE) == TRUE ) {
#        # print out colors that are 'bad'
#        cat( "Some colors are not properly formatted hex codes:",
#             paste0("\t",col[bad_cols == TRUE]), sep="\n")
#        stop()
#      }
#    }
#  
#    if (is.character(col_alpha) == FALSE) stop("col_alpha must be a character string")
#    if (is.numeric(offset) == FALSE) stop("offset must be a numeric value")
#    xaxt <- match.arg(xaxt, choices = c("n", "s"))
#    yaxt <- match.arg(yaxt, choices = c("n", "s"))
#  
#  
#    # Check that fig type is valid
#    valid_fig_types <- c("speciation rate",
#                         "extinction rate",
#                         "net-diversification rate",
#                         "relative-extinction rate")
#  
#    invalid_fig_types <- fig_types[!fig_types %in% valid_fig_types]
#  
#    if ( length(invalid_fig_types) > 0 ) {
#      stop("\nThe following figure types are invalid: ", paste(invalid_fig_types, collapse = ", "), ".",
#           "\nValid options are: ", paste(valid_fig_types, collapse = ", ") ,".")
#    }
#  
#  
#    # Make color vector
#    if ( is.null(col) ) {
#      col <- c("speciation rate"="#984EA3",
#               "speciation shift times"="#984EA3",
#               "speciation Bayes factors"="#984EA3",
#               "extinction rate"="#E41A1C",
#               "extinction shift times"="#E41A1C",
#               "extinction Bayes factors"="#E41A1C",
#               "fossilization rate"="#ffa31a",
#               "fossilization shift times"="#ffa31a",
#               "fossilization Bayes factors"="#ffa31a",
#               "net-diversification rate"="#377EB8",
#               "relative-extinction rate"="#FF7F00",
#               "mass extinction times"="#4DAF4A",
#               "mass extinction Bayes factors"="#4DAF4A")
#    } else {
#      names(col) <- fig_types
#    }
#  
#    # Compute the axes
#    intervals <- output$`speciation time`[1,grep("interval", colnames(output$`speciation time`))]
#    tree_age <- max(intervals)
#    num_intervals <- length(intervals) - 1
#    plot_at <- 0:num_intervals
#    interval_size <- tree_age/num_intervals
#    labels <- pretty( c(0,tree_age) ) + offset
#    labels_at <- num_intervals - ( pretty(c(0,tree_age)) / interval_size )
#    ages <- seq( 0,tree_age, length.out = num_intervals + 1 )
#  
#    for( type in fig_types ) {
#  
#        this_output <- output[[type]][ ,grep("[0-9]", colnames(output[[type]])) ]
#        mean_this_output <- colMeans(this_output)
#        quantiles_this_output <- apply(this_output, 2, quantile, prob = c(0.025,0.975))
#  
#        # find the limits of the y-axis
#        # we always want the speciation and extinction rates to be on the same scale
#        if ( type %in% c("speciation rate","extinction rate")) {
#  
#          quantiles_speciation <-
#            apply(output[["speciation rate"]][ ,grep("[0-9]", colnames(output[["speciation rate"]])) ],
#                  2, quantile, prob=c(0.025,0.975))
#          quantiles_extinction <-
#            apply(output[["extinction rate"]][ ,grep("[0-9]", colnames(output[["extinction rate"]])) ],
#                  2, quantile, prob=c(0.025,0.975))
#          ylim <- c( min(0, quantiles_speciation, quantiles_extinction),
#                     max(quantiles_speciation, quantiles_extinction) )
#  
#        } else {
#          ylim <- c( min(0, quantiles_this_output),
#                     max(quantiles_this_output) )
#        }
#  
#        #### plot
#          plot( x = plot_at, y = c( rev(mean_this_output)),
#                type = "l", ylim = ylim, xaxt = xaxt, col = col[type],
#                ylab = "rate", main = type, xlab = xlab,...)
#          polygon( x = c( 0:ncol(quantiles_this_output),
#                          ncol(quantiles_this_output):0 ),
#                   y = c( rev(c(quantiles_this_output[1,1],
#                            quantiles_this_output[1,])),
#                          c(quantiles_this_output[2,1],
#                                quantiles_this_output[2,])),
#                   border = NA, col = paste(col[type], col_alpha, sep="") )
#          axis(1, at = labels_at, labels = labels)
#        }
#  
#      }

## ---- eval = FALSE-------------------------------------------------------
#  Encoding: UTF-8
#  LazyData: true
#  RoxygenNote: 6.1.0
#  Imports: ape, phytools, ggtree, treeio, dplyr, tidytree, methods
#  Suggests: testthat,
#      knitr,
#      rmarkdown

## ---- eval = FALSE-------------------------------------------------------
#  #' Plot Diversification Rates
#  #'
#  #' Plots the output of a episodic diversification rate analysis
#  #'
#  #' Plots the output of diversification rate analyses. Takes as input the
#  #' output of processDivRates() and plotting parameters. Does not return an
#  #' object. For now, valid figure types include: "speciation rate",
#  #' "extinction rate","net-diversification rate", and "relative-extinction rate".
#  #' If colors are not provided (parameter col), the plot defaults to preset colors.
#  #'
#  #'
#  #' @param output (list; no default) The processed output for plotting (output of processDivRates()).
#  #' @param fig_types (character vector, c("speciation rate", "extinction rate", "net-diversification rate", "relative-extinction rate")) Which aspects of the model to visualize. See details for a complete description.
#  #' @param xlab (character; "million years ago") The label of the x-axis.
#  #' @param col (character; NULL) Colors used for printing in hex code. Must be of same length as fig_types.
#  #' @param col_alpha (character; "50") Alpha channel parameter for credible intervals plotting. May range from 00 (completely transparent) to FF (completely opaque).
#  #' @param offset (numeric; 0) Value to offset the x-axis by.
#  #' @param xaxt (character; "n") The type of x-axis to plot. By default, no x-axis is plotted (recommended).
#  #' @param yaxt (character; "s") The type of y-axis to plot.
#  #' @param ... Parameters delegated to various plotting functions.
#  #'
#  #' @return Plots diversification rates, does not return an object.
#  #'
#  #' @examples
#  #'
#  #' \dontrun{
#  #' # first run processDivRates()
#  #' speciation_time_file <- system.file("extdata",
#  #'     "epi_bd/primates_EBD_speciation_times.p", package="RevGadgets")
#  #' speciation_rate_file <- system.file("extdata",
#  #'     "epi_bd/primates_EBD_speciation_rates.p", package="RevGadgets")
#  #' extinction_time_file <- system.file("extdata",
#  #'     "epi_bd/primates_EBD_extinction_times.p", package="RevGadgets")
#  #' extinction_rate_file <- system.file("extdata",
#  #'     "epi_bd/primates_EBD_extinction_rates.p", package="RevGadgets")
#  #'
#  #' primates <- processDivRates(speciation_time_log = speciation_time_file,
#  #'                             speciation_rate_log = speciation_rate_file,
#  #'                             extinction_time_log = extinction_time_file,
#  #'                             extinction_rate_log = extinction_rate_file,
#  #'                             burnin = 0.25)
#  #' # then plot results:
#  #' plotDivRates(output = primates)
#  #' }
#  #'
#  #' @export
#  #' @importFrom graphics plot polygon
#  #' @importFrom stats quantile

## ---- eval = FALSE-------------------------------------------------------
#  #'@param output (list; no default) The processed output for plotting (output of processDivRates()).

## ----  eval = FALSE------------------------------------------------------
#  #' @export
#  #' @importFrom graphics plot polygon
#  #' @importFrom stats quantile

## ---- eval = FALSE-------------------------------------------------------
#  
#  plotDivRates <- function(output,
#                           fig_types = c("speciation rate",
#                                         "extinction rate",
#                                         "net-diversification rate",
#                                         "relative-extinction rate"),
#                           xlab = "million years ago",
#                           col = NULL,
#                           col_alpha = "50",
#                           offset = 0,
#                           xaxt = "n",
#                           yaxt = "s",
#                           ...){
#    # Enforce argument matching
#    if (is.list(output) == FALSE) stop("output must be list of processed rates")
#    if (is.character(fig_types) == FALSE) stop("fig_types must be a character string or vector of strings")
#    if (is.character(xlab) == FALSE) stop("xlab must be a single character string")
#  
#    # check length of color vector and ensure colors are formatted as hex codes
#    if (is.null(col) == FALSE) {
#      if (length(col) != length(fig_types))
#        stop("col and fig_types must be of the same length")
#      bad_cols <- rep("FALSE", times = length(col))
#      for (i in 1:length(col)) {
#        col_split <- unlist(strsplit(col[i], split = ""))
#        if (length(col_split) != 7 | col_split[1] != "#") {
#          bad_cols[i] <- TRUE
#        }
#      }
#      if ( any(bad_cols == TRUE) == TRUE ) {
#        # print out colors that are 'bad'
#        cat( "Some colors are not properly formatted hex codes:",
#             paste0("\t",col[bad_cols == TRUE]), sep="\n")
#        stop()
#      }
#    }
#  
#    if (is.character(col_alpha) == FALSE) stop("col_alpha must be a character string")
#    if (is.numeric(offset) == FALSE) stop("offset must be a numeric value")
#    xaxt <- match.arg(xaxt, choices = c("n", "s"))
#    yaxt <- match.arg(yaxt, choices = c("n", "s"))
#  
#  
#    # Check that fig type is valid
#    valid_fig_types <- c("speciation rate",
#                         "extinction rate",
#                         "net-diversification rate",
#                         "relative-extinction rate")
#  
#    invalid_fig_types <- fig_types[!fig_types %in% valid_fig_types]
#  
#    if ( length(invalid_fig_types) > 0 ) {
#      stop("\nThe following figure types are invalid: ", paste(invalid_fig_types, collapse = ", "), ".",
#           "\nValid options are: ", paste(valid_fig_types, collapse = ", ") ,".")
#    }
#  
#  
#    # Make color vector
#    if ( is.null(col) ) {
#      col <- c("speciation rate"="#984EA3",
#               "speciation shift times"="#984EA3",
#               "speciation Bayes factors"="#984EA3",
#               "extinction rate"="#E41A1C",
#               "extinction shift times"="#E41A1C",
#               "extinction Bayes factors"="#E41A1C",
#               "fossilization rate"="#ffa31a",
#               "fossilization shift times"="#ffa31a",
#               "fossilization Bayes factors"="#ffa31a",
#               "net-diversification rate"="#377EB8",
#               "relative-extinction rate"="#FF7F00",
#               "mass extinction times"="#4DAF4A",
#               "mass extinction Bayes factors"="#4DAF4A")
#    } else {
#      names(col) <- fig_types
#    }
#  
#    # Compute the axes
#    intervals <- output$`speciation time`[1,grep("interval", colnames(output$`speciation time`))]
#    tree_age <- max(intervals)
#    num_intervals <- length(intervals) - 1
#    plot_at <- 0:num_intervals
#    interval_size <- tree_age/num_intervals
#    labels <- pretty( c(0,tree_age) ) + offset
#    labels_at <- num_intervals - ( pretty(c(0,tree_age)) / interval_size )
#    ages <- seq( 0,tree_age, length.out = num_intervals + 1 )
#  
#    for( type in fig_types ) {
#  
#        this_output <- output[[type]][ ,grep("[0-9]", colnames(output[[type]])) ]
#        mean_this_output <- colMeans(this_output)
#        quantiles_this_output <- apply(this_output, 2, quantile, prob = c(0.025,0.975))
#  
#        # find the limits of the y-axis
#        # we always want the speciation and extinction rates to be on the same scale
#        if ( type %in% c("speciation rate","extinction rate")) {
#  
#          quantiles_speciation <-
#            apply(output[["speciation rate"]][ ,grep("[0-9]", colnames(output[["speciation rate"]])) ],
#                  2, quantile, prob=c(0.025,0.975))
#          quantiles_extinction <-
#            apply(output[["extinction rate"]][ ,grep("[0-9]", colnames(output[["extinction rate"]])) ],
#                  2, quantile, prob=c(0.025,0.975))
#          ylim <- c( min(0, quantiles_speciation, quantiles_extinction),
#                     max(quantiles_speciation, quantiles_extinction) )
#  
#        } else {
#          ylim <- c( min(0, quantiles_this_output),
#                     max(quantiles_this_output) )
#        }
#  
#        #### plot
#          plot( x = plot_at, y = c( rev(mean_this_output)),
#                type = "l", ylim = ylim, xaxt = xaxt, col = col[type],
#                ylab = "rate", main = type, xlab = xlab,...)
#          polygon( x = c( 0:ncol(quantiles_this_output),
#                          ncol(quantiles_this_output):0 ),
#                   y = c( rev(c(quantiles_this_output[1,1],
#                            quantiles_this_output[1,])),
#                          c(quantiles_this_output[2,1],
#                                quantiles_this_output[2,])),
#                   border = NA, col = paste(col[type], col_alpha, sep="") )
#          axis(1, at = labels_at, labels = labels)
#        }
#  
#      }

## ---- eval = FALSE-------------------------------------------------------
#  plotDivRates <- function(output,
#                           fig_types = c("speciation rate",
#                                         "extinction rate",
#                                         "net-diversification rate",
#                                         "relative-extinction rate"),
#                           xlab = "million years ago",
#                           col = NULL,
#                           col_alpha = "50",
#                           offset = 0,
#                           xaxt = "n",
#                           yaxt = "s",
#                           ...){

## ---- eval = FALSE-------------------------------------------------------
#    # Enforce argument matching
#    if (is.list(output) == FALSE) stop("output must be list of processed rates")
#    if (is.character(fig_types) == FALSE) stop("fig_types must be a character string or vector of strings")
#    if (is.character(xlab) == FALSE) stop("xlab must be a single character string")
#  
#    # check length of color vector and ensure colors are formatted as hex codes
#    if (is.null(col) == FALSE) {
#      if (length(col) != length(fig_types))
#        stop("col and fig_types must be of the same length")
#      bad_cols <- rep("FALSE", times = length(col))
#      for (i in 1:length(col)) {
#        col_split <- unlist(strsplit(col[i], split = ""))
#        if (length(col_split) != 7 | col_split[1] != "#") {
#          bad_cols[i] <- TRUE
#        }
#      }
#      if ( any(bad_cols == TRUE) == TRUE ) {
#        # print out colors that are 'bad'
#        cat( "Some colors are not properly formatted hex codes:",
#             paste0("\t",col[bad_cols == TRUE]), sep="\n")
#        stop()
#      }
#    }
#  
#    if (is.character(col_alpha) == FALSE) stop("col_alpha must be a character string")
#    if (is.numeric(offset) == FALSE) stop("offset must be a numeric value")
#    xaxt <- match.arg(xaxt, choices = c("n", "s"))
#    yaxt <- match.arg(yaxt, choices = c("n", "s"))
#  
#  
#    # Check that fig type is valid
#    valid_fig_types <- c("speciation rate",
#                         "extinction rate",
#                         "net-diversification rate",
#                         "relative-extinction rate")
#  
#    invalid_fig_types <- fig_types[!fig_types %in% valid_fig_types]
#  
#    if ( length(invalid_fig_types) > 0 ) {
#      stop("\nThe following figure types are invalid: ", paste(invalid_fig_types, collapse = ", "), ".",
#           "\nValid options are: ", paste(valid_fig_types, collapse = ", ") ,".")
#    }
#  

## ---- eval = FALSE-------------------------------------------------------
#  # read in tree(s) of type nexus or newick
#  
#  # Make color vector
#    if ( is.null(col) ) {
#      col <- c("speciation rate"="#984EA3",
#               "speciation shift times"="#984EA3",
#               "speciation Bayes factors"="#984EA3",
#               "extinction rate"="#E41A1C",
#               "extinction shift times"="#E41A1C",
#               "extinction Bayes factors"="#E41A1C",
#               "fossilization rate"="#ffa31a",
#               "fossilization shift times"="#ffa31a",
#               "fossilization Bayes factors"="#ffa31a",
#               "net-diversification rate"="#377EB8",
#               "relative-extinction rate"="#FF7F00",
#               "mass extinction times"="#4DAF4A",
#               "mass extinction Bayes factors"="#4DAF4A")
#    } else {
#      names(col) <- fig_types
#    }
#  
#    # Compute the axes
#    intervals <- output$`speciation time`[1,grep("interval", colnames(output$`speciation time`))]
#    tree_age <- max(intervals)
#    num_intervals <- length(intervals) - 1
#    plot_at <- 0:num_intervals
#    interval_size <- tree_age/num_intervals
#    labels <- pretty( c(0,tree_age) ) + offset
#    labels_at <- num_intervals - ( pretty(c(0,tree_age)) / interval_size )
#    ages <- seq( 0,tree_age, length.out = num_intervals + 1 )
#  
#    for( type in fig_types ) {
#  
#        this_output <- output[[type]][ ,grep("[0-9]", colnames(output[[type]])) ]
#        mean_this_output <- colMeans(this_output)
#        quantiles_this_output <- apply(this_output, 2, quantile, prob = c(0.025,0.975))
#  
#        # find the limits of the y-axis
#        # we always want the speciation and extinction rates to be on the same scale
#        if ( type %in% c("speciation rate","extinction rate")) {
#  
#          quantiles_speciation <-
#            apply(output[["speciation rate"]][ ,grep("[0-9]", colnames(output[["speciation rate"]])) ],
#                  2, quantile, prob=c(0.025,0.975))
#          quantiles_extinction <-
#            apply(output[["extinction rate"]][ ,grep("[0-9]", colnames(output[["extinction rate"]])) ],
#                  2, quantile, prob=c(0.025,0.975))
#          ylim <- c( min(0, quantiles_speciation, quantiles_extinction),
#                     max(quantiles_speciation, quantiles_extinction) )
#  
#        } else {
#          ylim <- c( min(0, quantiles_this_output),
#                     max(quantiles_this_output) )
#        }
#  
#        #### plot
#          plot( x = plot_at, y = c( rev(mean_this_output)),
#                type = "l", ylim = ylim, xaxt = xaxt, col = col[type],
#                ylab = "rate", main = type, xlab = xlab,...)
#          polygon( x = c( 0:ncol(quantiles_this_output),
#                          ncol(quantiles_this_output):0 ),
#                   y = c( rev(c(quantiles_this_output[1,1],
#                            quantiles_this_output[1,])),
#                          c(quantiles_this_output[2,1],
#                                quantiles_this_output[2,])),
#                   border = NA, col = paste(col[type], col_alpha, sep="") )
#          axis(1, at = labels_at, labels = labels)
#        }
#  }
#  

## ---- eval = FALSE-------------------------------------------------------
#  
#   # return a list of processed data frames
#   res <- list("speciation rate" = speciation_rate,
#               "extinction rate" = extinction_rate,
#               "net-diversification rate" = net_diversification_rate,
#               "relative-extinction rate" = relative_extinction_rate,
#               "speciation time" = speciation_time,
#               "extinction time" = extinction_time)
#   return(res)
#  

## ---- eval = FALSE-------------------------------------------------------
#  context("tests the processDivRates function")
#  
#  # relies heavily on readTrace(), so we only test for elements not
#  # included in the readTrace() testing.
#  
#  test_that("processes birth-death scripts", {
#    file_spectimes <- system.file("extdata", "epi_bd/primates_EBD_speciation_times.p", package="RevGadgets")
#    file_specrates <- system.file("extdata", "epi_bd/primates_EBD_speciation_rates.p", package="RevGadgets")
#    file_exttimes <- system.file("extdata", "epi_bd/primates_EBD_extinction_times.p", package="RevGadgets")
#    file_extrates <- system.file("extdata", "epi_bd/primates_EBD_extinction_rates.p", package="RevGadgets")
#  
#    primates <- processDivRates(speciation_time_log = file_spectimes,
#                                speciation_rate_log = file_specrates,
#                                extinction_time_log = file_exttimes,
#                                extinction_rate_log = file_extrates,
#                                burnin = 0.25)
#    expect_equal(length(primates), 6)
#    expect_equal(class(primates), "list")
#    expect_equal(nrow(primates[[1]]), 3750)
#  })
#  


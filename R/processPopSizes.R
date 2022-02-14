#' Process Population Sizes
#'
#'
#' Processing the output of a coalescent demographic analysis.
#'
#' For processing the output of a coalescent demographic analysis.
#' processPopSizes() can process output from different analyses,
#' either an event-based skyline analysis (method = "events"),
#' or an analysis with user-defined intervals (method = "specified"),
#' or an analysis with only one population size and no interval times
#' (method = "constant").
#' processPopSizes() assumes that the the first size parameter (i.e. population_size[1])
#' corresponds to the present. processPopSizes() partly
#' relies on readTrace and produces a list object that can be read by
#' plotPopSizes() to visualize the results. For now, only one log file per
#' parameter type is accepted (i.e. log files from multiple runs must be
#' combined before reading into the function).
#'
#' @param population_size_log (vector of character strings or
#' single character string; "") Path to population sizes log file(s)
#' @param interval_change_points_log (vector of character strings or
#' single character string; "") Path to interval change points log file(s)
#' @param method which method was chosen for the analysis, "events" - event-based coalescent process,
#' "specified" - coalescent process with user-defined interval times, "constant" - constant coalescent process
#' @param burnin (single numeric value; default = 0) Fraction of generations to
#'  discard (if value provided is between 0 and 1) or number of generations (if
#'  value provided is greater than 1).
#' @param probs (numeric vector; c(0.025, 0.975)) a vector of length two
#' containing the upper and lower bounds for the confidence intervals.
#' @param summary typically "mean" or "median"; the metric to summarize the
#' posterior distribution. Defaults to "mean"
#' @param num_grid_points only if method = "events", defines the number of grid points through time for which to 
#' evaluate the demographic functions
#' @return List object with processed rate and, if applicable, time parameters.
#'
#' @export
#' @importFrom tibble as_tibble tibble

processPopSizes <- function(population_size_log = "",
                            interval_change_points_log = "",
                            method = "events",
                            burnin = 0.25,
                            probs = c(0.025, 0.975),
                            summary = "mean",
                            num_grid_points = 100){
  if (method == "events"){
    times <- read_file(interval_change_points_log, burnin = burnin)
    pop_size <- read_file(population_size_log, burnin = burnin)
    
    orders <- lapply(times, order)
    
    pop_size_ordered <- list() 
    for (i in seq_along(pop_size)){
      res <- c(pop_size[[i]][1], pop_size[[i]][-1][orders[[i]]])
      pop_size_ordered[[i]] <- res
    }
    
    pop_size_trajectories <- list()
    for (i in seq_along(pop_size_ordered)){
      if(length(times[[i]]) > 0){
        f <- approxfun(sort(times[[i]]),
                       tail(pop_size_ordered[[i]], n = -1), 
                       yleft = pop_size_ordered[[i]][1],
                       yright = tail(pop_size_ordered[[i]], n = 1),
                       method = "constant")
      }else{
        f <- function(t) pop_size_ordered[[i]][1] + t*0
      }
      pop_size_trajectories[[i]] <- f
    }
    
    x <- seq(0, suppressWarnings(max(sapply(times, max))), length.out = num_grid_points)
    m <- sapply(pop_size_trajectories, function(e) e(x))
    quantiles <- apply(m, 1, function(x) quantile(x, probs = probs))
    
    plotdf <- as_tibble(apply(m, 1, summary))
    plotdf$lower <- quantiles[1,]
    plotdf$upper <- quantiles[2,]
    plotdf$time <- x
    
  } else if (method == "specified"){
    interval_time <- readTrace(paths = interval_change_points_log,
                               burnin = burnin)[[1]]
    pop_size <- readTrace(paths = population_size_log,
                          burnin = burnin)[[1]]

    interval_time$`interval_times[0]` <- rep(0, nrow(interval_time))

    rates <- list(
      "population size" = pop_size,
      "coalescent time" = interval_time
    )

    plotdf <- .makePlotData(rates = rates, probs = probs, summary = summary)
    
  } else if (method == "constant") {
    pop_size <- readTrace(paths = population_size_log)[[1]]
    interval_time <- c("interval_times[0]" = 0, "interval_times[1]" = Inf)
    
    pop_sizes <- pop_size[, ncol(pop_size)]
    
    summary_size <- apply(as.matrix(pop_sizes), 2, FUN = summary)
    quantiles <- apply(as.matrix(pop_sizes), 2,
                       quantile,
                       probs = probs)
    
    plotdf <- dplyr::tibble(.rows = length(summary_size))
    plotdf["value"] <- summary_size
    plotdf["lower"] <- quantiles[1, ]
    plotdf["upper"] <- quantiles[2, ]
    plotdf$time <- interval_time[1]
    plotdf$time_end <- interval_time[2]

  }
  
  return(plotdf)
}

read_file <- function(path, burnin = 0.25) {
  
  `%>%` <- dplyr::`%>%`
  
  res <- path %>% 
    readLines() %>%
    tail(n = -1)
    
  if (burnin >= length(res))
    stop("Burnin larger than provided trace file")

  if (burnin >= 1) {
    res <- res[(burnin + 1):length(res)]
  } else if (burnin < 1 & burnin > 0) {
    discard <- ceiling(burnin * length(res))
    res <- res[(discard + 1):length(res)]
  } else if (burnin == 0) {
    res <- res
  } else {
    stop("What have you done?")
  }
  
  res <- res %>%
    strsplit("\t") %>% 
    lapply(function(x) tail(x, n = -4)) %>%
    lapply(as.numeric)
  
  return(res)
}

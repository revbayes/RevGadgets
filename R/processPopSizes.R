#' Process Population Sizes
#'
#'
#' Processing the output of a coalescent demographic analysis.
#'
#' For processing the output of a coalescent demographic analysis.
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
#' single character string; "") Path to interval change points log file(s).
#' If not given, a constant process with only one population size is assumed.
#' @param burnin (single numeric value; default = 0) Fraction of generations to
#'  discard (if value provided is between 0 and 1) or number of generations (if
#'  value provided is greater than 1).
#' @param probs (numeric vector; c(0.025, 0.975)) a vector of length two
#' containing the upper and lower bounds for the confidence intervals.
#' @param summary typically "mean" or "median"; the metric to summarize the
#' posterior distribution. Defaults to "mean"
#' @param num_grid_points defines the number of grid points through time for which to 
#' evaluate the demographic functions
#' @param max_age defines the maximal age up to which the demographic functions should be evaluated.
#' If not provided, it will either be automatically set to 1e5 (in case of a constant process) or
#' to the maximal age provided with the interval_change_points_log.
#' @return List object with processed rate and, if applicable, time parameters.
#'
#' @export
#' @importFrom tibble as_tibble tibble

processPopSizes <- function(population_size_log = "",
                            interval_change_points_log = "",
                            burnin = 0.25,
                            probs = c(0.025, 0.975),
                            summary = "mean",
                            num_grid_points = 100,
                            max_age = NULL){
  constant_dem = FALSE
  
  if (interval_change_points_log == ""){
    constant_dem = TRUE
  }
  
  if (constant_dem == TRUE){
    pop_size <- readTrace(paths = population_size_log)[[1]]
    pop_sizes <- pop_size[, ncol(pop_size)]

    summary_size <- apply(as.matrix(pop_sizes), 2, FUN = summary)
    quantiles <- apply(as.matrix(pop_sizes), 2,
                       quantile,
                       probs = probs)
    
    if (is.null(max_age)){
      x <- seq(0, 1e5, length.out = num_grid_points)
    } else {
      x <- seq(0, max_age, length.out = num_grid_points)
    }

    plotdf <- dplyr::tibble(.rows = length(x))
    plotdf["value"] <- rep(summary_size, length(x))
    plotdf["lower"] <- rep(quantiles[1, ], length(x))
    plotdf["upper"] <- rep(quantiles[2, ], length(x))
    plotdf$time <- x

  } else {
    pop_size <- read_file(population_size_log, burnin = burnin)
    times <- read_file(interval_change_points_log, burnin = burnin)
    
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
    
    if (is.null(max_age)){
      x <- seq(0, suppressWarnings(max(sapply(times, max))), length.out = num_grid_points)
    } else {
      x <- seq(0, max_age, length.out = num_grid_points)
    }
    
    m <- sapply(pop_size_trajectories, function(e) e(x))
    quantiles <- apply(m, 1, function(x) quantile(x, probs = probs))
    
    plotdf <- as_tibble(apply(m, 1, summary))
    plotdf$lower <- quantiles[1,]
    plotdf$upper <- quantiles[2,]
    plotdf$time <- x
  }

  return(plotdf)
}

read_file <- function(path, burnin = 0.25) {
  
  `%>%` <- dplyr::`%>%`
  
  res <- path %>% 
    readLines() %>%
    tail(n = -1)
  
  names <- path %>% 
    readLines() %>%
    head(n = 1) %>%
    strsplit("\t")
    
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
  
  names_to_exclude = c("Iteration|Replicate_ID|Posterior|Likelihood|Prior")
  cols_to_exclude = length(grep(pattern = names_to_exclude, names[[1]]))
  
  res <- res %>%
    strsplit("\t") %>% 
    lapply(function(x) tail(x, n = -cols_to_exclude)) %>%
    lapply(as.numeric)
  
  return(res)
}

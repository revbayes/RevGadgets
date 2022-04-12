#' Remove Burnin
#'
#' Removes burnin from MCMC trace
#'
#' Removes burnin from an MCMC trace, such as the output of readTrace(). If
#' multiple traces are provided, this function will remove the burnin from each.
#'
#' @param trace (list of data frames; no default) Name of a list of data frames,
#' such as produced by readTrace().
#'
#' @param burnin (single numeric value; 0.1) Fraction of generations to
#' discard (if value provided is between 0 and 1) or number of generations (if
#' value provided is greater than 1).
#'
#' @return List of dataframes (of length 1 if only 1 log file provided).
#'
#' @examples
#'
#' \donttest{
#' # download the example dataset to working directory
#' url_gtr <-
#'     "https://revbayes.github.io/tutorials/intro/data/primates_cytb_GTR.log"
#' dest_path_gtr <- "primates_cytb_GTR.log"
#' download.file(url_gtr, dest_path_gtr)
#'
#' # to run on your own data, change this to the path to your data file
#' file_single <- dest_path_gtr
#'
#' one_trace <- readTrace(paths = file_single)
#' one_trace_burnin <- removeBurnin(trace = one_trace, burnin = 0.1)
#'
#' # remove file
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_gtr)
#' }
#'
#' @export
#'

removeBurnin <- function(trace, burnin) {
  if (is.list(trace) == FALSE)
    stop("trace should be a list of data frames")
  if (is.data.frame(trace[[1]]) == FALSE)
    stop("trace should be a list of data frames")
  if (is.numeric(burnin) == FALSE)
    stop("burnin must be a single numeric value")
  if (burnin < 0)
    stop("burnin must be a positive value")

  for (i in seq_len(length(trace))) {
    if (burnin >= nrow(trace[[i]]))
      stop("Burnin larger than provided trace file")

    if (burnin >= 1) {
      trace[[i]] <- trace[[i]][(burnin + 1):nrow(trace[[i]]),]
    } else if (burnin < 1 & burnin > 0) {
      discard <- ceiling(burnin * nrow(trace[[i]]))
      trace[[i]] <- trace[[i]][(discard + 1):nrow(trace[[i]]),]
    } else if (burnin == 0) {
      trace[[i]] <- trace[[i]]
    }
  }
  return(trace)
}

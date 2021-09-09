#' Combine traces
#'
#' Combine traces into one trace file
#'
#' Combines multiple traces from independent MCMC replicates
#' into one trace file.
#'
#' @param traces (list of data frames; no default) Name of a list of data
#' frames, such as produced by readTrace().
#'
#' @param burnin (single numeric value; default = 0.0) Fraction of generations
#' to discard (if value provided is between 0 and 1) or number of generations
#' to discard (if value provided is greater than 1) before combining the
#' samples.
#'
#' @return combineTraces() returns a list of data frames of length 1,
#' corresponding to the combination of the provided samples.
#'
#' @examples
#'
#' \donttest{
#'
#' #' # download the example dataset to working directory
#' url_1 <-
#' "https://revbayes.github.io/tutorials/intro/data/primates_cytb_GTR_run_1.log"
#' dest_path_1 <- "primates_cytb_GTR_run_1.log"
#' download.file(url_1, dest_path_1)
#'
#' url_2 <-
#' "https://revbayes.github.io/tutorials/intro/data/primates_cytb_GTR_run_2.log"
#' dest_path_2 <- "primates_cytb_GTR_run_2.log"
#' download.file(url_2, dest_path_2)
#'
#' # to run on your own data, change this to the path to your data file
#' file_1 <- dest_path_1
#' file_2 <- dest_path_2
#'
#' # read in the multiple trace files
#' multi_trace <- readTrace(path = c(file_1, file_2), burnin = 0.0)
#'
#' # combine samples after discarding 10% burnin
#' combined_trace <- combineTraces(trace = multi_trace,
#'                                 burnin = 0.1)
#'
#' # remove files
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_1, dest_path_2)
#' }
#'
#'
#' @export


combineTraces <- function(traces, burnin = 0.0) {
  # enforce argument matching
  if (is.list(traces) == FALSE)
    stop("trace should be a list of data frames")
  if (is.data.frame(traces[[1]]) == FALSE)
    stop("trace should be a list of data frames")
  if (burnin < 0)
    stop("burnin must be a positive value")

  # only combine if there are multiple traces
  num_traces <- length(traces)
  if (num_traces < 2) {
    stop("Provided a list of 1 trace; can only combine multiple traces")
  }

  # check that the file headings match for all traces
  header <- vector("list", num_traces)
  for (i in 1:num_traces) {
    header[[i]] <- colnames(traces[[i]])
  }
  all_headers <- unique(unlist(header))
  for (i in 1:num_traces) {
    if (setequal(all_headers, header[[i]]) == FALSE) {
      stop("Not all headers of trace files match")
    }
  }

  # order columns, discard burnin
  post_burnin_samples <- vector("list", num_traces)
  for (i in 1:num_traces) {
    # get the trace
    this_trace <- traces[[i]]

    # order the trace
    this_trace <- this_trace[, all_headers]

    # discard burnin
    if (burnin >= 1) {
      post_burnin_samples[[i]] <-
        this_trace[(burnin + 1):nrow(this_trace),]
    } else if (burnin < 1 & burnin > 0) {
      discard <- ceiling(burnin * nrow(this_trace))
      post_burnin_samples[[i]] <-
        this_trace[(discard + 1):nrow(this_trace),]
    } else if (burnin == 0) {
      post_burnin_samples[[i]] <- this_trace
    } else {
      stop("What have you done?")
    }

  }

  # concatenate the samples
  output <- do.call(rbind, post_burnin_samples)

  # re-index the generations
  rownames(output) <- seq_len(nrow(output))
  if ("Iteration" %in% all_headers) {
    output$Iteration <- seq_len(nrow(output))
  }

  # return
  return(list("combined" = output))

}

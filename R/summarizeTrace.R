#' Summarize trace
#'
#' Summarizes trace file(s) that have been read into memory
#'
#' Summarizes a trace file for continuous or discrete characters by
#' computing the mean and 95\% credible interval for quantitative
#' character and the 95\% credible set for discrete characters.
#'
#' @param trace (list of data frames; no default) Name of a list of data frames,
#' such as produced by readTrace(). If the readTrace() output
#' contains multiple traces (such as from multiple runs), summarizeTrace() will
#' provide summaries for each trace individually, as well as the combined trace.
#'
#' @param vars (character or character vector; no default) The name of the
#' variable(s) to be summarized.
#'
#'
#' @return summarizeTrace() returns a list of the length of provided variables.
#' For quantitative variables, it returns the mean and 95% credible interval.
#' For discrete variables, it returns the 95% credible set of states and their
#' associated probabilities.
#'
#' @examples
#'
#' \donttest{
#' # continuous character only example, one run
#'
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
#' trace_sum <- summarizeTrace(trace = one_trace,
#'                             vars = c("pi[1]","pi[2]","pi[3]","pi[4]"))
#' trace_sum[["pi[1]"]]
#'
#' # remove file
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_gtr)
#'
#' # continuous character example, multiple runs
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
#' trace_sum_multi <- summarizeTrace(trace = multi_trace,
#'                                   vars = c("pi[1]","pi[2]","pi[3]","pi[4]"))
#' trace_sum_multi[["pi[1]"]]
#'
#' # remove files
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_1, dest_path_2)
#'
#' # discrete character example
#'
#' # download the example dataset to working directory
#' url_rj <- "https://revbayes.github.io/tutorials/intro/data/freeK_RJ.log"
#' dest_path_rj <- "freeK_RJ.log"
#' download.file(url_rj, dest_path_rj)
#'
#' file <- dest_path_rj
#' trace <- readTrace(path = file)
#'
#' trace_sum_discrete <- summarizeTrace(trace = trace,
#'                                      vars = c("prob_rate_12",
#'                                               "prob_rate_13",
#'                                               "prob_rate_31",
#'                                               "prob_rate_32"))
#' trace_sum_discrete[["prob_rate_12"]]
#'
#' #' # remove file
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_rj)
#' }
#'
#' @export

summarizeTrace <- function(trace, vars) {
  # enforce argument matching
  if (is.list(trace) == FALSE)
    stop("trace should be a list of data frames")
  if (is.data.frame(trace[[1]]) == FALSE)
    stop("trace should be a list of data frames")
  if (is.character(vars) == FALSE)
    stop("vars should be a character vector")

  # ensure variable names present in data frame
  if (any(vars %in% colnames(trace[[1]]) == FALSE) == TRUE) {
    stop(
      paste0(
      "The following variables you provided are not present in trace file:",
      paste0("\t", vars[!vars %in% colnames(trace[[1]])]),
      sep = "\n"
    )
    )
  }

  # subset to desired characters
  output <- list()
  # pass through vars and summarize each
  for (i in seq_len(length(vars))) {
    output[[i]] <- list()
    for (j in seq_len(length(trace))) {
      col <- trace[[j]][, vars[i]]
      if (class(col) == "numeric") {
        q_2.5 <- quantile(col, prob = c(0.025, 0.975))[1]
        q_97.5 <- quantile(col, prob = c(0.025, 0.975))[2]
        names(q_2.5) <- NULL
        names(q_97.5) <- NULL
        output[[i]][[j]] <- c(
          mean = mean(col),
          median = stats::median(col),
          MAP = getMAP(col),
          quantile_2.5 = q_2.5,
          quantile_97.5 = q_97.5
        )
        if (is.null(names(trace)[[j]]) == FALSE) {
          names(output[[i]])[j] <- names(trace)[[j]]
        } else {
          names(output[[i]])[j] <- paste0("trace_", j)
        }
      } else if (class(col) == "integer" |
                 class(col) == "character") {
        credible_set <- col
        state_probs <-
          sort(table(credible_set) / length(credible_set),
               decreasing = TRUE)
        cred_set <-
          state_probs[1:min(which((cumsum(
            state_probs
          ) >= 0.95) == TRUE))]
        output[[i]][[j]] <- cred_set
        if (is.null(names(trace)[[j]]) == FALSE) {
          names(output[[i]])[j] <- names(trace)[[j]]
        } else {
          names(output[[i]])[j] <- paste0("trace_", j)
        }

      }
    }
    #for multiple traces, combine and then summarize
    if (length(trace) > 1) {
      combined_trace <- do.call("rbind", trace)
      col <- combined_trace[, vars[i]]
      num_traces <- length(trace)
      if (class(col) == "numeric") {
        q_2.5 <- quantile(col, prob = c(0.025, 0.975))[1]
        q_97.5 <- quantile(col, prob = c(0.025, 0.975))[2]
        names(q_2.5) <- NULL
        names(q_97.5) <- NULL
        output[[i]][[num_traces + 1]] <- c(
          mean = mean(col),
          median = stats::median(col),
          MAP = getMAP(col),
          quantile_2.5 = q_2.5,
          quantile_97.5 = q_97.5
        )
        names(output[[i]])[num_traces + 1] <- "Combined"
      } else if (class(col) == "integer" |
                 class(col) == "character") {
        credible_set <- col
        state_probs <-
          sort(table(credible_set) / length(credible_set),
               decreasing = TRUE)
        cred_set <-
          state_probs[1:min(which((cumsum(
            state_probs
          ) >= 0.95) == TRUE))]
        output[[i]][[num_traces + 1]] <- cred_set
        names(output[[i]])[num_traces + 1] <- "Combined"
      }
    }
  }
  names(output) <- vars
  return(output)
}

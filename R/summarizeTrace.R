#' Summarize trace
#'
#' Summarizes trace file(s) that have been read into memory
#'
#' Summarizes a trace file for continuous or discrete characters by
#' computing the mean and 95% credible interval for quantitative
#' character and the 95% credible set for discrete characters.
#'
#' @param trace (list of data frames; no default) Name of a list of data frames,
#' such as produced by readTrace(). If the readTrace() output
#' contains multiple traces (such as from multiple runs), summarizeTrace() will provide
#' summaries for each trace individually, as well as the combined trace.
#'
#' @param vars (character or character vector; no default) The name of the variable(s)
#' to be summarized.
#'
#'
#' @return summarizeTrace() returns a list of the length of provided variables. For quantitative
#' variables, it returns the mean and 95% credible interval. For discrete variables, it returns
#' the 95% credible set of states and their associated probabilities.
#'
#' @examples
#'
#' \dontrun{
#' # continuous character only example, one run
#' file <- system.file("extdata",
#'     "sub_models/primates_cytb_covariotide.p", package="RevGadgets")
#' one_trace <- readTrace(path = file)
#' trace_sum <- summarizeTrace(trace = one_trace,
#'                             vars = c("Likelihood",
#'                                      "lambda",
#'                                    "TL"))
#'
#' # discrete character example, multiple runs
#' file1 <- system.file("extdata",
#'     "comp_method_disc/mkstates_run_1.txt", package="RevGadgets")
#' file2 <- system.file("extdata",
#'     "comp_method_disc/mkstates_run_2.txt", package="RevGadgets")
#'
#' multi_trace <- readTrace(path = c(file1, file2))
#' trace_sum_multi <- summarizeTrace(trace = multi_trace,
#'                                   vars = c("end_680",
#'                                            "end_681",
#'                                            "end_682",
#'                                            "end_683"))

#'
#' }
#'
#' @export

summarizeTrace <- function(trace, vars) {

  # enforce argument matching
  if (is.list(trace) == FALSE) stop("trace should be a list of data frames")
  if (is.data.frame(trace[[1]]) == FALSE) stop("trace should be a list of data frames")
  if (is.character(vars) == FALSE) stop("vars should be a character vector")

  # ensure variable names present in data frame
  if (any(vars %in% colnames(trace[[1]]) == FALSE) == TRUE) {
    cat("The following variables you provided are not present in trace file:",
        paste0("\t", vars[!vars %in% colnames(trace[[1]])]), sep = "\n")
    stop("oops!")
  }

  # subset to desired characters
  output <- list()

  # pass through vars and summarize each
  for (i in 1:length(vars)) {
    output[[i]] <-list()
    for (j in 1:length(trace)) {
    col <- trace[[j]][,vars[i]]
    if (class(col) == "numeric"){
      output[[i]][[j]] <- data.frame(mean = mean(col),
                                quantile_2.5 = quantile(col, prob = c(0.025,0.975))[1],
                                quantile_97.5 = quantile(col, prob = c(0.025,0.975))[2])
      names(output[[i]])[j] <- paste0("trace_", j)
    } else if (class(col) == "integer" | class(col) == "character" ){
      state_probs <- sort(table(col)/length(col), decreasing = TRUE)
      cred_set <- state_probs[1:min(which((cumsum(state_probs) >= 0.95) == TRUE))]
      output[[i]][[j]] <- cred_set
      names(output[[i]])[j] <- paste0("trace_", j)

    }
    }
    #for multiple traces, combine and then summarize
    if (length(trace) > 1) {
      combined_trace <- do.call("rbind", trace)
      col <- combined_trace[,vars[i]]
      num_traces <- length(trace)
      if (class(col) == "numeric"){
        output[[i]][[num_traces + 1]] <- data.frame(mean = mean(col),
                                       quantile_2.5 = quantile(col, prob = c(0.025,0.975))[1],
                                       quantile_97.5 = quantile(col, prob = c(0.025,0.975))[2])
        names(output[[i]])[num_traces + 1] <- "Combined"
      } else if (class(col) == "integer" | class(col) == "character" ){
        state_probs <- sort(table(col)/length(col), decreasing = TRUE)
        cred_set <- state_probs[1:min(which((cumsum(state_probs) >= 0.95) == TRUE))]
        output[[i]][[num_traces + 1]] <- cred_set
        names(output[[i]])[num_traces + 1] <- "Combined"
      }
      }
  }
  names(output) <- vars
  return(output)
}

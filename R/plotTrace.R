#' Plot trace
#'
#' Plots the posterior distributions of variables from trace file.
#'
#' Plots the posterior distributions of continuous variables from one or
#' multiple traces (as in, from multiple runs). If multiple traces are provided,
#' plotTrace() will plot each run independently as well as plot the combined output.
#' When multiple variables are listed, the default behavior is to overlay
#' the variables, unless the user specifies combine = FALSE. Note that for variables
#' with very different distributions, overlaying the plots may result in illegible
#' figures. In these cases, we recommend combine = FALSE.
#'
#'
#'
#' @param trace (list of data frames; no default) Name of a list of data frames,
#' such as produced by readTrace(). If the readTrace() output
#' contains multiple traces (such as from multiple runs), summarizeTrace() will provide
#' summaries for each trace individually, as well as the combined trace.
#'
#' @param vars (character or character vector; NULL) The specific name(s) of the variable(s)
#' to be summarized.
#'
#' @param match (character; NULL) Plot all variables together?
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
#' plots <- plotTrace(trace = one_trace,
#'                     vars = c("pi[1]","pi[2]","pi[3]","pi[4]"))
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


plotTrace <- function(trace, vars = NULL, match = NULL) {

  # enforce argument matching
  if (is.list(trace) == FALSE) stop("trace should be a list of data frames")
  if (is.data.frame(trace[[1]]) == FALSE) stop("trace should be a list of data frames")
  if (is.null(vars) & is.null(match)) stop("Either vars or match should be specified")
  if (!is.null(vars) & !is.null(match)) stop("Only vars OR match should be specified")
  if (is.character(vars) == FALSE | !is.null(vars)) stop("vars should be a character vector")
  if (is.character(match) == FALSE | !is.null(match)) stop("vars should be a character vector")

  # ensure variable names present in data frame
  if (any(vars %in% colnames(trace[[1]]) == FALSE) == TRUE) {
    cat("The following variables you provided are not present in trace file:",
        paste0("\t", vars[!vars %in% colnames(trace[[1]])]), sep = "\n")
    stop("oops!")
  }

    plots <- list()
    for(i in 1:length(trace)){
      if (length(vars) > 1) {
        trace[[i]][,vars] %>%
          melt() %>%
          ggplot(aes(x = value, fill = variable)) +
          geom_density(alpha = 0.5) +
          scale_fill_manual(values = colFun(length(vars))) +
          #scale_fill_brewer(palette = "Set1") +
          theme_few() +
          ggtitle(label = paste("Trace",i, sep = " ") ) +
          theme(plot.title = element_text(hjust = 0.5)) ->
          plots[[i]]
      } else if (length(vars) == 1) {
        data.frame(value = trace[[i]][,vars]) %>%
          ggplot(aes(x = value)) +
          geom_density(fill = "#e41a1c") +
          theme_few() +
          ggtitle(label = paste0("Trace ",i,": ",vars)) +
          theme(plot.title = element_text(hjust = 0.5)) ->
          plots[[i]]
      }
    }
    # add combined trace plots if multiple traces provided
    if (length(trace) > 1){
      if (length(vars) > 1) {
        do.call("rbind", trace)[,vars] %>%
          melt() %>%
          ggplot(aes(x = value, fill = variable)) +
          geom_density(alpha = 0.5) +
          #scale_fill_brewer(palette = "Set1") +
          scale_fill_manual(values = colFun(length(vars))) +
          theme_few() +
          ggtitle(label = "Combined Trace:") +
          theme(plot.title = element_text(hjust = 0.5)) ->
          plots[[length(trace) + 1]]
      } else if (length(vars) == 1) {
        data.frame(value =  do.call("rbind", trace)[,vars]) %>%
          ggplot(aes(x = value)) +
          geom_density(fill = "#e41a1c") +
          theme_few() +
          ggtitle(label = paste("Combined Trace:", vars)) +
          theme(plot.title = element_text(hjust = 0.5)) ->
          plots[[length(trace) + 1]]
      }
    }

  return(plots)
}

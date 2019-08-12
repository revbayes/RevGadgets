#' Summarize trace
#'
#' Summarizes trace file(s) that have been read into memory
#'
#' Summarizes a trace file for continuous or discrete characters by
#' computing the mean and 95% credible interval for quantitative
#' character and the probability for each state of discrete characters.
#'
#' @param trace (data frame; no default) Name of a single data frame,
#' such as one data frame produced by readTrace(). If the readTrace() output
#' contains multiple traces (such as from multiple runs), you should
#' combine them into one data frame prior to summarizing or summarize each
#' trace separately.
#' @param quant_char (character vector; default = "") Vector of the names of
#' quantitative characters to summarize.
#' @param disc_char (character vector; default = "") Vector of the names  of
#' discrete characters to summarize.
#'
#'
#'
#' @return List of data frames of length two. The first element
#' contains the summaries of quantitative variables (if no quantitative
#' variables included, element is empty). The second element contains
#' the summaries of discrete variables (again, empty if none provided).
#'
#' @examples
#'
#' \dontrun{
#' file <- system.file("extdata",
#'     "sub_models/primates_cytb_covariotide.p", package="RevGadgets")
#' one_trace <- readTrace(path = file)
#'
#' # continuous character example
#' trace_sum <- summarizeTrace(trace = one_trace[[1]]
#'                             quant_char = c("er[1]","er[2]",
#'                                            "er[3]","er[4]",
#'                                            "er[5]","er[6]"))
#' }
#'
#' @export

summarizeTrace <- function(trace, quant_char = "", disc_char = "") {

  # enforce argument matching
  if (is.data.frame(trace) == FALSE) stop("trace should be a data frame")
  if (is.character(quant_char) == FALSE) stop("quant_char should be a character vector")
  if (is.character(disc_char) == FALSE) stop("disc_char should be a character vector")

  # ensure variable names present in data frame
  if (sum(quant_char != "") > 0) {
  if (any(quant_char %in% colnames(trace) == FALSE) == TRUE) {
    cat("The following quantitative characters you provided are not present in trace file:",
        paste0("\t", quant_char[!quant_char %in% colnames(trace)]), sep = "\n")
    stop("oops!")
  }
  }
  if (sum(disc_char != "") > 0) {
    if (any(disc_char %in% colnames(trace) == FALSE) == TRUE) {
      cat("The following discrete characters you provided are not present in trace file:",
          paste0("\t", disc_char[!disc_char %in% colnames(trace)]), sep = "\n")
      stop("oops!")
    }
  }
  output <- list()
  # process quantitative characters
  if (sum(quant_char != "") > 0){
      df_sub <- trace[,quant_char]
      means <- colMeans(df_sub)
      quantiles <- apply(df_sub, 2, quantile, prob = c(0.025,0.975))
      res_quant<- rbind(means,quantiles)
      output[[1]] <- res_quant
  } else if (sum(quant_char != "") == 0) { output[[1]] <- NA }

  # process discrete characters
  if (sum(disc_char != "") > 0){
    ngen <- nrow(trace)
    df_sub <- trace[,disc_char]
    for (i in 1:ncol(df_sub)) {
      t <- table(df_sub[,i])/ngen
    }

  } else if (sum(disc_char != "") == 0) { output[[2]] <- NA }

  names(output) <- c("quantitative", "discrete")
  return(output)
}

#' Title
#'
#' @param path (vector of character strings; no default) File path(s) to
#' trace file.
#' @param speciation (single character string; "speciation") RevBayes variable
#' name
#' @param extinction (single character string; "extinction") RevBayes variable
#' name
#' @param speciation_hidden (single character string; "speciation_hidden")
#' RevBayes variable name
#' @param rates (vector; c(speciation, extinction, "net-diversification"))
#' names of rates to be included in plot
#' @param ... additional arguments passed to readTrace()
#'
#' @return a data frame
#' @examples
#'
#' \donttest{
#' # download the example dataset to working directory
#'
#' url <-
#'   "https://revbayes.github.io/tutorials/intro/data/primates_BiSSE_activity_period.log"
#' dest_path <- "primates_BiSSE_activity_period.log"
#' download.file(url, dest_path)
#'
#' # to run on your own data, change this to the path to your data file
#' bisse_file <- dest_path
#'
#' pdata <- processSSE(bisse_file)
#'
#' # remove file
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path)
#' }
#'
#' @export
processSSE <- function(path,
                       speciation = "speciation",
                       extinction = "extinction",
                       speciation_hidden = "speciation_hidden",
                       rates = c(speciation, extinction, "net-diversification"),
                       ...) {
  # parameter compatibility checks
  if (!is.character(speciation))
    stop("speciation should be a single character string")
  if (length(speciation) != 1)
    stop("speciation should be a single character string")
  if (!is.character(extinction))
    stop("extinction should be a single character string")
  if (length(extinction) != 1)
    stop("extinction should be a single character string")
  if (!is.character(speciation_hidden))
    stop("speciation_hidden should be a single character string")
  if (length(speciation_hidden) != 1)
    stop("speciation_hidden should be a single character string")

  # read in trace
  tr <- readTrace(paths = path, ...)[[1]]

  # process trace
  n_hidden <-
    max(1, sum(grepl(
      paste0(speciation_hidden, "\\["), names(tr)
    )))
  n_states <-
    sum(grepl(paste0(speciation, "\\["), names(tr))) / n_hidden
  n_rates <- n_hidden * n_states

  for (index in 1:n_rates) {
    netdiv <-
      tr[[paste0(speciation,
                 "[",
                 index,
                 "]")]] - tr[[paste0(extinction,
                                     "[",
                                     index,
                                     "]")]]
    tr[[paste0("net-diversification[", index, "]")]] <- netdiv
  }

  dfs <- list()
  m <- 1
  for (k in seq_along(rates)) {
    for (i in 1:n_states) {
      for (j in 1:n_hidden) {
        hiddenletter <- LETTERS[j]
        index <- i + (j * n_states) - n_states

        value <- tr[[paste0(rates[k], "[", index, "]")]]

        df1 <- data.frame(
          "value" = value,
          "rate" = rates[k],
          "hidden_state" = hiddenletter,
          "label" = paste0(i - 1, hiddenletter),
          "observed_state" = as.factor(i - 1),
          "Iteration" = tr$Iteration
        )


        dfs[[m]] <- df1
        m <- m + 1
      }
    }
  }
  res <- dplyr::bind_rows(dfs)
  return(res)
}

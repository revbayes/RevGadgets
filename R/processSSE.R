#' Title
#'
#' @param tr a trace object
#' @param speciation revbayes variable name
#' @param extinction revbayes variable name
#' @param speciation_hidden revbayes variable name
#' @param rates names of rates to be included in plot
#'
#' @return a data frame
#' @examples bisse_file <- system.file("extdata", "sse/primates_BiSSE_activity_period.log", package="RevGadgets")
#'
#' tr <- readTrace(bisse_file)[[1]]
#' pdata <- processSSE(tr)
#' head(pdata)
#'
#' @export
processSSE <- function(tr,
                       speciation = "speciation",
                       extinction = "extinction",
                       speciation_hidden = "speciation_hidden",
                       rates = c(speciation, extinction, "net-diversification")
){
  n_hidden <- max(1, sum(grepl(paste0(speciation_hidden, "\\["), names(tr))))
  n_states <- sum(grepl(paste0(speciation, "\\["), names(tr))) / n_hidden
  n_rates <- n_hidden*n_states

  for (index in 1:n_rates){
    netdiv <- tr[[paste0(speciation, "[", index, "]")]] - tr[[paste0(extinction, "[", index, "]")]]
    tr[[paste0("net-diversification[",index,"]")]] <- netdiv
  }

  dfs <- list(); m <- 1
  for (k in seq_along(rates)){
    for (i in 1:n_states){
      for (j in 1:n_hidden){
        hiddenletter <- LETTERS[j]
        index <- i + (j*n_states) - n_states

        value <- tr[[paste0(rates[k], "[", index, "]")]]

        df1 <- data.frame("value" = value,
                          "rate" = rates[k],
                          "hidden_state" = hiddenletter,
                          "label" = paste0(i-1, hiddenletter),
                          "observed_state" = as.factor(i-1),
                          "Iteration" = tr$Iteration)


        dfs[[m]] <- df1
        m <- m + 1
      }
    }
  }
  res <- dplyr::bind_rows(dfs)
  return(res)
}


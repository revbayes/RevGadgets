#' plotMuSSE
#'
#' @param rates (data.frame; no default) a data frame containing columns "value", "rate",
#' "hidden_state", "observed_state" (such as the output of processSSE())
#'
#' @return a ggplot object
#' @examples
#' \dontrun{
#' bisse_file <- system.file("extdata", "sse/primates_BiSSE_activity_period.log", package="RevGadgets")
#'
#' pdata <- processSSE(bisse_file)
#' p <- plotMuSSE(pdata);p
#'
#' # change colors:
#' p + scale_fill_manual(values = c("red","green"))
#'
#' # change x-axis label
#' library(ggplot2)
#' p + xlab("Rate (events/Ma)")
#' }
#' @export

plotMuSSE <- function(rates){
  if (is.data.frame(rates) == FALSE) stop("rates should be a data frame")
  p <- ggplot2::ggplot(rates, ggplot2::aes(x = value, fill = observed_state)) +
       ggplot2::geom_density(alpha=0.8) +
       ggplot2::scale_fill_manual(values = colFun(length(unique(rates$observed_state))),
                                  name = "Observed state") +
       ggplot2::facet_wrap( ~ rate,
                            scales = "free",
                            ncol = 1,
                            labeller = ggplot2::labeller(rate = .titleFormatLabeller)) +
       ggplot2::xlab("Rate") +
       ggplot2::ylab("Posterior density") +
       ggplot2::theme_bw() +
       ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                      panel.grid.minor = ggplot2::element_blank(),
                      strip.background = ggplot2::element_blank())
  return(p)
}

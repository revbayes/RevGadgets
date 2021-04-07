#' plotHiSSE
#'
#' @inheritParams plotMuSSE
#'
#' @return a ggplot object
#' @examples
#' \dontrun{
#' hisse_file <- system.file("extdata", "sse/primates_HiSSE_2.p", package="RevGadgets")
#' pdata <- processSSE(hisse_file)
#' p <- plotHiSSE(pdata);p
#'
#' # change colors:
#' library(ggplot2)
#' p + scale_fill_manual(values = c("red","green"))
#'
#' # change x-axis label
#' library(ggplot2)
#' p + xlab("Rate (events/Ma)")
#' }
#' @export
plotHiSSE <- function(rates){
  if (is.data.frame(rates) == FALSE) stop("rates should be a data frame")
  p <- ggplot2::ggplot(rates, ggplot2::aes(x = value, fill = observed_state)) +
    ggplot2::geom_density(alpha=0.8) +
    ggplot2::facet_grid(rate ~ hidden_state,
                        scales = "free",
                        labeller = ggplot2::labeller(rate = .titleFormatLabeller)) +
    ggplot2::scale_fill_manual(values = colFun(length(unique(rates$observed_state))),
                               name = "Observed state") +
    ggplot2::xlab("Rate") +
    ggplot2::ylab("Posterior density") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank())
  return(p)
}

#' plotMuSSE
#'
#' @param x a data frame containing columns "value", "rate", "hidden_state", "observed_state"
#'
#' @return a ggplot object#'
#' @examples bisse_file <- system.file("extdata", "sse/primates_BiSSE_activity_period.log", package="RevGadgets")
#'
#' tr <- readTrace(bisse_file)[[1]]
#'
#' pdata <- processSSE(tr)
#' p <- plotMuSSE(pdata)
#'
#' @export
plotMuSSE <- function(x){
  p <- ggplot2::ggplot(x, ggplot2::aes(x = value, fill = observed_state)) +
    ggplot2::geom_density(alpha=0.8) +
    ggplot2::facet_wrap( ~ rate,
                         scales = "free",
                         ncol = 1) +
    ggplot2::xlab("Rate (events/time)") +
    ggplot2::ylab("Posterior density") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank())
  return(p)
}

#' plotHiSSE
#'
#' @inheritParams plotMuSSE
#'
#' @return a ggplot object
#' @examples hisse_file <- system.file("extdata", "sse/primates_HiSSE_2.log", package="RevGadgets")
#'
#' tr <- readTrace(hisse_file)[[1]]
#'
#' pdata <- processSSE(tr)
#' p <- plotHiSSE(pdata)
#' @export
plotHiSSE <- function(x){
  p <- ggplot2::ggplot(x, ggplot2::aes(x = value, fill = observed_state)) +
    ggplot2::geom_density(alpha=0.8) +
    ggplot2::facet_grid(rate ~ hidden_state,
                        scales = "free") +
    ggplot2::xlab("Rate (events/time)") +
    ggplot2::ylab("Posterior density") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank())
  return(p)
}

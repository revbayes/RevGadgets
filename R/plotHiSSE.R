#' plotHiSSE
#'
#' @inheritParams plotMuSSE
#'
#' @return a ggplot object
#' @examples
#' \donttest{
#' # download the example dataset to working directory
#'
#' url <- "https://revbayes.github.io/tutorials/intro/data/primates_HiSSE_2.log"
#' dest_path <- "primates_HiSSE_2.log"
#' download.file(url, dest_path)
#'
#' # to run on your own data, change this to the path to your data file
#' hisse_file <- dest_path
#'
#' pdata <- processSSE(hisse_file)
#' p <- plotHiSSE(pdata);p
#'
#' # change colors:
#' p + ggplot2::scale_fill_manual(values = c("red","green"))
#'
#' # change x-axis label
#' p + ggplot2::xlab("Rate (events/Ma)")
#'
#' # remove file
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path)
#'
#' }
#' @export
plotHiSSE <- function(rates) {
  if (is.data.frame(rates) == FALSE)
    stop("rates should be a data frame")
  p <-
    ggplot2::ggplot(rates, ggplot2::aes(x = value, fill = observed_state)) +
    ggplot2::geom_density(alpha = 0.8) +
    ggplot2::facet_grid(
      rate ~ hidden_state,
      scales = "free",
      labeller = ggplot2::labeller(rate = .titleFormatLabeller)
    ) +
    ggplot2::scale_fill_manual(values = colFun(length(unique(
      rates$observed_state
    ))),
    name = "Observed state") +
    ggplot2::xlab("Rate") +
    ggplot2::ylab("Posterior density") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    )
  return(p)
}

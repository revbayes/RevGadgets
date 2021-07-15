#' Plot Diversification Rates
#'
#' Plots the output of an episodic diversification rate analysis
#'
#' Plots the output of episodic diversification rate analyses. Takes as
#' input the output of processDivRates() and plotting parameters.
#' For now, only variable names (under "item") that contain the word "rate" are
#' included in the plot.
#'
#' The return object can be manipulated. For example, you can change the
#' axis labels, the color palette, whether the axes are to be linked, or the
#' overall plotting style/theme, just as with any ggplot object.
#'
#'
#' @param rates (list of dataframes; no default) A list of dataframes,
#' such as produced by processDivRates(), containing the data on rates
#' and interval times for each type of rate to be plotted (e.g.
#' speciation rate, etc.).
#'
#' @param facet (logical; TRUE) plot rates in separate facets.
#'
#'
#' @return A ggplot object
#'
#' @examples
#'
#' \dontrun{
#' # first run processDivRates()
#' speciation_time_file <- system.file("extdata",
#'     "epi_bd/primates_EBD_speciation_times.p", package="RevGadgets")
#' speciation_rate_file <- system.file("extdata",
#'     "epi_bd/primates_EBD_speciation_rates.p", package="RevGadgets")
#' extinction_time_file <- system.file("extdata",
#'     "epi_bd/primates_EBD_extinction_times.p", package="RevGadgets")
#' extinction_rate_file <- system.file("extdata",
#'     "epi_bd/primates_EBD_extinction_rates.p", package="RevGadgets")
#'
#' rates <- processDivRates(speciation_time_log = speciation_time_file,
#'                          speciation_rate_log = speciation_rate_file,
#'                          extinction_time_log = extinction_time_file,
#'                          extinction_rate_log = extinction_rate_file,
#'                          burnin = 0.25)
#'
#' # then plot results:
#' p <- plotDivRates(rates = rates);p
#'
#' # change the x-axis
#' library(ggplot2)
#' p <- p + xlab("Thousands of years ago");p
#'
#' # change the colors
#' library(ggplot2)
#' p <- p + scale_fill_manual(values = c("red", "green", "yellow", "purple")) +
#'   scale_color_manual(values = c("red", "green", "yellow", "purple"));p
#'
#' # let's say we don't want to plot relative-extinction rate,
#' # and use the same y-axis for all three rates
#' library(ggplot2)
#' rates[grep("relative-extinction", names(rates))] <- NULL
#' p2 <- plotDivRates(rates)
#' p2 <- p2 + facet_wrap(vars(item), scale = "fixed");p2
#' }
#'
#' @export

plotDivRates <- function(rates, facet = TRUE){
  message("Using default time units in x-axis label: Age (Ma)")
  rates_to_plot <- unique(rates$item)[grep("rate", unique(rates$item))]
  `%>%` <- dplyr::`%>%`

    p <- rates %>%
    subset(grepl("rate", item)) %>%
    ggplot2::ggplot(ggplot2::aes(time, value, color = item))  +
    ggplot2::geom_step(ggplot2::aes(time, value),
                       direction = "vh") +
    geom_stepribbon(ggplot2::aes(x = time,
                                 ymin = lower,
                                 ymax = upper,
                                 fill = item),
                    direction = "vh",
                    alpha = 0.4,
                    color = NA) +
    ggplot2::scale_x_reverse() +
    ggplot2::xlab("Age (Ma)") +
    ggplot2::ylab("Rate") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = colFun(length(rates_to_plot))) +
    ggplot2::scale_fill_manual(values = colFun(length(rates_to_plot)))

    if (facet){
      p <- p +
        ggplot2::facet_wrap(dplyr::vars(item),
                            scales = "free_y",
                            labeller =
                              ggplot2::labeller(item = .titleFormatLabeller))
    }


  return(p)
}

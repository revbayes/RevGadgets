#' Plot Population Sizes
#'
#' Plots the output of a coalescent demographic analysis.
#'
#' Plots the output of coalescent demographic analyses. Takes as
#' input the output of processPopSizes() and plotting parameters.
#'
#' The return object can be manipulated. For example, you can change the
#' axis labels, the color palette, whether the axes are to be linked, or the
#' overall plotting style/theme, just as with any ggplot object.
#'
#' @param df (data frame) such as produced by processPopSizes(), containing
#' the data on population sizes and corresponding grid points (points in time for population size evaluation)
#' @param plot_CIs (boolean; default: TRUE) specifies whether the credible intervals should be plotted.
#' @param add (boolean; default: FALSE) specifies whether the new plot should be added to an existing ggplot2 object. If TRUE,
#' the existing_plot has to be given.
#' @param existing_plot (ggplot2 object; default: NULL) a ggplot2 object to which the new plot should be added.
#' @param col (string; default: "#00883a") color for the trajectories
#'
#'
#' @return a ggplot object
#'
#' @examples
#' df <- dplyr::tibble("time" = c(0.0, 1.0, 2.0, 3.0, 4.0),
#'                     "value" = c(1.0, 1.5, 2.0, 1.5, 1.5),
#'                     "upper" = c(3.5, 7.0, 6.5, 5.0, 5.0),
#'                     "lower" = c(0.5, 0.1, 0.5, 0.5, 0.8))
#'
#' plotPopSizes(df)
#'
#' @export

plotPopSizes <- function(df,
                         plot_CIs = TRUE,
                         add = FALSE,
                         existing_plot = NULL,
                         col = "#00883a"){
  if (add == TRUE && is.null(existing_plot)){
    stop("Please provide an existing plot if you want to add this one.")
  }

  `%>%` <- dplyr::`%>%`

  if (add == FALSE){
    message("Using default time units in x-axis label: Age (years)")

    p <- df %>%
      ggplot2::ggplot(ggplot2::aes(x = time, y = value)) +
      ggplot2::geom_line(color = col, linewidth = 0.8) +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_reverse() +
      ggplot2::xlab("Age (years)") +
      ggplot2::ylab("Population Size")

  } else {
    p <- existing_plot +
      ggplot2::geom_line(data = df, 
                         ggplot2::aes(x = time, y = value), 
                         color = col, 
                         linewidth = 0.8)
  }

  if (plot_CIs == TRUE){
    p <- p + 
      ggplot2::geom_ribbon(data = df, 
                           ggplot2::aes(ymin = lower, ymax = upper), 
                           fill = col, 
                           alpha = 0.4)
  }

  p <- p +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank())

  return(p)
}

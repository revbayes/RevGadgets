#' Plot Population Sizes
#' 
#' Plots the output of a coalescent demographic analysis.
#'
#' Plots the output of coalescent demographic analyses. Takes as
#' input the output of processPopSizes() and plotting parameters.
#' For now, only variable names that contain the word "size" are
#' included in the plot.
#'
#' The return object can be manipulated. For example, you can change the
#' axis labels, the color palette, whether the axes are to be linked, or the
#' overall plotting style/theme, just as with any ggplot object.
#'
#' @param df (data frame) such as produced by processPopSizes(), containing 
#' the data on population sizes and, if applicable, interval times
#' @param plot_var (string, default: "size") only include variables that contain the string, default "size" in the name
#' @param add (boolean; default: TRUE) specifies whether the new plot should be added to an existing ggplot2 object. If TRUE,
#' the existing_plot has to be given.
#' @param existing_plot (ggplot2 object; default: NULL) a ggplot2 object to which the new plot should be added.
#' @param col (string; default: "#00883a") color for the trajectories
#'
#'
#' @return a ggplot object
#' 
#' @examples
#' library(tibble)
#' 
#' df <- tibble("time" = c(0.0, 1.0, 2.0, 3.0, 4.0),
#'              "time_end" = c(1.0, 2.0, 3.0, 4.0, 5.0),
#'              "value" = c(1.0, 1.5, 2.0, 1.5, 1.5),
#'              "upper" = c(3.5, 7.0, 6.5, 5.0, 5.0),
#'              "lower" = c(0.5, 0.0, 0.5, 0.5, 0.8),
#'              "item" = "population size")
#'              
#' plotPopSizes(df)
#'
#' @export
#' @importFrom ggplot2 aes ggplot theme xlab ylab theme_bw scale_color_manual scale_fill_manual scale_x_reverse labeller facet_wrap element_blank geom_segment geom_rect geom_line geom_ribbon scale_y_log10 coord_cartesian xlim
#' @importFrom dplyr bind_rows
#' @importFrom utils head tail
#' 
plotPopSizes <- function(df, plot_var = "size", add = FALSE, existing_plot = NULL, col = "#00883a"){
  if (add == TRUE && is.null(existing_plot)){
    stop("Please provide an existing plot if you want to add this one.")
  }
  
  `%>%` <- dplyr::`%>%`
  
  if (add == FALSE){
    message("Using default time units in x-axis label: Age (years)")
    
    p <- df %>%
      ggplot(aes(x = time, y = value)) +
      geom_line(color = col, size = 0.8) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = col, alpha = 0.4) +
      scale_y_log10() +
      scale_x_reverse() +
      xlab("Age (years)") +
      ylab("Population Size")
  } else {
    p <- existing_plot +
      geom_line(data = df, aes(x = time, y = value), color = col, size = 0.8) +
      geom_ribbon(data = df, aes(ymin = lower, ymax = upper), fill = col, alpha = 0.4)        
  }

  
  p <- p +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank())
  
  return(p)
}


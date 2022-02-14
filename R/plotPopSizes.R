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
#' @param df a data frame, such as produced by processPopSizes(), containing 
#' the data on population sizes and, if applicable, interval times
#' @param plot_var only include variables that contain the string, default "size" in the name
#' @param method which method was chosen for the analysis, "events" - event-based coalescent process,
#' "specified" - coalescent process with user-defined interval times, "constant" - constant coalescent process
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
#' plotPopSizes(df, method = "specified")
#'
#' @export
#' @importFrom ggplot2 aes ggplot theme xlab ylab theme_bw scale_color_manual scale_fill_manual scale_x_reverse labeller facet_wrap element_blank geom_segment geom_rect geom_line geom_ribbon scale_y_log10 coord_cartesian xlim
#' @importFrom dplyr bind_rows
#' @importFrom utils head tail
#' 
plotPopSizes <- function(df, plot_var = "size", method = "events"){
  message("Using default time units in x-axis label: Age (years)")
  `%>%` <- dplyr::`%>%`
  
  if (method == "events"){
    p <- df %>%
      ggplot(aes(x = time, y = value)) +
      geom_line(color = "#00883a", size = 1) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#00883a", alpha = 0.4)
    
    
  } else if (method == "specified"){
    vert_lines <- lapply(grep(plot_var, unique(df$item), value = TRUE),
                         function(item) .makeVerticalLines(df, item)) %>%
      bind_rows()
    
    pdata <- df %>%
      subset(grepl(plot_var, item)) 
    
    rates_to_plot <- unique(pdata$item)
    
    p <- pdata %>%
      ggplot(aes(x = time, 
                 y = value, 
                 yend = value, 
                 xend = time_end))  +
      geom_segment(color = "#00883a", size = 1) + ## plot horizontal segments
      geom_segment(data = vert_lines, 
                   aes(y = y, x = x, yend = yend, xend = xend), color = "#00883a", size = 1.2) + ## plot the vertical segments
      geom_rect(aes(xmin = time, xmax = time_end, ymin = lower, ymax = upper), fill = "#00883a",
                alpha = 0.4)
    
  } else if (method == "constant"){

    p <- df %>%
      ggplot(aes(x = time, 
                 y = value, 
                 yend = value, 
                 xend = time_end))  +
      geom_segment(color = "#00883a", size = 1) + ## plot horizontal segment
      geom_rect(aes(xmin = time, xmax = time_end, ymin = lower, ymax = upper), fill = "#00883a",
                alpha = 0.4) +
      xlim(1e5, 0)
    
  } else {
    stop("Please choose as method either 'events', 'specified' or 'constant'")
  }
  
  p <- p +
    scale_y_log10() +
    scale_x_reverse() +
    xlab("Age (years)") +
    ylab("Population Size") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank())
  
  return(p)
}
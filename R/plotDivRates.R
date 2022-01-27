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
#' \donttest{
#' # download the example datasets to working directory
#'
#' url_ex_times <-
#'    "https://revbayes.github.io/tutorials/intro/data/primates_EBD_extinction_times.log"
#' dest_path_ex_times <- "primates_EBD_extinction_times.log"
#' download.file(url_ex_times, dest_path_ex_times)
#'
#' url_ex_rates <-
#'    "https://revbayes.github.io/tutorials/intro/data/primates_EBD_extinction_rates.log"
#' dest_path_ex_rates <- "primates_EBD_extinction_rates.log"
#' download.file(url_ex_rates, dest_path_ex_rates)
#'
#' url_sp_times <-
#'    "https://revbayes.github.io/tutorials/intro/data/primates_EBD_speciation_times.log"
#' dest_path_sp_times <- "primates_EBD_speciation_times.log"
#' download.file(url_sp_times, dest_path_sp_times)
#'
#' url_sp_rates <-
#'    "https://revbayes.github.io/tutorials/intro/data/primates_EBD_speciation_rates.log"
#' dest_path_sp_rates <- "primates_EBD_speciation_rates.log"
#' download.file(url_sp_rates, dest_path_sp_rates)
#'
#' # to run on your own data, change this to the path to your data file
#' speciation_time_file <- dest_path_sp_times
#' speciation_rate_file <- dest_path_sp_rates
#' extinction_time_file <- dest_path_ex_times
#' extinction_rate_file <- dest_path_ex_rates
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
#' p <- p + xlab("Thousands of years ago");p
#'
#' # change the colors
#' p <- p + scale_fill_manual(values = c("red",
#'                                                "green",
#'                                                "yellow",
#'                                                "purple")) +
#'   scale_color_manual(values = c("red",
#'                                          "green",
#'                                          "yellow",
#'                                          "purple"));p
#'
#' # let's say we don't want to plot relative-extinction rate,
#' # and use the same y-axis for all three rates
#' rates <- rates[!grepl("relative-extinction", rates$item),]
#' p2 <- plotDivRates(rates)
#' p2 <- p2 + facet_wrap(vars(item), scale = "fixed");p2
#'
#' # remove files
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_sp_times, dest_path_ex_times,
#'             dest_path_sp_rates, dest_path_ex_rates)
#'
#' }
#'
#' @export
#' @importFrom ggplot2 aes ggplot theme xlab ylab theme_bw scale_color_manual scale_fill_manual scale_x_reverse labeller facet_wrap element_blank geom_segment
#' @importFrom dplyr bind_rows

plotDivRates <- function(rates, facet = TRUE){
  message("Using default time units in x-axis label: Age (Ma)")
  rates_to_plot <- unique(rates$item)[grep("rate", unique(rates$item))]
  `%>%` <- dplyr::`%>%`
  
  vert_lines <- lapply(grep("rate", unique(rates$item), value = TRUE),
                       function(item) make_vertical_lines(rates, item)) %>%
    bind_rows()
  
  p <- rates %>%
    subset(grepl("rate", item)) %>%
    ggplot(aes(x = time, 
               y = value, 
               yend = value, 
               xend = time_end))  +
    geom_segment(aes(color = item)) + ## plot horizontal segments
    geom_segment(data = vert_lines, 
                 aes(y = y, x = x, yend = yend, xend = xend, color = item)) + ## plot the vertical segments
    geom_rect(aes(xmin = time, xmax = time_end, ymin = lower, ymax = upper, fill = item),
              alpha = 0.4) +
    scale_x_reverse() +
    xlab("Age (Ma)") +
    ylab("Rate") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()) +
    scale_color_manual(values = colFun(length(rates_to_plot))) +
    scale_fill_manual(values = colFun(length(rates_to_plot)))
  
  if (facet){
    p <- p +
      facet_wrap(dplyr::vars(item),
                 scales = "free_y",
                 labeller =
                   labeller(item = .titleFormatLabeller))
  }
  
  
  return(p)
}

make_vertical_lines <- function(df, item1){
  d <- subset(df, item == item1)
  
  res <- tibble::tibble(y = head(d$value, n = -1), 
                        yend = tail(d$value, n = -1), 
                        x =  tail(d$time, n = -1), 
                        xend =  tail(d$time, n = -1))
  res$item <- item1
  
  return(res)
}
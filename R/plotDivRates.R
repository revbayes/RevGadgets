#' Plot Diversification Rates
#'
#' Plots the output of a episodic diversification rate analysis
#'
#' Plots the output of episodic diversification rate analyses. Takes as input the
#' output of processDivRates() and plotting parameters. Does not return an
#' object. For now, valid figure types include: "speciation rate",
#' "extinction rate","net-diversification rate", "relative-extinction rate",
#' and "fossilization rate".
#' If colors are not provided (parameter col), the plot defaults to preset colors.
#'
#'
#' @param rates (list; no default) The processed output for plotting (output of processDivRates()).
#' @param fig_types (character vector, c("speciation rate", "extinction rate", "net-diversification rate", "relative-extinction rate")) Which aspects of the model to visualize. "fossilization rate" may also be included but is not plotted by default. See details for a complete description.
#' @param xlab (character; "million years ago") The label of the x-axis.
#' @param col (character; NULL) Colors used for printing in hex code. Must be of same length as fig_types.
#' @param col_alpha (character; "50") Alpha channel parameter for credible intervals plotting. May range from 00 (completely transparent) to FF (completely opaque).
#' @param offset (numeric; 0) Value to offset the x-axis by.
#' @param xaxt (character; "n") The type of x-axis to plot. By default, no x-axis is plotted (recommended).
#' @param yaxt (character; "s") The type of y-axis to plot.
#' @param ... Parameters delegated to various plotting functions.
#'
#' @return Plots diversification rates, does not return an object.
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
#' primates <- processDivRates(speciation_time_log = speciation_time_file,
#'                             speciation_rate_log = speciation_rate_file,
#'                             extinction_time_log = extinction_time_file,
#'                             extinction_rate_log = extinction_rate_file,
#'                             burnin = 0.25)
#' # then plot results:
#' plotDivRates(rates = primates)
#' }
#'
#' @export

plotDivRates <- function(rates,
                         fig_types = c("speciation rate",
                                       "extinction rate",
                                       "net-diversification rate",
                                       "relative-extinction rate"),
                         xlab = "million years ago",
                         col = NULL,
                         col_alpha = "50",
                         offset = 0,
                         xaxt = "n",
                         yaxt = "s",
                         ...){
  # Enforce argument matching
  if (is.list(rates) == FALSE) stop("rates must be list of processed rates")
  if (is.character(fig_types) == FALSE) stop("fig_types must be a character string or vector of strings")
  if (is.character(xlab) == FALSE) stop("xlab must be a single character string")

  # check length of color vector and ensure colors are formatted as hex codes
  if (is.null(col) == FALSE) {
    if (length(col) != length(fig_types))
      stop("col and fig_types must be of the same length")
    bad_cols <- rep("FALSE", times = length(col))
    for (i in 1:length(col)) {
      col_split <- unlist(strsplit(col[i], split = ""))
      if (length(col_split) != 7 | col_split[1] != "#") {
        bad_cols[i] <- TRUE
      }
    }
    if ( any(bad_cols == TRUE) == TRUE ) {
      # print out colors that are 'bad'
      cat( "Some colors are not properly formatted hex codes:",
           paste0("\t",col[bad_cols == TRUE]), sep="\n")
      stop()
    }
  }

  if (is.character(col_alpha) == FALSE) stop("col_alpha must be a character string")
  if (is.numeric(offset) == FALSE) stop("offset must be a numeric value")
  xaxt <- match.arg(xaxt, choices = c("n", "s"))
  yaxt <- match.arg(yaxt, choices = c("n", "s"))


  # Check that fig type is valid
  valid_fig_types <- c("speciation rate",
                       "extinction rate",
                       "net-diversification rate",
                       "relative-extinction rate",
                       "fossilization rate")

  invalid_fig_types <- fig_types[!fig_types %in% valid_fig_types]

  if ( length(invalid_fig_types) > 0 ) {
    stop("\nThe following figure types are invalid: ", paste(invalid_fig_types, collapse = ", "), ".",
         "\nValid options are: ", paste(valid_fig_types, collapse = ", ") ,".")
  }


  # Make color vector
  if ( is.null(col) ) {

    col <- .colFun(5)
    names(col) <- c("speciation rate",
                    "extinction rate",
                    "net-diversification rate",
                    "relative-extinction rate",
                    "fossilization rate")
  } else {
    names(col) <- fig_types
  }

  # Compute the axes
  intervals <- rates$`speciation time`[1,grep("interval", colnames(rates$`speciation time`))]
  tree_age <- max(intervals)
  num_intervals <- length(intervals) - 1
  plot_at <- 0:num_intervals
  interval_size <- tree_age/num_intervals
  labels <- pretty( c(0,tree_age) ) + offset
  labels_at <- num_intervals - ( pretty(c(0,tree_age)) / interval_size )
  ages <- seq( 0,tree_age, length.out = num_intervals + 1 )
  for( type in fig_types ) {

      these_rates <- rates[[type]][ ,grep("[0-9]", colnames(rates[[type]])) ]
      mean_these_rates <- colMeans(these_rates)
      quantiles_these_rates <- apply(these_rates, 2, quantile, prob = c(0.025,0.975))

      #find the limits of the y-axis
      #we always want the speciation and extinction rates to be on the same scale
      if ( type %in% c("speciation rate","extinction rate")) {

        quantiles_speciation <-
          apply(rates[["speciation rate"]][ ,grep("[0-9]", colnames(rates[["speciation rate"]])) ],
                2, quantile, prob=c(0.025,0.975))
        quantiles_extinction <-
          apply(rates[["extinction rate"]][ ,grep("[0-9]", colnames(rates[["extinction rate"]])) ],
                2, quantile, prob=c(0.025,0.975))
        ylim <- c( min(0, quantiles_speciation, quantiles_extinction),
                   max(quantiles_speciation, quantiles_extinction) )

      } else {
        ylim <- c( min(0, quantiles_these_rates),
                   max(quantiles_these_rates) )
      }

      #### plot
        plot( x = plot_at, y = c( rev(mean_these_rates)),
              type = "l", ylim = ylim, xaxt = xaxt, col = col[type],
              ylab = "rate", main = type, xlab = xlab,...)
        polygon( x = c( 0:ncol(quantiles_these_rates),
                        ncol(quantiles_these_rates):0 ),
                 y = c( rev(c(quantiles_these_rates[1,1],
                          quantiles_these_rates[1,])),
                        c(quantiles_these_rates[2,1],
                              quantiles_these_rates[2,])),
                 border = NA, col = paste(col[type], col_alpha, sep="") )
        axis(1, at = labels_at, labels = labels)
      }

    }


computeMeanInterval <- function(item, rates, probs){
  interval_times <- unlist(rates[["speciation time"]][1,grepl("interval_times", names(rates$`speciation time`))])
  interval_times <- sort(interval_times) # For some reason these are ordered differently than rate vectors

  rate <- rates[[item]]
  rate <- rate[, grep("[0-9]", colnames(rate))]

  mean_rate <- colMeans(rate)
  quantiles <- apply(rate, 2,
                     quantile,
                     probs = probs)

  df <- tibble(.rows = length(mean_rate))
  df["mean"] <- mean_rate
  df["lower"] <- quantiles[1,]
  df["upper"] <- quantiles[2,]
  df$time <- interval_times
  df$item <- item

  return(df)
}

removeNull <- function(x){
  res <- x[which(!sapply(x, is.null))]
}

makePlotData <- function(rates, probs = c(0.025, 0.975)){
  rates <- removeNull(rates)
  res <- lapply(names(rates), function(e) computeMeanInterval(e, rates = rates, probs = probs))
  plotdata <- do.call(rbind, res)
  plotdata$item <- factor(plotdata$item,
                          levels = c("speciation rate", "extinction rate", "speciation time", "extinction time",
                                     "net-diversification rate", "relative-extinction rate"))
  return(plotdata)
}


#' Plot Diversification Rates
#'
#' Plots the output of an episodic diversification rate analysis
#'
#' Plots the output of episodic diversification rate analyses. Takes as input the
#' output of processDivRates() and plotting parameters. For now, only variable names
#' (under "item") that contain the word "rate" are included in the plot.
#'
#' The return object can be manipulated. For example, you can change the axis labels,
#' the color palette, whether the axes are to be linked, or the overall plotting style/theme,
#' just as with any ggplot object.
#'
#'
#' @param plotdata (data frame) A table consisting of the following columns:
#' \itemize{
#' \item{mean - the mean of the variable},
#' \item{lower - the lower bounds of the credibility interval}
#' \item{upper - the upper bounds of the credibility interval}
#' \item{time - the time units on the x-axis}
#' \item{item - the variable name, e.g. "speciation rate" or "extinction rate"}
#' }
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
#' primates <- processDivRates(speciation_time_log = speciation_time_file,
#'                             speciation_rate_log = speciation_rate_file,
#'                             extinction_time_log = extinction_time_file,
#'                             extinction_rate_log = extinction_rate_file,
#'                             burnin = 0.25)
#'
#' # generate the plot data
#' plotdata <- makePlotData(primates, probs = c(0.05, 0.95)) # Specify credibility interval for the ribbon plot
#'
#' # then plot results:
#' p <- plotDivRates2(plotdata);p
#'
#' # let's say we want to change the x-axis
#' p <- p + xlab("Millions of years ago");p
#'
#' # let's say we don't want to plot relative-extinction rate,
#' # and use the same y-axis for all three rates
#' plotdata2 <- plotdata %>%
#'                filter(item != "relative-extinction rate")
#' p2 <- plotDivRates2(plotdata2)
#' p2 <- p2 +
#'        facet_wrap(vars(item),
#'        scale = "fixed");p2
#' }
#'
#' @export
plotDivRates2 <- function(plotdata){
  p <- plotdata %>%
    subset(grepl("rate", item)) %>%
    ggplot(aes(x, mean, color = item))  +
    geom_line(aes(time, mean)) +
    geom_ribbon(aes(x = time,
                    ymin = lower,
                    ymax = upper,
                    fill = item),
                alpha = 0.4) +
    scale_x_reverse() +
    xlab("time") +
    theme_bw() +
    theme(axis.title.y=element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    facet_wrap(vars(item), scales = "free_y")

  return(p)
}

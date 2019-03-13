#' Plot Diversification Rates
#'
#' Plots the output of a episodic diversification rate analysis
#'
#' Plots the output of diversification rate analyses. Takes as input the
#' output of processDivRates() and plotting parameters. Does not return an
#' object. For now, valid figure types include: "speciation rate",
#' "extinction rate","net-diversification rate", and "relative-extinction rate".
#' If colors are not provided (parameter col), the plot defaults to preset colors.
#'
#'
#' @param output (list; no default) The processed output for plotting (output of processDivRates()).
#' @param fig_types (character vector, c("speciation rate", "extinction rate", "net-diversification rate", "relative-extinction rate")) Which aspects of the model to visualize. See details for a complete description.
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
#' pdf("primates_example.pdf")
#' plotDivRates(output = primates)
#' dev.off()


plotDivRates <- function(output,
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
  if (is.list(output) == FALSE) stop("output must be list of processed rates")
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
                       "relative-extinction rate")

  invalid_fig_types <- fig_types[!fig_types %in% valid_fig_types]

  if ( length(invalid_fig_types) > 0 ) {
    stop("\nThe following figure types are invalid: ", paste(invalid_fig_types, collapse = ", "), ".",
         "\nValid options are: ", paste(valid_fig_types, collapse = ", ") ,".")
  }


  # Make color vector
  if ( is.null(col) ) {
    col <- c("speciation rate"="#984EA3",
             "speciation shift times"="#984EA3",
             "speciation Bayes factors"="#984EA3",
             "extinction rate"="#E41A1C",
             "extinction shift times"="#E41A1C",
             "extinction Bayes factors"="#E41A1C",
             "fossilization rate"="#ffa31a",
             "fossilization shift times"="#ffa31a",
             "fossilization Bayes factors"="#ffa31a",
             "net-diversification rate"="#377EB8",
             "relative-extinction rate"="#FF7F00",
             "mass extinction times"="#4DAF4A",
             "mass extinction Bayes factors"="#4DAF4A")
  } else {
    names(col) <- fig_types
  }

  # Compute the axes
  intervals <- output$`speciation time`[1,grep("interval", colnames(output$`speciation time`))]
  tree_age <- max(intervals)
  num_intervals <- length(intervals) - 1
  plot_at <- 0:num_intervals
  interval_size <- tree_age/num_intervals
  labels <- pretty( c(0,tree_age) ) + offset
  labels_at <- num_intervals - ( pretty(c(0,tree_age)) / interval_size )
  ages <- seq( 0,tree_age, length.out = num_intervals + 1 )

  for( type in fig_types ) {

      this_output <- output[[type]][ ,grep("[0-9]", colnames(output[[type]])) ]
      mean_this_output <- colMeans(this_output)
      quantiles_this_output <- apply(this_output, 2, quantile, prob = c(0.025,0.975))

      # find the limits of the y-axis
      # we always want the speciation and extinction rates to be on the same scale
      if ( type %in% c("speciation rate","extinction rate")) {

        quantiles_speciation <-
          apply(output[["speciation rate"]][ ,grep("[0-9]", colnames(output[["speciation rate"]])) ],
                2, quantile, prob=c(0.025,0.975))
        quantiles_extinction <-
          apply(output[["extinction rate"]][ ,grep("[0-9]", colnames(output[["extinction rate"]])) ],
                2, quantile, prob=c(0.025,0.975))
        ylim <- c( min(0, quantiles_speciation, quantiles_extinction),
                   max(quantiles_speciation, quantiles_extinction) )

      } else {
        ylim <- c( min(0, quantiles_this_output),
                   max(quantiles_this_output) )
      }

      #### plot
        plot( x = plot_at, y = c( rev(mean_this_output)),
              type = "l", ylim = ylim, xaxt = xaxt, col = col[type],
              ylab = "rate", main = type, xlab = xlab,...)
        polygon( x = c( 0:ncol(quantiles_this_output),
                        ncol(quantiles_this_output):0 ),
                 y = c( rev(c(quantiles_this_output[1,1],
                          quantiles_this_output[1,])),
                        c(quantiles_this_output[2,1],
                              quantiles_this_output[2,])),
                 border = NA, col = paste(col[type], col_alpha, sep="") )
        axis(1, at = labels_at, labels = labels)
      }

    }


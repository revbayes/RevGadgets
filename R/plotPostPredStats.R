#' plot Posterior Predictive Statistics
#'
#' Plots the posterior predictive statistics data
#'
#' Produces one ggplot object per metric. Intended
#' to plot the results of the RevBayes tutorial:
#' Assessing Phylogenetic Reliability Using RevBayes and P3
#' Model adequacy testing using posterior prediction (Data Version).
#'
#' @param data (list of data frames; no default) A list of data frames
#' of the empirical and simulated values, such as the output of
#' processPostPredStats.R
#' 
#' @param prob (vector of numerics; default c(0.9, 0.95)) The
#' posterior-predictive intervals to shade.
#' 
#' @param col (vector of colors; default NULL) The colors for each quantile.
#' Defaults to blue and red.
#' 
#' @param side (character; default "both") Whether the plotted/colored 
#' intervals are on "both" sides, the "left" side, or the "right" 
#' side of the distribution.
#' 
#' @param type (character; default "strict") Whether equal values are
#' considered as less extreme as the observed data ("strict") or half of the
#' equal values are considered to be higher and half to be lower ("midpoint")
#' 
#' @param PPES (boolean; default FALSE) Whether we provide the posterior
#' predictive effect size (PPES).
#' 
#' @param ... Additional arguments are passed to stats::density(). 
#'
#' @return A list of ggplot objects, where each plot contains a density
#' distribution of the predicted values and a dashed line of the empirical
#' value. The blue shaded region of the density plot corresponds to the 5\%
#' two-sided quantile and the orange corresponds to the 2\% two-sided quantile.
#'
#' @details Each plot shows the rejection region for the provided quantiles,
#' as well as a p-value for the observed statistic. If side="left" (or "right"),
#' then the p-value is the fraction of simulated statistics that are less than
#' ( or greater than) or equal to the observed statistic. If side="both", then
#' the p-value is calculated by first fitting a KDE to the samples, then
#' computing the fraction of simulated statistics with density lower than the
#' density of he observed statistic; in this sense, the "both" option computes
#' the size of HPD defined by the observed statistic.
#'
#' @examples
#'
#' \donttest{
#' # download the example datasets to working directory
#'
#' url_emp <-
#'    "https://revbayes.github.io/tutorials/intro/data/empirical_data_pps_example.csv"
#' dest_path_emp <- "empirical_data_pps_example.csv"
#' download.file(url_emp, dest_path_emp)
#'
#' url_sim <-
#'    "https://revbayes.github.io/tutorials/intro/data/simulated_data_pps_example.csv"
#' dest_path_sim <- "simulated_data_pps_example.csv"
#' download.file(url_sim, dest_path_sim)
#'
#' # to run on your own data, change this to the path to your data file
#' file_sim <- dest_path_sim
#' file_emp <- dest_path_emp
#'
#' t <- processPostPredStats(path_sim = file_sim,
#'                          path_emp = file_emp)
#' plots <- plotPostPredStats(data = t)
#' plots[[1]]
#'
#' # remove files
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_sim, dest_path_emp)
#'
#' }
#'
#' @export

plotPostPredStats <- function(data,
                              prob  = c(0.9, 0.95),
                              col   = NULL,
                              side  = "both",
                              type  = "strict",
                              PPES  = FALSE,
                              ...) {
  if (is.list(data) == FALSE)
    stop("Argument data must be a list.")
  if ("simulated" %in% names(data) == FALSE)
    stop("Argument data must be a contain an element called simulated.")
  if ("observed" %in% names(data) == FALSE)
    stop("Argument data must be a contain an element called observed.")
  if (is.data.frame(data$simulated) == FALSE)
    stop("data$simulated must be a data.frame.")
  if (is.data.frame(data$observed) == FALSE)
    stop("data$observed must be a data.frame.")
  if (side %in% c("both", "left", "right") == FALSE)
    stop("Invalid side argument.")
  if (type %in% c("strict", "midpoint") == FALSE)
    stop("Invalid type argument.")

  if (is.null(col)) {
    col <- grDevices::colorRampPalette(colFun(2))(length(prob))
  }
  if (length(col) != length(prob))
    stop("Number of colors does not match the number of quantiles.")

  # sort the probs
  prob <- sort(prob)

  # retrieve the data
  sim <- data$simulated
  obs <- data$observed

  # check names of statistics
  sim_stats <- colnames(sim)
  obs_stats <- colnames(obs)

  # get the intersect
  names <- intersect(sim_stats, obs_stats)
  if (length(names) == 0) {
    stop("data$simulated and data$observed do not contain the same statistics.")
  }
  if (length(setdiff(obs_stats, sim_stats)) > 0) {
    warning(
      "data$simulated and data$observed do not share all the same statistics.
      Only the shared statistics will be plotted."
    )
  }

  # containers for plots
  plots <- vector("list", length(names))
  for (i in seq_len(length(names))) {
    # xlim values
    min_value <- min(sim[, i], obs[[i]])
    max_value <- max(sim[, i], obs[[i]])
    spread_value <- max_value - min_value
#    spread_value <- ifelse( spread_value > 0, spread_value, min_value*0.01 )
#    spread_value <- ifelse( spread_value > 0, spread_value, 0.05 )

    # fit a kernel density
    kde <- density(sim[, i], ...)
    pdf <- approxfun(kde)

    # compute the p-value
    if (side == "both") {
      # compute the density at the sampled point
      dens <- pdf(obs[, i])

      # compute the amount of the distribution with equal or
      # lower density
      if (is.na(dens)) {
        p_value <- 0.0
      } else {
        p_value <- mean(pdf(sim[, i]) <= dens)
      }

    } else if (side == "left" && type == "strict") {
      # lower p-value
      p_value <- mean(sim[, i] < obs[, i])
    } else if (side == "left" && type == "midpoint") {
      # lower midpoint p-value
      p_value <- mean(sim[, i] < obs[, i]) + mean(sim[, i] == obs[, i]) / 2.0
    } else if (side == "right" && type == "strict" ) {
      # upper p-value
      p_value <- mean(sim[, i] > obs[, i])
    } else if (side == "right" && type == "midpoint" ) {
      # upper midpoint p-value
      p_value <- mean(sim[, i] > obs[, i]) + mean(sim[, i] == obs[, i]) / 2.0
    }

    # compute the posterior predictive effect size
    ppes <- abs( obs[, i] - stats::median(sim[, i]) ) / stats::sd(sim[, i])
    ppes <- ifelse( stats::sd(sim[, i]) == 0, 0, ppes )

    # make dataframe of plotting data
    df <- data.frame((kde)[c("x", "y")])

    # make the p-value label
    p_lab <- paste0("p=", sprintf("%.3f", p_value))
    p_x   <- max_value + 0.25 * spread_value
    p_y   <- max(df$y)

    # make the ppes label
    ppes_lab <- paste0("ppes=", sprintf("%.3f", ppes))
    ppes_x   <- min_value
    ppes_y   <- max(df$y)

    # plot
    p <- ggplot2::ggplot(df, ggplot2::aes(x, y))
    for (q in seq_len(length(prob))) {
      this_q <- prob[q]
      if (side == "left") {
        l <- 1 - this_q
        p <-
          p +
          ggplot2::geom_area(data = df[df$x <= quantile(sim[, i], prob = l), ],
                                 fill = col[q])
      } else if (side == "right") {
        u <- this_q
        p <-
          p +
          ggplot2::geom_area(data = df[df$x >= quantile(sim[, i], prob = u), ],
                                 fill = col[q])
      } else {
        # compute the quantiles
        l <- (1 - this_q) / 2
        u <- 1 - l
        p <-
          p +
          ggplot2::geom_area(data = df[df$x <= quantile(sim[, i], prob = l), ],
                                 fill = col[q])
        p <-
          p +
          ggplot2::geom_area(data = df[df$x >= quantile(sim[, i], prob = u), ],
                                 fill = col[q])
      }

    }
    p <- p + ggplot2::geom_line() +
      ggplot2::xlim(c(
        min_value - 0.25 * spread_value,
        max_value + 0.25 * spread_value
      )) +
      ggplot2::geom_vline(xintercept = obs[[i]],
                          linetype = "dashed") +
      ggplot2::xlab(names[i]) +
      ggplot2::ylab("Density") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::annotate(
        "text",
        x = p_x,
        y = p_y,
        label = p_lab,
        size = 3,
        hjust = 1
      )
      if ( PPES ) {
        p <- p + ggplot2::annotate(
                 "text",
                 x = ppes_x,
                 y = ppes_y,
                 label = ppes_lab,
                 size = 3,
                 hjust = 1
                 )
      }

    plots[[i]] <- p

  }


  # name each element of the list according to the statistic
  names(plots) <- colnames(data[[2]])

  return(plots)

}

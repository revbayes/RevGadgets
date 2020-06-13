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
#' @return A list of ggplot objects, where each plot contains a density distribution
#' of the predicted values and a dashed line of the empirical value. The light blue
#' shaded region of the density plot corresponds to the 5% two-sided quantile and
#' the dark blue corresponds to the 2% two-sided quantile.
#'
#' @examples
#'
#' \dontrun{
#' file_sim <- system.file("extdata",
#'     "PPS/simulated_data_pps_example.csv", package="RevGadgets")
#' file_emp <- system.file("extdata",
#'     "PPS/empirical_data_pps_example.csv", package="RevGadgets")
#' t <- processPostPredStats(path_sim = file_sim,
#'                          path_emp = file_emp)
#' plots <- plotPostPredStats(data = t)
#'
#' }
#'
#' @export

plotPostPredStats <- function(data) {
  sim <- data$simulated
  emp <- data$empirical
  names <- colnames( emp )

  min_value <- c()
  max_value <- c()
  spread_value <- c()
  plots <- list()
  for ( i in 1:length(names)) {
    # xlim values
    min_value[i] <- min(sim[,i], emp[[i]])
    max_value[i] <- max(sim[,i], emp[[i]])
    spread_value[i] <- max_value[i] - min_value[i]
    # quantile values for fill
    q01  <- quantile(sim[,i],prob = 0.01)
    q025 <- quantile(sim[,i],prob = 0.025)
    q975 <- quantile(sim[,i],prob = 0.975)
    q99  <- quantile(sim[,i],prob = 0.99)
    # make dataframe of plotting data
    df <- data.frame(density(sim[,i])[c("x", "y")])
    # plot
    p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
      ggplot2::geom_area(data = subset(df, x < q025),
                         fill =  "#14d2dc") +
      ggplot2::geom_area(data = subset(df, x < q01),
                         fill = "#005ac8") +
      ggplot2::geom_area(data = subset(df, x > q975),
                         fill = "#14d2dc") +
      ggplot2::geom_area(data = subset(df, x > q99),
                         fill = "#005ac8") +
      ggplot2::geom_line() +
      ggplot2::xlim(c(min_value[i]-0.25*spread_value[i],
                      max_value[i]+0.25*spread_value[i])) +
      ggplot2::geom_vline(xintercept = emp[[i]],
                          linetype = "dashed") +
      ggplot2::xlab(names[i]) +
      ggplot2::ylab("Density") +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major =
                       ggplot2::element_blank(),
                     panel.grid.minor =
                       ggplot2::element_blank())
    plots[[i]] <- p
  }
  return(plots)
}

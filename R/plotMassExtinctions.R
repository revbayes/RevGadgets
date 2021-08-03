#' Plot Mass Extinction Support
#'
#' Plots the support (as 2ln Bayes factors) for mass extinctions.
#'
#' Works only for analyses with a fixed grid where mass extinctions may occur.
#'
#' The return object can be manipulated. For example, you can change the axis
#' labels, the color palette, whether the axes are to be linked, or the overall
#' plotting style/theme, just as with any ggplot object.
#'
#' @param mass_extinction_trace (list; no default) The processed Rev output of
#' the mass extinction probabilities (output of readTrace()).
#' @param mass_extinction_times (numeric; no default) Vector of the fixed grid
#' of times at which mass extinctions were allowed to occur.
#' @param mass_extinction_name (character; no default) The name of the mass
#' extinction probability parameter (e.g. "mass_extinction_probabilities")
#' for which support is to be calculated/plotted.
#' @param prior_prob (numeric; no default) The per-interval prior probability
#' of a mass extinction (one minus the p parameter in RevBayes'
#' dnReversibleJumpMixture()).
#' @param return_2lnBF (logical; TRUE) Should the 2ln(BF) be returned (if TRUE)
#' or simply the BF (if FALSE)?
#'
#' @return A ggplot object
#'
#' @references
#'
#' Kass and Raftery (1995) Bayes Factors.
#' \emph{JASA}, \bold{90 (430)}, 773-795.
#'
#' @examples
#'\donttest{
#'
#' # download the example dataset to working directory
#' url <-
#'   "https://revbayes.github.io/tutorials/intro/data/crocs_mass_extinction_probabilities.log"
#' dest_path <- "crocs_mass_extinction_probabilities.log"
#' download.file(url, dest_path)
#'
#' # to run on your own data, change this to the path to your data file
#' mass_extinction_probability_file <- dest_path
#'
#' mass_extinction_probabilities <-
#'               readTrace(mass_extinction_probability_file,burnin = 0.25)
#'
#' # prior probability of mass extinction at any time
#' prior_n_expected <- 0.1
#' n_intervals <- 100
#' prior_prob <- prior_n_expected/(n_intervals-1)
#'
#' # times when mass extinctions were allowed
#' tree_age <- 243.5
#' interval_times <- tree_age * seq(1/n_intervals,(n_intervals-1) /
#'                                                    n_intervals,1/n_intervals)
#'
#' # then plot results:
#' p <- plotMassExtinctions(mass_extinction_trace=mass_extinction_probabilities,
#'                          mass_extinction_times=interval_times,
#'                          mass_extinction_name="mass_extinction_probabilities"
#'                          ,prior_prob);p
#'
#' # remove file
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path)
#' }
#'
#' @export

plotMassExtinctions <- function(mass_extinction_trace,
                                mass_extinction_times,
                                mass_extinction_name,
                                prior_prob,
                                return_2lnBF = TRUE) {
  # Condense log lists into data frames
  me_log <- do.call(rbind, mass_extinction_trace)

  # Find parameter and times, remove rest of trace
  is_me <-
    grepl(paste0(mass_extinction_name, "["), names(me_log), fixed = TRUE)

  if (sum(is_me) < 1) {
    stop(
      paste0(
        "Cannot find mass extinction probability parameter named \"",
        rate.name,
        "\" in mass_extinction_trace"
      )
    )
  }

  me_log <- me_log[, is_me]

  if (dim(me_log)[2] != length(mass_extinction_times)) {
    stop("Number of provided interval times does not match number
         of mass extinctions.")
  }

  niter <- dim(me_log)[1]

  me_prob <-
    1.0 - apply(me_log, 2, function(x) {
      sum(x < .Machine$double.eps)
    }) / niter

  # Avoid infinite Bayes Factors
  if (any(me_prob == 0.0)) {
    me_prob[me_prob == 0] <- 1 / (niter + 1)
  } else if (any(me_prob == 1.0)) {
    me_prob[me_prob == 1.0] <- niter / (niter + 1)
  }

  # Prior odds are 1 under a MRF model
  posterior_odds <- me_prob / (1 - me_prob)
  prior_odds <- prior_prob / (1 - prior_prob)
  BF <- posterior_odds / prior_odds

  if (return_2lnBF == TRUE) {
    BF <- 2 * log(BF)
  }

  # Data frame for ggplot
  bfdf <- data.frame(bf = BF, time = mass_extinction_times)

  p <-
    ggplot2::ggplot(data = bfdf,
                    mapping = ggplot2::aes(x = time, y = bf)) +
    ggplot2::geom_rect(mapping = ggplot2::aes(
      ymin = -Inf,
      ymax = 2,
      xmin = -Inf,
      xmax = Inf
    ),
    fill = "grey70") +
    ggplot2::geom_rect(mapping = ggplot2::aes(
      ymin = 2,
      ymax = 6,
      xmin = -Inf,
      xmax = Inf
    ),
    fill = "grey80") +
    ggplot2::geom_rect(mapping = ggplot2::aes(
      ymin = 6,
      ymax = 10,
      xmin = -Inf,
      xmax = Inf
    ),
    fill = "grey90") +
    #ggplot2::geom_col(ggplot2::aes(fill = "red")) +
    ggplot2::geom_col(fill = colFun(2)[2]) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_x_reverse() +
    ggplot2::xlab("million years ago") +
    ggplot2::ylab("2 log Bayes factor") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_rect(color = "black", fill = NA)
    )

  return(p)
}

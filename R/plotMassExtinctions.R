#' Plot Mass Extinction Support
#'
#' Plots the support (as 2ln Bayes factors) for mass extinctions.
#'
#' Works only for analyses with a fixed grid where mass extinctions may occur.
#'
#' The return object can be manipulated. For example, you can change the axis labels,
#' the color palette, whether the axes are to be linked, or the overall plotting style/theme,
#' just as with any ggplot object.
#'
#' @param mass.extinction.trace (list; no default) The processed Rev output of the mass extinction probabilities (output of readTrace()).
#' @param mass.extinction.times (numeric; no default) Vector of the fixed grid of times at which mass extinctions were allowed to occur.
#' @param mass.extinction.name (character; no default) The name of the mass extinction probability parameter (e.g. "mass_extinction_probabilities") for which support is to be calculated/plotted.
#' @param prior.prob (numeric; no default) The per-interval prior probability of a mass extinction (one minus the p parameter in RevBayes' dnReversibleJumpMixture()).
#' @param lnBF (logical; TRUE) Should the 2ln(BF) be returned (if TRUE) or simply the BF (if FALSE)?
#'
#' @return A ggplot object
#'
#' @references
#'
#' Kass and Raftery (1995) Bayes Factors.
#' \emph{JASA}, \bold{90 (430)}, 773-795.
#'
#' @examples
#'\dontrun{
#'
#' mass_extinction_probability_file <- system.file("extdata",
#'     "mass_extinction/crocs_mass_extinction_probabilities.p", package="RevGadgets")
#'
#' mass_extinction_probabilities <- readTrace(mass_extinction_probability_file,burnin = 0.25)
#'
#' # prior probability of mass extinction at any time
#' prior_n_expected <- 0.1
#' n_intervals <- 100
#' prior_prob <- prior_n_expected/(n_intervals-1)
#'
#' # times when mass extinctions were allowed
#' tree_age <- 243.5
#' interval_times <- tree_age * seq(1/n_intervals,(n_intervals-1)/n_intervals,1/n_intervals)
#'
#' # then plot results:
#' p <- plotMassExtinctions(mass.extinction.trace=mass_extinction_probabilities,mass.extinction.times=interval_times,mass.extinction.name="mass_extinction_probabilities",prior_prob);p
#'}
#'
#' @export

plotMassExtinctions <- function(mass.extinction.trace,mass.extinction.times,mass.extinction.name,prior.prob,return.2lnBF=TRUE) {

  # Condense log lists into data frames
  me_log <- do.call(rbind,mass.extinction.trace)

  # Find parameter and times, remove rest of trace
  is_me <- grepl(paste0(mass.extinction.name,"["),names(me_log),fixed=TRUE)

  if ( sum(is_me) < 1 ) {
    stop(paste0("Cannot find mass extinction probability parameter named \"",rate.name,"\" in mass.extinction.trace"))
  }

  me_log <- me_log[,is_me]

  if ( dim(me_log)[2] != length(mass.extinction.times) ) {
    stop("Number of provided interval times does not match number of mass extinctions.")
  }

  niter <- dim(me_log)[1]

  me_prob <- 1.0 - apply(me_log,2,function(x){sum(x < .Machine$double.eps)})/niter

  # Avoid infinite Bayes Factors
  if ( any(me_prob == 0.0) ) {
    me_prob[me_prob == 0] <- 1 / (niter + 1)
  } else if ( any(me_prob == 1.0) ) {
    me_prob[me_prob == 1.0] <- niter / (niter + 1)
  }

  # Prior odds are 1 under a MRF model
  posterior_odds <- me_prob / (1 - me_prob)
  prior_odds <- prior.prob / (1 - prior.prob)
  BF <- posterior_odds / prior_odds

  if ( return.2lnBF == TRUE ) {
    BF <- 2 * log(BF)
  }

  # Data frame for ggplot
  bfdf <- data.frame(bf=BF,time=mass.extinction.times)

  p <- ggplot2::ggplot(data=bfdf, mapping=aes(x=time,y=bf)) +
    ggplot2::geom_rect(mapping=aes(ymin=-Inf,ymax=2,xmin=-Inf,xmax=Inf),fill="grey70") +
    ggplot2::geom_rect(mapping=aes(ymin=2,ymax=6,xmin=-Inf,xmax=Inf),fill="grey80") +
    ggplot2::geom_rect(mapping=aes(ymin=6,ymax=10,xmin=-Inf,xmax=Inf),fill="grey90") +
    ggplot2::geom_col(aes(fill="red")) +
    ggplot2::geom_hline(yintercept=0,linetype="dashed") +
    ggplot2::scale_x_reverse() + xlab("million years ago") +
    ggplot2::ylab("2 log Bayes factor") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position="none",panel.border=element_rect(color="black",fill=NA))

  return(p)
}
#' Priors from MCMC samples
#'
#'
#' Turn posterior samples collected by MCMC into a parametric prior
#' distribution.
#'
#' The distributions are fit by the method of moments.
#' The function allows inflating the prior variance relative to the posterior
#' being supplied.
#'
#' @param samples (numeric vector; no default) MCMC samples for a single
#' parameter.
#' @param distribution (character; no default) The distribution to fit. Options
#' are gamma for strictly positive
#' parameters and normal for unbounded parameters.
#' @param variance_inflation_factor (single numeric value; default = 2.0) Makes
#' the prior variance larger than the variance of the posterior
#'
#' @return Numeric vector of parameters with names (to avoid rate/scale and
#' var/sd confusion).
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
#' primates <- processDivRates(speciation_time_log = speciation_time_file,
#'                             speciation_rate_log = speciation_rate_file,
#'                             extinction_time_log = extinction_time_file,
#'                             extinction_rate_log = extinction_rate_file,
#'                             burnin = 0.25)
#'
#' speciation_rates <-
#'        dplyr::pull(primates[which(primates$item == "speciation rate"),],
#'                   "value")
#' speciation_1_gamma_prior <-
#'           posteriorSamplesToParametricPrior(speciation_rates,"gamma")
#'
#' # remove files
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_sp_times, dest_path_ex_times,
#'             dest_path_sp_rates, dest_path_ex_rates)
#' }
#'
#' @export

posteriorSamplesToParametricPrior <- function(samples,
                                              distribution,
                                              variance_inflation_factor =
                                                2.0) {
  if (!is.numeric(samples)) {
    stop("This function requires numeric input data.")
  }
  if (!length(samples) >= 3) {
    stop("Must provide at least 3 samples to compute mean/variance
         (many more is recommended).")
  }
  if (variance_inflation_factor < 1.0) {
    stop("Not recommended to decrease variance from posterior.")
  }
  if (length(distribution) > 1) {
    stop("Only one distribution can be selected.")
  }

  par <- numeric()

  sample_mean <- mean(samples)
  sample_var <- variance_inflation_factor * stats::var(samples)

  if (grepl("gam", tolower(distribution))) {
    # Fit a gamma distribution
    b <- sample_mean / sample_var
    a <- sample_mean * b
    par <- c(a, b)
    names(par) <- c("gamma.shape", "gamma.rate")
  } else if (grepl("^norm", tolower(distribution))) {
    # Fit a normal distribution
    par <- c(sample_mean, sqrt(sample_var))
    names(par) <- c("normal.mean", "normal.sd")
  } else {
    stop("Unrecognized distribution choice.")
  }

  return(par)
}

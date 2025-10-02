#' Process Simulation-Based Calibration (SBC) Results
#'
#' @description
#' Reads and processes the output of a set of simulation replicates to assesss
#' calibration. For each specified parameter, it calculates the coverage
#' frequency across a range of Higest Posterior Density (HPD) interval widths.
#'
#'
#' @param path (character; no default) The path to the main output directory which contains
#'   the individual simulation subdirectories (e.g., "output_<analysis_name>").
#' @param results_path (character; optional) The path to a directory where the
#'   processed RDS file for each parameter will be saved. If NULL (the default),
#'   the results will not be saved to disk.
#' @param parameters (character vector; optional) A vector of parameter names to process.
#'   If NULL (the default), the function will automatically find all parameters in the
#'   output files.
#' @param n_reps (numeric; default: 1000) The total number of simulation replicates.
#' @param n_bins (numeric; default 50) The number of bins to use for the HPD widths, ranging
#'   from 0 to 1.
#' @param burnin (numeric; default: 0.25) The fraction of samples from the MCMC chain to
#'   discard as burn-in.
#' @param verbose (logical; default: TRUE) Whether to print progress bar and summary
#'   messages to the console.
#'
#' @return A named list of data frames. Each data frame corresponds to a parameter and contains
#' the HPD witdth, the count of simulations where the true value was within the HPD interval,
#' the total count and the coverage frequency.
#'
#' @examples
#' \dontrun{
#' # assuming simulation outputs are in "output_my_analysis"
#' sbc_results <- processSBC(
#'   path = "output_my_analysis",
#'   results_path = "results_my_analysis",
#'   n_reps = 1000,
#'   parameters = c("alpha", "sigma")
#' )
#'
#' # The results for the 'alpha' parameter can be accessed as:
#' alpha_results <- sbc_results$alpha
#' }
#'
#' @export
#' @importFrom coda HPDinterval as.mcmc
#' @importFrom utils txtProgressBar setTxtProgressBar
processSBC <- function(
    path,
    results_path = NULL,
    parameters = NULL,
    n_reps = 1000,
    n_bins = 50,
    burnin = 0.25,
    verbose = TRUE) {
  # check for the required 'path' argument
  if (missing(path)) {
    stop("Argument 'path' is required. Please provide the path to the simulation output directory.")
  }

  # if a results path is provided, create the directory
  if (!is.null(results_path)) {
    dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
  }

  # auto-detect parameters if not provided
  if (is.null(parameters)) {
    first_sim_file <- file.path(path, "Validation_Sim_0", "posterior_samples.var")
    if (!file.exists(first_sim_file)) {
      stop("Cannot find output files. Checked for: ", first_sim_file)
    }
    col_names <- colnames(read.table(first_sim_file, header = TRUE, sep = "\t", check.names = FALSE, nrow = 1))

    parameters <- setdiff(col_names, c("Iteration", "Posterior", "Likelihood", "Prior", "branch_rates"))
  }

  if (verbose) {
    cat("Processing parameters:\n", paste(parameters, collapse = ", "), "\n\n")
  }

  # initialize results list
  all_results <- list()

  # iterate over each parameter
  for (param in parameters) {

    # initialize the data frame to store coverage probabilities
    hpd_width <- seq(from = 0.0, to 1.0, length.out = n_bins + 1)
    coverage_probs <- data.frame(total_count = 0, in_count = 0, hpd_width = hpd_width)

    # initialize the progress bar if verbose
    if (verbose) {
      pb <- txtProgressBar(min = 0, max = n_reps, char = "*", style = 3)
    }

    # iterate over each simulation replicate
    for (i in 1:n_reps) {
      if (verbose) utils::setTxtProgressBar(pb, i)

      sim_dir <- file.path(path, path0("Validation_Sim_", i - 1))
      posterior_file <- file.path(sim_dir, "posterior_samples.var")

      if (!file.exists(posterior_file)) next

      data <- read.table(posterior_file, header = TRUE, sep = "\t", check.names = FALSE)

      # extract samples and apply burn-in
      num_samples <- nrow(data)
      start_index <- floor(burnin * num_samples) + 1
      x <- coda::as.mcmc(data[start_index:num_samples, param])

      # read the true value
      true_val_ext <- ifelse(param == "branch_rates", ".out", ".txt")
      true_val_file <- file.path(sim_dir, paste0(param, true_val_ext))

      if (!file.exists(true_val_file)) next
      true_val <- read.table(true_val_file)[1, 1]

      # calculate coverage probabilities for each HPD width
      for (k in 1:(n_bins + 1)) {
        hpd <- coda::HPDinterval(x, prob = hpd_width[k])
        if (true_val >= hpd[1, 1] && true_val <= hpd[1, 2]) {
          coverage_probs$in_count[k] <- coverage_probs$in_count[k] + 1
        }
        coverage_probs$total_count[k] <- coverage_probs$total_count[k] + 1
      }
    }

    if (verbose) close(pb)

    # calculate coverage frequency
    # avoid division by zero if total_count is zero
    coverage_probs$freq <- ifelse(coverage_probs$total_count > 0,
      coverage_probs$in_count / coverage_probs$total_count,
      0
    )

    # save the parameter's results into the main list
    all_results[[param]] <- coverage_probs

    # save to RDS if results_path is provided
    if (!is.null(results_path)) {
      saveRDS(coverage_probs, file = file.path(results_path, paste0(param, ".rds")))
    }

    if (verbose) {
      cat("\nResults for parameter:", param, "\n")
      # A more readable printout
      print(coverage_probs, row.names = FALSE)
      cat("\n")
    }
  } # end of parameters loop

  # return the list containing all results
  return(all_results)
}

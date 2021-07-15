#' process Posterior Predictive Statistics
#'
#' Reads in and processes posterior-predictive statistics
#'
#' @param path_sim (character string; no default) Path to the .csv file
#' containing the simulated data results
#' @param path_emp (character string; no default) Path to the .csv file
#' containing the empirical values
#'
#' @return A list of data frames
#'
#' @examples
#'
#' \dontrun{
#' file_sim <- system.file("extdata",
#'     "PPS/simulated_data_pps_example.csv", package="RevGadgets")
#' file_emp <- system.file("extdata",
#'     "PPS/empirical_data_pps_example.csv", package="RevGadgets")
#' t <- processPostPredStats(path_sim = file_sim,
#'                                           path_emp = file_emp)
#' }
#'
#' @export

processPostPredStats <- function(path_sim, path_emp) {
  # parameter checks
  paths <- c(path_sim, path_emp)

  character_paths_are_strings <- is.character(paths)
  if (any(character_paths_are_strings == FALSE) == TRUE) {
    # print out the ones that are not character strings
    cat("Some paths are not character strings:",
        paste0("\t", paths[character_paths_are_strings == FALSE]),
        sep = "\n")
    stop()
  }

  do_files_exist <- file.exists(paths)
  if (any(do_files_exist == FALSE) == TRUE) {
    # print out paths to files that don't exist
    cat("Some files do not exist:",
        paste0("\t", paths[do_files_exist == FALSE]), sep = "\n")
    stop()
  }

  # read in data
  posterior_predictive_statistics <- read.table(path_sim,
                                                header = TRUE,
                                                sep = ",",
                                                check.names = FALSE)

  observed_statistics <- read.table(path_emp,
                                    header = TRUE,
                                    sep = ",",
                                    check.names = FALSE)

  # check that the statistics match
  if (length(setdiff(
    colnames(posterior_predictive_statistics),
    colnames(observed_statistics)
  )) > 0) {
    cat("Simulated and observed files do not have the same statistics.\n")
    stop()
  }

  # return list
  return(list(simulated = posterior_predictive_statistics,
              observed  = observed_statistics))
}


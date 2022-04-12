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
#'
#' # remove files
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_sim, dest_path_emp)
#' }
#'
#'
#' @export

processPostPredStats <- function(path_sim, path_emp) {
  # parameter checks
  paths <- c(path_sim, path_emp)

  character_paths_are_strings <- is.character(paths)
  if (any(character_paths_are_strings == FALSE) == TRUE) {
    # print out the ones that are not character strings
    stop(
      paste0("Some paths are not character strings:",
        paste0("\t", paths[character_paths_are_strings == FALSE]),
        sep = "\n")
    )
  }

  do_files_exist <- file.exists(paths)
  if (any(do_files_exist == FALSE) == TRUE) {
    # print out paths to files that don't exist
    stop(
      paste0("Some files do not exist:",
        paste0("\t", paths[do_files_exist == FALSE]), sep = "\n")
    )
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
    stop("Simulated and observed files do not have the same statistics.\n")
  }

  # return list
  return(list(simulated = posterior_predictive_statistics,
              observed  = observed_statistics))
}


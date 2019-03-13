#' Read trace
#'
#' Reads in MCMC log files
#'
#' Reads in one or multiple MCMC log files from the same analysis
#' and discards a user-specified burn-in, compatible with multiple monitor types.
#' If the trace contains vectors of vectors and the user does not specify
#' format = "complex", readTrace() will read in those columns as factors
#' rather than as numeric vectors.
#'
#' @param path (vector of character strings; no default) File path(s) to trace file.
#' @param format (single character string; default = simple) Indicates type of
#' MCMC trace, complex indicates cases where trace contains vectors of vectors/
#' matrices - mnStochasticVariable monitor will sometimes be of this type.
#' @param delim (single character string; default = "\\t") Delimiter of file.
#' @param burnin (single numeric value; default = 0) Fraction of generations to
#' discard (if value provided is between 0 and 1) or number of generations (if
#' value provided is greater than 1).
#' @param check.names (logical; default = FALSE) Passed to read.table(); indicates
#' if read.table() should check column names and replace syntactically invalid
#' characters.
#' @param ... (various) Additional arguments passed to read.table().
#'
#' @return List of dataframes (of length 1 if only 1 log file provided).
#'
#' @examples
#' one_trace <- readTrace(path = "data/sub_models/primates_cytb_covariotide.log")
#' multi_trace <- readTrace(path = c("data/sub_models/primates_cytb_covariotide_run_1.log","data/sub_models/primates_cytb_covariotide_run_1.log"))

readTrace <- function(path, format = "simple",
                      delim="\t", burnin = 0, check.names = FALSE, ...){

  # enforce argument matching

  character_paths_are_strings <- is.character(path)
  if ( any(character_paths_are_strings == FALSE) == TRUE ) {
    # print out the ones that are not character strings
    cat( "Some paths are not character strings:",
         paste0("\t",path[character_paths_are_strings == FALSE]), sep="\n")
    stop()
  }

  do_files_exist <- file.exists(path)
  if ( any(do_files_exist == FALSE) == TRUE ) {
    # print out paths to files that don't exist
    cat( "Some files do not exist:",
         paste0("\t",path[do_files_exist == FALSE]), sep="\n")
    stop()
  }

  format <- match.arg(format, choices = c("simple", "complex"))

  if (is.character(delim) == FALSE) stop("delim must be a single character string")

  if (is.numeric(burnin) == FALSE) stop("burnin must be a single numeric value")
  if (burnin < 0) stop("burnin must be a positive value")

  num_paths <- length(path)

  # check that the file headings match for all traces

  header <- vector("list", num_paths)
  for (i in 1:num_paths) {
    header[[i]] <- colnames(read.table(file = path[i], header = TRUE, sep = delim, check.names = check.names, nrows=0))
  }

  all_headers <- unique(unlist(header))
  for (i in 1:length(header)) {
    if (setequal(all_headers, header[[i]]) == FALSE) {
      stop("Not all headers of trace files match")
    }
  }


  # read in the traces

  if (format == "simple") {
    output <- vector("list", num_paths)
    for (i in 1:num_paths) {

      cat(paste0("Reading in log file ",i),"\n",sep="")

      out <- read.table(file = path[i], header = TRUE,
                        sep = delim, check.names = check.names, ...)

      if (burnin >= nrow(out)) stop("Burnin larger than provided trace file")

      if (burnin >= 1) {
        output[[i]] <- out[(burnin+1):nrow(out), ]
      } else if (burnin < 1 & burnin > 0) {
        discard <- ceiling(burnin*nrow(out))
        output[[i]] <- out[(discard+1):nrow(out), ]
      } else if (burnin == 0) {
        output[[i]] <- out
      } else {
        stop("What have you done?")
        }
    }
  } else if (format == "complex") {
    stop("Complex trace type currently not supported")
  } else {
    stop("Format is not of type simple or complex")
    }

  return(output)
}

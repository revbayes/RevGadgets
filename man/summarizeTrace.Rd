% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarizeTrace.R
\name{summarizeTrace}
\alias{summarizeTrace}
\title{Summarize trace}
\usage{
summarizeTrace(trace, vars)
}
\arguments{
\item{trace}{(list of data frames; no default) Name of a list of data frames,
such as produced by readTrace(). If the readTrace() output
contains multiple traces (such as from multiple runs), summarizeTrace() will
provide summaries for each trace individually, as well as the combined trace.}

\item{vars}{(character or character vector; no default) The name of the
variable(s) to be summarized.}
}
\value{
summarizeTrace() returns a list of the length of provided variables.
For quantitative variables, it returns the mean and 95% credible interval.
For discrete variables, it returns the 95% credible set of states and their
associated probabilities.
}
\description{
Summarizes trace file(s) that have been read into memory
}
\details{
Summarizes a trace file for continuous or discrete characters by
computing the mean and 95\% credible interval for quantitative
character and the 95\% credible set for discrete characters.
}
\examples{

\donttest{
# continuous character only example, one run

# download the example dataset to working directory
url_gtr <-
    "https://revbayes.github.io/tutorials/intro/data/primates_cytb_GTR.log"
dest_path_gtr <- "primates_cytb_GTR.log"
download.file(url_gtr, dest_path_gtr)

# to run on your own data, change this to the path to your data file
file_single <- dest_path_gtr

one_trace <- readTrace(paths = file_single)
trace_sum <- summarizeTrace(trace = one_trace,
                            vars = c("pi[1]","pi[2]","pi[3]","pi[4]"))
trace_sum[["pi[1]"]]

# remove file
# WARNING: only run for example dataset!
# otherwise you might delete your data!
file.remove(dest_path_gtr)

# continuous character example, multiple runs

#' # download the example dataset to working directory
url_1 <-
"https://revbayes.github.io/tutorials/intro/data/primates_cytb_GTR_run_1.log"
dest_path_1 <- "primates_cytb_GTR_run_1.log"
download.file(url_1, dest_path_1)

url_2 <-
"https://revbayes.github.io/tutorials/intro/data/primates_cytb_GTR_run_2.log"
dest_path_2 <- "primates_cytb_GTR_run_2.log"
download.file(url_2, dest_path_2)

# to run on your own data, change this to the path to your data file
file_1 <- dest_path_1
file_2 <- dest_path_2

# read in the multiple trace files
multi_trace <- readTrace(path = c(file_1, file_2), burnin = 0.0)

trace_sum_multi <- summarizeTrace(trace = multi_trace,
                                  vars = c("pi[1]","pi[2]","pi[3]","pi[4]"))
trace_sum_multi[["pi[1]"]]

# remove files
# WARNING: only run for example dataset!
# otherwise you might delete your data!
file.remove(dest_path_1, dest_path_2)

# discrete character example

# download the example dataset to working directory
url_rj <- "https://revbayes.github.io/tutorials/intro/data/freeK_RJ.log"
dest_path_rj <- "freeK_RJ.log"
download.file(url_rj, dest_path_rj)

file <- dest_path_rj
trace <- readTrace(path = file)

trace_sum_discrete <- summarizeTrace(trace = trace,
                                     vars = c("prob_rate_12",
                                              "prob_rate_13",
                                              "prob_rate_31",
                                              "prob_rate_32"))
trace_sum_discrete[["prob_rate_12"]]

#' # remove file
# WARNING: only run for example dataset!
# otherwise you might delete your data!
file.remove(dest_path_rj)
}

}

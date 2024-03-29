% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMuSSE.R
\name{plotMuSSE}
\alias{plotMuSSE}
\title{plotMuSSE}
\usage{
plotMuSSE(rates)
}
\arguments{
\item{rates}{(data.frame; no default) a data frame containing columns
"value", "rate", "hidden_state", "observed_state" (such as the output
of processSSE())}
}
\value{
a ggplot object
}
\description{
plotMuSSE
}
\examples{
\donttest{

# download the example dataset to working directory

url <-
  "https://revbayes.github.io/tutorials/intro/data/primates_BiSSE_activity_period.log"
dest_path <- "primates_BiSSE_activity_period.log"
download.file(url, dest_path)

# to run on your own data, change this to the path to your data file
bisse_file <- dest_path

pdata <- processSSE(bisse_file)
p <- plotMuSSE(pdata);p

# change colors:
p + ggplot2::scale_fill_manual(values = c("red","green"))

# change x-axis label
p + ggplot2::xlab("Rate (events/Ma)")

# remove file
# WARNING: only run for example dataset!
# otherwise you might delete your data!
file.remove(dest_path)
}
}

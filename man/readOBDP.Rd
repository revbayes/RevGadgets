% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readOBDP.R
\name{readOBDP}
\alias{readOBDP}
\title{Read OBDP Outputs}
\usage{
readOBDP(
  start_time_trace_file,
  popSize_distribution_matrices_file,
  trees_trace_file
)
}
\arguments{
\item{start_time_trace_file}{(character; no default)  Trace of the starting times along the MCMC chain.}

\item{popSize_distribution_matrices_file}{(character; no default)  Kt matrices computed with `fnInferAncestralPopSize` in RevBayes.}

\item{trees_trace_file}{(character; no default)  Trace of the trees.}
}
\value{
A data.frame
}
\description{
Reads and formats the outputs of an analysis with the Occurrence Birth Death Process (MCMC parameter 
inference + diversity estimation)
}
\examples{

\dontrun{
# first run readOBDP()
start_time_trace_file <- 
     system.file("extdata", "obdp/start_time_trace.p", package="RevGadgets")
popSize_distribution_matrices_file <- 
     system.file("extdata", "obdp/Kt_trace.p", package="RevGadgets")
trees_trace_file <- 
     system.file("extdata", "obdp/mcmc_OBDP_trees.p", package="RevGadgets")
    
Kt_mean <- readOBDP( start_time_trace_file=start_time_trace_file, 
                     popSize_distribution_matrices_file=popSize_distribution_matrices_file, 
                     trees_trace_file=trees_trace_file )

# then get the customized ggplot object with plotDiversityOBDP()
p <- plotDiversityOBDP( Kt_mean,
                        xlab="Time (My)",
                        ylab="Number of lineages",
                        xticks_n_breaks=21,
                        col_Hidden="dodgerblue3",
                        col_LTT="gray25",
                        col_Total="forestgreen",
                        col_Hidden_interval="dodgerblue2",
                        col_Total_interval="darkolivegreen4",
                        palette_Hidden=c("transparent", "dodgerblue2", "dodgerblue3", 
                                         "dodgerblue4", "black"),
                        palette_Total=c("transparent", "green4", "forestgreen", "black"),
                        line_size=0.7,
                        interval_line_size=0.5,
                        show_Hidden=TRUE,
                        show_LTT=TRUE,
                        show_Total=TRUE,
                        show_intervals=TRUE,
                        show_densities=TRUE,
                        show_expectations=TRUE,
                        use_interpolate=TRUE )

# basic plot
p

# option: add a stratigraphic scale
library(deeptime)
library(ggplot2)
q <- gggeo_scale(p, dat="periods", height=unit(1.3, "line"), abbrv=F, size=4.5, neg=T)
r <- gggeo_scale(q, dat="epochs", height=unit(1.1, "line"), abbrv=F, size=3.5, neg=T, 
                    skip=c("Paleocene", "Pliocene", "Pleistocene", "Holocene"))
s <- gggeo_scale(r, dat="stages", height=unit(1, "line"), abbrv=T, size=2.5, neg=T)
s
}

}

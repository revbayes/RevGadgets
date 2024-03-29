% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processAncStates.R
\name{processAncStates}
\alias{processAncStates}
\title{Process Ancestral States}
\usage{
processAncStates(
  path,
  state_labels = NULL,
  labels_as_numbers = FALSE,
  missing_to_NA = TRUE
)
}
\arguments{
\item{path}{(character string; no default) File path to annotated tree.}

\item{state_labels}{(character vector; NULL) Vector of labels for ancestral
states named with the current state labels in annotated tree file
(as characters).}

\item{labels_as_numbers}{(logical; FALSE) Should the state labels be treated
as integers (for example, as chromosome numbers)?}

\item{missing_to_NA}{(logical; TRUE) Should missing data, coded as "?", be
coded to NA? If TRUE, the state will not be plotted. If FALSE, it will be
considered an additional state  when plotting.}
}
\value{
A treedata object
}
\description{
Process data for ancestral states plotting
}
\examples{

# standard ancestral state estimation example
file <- system.file("extdata",
                    "comp_method_disc/ase_freeK.tree",
                    package="RevGadgets")
example <- processAncStates(file,
                            state_labels = c("1" = "Awesome",
                                             "2" = "Beautiful",
                                             "3" = "Cool!"))

#chromosome evolution example
file <- system.file("extdata",
                    "chromo/ChromEvol_simple_final.tree",
                    package="RevGadgets")
chromo_example <- processAncStates(file, labels_as_numbers = TRUE)

}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeNetFromReference.R
\name{makeNetFromReference}
\alias{makeNetFromReference}
\title{Creates full regulatory network from reference files}
\usage{
makeNetFromReference(
  x = NULL,
  species = "Mouse",
  suffix = "_TF_Reference.tsv",
  reference_dir = "?",
  verbose = T
)
}
\arguments{
\item{x}{Use only these Transcription Factors. \code{NULL} (the default) creates
the biggest possible network.}

\item{species}{Organism. Respective reference file must exist.}

\item{suffix}{Reference file name is constructed by pasting \code{species} and \code{suffix}.}

\item{reference_dir}{Directory containing the reference file. \verb{?} (the default) uses
the reference directory of the package.}

\item{verbose}{Shows some information during process.}
}
\value{
A 2-column matrix representing directed edges of the network
}
\description{
Creates full regulatory network from reference files
}

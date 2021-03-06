% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/produce_reference.R
\name{produce_reference}
\alias{produce_reference}
\title{Produce a list of Transcription Factors interactions}
\usage{
produce_reference(
  t,
  tf,
  biomart,
  org_id = "",
  thresholdCols = c(8, 10),
  thresholds = c(),
  thresholdMethod = c("and", "or"),
  as.data.frame = F,
  verbose = T
)
}
\arguments{
\item{t}{Protein interactions from STRING database.}

\item{tf}{Transcription factors (gene names).}

\item{biomart}{Protein-Gene matching table. Normally from BioMart database.}

\item{org_id}{STRING organism ID. If specified, it is assumed that protein names
in the interactions table will have the form \emph{<org_id>.<protein>}.}

\item{thresholdCols}{Column indices in protein interactions dataset,
to be checked against thresholds}

\item{thresholds}{Respective thresholds for thresholdCols.
Only greater values than the threshold are accepted.}

\item{thresholdMethod}{The method to be used for threshold checking. See 'Details'.}

\item{as.data.frame}{Return data frame instead of list}

\item{verbose}{Print stats during processing}
}
\value{
List or data frame containing the TF interactions. NULL in case of error.
Format is: Count of interactions, TF, Interactors
}
\description{
Produces and returns a list of Transcription Factors interactions,
by combining provided information.
}
\details{
Two methods are supported for threshold checking. The default is \code{and}
where the thresholds are checked one after another, and the other is \code{or},
where each condition is checked separately and the results are combined.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgp.R
\name{ValuationSampleAnalytic}
\alias{ValuationSampleAnalytic}
\title{Data Generating Processes (DGPs)}
\usage{
ValuationSampleAnalytic(
  NumSamples,
  DistrType,
  Location = NULL,
  DistrScale = NULL,
  ShapeA = NULL,
  ShapeB = NULL
)
}
\arguments{
\item{NumSamples}{Number of valuations to sample.}

\item{DistrType}{The distribution type ('beta', 'gumbel', 'frechet', 'gev').}

\item{Location}{Distribution parameter (optional).}

\item{DistrScale}{Distribution parameter (optional).}

\item{ShapeA}{Distribution parameter (optional).}

\item{ShapeB}{Distribution parameter (optional).}
}
\value{
A vector of length NumSamples containing the sampled valuations.
}
\description{
Functions to generate valuation samples using static, time-varying, and left-digit discontinuity approaches.

Generates a sample from a static distribution of valuations.
}

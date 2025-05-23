% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expectedrewards.R
\name{GetCDF_Analytic}
\alias{GetCDF_Analytic}
\title{Generate CDF from Analytic Distribution}
\usage{
GetCDF_Analytic(
  DistrType,
  Location = NULL,
  DistrScale = NULL,
  ShapeA = NULL,
  ShapeB = NULL,
  granularity = 1e-06
)
}
\arguments{
\item{DistrType}{The distribution type. One of 'beta', 'gumbel', 'frechet', or 'gev'.}

\item{Location}{The location parameter for the distribution (optional).}

\item{DistrScale}{The scale parameter for the distribution (optional).}

\item{ShapeA}{Beta distribution parameter (optional).}

\item{ShapeB}{Beta distribution parameter (optional).}

\item{granularity}{The spacing between x-values in the CDF support (default is 1E-6).}
}
\value{
A list with two elements: \code{x}, the input values (support), and \code{y}, the corresponding CDF values.
}
\description{
Computes the cumulative distribution function (CDF) values for a specified distribution using analytic formulas.
}

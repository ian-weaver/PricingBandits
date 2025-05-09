% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expectedrewards.R
\name{ExpectedRewards}
\alias{ExpectedRewards}
\title{Expected Rewards Across Price Sets}
\usage{
ExpectedRewards(PriceSets, DistrClass, DistrType = NULL, ...)
}
\arguments{
\item{PriceSets}{A list of numeric vectors, where each vector contains prices to evaluate.}

\item{DistrClass}{The type of distribution ('analytic', 'empirical', 'leftdigit', 'timevarying').}

\item{DistrType}{The distribution type. One of 'beta', 'gumbel', 'frechet', or 'gev'.}

\item{...}{Additional parameters passed to the corresponding CDF-generating function based on the specified \code{DistrClass}.}
}
\value{
A list containing:
\describe{
  \item{\code{TrueOptimal}}{The maximum expected reward across all possible prices in the CDF.}
  \item{\code{Ratios}}{A named vector with the maximum reward in each price set divided by the true optimal reward.}
  \item{\code{Rewards}}{A list of expected rewards for each price set.}
}
}
\description{
Computes expected rewards for multiple sets of prices under a specified distribution class. Supports analytic, left-digit adjusted, and empirical distributions.
}

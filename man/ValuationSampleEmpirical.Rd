% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgp.R
\name{ValuationSampleEmpirical}
\alias{ValuationSampleEmpirical}
\title{Empirical Data Valuation Sample}
\usage{
ValuationSampleEmpirical(NumSamples, CDF, Scale)
}
\arguments{
\item{NumSamples}{Integer; the number of valuations to sample.}

\item{CDF}{A data frame representing the empirical CDF. The first column 
contains WTP, and the second column contains cumulative probabilities.}

\item{Scale}{Numeric; the scaling factor applied to the WTP data to ensure 
prices are normalized (e.g., between 0 and 1).}
}
\value{
A numeric vector of sampled valuations.
}
\description{
Generates samples from an empirical CDF of valuations.
}
\examples{
# Example CDF (prices and cumulative probabilities)
CDF <- data.frame(
  WTP = c(0.2, 0.5, 0.8),
  Prob = c(0.3, 0.6, 0.9)
)
}

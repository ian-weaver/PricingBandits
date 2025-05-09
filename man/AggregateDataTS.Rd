% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datapreprocessing.R
\name{AggregateDataTS}
\alias{AggregateDataTS}
\title{Aggregate Data for Thompson Sampling (TS) Algorithm}
\usage{
AggregateDataTS(PricesTested, PurchaseDecisions, TestX)
}
\arguments{
\item{PricesTested}{A vector of prices (values between 0 and 1) that have been tested so far in the experiment.}

\item{PurchaseDecisions}{A vector of associated purchase decisions for the prices tested.}

\item{TestX}{A set of prices (values between 0 and 1) to test.}
}
\value{
A list containing:
  \describe{
    \item{s_t}{The sum of purchase decisions for each price.}
    \item{n_t}{The number of observations for each price.}
  }
}
\description{
This function processes price testing data for the Thompson Sampling algorithm. It adds prior data
for untested prices and aggregates purchase decisions by price. Prices should be in the range [0, 1].
}
\examples{
PricesTested <- c(0.1, 0.2, 0.2, 0.9)
PurchaseDecisions <- c(1, 0, 1, 0)
TestX <- seq(10)/10
AggregateDataTS(PricesTested, PurchaseDecisions, TestX)
}

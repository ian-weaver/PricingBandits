% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datapostprocessing.R
\name{CalculateMaxRepeats}
\alias{CalculateMaxRepeats}
\title{Calculate Repeated Maximum Rewards Across Simulation Blocks}
\usage{
CalculateMaxRepeats(SimData, Prices, SeasonLength)
}
\arguments{
\item{SimData}{A data frame containing simulation results.}

\item{Prices}{A vector of prices in the test set.}

\item{SeasonLength}{Number of consumers per season; determines the frequency of seasonal demand shifts.}
}
\value{
A numeric vector containing the repeated maximum values for each block in each simulation.
}
\description{
For each row in a simulation matrix, this function finds the maximum reward 
within each block of prices and repeats it a fixed number of times (SeasonLength) to 
construct a vector of repeated maximums.
}

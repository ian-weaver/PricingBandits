% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gpregression.R
\name{NLML}
\alias{NLML}
\title{Gaussian Process Regression}
\usage{
NLML(Theta, TrainX, TrainY, sigma2_y)
}
\arguments{
\item{Theta}{Vector of hyperparameters (e.g., sigma_f and l).}

\item{TrainX}{Matrix of training points (m x d).}

\item{TrainY}{Vector of training targets (m x 1).}

\item{sigma2_y}{Observation noise variance (sigma_y^2).}
}
\value{
Negative log marginal likelihood value.
}
\description{
Functions to perform Gaussian Process regression, optimize hyperparameters, and compute posterior predictions.

Computes the negative log marginal likelihood for Gaussian Process regression.
}

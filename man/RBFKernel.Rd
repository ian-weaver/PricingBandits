% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernels.R
\name{RBFKernel}
\alias{RBFKernel}
\title{Kernel Functions}
\usage{
RBFKernel(x_i, x_j, sigma_f, l)
}
\arguments{
\item{x_i}{A point with d dimensions.}

\item{x_j}{A point with d dimensions.}

\item{sigma_f}{Hyperparameter defining the vertical scale.}

\item{l}{Hyperparameter defining the horizontal scale.}
}
\value{
Gaussian kernel value between two points.
}
\description{
Collection of functions to compute RBF kernels and covariance matrices.

Calculates the RBF kernel between two points.
}

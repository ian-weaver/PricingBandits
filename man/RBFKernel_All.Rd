% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernels.R
\name{RBFKernel_All}
\alias{RBFKernel_All}
\title{Generalized RBF Kernel}
\usage{
RBFKernel_All(x_i, x_j, d_i, d_j, sigma_f, l)
}
\arguments{
\item{x_i}{A point with d dimensions.}

\item{x_j}{A point with d dimensions.}

\item{d_i}{Order of derivative for point x_i (0 or 1).}

\item{d_j}{Order of derivative for point x_j (0 or 1).}

\item{sigma_f}{Hyperparameter defining the vertical scale.}

\item{l}{Hyperparameter defining the horizontal scale.}
}
\value{
Kernel value based on the derivative orders of the input points.
}
\description{
Computes RBF kernels for combinations of points and derivatives.
}

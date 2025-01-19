#' @title Kernel Functions
#' @description Collection of functions to compute RBF kernels and covariance matrices.

### Kernels ### ---------------------------------------------

#' @title Radial Basis Function (RBF) Kernel
#' @description Calculates the RBF kernel between two points.
#' @param x_i A point with d dimensions.
#' @param x_j A point with d dimensions.
#' @param sigma_f Hyperparameter defining the vertical scale.
#' @param l Hyperparameter defining the horizontal scale.
#' @return Gaussian kernel value between two points.
#' @export
RBFKernel <- function(x_i, x_j, sigma_f, l) {
  sigma_f^2 * exp(-(x_i - x_j)^2 / (2 * l^2))
}

#' @title RBF Kernel (Point to Derivative)
#' @description Calculates the RBF kernel between a point and derivative.
#' @param x_i A point with d dimensions.
#' @param x_j A point (derivative) with d dimensions.
#' @param sigma_f Hyperparameter defining the vertical scale.
#' @param l Hyperparameter defining the horizontal scale.
#' @return Kernel value between the point and derivative.
#' @export
RBFKernel_01 <- function(x_i, x_j, sigma_f, l) {
  sigma_f^2 / l^2 * (x_i - x_j) * exp(-(x_i - x_j)^2 / (2*l^2))
}

#' @title RBF Kernel (Derivative to Derivative)
#' @description Calculates the RBF kernel between two derivatives.
#' @param x_i A point (derivative) with d dimensions.
#' @param x_j A point (derivative) with d dimensions.
#' @param sigma_f Hyperparameter defining the vertical scale.
#' @param l Hyperparameter defining the horizontal scale.
#' @return Kernel value between the two derivatives.
#' @export
RBFKernel_11 <- function(x_i, x_j, sigma_f, l) {
  sigma_f^2 / l^4 * (l^2 - (x_i - x_j)^2) * exp(-(x_i - x_j)^2 / (2*l^2))
}

#' @title Generalized RBF Kernel
#' @description Computes RBF kernels for combinations of points and derivatives.
#' @param x_i A point with d dimensions.
#' @param x_j A point with d dimensions.
#' @param d_i Order of derivative for point x_i (0 or 1).
#' @param d_j Order of derivative for point x_j (0 or 1).
#' @param sigma_f Hyperparameter defining the vertical scale.
#' @param l Hyperparameter defining the horizontal scale.
#' @return Kernel value based on the derivative orders of the input points.
#' @export
RBFKernel_All <- function(x_i, x_j, d_i, d_j, sigma_f, l) {
  case_when(
    d_i == 0 & d_j == 0 ~ RBFKernel(x_i, x_j, sigma_f, l),
    d_i == 0 & d_j == 1 ~ RBFKernel_01(x_i, x_j, sigma_f, l),
    d_i == 1 & d_j == 0 ~ RBFKernel_01(x_j, x_i, sigma_f, l),
    d_i == 1 & d_j == 1 ~ RBFKernel_11(x_i, x_j, sigma_f, l)
  )
}

### Covariance Matrices ### ------------------------------

#' @title Covariance Matrix from Kernel
#' @description Computes the covariance matrix between two sets of points.
#' @param X1 Matrix of m points (m x d).
#' @param X2 Matrix of n points (n x d).
#' @param kernel Kernel function to compute covariance.
#' @param sigma_f Hyperparameter defining the vertical scale.
#' @param l Hyperparameter defining the horizontal scale.
#' @return Covariance matrix of size m x n.
#' @export
CovarianceFromKernel <- function(X1, X2, kernel, sigma_f, l) {
  outer(X1, X2, kernel, sigma_f, l)
}

#' @title Joint Covariance Matrix from Kernel
#' @description Computes the joint covariance matrix for a set of points and derivatives.
#' @param X Vector of points.
#' @param Index Vector indicating the derivative order at each point (0 or 1).
#' @param kernel Kernel function to compute covariance.
#' @param sigma_f Hyperparameter defining the vertical scale.
#' @param l Hyperparameter defining the horizontal scale.
#' @return Joint covariance matrix of size length(X) x length(X).
#' @export
JointCovFromKernel <- function(X, Index, kernel, sigma_f, l) {
  outer(1:length(X), 1:length(X),
        function(i, j) kernel(X[i], X[j], Index[i], Index[j], sigma_f, l))
}

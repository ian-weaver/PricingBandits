#' @title Kernel Functions
#' @description Functions to compute the RBF kernel (including its derivative
#' cross-covariances) and the covariance matrices built from it.

### Kernel ### ---------------------------------------------

#' @title Radial Basis Function (RBF) Kernel
#' @description Calculates the RBF kernel between two points, either of which may
#' represent the first derivative of the underlying function rather than a
#' function value. The derivative orders `d_i` and `d_j` (0 = function value,
#' 1 = first derivative; default 0) select the appropriate covariance:
#' value-value, value-derivative, derivative-value, or derivative-derivative.
#' @param x_i A point with d dimensions.
#' @param x_j A point with d dimensions.
#' @param sigma_f Hyperparameter defining the vertical scale.
#' @param l Hyperparameter defining the horizontal scale.
#' @param d_i Derivative order of `x_i` (0 or 1; default 0).
#' @param d_j Derivative order of `x_j` (0 or 1; default 0).
#' @return Kernel value between the two (possibly derivative) points.
#' @examples
#' RBFKernel(0.2, 0.5, 0.7, 0.2)                    # value-value
#' RBFKernel(0.2, 0.5, 0.7, 0.2, d_j = 1)           # value-derivative
#' RBFKernel(0.2, 0.5, 0.7, 0.2, d_i = 1, d_j = 1)  # derivative-derivative
#' @export
RBFKernel <- function(x_i, x_j, sigma_f, l, d_i = 0, d_j = 0) {
  # fast path for the common all-values case (e.g. likelihood evaluations)
  if (all(d_i == 0) && all(d_j == 0)) {
    return(sigma_f^2 * exp(-(x_i - x_j)^2 / (2 * l^2)))
  }
  case_when(
    d_i == 0 & d_j == 0 ~ sigma_f^2 * exp(-(x_i - x_j)^2 / (2 * l^2)),
    d_i == 0 & d_j == 1 ~ sigma_f^2 / l^2 * (x_i - x_j) * exp(-(x_i - x_j)^2 / (2*l^2)),
    d_i == 1 & d_j == 0 ~ sigma_f^2 / l^2 * (x_j - x_i) * exp(-(x_j - x_i)^2 / (2*l^2)),
    d_i == 1 & d_j == 1 ~ sigma_f^2 / l^4 * (l^2 - (x_i - x_j)^2) * exp(-(x_i - x_j)^2 / (2*l^2))
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
        function(i, j) kernel(X[i], X[j], sigma_f, l, Index[i], Index[j]))
}

#' @title Gaussian Process Regression
#' @description Functions to perform Gaussian Process regression, optimize hyperparameters, and compute posterior predictions.

### Obtain Optimal Hyperparameters ### ---------------------------------------------

#' @title Negative Log Marginal Likelihood (NLML)
#' @description Computes the negative log marginal likelihood for Gaussian Process regression.
#' @param Theta Vector of hyperparameters (e.g., sigma_f and l).
#' @param TrainX Matrix of training points (m x d).
#' @param TrainY Vector of training targets (m x 1).
#' @param sigma2_y Observation noise variance (sigma_y^2).
#' @return Negative log marginal likelihood value.
#' @export
NLML <- function(Theta, TrainX, TrainY, sigma2_y){
  K      = Matrix(CovarianceFromKernel(TrainX, TrainX, RBFKernel, Theta[1], Theta[2]) 
                  + (sigma2_y)*diag(length(TrainX)))
  Output = as.vector(0.5*determinant(K)$modulus[1] + TrainY %*% (solve(K) %*% TrainY)) 
                  + length(TrainX)*log(2*pi)
  return(Output)
}

#' @title Optimal Hyperparameters
#' @description Optimizes hyperparameters for Gaussian Process regression using NLopt.
#' @param TrainX Matrix of training points (m x d).
#' @param TrainY Vector of training targets (m x 1).
#' @param sigma2_y Observation noise variance (sigma_y^2).
#' @return Vector of optimized hyperparameters (sigma_f and l).
#' @export
OptimalHyperparameters <- function(TrainX, TrainY, sigma2_y){
  
  sigma_f    = 0.7
  l          = 0.2
  Parameters = c(sigma_f, l)
  Upper      = c(1.5, 0.4)
  Lower      = c(0.15, 0.1)
  
  tryCatch(
    expr = { 
      Results = nloptr(x0       = Parameters,
                       eval_f   = NLML,
                       lb       = Lower,
                       ub       = Upper,
                       opts     = list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel" = 1.0e-8),
                       TrainX   = TrainX,
                       TrainY   = TrainY,
                       sigma2_y = sigma2_y)
      return(Results$solution)
    },
    error = function(e){
      message("Warning: Error in optim (non-finite value supplied) - using priors instead")
      return(c(sigma_f, l))
    }
  )
}

### Obtain Mean and Covariance ### ---------------------------------------------

#' @title Posterior Prediction (Joint GP with Derivatives)
#' @description Computes the posterior mean and covariance for Gaussian Process regression with derivatives.
#' @param TrainX Vector of training points.
#' @param TrainY Vector of training targets.
#' @param TestX Vector of test points.
#' @param TestD Vector of test points for derivative prediction.
#' @param sigma_f Hyperparameter for kernel vertical scale.
#' @param l Hyperparameter for kernel horizontal scale.
#' @param sigma2_y Observation noise variance (sigma_y^2).
#' @return List containing posterior mean vector and covariance matrix.
#' @export
PosteriorPrediction <- function(TrainX, TrainY, TestX, TestD, sigma_f, l, sigma2_y){
  
  TrainIndex = rep(0, length(TrainX))
  TestIndex  = c(rep(0, length(TestX)), rep(1, length(TestD)))
  TestAll    = c(TestX, TestD)
  
  K_00 = JointCovFromKernel(TrainX, TrainIndex, RBFKernel_All, sigma_f, l)
  K_11 = JointCovFromKernel(TestAll, TestIndex, RBFKernel_All, sigma_f, l)
  K_01 = outer(1:length(TrainX), 1:length(TestAll),
               function(i, j) {RBFKernel_All(TrainX[i], TestAll[j], TrainIndex[i], TestIndex[j], sigma_f, l)})
  
  MeanPred = t(K_01) %*% (chol2inv(chol(K_00 + sigma2_y*diag(nrow(K_00))))) %*% TrainY
  CovPred  = K_11 - t(K_01) %*% (chol2inv(chol(K_00 + sigma2_y*diag(nrow(K_00))))) %*% K_01
  
  return (list(MeanPred, CovPred))
}

### Fix Improper Covariance Matrices ### ---------------------------------------------

#' @title Positive Definite Covariance Matrix
#' @description Ensures a covariance matrix is positive definite by adjusting negative eigenvalues.
#' @param CovMatrix Matrix to be corrected.
#' @param zilon Small value to replace negative eigenvalues.
#' @return Positive definite covariance matrix.
#' @export
MakePosDefinitive <- function(CovMatrix, zilon){
  
  NewMat   = as.matrix(forceSymmetric(CovMatrix))
  MinEigen = min(eigen(NewMat)$values)
  NegError = ifelse (MinEigen < zilon, TRUE, FALSE)
  
  while (NegError) {
    NewEig   = eigen(NewMat)
    NewEig2  = pmax(zilon, NewEig$values)
    NewMat   = NewEig$vectors %*% diag(NewEig2) %*% t(NewEig$vectors)
    NewMat   = as.matrix(forceSymmetric(NewMat))
    NegError = ifelse(min(eigen(NewMat)$values) < 0, TRUE, FALSE)
  }
  
  return(NewMat)
}

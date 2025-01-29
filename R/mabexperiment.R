#' @title Multi-Armed Bandit Experiment Framework
#' @description This file contains functions for running multi-armed bandit (MAB) experiments, supporting various policies, Gaussian Process-based methods, and monotonic variants. It allows simulation of consumer purchase decisions under different pricing strategies.

#' Run a Multi-Armed Bandit (MAB) Experiment
#'
#' This function simulates a Multi-Armed Bandit experiment using a specified policy and test settings.
#' It supports both Gaussian Process-based and monotonic variants, with optional heterogeneous noise modeling.
#'
#' @param Valuations A vector of consumer valuations for the product, used to simulate purchase decisions.
#' @param Policy The policy function to be used in the experiment (e.g., UCB, TS, GPTS, GPUCB).
#' @param TestX A vector of prices (values between 0 and 1) to test.
#' @param NumIter Number of iterations (time steps) to run the experiment.
#' @param BatchSize Number of consumers tested before the policy updates its action scores.
#' @param NumKnots Number of knots for monotonic Gaussian Process variants (optional, default is `NULL`).
#' @param Knots A vector of knot locations for monotonic Gaussian Process regression (optional, default is `NULL`).
#' @param BasisFunctions A matrix of basis functions for monotonic Gaussian Process regression (optional, default is `NULL`).
#' @param HeteroNoise Logical; if `TRUE`, allows for heterogeneous noise modeling (optional, default is `FALSE`).
#' @param Reset Number of consumers before the experiment history is wiped (optional, default is `NULL`).
#' @return A data frame containing two columns:
#'   \describe{
#'     \item{PricesTested}{A vector of prices tested over the experiment.}
#'     \item{PurchaseDecisions}{A vector of binary purchase decisions corresponding to each price tested.}
#'   }
#' @examples
#' Valuations <- runif(100, min = 0.1, max = 0.9)
#' TestX <- seq(10)/10
#' NumIter <- 100
#' BatchSize <- 10
#' MABExperiment(Valuations, UCB, TestX, NumIter, BatchSize)
#' @export
MABExperiment <- function(Valuations, Policy, TestX, NumIter, BatchSize,
                          NumKnots = NULL, Knots = NULL, BasisFunctions = NULL, 
                          HeteroNoise = FALSE, Reset = NULL) {
  
  # Setup for all GP Variants
  PricesTested      = rep(NA, NumIter)
  PurchaseDecisions = rep(NA, NumIter)
  Index2ActionHash  = hash(seq(length(TestX)), TestX)
  sigma2_y          = rep(0.25, length(TestX) + 1)
  ActionScores      = rep(0.5, length(TestX))
  GP                = any(sapply(c(GPTS, GPUCB, GPTS_Mono, GPUCB_Mono), identical, Policy))
  Mono              = any(sapply(c(GPTS_Mono, GPUCB_Mono), identical, Policy))
  
  # Additional setup for Monotonic Variants
  if (Mono){
    Knots          = c(0, seq(NumKnots - 1)/(NumKnots - 1))
    BasisFunctions = matrix(0, length(TestX), NumKnots)
    for(j in 1:length(TestX)){
      BasisFunctions[j, 1:NumKnots] = BasisFunction(TestX[j], NumKnots)
    }
  }
  
  # Run the Bandit
  for (t in 1:NumIter){
    
    # Update depending on batch size
    if ((t-1)%%BatchSize == 0){
      
      # Amend training data if algorithm is resetting
      if (!is.null(Reset)) {
        ToKeep = ifelse(t %% Reset == 1, 0, t %% Reset - 1)
        if (ToKeep == 0) {
          TrainX = TrainY = numeric(0)  
        } else {
          TrainX = na.omit(PricesTested[(t - ToKeep):(t-1)])
          TrainY = na.omit(PurchaseDecisions[(t - ToKeep):(t-1)])
        }
      } else {
        TrainX = na.omit(PricesTested)
        TrainY = na.omit(PurchaseDecisions)
      }
      
      # Pick first price randomly for GP variants
      if (length(TrainX) == 0 && GP){
      PriceToTest = sample(TestX)[1]
      } else {
        
        # If heterogeneous noise 
        if (HeteroNoise){
          sigma2_y = NoiseSample(TrainX, TrainY, TestX, BatchSize)
        }
        
        ActionScores = PolicyEvaluation(Policy, TrainX, TrainY, TestX, Knots, sigma2_y, 
                                        BatchSize, BasisFunctions, GP, Mono, ActionScores)
        a            = which.max(ActionScores)
        PriceToTest  = Index2ActionHash[[as.character(a)]]
      }
    }
    PurchaseDecision     = ifelse(PriceToTest < Valuations[t] + 1e-10, 1, 0)
    PricesTested[t]      = PriceToTest
    PurchaseDecisions[t] = PurchaseDecision
  }
  Output = data.frame(PricesTested, PurchaseDecisions)
  return(Output)
}
  
  
  
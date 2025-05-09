#' @title Multi-Armed Bandit Experiment Framework
#' @description Runs multi-armed bandit (MAB) experiments, supporting TS, UCB, Gaussian Process-based variants, and monotonic and heterogeneous variants. 
#' It allows simulation of consumer purchase decisions under different MAB policies and underlying WTP distributions.
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


#' @title Multi-Armed Bandit Experiments
#' @description Allows for multiple MAB experiment simulations to be run
#' @param Seeds A vector of the length of the number of simulations allowing reproducibility.
#' @param DGPType The type of dgp ('static', 'timevarying', 'leftdigit', 'empirical').
#' @param Policy The policy function to be used in the experiment (e.g., UCB, TS, GPTS, GPUCB).
#' @param TestX A vector of prices (values between 0 and 1) to test.
#' @param NumIter Number of iterations (time steps) to run the experiment.
#' @param BatchSize Number of consumers tested before the policy updates its action scores.
#' @param NumKnots Number of knots for monotonic Gaussian Process variants (monotonic variants).
#' @param Knots A vector of knot locations for monotonic Gaussian Process regression (monotonic variants).
#' @param BasisFunctions A matrix of basis functions for monotonic Gaussian Process regression (monotonic variants).
#' @param HeteroNoise Logical; if `TRUE`, allows for heterogeneous noise modeling (optional, default is `FALSE`).
#' @param Reset Number of consumers before the experiment history is wiped (optional - time varying).
#' @param DistrType The distribution type ('beta', 'gumbel', 'frechet', 'gev') (optional).
#' @param Location Distribution parameter (optional).
#' @param DistrScale Distribution parameter (optional).
#' @param ShapeA Distribution parameter (optional).
#' @param ShapeB Distribution parameter (optional).
#' @param SeasonLength Length of each season. (timevarying)
#' @param Rho Maximum seasonal shock (timevarying).
#' @param CDFGranularity Step size for generating the x-axis points of the CDF (leftdigit).
#' @param OneCentScaled Value of one cent when scaled to 0-1 range (leftdigit).
#' @param DP Vector of discontinuity points (leftdigit).
#' @param Zeta Scaling factor for the gap aka jump parameter (leftdigit).
#' @param CDF A data frame representing the empirical CDF. The first column 
#' contains WTP, and the second column contains cumulative probabilities (empirical).
#' @param Scale Numeric; the scaling factor applied to the WTP data to ensure 
#' prices are normalized (e.g., between 0 and 1) (empirical).
#' @return A data frame containing three columns:
#'   \describe{
#'     \item{PricesTested}{A vector of prices tested over the experiment.}
#'     \item{PurchaseDecisions}{A vector of binary purchase decisions corresponding to each price tested.}
#'     \item{SimulationID}{A vector containing the Simulation ID corresponding to the data}
#'   }
#' @export
MultipleMABExperiments <- function(Seeds, DGPType, Policy, TestX, NumIter, BatchSize,
                                   NumKnots = NULL, Knots = NULL, BasisFunctions = NULL, 
                                   HeteroNoise = FALSE, Reset = NULL, DistrType = NULL, Location = NULL, 
                                   DistrScale = NULL, ShapeA = NULL, ShapeB = NULL, SeasonLength = NULL, Rho = NULL, 
                                   CDFGranularity = NULL, OneCentScaled = NULL, DP = NULL, Zeta = NULL, 
                                   CDF = NULL, Scale = NULL){
  
  ResultsList = vector("list", length = length(Seeds))
  for (i in 1:length(Seeds)){
    set.seed(Seeds[i])
    Valuations = ValuationSample(DGPType, NumIter, DistrType, Location, DistrScale, ShapeA, ShapeB, 
                                 SeasonLength, Rho, CDFGranularity, OneCentScaled, DP, Zeta, CDF, Scale)
    ResultsList[[i]] = MABExperiment(Valuations, Policy, TestX, NumIter, BatchSize, NumKnots, Knots, 
                                     BasisFunctions, HeteroNoise, Reset)
    print(paste("Simulation", i, "complete."))
  }
  Results_df = do.call(rbind, ResultsList)
  Results_df$SimulationID = rep(1:length(Seeds), each = NumIter)
  return (Results_df)
}
  

  
  
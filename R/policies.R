#' @title Bandit Policies for Pricing Experiments
#' @description Collection of functions implementing various policies for pricing experiments, including Upper Confidence Bound (UCB), Thompson Sampling (TS), Gaussian Process-based policies, and their monotonic extensions.

### Baseline Policies ----------------------------------------------------------

#' Upper Confidence Bound (UCB) Policy
#'
#' This function implements the UCB policy for pricing experiments. It calculates action scores
#' based on aggregated data and selects prices to test using an upper confidence bound.
#'
#' @param PricesTested A vector of prices (values between 0 and 1) that have been tested so far in the experiment.
#' @param PurchaseDecisions A vector of associated purchase decisions for the prices tested.
#' @param TestX A set of prices (values between 0 and 1) to test.
#' @return A vector of action scores for the prices in `TestX`.
#' @examples
#' PricesTested <- c(0.1, 0.2, 0.2, 0.9)
#' PurchaseDecisions <- c(1, 0, 1, 0)
#' TestX <- seq(10)/10
#' UCB(PricesTested, PurchaseDecisions, TestX)
#' @export
UCB <- function(PricesTested, PurchaseDecisions, TestX){
  # Pre-processing
  AggregatedData = AggregateDataUCB(PricesTested, PurchaseDecisions, TestX)
  s_t            = AggregatedData[[1]]
  n_t            = AggregatedData[[2]]
  t              = sum(n_t)
  # Obtain Action Scores
  V_t          = (s_t*TestX^2)/n_t - (TestX*(s_t/n_t))^2 + (2*log(t)/n_t)^0.5
  ActionScores = TestX*((s_t/n_t) + ((log(t)/n_t) * pmin(rep(0.25, length(V_t)), V_t))^0.5)
  return(ActionScores)
}

#' Thompson Sampling (TS) Policy
#'
#' This function implements the Thompson Sampling policy for pricing experiments. It calculates
#' action scores by drawing from a Beta distribution based on aggregated data.
#'
#' @param PricesTested A vector of prices (values between 0 and 1) that have been tested so far in the experiment.
#' @param PurchaseDecisions A vector of associated purchase decisions for the prices tested.
#' @param TestX A set of prices (values between 0 and 1) to test.
#' @return A vector of action scores for the prices in `TestX`.
#' @examples
#' PricesTested <- c(0.1, 0.2, 0.2, 0.9)
#' PurchaseDecisions <- c(1, 0, 1, 0)
#' TestX <- seq(10)/10
#' TS(PricesTested, PurchaseDecisions, TestX)
#' @export
TS <- function(PricesTested, PurchaseDecisions, TestX){
  # Pre-processing
  AggregatedData = AggregateDataTS(PricesTested, PurchaseDecisions, TestX)
  s_t            = AggregatedData[[1]]
  n_t            = AggregatedData[[2]]
  # Obtain Action Scores
  RandomDraws  = unlist(Map(rbeta, 1, s_t, n_t - s_t))
  ActionScores = TestX*RandomDraws
  return(ActionScores)
}


### GP Policies ----------------------------------------------------------------

#' Gaussian Process Upper Confidence Bound (GPUCB) Policy
#'
#' This function implements the GPUCB policy for pricing experiments. It uses Gaussian Process regression
#' to calculate action scores based on posterior predictions and confidence intervals.
#'
#' @param PricesTested A vector of prices (values between 0 and 1) that have been tested so far in the experiment.
#' @param PurchaseDecisions A vector of associated purchase decisions for the prices tested.
#' @param TestX A set of prices (values between 0 and 1) to test.
#' @param sigma2_y Observation noise variance (must be less than 0.25).
#' @param BatchSize Number of consumers tested before an algorithm can update.
#' @return A vector of action scores for the prices in `TestX`.
#' @examples
#' PricesTested <- c(0.1, 0.2, 0.2, 0.9)
#' PurchaseDecisions <- c(1, 0, 1, 0)
#' TestX <- seq(10)/10
#' sigma2_y <- rep(0.25, 10)
#' BatchSize <- 10
#' GPUCB(PricesTested, PurchaseDecisions, TestX, sigma2_y, BatchSize)
#' @export
GPUCB <- function(PricesTested, PurchaseDecisions, TestX, sigma2_y, BatchSize){
  # Pre-processing
  AggregatedData = AggregateDataGP(PricesTested, PurchaseDecisions, TestX, sigma2_y, BatchSize)
  TrainX         = AggregatedData$PricesTested
  TrainY         = AggregatedData$PurchaseRates
  sigma2_y       = AggregatedData$sigma2_y
  OptParams      = OptimalHyperparameters(TrainX, TrainY, sigma2_y)
  t              = length(na.omit(PricesTested)) + BatchSize
  # Obtain Action Scores
  GP           = PosteriorPrediction(TrainX, TrainY, TestX, NULL, OptParams[1], OptParams[2], sigma2_y)
  GP[[2]]      = MakePosDefinitive(GP[[2]], 1e-10)
  SD           = abs(diag(GP[[2]]))^0.5
  D            = length(TestX)
  Delta        = 0.1
  Beta         = 2*log((D*t^2*pi^2)/(6*Delta))/5
  ActionScores = TestX*(GP[[1]] + sqrt(Beta)*SD)
  return(ActionScores)
}

#' Gaussian Process Thompson Sampling (GPTS) Policy
#'
#' This function implements the GPTS policy for pricing experiments. It uses Gaussian Process regression
#' to generate posterior predictions and action scores.
#'
#' @param PricesTested A vector of prices (values between 0 and 1) that have been tested so far in the experiment.
#' @param PurchaseDecisions A vector of associated purchase decisions for the prices tested.
#' @param TestX A set of prices (values between 0 and 1) to test.
#' @param sigma2_y Observation noise variance (must be less than 0.25).
#' @param BatchSize Number of consumers tested before an algorithm can update.
#' @return A vector of action scores for the prices in `TestX`.
#' @examples
#' PricesTested <- c(0.1, 0.2, 0.2, 0.9)
#' PurchaseDecisions <- c(1, 0, 1, 0)
#' TestX <- seq(10)/10
#' sigma2_y <- rep(0.25, 10)
#' BatchSize <- 10
#' GPTS(PricesTested, PurchaseDecisions, TestX, sigma2_y, BatchSize)
#' @export
GPTS <- function(PricesTested, PurchaseDecisions, TestX, sigma2_y, BatchSize){
  # Pre-processing
  AggregatedData = AggregateDataGP(PricesTested, PurchaseDecisions, TestX, sigma2_y, BatchSize)
  TrainX         = AggregatedData$PricesTested
  TrainY         = AggregatedData$PurchaseRates
  sigma2_y       = AggregatedData$sigma2_y
  OptParams      = OptimalHyperparameters(TrainX, TrainY, sigma2_y)
  # Obtain Action Scores
  GP           = PosteriorPrediction(TrainX, TrainY, TestX, NULL, OptParams[1], OptParams[2], sigma2_y)
  GP[[2]]      = MakePosDefinitive(GP[[2]], 1e-10)
  RandomDraws  = mvrnorm(1, GP[[1]], GP[[2]])
  ActionScores = RandomDraws*TestX
  return (ActionScores)
}

### GP-Mono policies -----------------------------------------------------------

#' Gaussian Process Upper Confidence Bound Monotonic (GPUCB_Mono) Policy
#'
#' This function implements the monotonic Gaussian Process Upper Confidence Bound (GPUCB_Mono) policy for pricing experiments.
#' It uses Gaussian Process regression to calculate action scores based on posterior predictions and confidence intervals 
#' while ensuring monotonicity.
#'
#' @param PricesTested A vector of prices (values between 0 and 1) that have been tested so far in the experiment.
#' @param PurchaseDecisions A vector of associated purchase decisions for the prices tested.
#' @param TestX A set of prices (values between 0 and 1) to test.
#' @param Knots Knots for the basis functions used in the regression.
#' @param sigma2_y Observation noise variance (must be less than 0.25).
#' @param BatchSize Number of consumers tested before an algorithm can update.
#' @param BasisFunctions Basis functions for the monotonic Gaussian Process regression.
#' @param LB1 Lower bound for the test prices.
#' @param LB2 Lower bound for the knots.
#' @param UB Upper bound for both test prices and knots.
#' @return A vector of action scores for the prices in `TestX`.
#' @export
GPUCB_Mono <- function(PricesTested, PurchaseDecisions, TestX, Knots, sigma2_y, BatchSize, 
                       BasisFunctions, LB1, LB2, UB){
  # Pre-processing
  AggregatedData = AggregateDataGP(PricesTested, PurchaseDecisions, TestX, sigma2_y, BatchSize)
  TrainX         = AggregatedData$PricesTested
  TrainY         = AggregatedData$PurchaseRates
  sigma2_y       = AggregatedData$sigma2_y
  OptParams      = OptimalHyperparameters(TrainX, TrainY, sigma2_y)
  t              = length(na.omit(PricesTested)) + BatchSize
  
  # Obtain Action Scores
  GP           = PosteriorPrediction(TrainX, TrainY, 0, Knots, OptParams[1], OptParams[2], sigma2_y)
  GP[[2]]      = MakePosDefinitive(GP[[2]], 1e-10)
  RandomDraws  = suppressWarnings(
    rtmvnorm(n     = 1E4,
             mu    = as.vector(GP[[1]]),
             sigma = forceSymmetric(GP[[2]]),
             lb    = c(LB1, rep(LB2, length(Knots))),
             ub    = c(UB, rep(0, length(Knots))))
  )
  DemandDraws  = matrix(0, nrow = nrow(RandomDraws), ncol = length(TestX))
  for (i in 1:nrow(RandomDraws)){
    DemandDraws[i, ] = as.vector(BasisFunctions %*% RandomDraws[i,][-1]) + RandomDraws[i,][1]
  }
  Mean         = apply(DemandDraws, 2, mean)
  SD           = apply(DemandDraws, 2, sd)
  D            = length(TestX)
  Delta        = 0.1
  Beta         = 2*log((D*t^2*pi^2)/(6*Delta))/5
  ActionScores = TestX*(Mean + sqrt(Beta)*SD)
  return(ActionScores)
}

#' Gaussian Process Thompson Sampling Monotonic (GPTS_Mono) Policy
#'
#' This function implements the monotonic Gaussian Process Thompson Sampling (GPTS_Mono) policy for pricing experiments.
#' It uses Gaussian Process regression to generate posterior predictions and action scores while ensuring monotonicity.
#'
#' @param PricesTested A vector of prices (values between 0 and 1) that have been tested so far in the experiment.
#' @param PurchaseDecisions A vector of associated purchase decisions for the prices tested.
#' @param TestX A set of prices (values between 0 and 1) to test.
#' @param Knots Knots for the basis functions used in the regression.
#' @param sigma2_y Observation noise variance (must be less than 0.25).
#' @param BatchSize Number of consumers tested before an algorithm can update.
#' @param BasisFunctions Basis functions for the monotonic Gaussian Process regression.
#' @param LB1 Lower bound for the test prices.
#' @param LB2 Lower bound for the knots.
#' @param UB Upper bound for both test prices and knots.
#' @return A vector of action scores for the prices in `TestX`.
#' @export
GPTS_Mono <- function(PricesTested, PurchaseDecisions, TestX, Knots, sigma2_y, BatchSize, 
                      BasisFunctions, LB1, LB2, UB){
  # Pre-processing
  AggregatedData = AggregateDataGP(PricesTested, PurchaseDecisions, TestX, sigma2_y, BatchSize)
  TrainX         = AggregatedData$PricesTested
  TrainY         = AggregatedData$PurchaseRates
  sigma2_y       = AggregatedData$sigma2_y
  OptParams      = OptimalHyperparameters(TrainX, TrainY, sigma2_y)
  
  # Obtain Action Scores
  GP           = PosteriorPrediction(TrainX, TrainY, 0, Knots, OptParams[1], OptParams[2], sigma2_y)
  GP[[2]]      = MakePosDefinitive(GP[[2]], 1e-10)
  RandomDraws  = suppressWarnings(
    rtmvnorm(n     = 1,
             mu    = as.vector(GP[[1]]),
             sigma = forceSymmetric(GP[[2]]),
             lb    = c(LB1, rep(LB2, length(Knots))),
             ub    = c(UB, rep(0, length(Knots))))
  )
  DemandDraw   = as.vector(BasisFunctions %*% RandomDraws[-1]) + RandomDraws[1]
  ActionScores = DemandDraw*TestX
  return(ActionScores)
}

### Wrappers for evaluating the policies ---------------------------------------

#' Evaluate Non-Monotonic Policies
#'
#' This function evaluates a non-monotonic policy by executing it with the provided parameters.
#' If the policy relies on Gaussian Processes (GP) and encounters an error, it falls back to using
#' the action scores from the previous round.
#'
#' @param Policy The policy function to evaluate.
#' @param PricesTested A vector of prices (values between 0 and 1) that have been tested so far in the experiment.
#' @param PurchaseDecisions A vector of associated purchase decisions for the prices tested.
#' @param TestX A set of prices (values between 0 and 1) to test.
#' @param sigma2_y Observation noise variance (must be less than 0.25).
#' @param BatchSize Number of consumers tested before an algorithm can update.
#' @param GP Logical; if TRUE, the policy relies on Gaussian Processes.
#' @param PrevAS A vector of previous action scores to use as a fallback in case of errors.
#' @return A vector of action scores for the prices in `TestX`.
#' @export
NonMonoPolicyEval <- function(Policy, PricesTested, PurchaseDecisions, TestX, sigma2_y, BatchSize, GP, PrevAS) {
  
  # Try to evaluate the policy
  ActionScores <- tryCatch(
    expr = {
      if (GP) {
        Policy(PricesTested, PurchaseDecisions, TestX, sigma2_y, BatchSize)
      } else {
        Policy(PricesTested, PurchaseDecisions, TestX)
      }
    },
    error = function(e) {
      message("Warning: Error in GP - use last round's information")
      PrevAS
    }
  )
  
  return(ActionScores)
}


#' Evaluate Policies for Pricing Experiments
#'
#' This function evaluates either monotonic or non-monotonic policies by executing them with
#' the provided parameters. For monotonic policies, it includes mechanisms to handle timeout
#' scenarios and fallback options.
#'
#' @param Policy The policy function to evaluate.
#' @param PricesTested A vector of prices (values between 0 and 1) that have been tested so far in the experiment.
#' @param PurchaseDecisions A vector of associated purchase decisions for the prices tested.
#' @param TestX A set of prices (values between 0 and 1) to test.
#' @param Knots Knots for the basis functions used in monotonic policies.
#' @param sigma2_y Observation noise variance (must be less than 0.25).
#' @param BatchSize Number of consumers tested before an algorithm can update.
#' @param BasisFunctions Basis functions for monotonic Gaussian Process regression.
#' @param GP Logical; if TRUE, the policy relies on Gaussian Processes.
#' @param Mono Logical; if TRUE, evaluates a monotonic policy.
#' @param PrevAS A vector of previous action scores to use as a fallback in case of errors.
#' @return A vector of action scores for the prices in `TestX`.
#' @export
PolicyEvaluation <- function(Policy, PricesTested, PurchaseDecisions, TestX, Knots, sigma2_y, 
                             BatchSize, BasisFunctions, GP, Mono, PrevAS) {
  
  if (!Mono) {
    # Non-Monotonic Policy Evaluation
    ActionScores <- NonMonoPolicyEval(Policy, PricesTested, PurchaseDecisions, TestX, sigma2_y, BatchSize, GP, PrevAS)
  } else {
    # Monotonic Policy Evaluation
    ActionScores <- tryCatch(
      expr = {
        # First attempt with original settings and a timeout
        withTimeout({
          Policy(PricesTested, PurchaseDecisions, TestX, Knots, sigma2_y, BatchSize, 
                 BasisFunctions, LB1 = -Inf, LB2 = -Inf, UB = 2)
        }, timeout = 5)
      },
      TimeoutException = function(ex) {
        # Second attempt with modified bounds and a timeout
        tryCatch(
          expr = {
            withTimeout({
              Policy(PricesTested, PurchaseDecisions, TestX, Knots, sigma2_y, BatchSize, 
                     BasisFunctions, LB1 = -Inf, LB2 = -Inf, UB = 2)
            }, timeout = 5)
          },
          error = function(e) {
            # Last resort: use Non-Monotonic method
            NonMonoPolicyEval(Policy, PricesTested, PurchaseDecisions, TestX, sigma2_y, BatchSize, GP, PrevAS)
          }
        )
      }
    )
  }
  return(ActionScores)
}


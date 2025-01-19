#' @title Noise Sampling for Heteroscedastic Gaussian Processes
#' @description This function generates a sample of heteroscedastic noise variance (\eqn{\sigma^2_y}) for a Gaussian Process.
#' The function processes price testing data and uses posterior predictions to generate a random draw.
#' Prices should be in the range [0, 1], and observation noise variance (\eqn{\sigma^2_y}) must be less than 0.25.
#'
#' @param PricesTested A vector of prices (values between 0 and 1) that have been tested so far in the experiment.
#' @param PurchaseDecisions A vector of associated purchase decisions for the prices tested.
#' @param TestX A set of prices (values between 0 and 1) to test.
#' @param BatchSize Number of consumers tested before an algorithm can update.
#' @return A vector containing the sampled heteroscedastic noise variances (\eqn{\sigma^2_y}).
#' @examples
#' PricesTested <- c(0.1, 0.2, 0.2, 0.9)
#' PurchaseDecisions <- c(1, 0, 1, 0)
#' TestX <- seq(10)/10
#' sigma2_y <- rep(0.25, 10)
#' BatchSize <- 10
#' NoiseSample(PricesTested, PurchaseDecisions, TestX, BatchSize)
#' @export
NoiseSample <- function(PricesTested, PurchaseDecisions, TestX, BatchSize){
  # Pre-processing
  sigma2_y       = rep(0.25, length(TestX) + 1)
  AggregatedData = AggregateDataGP(PricesTested, PurchaseDecisions, TestX, sigma2_y, BatchSize)
  TrainX         = AggregatedData$PricesTested
  TrainY         = AggregatedData$PurchaseRates
  sigma2_y       = AggregatedData$sigma2_y
  OptParams      = OptimalHyperparameters(TrainX, TrainY, sigma2_y)
  # Obtain sigma2_y sample
  GP           = PosteriorPrediction(TrainX, TrainY, c(0, TestX), NULL, OptParams[1], OptParams[2], sigma2_y)
  GP[[2]]      = MakePosDefinitive(GP[[2]], 1e-10)
  RandomDraws  = mvrnorm(1, GP[[1]], GP[[2]])
  RandomDraws[RandomDraws > 0.99] = 0.99
  RandomDraws[RandomDraws < 0.01] = 0.01
  sigma2_y     = RandomDraws*(1 - RandomDraws)
  return (sigma2_y)
}






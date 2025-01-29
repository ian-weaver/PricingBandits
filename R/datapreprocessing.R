#' @title GP Variants
#' @description Collection of functions to aggregate data for various algorithms, including UCB, Thompson Sampling, and Gaussian Processes.

#' Aggregate Data for Upper Confidence Bound (UCB) Algorithm
#'
#' This function processes price testing data for the UCB algorithm. It aggregates purchase decisions
#' by price, includes priors for prices that have not yet been tested, and returns summarized data.
#' Prices should be in the range [0, 1].
#'
#' @param PricesTested A vector of prices (values between 0 and 1) that have been tested so far in the experiment.
#' @param PurchaseDecisions A vector of associated purchase decisions for the prices tested.
#' @param TestX A set of prices (values between 0 and 1) to test.
#' @return A list containing:
#'   \describe{
#'     \item{s_t}{The sum of purchase decisions for each price.}
#'     \item{n_t}{The number of observations for each price.}
#'   }
#' @examples
#' PricesTested <- c(0.1, 0.2, 0.2, 0.9)
#' PurchaseDecisions <- c(1, 0, 1, 0)
#' TestX <- seq(10)/10
#' AggregateDataUCB(PricesTested, PurchaseDecisions, TestX)
#' @export
AggregateDataUCB <- function(PricesTested, PurchaseDecisions, TestX) {
  
  RawData = data.frame(PricesTested, PurchaseDecisions)
  SumData = RawData %>% 
    group_by(PricesTested) %>% 
    summarise(s_t = sum(PurchaseDecisions), n_t = n(), .groups = 'drop')
  
  MissingData = data.frame(PricesTested = TestX, s_t = 1, n_t = 1)
  Result = MissingData %>% 
    left_join(SumData, by = 'PricesTested', suffix = c("", ".y")) %>% 
    mutate(s_t = coalesce(s_t.y, s_t), 
           n_t = coalesce(n_t.y, n_t)) %>% 
    dplyr::select(-s_t.y, -n_t.y)
  
  return(list(s_t = Result$s_t, n_t = Result$n_t))
}

#' Aggregate Data for Thompson Sampling (TS) Algorithm
#'
#' This function processes price testing data for the Thompson Sampling algorithm. It adds prior data
#' for untested prices and aggregates purchase decisions by price. Prices should be in the range [0, 1].
#'
#' @param PricesTested A vector of prices (values between 0 and 1) that have been tested so far in the experiment.
#' @param PurchaseDecisions A vector of associated purchase decisions for the prices tested.
#' @param TestX A set of prices (values between 0 and 1) to test.
#' @return A list containing:
#'   \describe{
#'     \item{s_t}{The sum of purchase decisions for each price.}
#'     \item{n_t}{The number of observations for each price.}
#'   }
#' @examples
#' PricesTested <- c(0.1, 0.2, 0.2, 0.9)
#' PurchaseDecisions <- c(1, 0, 1, 0)
#' TestX <- seq(10)/10
#' AggregateDataTS(PricesTested, PurchaseDecisions, TestX)
#' @export
AggregateDataTS <- function(PricesTested, PurchaseDecisions, TestX){
  
  PricesTested       = c(PricesTested, rep(TestX, 2))
  PurchaseDecisions  = c(PurchaseDecisions, rep(0, length(TestX)), rep(1, length(TestX)))
  
  RawData = data.frame(PricesTested, PurchaseDecisions)
  SumData = RawData %>% 
    group_by(PricesTested) %>% 
    summarise(s_t = sum(PurchaseDecisions), n_t = n(), .groups = 'drop')
  
  return(list(s_t = SumData$s_t, n_t = SumData$n_t))
}

#' Aggregate Data for Gaussian Process (GP) Algorithms
#'
#' This function processes price testing data for Gaussian Process algorithms. It computes purchase rates
#' and adjusts observation noise variance for each price. Prices should be in the range [0, 1], and
#' observation noise variance (sigma2_y) must be less than 0.25.
#'
#' @param PricesTested A vector of prices (values between 0 and 1) that have been tested so far in the experiment.
#' @param PurchaseDecisions A vector of associated purchase decisions for the prices tested.
#' @param TestX A set of prices (values between 0 and 1) to test.
#' @param sigma2_y Observation noise variance (sigma_y^2), must be less than 0.25.
#' @param BatchSize Number of consumers tested before an algorithm can update.
#' @return A data frame containing:
#'   \describe{
#'     \item{PricesTested}{The prices included in the dataset.}
#'     \item{PurchaseRates}{The mean purchase rate for each price.}
#'     \item{sigma_2y}{The adjusted observation noise variance for each price.}
#'   }
#' @examples
#' PricesTested <- c(0.1, 0.2, 0.2, 0.9)
#' PurchaseDecisions <- c(1, 0, 1, 0)
#' TestX <- seq(10)/10
#' sigma2_y <- rep(0.25, 10)
#' BatchSize <- 10
#' AggregateDataGP(PricesTested, PurchaseDecisions, TestX, sigma2_y, BatchSize)
#' @export
AggregateDataGP <- function(PricesTested, PurchaseDecisions, TestX, sigma2_y, BatchSize){
  
  PricesTested      = c(PricesTested, rep(0, BatchSize))
  PurchaseDecisions = c(PurchaseDecisions, rep(1, BatchSize))
  
  RawData  = data.frame(PricesTested, PurchaseDecisions)
  MeanData = RawData %>% 
    group_by(PricesTested) %>% 
    summarise(PurchaseRates = mean(PurchaseDecisions), .groups = 'drop')
  
  n_t               = as.vector(table(PricesTested))
  sigma2_y          = sigma2_y[which(round(c(0, TestX), 2) %in% round(MeanData$PricesTested, 2))]
  MeanData$sigma2_y = sigma2_y/n_t
  return (MeanData)
}


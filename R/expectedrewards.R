### GET CDF ### ----------------------------------------------------------------

#' @title Generate CDF from Analytic Distribution
#' @description Computes the cumulative distribution function (CDF) values for a specified distribution using analytic formulas.
#' @param DistrType The distribution type. One of 'beta', 'gumbel', 'frechet', or 'gev'.
#' @param Location The location parameter for the distribution (optional).
#' @param DistrScale The scale parameter for the distribution (optional).
#' @param ShapeA Beta distribution parameter (optional).
#' @param ShapeB Beta distribution parameter (optional).
#' @param granularity The spacing between x-values in the CDF support (default is 1E-6).
#' @return A list with two elements: \code{x}, the input values (support), and \code{y}, the corresponding CDF values.
#' @export
GetCDF_Analytic <- function(DistrType, Location = NULL, DistrScale = NULL, ShapeA = NULL, ShapeB = NULL, granularity = 1E-6) {
  CDF_x <- seq(granularity, 1, by = granularity)
  CDF_y <- switch(DistrType,
                  beta    = pbeta(CDF_x, ShapeA, ShapeB),
                  gumbel  = pgumbel(CDF_x, Location, DistrScale),
                  frechet = pfrechet(CDF_x, Location, DistrScale, ShapeA),
                  gev     = pgev(CDF_x, Location, DistrScale, ShapeA))
  list(x = CDF_x, y = CDF_y)
}

#' @title Left-Digit CDF Generator
#' @description Computes the cumulative distribution function (CDF) for a left-digit adjusted Beta distribution, incorporating valuation discontinuities.
#' @param CDFGranularity Step size for generating the x-axis points of the CDF.
#' @param OneCentScaled Value of one cent when scaled to the 0–1 range.
#' @param ShapeA Beta distribution parameter.
#' @param ShapeB Beta distribution parameter.
#' @param DP Vector of discontinuity points.
#' @param Zeta Scaling factor for the discontinuity gap (jump parameter).
#' @return A list with two elements: \code{x}, the support points for the CDF, and \code{y}, the corresponding CDF values adjusted for left-digit discontinuities.
#' @export
GetCDF_LeftDigit <- function(CDFGranularity, OneCentScaled, ShapeA, ShapeB, DP, Zeta) {
  DG     <- ObtainDiscontinuityGaps(OneCentScaled, ShapeA, ShapeB, DP, Zeta)
  CDF_x  <- seq(0, 1, CDFGranularity)
  CDF_y  <- LeftDigitCDF(CDFGranularity, DP, DG, ShapeA, ShapeB)
  list(x = CDF_x, y = CDF_y)
}

#' @title Rescale Empirical CDF
#' @description Takes an empirical CDF and rescales the x-axis to fit within the 0–1 range.
#' @param PreScaledCDF A two-column matrix or data frame where the first column contains original x-values (e.g., prices) and the second column contains corresponding CDF values.
#' @param Scale A scalar used to rescale the x-axis values (e.g., 0.001 to convert prices up to 1000 into the [0,1] range).
#' @return A list with two elements: \code{x}, the rescaled support, and \code{y}, the original CDF values.
#' @export
GetCDF_Empirical <- function(PreScaledCDF, Scale) {
  list(x = PreScaledCDF[,1] * Scale, y = PreScaledCDF[,2])
}

### Expected Rewards ### -------------------------------------------------------

#' @title Expected Rewards from CDF
#' @description Computes the expected reward at each price using a given CDF. The expected reward is calculated as the price multiplied by the probability of purchase (1 - CDF at that price).
#' @param Prices A vector of prices in the test set.
#' @param CDF_x Vector of x-values (support) for the CDF.
#' @param CDF_y Vector of CDF values corresponding to \code{CDF_x}.
#' @return A numeric vector of expected rewards, one for each price in \code{Prices}.
#' @export
ComputeExpectedRewardsFromCDF <- function(Prices, CDF_x, CDF_y) {
  RoundPlaces = ceiling(log10(length(CDF_x)))
  Indices     = match(round(as.numeric(as.character(Prices)), RoundPlaces), round(as.numeric(as.character(CDF_x)), RoundPlaces))
  Prices * (1 - CDF_y[Indices])
}

#' @title Expected Rewards Across Price Sets
#' @description Computes expected rewards for multiple sets of prices under a specified distribution class. Supports analytic, left-digit adjusted, and empirical distributions.
#' @param PriceSets A list of numeric vectors, where each vector contains prices to evaluate.
#' @param DistrClass The type of distribution ('analytic', 'empirical', 'leftdigit', 'timevarying').
#' @param DistrType The distribution type. One of 'beta', 'gumbel', 'frechet', or 'gev'.
#' @param ... Additional parameters passed to the corresponding CDF-generating function based on the specified \code{DistrClass}.
#' @return A list containing:
#' \describe{
#'   \item{\code{TrueOptimal}}{The maximum expected reward across all possible prices in the CDF.}
#'   \item{\code{Ratios}}{A named vector with the maximum reward in each price set divided by the true optimal reward.}
#'   \item{\code{Rewards}}{A list of expected rewards for each price set.}
#' }
#' @export
ExpectedRewards <- function(PriceSets, DistrClass, DistrType = NULL, ...) {
  
  CDF <- switch(DistrClass,
                analytic  = GetCDF_Analytic(DistrType, ...),
                leftdigit = GetCDF_LeftDigit(...),
                empirical = GetCDF_Empirical(...)
  )
  
  # Compute expected rewards for each price set
  ExRewards = lapply(PriceSets, function(prices) {
    ComputeExpectedRewardsFromCDF(prices, CDF$x, CDF$y)
  })
  names(ExRewards) = paste0("A", lengths(PriceSets))
  
  # Compute the true optimal
  TrueOptimal = max(CDF$x * (1 - CDF$y))
  
  # Ratios relative to true optimal
  Ratios = sapply(ExRewards, function(reward_vec) max(reward_vec) / TrueOptimal)
  
  list(TrueOptimal = TrueOptimal, Ratios = Ratios, Rewards = ExRewards)
}

### TIME-VARYING CASE ### ------------------------------------------------------

#' @title Expected Rewards in Time-Varying Environment
#' @description Computes expected rewards over multiple pricing periods in a time-varying demand setting, where consumer valuations shift according to a seasonal shock process.
#' @param Prices A vector of prices in the test set.
#' @param Seeds A vector of random seeds used to generate different demand realizations.
#' @param NumIter Integer. Number of iterations (consumers) in each simulation run.
#' @param SeasonLength Number of consumers per season; determines the frequency of seasonal demand shifts.
#' @param Rho Maximum seasonal shock; actual shocks are drawn uniformly from \code{[-Rho, Rho]}.
#' @param DistrType Distribution type for baseline valuations (e.g., 'beta').
#' @param ShapeA Beta distribution parameter.
#' @param ShapeB Beta distribution parameter.
#' @return A matrix of expected rewards with rows corresponding to each seed and columns corresponding to price–season combinations.
#' @export
ExpectedRewardsTV <- function(Prices, Seeds, NumIter, SeasonLength, Rho, DistrType, ShapeA, ShapeB){
  ExpectedRewards = matrix(0, nrow=length(Seeds), ncol=length(Prices)*(NumIter/SeasonLength))
  for (i in Seeds){
    set.seed(i)
    Samples = ValuationSampleTimeVarying(NumIter, SeasonLength, Rho, DistrType, ShapeA, ShapeB)
    Valuations = Samples[[1]]
    Shocks     = unique(Samples[[2]])
    ExRewards_i = c()
    for (j in Shocks){
      Demands = (1 - pbeta(Prices - j, ShapeA, ShapeB))
      Rewards = Prices*Demands
      ExRewards_i = c(ExRewards_i, Rewards)
    }
    ExpectedRewards[i, ] = ExRewards_i
  }
  return(ExpectedRewards)
}
#' @title Data Generating Processes (DGPs)
#' @description Functions to generate valuation samples using static, time-varying, and left-digit discontinuity approaches.
#' @import dplyr Matrix
#' @importFrom stats na.omit pbeta rbeta runif sd
#' @importFrom evd rgumbel rfrechet rgev
#' @importFrom nloptr nloptr
#' @importFrom MASS mvrnorm
#' @importFrom TruncatedNormal rtmvnorm
#' @importFrom R.utils withTimeout
#' @importFrom hash hash

### STATIC CASE ### ---------------------------------------------------

#' @title Static Valuation Sample
#' @description Generates a sample from a static distribution of valuations.
#' @param NumSamples Number of valuations to sample.
#' @param Type The distribution type ('beta', 'gumbel', 'frechet', 'gev').
#' @param Location Distribution parameter (optional).
#' @param DistrScale Distribution parameter (optional).
#' @param ShapeA Distribution parameter (optional).
#' @param ShapeB Distribution parameter (optional).
#' @return A vector of length NumSamples containing the sampled valuations.
#' @export
ValuationSampleStatic <- function(NumSamples, Type, Location = NULL, DistrScale = NULL, 
                                  ShapeA = NULL, ShapeB = NULL){
  switch(Type,
         beta    = rbeta(NumSamples, ShapeA, ShapeB),
         gumbel  = rgumbel(NumSamples, Location, DistrScale),
         frechet = rfrechet(NumSamples, Location, DistrScale, ShapeA),
         gev     = rgev(NumSamples, Location, DistrScale, ShapeA))
}


### TIME-VARYING CASE ### --------------------------------------------

#' @title Time-Varying Valuation Sample
#' @description Generates a time-varying valuation sample by introducing seasonal shocks.
#' @param NumSamples Number of valuations to sample.
#' @param SeasonLength Length of each season.
#' @param Rho Maximum seasonal shock (uniform distribution between -Rho and Rho).
#' @param Type The distribution type ('beta', 'gumbel', 'frechet', 'weibull', 'gev').
#' @param ShapeA Distribution parameter (optional).
#' @param ShapeB Distribution parameter (optional).
#' @param Location Distribution parameter (optional).
#' @param DistrScale Distribution parameter (optional).
#' @return A list containing seasonal valuations and associated shocks.
#' @export
ValuationSampleTimeVarying <- function(NumSamples, SeasonLength, Rho, Type, ShapeA = NULL, 
                                       ShapeB = NULL, Location = NULL, DistrScale = NULL){
  StaticVals     = ValuationSampleStatic(NumSamples, Type, Location, DistrScale, ShapeA, ShapeB)
  Seasons        = unname(lengths(split(StaticVals, rep(1:ceiling(length(StaticVals)/SeasonLength), 
                                                        each = SeasonLength, length.out = length(StaticVals)))))
  SeasonalShocks = runif(length(Seasons), min = -Rho, max = Rho)
  Shocks         = rep(SeasonalShocks, Seasons)
  SeasonalVals   = StaticVals + Shocks
  return (list(SeasonalVals, Shocks))
}

### LEFT DIGIT CASE ### ----------------------------------------------

#' @title Discontinuity Gap Calculation
#' @description Calculates discontinuity gaps in demand at specified price points.
#' @param OneCentScaled Value of one cent when scaled to 0-1 range.
#' @param ShapeA Beta distribution parameter.
#' @param ShapeB Beta distribution parameter.
#' @param DP Vector of discontinuity points.
#' @param Zeta Scaling factor for the gap (jump parameter).
#' @return A vector of discontinuity gaps.
#' @export
ObtainDiscontinuityGaps <- function(OneCentScaled, ShapeA, ShapeB, DP, Zeta) {
  Prices                = seq(0, 1 - OneCentScaled, OneCentScaled)
  ContinuousDemandCurve = 1 - pbeta(Prices, ShapeA, ShapeB)
  UsualGap              = numeric(length(DP))
  
  for (i in seq_along(DP)) {
    IndexBefore = which.min(abs(Prices - (DP[i] - OneCentScaled)))
    IndexAfter  = which.min(abs(Prices - DP[i]))
    UsualGap[i] = ContinuousDemandCurve[IndexBefore] - ContinuousDemandCurve[IndexAfter]
  }
  
  DG = UsualGap * Zeta
  return (DG)
}


#' @title Left-Digit CDF
#' @description Creates a CDF with discontinuity gaps at specified points.
#' @param CDFGranularity Step size for generating the x-axis points of the CDF.
#' @param DP Vector of discontinuity points.
#' @param DG Vector of discontinuity gaps.
#' @param ShapeA Beta distribution parameter.
#' @param ShapeB Beta distribution parameter.
#' @return A vector of discrete CDF values for the left-digit case.
#' @export
LeftDigitCDF <- function(CDFGranularity, DP, DG, ShapeA, ShapeB){
  x    = seq(0, 1, CDFGranularity)
  z    = 1e-5
  DP   = c(0, DP, 1)
  DG   = c(0, DG)
  LCDF = rep(0, length(x))
  for (i in 1:(length(DP)-1)){
    LCDF[x > (DP[i]-z) & x < (DP[i+1]-z)] <- pbeta(x[x > (DP[i]-z) & x < (DP[i+1]-z)], ShapeA, ShapeB)*(1-sum(DG)) + sum(DG[1:i])
  }
  LCDF[length(x)] = 1
  return(LCDF)
}

#' @title Inverse CDF
#' @description Computes the inverse of the CDF at a given point.
#' @param u Point at which to find the inverse CDF.
#' @param CDF_x Vector of x-axis points of the CDF.
#' @param CDF_y Vector of the CDF values.
#' @return The inverse CDF value at point u.
#' @export
Inverse <- function(u, CDF_x, CDF_y){
  Index  = which(CDF_y > u)[1] - 1
  Sample = CDF_x[Index]
  return (Sample)
}

#' @title Left-Digit Valuation Sample
#' @description Generates samples from a left-digit adjusted distribution of valuations.
#' @param NumSamples Number of valuations to sample.
#' @param CDFGranularity Step size for generating the x-axis points of the CDF.
#' @param OneCentScaled Value of one cent when scaled to 0-1 range.
#' @param ShapeA Beta distribution parameter.
#' @param ShapeB Beta distribution parameter.
#' @param DP Vector of discontinuity points.
#' @param Zeta Scaling factor for the gap (jump parameter).
#' @return A vector of sampled valuations.
#' @export
ValuationSampleLD <- function(NumSamples, CDFGranularity, OneCentScaled, ShapeA, ShapeB, DP, Zeta){
  DG         = ObtainDiscontinuityGaps(OneCentScaled, ShapeA, ShapeB, DP, Zeta)
  CDF_x      = seq(0, 1, CDFGranularity)
  CDF_y      = LeftDigitCDF(CDFGranularity, DP, DG, ShapeA, ShapeB)
  Draws      = runif(NumSamples)
  Valuations = sapply(Draws, Inverse, CDF_x, CDF_y)
  return(Valuations)
}

### Field Data ### ----------------------------------------------

#' @title Field Data Valuation Sample
#' @description Generates samples from an empirical CDF of valuations.
#' @param NumSamples Integer; the number of valuations to sample.
#' @param CDF A data frame representing the empirical CDF. The first column 
#' contains WTP, and the second column contains cumulative probabilities.
#' @param Scale Numeric; the scaling factor applied to the WTP data to ensure 
#' prices are normalized (e.g., between 0 and 1).
#' @return A numeric vector of sampled valuations.
#' @examples
#' # Example CDF (prices and cumulative probabilities)
#' CDF <- data.frame(
#'   WTP = c(0.2, 0.5, 0.8),
#'   Prob = c(0.3, 0.6, 0.9)
#' )
#' @export
ValuationSampleField <- function(NumSamples, CDF, Scale){
  CDF[,1] = CDF[,1]*Scale
  # Complete CDF bounds (can only test prices within original empirical range)
  if (min(CDF[,2]) > 0){
    CDF = rbind(c(min(CDF[,1] - 0.01), 0), CDF) 
  }
  if (max(CDF[,2]) < 1){
    CDF = rbind(CDF, c(max(CDF[,1] + 0.01), 1))
  }
  Draws      = runif(NumSamples)
  Valuations = vapply(Draws, Inverse, CDF_x = CDF[,1], CDF_y = CDF[,2], FUN.VALUE = numeric(1))
  return(Valuations)
}


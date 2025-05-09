### Process Expected Rewards from Experiments ### ------------------------------

#' @title Prepare Simulation Data with Expected Rewards
#' @description
#' Helper function to process simulation data by computing expected rewards from tested prices,
#' assigning simulation IDs, and determining the optimal reward.
#' @param SimData A data frame containing simulation results, with a \code{PricesTested} column.
#' @param ExRewards A vector of expected rewards for each price in \code{Prices}.
#' @param Prices A numeric vector of all candidate prices.
#' @param NumIter Integer. Number of iterations per simulation.
#' @param NumSimulations Integer. Total number of simulations.
#' @param PriceGranularity Numeric. Step size of the price grid.
#' @param TrueOptimal Optional numeric. The true optimal reward value.
#' @return A list with:
#' \describe{
#'   \item{\code{SimData}}{The input \code{SimData} augmented with expected rewards, simulation IDs, and iterations.}
#'   \item{\code{OptimalReward}}{The optimal reward used for normalization.}
#' }
PrepareSimData <- function(SimData, ExRewards, Prices, NumIter, NumSimulations,
                                              PriceGranularity, TrueOptimal = NULL) {
  SimData$PricesTested = round(as.numeric(as.character(SimData$PricesTested)), -log10(PriceGranularity))
  Prices = round(as.numeric(as.character(Prices)), -log10(PriceGranularity))
  Index  = match(SimData$PricesTested, Prices)
  SimData$ExpectedReward = ExRewards[Index]
  SimData$SimID          = rep(1:NumSimulations, each = NumIter)
  SimData$Iteration      = rep(1:NumIter, times = NumSimulations)
  return(SimData)
}

#' @title Compute Mean Expected Reward as a Percentage of Optimal
#' @description
#' Computes the mean expected reward at each iteration of a bandit simulation,
#' expressed as a percentage of the true optimal reward. The function accounts
#' for multiple simulations and optionally computes cumulative performance over time.
#' @param SimData A data frame containing simulation results. Must include a \code{PricesTested}
#' column representing the price chosen at each iteration across simulations.
#' @param ExRewards A list containing the true optimal price, ratio of maximum possible to true optimal,
#' and the expected rewards at every price in test set.
#' @param Prices A vector of prices in the test set.
#' @param NumIter Integer. Number of iterations (consumers) in each simulation run.
#' @param NumSimulations Integer. Number of independent simulation runs.
#' @param PriceGranularity Numeric. Granularity (step size) of the price grid; used to match
#' observed prices to expected reward values via rounding.
#' @param TrueOptimal Optional numeric. If provided, used as the benchmark optimal reward.
#' If \code{NULL}, the optimal is computed as the maximum of \code{ExRewards}.
#' @param Cumulative Logical. If \code{TRUE}, computes the cumulative mean expected reward
#' as a percentage of cumulative optimal reward.
#' @return A numeric vector of length \code{NumIter}:
#' \itemize{
#'   \item If \code{Cumulative = FALSE}, each entry is the average expected reward at that iteration,
#'         normalized by the optimal reward.
#'   \item If \code{Cumulative = TRUE}, each entry reflects the cumulative average expected reward
#'         up to that iteration, normalized by the cumulative optimal reward.
#' }
#' @export
ComputeMeanExpectedReward <- function(SimData, ExRewards, Prices, NumIter, NumSimulations, 
                                      PriceGranularity, TrueOptimal = NULL, Cumulative = FALSE) {
  OptimalReward = if (!is.null(TrueOptimal)) TrueOptimal else max(ExRewards, na.rm = TRUE)
  SimData = PrepareSimData(SimData, ExRewards, Prices, NumIter, NumSimulations, PriceGranularity, OptimalReward)
  MeanExRewards = tapply(SimData$ExpectedReward, SimData$Iteration, mean, na.rm = TRUE)
  MeanExpectedRewardPercentage = (MeanExRewards / OptimalReward) * 100
  if (Cumulative) {
    CumExRewards = cumsum(MeanExRewards)
    CumOptimalRewards = cumsum(rep(OptimalReward, length(MeanExRewards)))
    MeanExpectedRewardPercentage = (CumExRewards / CumOptimalRewards) * 100
  }
  return(unname(MeanExpectedRewardPercentage))
}

#' @title Compute Expected Reward Matrix as a Percentage of Optimal
#' @description
#' Computes the expected reward matrix at each iteration of a bandit simulation,
#' expressed as a percentage of the true optimal reward. The function accounts
#' for multiple simulations and optionally computes cumulative performance over time.
#' @param SimData A data frame containing simulation results. Must include a \code{PricesTested}
#' column representing the price chosen at each iteration across simulations.
#' @param ExRewards A list containing the true optimal price, ratio of maximum possible to true optimal,
#' and the expected rewards at every price in test set.
#' @param Prices A vector of prices in the test set.
#' @param NumIter Integer. Number of iterations (consumers) in each simulation run.
#' @param NumSimulations Integer. Number of independent simulation runs.
#' @param PriceGranularity Numeric. Granularity (step size) of the price grid; used to match
#' observed prices to expected reward values via rounding.
#' @param TrueOptimal Optional numeric. If provided, used as the benchmark optimal reward.
#' If \code{NULL}, the optimal is computed as the maximum of \code{ExRewards}.
#' @param Cumulative Logical. If \code{TRUE}, computes the cumulative mean expected reward
#' as a percentage of cumulative optimal reward.
#' @return A numeric matrix of dimension \code{NumIter} × \code{NumSimulations}, where each column corresponds to a simulation run and each row to an iteration:
#' \itemize{
#'   \item If \code{Cumulative = FALSE}, each entry represents the expected reward at a given iteration in a given simulation, normalized by the true optimal reward.
#'   \item If \code{Cumulative = TRUE}, each entry represents the cumulative expected reward up to that iteration in a given simulation, normalized by the cumulative optimal reward.
#' }
#' 
#' @export
ComputeExpectedRewardMatrix <- function(SimData, ExRewards, Prices, NumIter, NumSimulations, 
                                        PriceGranularity, TrueOptimal = NULL, Cumulative = FALSE) {
  
  OptimalReward = if (!is.null(TrueOptimal)) TrueOptimal else max(ExRewards, na.rm = TRUE)
  SimData = PrepareSimData(SimData, ExRewards, Prices, NumIter, NumSimulations, PriceGranularity, OptimalReward)
  
  ExpectedRewardMatrix = matrix(NA, nrow = NumIter, ncol = NumSimulations)
  for (sim in 1:NumSimulations) {
    sim_data = SimData[SimData$Simulation == sim, ]
    ExpectedRewardMatrix[sim_data$Iteration, sim] = sim_data$ExpectedReward
  }
  
  ExpectedRewardPercentageMatrix = (ExpectedRewardMatrix / OptimalReward) * 100
  if (Cumulative) {
    CumExRewardsMatrix     = apply(ExpectedRewardMatrix, 2, cumsum)
    CumOptimalRewardMatrix = matrix(rep(cumsum(rep(OptimalReward, NumIter)), NumSimulations),
                                     nrow = NumIter, ncol = NumSimulations)
    ExpectedRewardPercentageMatrix = (CumExRewardsMatrix / CumOptimalRewardMatrix) * 100
  }
  return(ExpectedRewardPercentageMatrix)
}

### Histogram of Arms Played ### -----------------------------------------------

#' @title Compute Normalized Histogram of Arms Played
#' @description
#' This function computes a normalized histogram of the prices (arms) tested during an experiment,
#' comparing them against the full set of prices. It returns the relative frequency of each price being played.
#' @param PricesTested A vector of prices that were played during the experiment.
#' @param NumIter Integer. Number of iterations (consumers) in each simulation run.
#' @param Prices A vector of prices in the test set.
#' @return A numeric vector of relative frequencies, one for each price in the test set.
#' @export
Histogram <- function(PricesTested, NumIter, Prices){
  Tested = c(PricesTested, Prices)
  Table  = as.vector(table(Tested))
  Table  = Table - 1
  return(Table/(NumIter*10))
}

### Process Time Varying Results ### -------------------------------------------

#' @title Calculate Repeated Maximum Rewards Across Simulation Blocks
#' @description For each row in a simulation matrix, this function finds the maximum reward 
#' within each block of prices and repeats it a fixed number of times (SeasonLength) to 
#' construct a vector of repeated maximums.
#' @param SimData A data frame containing simulation results.
#' @param Prices A vector of prices in the test set.
#' @param SeasonLength Number of consumers per season; determines the frequency of seasonal demand shifts.
#' @return A numeric vector containing the repeated maximum values for each block in each simulation.
#' @export
CalculateMaxRepeats <- function(SimData, Prices, SeasonLength) {
  NumPrices = length(Prices)
  MaxValues = numeric()
  for (i in 1:nrow(SimData)) {
    for (j in seq(1, ncol(SimData), by = NumPrices)) {
      CurrentMax = max(SimData[i, j:min(j+NumPrices-1, ncol(SimData))])
      MaxValues  = c(MaxValues, rep(CurrentMax, SeasonLength))
    }
  }
  return(MaxValues)
}

#' @title Append Expected and Optimal Rewards to Simulation Data
#' @description
#' This function augments a simulation data frame by computing and appending columns 
#' for realized rewards, expected rewards, and optimal rewards for each observation, 
#' using information on prices, expected reward matrices, and seasonal structure.
#' @param SimData A data frame containing simulation results.
#' @param Prices A vector of prices in the test set.
#' @param NumIter Integer. Number of iterations (consumers) in each simulation run.
#' @param NumSimulations Integer. Number of independent simulation runs.
#' @param SeasonLength Number of consumers per season; determines the frequency of seasonal demand shifts.
#' @param ExRewardsTV A matrix containing the expected rewards at every price in test set for every simulation (time varying).
#' @param SeasonOpt A numeric vector containing the repeated optimal values for each block in each simulation.
#' @return A data frame identical to \code{SimData} but with additional columns for realized reward, 
#' action ID, iteration ID, season index, expected reward, and optimal reward.
#' @export
AppendExRewards <- function(SimData, Prices, NumIter, NumSimulations, SeasonLength, ExRewardsTV, SeasonOpt){
  SimData$Reward         = SimData$PricesTested * SimData$PurchaseDecisions
  SimData$ActionID       = as.integer(SimData$PricesTested*length(Prices))
  SimData$IterationID    = rep(seq(NumIter), NumSimulations)
  SimData$Season         = ((SimData$IterationID - 1) %/% SeasonLength) + 1
  SimData$ColumnID       = (SimData$Season*length(Prices) - length(Prices)) + SimData$ActionID
  SimData$ExpectedReward = apply(SimData, 1, function(x) ExRewardsTV[as.integer(x['SimulationID']), as.integer(x['ColumnID'])])
  SimData$OptReward      = SeasonOpt
  return(SimData)
}

#' @title Compute Cumulative Percentage of Optimal Rewards Over Time
#' @description
#' This function calculates the cumulative percentage of the optimal reward achieved 
#' across all simulations over time, based on expected and optimal rewards in the simulation data.
#' @param SimData A data frame containing simulation results.
#' @param NumIter Integer. Number of iterations (consumers) in each simulation run.
#' @param NumSimulations Integer. Number of independent simulation runs.
#' @return A numeric vector of length \code{NumIter} giving the cumulative percentage of the optimal reward 
#' obtained at each time step.
#' @export
CumulativeRewardTV <- function(SimData, NumIter, NumSimulations) {
  rsEx         = matrix(df$ExpectedReward, nrow = NumSimulations, ncol = NumIter, byrow = TRUE)
  ExColSums    = colSums(rsEx)
  rsOpt        = matrix(df$OptReward, nrow = NumSimulations, ncol = NumIter, byrow = TRUE)
  OptColSums   = colSums(rsOpt)
  PctOfOptimal = (cumsum(ExColSums)/cumsum(OptColSums))*100
  return(PctOfOptimal)
}
#' @title Subset Bandit Simulation Results for Plotting
#' @description 
#' This function extracts a subset of simulation results from a larger list based on the
#' specified distribution class, bandit policy, and identifying parameters such as the 
#' number of arms, beta shape, gap parameter (zeta), or seasonal shock (rho).
#' @param ResultsList A named list containing simulation results from various bandit experiments
#' @param DistrClass The type of distribution ('analytic', 'empirical', 'leftdigit', 'timevarying').
#' @param policies A vector of policy names (e.g., "TS", "GPTS") to filter results.
#' @param arm The number of arms as a character string (e.g., "10").
#' @param beta The beta distribution type as a character string (e.g., "29" for Beta(2,9)).
#' @param zeta Scaling factor for the discontinuity gap as a character string. (e.g., "20")
#' @param rho Maximum seasonal shock as a character string. (e.g., "01")
#' @return A named list containing only the simulation results matching the specified criteria.
#' @export
SubsetResults <- function(ResultsList, DistrClass, policies, arm = NULL, beta = NULL, zeta = NULL, rho = NULL) {
  Pattern <- switch(DistrClass,
                    "analytic" = {
                      if (is.null(beta)) stop("beta must be provided for analytic DistrClass")
                      paste0("Beta", beta, "_(", paste(policies, collapse = "|"), ")_", arm, "$")
                    },
                    "empirical" = {
                      paste0("Field_(", paste(policies, collapse = "|"), ")_", arm, "$")
                    },
                    "leftdigit" = {
                      if (is.null(beta) || is.null(zeta)) stop("Both beta and zeta must be provided for leftdigit DistrClass")
                      paste0("Beta", beta, "Zeta", zeta, "(", paste(policies, collapse = "|"), ")_", arm, "$")
                    },
                    "timevarying" = {
                      if (is.null(beta) || is.null(rho)) stop("Both beta and rho must be provided for timevarying DistrClass")
                      paste0("Beta", beta, "Rho", rho, "(", paste(policies, collapse = "|"), ")_10$")
                    },
                    stop("Unknown DistrClass: must be 'analytic', 'empirical', 'leftdigit', or 'timevarying'")
  )
  
  SelectedResults <- ResultsList[grep(Pattern, names(ResultsList))]
  return(SelectedResults)
}

#' @title Plot Cumulative Percentage of Optimal Rewards Over Time
#' @description
#' This function generates a line plot showing the cumulative percentage of optimal rewards
#' achieved over time for different bandit policies. It visualizes performance over the course
#' of consumer interactions in a simulation, with optional customization for appearance.
#' @param SelectedResults A list containing the relevant subset of the results for the plot.
#' @param PolicyNames A vector containing the policy names in the plot.
#' @param NumIter Integer. Number of iterations (consumers) in each simulation run.
#' @param Title A character string to be used as the title of the plot.
#' @param LineColors A character vector of colors for each policy line.
#' @param LineTypes A vector of line types corresponding to each policy.
#' @param RatioMaxTrue Ratio between max reward possible from price set and true optimal (value between 0 and 1)
#' @param legend_position Character string specifying the position of the legend ("none", "bottom").
#' @return A ggplot2 object representing the cumulative reward comparison plot.
#' @export
CumulativeRewardsPlot <- function(SelectedResults, PolicyNames, NumIter, Title,
                                LineColors, LineTypes, RatioMaxTrue, legend_position = "none"){
  
  Iteration            = seq(NumIter)
  Data                 = as.data.frame(do.call(cbind, SelectedResults))
  Data                 = cbind(seq(NumIter), Data)
  colnames(Data)       = c("Iteration", PolicyNames)
  DataMelted           = reshape2::melt(Data, id.var='Iteration')
  colnames(DataMelted) = c("Iteration", "Algorithm", "Value")
  DataMelted$Algorithm = factor(DataMelted$Algorithm, levels = PolicyNames)
  
  p <- ggplot(DataMelted, aes(x = Iteration, y = Value, linetype= Algorithm, colour=Algorithm)) +
    geom_line(linewidth=1.2) +
    geom_hline(yintercept = 100*(RatioMaxTrue), color = "black", linetype = "longdash", linewidth = 0.7) +
    scale_x_continuous(limits = c(0, NumIter)) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_minimal() +
    scale_color_manual(values = LineColors) +
    scale_linetype_manual(values = LineTypes) +
    labs(title = Title, x ="# of consumers tested", y = "% of optimal rewards") +
    theme(plot.title        = element_text(hjust = 0.5)) +
    theme(plot.title        = element_text(size = 24),
          axis.title.x      = element_text(size = 22),
          axis.title.y      = element_text(size = 22),
          axis.text.x       = element_text(size = 18),
          axis.text.y       = element_text(size = 18),
          legend.text       = element_text(size = 26),
          legend.title      = element_text(size = 26),
          legend.key.width  = unit(2.2, "cm"),
          legend.key.height = unit(0.5, "cm"))
  guide_args <- list(override.aes = list(linewidth = 2))
  if (legend_position != "none") {
    guide_args <- utils::modifyList(guide_args, list(nrow = 2, byrow = TRUE))
  }
  p <- p + guides(color = do.call(guide_legend, guide_args))
  return (p)
}

#' @title Plot Histogram of Arms Played with Relative Reward Benchmarks
#' @description
#' This function generates a bar plot showing the percentage of times each price (arm) was chosen 
#' by different bandit policies. A red benchmark line is overlaid to indicate how good each price is 
#' relative to the true optimal in terms of expected reward. The plot allows for visual comparison 
#' of policy behavior against the relative reward quality of each price.
#' @param SelectedResults A list containing the relevant subset of the results for the plot.
#' @param OptVec A numeric vector indicating the expected reward of each price relative to the true optimal.
#' @param Title A character string to be used as the title of the plot.
#' @param PolicyNames A vector containing the policy names in the plot.
#' @param Colors A vector of colors for the plot.
#' @param Prices A vector of prices in the test set.
#' @param PriceGranularity A numeric value indicating the spacing between price levels.
#' @param y_max A numeric value specifying the maximum value of the y-axis (default is 100).
#' @return A ggplot2 object representing the histogram of arms played, with overlaid benchmarks for relative expected rewards.
#' @export
ArmsPlayedPlot <- function(SelectedResults, OptVec, Title, PolicyNames, Colors, 
                             Prices, PriceGranularity, y_max = 100){
  
  df1           = cbind(Prices, data.frame(SelectedResults))
  colnames(df1) = c("Price", PolicyNames)
  DataMelted    = reshape2::melt(df1, id.var='Price') 
  df2           = data.frame(Price = Prices, Optimal = OptVec * y_max)
  max_optimal   = max(df2$Optimal)
  
  ggplot() +
    geom_bar(data = DataMelted, aes(x = Price, y = value, fill = variable), 
             stat = "identity", position = "dodge", width = 0.8 / length(Prices), 
             alpha = 1) + 
    geom_segment(data = df2, aes(x = Price - PriceGranularity / 2, 
                                 xend = Price + PriceGranularity / 2, 
                                 y = Optimal, yend = Optimal), 
                 color = "red", linewidth = 0.8) +
    geom_point(data = subset(df2, Optimal == max_optimal), aes(x = Price, y = Optimal),
               color = "black", shape = 16, size = 4) +
    labs(fill = "Algorithm") +
    xlab(paste0("price (granularity = ", PriceGranularity, ")")) +
    ylab("% of time arm is chosen") +
    ggtitle(Title) +
    theme_minimal() +
    scale_fill_manual(values = Colors) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = unique(pretty(Prices, n = min(length(Prices), 10)))) +
    scale_y_continuous(
      limits   = c(0, y_max),
      sec.axis = sec_axis(~ . * (100 / y_max), name = "% of optimal", breaks = seq(0, 100, by = 10))
    ) +
    theme(plot.title   = element_text(size = 24),
          axis.title.x = element_text(size = 22),
          axis.title.y = element_text(size = 22),
          axis.text.x  = element_text(size = 18),
          axis.text.y  = element_text(size = 18),
          legend.text  = element_text(size = 22),
          legend.title = element_text(size = 22))
}
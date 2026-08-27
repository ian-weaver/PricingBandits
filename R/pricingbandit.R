#' @title Run a Pricing Bandit Experiment
#'
#' @description
#' Wrapper for running a single multi-armed bandit pricing experiment with any of the
#' supported policies, selected by name. The willingness-to-pay (WTP) distribution is
#' entirely user-specified: pass one draw per consumer via `valuations` (e.g.
#' `rbeta(2500, 2, 9)`, draws from an empirical CDF, or any other process). At each
#' round the policy posts a price from `prices`; the consumer purchases if and only if
#' their valuation exceeds the price; the policy re-optimizes every `batch_size`
#' consumers.
#'
#' Available policies:
#' \describe{
#'   \item{"UCB"}{Upper Confidence Bound on independent arms.}
#'   \item{"TS"}{Thompson Sampling with Beta posteriors on independent arms.}
#'   \item{"GP-UCB"}{Gaussian Process UCB: correlated arms through a GP demand curve.}
#'   \item{"GP-TS"}{Gaussian Process Thompson Sampling.}
#'   \item{"GP-UCB-M"}{Monotonic GP-UCB: demand draws are monotone everywhere by construction (basis-function reconstruction from derivatives constrained to be non-positive).}
#'   \item{"GP-TS-M"}{Monotonic GP-TS (same construction).}
#' }
#'
#' @param valuations Numeric vector of consumer valuations (WTP), one per consumer.
#'   Its length determines the number of iterations.
#' @param prices Numeric vector of candidate prices (values between 0 and 1) - the arms.
#' @param policy Character; one of `"UCB"`, `"TS"`, `"GP-UCB"`, `"GP-TS"`,
#'   `"GP-UCB-M"`, `"GP-TS-M"`.
#' @param batch_size Number of consumers tested before the policy updates (default 10).
#' @param num_knots Number of knots for the monotonic ("-M") variants (default 11,
#'   the value used throughout the paper). The default is deliberately independent of
#'   the number of arms: with many more knots the truncated sampler operates in a much
#'   higher-dimensional, near-degenerate space and can break down, so 11 is recommended
#'   even for dense price grids (e.g. 100 arms).
#' @param timeout Seconds allowed for each truncated-sampling attempt in the
#'   monotonic fallback chain before the next fallback is tried (default 5).
#'   Increase on slow machines to give the exact sampler more time; decrease to
#'   fail over to the cheaper approximations sooner.
#' @param hetero Logical; if `TRUE`, use heterogeneous observation noise sampled each
#'   update via \code{\link{NoiseSample}} (default `FALSE`).
#' @param reset Optional; number of consumers after which the experiment history is
#'   wiped (for, e.g., time-varying demand). Default `NULL` (never reset).
#' @return A data frame with columns `PricesTested` and `PurchaseDecisions` (one row
#'   per consumer), with attributes:
#'   \describe{
#'     \item{policy}{The policy name used.}
#'     \item{diagnostics}{A list of fallback counters (see \code{\link{GetDiagnostics}}).}
#'   }
#' @examples
#' set.seed(1)
#' valuations <- rbeta(200, 2, 9)
#' out <- PricingBandit(valuations, prices = seq(10)/10, policy = "TS")
#' table(out$PricesTested)
#' @export
PricingBandit <- function(valuations, prices, policy = "GP-TS-M", batch_size = 10,
                          num_knots = 11, hetero = FALSE, reset = NULL, timeout = 5){

  PolicyMap = list("UCB"      = UCB,
                   "TS"       = TS,
                   "GP-UCB"   = GPUCB,
                   "GP-TS"    = GPTS,
                   "GP-UCB-M" = GPUCB_Mono,
                   "GP-TS-M"  = GPTS_Mono)

  policy = match.arg(policy, names(PolicyMap))

  if (!is.numeric(valuations) || length(valuations) < 1)
    stop("`valuations` must be a non-empty numeric vector (one WTP draw per consumer).")
  if (!is.numeric(prices) || length(prices) < 2 || any(prices <= 0) || any(prices > 1))
    stop("`prices` must be a numeric vector of at least two prices in (0, 1].")
  if (!is.numeric(batch_size) || batch_size < 1)
    stop("`batch_size` must be a positive integer.")
  if (!is.numeric(timeout) || timeout <= 0)
    stop("`timeout` must be a positive number of seconds.")

  Output = MABExperiment(Valuations  = valuations,
                         Policy      = PolicyMap[[policy]],
                         TestX       = prices,
                         NumIter     = length(valuations),
                         BatchSize   = batch_size,
                         NumKnots    = num_knots,
                         HeteroNoise = hetero,
                         Reset       = reset,
                         Timeout     = timeout)

  attr(Output, "policy")      = policy
  attr(Output, "diagnostics") = GetDiagnostics()
  return(Output)
}

#' @title Experiment Diagnostics
#' @description Internal infrastructure for tracking fallback/error-handling
#' events and monotonicity diagnostics during an experiment. Counters are
#' RNG-neutral: they never draw random numbers, so instrumented runs remain
#' bit-identical to uninstrumented ones.
#' @keywords internal
#' @noRd
.pb_diag <- new.env(parent = emptyenv())

#' @title Reset Experiment Diagnostics
#' @description Resets all internal fallback counters and the monotonicity log.
#' Called automatically at the start of \code{\link{MABExperiment}}; can be
#' called manually before using policy functions directly.
#' @return Invisibly, \code{NULL}.
#' @export
ResetDiagnostics <- function(){
  assign("HyperoptPrior",      0L, envir = .pb_diag)  # OptimalHyperparameters failed -> priors used
  assign("GPErrorPrevAS",     0L, envir = .pb_diag)  # GP policy errored -> previous action scores reused
  assign("MonoTimeout1",      0L, envir = .pb_diag)  # basis-mono: first attempt timed out
  assign("MonoRetryError",    0L, envir = .pb_diag)  # basis-mono: retry errored -> last resort
  assign("PolicyUpdates",     0L, envir = .pb_diag)  # number of policy evaluations performed
  invisible(NULL)
}

# Increment a counter (internal, RNG-neutral)
BumpDiag <- function(name){
  if (!exists(name, envir = .pb_diag, inherits = FALSE)) {
    assign(name, 0L, envir = .pb_diag)
  }
  assign(name, get(name, envir = .pb_diag) + 1L, envir = .pb_diag)
  invisible(NULL)
}

#' @title Get Experiment Diagnostics
#' @description Returns the internal fallback counters and monotonicity log
#' accumulated since the last \code{\link{ResetDiagnostics}} call. The same
#' information is attached as the \code{"diagnostics"} attribute of
#' \code{\link{PricingBandit}} output.
#' @return A named list of counters.
#' @export
GetDiagnostics <- function(){
  as.list(.pb_diag)
}

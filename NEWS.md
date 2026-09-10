# PricingBandits 2.1.0

* The four kernel functions are consolidated into one: `RBFKernel()` now takes
  derivative orders `d_i` and `d_j` (0 = function value, 1 = first derivative;
  default 0) and covers the value-value, value-derivative, derivative-value,
  and derivative-derivative cases itself. `RBFKernel_01()`, `RBFKernel_11()`,
  and `RBFKernel_All()` are removed. Results are unchanged: the new function
  performs the identical arithmetic, verified bit-for-bit against 2.0.0 at the
  kernel, posterior, and full seed-matched experiment level.
* DESCRIPTION reference formatted per CRAN guidance (year added, quotes
  removed from acronyms).

# PricingBandits 2.0.0

First public release of the streamlined package, extracted from the replication
code for Weaver, Kumar, and Jain, "Nonparametric Pricing Bandits Leveraging
Informational Externalities to Learn the Demand Curve" (Marketing Science).

* `PricingBandit()`: single entry point running one pricing experiment with any
  of six policies ("UCB", "TS", "GP-UCB", "GP-TS", "GP-UCB-M", "GP-TS-M"),
  a user-supplied willingness-to-pay vector, and an arbitrary price grid.
* Heterogeneous-noise extension (`hetero = TRUE`) available for all Gaussian
  process variants, and periodic history resets (`reset`) for time-varying
  demand.
* Lightweight run diagnostics (`GetDiagnostics()`) counting how often numerical
  fallback paths fired; the counters draw no random numbers, so results are
  unaffected.
* `num_knots` defaults to 11 for the monotonic variants regardless of the
  number of arms (larger values can destabilize the truncated sampler).
* The truncated-sampling `timeout` (default 5 seconds per attempt) is a
  user-settable argument of `PricingBandit()`.
* Results validated seed-for-seed bit-identical against the paper's
  replication code for all policies.

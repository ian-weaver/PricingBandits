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

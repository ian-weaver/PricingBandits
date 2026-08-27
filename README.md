# PricingBandits

Multi-armed bandit approaches to pricing experiments with an unknown demand
curve, from:

> Weaver, I. N., Kumar, V., & Jain, L. *Nonparametric Pricing Bandits Leveraging
> Informational Externalities to Learn the Demand Curve*. Marketing Science.

This package contains the streamlined, core bandit machinery. The full
replication code for the paper (WTP construction, expected-reward analysis,
figures) lives in the separate replication repository; results from this
package are validated seed-for-seed bit-identical against it.

> Looking for the **Python version**? See
> [pricingbandits-py](https://github.com/ian-weaver/pricingbandits-py).

## Installation

```r
# install.packages("devtools")
devtools::install_github("ian-weaver/PricingBandits")
```

## Quick start

You supply one consumer valuation (willingness to pay) per round — drawn from
*any* distribution you like — plus a set of candidate prices. The chosen policy
prices each arriving consumer and learns from buy/no-buy feedback alone.

```r
library(PricingBandits)

set.seed(1)
valuations <- rbeta(2500, 2, 9)          # any WTP process you like

out <- PricingBandit(valuations,
                     prices     = seq(10)/10,
                     policy     = "GP-TS-M",
                     batch_size = 10)

table(out$PricesTested)
attr(out, "diagnostics")
```

## Policies

| Name | Description |
|------|-------------|
| `"UCB"` | Upper Confidence Bound on independent arms |
| `"TS"` | Thompson Sampling with Beta posteriors per arm |
| `"GP-UCB"`, `"GP-TS"` | Gaussian-process variants — arms correlated through a GP demand curve (first informational externality) |
| `"GP-UCB-M"`, `"GP-TS-M"` | Monotonic variants — demand draws weakly decreasing everywhere by construction, via basis-function reconstruction from derivative-constrained draws (second informational externality) |

Options: `hetero = TRUE` enables the heterogeneous-noise extension for GP
variants; `reset = n` wipes the history every `n` consumers (time-varying
demand); `num_knots` controls the monotonic variants' constraint grid and
defaults to 11 **regardless of the number of arms** — with many more knots the
truncated sampler operates in a near-degenerate high-dimensional space and can
break down, so keep it modest even for dense price grids. `timeout` (default
5 seconds) bounds each truncated-sampling attempt before the fallback chain
advances; raise it on slow machines, lower it to fail over to cheaper
approximations sooner.

## Diagnostics

Every run counts the numerical fallback paths (hyperparameter-optimization
failures, truncated-sampler timeouts, last-resort samplers); in normal
operation these are all zero. See `GetDiagnostics()`.

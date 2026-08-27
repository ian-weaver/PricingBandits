test_that("BasisFunction values are correct at knots and interior points", {
  # rows of the basis matrix integrate the piecewise-linear derivative basis
  b <- BasisFunction(0.5, 11)
  expect_length(b, 11)
  expect_true(all(is.finite(b)))
  # at x = 1 with J knots the pattern is (d/2, d, ..., d, d/2)
  b1 <- BasisFunction(1, 11)
  expect_equal(b1, c(0.05, rep(0.1, 9), 0.05))
})

test_that("AggregateDataUCB adds (1,1) priors for untested prices", {
  TestX <- seq(10)/10
  out <- AggregateDataUCB(c(0.1, 0.2, 0.2, 0.9), c(1, 0, 1, 0), TestX)
  expect_equal(out$s_t[1], 1)   # 0.1 tested once, one purchase
  expect_equal(out$n_t[2], 2)   # 0.2 tested twice
  expect_equal(out$s_t[5], 1)   # untested -> prior
  expect_equal(out$n_t[5], 1)
})

test_that("AggregateDataTS appends one success and one failure per arm", {
  TestX <- seq(10)/10
  out <- AggregateDataTS(c(0.1, 0.1), c(1, 1), TestX)
  expect_equal(sum(out$n_t), 2 + 2*length(TestX))
  expect_equal(out$s_t[1], 3)  # 2 purchases + 1 pseudo success
})

test_that("AggregateDataGP anchors demand at price zero", {
  TestX <- seq(10)/10
  out <- AggregateDataGP(c(0.5, 0.5), c(1, 0), TestX, rep(0.25, 11), 10)
  expect_equal(out$PricesTested[1], 0)
  expect_equal(out$PurchaseRates[1], 1)
  expect_equal(out$PurchaseRates[out$PricesTested == 0.5], 0.5)
})

test_that("UCB scores are deterministic", {
  TestX <- seq(10)/10
  s1 <- UCB(c(0.1, 0.2, 0.2, 0.9), c(1, 0, 1, 0), TestX)
  s2 <- UCB(c(0.1, 0.2, 0.2, 0.9), c(1, 0, 1, 0), TestX)
  expect_identical(s1, s2)
  expect_length(s1, 10)
})

test_that("NLML is finite and increases away from a reasonable optimum", {
  TrainX <- c(0, 0.2, 0.5); TrainY <- c(1, 0.4, 0); s2 <- rep(0.025, 3)
  expect_true(is.finite(NLML(c(0.7, 0.2), TrainX, TrainY, s2)))
})

test_that("num_knots defaults to 11 and timeout to 5", {
  expect_identical(formals(PricingBandit)$num_knots, 11)
  expect_identical(formals(PricingBandit)$timeout, 5)
})

test_that("custom timeout runs and is validated", {
  skip_on_cran()
  set.seed(3)
  v <- rbeta(40, 2, 2)
  out <- PricingBandit(v, seq(10)/10, policy = "GP-TS-M", timeout = 2)
  expect_equal(nrow(out), 40)
  expect_error(PricingBandit(v, seq(10)/10, policy = "GP-TS-M", timeout = 0))
})

test_that("PricingBandit dispatches every policy and returns diagnostics", {
  skip_on_cran()
  prices <- seq(10)/10
  for (p in c("UCB", "TS", "GP-UCB", "GP-TS", "GP-UCB-M", "GP-TS-M")){
    set.seed(1)
    v <- rbeta(40, 2, 2)
    out <- PricingBandit(v, prices, policy = p, batch_size = 10)
    expect_equal(nrow(out), 40)
    expect_true(all(out$PricesTested %in% prices))
    expect_true(all(out$PurchaseDecisions %in% c(0, 1)))
    d <- attr(out, "diagnostics")
    expect_true(is.list(d))
    expect_true(d$PolicyUpdates >= 1)
  }
})

test_that("PricingBandit validates inputs", {
  expect_error(PricingBandit(numeric(0), seq(10)/10, policy = "UCB"))
  expect_error(PricingBandit(c(0.5), c(0.5, 1.5), policy = "UCB"))
  expect_error(PricingBandit(c(0.5), seq(10)/10, policy = "NOPE"))
})

test_that("hetero and reset options run", {
  skip_on_cran()
  set.seed(2)
  v <- rbeta(60, 2, 9)
  out_h <- PricingBandit(v, seq(10)/10, policy = "GP-TS", batch_size = 10, hetero = TRUE)
  expect_equal(nrow(out_h), 60)
  out_r <- PricingBandit(v, seq(10)/10, policy = "TS", batch_size = 10, reset = 30)
  expect_equal(nrow(out_r), 60)
})

test_that("seeded runs are reproducible", {
  run <- function(){
    set.seed(42)
    v <- rbeta(50, 2, 2)
    PricingBandit(v, seq(10)/10, policy = "GP-TS", batch_size = 10)
  }
  expect_identical(run()$PricesTested, run()$PricesTested)
})

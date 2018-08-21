context("Prior specification")

test_that("Normal prior constructed correctly", {
  set.seed(584034)
  n <- 100
  prior <- make_prior(rnorm, 0, 1, TRUE)
  prior_samples <- replicate(n, prior$rfn())
  expect_equal(prior$pars, c(0, 1))
  expect_equal(mean(prior_samples), 0, tol = 2 / sqrt(n))
  expect_equal(sd(prior_samples), 1, tol = 2 / sqrt(n))
  expect_true(prior$use)
  expect_equal(prior$label, "Normal(0, 1)")
})

test_that("Log-normal prior constructed correctly", {
  set.seed(584034)
  n <- 100
  prior <- make_prior(rlnorm, 0, 1, TRUE, "LogNormal(0, 1)")
  prior_samples <- replicate(n, prior$rfn())
  expect_equal(prior$pars, c(0, 1))
  expect_equal(mean(prior_samples), exp(0) + 1^2 / 2, tol = 2 / sqrt(n))
  expect_equal(median(prior_samples), exp(0), tol = 2 / sqrt(n))
  expect_true(prior$use)
  expect_equal(prior$label, "LogNormal(0, 1)")
})

test_that("Uniform prior constructed correctly", {
  set.seed(584034)
  n <- 100
  prior <- make_prior(runif, 0, 1, TRUE)
  prior_samples <- replicate(n, prior$rfn())
  expect_equal(prior$pars, c(0, 1))
  expect_equal(mean(prior_samples), 0.5, tol = 2 / sqrt(n))
  expect_equal(sd(prior_samples), sqrt(1/12), tol = 2 / sqrt(n))
  expect_true(all((prior_samples > 0) & (prior_samples < 1)))
  expect_true(prior$use)
  expect_equal(prior$label, "Uniform(0, 1)")
})

test_that("Log-uniform sampler correct", {
  set.seed(484948)
  n <- 100
  rlu_min <- 0.01
  rlu_max <- 0.2
  samp <- rlunif(n, rlu_min, rlu_max)
  analytic_mean <- (rlu_max - rlu_min) / (log(rlu_max) - log(rlu_min))
  expect_equal(mean(samp), analytic_mean, tol = 2 / sqrt(n))
  expect_true(all((samp > 0.01) & (samp < 0.2)))
})

test_that("Log-uniform prior constructed correctly", {
  set.seed(584034)
  n <- 100
  rlu_min <- 0.01
  rlu_max <- 0.2
  prior <- make_prior(rlunif, rlu_min, rlu_max, TRUE)
  prior_samples <- replicate(n, prior$rfn())
  expect_equal(prior$pars, c(rlu_min, rlu_max))
  expect_true(prior$use)
  expect_equal(prior$label, "Log-uniform(0.01, 0.2)")
})

test_that("Specify constant \"prior\" constructed correctly", {
  n <- 100
  z <- 2.39
  prior <- make_prior(z)
  prior_samples <- replicate(n, prior$rfn())
  expect_equal(prior_samples, rep(z, n))
  expect_equal(prior$pars, z)
  expect_true(prior$use)
  expect_equal(prior$label, "Constant(2.39)")
})

test_that("Unused priors constructed correctly", {
  prior <- make_prior(use = FALSE)
  expect_true(is.na(prior$rfn()))
  expect_false(prior$use)
})

context("Test analytic parameter estimation")

test_that("K loss function functioning", {
  K <- 600
  param_sample <- list(r_max = 0.08, z = 2.39)
  target_N <- 750
  tspan <- c(1990, 2000)
  data <- list(catch = c(20, 19, 46, 45, 43, 27, 40, 48, 44, 49))
  loss <- target_K(K = K,
                   param_sample = param_sample,
                   target_N = target_N,
                   data = data,
                   tspan = tspan)
  expect_equal(loss, -401.153651630282)
})

test_that("We can find K analytically", {
  param_sample <- list(r_max = 0.08, z = 2.39)
  target_N <- 350
  tspan <- c(1990, 2000)
  data <- list(catch = c(20, 19, 46, 45, 43, 27, 40, 48, 44, 49))
  data0 <- list(catch = rep(0, diff(tspan)))
  control <- sir_control()
  K_est <- find_K(param_sample, target_N, tspan, data, control)
  K_est0 <- find_K(param_sample, target_N, tspan, data0, control)
  expect_equal(K_est, 601.038771457972)
  expect_equal(K_est0, target_N)
})

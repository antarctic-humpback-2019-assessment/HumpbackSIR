context("Population dynamics")

tspan <- c(1990, 2000)
param <- list(r_max = 0.2,
              K = 1000,
              z = 2.39)
data0 <- list(catch = rep(0, diff(tspan)))
data_sust <- list(catch_sust = rep(1 - 0.6 * param$K, diff(tspan)))
num_years <- diff(tspan) + 1

test_that("K steady state", {
  P_Ksteady <- project_population(param_sample = param,
                                  data = data0,
                                  tspan = tspan)
  Pmin_Ksteady <- minimum_population(P_Ksteady)
  expect_equal(length(P_Ksteady$N), num_years)
  expect_equal(Pmin_Ksteady$N, rep(param$K, num_years))
  expect_equal(P_Ksteady$N[10], param$K)
  expect_false(check_mvp_violated(P_Ksteady, 0))
})

test_that("Zero steady state", {
  skip("Minimum population currently set to 1")
  P_0steady <- project_population(param_sample = param,
                                  data = data0,
                                  tspan = tspan)
  Pmin_0steady <- minimum_population(P_0steady)
  expect_equal(length(P_0steady$N), num_years)
  expect_equal(Pmin_0steady$N, 0)
  expect_equal(P_0steady$N[0], 0)
  expect_false(check_mvp_violated(P_0steady, 0))
})

test_that("Gen. logistic with catches", {
  tspan <- c(1990, 2000)
  data <- list(catch = c(10, 20, 40, 80, 160, 320, 40, 40, 20, 10))
  P_catches <- project_population(param_sample = param,
                                  data = data,
                                  tspan = tspan)
  Pmin_catches <- minimum_population(P_catches)
  ## Predicted values generated 2018-08-19; testing this will only tell us if
  ## values change
  pred <- c(1000, 990, 974.69935403084, 946.280493711608, 889.67821393339,
            773.049600240948, 524.089030824881, 566.52925365989,
            610.697507599572, 675.254944423984, 747.470442533935)
  expect_equal(P_catches$N, pred)
  expect_equal(Pmin_catches$year, 1996)
  expect_false(check_mvp_violated(P_catches, 400))
  expect_true(check_mvp_violated(P_catches, 600))
})


test_that("Gen. logistic with catches and violates MVP", {
  data <- list(catches = c(300, 300, 200, 200, 100, 100, 100, 0, 0, 0))
  P_catches_mvp <- project_population(param_sample = param,
                                      data = data,
                                      tspan = tspan)
  Pmin_catches_mvp <- minimum_population(P_catches_mvp)
  ## Predicted values generated 2018-08-19; testing this will only tell us if
  ## values change
  pred <- c(1000, 700, 480.308519406096, 359.721390994838, 225.417459001219,
            169.219630256064, 102.578850737432, 23.0057998076636,
            27.6064004514929, 33.126642886365, 39.750046460377)
  expect_equal(P_catches_mvp$N, pred)
  expect_equal(Pmin_catches_mvp$year, 1997)
  expect_true(check_mvp_violated(P_catches_mvp, 100))
})

context("Population dynamics")

num_yrs <- 10
start_yr <- 1
r_max <- 0.2
z <- 2.39
K <- 1000
N1 <- K
catch0 <- rep(0, num_yrs)
catch_sust <- rep(1 - 0.6 * K, num_yrs)

test_that("K steady state", {
  P_Ksteady <- GENERALIZED_LOGISTIC(r_max, K, N1, z,
                                    start_yr, num_yrs,
                                    catches = catch0,
                                    MVP = 0)
  expect_equal(length(P_Ksteady$Min_Yr), num_yrs)
  expect_equal(length(P_Ksteady$Pred_N), num_yrs)
  expect_equal(P_Ksteady$Min_Pop, K)
  expect_equal(P_Ksteady$Pred_N[10], K)
  expect_false(P_Ksteady$Violate_Min_Viable_Pop)
})

test_that("Zero steady state", {
  skip("Minimum population currently set to 1")
  P_0steady <- GENERALIZED_LOGISTIC(r_max, K, N1 = 0, z,
                                    start_yr, num_yrs,
                                    catches = catch0,
                                    MVP = 0)
  expect_equal(length(P_Ksteady$Min_Yr), num_yrs)
  expect_equal(length(P_Ksteady$Pred_N), num_yrs)
  expect_equal(P_Ksteady$Min_Pop, 0)
  expect_equal(P_Ksteady$Pred_N[10], 0)
  expect_false(P_Ksteady$Violate_Min_Viable_Pop)
})

test_that("Gen. logistic with MYS catch", {
  ## Need to start near steady state; values found with `optimize`
  num_yrs_sust <- 20
  P_sust <- GENERALIZED_LOGISTIC(r_max, K, N1 = 600.0109, z,
                                 start_yr, num_yrs_sust,
                                 catches = rep(84.6033, num_yrs_sust),
                                 MVP = 600)
  expect_equal(length(P_sust$Pred_N), num_yrs_sust)
  ## Increased tolerance; approaches steady state slowly
  expect_equal(P_sust$Pred_N[num_yrs_sust], 600, tol = 5)
  ## expect_equal(P_sust$Min_Yr, num_yrs_sust)
  ## expect_equal(P_sust$Min_Pop, P_sust$Pred_N[num_yrs_sust])
  expect_false(P_sust$Violate_Min_Viable_Pop)
})

test_that("Gen. logistic with catches", {
  catches <- c(10, 20, 40, 80, 160,
               320, 40, 40, 20, 10)
  P_catches <- GENERALIZED_LOGISTIC(r_max, K, N1, z,
                                    start_yr, num_yrs,
                                    catches = catches,
                                    MVP = 100)
  ## Predicted values generated 2018-08-19; testing this will only tell us if
  ## values change
  pred <- c(1000, 990, 974.69935403084, 946.280493711608, 889.67821393339,
            773.049600240948, 524.089030824881, 566.52925365989,
            610.697507599572, 675.254944423984)
  expect_equal(P_catches$Pred_N, pred)
  expect_equal(P_catches$Min_Yr, 7)
  expect_false(P_catches$Violate_Min_Viable_Pop)
})


test_that("Gen. logistic with catches and violates MVP", {
  catches <- c(300, 300, 200, 200, 100,
               100, 100, 0, 0, 0)
  P_catches_mvp <- GENERALIZED_LOGISTIC(r_max, K, N1, z,
                                        start_yr, num_yrs,
                                        catches = catches,
                                        MVP = 100)
  ## Predicted values generated 2018-08-19; testing this will only tell us if
  ## values change
  pred <- c(1000, 700, 480.308519406096, 359.721390994838, 225.417459001219,
            169.219630256064, 102.578850737432, 23.0057998076636,
            27.6064004514929, 33.126642886365)
  expect_equal(P_catches_mvp$Pred_N, pred)
  expect_equal(P_catches_mvp$Min_Yr, 8)
  expect_true(P_catches_mvp$Violate_Min_Viable_Pop)
})

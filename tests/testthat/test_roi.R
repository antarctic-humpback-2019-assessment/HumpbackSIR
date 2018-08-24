context("Rate of increase function")

test_that("ROI returns correctly", {
  param <- list(r_max = 0.05, z = 2.39, K = 30000)
  data = list(catch = Catch.data)
  names(data$catch) <- c("year", "catch")
  tspan <- c(min(data$catch$year), max(data$catch$year))
  Pred_N <- project_population(param, data, tspan)
  COMPUTING.ROI(data = Rel.Abundance,
                Pred_N = Pred_N,
                start_Yr = tspan[1])
})

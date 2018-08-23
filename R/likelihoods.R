
construct_ia_lik <- function(datname, dfn, data) {
  function(trajectory, param_sample, log = FALSE) {
    q <- param_sample[[datname]]
    obs <- data[[datname]]
    pred <- q * trajectory[trajectory$year == data[[datname]]$year]
    sum(mapply(dfn, obs, pred, data[[datname]]$sd))
  }
}

construct_abs_lik <- function(datname, dfn, data) {
  function(trajectory, param_sample, log = FALSE) {
    obs <- data[[datname]]
    pred <- trajectory[trajectory$year == data[[datname]]$year]
    sum(mapply(dfn, obs, pred, data[[datname]]$sd))
  }
}

## Data: list with year range (vector), growth rate, and sd. This one doesn't
## have to be a data.frame).
construct_growth_lik <- function(datname, dfn, data) {
  function(trajectory, param_sample, log = FALSE) {
    obs <- data[[datname]]$growth_rate
    pred_pop <- trajectory[trajectory$year == data[[datname]]$year]
    pred <- PRED.GROWTH.RATE(data[[datname]]$years,
                             pred_pop, start_Yr = 0L)
  }

}

dgamma_meansd <- function(x, mean, sd) {
  rate <- mean / sd ^ 2
  shape <- mean * rate
  dgamma(x, shape = shape, rate = rate)
}

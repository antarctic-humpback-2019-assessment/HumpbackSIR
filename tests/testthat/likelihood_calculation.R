context("Likelihood calculation")

test_that("Growth rate likelihood calculated correctly", {
    nll_gr <- LNLIKE.GR(0, 1, 2, log = T)
    expected_nll <- -dnorm( x = 0 , mean = 1 , sd = 2, log = TRUE )
    expect_equal(nll_gr, expected_nll)
})


LNLIKE.Ns <- function(Obs.N, Pred_N, start_yr, add_cv, log = TRUE) {
    N.yrs <- Obs.N$Year-start_yr+1
    nll_n <- -sum(
        dlnorm( # NOTE: can be changed to dlnorm_zerb
            x = Obs.N$N.obs,
            meanlog = log( Pred_N[N.yrs] ),
            sdlog = Obs.N$Sigma + add_cv,
            log))  ## Years for which Ns are available
    nll_n
}

test_that("Numbers likelihood calculated correctly", {
    set.seed(13212)
    n_obs <- 3
    predicted <- c( 3000, 4123, 5412)
    sigma <- 1
    N.obs = rlnorm(n_obs, meanlog = log(predicted), sdlog = sigma )
    observed <- data.frame(Year = c(1:n_obs), N.obs = N.obs, Sigma = rep(sigma, n_obs) )
    nll_N <- LNLIKE.Ns(observed, predicted, start_yr = 1, add_cv = 0, log = TRUE)
    expected_nll <- -sum(dlnorm( x = N.obs, meanlog = log(predicted), sdlog = sigma, log = TRUE ))
    expect_equal(nll_N, expected_nll)
})
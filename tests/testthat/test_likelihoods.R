context("Test likelihood constructors")

## Construct data for abundance likelihoods
data <- list(catch = Catch.data,
             index1 = Rel.Abundance[Rel.Abundance$Index == 1, ],
             index2 = Rel.Abundance[Rel.Abundance$Index == 2, ],
             abs_abund = Abs.Abundance.2005)
names(data$catch) <- c("year", "catch")
names(data$index1) <- c("index", "year", "obs", "cv")
data$index1$sd <- cv_to_sd(data$index1$cv)
names(data$index2) <- c("index", "year", "obs", "cv")
data$index2$sd <- cv_to_sd(data$index2$cv)
names(data$abs_abund) <- c("year", "obs", "cv")
data$abs_abund$sd <- cv_to_sd(data$abs_abund$cv)
data$missing <- data$abs_abund
data$missing$obs[1] <- NA


tspan <- c(min(data$catch$year), max(data$catch$year))
param <- list(r_max = 0.10,
              z = 2.39,
              K = 22000,
              q1 = 2,
              q2 = 1.2)
pred_N <- project_population(param, data, tspan)

test_that("Likelihood data checks are correct", {
  expect_true(is.null(check_lik_data(data, "index1")))
  expect_error(check_lik_data(data, "catch"), "year, obs, sd")
  expect_error(check_lik_data(data, "asdf"), "does not contain")
  expect_error(check_lik_data(data, "missing"), "Missing data not allowed")
})

ia1_likfun <- construct_ia_loglik(data_name = "index1",
                                  par_name = "q1",
                                  dfn = dlnorm,
                                  data = data,
                                  mean_link = log)
ia1_lik <- ia1_likfun(pred_N, param)
expect_equal(ia1_lik, -50.0604955045966)
ia2_likfun <- construct_ia_loglik(data_name = "index2",
                                  par_name = "q2",
                                  dfn = dlnorm,
                                  data = data,
                                  mean_link = log)
ia2_lik <- ia2_likfun(pred_N, param)

test_that("Index likelihoods construct", {
  expect_equal(ia2_lik, -2391.62382722937)
  expect_equal(ia1_lik + ia2_lik, -2441.68432273397)
})

abs_likfun_lnorm <- construct_abs_loglik(data_name = "abs_abund",
                                         dfn = dlnorm,
                                         data = data,
                                         mean_link = log)
abs_lnorm_loglik <- abs_likfun_lnorm(pred_N, param)
abs_likfun_norm <- construct_abs_loglik(data_name = "abs_abund",
                                        dfn = dnorm,
                                        data = data,
                                        mean_link = identity)
abs_norm_loglik <- abs_likfun_norm(pred_N, param)

test_that("Absolute likelihoods construct", {
  expect_equal(abs_lnorm_loglik, -1149.6769432598)
  expect_equal(abs_norm_loglik, -1424787306.52207)
})

## Add growth rate data
data$growth1 <- data.frame(start_year = 1995, end_year = 1998,
                           obs = 0.074, sd = 0.033)
data$growth2 <- data.frame(start_year = c(1995, 1998),
                           end_year = c(1998, 2001),
                           obs = c(0.074, 0.08),
                           sd = c(0.033, 0.03))
data$growth_names <- data$growth1
names(data$growth_names) <- paste0("a", 1:4)
data$growthNA <- data$growth2
data$growthNA$growth_rate[2] <- NA

test_that("Growth likelihood data checks check", {
  expect_true(is.null(check_growth_data(data, "growth1")))
  expect_true(is.null(check_growth_data(data, "growth2")))
  expect_error(check_growth_data(data, "growth_names"), "start_year, end_year")
  expect_error(check_growth_data(data, "asdf"), "does not contain")
  expect_error(check_growth_data(data, "growthNA"), "Missing data")
})

growth1_loglik <- construct_growth_loglik("growth1", dnorm, data, identity)
gr1 <- growth1_loglik(pred_N, param)
growth2_loglik <- construct_growth_loglik("growth2", dnorm, data, identity)
gr2 <- growth2_loglik(pred_N, param)

test_that("Growth rate likelihoods construct", {
  expect_equal(gr1, 2.28380425676316)
  expect_equal(gr2, 4.74120053121326)
})

liklist_ia1 <- list(ia1_likfun = ia1_likfun)
liklist_abs <- list(abs_likfun = abs_likfun_lnorm)
liklist_growth <- list(growth1_loglik = growth1_loglik)
liklist <- list(ia1_likfun = ia1_likfun,
                ia2_likfun = ia2_likfun,
                abs_likfun = abs_likfun_lnorm,
                growth_lik = growth1_loglik)

test_that("Lists of likelihoods return expected values", {
  expect_equal(calc_lik(pred_N, param, liklist_ia1, log = TRUE),
               ia1_lik)
  expect_equal(calc_lik(pred_N, param, liklist_abs, log = TRUE),
               abs_lnorm_loglik)
  expect_equal(calc_lik(pred_N, param, liklist_growth, log = TRUE),
               gr1)
  expect_equal(calc_lik(pred_N, param, liklist, log = TRUE),
               ia1_lik + ia2_lik + abs_lnorm_loglik + gr1)
  expect_equal(calc_lik(pred_N, param, liklist, log = FALSE),
               exp(ia1_lik + ia2_lik + abs_lnorm_loglik + gr1))
})


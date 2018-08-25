context("Test likelihood constructors")

data <- list(catch = Catch.data,
             index1 = Rel.Abundance[Rel.Abundance$Index == 1, ],
             index2 = Rel.Abundance[Rel.Abundance$Index == 2, ])
names(data$catch) <- c("year", "catch")
names(data$index1) <- c("index", "year", "obs", "cv")
data$index1$sd <- cv_to_sd(data$index1$cv)
names(data$index2) <- c("index", "year", "obs", "cv")
data$index2$sd <- cv_to_sd(data$index2$cv)

tspan <- c(min(data$catch$year), max(data$catch$year))
param <- list(r_max = 0.10,
              z = 2.39,
              K = 22000,
              q1 = 2,
              q2 = 1.2)
pred_N <- project_population(param, data, tspan)

test_that("Index likelihoods construct", {
  ia1_likfun <- construct_ia_loglik(data_name = "index1",
                                    par_name = "q1",
                                    dfn = dlnorm,
                                    data,
                                    mean_link = log)
  ia1_lik <- ia1_likfun(pred_N, param)
  expect_equal(ia1_lik, -50.0604955045966)
  ia2_likfun <- construct_ia_loglik("index2",
                                    "q2",
                                    dlnorm,
                                    data,
                                    mean_link = log)
  ia2_lik <- ia2_likfun(pred_N, param)
  expect_equal(ia2_lik, -2391.62382722937)
  old_lik <- LNLIKE.IAs(Rel.Abundance, Pred_N, tspan[1],
                        c(param$q1, param$q2), add.CV = 0,
                        log = TRUE)
  expect_equal(ia1_lik + ia2_lik, -2441.68432273397)
})

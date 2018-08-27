context("SIR")

test_that("Example runs", {
  set.seed(48448)

  ## Prepare data
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
  data$growth_rate <- data.frame(start_year = 1995, end_year = 1998,
                                 obs = 0.074, sd = 0.033)
  data$n_haplo <- 0

  ## Declare parameters
  priors <- make_prior_list(r_max = make_prior(runif, 0, 0.106),
                            K = make_prior(use = FALSE),
                            z = make_prior(2.39),
                            N_obs = make_prior(runif, 500, 200000, year = 2005),
                            add_CV = make_prior(use = FALSE),
                            q_ia1 = make_prior(use = FALSE),
                            q_ia2 = make_prior(use = FALSE))

  ## Declare likelihoods
  liklist <- list(construct_ia_loglik("index1", "q_ia1", dlnorm,
                                      data, mean_link = log),
                  construct_ia_loglik("index2", "q_ia2", dlnorm,
                                      data, mean_link = log),
                  construct_abs_loglik("abs_abund", dlnorm,
                                       data, mean_link = log),
                  construct_growth_loglik("growth_rate", dnorm,
                                          data, mean_link = "identity"),
                  construct_mvp_loglik("n_haplo", data))


  sir <- HUMPBACK.SIR(file.name = "test.N2005",
                      n.resamples = 100,
                      data = data,
                      priors = priors,
                      liklist = liklist,
                      output.Yrs = c(2005, 2006),
                      control = sir_control())
  ## Results generated 2018-08-19
  ## Only using 100 samples so tests don't take forever
  ## Only checking a few of these; not clear which should be targeted for tests
  ## yet. Maybe use a fake data set for testing? At least this shows that the
  ## function runs?
  r_max_summ <- c(0.0747308365546679, 0.0785385308959521, 0.0325653742793947,
                  0.103600497766363, 0.0427498911390547, 0.101975065643247,
                  0.0246829054611735, 0.104749660704751, 100)
  expect_equal(resample_summary$output_summary$r_max, r_max_summ)
  K_summ <- c(24341.0440866569, 24013.5352974968, 22809.4671260726,
              27420.3612174263, 22881.3564593546, 26393.1754034563,
              22759.0720310569, 28434.3853907956, 100)
  expect_equal(resample_summary$output_summary$K, K_summ)
  Nmin_summ <- c(523.562366839516, 387.891026320824, 176.885416772425,
                 1724.51699097789, 183.140884857596, 1266.76707418949,
                 163.174008013952, 2201.98115521626, 100)
  expect_equal(resample_summary$output_summary$Nmin, Nmin_summ)
  ## make sure that relevant files are created (and then removed)
  wd <- getwd()
  outfiles <- c(paste0(wd, "/test.N2005_resamples_trajectories.csv"),
                paste0(wd, "/test.N2005_resamples_output.csv"))
  expect_true(file.exists(outfiles[1]))
  file.remove(outfiles[1])
  expect_true(file.exists(outfiles[2]))
  file.remove(outfiles[2])
})

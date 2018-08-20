context("SIR")

test_that("Example runs", {
  set.seed(48448)
  sir <- HUMPBACK.SIR(file.name = "test.N2005",
                      n.samples = NULL,
                      n.resamples = 100,
                      prior.K = c(NA, NA, NA),
                      prior.r_max = c("uniform", 0, 0.106),
                      r_max.bound = c(0, 0.106),
                      prior.N.obs = c("uniform", 500, 20000),
                      prior.add.CV = c("uniform", 0, 1, FALSE),
                      prior.z = c(NA, 2.39, NA),
                      q.prior.IA = c("uniform", 0, 1, FALSE),
                      q.prior.Count = c("uniform", 0, 1, FALSE),
                      Klim = c(1, 500000),
                      target.Yr = 2005,
                      num.haplotypes = 0,
                      tolerance.for.bisection = 0.0001,
                      output.Yrs = c(2005, 2006),
                      abs.abundance = Abs.Abundance.2005,
                      rel.abundance = Rel.Abundance,
                      rel.abundance.key = TRUE,
                      count.data = Count.Data,
                      count.data.key = FALSE,
                      growth.rate.obs = c(0.074, 0.033, TRUE),
                      growth.rate.Yrs = c(1995, 1996, 1997, 1998),
                      catch.data = Catch.data,
                      Threshold = 1e-17,
                      Print = 0)
  ## Results generated 2018-08-19
  ## Only using 100 samples so tests don't take forever
  ## Only checking a few of these; not clear which should be targeted for tests
  ## yet. Maybe use a fake data set for testing? At least this shows that the
  ## function runs?
  r_max_summ <- c(0.0747308365546679, 0.0785385308959521, 0.0325653742793947,
                  0.103600497766363, 0.0427498911390547, 0.101975065643247,
                  0.0246829054611735, 0.104749660704751, 100)
  expect_equal(sir$resamples.output.summary$r_max, r_max_summ)
  K_summ <- c(24341.0440866569, 24013.5352974968, 22809.4671260726,
              27420.3612174263, 22881.3564593546, 26393.1754034563,
              22759.0720310569, 28434.3853907956, 100)
  expect_equal(sir$resamples.output.summary$K, K_summ)
  Nmin_summ <- c(523.562366839516, 387.891026320824, 176.885416772425,
                 1724.51699097789, 183.140884857596, 1266.76707418949,
                 163.174008013952, 2201.98115521626, 100)
  expect_equal(sir$resamples.output.summary$Nmin, Nmin_summ)
  ## make sure that relevant files are created (and then removed)
  wd <- getwd()
  outfiles <- c(paste0(wd, "/test.N2005_samples.output.csv"),
                paste0(wd, "/test.N2005_resample.trajectories.csv"),
                paste0(wd, "/test.N2005_resamples.output.csv"))
  expect_true(file.exists(outfiles[1]))
  file.remove(outfiles[1])
  expect_true(file.exists(outfiles[2]))
  file.remove(outfiles[2])
  expect_true(file.exists(outfiles[3]))
  file.remove(outfiles[3])
})

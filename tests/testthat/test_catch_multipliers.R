context("Catch multipliers")

test_that("Default catch multiplier is equal to 1", {
    ml <- make_multiplier_list()
    expect_equal(length(ml), 1)
    expect_equal(ml$c_mult_1$rfn(), 1)
    expect_equal(ml$c_mult_1$label, "Constant(1)")
})

test_that("Number of catch mutipliers and sampling periods mismatch", {
    expect_error(
        HUMPBACK.SIR(file_name = "test.N2005",
                     n_resamples = 100,
                     priors = make_prior_list(),
                     catch_multipliers = make_multiplier_list(make_prior(1),
                                                              make_prior(1)),
                     target.Yr = 2005,
                     num.haplotypes = 0,
                     output.Yrs = c(2005, 2006),
                     abs.abundance = Abs.Abundance.2005,
                     rel.abundance = Rel.Abundance,
                     rel.abundance.key = TRUE,
                     count.data = Count.Data,
                     count.data.key = FALSE,
                     growth.rate.obs = c(0.074, 0.033, TRUE),
                     growth.rate.Yrs = c(1995, 1996, 1997, 1998),
                     catch.data = Catch.data,
                     control = sir_control(threshold = 5e-19)))
})

test_that("Make catch multiplier can take multiple arguments", {
    catch_multipliers <- make_multiplier_list(make_prior(1),
                                              make_prior(2),
                                              make_prior(3),
                                              make_prior(4),
                                              make_prior(5),
                                              make_prior(6))

    expect_equal(length(catch_multipliers), 6)
    expect_setequal(names(catch_multipliers), paste0("c_mult_", 1:6))
    expect_setequal(sapply(catch_multipliers, function(x) x$rfn()), 1:6)
})

test_that("Multiple catch time series annd multipliers work", {
  skip("Move into integrated SIR test")
    Catch.data$Premodern <- rep(100, nrow(Catch.data))
    sir <- HUMPBACK.SIR(file_name = "test.N2005",
                        n_resamples = 100,
                        priors = make_prior_list(),
                        catch_multipliers = make_multiplier_list(make_prior(1),
                                                                 make_prior(2)),
                        target.Yr = 2005,
                        num.haplotypes = 0,
                        output.Yrs = c(2005, 2006),
                        abs.abundance = Abs.Abundance.2005,
                        rel.abundance = Rel.Abundance,
                        rel.abundance.key = TRUE,
                        count.data = Count.Data,
                        count.data.key = FALSE,
                        growth.rate.obs = c(0.074, 0.033, TRUE),
                        growth.rate.Yrs = c(1995, 1996, 1997, 1998),
                        catch.data = Catch.data,
                        control = sir_control(threshold = 1e-18,
                                              progress_bar = TRUE))
    expect_equal(sir$resamples_output$catch_multiplier_1, rep(1, 100))
    expect_equal(sir$resamples_output$catch_multiplier_2, rep(2, 100))
})

context("Test control function")

test_that("SIR control working", {
  expect_warning(sir_control(progress_bar = TRUE,
                             verbose = 1))
})

context("Test growth rate functions")

test_that("Growth rate returned correctly", {
  years  <- c(1995, 1998)
  pop <- c(1000, 2000)
  gr <- calc_growth_rate(years[1], years[2], pop[1], pop[2])
  expect_equal(gr, 0.231049060186648)
})

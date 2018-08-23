context("Test growth rate functions")

test_that("Growth rate returned correctly", {
  years  <- c(1995, 1998)
  pop <- c(1000, 2000)
  gr <- calc_growth_rate(years, pop)
  expect_equal(gr, 0.231049060186648)
})

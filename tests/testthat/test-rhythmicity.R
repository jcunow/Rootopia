# Sine fitting + rhythmicity. Need minpack.lm for the non-linear fit.

test_that("fit_sine_curve recovers a known sine", {
  skip_if_not_installed("minpack.lm")
  x <- seq(0, 48, by = 1)
  y <- 5 * sin(2 * pi * x / 24 + 0.5) + 10
  fit <- fit_sine_curve(x, y,
                        par_start = list(amp = 3, phase = 0, offset = 0, period = 24))
  expect_type(fit, "list")
  expect_equal(unname(abs(fit$amp)), 5, tolerance = 0.5)
  expect_equal(unname(fit$offset), 10, tolerance = 0.5)
})

test_that("rhythmicity returns its documented list", {
  skip_if_not_installed("minpack.lm")
  set.seed(1)
  x <- seq(0, 48, by = 1)
  y <- 5 * sin(2 * pi * x / 24) + 10 + rnorm(length(x), 0, 0.2)
  res <- rhythmicity(x, y, fix_period = 24)
  expect_type(res, "list")
})

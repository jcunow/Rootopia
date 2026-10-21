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

# A sine with known amplitude, period, phase and offset is the oracle: the fit
# has to give those numbers back. Nothing below is read off the function.

.rhy_signal <- function(amp = 1.5, period = 24, phase = 4, offset = 4,
                        sd = 0.1, n = 200, seed = 123) {
  set.seed(seed)
  x <- seq(0, 2 * period, length.out = n)
  list(x = x, y = amp * sin(2 * pi / period * (x + phase)) + offset +
                  stats::rnorm(n, 0, sd),
       amp = amp, period = period, phase = phase, offset = offset)
}

test_that("rhythmicity recovers the parameters of a known sine", {
  skip_if_not_installed("minpack.lm")
  s <- .rhy_signal()
  r <- rhythmicity(s$x, s$y, fix_period = s$period, method = "F")

  expect_equal(unname(r$amplitude), s$amp,    tolerance = 0.05)
  expect_equal(unname(r$phase),     s$phase,  tolerance = 0.05)
  expect_equal(unname(r$offset),    s$offset, tolerance = 0.05)
  expect_equal(unname(r$period),    s$period)

  # amp * sin(2*pi/P * (x + phase)) peaks where the argument is pi/2, i.e. at
  # x = P/4 - phase, so peakTime is fixed by the fitted parameters.
  expect_equal(unname(r$peakTime), (s$period / 4 - s$phase) %% s$period,
               tolerance = 0.05)
  expect_gt(r$R2, 0.98)
})

test_that("a free period is recovered too", {
  skip_if_not_installed("minpack.lm")
  s <- .rhy_signal()
  r <- rhythmicity(s$x, s$y, fix_period = NULL, method = "F",
                   par_start = list(amp = 1, phase = 0, offset = 0, period = 20))
  expect_equal(unname(r$period), s$period, tolerance = 0.1)
  expect_equal(unname(r$amplitude), s$amp, tolerance = 0.05)
})

test_that("the test statistics follow their definitions", {
  skip_if_not_installed("minpack.lm")
  s <- .rhy_signal()
  n <- length(s$y)

  # F = ((rss_null - rss_model)/(k-1)) / (rss_model/(n-k)), and R2 is
  # 1 - rss_model/rss_null, so F is fixed by R2 and the two df. Written this
  # way the check uses only returned values and the definition of an F test,
  # so it cannot be satisfied by copying the implementation.
  f <- rhythmicity(s$x, s$y, fix_period = 24, method = "F")
  k <- 3L                                    # amp, phase, offset
  expect_equal(unname(f$stat$F_stat),
               (f$R2 / (1 - f$R2)) * (n - k) / (k - 1))
  expect_equal(unname(f$df$df1), k - 1)
  expect_equal(unname(f$df$df2), n - k)

  # The likelihood ratio for Gaussian errors is n * log(rss_null/rss_model),
  # which is -n * log(1 - R2).
  l <- rhythmicity(s$x, s$y, fix_period = 24, method = "LR")
  expect_equal(unname(l$stat$LR_stat), -n * log(1 - l$R2))
  expect_equal(unname(l$p_value),
               stats::pchisq(unname(l$stat$LR_stat), k - 1, lower.tail = FALSE))
})

test_that("FLR is the F test for these errors, and LR is not", {
  skip_if_not_installed("minpack.lm")
  # A moderate signal, so the statistics are ordinary numbers. Comparing
  # p-values would prove nothing here: all.equal() falls back to an absolute
  # difference once the values are tiny, so two p-values of 1e-211 and 1e-208
  # compare equal.
  s <- .rhy_signal(sd = 1.5, n = 60)
  f   <- rhythmicity(s$x, s$y, fix_period = 24, method = "F")
  flr <- rhythmicity(s$x, s$y, fix_period = 24, method = "FLR")
  lr  <- rhythmicity(s$x, s$y, fix_period = 24, method = "LR")

  # (exp(LR/n) - 1) * (n-k)/r expands to ((rss_null - rss_model)/rss_model) *
  # (n-k)/r, which is the F statistic itself, on the same two df. The finite
  # sample correction therefore coincides exactly with the F test for these
  # errors; it is not an independent third answer.
  expect_equal(unname(flr$stat$FLR_stat), unname(f$stat$F_stat))
  expect_equal(unname(flr$df$df1), unname(f$df$df1))
  expect_equal(unname(flr$df$df2), unname(f$df$df2))

  # The large sample test is a different statistic against a different
  # reference distribution, and reports the df under different names.
  expect_gt(abs(unname(lr$stat$LR_stat) - unname(f$stat$F_stat)), 1)
  expect_named(lr$df, c("n", "dfdiff"))
  expect_named(f$df,  c("df1", "df2"))
})

test_that("noise is not reported as rhythmic", {
  skip_if_not_installed("minpack.lm")
  # The other half of having power: a test that only ever says "significant"
  # would pass every assertion above.
  set.seed(9)
  x <- seq(0, 48, length.out = 200)
  y <- stats::rnorm(200, 0, 1)
  for (m in c("F", "FLR", "LR")) {
    r <- rhythmicity(x, y, fix_period = 24, method = m)
    expect_lt(r$R2, 0.05)
    expect_gt(r$p_value, 0.1)
  }
})

test_that("rhythmicity and fit_sine_curve reject malformed input", {
  skip_if_not_installed("minpack.lm")
  expect_error(fit_sine_curve(1:3, 1:4), "equal length")
  expect_error(rhythmicity(1:10, 1:10, fix_period = "24"), "must be numeric")
  expect_error(rhythmicity(1:10, 1:10, method = "wald"), "must be one of")
})

# modal_peaks: peak detection on a distribution (display_type = "none" to avoid
# opening a graphics device during tests).

test_that("modal_peaks finds two peaks in a clear bimodal sample", {
  set.seed(2)
  x <- c(rnorm(300, -3, 0.5), rnorm(300, 3, 0.5))
  res <- modal_peaks(x, prominence_threshold = 0.001, display_type = "none")
  expect_type(res, "list")
  expect_true("peak_x" %in% names(res))
  expect_gte(length(res$peak_x), 2)
})

test_that("modal_peaks rejects an invalid display_type", {
  expect_error(modal_peaks(rnorm(50), display_type = "bogus"))
})

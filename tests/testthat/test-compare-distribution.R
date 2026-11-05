# compare_depth_distribution: distance between two depth profiles (numeric in,
# numeric out). Pure-numeric, so precise property checks are possible.

test_that("identical profiles give zero distance, different ones give > 0", {
  P <- c(1, 5, 10, 5, 1)
  expect_equal(compare_depth_distribution(P, P, metric = "wasserstein"),
               0, tolerance = 1e-6)
  Q <- c(1, 2, 3, 8, 12)
  expect_gt(compare_depth_distribution(P, Q, metric = "wasserstein"), 0)
})

test_that("js and kl metrics return a single numeric value", {
  P <- c(2, 4, 8, 4, 2); Q <- c(1, 3, 9, 5, 2)
  expect_true(is.numeric(compare_depth_distribution(P, Q, metric = "js")))
  expect_true(is.numeric(compare_depth_distribution(P, Q, metric = "kl")))
})

test_that("unequal lengths are padded with zeros, the same way for every metric", {
  P <- c(0.2, 0.3, 0.5); Q <- c(0.5, 0.5)
  for (m in c("js", "kl")) {
    expect_equal(compare_depth_distribution(P, Q, metric = m),
                 compare_depth_distribution(P, c(Q, 0), metric = m))
    expect_equal(compare_depth_distribution(P, Q, metric = m, cut = TRUE),
                 compare_depth_distribution(P[1:2], Q, metric = m))
  }
  # P has mass where the padded Q has none, so KL is infinite; JS is not.
  expect_equal(compare_depth_distribution(P, Q, metric = "kl"), Inf)
  expect_true(is.finite(compare_depth_distribution(P, Q, metric = "js")))
  expect_error(compare_depth_distribution(P, Q, alignPQ = FALSE), "differ in length")
  expect_error(compare_depth_distribution(P, Q, index = 1:2), "one position per bin")

  expect_equal(compare_depth_distribution(P, Q, metric = "wasserstein"),
               compare_depth_distribution(P, c(Q, 0), metric = "wasserstein"))
  # All mass moved one bin deeper costs exactly one bin.
  expect_equal(compare_depth_distribution(c(1, 0), c(0, 1)), 1)
})

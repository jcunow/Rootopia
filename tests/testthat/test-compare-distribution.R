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

# Pure-numeric distribution metrics: precise value + edge-case tests.

test_that("MRD is the root-weighted mean depth", {
  expect_equal(MRD(w = c(10, 20, 30), roots = c(0, 0, 1)), 30, tolerance = 1e-6)
  expect_equal(MRD(w = c(10, 20),     roots = c(1, 1)),     15, tolerance = 1e-6)
})

test_that("MRD validates inputs", {
  expect_error(MRD(c(1, 2), c(1, 2, 3)))      # length mismatch
  expect_error(MRD(c(1, 2), c(0, 0)))         # zero roots
  expect_error(MRD(c(1, 2), c(-1, 2)))        # negative roots
  expect_error(MRD(c(1, NA), c(1, 2)))        # NA
})

test_that("root_accumulation returns a cumulative vector aligned to input rows", {
  df <- data.frame(
    depth  = c(seq(0, 80, 20), seq(0, 80, 20)),
    Plot   = c(rep("a", 5), rep("b", 5)),
    rootpx = c(5, 50, 20, 15, 5, 10, 40, 30, 10, 5)
  )
  acc <- root_accumulation(df, group = "Plot", depth = "depth", variable = "rootpx")
  expect_length(acc, nrow(df))
  expect_equal(acc[5], 95)                     # group a cumulative total
  expect_true(all(diff(acc[1:5]) >= 0))        # monotonic within a group
})

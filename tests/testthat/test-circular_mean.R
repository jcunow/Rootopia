# Pure-numeric leaf function: demonstrates parameter + edge-case testing with
# no image dependencies. Good template for the other math helpers
# (MRD, RPI, root_accumulation, count_pixels, ...).

test_that("known circular means are correct", {
  expect_equal(circular_mean(c(0, 90), input_units = "degrees",
                             output_units = "degrees"), 45)
  # wrap-around: mean of 350 and 10 degrees is the 0/360 direction -> 0
  expect_equal(circular_mean(c(350, 10)), 0, tolerance = 1e-6)
})

test_that("unit conversions round-trip", {
  deg <- circular_mean(c(0, 90), output_units = "degrees")
  rad <- circular_mean(c(0, 90), output_units = "radians")
  expect_equal(rad, deg * pi / 180, tolerance = 1e-8)
})

test_that("invalid arguments error", {
  expect_error(circular_mean("a"))
  expect_error(circular_mean(numeric(0)))
  expect_error(circular_mean(c(1, 2), input_units = "grads"))
})

test_that("NA / Inf are dropped with a warning", {
  expect_warning(out <- circular_mean(c(0, 90, NA, Inf)))
  expect_equal(out, 45)
})

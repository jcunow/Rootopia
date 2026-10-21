# Regression test for the tryCatch(warning=) fix (#1): an internal warning()
# must NOT abort the function and return NULL. Before the fix, an Inf in the
# depthmap raised a warning that the top-level handler turned into an early
# NULL return; now it warns and still returns the binned raster.

test_that("binning warns on Inf but still returns a SpatRaster", {
  skip_if_not_installed("terra")
  dm <- terra::rast(matrix(c(1, 2, Inf, 4, 5, 6), nrow = 2))
  expect_warning(res <- binning(dm, nn = 1))
  expect_s4_class(res, "SpatRaster")          # would be NULL before the #1 fix
})

test_that("binning produces multiples of the bin width", {
  skip_if_not_installed("terra")
  dm  <- terra::rast(matrix(seq(0, 19), nrow = 4))
  res <- binning(dm, nn = 5, round_option = "floor")
  vals <- terra::values(res)
  expect_true(all(vals %% 5 == 0 | is.na(vals)))
})

test_that("binning validates its arguments", {
  skip_if_not_installed("terra")
  dm <- terra::rast(matrix(1:4, 2))
  expect_error(binning(dm))                      # missing nn
  expect_error(binning(dm, nn = -1))             # non-positive
  expect_error(binning(dm, nn = 5, round_option = "nope"))
})


test_that("the three rounding options place depths in the bins they name", {
  skip_if_not_installed("terra")
  # nn * f(depth / nn) for f = round, floor, ceiling, worked out by hand for
  # bin width 5. No depth sits on a .5 boundary, so R's banker's rounding at
  # exactly .5 never decides a case here.
  dm <- terra::rast(matrix(c(0, 2, 3, 7, 8, 12), nrow = 1))

  expect_equal(as.vector(terra::values(binning(dm, nn = 5, round_option = "rounding"))),
               c(0, 0, 5, 5, 10, 10))
  expect_equal(as.vector(terra::values(binning(dm, nn = 5, round_option = "floor"))),
               c(0, 0, 0, 5, 5, 10))
  expect_equal(as.vector(terra::values(binning(dm, nn = 5, round_option = "ceiling"))),
               c(0, 5, 5, 10, 10, 15))
})

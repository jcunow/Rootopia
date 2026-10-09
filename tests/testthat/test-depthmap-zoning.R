# Depth map construction + zoning / slicing utilities.

test_that("create_depthmap returns a raster matching the input grid", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  dm  <- create_depthmap(img, sinoid = TRUE, tube_thicc = 7, tilt = 45, dpi = 300)
  expect_s4_class(dm, "SpatRaster")
  # the depth map carries the input grid: same dims, same extent/res/CRS
  expect_equal(dim(dm)[1:2], dim(img)[1:2])
  expect_true(terra::compareGeom(img, dm, stopOnError = FALSE))
})

test_that("create_depthmap validates parameters", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  expect_error(create_depthmap(img, tilt = 0))
  expect_error(create_depthmap(img, tube_thicc = -1))
  expect_error(create_depthmap(img, center_offset = 2))   # must be in [0,1]
})

test_that("center_offset sweeps the tube top over the full image height", {
  skip_if_not_installed("terra")
  # One period of the curvature wave spans the tube circumference,
  # round(tube_thicc * pi * dpi / 2.54) - 1 = 129 px here, so an image 129 rows
  # tall holds exactly one period and the trough row reads the phase directly.
  # The trough is the shallowest point, i.e. the top of the tube, and it sits on
  # row 1 + (nr - 1) * center_offset.
  img <- make_binary_spatraster(nr = 129)
  tube_top_row <- function(offset) {
    dm <- create_depthmap(img, sinoid = TRUE, tube_thicc = 3.5, tilt = 45,
                          dpi = 30, center_offset = offset)
    which.min(terra::as.matrix(dm, wide = TRUE)[, 1])
  }
  expect_equal(tube_top_row(0),    1L)
  expect_equal(tube_top_row(0.25), 33L)
  expect_equal(tube_top_row(0.50), 65L)
  expect_equal(tube_top_row(0.75), 97L)
  # Monotone over the whole range, not just the first half.
  offsets <- seq(0, 0.99, by = 0.01)
  expect_true(!is.unsorted(vapply(offsets, tube_top_row, integer(1))))
  expect_equal(tube_top_row(0.99), 128L)
})

test_that("center_offset 0 and 1 are the same point on the tube", {
  skip_if_not_installed("terra")
  # A full turn returns to the image seam, so the two ends agree and the trough
  # sits on the first and last row at once.
  img <- make_binary_spatraster(nr = 129)
  dm <- function(offset) {
    terra::as.matrix(create_depthmap(img, sinoid = TRUE, tube_thicc = 3.5,
                                     tilt = 45, dpi = 30,
                                     center_offset = offset), wide = TRUE)
  }
  expect_equal(dm(0), dm(1))
  col1 <- dm(0)[, 1]
  expect_equal(col1[1], col1[129])
  expect_equal(which.min(col1), 1L)
})

test_that("slice_rotation splits into n contiguous slices", {
  skip_if_not_installed("terra")
  img    <- make_binary_spatraster()
  slices <- slice_rotation(img, n = 3)
  expect_length(slices, 3)
  expect_true(all(vapply(slices, inherits, logical(1), "SpatRaster")))
  total_rows <- sum(vapply(slices, terra::nrow, numeric(1)))  # terra::nrow is double
  expect_equal(total_rows, terra::nrow(img))              # exact, no overlap/gap
})

test_that("slice_rotation validates n", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  expect_error(slice_rotation(img, n = 0))
  expect_error(slice_rotation(img, n = terra::nrow(img) + 1))
})

test_that("depth_zoning masks pixels outside the selected depth bin", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster(24, 28)
  dmat <- matrix(10, 24, 28); dmat[1:12, ] <- 5      # top half = 5, bottom = 10
  dm <- terra::rast(dmat)
  z <- depth_zoning(img, depth_map = dm, depth = 5)
  expect_s4_class(z, "SpatRaster")
  expect_equal(terra::nrow(z), terra::nrow(img))      # extent unchanged
  zmat <- matrix(terra::values(z), nrow = 24, byrow = TRUE)
  expect_true(all(is.na(zmat[13:24, ])))              # depth-10 half is masked out
  expect_false(all(is.na(zmat[1:12, ])))              # depth-5 half retained
})

test_that("depth_zoning requires depth_map and depth", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  expect_error(depth_zoning(img, depth = 5))
  expect_error(depth_zoning(img, depth_map = make_binary_spatraster()))
})

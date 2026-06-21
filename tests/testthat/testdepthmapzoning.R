# Depth map construction + zoning / slicing utilities.

test_that("create_depthmap returns a raster matching the input grid", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  dm  <- create_depthmap(img, sinoid = TRUE, tube.thicc = 7, tilt = 45, dpi = 300)
  expect_s4_class(dm, "SpatRaster")
})

test_that("create_depthmap validates parameters", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  expect_error(create_depthmap(img, tilt = 0))
  expect_error(create_depthmap(img, tube.thicc = -1))
  expect_error(create_depthmap(img, center.offset = 2))   # must be in [0,1]
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

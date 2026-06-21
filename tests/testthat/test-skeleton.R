# Skeletonisation + skeleton point detection.

test_that("skeletonize_image returns a binary single-layer raster", {
  skip_if_not_installed("terra")
  sk <- skeletonize_image(make_binary_spatraster(), verbose = FALSE)
  expect_s4_class(sk, "SpatRaster")
  vals <- terra::values(sk)
  expect_true(all(vals %in% c(0, 1) | is.na(vals)))
  expect_gt(sum(vals == 1, na.rm = TRUE), 0)     # something survived
})

test_that("skeletonize_image errors on an empty image", {
  skip_if_not_installed("terra")
  expect_error(skeletonize_image(terra::rast(matrix(0, 24, 28)), verbose = FALSE))
})

test_that("detect_skeleton_points returns a list", {
  skip_if_not_installed("terra")
  res <- detect_skeleton_points(make_binary_spatraster(), skeletonize = TRUE)
  expect_type(res, "list")
})

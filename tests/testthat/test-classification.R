# RGB material classification. Uses the soil_* names (post-rename); adjust if
# your package still uses the older names.

test_that("classify_soil_rgb returns map + metrics on an RGB raster", {
  skip_if_not_installed("terra")
  rgb <- terra::rast(make_rgb_array() * 255)
  res <- classify_soil_rgb(rgb, verbose = FALSE)
  expect_type(res, "list")
  expect_true(all(c("map", "metrics", "centroids") %in% names(res)))
  expect_s4_class(res$map, "SpatRaster")
})

test_that("build_soil_centroids returns a centroid table from picks", {
  picks <- list(
    dark = matrix(c(28, 22, 18, 32, 26, 21), ncol = 3, byrow = TRUE),
    root = matrix(c(180, 160, 130, 175, 155, 125), ncol = 3, byrow = TRUE)
  )
  max_dist <- c(dark = 14, root = 26)
  cents <- build_soil_centroids(picks, max_dist, verbose = FALSE)
  expect_s3_class(cents, "data.frame")
  expect_true(all(c("class", "L", "A", "B", "MAX_DIST") %in% names(cents)))
  expect_equal(nrow(cents), 2L)
})

test_that("build_soil_centroids errors when max_dist is missing a class", {
  picks <- list(dark = matrix(c(28, 22, 18), ncol = 3))
  expect_error(build_soil_centroids(picks, max_dist = c(other = 14), verbose = FALSE))
})

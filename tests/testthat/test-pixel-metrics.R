# Pixel / color metrics. These delegate input handling to load_flexible_image,
# so one representative input each, plus the property each is supposed to compute.

test_that("count_pixels counts foreground pixels", {
  skip_if_not_installed("terra")
  m <- make_binary_matrix()
  expect_equal(count_pixels(terra::rast(m)), sum(m))
})

test_that("rgb2gray returns a single-layer raster", {
  skip_if_not_installed("terra")
  rgb <- terra::rast(make_rgb_array())
  g <- rgb2gray(rgb)
  expect_s4_class(g, "SpatRaster")
  expect_equal(terra::nlyr(g), 1L)
})

test_that("rgb2gray requires exactly 3 layers", {
  skip_if_not_installed("terra")
  expect_error(rgb2gray(make_binary_spatraster()))
})

test_that("tube_coloration returns the documented color data.frame", {
  skip_if_not_installed("terra")
  rgb <- terra::rast(make_rgb_array() * 255)
  cv <- tube_coloration(rgb)
  expect_s3_class(cv, "data.frame")
  expect_true(all(c("rcc", "gcc", "bcc", "hue", "saturation", "luminosity")
                  %in% names(cv)))
})

# Shape only: that the output is a 0/1 raster is true even of a threshold that
# marks every pixel, so whether the two methods actually separate anything is
# tested against a known mask in test-image-threshold.R.
test_that("image_threshold returns a binarised raster across methods", {
  skip_if_not_installed("terra")
  set.seed(1)
  r <- terra::rast(matrix(runif(24 * 28), 24, 28))
  outs <- list(
    global   = image_threshold(r, threshold = 0.4, method = "global",
                               select_layer = 1, mask_layer = NULL,
                               binary_01 = TRUE),
    adaptive = image_threshold(r, method = "adaptive", adaptive_k = 1.15,
                               window_size = 21, select_layer = 1,
                               mask_layer = NULL, binary_01 = TRUE))
  for (out in outs) {
    expect_s4_class(out, "SpatRaster")
    vals <- terra::values(out)
    expect_true(all(vals %in% c(0, 1) | is.na(vals)))
  }
})

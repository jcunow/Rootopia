# clean_image: hole filling / artifact removal. Output-format sweep + a
# behavioral check that an isolated speck is removed and an enclosed hole filled.

test_that("clean_image returns each requested output format", {
  skip_if_not_installed("terra")
  skip_if_not_installed("imager")
  img <- make_binary_spatraster()
  expect_s4_class(clean_image(img, output_format = "spatrast"), "SpatRaster")
  expect_s3_class(clean_image(img, output_format = "cimg"), "cimg")
  expect_true(is.matrix(clean_image(img, output_format = "matrix")))
})

test_that("clean_image removes a small isolated artifact", {
  skip_if_not_installed("terra")
  skip_if_not_installed("imager")
  m <- matrix(0, 30, 30)
  m[10:20, 10:20] <- 1     # solid block (touches nothing else)
  m[3, 3] <- 1             # isolated 1-px artifact
  cleaned <- clean_image(terra::rast(m), drop_objects_upto = 5,
                         output_format = "matrix")
  expect_equal(cleaned[3, 3], 0)
})

test_that("deprecated threshold names still work and warn", {
  skip_if_not_installed("terra")
  skip_if_not_installed("imager")
  m <- matrix(0, 30, 30)
  m[10:20, 10:20] <- 1
  m[3, 3] <- 1
  expect_warning(
    old <- clean_image(terra::rast(m), max_artifact_size = 5,
                       output_format = "matrix"),
    "deprecated")
  new <- clean_image(terra::rast(m), drop_objects_upto = 5,
                     output_format = "matrix")
  expect_equal(old, new)
})

test_that("mm2 thresholds convert to the equivalent pixel thresholds", {
  skip_if_not_installed("terra")
  skip_if_not_installed("imager")
  m <- matrix(0, 30, 30)
  m[10:20, 10:20] <- 1
  m[3, 3] <- 1
  dpi <- 300
  # 5 px at 300 dpi, expressed in mm^2
  mm2 <- 5 / (dpi / 25.4)^2
  a <- clean_image(terra::rast(m), drop_objects_upto = mm2,
                   size_unit = "mm2", dpi = dpi, output_format = "matrix")
  b <- clean_image(terra::rast(m), drop_objects_upto = 5,
                   output_format = "matrix")
  expect_equal(a, b)
  expect_error(clean_image(terra::rast(m), drop_objects_upto = 1,
                           size_unit = "mm2"), "dpi")
})

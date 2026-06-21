# Smoke tests for the remaining exported functions + buffer/zone helpers.

test_that("deep_drive returns a single proportion (or NA)", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  dm  <- terra::t(create_depthmap(img, center_offset = 0, tube_thicc = 3.5))
  res <- suppressWarnings(deep_drive(DepthMap = dm, RootMap = img,
                                     select_layer_rm = 1))
  expect_true(length(res) == 1L && (is.na(res) || (res >= 0 && res <= 1)))
})

test_that("deep_drive 'all' returns spatial outputs", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  dm  <- terra::t(create_depthmap(img, center_offset = 0, tube_thicc = 3.5))
  res <- suppressWarnings(deep_drive(DepthMap = dm, RootMap = img,
                                     select_layer_rm = 1, return = "all"))
  expect_type(res, "list")
  expect_true("deep_drive" %in% names(res))
})

test_that("create_root_buffer returns a raster", {
  skip_if_not_installed("terra")
  skip_if_not_installed("imager")
  img <- make_binary_spatraster()
  buf <- create_root_buffer(img, width = 2, halo_only = TRUE)
  expect_s4_class(buf, "SpatRaster")
})

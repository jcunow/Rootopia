# The input "seam": test the full input-type x output-format matrix HERE, once.
# Downstream functions that just delegate to load_flexible_image() then only
# need one representative input each.

test_that("every input type converts to a SpatRaster", {
  skip_if_not_installed("terra")
  expect_runs_on_all_inputs(
    load_flexible_image,
    output_format = "spatrast", scale = "none",
    check = function(res, nm) expect_s4_class(res, "SpatRaster")
  )
})

test_that("output_format is case-insensitive and accepts aliases", {
  skip_if_not_installed("terra")
  arr <- make_rgb_array()
  for (fmt in c("SpatRaster", "spatrast", "SPATRASTER")) {
    expect_s4_class(load_flexible_image(arr, output_format = fmt, scale = "none"),
                    "SpatRaster")
  }
})

test_that("scale transforms land in the expected range", {
  skip_if_not_installed("terra")
  arr255 <- make_rgb_array() * 255

  r01 <- load_flexible_image(arr255, output_format = "spatrast", scale = "to_01")
  expect_lte(max(terra::global(r01, "max", na.rm = TRUE)[[1]]), 1)  # per-layer maxima

  rbin <- load_flexible_image(arr255, output_format = "spatrast", scale = "binary")
  vals <- terra::values(rbin)
  expect_true(all(vals %in% c(0, 1) | is.na(vals)))
})

test_that("invalid inputs error cleanly (not silently)", {
  expect_error(load_flexible_image("does_not_exist.tif"))
  expect_error(load_flexible_image(make_rgb_array(), output_format = "not_a_format",
                                   scale = "none"))
})

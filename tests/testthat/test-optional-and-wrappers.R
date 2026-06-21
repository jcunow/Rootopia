# Functions that need optional packages or real on-disk fixtures.
# These are guarded so the suite stays green in a minimal environment; fill in
# the fixtures to turn them into full tests.

# ---- optional-dependency functions (skip if the suggested pkg is absent) ----

test_that("root_scape_metrics runs when landscapemetrics is available", {
  skip_if_not_installed("terra")
  skip_if_not_installed("landscapemetrics")
  res <- root_scape_metrics(make_binary_spatraster(), select.layer = 1)
  expect_s3_class(res, "data.frame")
})

test_that("analyze_soil_texture runs on an RGB raster", {
  skip_if_not_installed("terra")
  rgb <- terra::rast(make_rgb_array() * 255)
  expect_error(analyze_soil_texture(rgb), regexp = NA)   # expect no error
})

test_that("estimate_rotation_center runs when RStoolbox is available", {
  skip_if_not_installed("terra")
  skip_if_not_installed("RStoolbox")
  rgb <- terra::rast(make_rgb_array() * 255)
  expect_silent_or_warning <- suppressWarnings(estimate_rotation_center(rgb))
  expect_true(is.numeric(expect_silent_or_warning) || is.na(expect_silent_or_warning))
})

test_that("estimate_rotation_shift runs when imagefx is available", {
  skip_if_not_installed("terra")
  skip_if_not_installed("imagefx")
  rgb <- terra::rast(make_rgb_array() * 255)
  res <- estimate_rotation_shift(rgb, rgb, select.layer = 1)
  expect_true(all(c("depth", "rotation", "peak") %in% names(res)))
})

# ---- wrappers that need real image directories / multi-session fixtures ----
# Provide a folder of segmented (and optionally RGB + skeleton) scans named per
# `group_regex`, then replace skip() with a real end-to-end assertion.

test_that("root_depth_metrics end-to-end (needs an image directory fixture)", {
  skip("Provide a directory of segmented scans + skeletons to enable this test.")
})

test_that("stitch_root_scans end-to-end (needs a scan-sequence fixture)", {
  skip("Provide a directory of overlapping scans to enable this test.")
})

test_that("root_turnover end-to-end (needs RootDetector-style images)", {
  skip("Provide two registered RootDetector difference images to enable this test.")
})

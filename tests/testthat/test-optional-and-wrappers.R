# Functions that need optional packages or real on-disk fixtures.
# These are guarded so the suite stays green in a minimal environment; fill in
# the fixtures to turn them into full tests.

# ---- optional-dependency functions (skip if the suggested pkg is absent) ----

test_that("root_scape_metrics runs when landscapemetrics is available", {
  skip_if_not_installed("terra")
  skip_if_not_installed("landscapemetrics")
  # landscapemetrics emits a static "use check_landscape()" advisory for any
  # non-geographic raster; muffle ONLY that message so genuine warnings surface.
  res <- withCallingHandlers(
    root_scape_metrics(make_binary_spatraster(), select_layer = 1),
    warning = function(w) {
      if (grepl("check_landscape", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
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
  res <- estimate_rotation_shift(rgb, rgb, select_layer = 1)
  expect_true(all(c("depth", "rotation", "peak") %in% names(res)))
})

# ---- wrappers exercised with the bundled example data ----

test_that("root_turnover (tc) compares two timepoints", {
  skip_if_not_installed("terra")
  data(skl_Oulanka2023_Session01_T067)
  data(skl_Oulanka2023_Session03_T067)
  t1 <- terra::rast(skl_Oulanka2023_Session01_T067)
  t2 <- terra::rast(skl_Oulanka2023_Session03_T067)
  res <- suppressWarnings(root_turnover(t1, t2, method = "tc",
                                        tc_method = "kimura", dpi = 150, unit = "cm"))
  expect_true(is.data.frame(res) || is.list(res))
})

test_that("root_turnover (dpc) runs on a single multi-layer image", {
  skip_if_not_installed("terra")
  data(TurnoverDPC_data)
  img <- terra::rast(TurnoverDPC_data)
  res <- suppressWarnings(root_turnover(img, method = "dpc", dpi = 150, unit = "cm"))
  expect_true(is.data.frame(res) || is.list(res))
})

test_that("root_depth_metrics runs on a one-image directory", {
  skip_if_not_installed("terra")
  data(seg_Oulanka2023_Session01_T067)
  dir <- tempfile("seg"); dir.create(dir)
  terra::writeRaster(terra::rast(seg_Oulanka2023_Session01_T067),
                     file.path(dir, "T067.tif"), overwrite = TRUE)
  res <- suppressWarnings(root_depth_metrics(
    path_seg         = paste0(dir, "/"),   # root_depth_metrics expects a trailing slash
    tube_names       = "T067",
    session          = "test",
    insertion_angles = 45,        # non-zero -> valid depth-map geometry
    soil_starts      = 0,
    dpi              = 150,
    calc_root_length         = FALSE,
    calc_diameter_stats      = FALSE,
    calc_density_metrics     = FALSE,
    calc_distribution_indices = FALSE,
    calc_advanced_metrics    = FALSE,
    verbose          = FALSE
  ))
  expect_s3_class(res, "data.frame")
  expect_gt(nrow(res), 0)
})

# stitch_root_scans needs an overlapping scan SEQUENCE, which the bundled
# single-tube data cannot provide. Supply a directory of overlapping frames to
# enable this test.
test_that("stitch_root_scans end-to-end (needs a scan-sequence fixture)", {
  skip("Provide a directory of overlapping scans to enable this test.")
})

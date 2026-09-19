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
    path_seg         = dir,       # a plain directory, no trailing separator
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
  expect_true(any(res$rootpx > 0))          # images were actually read, not skipped
})


test_that("the default-on derived metric groups produce values", {
  # Both blocks used dplyr::.data$col, which raises "Can't subset .data outside
  # of a data mask context". The wrapper's own fault tolerance caught it and
  # NA'd the columns, so every user silently got an empty mrd and
  # total.length.density from a run that reported success.
  skip_if_not_installed("terra")
  data(seg_Oulanka2023_Session01_T067)
  # Downsampled 4x purely for speed; the assertions are about which columns get
  # computed, which does not depend on image size.
  small <- terra::aggregate(terra::rast(seg_Oulanka2023_Session01_T067),
                            fact = 4, fun = "max")
  dir <- tempfile("seg"); dir.create(dir)
  terra::writeRaster(small, file.path(dir, "T067.tif"), overwrite = TRUE)
  res <- suppressWarnings(root_depth_metrics(
    path_seg            = dir,
    tube_names          = "T067",
    insertion_angles    = 45,
    dpi                 = 150,
    depth_interval_cm   = 10,
    calc_diameter_stats = FALSE,   # needs imager
    verbose             = FALSE
  ))
  expect_true(all(c("mrd", "total.length.density") %in% names(res)))
  expect_false(all(is.na(res$mrd)))
  expect_false(all(is.na(res$total.length.density)))
})


test_that("depth_interval_cm = NULL collapses each scan to a single row", {
  # A tray of washed roots has no depth axis: the whole scan is one bin, and
  # the densities become whole-scan numbers (how densely the tray was packed).
  skip_if_not_installed("terra")
  data(flatbed_scan_example)
  # Downsampled 4x purely for speed; what is asserted does not depend on size.
  small <- terra::aggregate(terra::rast(flatbed_scan_example)[[1]],
                            fact = 4, fun = "max")
  dir <- tempfile("seg"); dir.create(dir)
  terra::writeRaster(small, file.path(dir, "tray_01.tif"), overwrite = TRUE)

  run <- function(...) suppressWarnings(root_depth_metrics(
    path_seg = dir, tube_names = "tray_01", dpi = 150,
    calc_diameter_stats = FALSE,        # needs imager
    verbose = FALSE, ...))

  whole <- run(depth_interval_cm = NULL)
  expect_equal(nrow(whole), 1L)
  expect_equal(whole$depth, 0)
  expect_gt(whole$rootlength.density, 0)
  # With one bin there is no profile shape left for these to describe.
  expect_false(any(c("mrd", "total.length.density", "rootlength.fraction")
                   %in% names(whole)))

  # The whole-scan densities are the profile's, pooled over its bins.
  profile <- run(depth_interval_cm = 5)
  expect_gt(nrow(profile), 1L)

  cover <- sum(profile$rootpx, na.rm = TRUE) /
    sum(profile$rootpx + profile$voidpx, na.rm = TRUE) * 100
  expect_equal(whole$rootpx.density, cover, tolerance = 1e-8)

  # Length gets a loose tolerance on purpose: terra::terrain(v = "flowdir")
  # breaks ties at random, so root length moves by a few tenths of a percent
  # between identical runs. The bin width must not move it by more than that.
  pooled <- sum(profile$rootlength, na.rm = TRUE) /
    (sum(profile$rootpx + profile$voidpx, na.rm = TRUE) / (150 / 2.54)^2)
  expect_equal(whole$rootlength.density, pooled, tolerance = 0.02)
})


test_that("duplicate tube names are refused", {
  # Non-unique names key the tube-level joins onto one tube and multiply the
  # rows out, so they are caught before any image is read.
  skip_if_not_installed("terra")
  data(flatbed_scan_example)
  small <- terra::aggregate(terra::rast(flatbed_scan_example)[[1]],
                            fact = 8, fun = "max")
  dir <- tempfile("seg"); dir.create(dir)
  # Both names end in "fine.tif", so the derived default collapses to "Tfin".
  for (f in c("CLAS1A10-15fine.tif", "CLAS1B25-30fine.tif"))
    terra::writeRaster(small, file.path(dir, f), overwrite = TRUE)

  run <- function(...) suppressWarnings(root_depth_metrics(
    path_seg = dir, dpi = 150, depth_interval_cm = NULL,
    calc_diameter_stats = FALSE, verbose = FALSE, ...))

  expect_error(run(), "unique")                          # derived names collide
  expect_error(run(tube_names = c("A", "A")), "unique")  # passed in twice
  expect_error(run(tube_names = "A"), "unique")          # one name recycled
  expect_equal(run(tube_names = c("A", "B"))$Tube, c("A", "B"))
})


test_that("order_scheme names its own columns, and spurs can be pruned", {
  skip_if_not_installed("terra")
  data(flatbed_scan_example)
  small <- terra::aggregate(terra::rast(flatbed_scan_example)[[1]],
                            fact = 4, fun = "max")
  dir <- tempfile("seg"); dir.create(dir)
  terra::writeRaster(small, file.path(dir, "tray.tif"), overwrite = TRUE)

  run <- function(...) suppressWarnings(root_depth_metrics(
    path_seg = dir, tube_names = "tray", dpi = 150, depth_interval_cm = NULL,
    calc_diameter_stats = FALSE,        # needs imager
    calc_root_order_metrics = TRUE, verbose = FALSE, ...))

  strahler <- run()                                  # the default
  branch   <- run(order_scheme = "branch_order")

  # A column can never be read under the wrong definition.
  expect_true(all(c("mean.strahler_order", "max.strahler_order") %in% names(strahler)))
  expect_false("mean.branch_order" %in% names(strahler))
  expect_true("mean.branch_order" %in% names(branch))

  expect_gte(strahler$lateral_root_fraction, 0)
  expect_lte(strahler$lateral_root_fraction, 1)

  # Pruning works on the skeleton every metric is measured from, so it takes
  # both tips and length away.
  pruned <- run(prune_spur_length_cm = 0.3, prune_spur_iter = 2)
  tips <- function(r) r$main_root.n_tips + r$lateral_roots.n_tips
  expect_lt(tips(pruned), tips(strahler))
  expect_lt(pruned$rootlength, strahler$rootlength)
})


test_that("rotation_fixed_width controls the rotation-axis crop", {
  # It used to be hardcoded to 1800, which is wider than any tube in the
  # bundled data, so the crop silently clamped to the full image every time.
  skip_if_not_installed("terra")
  data(seg_Oulanka2023_Session01_T067)
  # Downsampled 4x: this test runs the whole wrapper twice, and at full size
  # that alone doubled the runtime of this file. The crop is a row count, so a
  # smaller image tests it just as well.
  small <- terra::aggregate(terra::rast(seg_Oulanka2023_Session01_T067),
                            fact = 4, fun = "max")
  dir <- tempfile("seg"); dir.create(dir)
  terra::writeRaster(small, file.path(dir, "T067.tif"), overwrite = TRUE)
  run <- function(...) suppressWarnings(root_depth_metrics(
    path_seg = dir, dpi = 150, insertion_angles = 45, depth_interval_cm = 10,
    calc_diameter_stats = FALSE, verbose = FALSE, ...))

  wide   <- run()                             # default 1800 > image height -> clamps
  narrow <- run(rotation_fixed_width = 150)   # fits, so it really crops

  expect_lt(sum(narrow$rootpx, na.rm = TRUE), sum(wide$rootpx, na.rm = TRUE))
})

# stitch_root_scans needs an overlapping scan SEQUENCE, which the bundled
# single-tube data cannot provide. Supply a directory of overlapping frames to
# enable this test.
test_that("stitch_root_scans end-to-end (needs a scan-sequence fixture)", {
  skip("Provide a directory of overlapping scans to enable this test.")
})

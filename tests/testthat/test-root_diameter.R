# Leaf compute function: ONE representative input (it delegates input handling
# to load_flexible_image), but a parameter sweep over its own `unit` argument
# plus the important edge case (empty image must error cleanly, not return junk).

test_that("root_diameter returns the documented structure", {
  skip_if_not_installed("terra")
  skip_if_not_installed("imager")
  img <- make_binary_spatraster()
  res <- root_diameter(img, unit = "px", select_layer = 1)

  expect_type(res, "list")
  expect_true(all(c("mean_diameter", "median_diameter", "diameters",
                    "root_volume", "root_surface_area") %in% names(res)))
  expect_true(is.numeric(res$mean_diameter) && length(res$mean_diameter) == 1)
})

test_that("root_diameter runs across all unit settings", {
  skip_if_not_installed("terra")
  skip_if_not_installed("imager")
  img  <- make_binary_spatraster()
  grid <- expand.grid(unit = c("px", "cm", "inch"), stringsAsFactors = FALSE)
  expect_runs_on_grid(
    function(unit) root_diameter(img, unit = unit, select_layer = 1, dpi = 300),
    grid
  )
})

test_that("unit scaling is internally consistent (cm == px * 2.54/dpi)", {
  skip_if_not_installed("terra")
  skip_if_not_installed("imager")
  img <- make_binary_spatraster()
  px  <- root_diameter(img, unit = "px",  select_layer = 1)$mean_diameter
  cm  <- root_diameter(img, unit = "cm",  select_layer = 1, dpi = 300)$mean_diameter
  expect_equal(cm, px * 2.54 / 300, tolerance = 1e-6)
})

test_that("an empty (all-zero) image errors cleanly", {
  skip_if_not_installed("terra")
  skip_if_not_installed("imager")
  empty <- terra::rast(matrix(0, 24, 28))
  expect_error(root_diameter(empty, unit = "px", select_layer = 1))
})

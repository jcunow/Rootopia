# estimate_soil_surface: detects the soil/tape boundary from an RGB scan.
# Needs RStoolbox; uses the bundled RGB example for realistic tape markers.

test_that("estimate_soil_surface returns soil/tape positions", {
  skip_if_not_installed("terra")
  skip_if_not_installed("RStoolbox")
  data(rgb_Oulanka2023_Session03_T067)
  rgb <- terra::rast(rgb_Oulanka2023_Session03_T067)
  res <- suppressWarnings(estimate_soil_surface(rgb, dpi = 150))
  expect_true(is.data.frame(res) || is.numeric(res))
})

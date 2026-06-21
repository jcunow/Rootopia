# Smoke tests for the plotting / file-output functions: build valid inputs,
# render to a temp PNG, and assert the file is produced without error. We do not
# check the visual content (that is impractical) -- only that it runs end-to-end.

make_order_map_for_plots <- function() {
  skip_if_not_installed("terra")
  skel <- make_binary_spatraster()
  branch_order_map(skel = skel, mask = skel, order = "branch_order", unit = "px")
}

test_that("render_order_overlay writes a PNG", {
  res  <- make_order_map_for_plots()
  segs <- attr(res$edges, "segments")
  dims <- attr(res$edges, "dims")
  f <- tempfile(fileext = ".png")
  # tiny synthetic graph -> harmless size/empty-layer warnings; this is a smoke test
  suppressWarnings(render_order_overlay(segs, res$edges$tip_order, dims, file = f))
  expect_true(file.exists(f))
})

test_that("plot_order_window writes a PNG", {
  res <- make_order_map_for_plots()
  f <- tempfile(fileext = ".png")
  suppressWarnings(
    plot_order_window(res$edges, make_binary_spatraster(),
                      r_range = c(1, 24), c_range = c(1, 28), file = f)
  )
  expect_true(file.exists(f))
})

test_that("plot_soil_classification writes a PNG", {
  skip_if_not_installed("terra")
  rgb <- terra::rast(make_rgb_array() * 255)
  cls <- classify_soil_rgb(rgb, verbose = FALSE)
  f <- tempfile(fileext = ".png")
  suppressWarnings(plot_soil_classification(cls, save_png = f))
  expect_true(file.exists(f))
})

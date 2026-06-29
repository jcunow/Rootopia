# load_flexible_image() must return a canonical (H, W, C) array for every input
# source, so terra/imager/file paths no longer disagree on orientation.

# Asymmetric reference (H != W, distinct features on the first row vs first
# column) so any transpose is detectable by value, not just by dim.
make_ref <- function(H = 6L, W = 10L) {
  ref <- array(0, c(H, W, 3))
  ref[1, , 1] <- 1      # bright top row in the red channel
  ref[, 1, 2] <- 1      # bright left column in the green channel
  ref[H, W, 3] <- 1     # single marker in the bottom-right, blue channel
  ref
}

test_that("array and SpatRaster inputs keep (H, W, C)", {
  ref <- make_ref()
  a_arr <- load_flexible_image(ref, output_format = "array", scale = "none")
  expect_equal(dim(a_arr), dim(ref))
  expect_equal(a_arr, ref)

  sr   <- terra::rast(ref)
  a_sr <- load_flexible_image(sr, output_format = "array", scale = "none")
  expect_equal(dim(a_sr)[1:2], dim(ref)[1:2])
  expect_equal(a_sr, ref)
})

test_that("array -> cimg -> array round-trips without transposing", {
  ref  <- make_ref()
  cm   <- load_flexible_image(ref, output_format = "cimg",  scale = "none")
  back <- load_flexible_image(cm,  output_format = "array", scale = "none")
  expect_equal(dim(back)[1:2], dim(ref)[1:2])
  expect_equal(back, ref)
})

test_that("PNG and TIFF files load in the same (H, W, C) orientation", {
  skip_if_not_installed("png")
  ref <- make_ref()

  png_f <- tempfile(fileext = ".png")
  png::writePNG(ref, png_f)
  a_png <- load_flexible_image(png_f, output_format = "array", scale = "none")
  expect_equal(dim(a_png)[1:2], dim(ref)[1:2])
  expect_equal(round(a_png, 4), ref)            # 8-bit PNG stores 0/1 exactly

  skip_if_not_installed("tiff")
  tif_f <- tempfile(fileext = ".tif")
  tiff::writeTIFF(ref, tif_f)
  a_tif <- load_flexible_image(tif_f, output_format = "array", scale = "none")
  expect_equal(dim(a_tif)[1:2], dim(ref)[1:2])
  expect_equal(a_png, a_tif)                    # PNG and TIFF now agree
})

test_that("a PNG loaded to SpatRaster matches the array reference (terra view)", {
  skip_if_not_installed("png")
  ref   <- make_ref()
  png_f <- tempfile(fileext = ".png"); png::writePNG(ref, png_f)
  r     <- load_flexible_image(png_f, output_format = "spatrast", scale = "none")
  expect_equal(c(terra::nrow(r), terra::ncol(r)), dim(ref)[1:2])
  expect_equal(round(terra::as.array(r), 4), ref)
})

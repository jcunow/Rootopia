# Orientation check against real scanner files, supplied by whoever runs it.
#
# test-load-orientation.R proves the contract on synthetic images this package
# writes itself. That cannot catch what real files do differently: EXIF
# orientation flags on camera and phone JPEGs, the TIFF Orientation tag some
# scanners set, palette and alpha channels, 16-bit depth. Those are decided by
# the decoder, not by Rootopia, and they vary between machines.
#
# To run it, point the variable at a directory of your own images and run the
# tests. Without it the whole file skips, so `R CMD check` and CI are unaffected.
#
#   Sys.setenv(ROOTOPIA_IMAGE_FIXTURES = "~/rootopia-fixtures")
#   testthat::test_local(filter = "load-orientation-local")
#
# What to put in the directory: the same scene saved several ways -- PNG, TIFF,
# JPEG, 8-bit and 16-bit, grayscale and RGB, and if you can, one JPEG straight
# off a device that writes an EXIF orientation flag. Make every image
# NON-SQUARE. A transpose is invisible in the dimensions of a square image, and
# the dimension check is the one that catches it on any format.

fixture_dir <- Sys.getenv("ROOTOPIA_IMAGE_FIXTURES", "")

fixture_files <- if (nzchar(fixture_dir) && dir.exists(path.expand(fixture_dir))) {
  list.files(path.expand(fixture_dir), full.names = TRUE,
             pattern = "\\.(png|tif|tiff|jpg|jpeg)$", ignore.case = TRUE)
} else {
  character(0)
}

# An independent decode of the same file, in (H, W[, C]). png and tiff both
# return row-major arrays with row 1 at the top of the image, which is exactly
# the layout load_flexible_image() promises -- that makes them a reference the
# loader had no hand in producing. Other formats have no such reference here,
# so they get the weaker cross-format checks only.
reference_decode <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "png" && requireNamespace("png", quietly = TRUE)) {
    return(png::readPNG(path))
  }
  if (ext %in% c("tif", "tiff") && requireNamespace("tiff", quietly = TRUE)) {
    return(as.array(tiff::readTIFF(path)))
  }
  NULL
}

hw <- function(x) {
  if (inherits(x, "SpatRaster")) return(c(terra::nrow(x), terra::ncol(x)))
  dim(x)[1:2]
}


test_that("fixture directory is usable", {
  if (!nzchar(fixture_dir)) {
    skip("set ROOTOPIA_IMAGE_FIXTURES to a directory of your own scans")
  }
  expect_true(dir.exists(path.expand(fixture_dir)))
  expect_gt(length(fixture_files), 0)

  # A square fixture cannot fail the dimension check even when the loader
  # transposes it, so say which ones are carrying less weight than they look.
  for (f in fixture_files) {
    ref <- reference_decode(f)
    if (!is.null(ref) && dim(ref)[1] == dim(ref)[2]) {
      message("[fixture] ", basename(f),
              " is square -- a transpose is undetectable by dimension here.")
    }
  }
})


for (f in fixture_files) {

  test_that(paste0(basename(f), ": loads in (H, W, C) like an independent decoder"), {
    skip_if_not_installed("terra")
    ref <- reference_decode(f)
    if (is.null(ref)) {
      skip(paste0("no independent decoder for .", tolower(tools::file_ext(f))))
    }

    got <- load_flexible_image(f, output_format = "array", scale = "none")

    # The orientation assertion. A failure here on a non-square image means the
    # loader transposed it: height and width came back swapped.
    expect_equal(hw(got), dim(ref)[1:2])

    # Values, once the shape is known to be right. A failure here with the
    # dimensions passing is a decoder difference -- gamma, alpha handling,
    # bit depth -- not an orientation bug. Channel counts may legitimately
    # differ (alpha dropped), so compare the channels both sources share.
    nc <- min(dim(got)[3], dim(ref)[3], na.rm = TRUE)
    if (!is.na(nc) && nc >= 1) {
      expect_equal(as.numeric(got[, , seq_len(nc)]),
                   as.numeric(ref[, , seq_len(nc)]),
                   tolerance = 1e-3)
    }
  })

  test_that(paste0(basename(f), ": every output_format agrees on height and width"), {
    skip_if_not_installed("terra")
    skip_if_not_installed("imager")

    a  <- load_flexible_image(f, output_format = "array",    scale = "none")
    r  <- load_flexible_image(f, output_format = "spatrast", scale = "none")
    ci <- load_flexible_image(f, output_format = "cimg",     scale = "none")

    # cimg is (W, H, D, C) by imager's convention, so its own dims are
    # transposed on purpose. Feed it back through the loader: the round trip
    # has to return to the same (H, W) the array branch produced.
    back <- load_flexible_image(ci, output_format = "array", scale = "none")

    expect_equal(hw(r),    hw(a))
    expect_equal(hw(back), hw(a))
    expect_equal(as.numeric(back), as.numeric(a), tolerance = 1e-3)
  })
}

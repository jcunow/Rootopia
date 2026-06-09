# ============================================================
# PATCH for load_flexible_image  (in R/ImageLoader2.R or wherever it lives)
#
# The bug
# -------
# When a cimg or multi-layer array is converted to a SpatRaster, the RGB
# channel metadata is never registered.  terra requires you to call
#   terra::RGB(r) <- 1:3    # for a 3-band image
# before terra::plotRGB() will work.  Without this call, plotRGB() throws:
#   Error: [set.RGB] value must have length 3 or 4
# even when the raster genuinely has 3 layers.
#
# Fix location
# ------------
# Find the section of load_flexible_image that builds/returns a SpatRaster
# (likely around the lines that do terra::rast(...) or as(img, "SpatRaster"))
# and add the terra::RGB() call immediately after the raster is constructed.
#
# Drop-in helper (can be called anywhere in the package)
# -------------------------------------------------------


#' Convert a cimg object to a SpatRaster, preserving RGB metadata
#'
#' Internal helper that replaces the bare `terra::rast(as.array(img))` pattern
#' used throughout the package.  The key addition is `terra::RGB(r) <- ...`,
#' which registers the channel indices so that `terra::plotRGB()` works without
#' the user having to know about this terra requirement.
#'
#' @param img A `cimg` object (any number of channels).
#' @param normalize Logical.  If `TRUE`, rescale values to [0, 1].
#' @return A `SpatRaster` with RGB metadata set when the image has 3 or 4
#'   channels.
#' @keywords internal
cimg_to_spatrast <- function(img, normalize = FALSE) {
  if (!inherits(img, "cimg"))
    stop("cimg_to_spatrast: input must be a cimg object")

  n_channels <- imager::spectrum(img)   # 4th dimension of the cimg array

  # imager stores pixels as [width, height, depth, channels].
  # terra::rast() expects a matrix (1-layer) or a 3D array [nrow, ncol, nlyr]
  # in row-major order.  We must permute axes.
  arr <- as.array(img)   # [width, height, depth=1, channels]

  if (n_channels == 1L) {
    # Single-channel: produce a plain matrix for terra
    mat <- arr[, , 1, 1]
    r   <- terra::rast(t(mat))         # terra::rast on a matrix = 1-layer
  } else {
    # Multi-channel: stack channels as separate layers
    # arr dims: [W, H, 1, C] -> we want [H, W, C] for terra
    layers <- vector("list", n_channels)
    for (ch in seq_len(n_channels)) {
      layers[[ch]] <- terra::rast(t(arr[, , 1, ch]))
    }
    r <- terra::rast(layers)
  }

  if (normalize) {
    mx <- max(terra::values(r), na.rm = TRUE)
    if (!is.na(mx) && mx > 0) r <- r / mx
  }

  # Register RGB channels so terra::plotRGB() works without error.
  # This is the fix for: Error [set.RGB] value must have length 3 or 4
  if (n_channels >= 3L)
    terra::RGB(r) <- 1:3
  if (n_channels == 4L)
    terra::RGB(r) <- 1:4

  r
}


# ============================================================
# HOW TO APPLY THIS PATCH TO load_flexible_image
# ============================================================
#
# Search load_flexible_image for the spatrast / SpatRaster conversion block.
# It will look something like one of these patterns:
#
#   Pattern A (from cimg):
#     result <- terra::rast(as.array(img))
#
#   Pattern B (from array/matrix):
#     result <- terra::rast(img_array)
#
#   Pattern C (from RasterBrick):
#     result <- terra::rast(img)
#
# Replace each with a call to cimg_to_spatrast() OR add the RGB setter inline:
#
#   # Pattern A replacement:
#   result <- cimg_to_spatrast(img, normalize = normalize)
#
#   # Inline fix if you prefer not to use the helper:
#   result <- terra::rast(as.array(img))
#   if (terra::nlyr(result) >= 3) terra::RGB(result) <- 1:3
#
# The terra::RGB() call costs nothing and is safe to add unconditionally for
# any raster with 3+ layers.

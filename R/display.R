# ---------------------------------------------------------------------------
# Faithful image display helper
# ---------------------------------------------------------------------------

#' Display a root image without resampling artifacts
#'
#' \code{terra::plot()} subsamples any raster larger than its \code{maxcell}
#' limit (500,000 cells by default) by regular decimation -- it keeps every
#' n-th row and column.  For root images this silently drops one-pixel-wide
#' features such as skeletons, fine roots, and thin branch-order lines, so the
#' preview no longer reflects the data the metrics are computed on.
#'
#' This helper avoids that.  When the image is larger than \code{max_dim} it is
#' reduced with \emph{block aggregation} instead of decimation: each output
#' cell summarises a block of input cells, so a feature present anywhere in a
#' block is retained.  The reduced image is then drawn at full cell resolution
#' with no interpolation, giving a faithful preview at a controlled size.
#'
#' @param x Image to display: anything accepted by \code{load_flexible_image()}
#'   (\code{SpatRaster}, \code{RasterLayer}/\code{RasterBrick}, array, matrix,
#'   or a file path).
#' @param fun Block-aggregation function applied when the image is downsampled
#'   for display.  \code{"max"} (default) preserves thin bright features and is
#'   the right choice for binary masks, skeletons, and integer label maps (a
#'   root pixel anywhere in a block keeps the block "root").  Use \code{"mean"}
#'   for continuous maps (depth, diameter) and for RGB photographs.
#' @param max_dim Integer. Target length of the longer image side, in pixels,
#'   after downsampling.  The image is block-aggregated by an integer factor
#'   until its larger dimension is at or below this value.  Use \code{Inf} to
#'   draw at native resolution (no aggregation).  Default \code{1000}.
#' @param col,legend,main Passed to \code{terra::plot()} for single-layer
#'   images.  Ignored for images with 3 or more bands, which are drawn with
#'   \code{terra::plotRGB()}.
#' @param ... Further graphical arguments passed to \code{terra::plot()}.
#' @return Invisibly, the (possibly aggregated) \code{SpatRaster} that was drawn.
#' @keywords internal
#' @examples
#' \dontrun{
#'   # Binary mask or skeleton: keep thin lines with block-max (the default)
#'   show_root_image(seg_Oulanka2023_Session03_T067)
#'
#'   # Continuous map: average instead, so values are not biased upward
#'   dm <- create_depthmap(seg_Oulanka2023_Session03_T067, select_layer = 2)
#'   show_root_image(dm, fun = "mean")
#' }
show_root_image <- function(x, fun = "max", max_dim = 1000,
                            col = NULL, legend = TRUE, main = "", ...) {
  r <- load_flexible_image(x, output_format = "spatrast",
                           scale = "none", select_layer = NULL)

  # Block-aggregate (not decimate) only when the image exceeds max_dim, so thin
  # features survive the size reduction.
  big <- max(dim(r)[1:2])
  if (is.finite(max_dim) && big > max_dim) {
    fact <- as.integer(ceiling(big / max_dim))
    if (fact > 1L) r <- terra::aggregate(r, fact = fact, fun = fun, na.rm = TRUE)
  }

  if (terra::nlyr(r) >= 3) {
    # RGB(A): plotRGB needs the value range; 0-1 floats vs 0-255 integers.
    mx <- suppressWarnings(max(terra::minmax(r), na.rm = TRUE))
    terra::plotRGB(r, r = 1, g = 2, b = 3,
                   scale   = if (is.finite(mx) && mx > 1) 255 else 1,
                   maxcell = terra::ncell(r))
  } else {
    # maxcell = ncell disables terra's own subsampling; smooth = FALSE keeps
    # nearest-neighbour rendering (no interpolation across cell edges).
    # Forward `col` only when supplied -- passing col = NULL would override
    # terra's default palette with NULL and error inside terra::plot().
    args <- list(r, maxcell = terra::ncell(r), smooth = FALSE,
                 legend = legend, main = main, ...)
    if (!is.null(col)) args$col <- col
    do.call(terra::plot, args)
  }

  invisible(r)
}

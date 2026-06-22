# Native-resolution plotting for SpatRasters
#
# Faithful display of high-resolution scans. When a several-thousand-pixel scan
# is shrunk to fit a figure, the graphics device drops one-pixel-wide roots and
# skeletons. zoom_plot() draws a full-image overview (with the zoom box marked)
# plus a native-resolution magnified inset, so thin features survive at ~1:1.

#' Plot a SpatRaster with a magnified native-resolution inset
#'
#' Displays a high-resolution `SpatRaster` without losing thin features to
#' device down-sampling. By default it draws two panels: a full-image
#' **overview** with the zoomed region outlined, and a **magnified inset** of
#' that region rendered at (near) 1:1 (`maxcell = Inf`).
#'
#' @section Why this exists:
#' Root masks and skeletons are often one pixel wide. Fitting a multi-thousand
#' pixel scan into a small figure forces the graphics device to drop those
#' pixels, so roots vanish or look broken. The magnified inset keeps a fraction
#' of the image at full pixel resolution; the effective magnification is
#' `1 / frac` (e.g. `frac = 0.25` magnifies 4x).
#'
#' @param x A `SpatRaster`.
#' @param frac Fraction of each axis kept in the magnified inset (0-1).
#'   Magnification is `1 / frac`. Default `0.3`.
#' @param center Where to centre the inset. One of:
#'   \itemize{
#'     \item `"center"` (default) — geometric centre of the image;
#'     \item `"densest"` — the `frac`-sized window containing the most
#'           non-zero / non-`NA` pixels (i.e. the most root material);
#'     \item a length-2 numeric `c(fx, fy)` of *relative* coordinates in
#'           `[0, 1]`, where `c(0, 0)` is bottom-left and `c(1, 1)` top-right.
#'   }
#' @param overview Logical. Draw the full-image overview panel with the zoom
#'   box outlined. Default `TRUE`. Set `FALSE` to draw only the inset.
#' @param layer Integer or `NULL`. For multi-layer rasters, which layer the
#'   `"densest"` search uses. Default `NULL` (first layer). Does not subset what
#'   is plotted — RGB rasters still plot as RGB.
#' @param box_col Colour of the zoom-box rectangle on the overview. Default
#'   `"red"`.
#' @param main Title stem. Panels are titled `"<main> - overview"` and
#'   `"<main> - inset (Nx)"`.
#' @param maxcell Passed to [terra::plot()]/[terra::plotRGB()]. Default `Inf`
#'   (no down-sampling) so thin features survive.
#' @param ... Further arguments passed to the underlying terra plot call.
#' @return Invisibly, the `terra::ext()` of the inset region.
#' @export
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' r <- terra::rast(seg_Oulanka2023_Session01_T067)
#' \dontrun{
#' zoom_plot(r, frac = 0.25)                 # 4x inset, centred
#' zoom_plot(r, center = "densest")          # zoom where the roots are
#' zoom_plot(r, center = c(0.2, 0.8))        # upper-left region
#' zoom_plot(r, overview = FALSE)            # inset only
#' }
zoom_plot <- function(x,
                      frac     = 0.3,
                      center   = "center",
                      overview = TRUE,
                      layer    = NULL,
                      box_col  = "red",
                      main     = "",
                      maxcell  = Inf,
                      ...) {

  if (!inherits(x, "SpatRaster"))
    stop("zoom_plot() expects a SpatRaster. Convert with load_flexible_image(..., output_format = 'spatrast').",
         call. = FALSE)
  if (!is.numeric(frac) || frac <= 0 || frac > 1)
    stop("`frac` must be a number in (0, 1].", call. = FALSE)

  ev <- as.vector(terra::ext(x))                 # xmin, xmax, ymin, ymax
  w  <- ev[2] - ev[1]; h <- ev[4] - ev[3]
  hw <- frac * w / 2;  hh <- frac * h / 2

  # Resolve the inset centre (cx, cy) in map coordinates.
  if (is.character(center)) {
    center <- match.arg(center, c("center", "densest"))
    if (center == "center") {
      cx <- (ev[1] + ev[2]) / 2; cy <- (ev[3] + ev[4]) / 2
    } else {
      cen <- .densest_window(x, frac, layer)
      cx <- cen[1]; cy <- cen[2]
    }
  } else {
    if (!is.numeric(center) || length(center) != 2)
      stop("`center` must be \"center\", \"densest\", or numeric c(fx, fy).",
           call. = FALSE)
    cx <- ev[1] + center[1] * w
    cy <- ev[3] + center[2] * h
  }

  # Clamp so the inset stays inside the image.
  cx <- min(max(cx, ev[1] + hw), ev[2] - hw)
  cy <- min(max(cy, ev[3] + hh), ev[4] - hh)
  box <- terra::ext(cx - hw, cx + hw, cy - hh, cy + hh)

  is_rgb <- terra::nlyr(x) >= 3L && isTRUE(terra::has.RGB(x))
  draw <- function(obj, title) {
    if (is_rgb) terra::plotRGB(obj, maxcell = maxcell, ...)
    else        terra::plot(obj, maxcell = maxcell, main = title, ...)
    if (is_rgb) graphics::title(main = title)
  }

  if (overview) {
    draw(x, paste0(main, " - overview"))
    terra::lines(box, col = box_col, lwd = 2)
  }
  draw(terra::crop(x, box),
       sprintf("%s - inset (%gx)", main, round(1 / frac, 1)))

  invisible(box)
}


#' Locate the densest frac-sized window in a raster
#'
#' Returns the map-coordinate centre `c(cx, cy)` of the `frac`-sized window
#' containing the most non-zero / non-`NA` cells. Uses a coarse block scan so
#' it stays cheap on large rasters.
#' @keywords internal
.densest_window <- function(x, frac, layer = NULL) {
  r <- if (!is.null(layer)) x[[layer]] else x[[1]]
  m <- terra::as.matrix(r, wide = TRUE)          # rows = y (top-down), cols = x
  m <- !is.na(m) & m != 0                        # presence mask
  nr <- nrow(m); nc <- ncol(m)

  wr <- max(1L, round(frac * nr))                # window size in cells
  wc <- max(1L, round(frac * nc))
  # Step on a grid of ~10 positions per axis to keep it fast.
  rs <- unique(pmax(1L, round(seq(1, nr - wr + 1, length.out = 10))))
  cs <- unique(pmax(1L, round(seq(1, nc - wc + 1, length.out = 10))))

  best <- -1; br <- 1L; bc <- 1L
  cs_sum <- rbind(0, apply(m, 2, cumsum))        # column-cumulative for fast sums
  for (r0 in rs) for (c0 in cs) {
    r1 <- min(nr, r0 + wr - 1L); c1 <- min(nc, c0 + wc - 1L)
    s  <- sum(cs_sum[r1 + 1L, c0:c1] - cs_sum[r0, c0:c1])
    if (s > best) { best <- s; br <- r0; bc <- c0 }
  }

  # Window-centre cell -> map coordinates.
  ctr_row <- br + wr / 2; ctr_col <- bc + wc / 2
  ev <- as.vector(terra::ext(x))
  cx <- ev[1] + (ctr_col / nc) * (ev[2] - ev[1])
  cy <- ev[4] - (ctr_row / nr) * (ev[4] - ev[3])   # row 1 = top = ymax
  c(cx, cy)
}

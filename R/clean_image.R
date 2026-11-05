# Image cleaning functions
#
# Hole filling and speck removal work on 0/1 matrices and label connected
# components with imager::label(). Edge smoothing and cimg output use imager.
# load_flexible_image() handles conversion from SpatRaster / matrix / file path.


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Label the connected non-zero regions of a 0/1 matrix
#'
#' Only non-zero cells get a label; zero and \code{NA} cells are \code{NA}.
#'
#' @param m A 0/1 matrix.
#' @param directions \code{8} (diagonal neighbors connect) or \code{4}.
#' @return Integer matrix of region labels, \code{NA} outside the regions.
#' @keywords internal
#' @noRd
.label_regions <- function(m, directions) {
  inside <- !is.na(m) & m != 0
  # imager::label() numbers every connected region of equal value, the 0-valued
  # background included, so the cells outside the regions are blanked after.
  # as.cimg() maps matrix rows to x and columns to y, so the labels unroll in
  # the same cell order as m.
  lab <- imager::label(imager::as.cimg(inside * 1),
                       high_connectivity = directions == 8)
  lab <- matrix(as.integer(lab), nrow(m), ncol(m))
  lab[!inside] <- NA
  lab
}

#' Labels of the regions that touch the matrix border
#' @keywords internal
#' @noRd
.border_labels <- function(lab) {
  b <- c(lab[1, ], lab[nrow(lab), ], lab[, 1], lab[, ncol(lab)])
  unique(b[!is.na(b)])
}


#' Fill internal black holes in a binary image
#'
#' Sets to 1 every black region (value 0) of at most \code{max_size} pixels
#' that is enclosed by white (value 1), i.e. does not touch the image border.
#' Holes are 4-connected: two background pixels that touch only at a corner
#' belong to different holes, because the 8-connected foreground passes
#' between them.
#'
#' @param img A 0/1 matrix.
#' @param max_size Maximum hole size in pixels to fill.  `NULL` (default)
#'   fills nothing; `Inf` fills every hole.
#' @return The matrix with holes filled.
#' @keywords internal
fill_holes <- function(img, max_size = NULL) {
  if (is.null(max_size)) return(img)
  lab <- .label_regions((img == 0) * 1, directions = 4)
  lab[lab %in% .border_labels(lab)] <- NA        # touches the border: background
  sizes <- table(lab)
  fill  <- as.integer(names(sizes)[sizes <= max_size])
  img[!is.na(lab) & lab %in% fill] <- 1
  img
}


#' Remove small white objects from a binary image
#'
#' Sets to 0 every white region (value 1) of at most \code{max_size} pixels.
#' Regions are 8-connected: pixels that touch at a corner belong to the same
#' object, so a diagonal run of root pixels is one object. With
#' \code{protect_border = TRUE}, regions touching the image border are kept
#' whatever their size.
#'
#' @param img A 0/1 matrix.
#' @param max_size Maximum object size in pixels to remove.  `NULL` (default)
#'   removes nothing; `Inf` removes every candidate region.
#' @param protect_border Logical.  If `TRUE`, white regions touching the image
#'   border are never removed (they are assumed to be roots leaving the frame).
#'   If `FALSE` (default), border-touching regions are subject to the same
#'   `max_size` test as any other region, so edge specks are removed too.
#' @return The matrix with objects removed.
#' @keywords internal
remove_small_objects <- function(img, max_size = NULL, protect_border = FALSE) {
  if (is.null(max_size)) return(img)
  lab <- .label_regions((img > 0) * 1, directions = 8)
  if (protect_border) lab[lab %in% .border_labels(lab)] <- NA
  sizes <- table(lab)
  drop  <- as.integer(names(sizes)[sizes <= max_size])
  img[!is.na(lab) & lab %in% drop] <- 0
  img
}


#' Create a morphological structuring element
#'
#' Returns a `cimg` kernel of a given shape for use with
#' `imager::dilate()` / `imager::erode()`.
#'
#' @param shape One of `"disk"` (default), `"square"`, or `"diamond"`.
#' @param size Kernel size in pixels (odd integer).  Even values are silently
#'   incremented by 1.
#' @return A `cimg` kernel.
#' @keywords internal
create_kernel <- function(shape = "disk", size = 3) {
  valid <- c("square", "diamond", "disk")
  if (!shape %in% valid)
    stop(sprintf("Invalid kernel shape '%s'. Choose from: %s",
                 shape, paste(valid, collapse = ", ")),
         call. = FALSE)

  if (size %% 2 == 0) size <- size + 1L

  kern   <- matrix(0, size, size)
  center <- floor(size / 2) + 1L
  radius <- floor(size / 2)

  for (i in seq_len(size)) {
    for (j in seq_len(size)) {
      kern[i, j] <- switch(shape,
        square  = 1,
        diamond = as.integer(abs(i - center) + abs(j - center) <= radius),
        disk    = as.integer(sqrt((i - center)^2 + (j - center)^2) <= radius)
      )
    }
  }

  imager::as.cimg(kern)
}


# ---------------------------------------------------------------------------
# Exported diagnostic helper
# ---------------------------------------------------------------------------

#' Summarize the component-size distribution of a binary image
#'
#' Reports the size distribution of the two component types that \code{\link{clean_image}}
#' acts on, to help you choose `fill_holes_upto` and `drop_objects_upto`:
#' \itemize{
#'   \item **holes** — `background (0)` regions fully enclosed by `root (1)`
#'         (segmentation gaps inside roots). Regions touching the image border
#'         are background, not holes, and are excluded.
#'   \item **root components** — connected `root (1)` regions. The whole root
#'         system is usually one large component; genuine artifacts are the
#'         small components at the low end of this distribution.
#' }
#'
#' @section No size threshold here:
#' This function applies **no** size cutoff — it characterizes *every*
#' component. There is therefore nothing special separating a "big root" from an
#' "artifact": both are root components, distinguished only by size. Use the
#' printed summary and the histograms to pick the `fill_holes_upto` and
#' `drop_objects_upto` values you then pass to \code{\link{clean_image}}.
#' Counts here are in pixels, so they match `clean_image()`'s default
#' `size_unit = "px"` directly; divide by \eqn{(dpi/25.4)^2} to express one in
#' mm\eqn{^2}.
#'
#' @param img A binary image (`SpatRaster`, `cimg`, matrix, or file path); any
#'   format accepted by \code{\link{load_flexible_image}}.
#' @param plot Logical. Draw size-distribution histograms (log10 x-axis, with
#'   mean and median marked). Default `TRUE`.
#' @param breaks Number of histogram breaks. Default `30`.
#' @param select_layer Integer or `NULL`. Layer to use for multi-layer inputs.
#' @return Invisibly, a list with numeric vectors `holes` and `objects` giving
#'   the size (in pixels) of every enclosed hole and every root component.
#' @export
#' @examples
#' \dontrun{
#' data(seg_Oulanka2023_Session01_T067)
#' sizes <- report_image_components(seg_Oulanka2023_Session01_T067)
#' quantile(sizes$objects)   # inspect the small end to set drop_objects_upto
#' }
report_image_components <- function(img, plot = TRUE, breaks = 30,
                                    select_layer = NULL) {
  arr <- load_flexible_image(img, output_format = "array", scale = "binary",
                             select_layer = select_layer)
  if (length(dim(arr)) > 2 && dim(arr)[3] > 1) {
    warning(sprintf(paste("img has %d layers and select_layer is not set;",
                          "reporting layer 1."), dim(arr)[3]), call. = FALSE)
  }
  m <- if (length(dim(arr)) > 2) arr[, , 1] else arr

  holes   <- .component_sizes(m == 0, directions = 4, exclude_border = TRUE)
  objects <- .component_sizes(m > 0,  directions = 8, exclude_border = FALSE)

  .print_size_summary("Enclosed holes  (background 0 surrounded by root 1)", holes)
  .print_size_summary("Root components (connected root-1 regions)",          objects)

  if (isTRUE(plot) && (length(holes) > 0 || length(objects) > 0)) {
    op <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(op))
    .plot_size_hist(holes,   "Enclosed hole sizes",   breaks)
    .plot_size_hist(objects, "Root component sizes",  breaks)
  }

  invisible(list(holes = holes, objects = objects))
}

#' Connected-component sizes (px) of a logical matrix
#'
#' Labels the connected `TRUE` regions of `x` with the same connectivity
#' clean_image() uses -- 8 for root objects, 4 for holes -- and returns one
#' size per region. With `exclude_border = TRUE`, regions touching the image
#' edge are dropped (used for enclosed holes).
#' @keywords internal
.component_sizes <- function(x, directions, exclude_border) {
  lab <- .label_regions(x * 1, directions = directions)
  if (exclude_border) lab[lab %in% .border_labels(lab)] <- NA
  if (all(is.na(lab))) return(integer(0))
  as.integer(table(lab))
}

#' Print a one-line size summary (n / min / median / mean / max)
#' @keywords internal
.print_size_summary <- function(title, sizes) {
  cat(sprintf("=== %s ===\n", title))
  if (length(sizes) == 0) { cat("  none\n\n"); return(invisible()) }
  cat(sprintf("  n = %d | min = %d | median = %g | mean = %.1f | max = %d  (px)\n\n",
              length(sizes), min(sizes), stats::median(sizes),
              mean(sizes), max(sizes)))
  invisible()
}

#' Histogram of component sizes on a log10 axis, with mean/median marked
#' @keywords internal
.plot_size_hist <- function(sizes, main, breaks) {
  if (length(sizes) == 0) {
    graphics::plot.new(); graphics::title(main = paste0(main, " (none)"))
    return(invisible())
  }
  lg <- log10(pmax(sizes, 1))
  graphics::hist(lg, breaks = breaks, main = main,
                 xlab = "component size (px, log10)",
                 col = "grey80", border = "white")
  graphics::abline(v = log10(mean(sizes)),          col = "red",  lwd = 2)
  graphics::abline(v = log10(stats::median(sizes)), col = "blue", lwd = 2, lty = 2)
  graphics::legend("topright", bty = "n",
                   legend = c(sprintf("mean = %.1f px",   mean(sizes)),
                              sprintf("median = %g px", stats::median(sizes))),
                   col = c("red", "blue"), lwd = 2, lty = c(1, 2))
  invisible()
}


# ---------------------------------------------------------------------------
# Edge smoothing (internal -- not exported)
# ---------------------------------------------------------------------------

#' Smooth object edges with morphological closing
#'
#' Applies a morphological closing (dilation then erosion) to smooth the
#' edges of binary objects.  Intended as a post-cleaning step before
#' skeletonization; do **not** apply after skeletonization.
#'
#' @param img A `cimg` binary image or any format accepted by
#'   \code{\link{load_flexible_image}}.
#' @param kernel_shape One of `"disk"` (default), `"square"`, `"diamond"`.
#' @param kernel_size Structuring element size (odd integer).  At 300 DPI,
#'   `kernel_size = 3` is a good starting point.
#' @param iterations Number of closing iterations.
#' @return A `cimg` binary image.
#' @importFrom imager dilate erode grayscale
#' @keywords internal
smooth_root_edges <- function(img,
                               kernel_shape = "disk",
                               kernel_size  = 3,
                               iterations   = 1) {
  img <- load_flexible_image(img, output_format = "cimg", scale = "binary")

  if (imager::spectrum(img) > 1)
    img <- imager::grayscale(img)

  kern <- create_kernel(shape = kernel_shape, size = kernel_size)

  for (i in seq_len(iterations)) {
    img <- imager::dilate(img, kern)
    img <- imager::erode(img, kern)
  }

  img
}


# ---------------------------------------------------------------------------
# Main exported function
# ---------------------------------------------------------------------------

#' Clean a binary root image
#'
#' Performs three sequential cleaning operations on a binary segmented image:
#' 1. **Hole filling** -- fills black regions enclosed by white (segmentation
#'    gaps inside roots) up to `fill_holes_upto` pixels.
#' 2. **Speck removal** -- deletes foreground blobs (false-positive root
#'    detections) up to `drop_objects_upto` pixels; blobs touching the image
#'    border are included unless `protect_border = TRUE`.
#' 3. **Edge smoothing** *(optional, off by default)* -- applies morphological
#'    closing to smooth jagged root edges.
#'
#' Each step runs only when its argument is set: with the defaults
#' (`NULL`, `NULL`, `edge_smooth = FALSE`) the image is binarized and returned
#' otherwise unchanged.
#'
#' @section Why clean before skeletonization:
#' `skeletonize_image()` thins the mask to a one-pixel line (Zhang-Suen
#' thinning). Thinning keeps the topology of the mask, so a small hole inside a
#' root becomes a loop in the skeleton, which reads as two spurious branch
#' points. An isolated artifact pixel becomes a phantom skeleton segment.
#' Cleaning first yields a much cleaner skeleton.
#'
#' @section Choosing thresholds:
#' Call \code{\link{report_image_components}} on your image first to see the actual
#' pixel counts of all holes and specks.  At 300 dpi, sensible starting values
#' are `fill_holes_upto = 50` and `drop_objects_upto = 10` pixels.
#'
#' A pixel threshold means a different physical size at every scan resolution:
#' 50 px is 0.36 mm\eqn{^2} at 300 dpi but 0.022 mm\eqn{^2} at 1200 dpi, so a
#' call tuned on one dataset cleans differently on another.  Set
#' `size_unit = "mm2"` and supply `dpi` to give the thresholds a fixed physical
#' meaning that carries across resolutions.
#'
#' @section Edge smoothing caution:
#' `edge_smooth = TRUE` applies a morphological closing that slightly dilates
#' then erodes root edges.  This can merge closely adjacent roots and alter
#' root diameter measurements.  Only use it when the segmentation output has
#' very jagged edges; leave it off (`FALSE`, the default) otherwise.
#'
#' @section Input must be binary:
#' `clean_image()` operates on a binary mask (root = 1, background = 0). If your
#' input is a raw probability / grayscale image from a segmentation model,
#' binarize it first with \code{\link{image_threshold}} and pass the result here.
#'
#' @param img A `cimg` object, `SpatRaster`, matrix, or file path (binary mask).
#' @param fill_holes_upto Largest enclosed hole that gets filled; anything
#'   bigger is left alone.  `NULL` (default) fills nothing; `Inf` fills every
#'   enclosed hole.  Units are set by `size_unit`.  See **Choosing
#'   thresholds** above.
#' @param drop_objects_upto Largest foreground blob that gets deleted;
#'   anything bigger is kept.  `NULL` (default) deletes nothing; `Inf` deletes
#'   every blob, roots included.  Units are set by `size_unit`.
#' @param size_unit Unit of the two thresholds above: `"px"` (default, pixel
#'   count) or `"mm2"` (square millimeters, which requires `dpi`).
#' @param dpi Scan resolution, required when `size_unit = "mm2"` and ignored
#'   otherwise.
#' @param max_hole_size,max_artifact_size Deprecated names for
#'   `fill_holes_upto` and `drop_objects_upto`.  Supplying either warns and
#'   sets the corresponding new argument.
#' @param protect_border Logical.  If `TRUE`, white regions touching the image
#'   border are kept regardless of size (roots leaving the frame).  If `FALSE`
#'   (default), border-touching specks are removed by the same
#'   `drop_objects_upto` test as anything else -- being at the edge does not
#'   protect an artifact.
#' @param edge_smooth Logical.  Apply morphological closing after hole/artifact
#'   cleaning.  Default `FALSE`.
#' @param kernel_shape Structuring element shape for edge smoothing:
#'   `"disk"` (default), `"square"`, or `"diamond"`.
#' @param kernel_size Structuring element size (odd integer).  Default `3`.
#' @param iterations Number of closing iterations for edge smoothing.
#'   Default `1`.
#' @param select_layer Integer or `NULL`.  Which layer to use for multi-layer
#'   inputs.
#' @param output_format Character.  Format of the returned object.  One of
#'   `"spatrast"` (default), `"cimg"`, or `"matrix"`.  Using `"spatrast"`
#'   means the result can be passed directly to `terra::plot()`,
#'   `skeletonize_image()`, `root_length()`, etc. without any further
#'   conversion.
#' @param report Logical.  If `TRUE`, also calls \code{\link{report_image_components}}
#'   on the *original* (uncleaned) image before cleaning.  When
#'   `output_format = "spatrast"` (default), the cleaned raster is returned
#'   directly even when `report = TRUE`; the report is printed as a side
#'   effect.  Default `FALSE`.
#' @return A cleaned image in the format specified by `output_format`
#'   (`SpatRaster`, `cimg`, or matrix).
#' @export
#' @seealso \code{\link{report_image_components}}, \code{\link{skeletonize_image}}, \code{\link{image_threshold}}
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img <- terra::rast(seg_Oulanka2023_Session01_T067)
#'
#'
#' # Clean: fill small holes, delete tiny specks -- returns SpatRaster
#' cleaned <- clean_image(img,
#'                        fill_holes_upto   = 50,
#'                        drop_objects_upto = 10,
#'                        select_layer      = 2)
#'
#' # The same thresholds expressed physically, so they mean the same thing
#' # whatever the scan resolution
#' cleaned_mm <- clean_image(img,
#'                           fill_holes_upto   = 0.36,
#'                           drop_objects_upto = 0.07,
#'                           size_unit         = "mm2",
#'                           dpi               = 300,
#'                           select_layer      = 2)
#'
#' # If you need a cimg for further imager operations:
#' cleaned_cimg <- clean_image(img, fill_holes_upto = 50,
#'                             output_format = "cimg", select_layer = 2)
clean_image <- function(img,
                         fill_holes_upto   = NULL,
                         drop_objects_upto = NULL,
                         size_unit         = c("px", "mm2"),
                         dpi               = NULL,
                         protect_border    = FALSE,
                         edge_smooth       = FALSE,
                         kernel_shape      = "disk",
                         kernel_size       = 3,
                         iterations        = 1,
                         select_layer      = NULL,
                         output_format     = "spatrast",
                         report            = FALSE,
                         max_hole_size     = NULL,
                         max_artifact_size = NULL) {

  output_format <- match.arg(output_format, c("spatrast", "cimg", "matrix"))
  size_unit     <- match.arg(size_unit)

  # Former argument names. They said how big a thing was but not what happened
  # to it, so a reader could not tell whether the value was the largest thing
  # kept or the largest thing removed.
  if (!is.null(max_hole_size)) {
    warning("`max_hole_size` is deprecated; use `fill_holes_upto`.", call. = FALSE)
    if (is.null(fill_holes_upto)) fill_holes_upto <- max_hole_size
  }
  if (!is.null(max_artifact_size)) {
    warning("`max_artifact_size` is deprecated; use `drop_objects_upto`.", call. = FALSE)
    if (is.null(drop_objects_upto)) drop_objects_upto <- max_artifact_size
  }

  # A threshold in pixels means a different physical size at every scan
  # resolution, so the same call cleans differently at 300 and at 1200 dpi.
  # Giving the size in mm^2 fixes the physical meaning and converts here.
  if (size_unit == "mm2") {
    if (is.null(dpi) || !is.numeric(dpi) || dpi <= 0)
      stop("`dpi` is required, and must be positive, when size_unit = \"mm2\".")
    px_per_mm2 <- (dpi / 25.4)^2
    if (!is.null(fill_holes_upto))   fill_holes_upto   <- fill_holes_upto   * px_per_mm2
    if (!is.null(drop_objects_upto)) drop_objects_upto <- drop_objects_upto * px_per_mm2
  }

  # (H, W) for one layer, (H, W, C) for several; each layer is cleaned on its own.
  arr <- load_flexible_image(img,
                             output_format = "array",
                             select_layer  = select_layer,
                             scale         = "binary")

  if (report) report_image_components(arr)

  clean_layer <- function(m) {
    m <- fill_holes(m, fill_holes_upto)
    remove_small_objects(m, drop_objects_upto, protect_border = protect_border)
  }
  if (length(dim(arr)) > 2) {
    for (k in seq_len(dim(arr)[3])) arr[, , k] <- clean_layer(arr[, , k])
  } else {
    arr <- clean_layer(arr)
  }

  if (edge_smooth) {
    arr <- load_flexible_image(
      smooth_root_edges(arr, kernel_shape, kernel_size, iterations),
      output_format = "array", scale = "none")
  }

  # Convert to the requested output format
  switch(output_format,
    spatrast = {
      r <- load_flexible_image(arr, output_format = "spatrast", scale = "none")
      # RGB metadata lets terra::plotRGB() draw a multi-layer result; for a
      # single-layer binary image, use terra::plot().
      if (terra::nlyr(r) >= 3L) terra::RGB(r) <- 1:3
      r
    },
    cimg   = load_flexible_image(arr, output_format = "cimg",   scale = "none"),
    matrix = load_flexible_image(arr, output_format = "matrix", scale = "none")
  )
}










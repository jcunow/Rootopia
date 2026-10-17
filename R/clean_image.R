# Image cleaning functions
#
# Dependencies: imager (for connected-component labeling, morphological ops)
# All functions work on `cimg` objects internally; load_flexible_image() handles
# conversion from SpatRaster / matrix / file path.


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Fill internal black holes in a binary image
#'
#' Sets to 1 all black regions (value = 0) that are completely surrounded by
#' white (value = 1) -- i.e. not connected to the image border.
#'
#' @param img A `cimg` binary image (values 0 and 1).
#' @param max_size Maximum hole size in pixels to fill.  If `NULL` (default),
#'   all holes are filled regardless of size.
#' @return A `cimg` with holes filled.
#' @keywords internal
fill_holes <- function(img, max_size = NULL) {
  inv  <- 1 - img
  lbl  <- imager::label(inv)
  nd   <- dim(lbl)

  border_labels <- unique(c(
    lbl[1,    , 1, 1],
    lbl[nd[1],, 1, 1],
    lbl[, 1,   1, 1],
    lbl[, nd[2], 1, 1]
  ))
  border_labels <- border_labels[border_labels > 0]

  lblv <- as.numeric(lbl)
  internal_labels <- lblv * as.numeric(!(lblv %in% border_labels))

  if (!is.null(max_size)) {
    sizes          <- table(internal_labels[internal_labels > 0])
    allowed        <- as.integer(names(sizes[sizes <= max_size]))
    fill_mask      <- internal_labels %in% allowed
  } else {
    fill_mask <- internal_labels > 0
  }

  img[fill_mask] <- 1
  img
}


#' Remove small isolated white artifacts from a binary image
#'
#' Sets to 0 all white regions (value = 1) that are not connected to the
#' image border.
#'
#' @param img A `cimg` binary image (values 0 and 1).
#' @param max_size Maximum artifact size in pixels to remove.  If `NULL`,
#'   all candidate white regions are removed.
#' @param protect_border Logical.  If `TRUE`, white regions touching the image
#'   border are never removed (they are assumed to be roots leaving the frame).
#'   If `FALSE` (default), border-touching regions are subject to the same
#'   `max_size` test as any other region, so edge specks are removed too.
#'   Border location confers no protection on its own.
#' @return A `cimg` with artifacts removed.
#' @keywords internal
remove_small_objects <- function(img, max_size = NULL, protect_border = FALSE) {
  lbl <- imager::label(img)
  nd  <- dim(lbl)

  lblv <- as.numeric(lbl)

  if (protect_border) {
    border_labels <- unique(c(
      lbl[1,    , 1, 1],
      lbl[nd[1],, 1, 1],
      lbl[, 1,   1, 1],
      lbl[, nd[2], 1, 1]
    ))
    border_labels <- border_labels[border_labels > 0]
    candidate <- lblv * as.numeric(!(lblv %in% border_labels))
  } else {
    candidate <- lblv                       # every white region is a candidate
  }

  if (!is.null(max_size)) {
    sizes        <- table(candidate[candidate > 0])
    small        <- as.integer(names(sizes[sizes <= max_size]))
    remove_mask  <- candidate %in% small
  } else {
    remove_mask <- candidate > 0
  }

  img[remove_mask] <- 0
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
  img <- load_flexible_image(img, output_format = "cimg", scale = "binary",
                             select_layer = select_layer)

  holes   <- .component_sizes(1 - img, exclude_border = TRUE)   # enclosed background (0)
  objects <- .component_sizes(img,     exclude_border = FALSE)  # root (1) components

  .print_size_summary("Enclosed holes  (background 0 surrounded by root 1)", holes)
  .print_size_summary("Root components (connected root-1 regions)",          objects)

  if (isTRUE(plot) && (length(holes) > 0 || length(objects) > 0)) {
    op <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(op))
    .plot_size_hist(holes,   "Enclosed hole sizes",   breaks)
    .plot_size_hist(objects, "Root component sizes",  breaks)
  }

  invisible(list(holes = holes, objects = objects))
}

#' Connected-component sizes (px) of a binary cimg
#'
#' Labels the connected `1`-regions of `bin_cimg` and returns one size per
#' component. With `exclude_border = TRUE`, components touching the image edge
#' are dropped (used for enclosed holes).
#' @keywords internal
.component_sizes <- function(bin_cimg, exclude_border) {
  lbl <- imager::label(bin_cimg)
  v   <- as.integer(lbl)
  keep <- v > 0
  if (exclude_border) {
    nd <- dim(lbl)
    bl <- unique(c(lbl[1, , 1, 1], lbl[nd[1], , 1, 1],
                   lbl[, 1, 1, 1], lbl[, nd[2], 1, 1]))
    bl <- bl[bl > 0]
    keep <- keep & !(v %in% bl)
  }
  if (!any(keep)) return(integer(0))
  as.integer(table(v[keep]))
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
#' skeletonisation; do **not** apply after skeletonisation.
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
#'    gaps inside roots).
#' 2. **Speck removal** -- deletes isolated foreground blobs not connected to
#'    the image border (false-positive root detections).
#' 3. **Edge smoothing** *(optional, off by default)* -- applies morphological
#'    closing to smooth jagged root edges.
#'
#' @section Why clean before skeletonisation:
#' `skeletonize_image()` uses the Medial Axis Transform, which is driven by
#' the distance transform.  Small holes inside a root inflate the local
#' distance values and force the medial axis to bifurcate around the hole,
#' producing spurious branching points.  Isolated artifact pixels produce
#' phantom skeleton segments.  Cleaning first yields a much cleaner skeleton.
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
#'   bigger is left alone.  `NULL` (default) fills every enclosed hole.  Units
#'   are set by `size_unit`.  See **Choosing thresholds** above.
#' @param drop_objects_upto Largest disconnected foreground blob that gets
#'   deleted; anything bigger is kept.  `NULL` (default) deletes every
#'   candidate blob.  Units are set by `size_unit`.
#' @param size_unit Unit of the two thresholds above: `"px"` (default, pixel
#'   count) or `"mm2"` (square millimetres, which requires `dpi`).
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

  img_cimg <- load_flexible_image(img,
                                   output_format = "cimg",
                                   select_layer  = select_layer,
                                   scale         = "binary")

  if (report) report_image_components(img_cimg)

  img_filled  <- fill_holes(img_cimg, fill_holes_upto)
  img_cleaned <- remove_small_objects(img_filled, drop_objects_upto,
                                      protect_border = protect_border)

  img_out <- if (edge_smooth) {
    smooth_root_edges(img_cleaned, kernel_shape, kernel_size, iterations)
  } else {
    img_cleaned
  }

  # Convert to the requested output format
  switch(output_format,
    spatrast = {
      r <- load_flexible_image(img_out,
                                output_format = "spatrast",
                                scale         = "none")
      # Ensure RGB metadata is set so terra::plotRGB() works on multi-layer results.
      # For a single-layer binary image, use terra::plot() rather than plotRGB().
      if (terra::nlyr(r) >= 3L) terra::RGB(r) <- 1:3
      r
    },
    cimg   = img_out,
    matrix = as.matrix(img_out)
  )
}










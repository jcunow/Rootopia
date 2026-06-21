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

  internal_labels <- lbl * as.numeric(!(lbl %in% border_labels))

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
#' @param max_size Maximum artifact size in pixels to remove.  If `NULL`
#'   (default), all isolated white regions are removed.
#' @return A `cimg` with artifacts removed.
#' @keywords internal
remove_small_objects <- function(img, max_size = NULL) {
  lbl <- imager::label(img)
  nd  <- dim(lbl)

  border_labels <- unique(c(
    lbl[1,    , 1, 1],
    lbl[nd[1],, 1, 1],
    lbl[, 1,   1, 1],
    lbl[, nd[2], 1, 1]
  ))
  border_labels <- border_labels[border_labels > 0]

  internal_labels <- lbl * as.numeric(!(lbl %in% border_labels))

  if (!is.null(max_size)) {
    sizes        <- table(internal_labels[internal_labels > 0])
    small        <- as.integer(names(sizes[sizes <= max_size]))
    remove_mask  <- internal_labels %in% small
  } else {
    remove_mask <- internal_labels > 0
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

#' Report sizes of holes and isolated artifacts in a binary image
#'
#' Prints a human-readable summary of all internal black holes and isolated
#' white artifacts in `img`, together with their pixel counts.  Use this
#' **before** calling [clean_image()] to decide on appropriate values for
#' `max_hole_size` and `max_artifact_size`.
#'
#' @section Choosing thresholds:
#' At 300 DPI a single root cross-section in a minirhizotron scan is
#' roughly 5–150 px^2 depending on root diameter.  As a starting point:
#' \itemize{
#'   \item `max_artifact_size = 10` removes single-pixel noise and tiny
#'         segmentation specks while preserving fine roots.
#'   \item `max_hole_size = 50` fills small gaps inside roots without
#'         merging genuinely separate objects.
#' }
#' Scale these linearly if your scanner DPI differs (e.g. at 150 DPI,
#' halve both values).
#'
#' @param img A `cimg` binary image (values 0 and 1), or any format accepted
#'   by [load_flexible_image()].
#' @return Invisibly `NULL`.  Prints to the console.
#' @keywords internal
#' @examples
#' \dontrun{
#' img <- imager::as.cimg(matrix(0, 50, 50))
#' img[10:20, 10:20] <- 1   # white square
#' img[13:15, 13:15] <- 0   # hole inside it
#' img[40, 40]        <- 1  # isolated artifact
#' report_image_components(img)
#' }
report_image_components <- function(img) {
  img <- load_flexible_image(img, output_format = "cimg", scale = "binary")

  cat("=== HOLES (black regions enclosed by white) ===\n")
  lbl_h   <- imager::label(1 - img)
  nd      <- dim(lbl_h)
  bl_h    <- unique(c(lbl_h[1,,1,1], lbl_h[nd[1],,1,1],
                       lbl_h[,1,1,1], lbl_h[,nd[2],1,1]))
  bl_h    <- bl_h[bl_h > 0]
  int_h   <- lbl_h * as.numeric(!(lbl_h %in% bl_h))
  ids_h   <- unique(as.numeric(int_h[int_h > 0]))

  if (length(ids_h) == 0) {
    cat("  No holes found.\n")
  } else {
    for (id in ids_h)
      cat(sprintf("  Hole %d: %d px\n", id,
                  sum(img[int_h == id] == 0, na.rm = TRUE)))
  }

  cat("\n=== ARTIFACTS (isolated white regions not touching border) ===\n")
  lbl_a   <- imager::label(img)
  bl_a    <- unique(c(lbl_a[1,,1,1], lbl_a[nd[1],,1,1],
                       lbl_a[,1,1,1], lbl_a[,nd[2],1,1]))
  bl_a    <- bl_a[bl_a > 0]
  int_a   <- lbl_a * as.numeric(!(lbl_a %in% bl_a))
  ids_a   <- unique(as.numeric(int_a[int_a > 0]))

  if (length(ids_a) == 0) {
    cat("  No isolated artifacts found.\n")
  } else {
    for (id in ids_a)
      cat(sprintf("  Artifact %d: %d px\n", id,
                  sum(img[int_a == id] == 1, na.rm = TRUE)))
  }

  invisible(NULL)
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
#'   [load_flexible_image()].
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
#' 2. **Artifact removal** — removes isolated white specks not connected to
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
#' Call [report_image_components()] on your image first to see the actual
#' pixel counts of all holes and artifacts.  At 300 DPI, sensible starting
#' values are `max_hole_size = 50` and `max_artifact_size = 10`.
#'
#' @section Edge smoothing caution:
#' `edge_smooth = TRUE` applies a morphological closing that slightly dilates
#' then erodes root edges.  This can merge closely adjacent roots and alter
#' root diameter measurements.  Only use it when the segmentation output has
#' very jagged edges; leave it off (`FALSE`, the default) otherwise.
#'
#' @section Optional pre-thresholding:
#' If the input is not yet a clean binary mask (e.g. a raw probability /
#' grayscale image from a segmentation model), set `pre_threshold` to binarize
#' it via [image_threshold()] *before* hole-filling and artifact removal. This
#' runs first because [image_threshold()] expects a non-binary `SpatRaster`;
#' once thresholded, the result feeds into the same hole-filling /
#' artifact-removal / edge-smoothing steps as a normal binary mask.
#'
#' @param img A `cimg` object, `SpatRaster`, matrix, or file path.
#' @param pre_threshold Numeric (0-1) or `NULL` (default). If not `NULL`,
#'   [image_threshold()] is applied to `img` first, using this value as its
#'   `threshold` argument, before any hole-filling or artifact removal.
#' @param pre_threshold_method Thresholding method passed to
#'   [image_threshold()] when `pre_threshold` is set: `"global"` (default) or
#'   `"adaptive"`.
#' @param pre_threshold_window_size Window size passed to [image_threshold()]
#'   when `pre_threshold_method = "adaptive"`. Default `15`.
#' @param max_hole_size Maximum hole size in pixels to fill.  If `NULL`, all
#'   enclosed holes are filled.  See **Choosing thresholds** above.
#' @param max_artifact_size Maximum artifact size in pixels to remove.  If
#'   `NULL`, all isolated white regions are removed.
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
#' @param report Logical.  If `TRUE`, also calls [report_image_components()]
#'   on the *original* (uncleaned) image before cleaning.  When
#'   `output_format = "spatrast"` (default), the cleaned raster is returned
#'   directly even when `report = TRUE`; the report is printed as a side
#'   effect.  Default `FALSE`.
#' @return A cleaned image in the format specified by `output_format`
#'   (`SpatRaster`, `cimg`, or matrix).
#' @export
#' @seealso [report_image_components()], [skeletonize_image()], [image_threshold()]
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img <- terra::rast(seg_Oulanka2023_Session01_T067)
#'
#'
#' # Clean: fill small holes, remove tiny artifacts -- returns SpatRaster
#' cleaned <- clean_image(img,
#'                        max_hole_size     = 50,
#'                        max_artifact_size = 10,
#'                        select_layer      = 2)
#'
#' # If you need a cimg for further imager operations:
#' cleaned_cimg <- clean_image(img, max_hole_size = 50,
#'                             output_format = "cimg", select_layer = 2)
clean_image <- function(img,
                         pre_threshold              = NULL,
                         pre_threshold_method       = "global",
                         pre_threshold_window_size  = 15,
                         max_hole_size     = NULL,
                         max_artifact_size = NULL,
                         edge_smooth       = FALSE,
                         kernel_shape      = "disk",
                         kernel_size       = 3,
                         iterations        = 1,
                         select_layer      = NULL,
                         output_format     = "spatrast",
                         report            = FALSE) {

  output_format <- match.arg(output_format, c("spatrast", "cimg", "matrix"))

  if (!is.null(pre_threshold)) {
    img <- load_flexible_image(img,
                                output_format = "spatrast",
                                select_layer  = select_layer,
                                scale         = "none")
    img <- image_threshold(img,
                            threshold   = pre_threshold,
                            method      = pre_threshold_method,
                            window_size = pre_threshold_window_size,
                            select_layer = NULL,
                            mask_layer   = NULL,
                            binary_01    = TRUE,
                            deblur       = FALSE)
    select_layer <- NULL
  }

  img_cimg <- load_flexible_image(img,
                                   output_format = "cimg",
                                   select_layer  = select_layer,
                                   scale         = "binary")

  if (report) report_image_components(img_cimg)

  img_filled  <- fill_holes(img_cimg, max_hole_size)
  img_cleaned <- remove_small_objects(img_filled, max_artifact_size)

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










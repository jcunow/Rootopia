####### create the depthmap

#' Assign a soil depth in centimeters to every pixel
#'
#' Builds a per-pixel depth map from the imaging geometry. \code{geometry}
#' selects the shape being imaged, and the shape decides which of the two
#' corrections apply: foreshortening of the depth axis by the insertion angle,
#' and the depth spread across a curved surface.
#'
#' @section Geometries:
#' \describe{
#'   \item{\code{"tube"}}{A cylindrical minirhizotron tube, unrolled by the
#'     scanner. Depth increases along the image width and is scaled by
#'     \code{cos(angle_from_vertical)}; a curvature wave spanning one
#'     circumference is added across the image height, scaled by
#'     \code{sin(angle_from_vertical)}. Needs \code{tube_thicc},
#'     \code{angle_from_vertical} and \code{center_offset}.}
#'   \item{\code{"window"}}{A planar rhizotron window. Depth increases along
#'     the image width and is scaled by \code{cos(angle_from_vertical)}, with
#'     no curvature term. \code{tube_thicc} and \code{center_offset} are
#'     unused.}
#'   \item{\code{"flat"}}{A surface with no depth axis at all -- a flatbed
#'     scan of washed roots on a tray. Every pixel is at the same depth,
#'     \code{-start_soil}, so a profile over this map collapses to one depth
#'     bin. \code{tube_thicc}, \code{angle_from_vertical} and
#'     \code{center_offset} are unused.}
#' }
#'
#' @param img Input image (accepts terra SpatRaster, matrix, array, or file
#'   path). Only its size is used, so any layer of a multi-layer image will do.
#'   A SpatRaster passes its extent, resolution and CRS on to the depth map.
#' @param mask Raster mask on the same grid as \code{img}, marking foreign
#'   objects such as tape with 1 (or TRUE); those cells are \code{NA} in the
#'   depth map. 0 and NA cells are kept.
#' @param sinoid Logical. Superseded by \code{geometry}, and kept so existing
#'   code runs unchanged: \code{TRUE} means \code{geometry = "tube"},
#'   \code{FALSE} means \code{geometry = "window"}. Ignored when
#'   \code{geometry} is given.
#' @param tube_thicc Numeric; diameter of minirhizotron tube in cm. Used by
#'   \code{geometry = "tube"} only
#' @param angle_from_vertical Numeric; insertion angle in degrees measured
#'   \strong{from vertical}, in \code{[0, 90]}. \code{0} is an upright tube and
#'   \code{90} a horizontal one, so the number is how far the tube leans over --
#'   the convention field notes usually record ("installed at 45 degrees").
#'   Used by \code{geometry = "tube"} and \code{"window"}.
#'
#'   The angle splits a tube's depth between its two axes, and the two terms are
#'   complements. The depth gained along the tube's length scales by
#'   \code{cos(angle_from_vertical)}: an upright tube gains a full centimeter of
#'   depth per centimeter of tube, a horizontal one gains none. The depth spread
#'   around its circumference scales by \code{sin(angle_from_vertical)}: an
#'   upright tube's cross-section is horizontal, so its whole circumference sits
#'   at one depth, while a horizontal tube's upper side is a full diameter above
#'   its lower side.
#' @param tilt Numeric. Superseded by \code{angle_from_vertical}, which states
#'   the convention in its name. Degrees from \emph{horizontal}, the opposite
#'   reference, so \code{tilt = t} is \code{angle_from_vertical = 90 - t}. Given
#'   it, the conversion is applied and a warning says so. \code{45} is the same
#'   angle either way
#' @param dpi Numeric; image resolution in dots per inch
#' @param start_soil Numeric; soil surface boundary in cm (0 = surface)
#' @param center_offset Numeric; the image row of the top of the tube -- the
#'   shallowest line of the curvature wave -- given in one of two forms:
#'   \itemize{
#'     \item \strong{Fraction} - a number in \code{[0, 1]}: the tube top is on
#'       row \code{1 + center_offset * (nrow - 1)}, so \code{0} puts it on the
#'       first row, \code{0.5} halfway down and \code{1} on the last row.
#'     \item \strong{Absolute row} - a number \code{> 1}: the row itself, at
#'       most \code{nrow}, e.g. from \code{estimate_rotation_center()}.
#'   }
#'   The wave itself
#'   always spans one full tube circumference (\code{tube_thicc * pi} cm at
#'   \code{dpi}); an image that covers only part of the circumference shows the
#'   matching part of the wave. Used by \code{geometry = "tube"} only
#' @param geometry Character; the shape being imaged. One of \code{"tube"}
#'   (default), \code{"window"} or \code{"flat"} -- see \strong{Geometries}.
#'
#' @return terra raster object containing the depth map
#' @export
#'
#' @author Johannes Cunow, Robert Weigel
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' seg_Oulanka2023_Session01_T067 = terra::rast(seg_Oulanka2023_Session01_T067)
#' img = seg_Oulanka2023_Session01_T067
#' mask = seg_Oulanka2023_Session01_T067[[1]] > 0   # layer 1 marks the tape
#' map = create_depthmap(img, mask, start_soil = 0.1,
#'   geometry = "tube",
#'   tube_thicc = 7,
#'   angle_from_vertical = 45,
#'   dpi = 300,
#'   center_offset = 0.1 )
#'
#' # A tray of washed roots: no depth axis, so every pixel gets one depth
#' flat = create_depthmap(img, geometry = "flat", dpi = 300)
#' terra::minmax(flat)
create_depthmap = function(img, mask = NULL, sinoid = TRUE,
                           tube_thicc = 7, tilt = NULL, dpi = 300,
                           start_soil = 0, center_offset = 0.5,
                           geometry = NULL, angle_from_vertical = 45) {

  # `missing()` reads the function's own frame, so both are resolved here rather
  # than inside the tryCatch() below, where the promise is evaluated.
  sinoid_given <- !missing(sinoid)
  angle_given  <- !missing(angle_from_vertical)

  tryCatch({
    # `geometry` names the shape; `sinoid` named only whether the curvature term
    # was added, so it cannot express the no-depth case. When both arrive,
    # `geometry` is the more specific statement and wins.
    if (is.null(geometry)) {
      if (!is.logical(sinoid) || length(sinoid) != 1L || is.na(sinoid))
        stop("sinoid must be TRUE or FALSE")
      geometry <- if (sinoid) "tube" else "window"
    } else {
      geometry <- match.arg(geometry, c("tube", "window", "flat"))
      if (sinoid_given) {
        implied <- if (isTRUE(sinoid)) "tube" else "window"
        if (!identical(implied, geometry))
          warning(sprintf(
            "Both `geometry = \"%s\"` and `sinoid = %s` were given; `geometry` is used.",
            geometry, sinoid), call. = FALSE)
      }
    }

    # `tilt` counted from horizontal, which is the opposite reference, so the
    # same number meant a different tube under each name. Converting here keeps
    # one angle convention past this point.
    if (!is.null(tilt)) {
      if (!is.numeric(tilt) || length(tilt) != 1L || is.na(tilt))
        stop("tilt must be a single number, in degrees from horizontal")
      if (angle_given && !isTRUE(all.equal(90 - tilt, angle_from_vertical)))
        stop("Both `tilt` and `angle_from_vertical` were given and they ",
             "describe different angles. Pass `angle_from_vertical` alone.",
             call. = FALSE)
      warning(sprintf(
        "`tilt` is superseded by `angle_from_vertical`, measured from vertical: use angle_from_vertical = %g.",
        90 - tilt), call. = FALSE)
      angle_from_vertical <- 90 - tilt
    }

    needs_tilt <- geometry %in% c("tube", "window")

    # Validate numeric parameters
    if (geometry == "tube" && (!is.numeric(tube_thicc) || tube_thicc <= 0))
      stop("tube_thicc must be a positive number")
    # The endpoints are both real cases: 0 is an upright tube, which gains a full
    # centimeter of depth per centimeter of tube and has its whole circumference
    # at one depth, and 90 is a horizontal tube, which is the reverse.
    if (needs_tilt && (!is.numeric(angle_from_vertical) ||
                       length(angle_from_vertical) != 1L ||
                       is.na(angle_from_vertical) ||
                       angle_from_vertical < 0 || angle_from_vertical > 90))
      stop("angle_from_vertical must be a single number between 0 and 90 degrees")
    if (!is.numeric(dpi) || dpi <= 0)
      stop("dpi must be a positive number")
    if (!is.numeric(start_soil))
      stop("start_soil must be a numeric value")
    if (geometry == "tube" &&
        (!is.numeric(center_offset) || length(center_offset) != 1L ||
         is.na(center_offset) || center_offset < 0))
      stop("center_offset must be a fraction in [0, 1] or a row number above 1")

    # Image validation
    if (is.null(img))
      stop("Input image cannot be NULL")

    # Only the size of the image is used, so any layer will do. A SpatRaster's
    # geometry is kept for the output.
    template <- if (inherits(img, "SpatRaster")) img[[1]] else NULL
    img <- load_flexible_image(img, select_layer = 1,
                               output_format = "spatrast", scale = "none")

    # Validate image dimensions
    if (any(dim(img)[1:2] <= 1))
      stop("Image dimensions must be greater than 1x1")

    # Mask validation and processing
    if (length(mask) == 0) {
      mask = img
      terra::values(mask) <- 0
    } else {
      # Check if mask dimensions match image
      if (!all(dim(mask)[1:2] == dim(img)[1:2]))
        stop("Mask dimensions must match image dimensions")
    }


    # Core computation with additional safety checks
    radiant = pi/180
    # The angle splits a tube's depth between its length and its circumference,
    # and the two factors are complements of each other. Along the length, an
    # upright tube (0) gains the full depth and a horizontal one (90) gains
    # none, which is cos. Across the circumference it is the reverse, which is
    # sin -- see `curve.factor` in the curvature branch below. "flat" has no
    # depth axis at all, so both factors are 0 whatever the angle says.
    tilt.factor  = if (needs_tilt) cos(angle_from_vertical * radiant) else 0
    curve.factor = if (needs_tilt) sin(angle_from_vertical * radiant) else 0
    px.to.cm.h = (2.54/dpi)

    # Image geometry. Depth increases along the image WIDTH (columns); the
    # sinusoidal tube-curvature offset varies along the HEIGHT (rows). The map
    # is built directly in `img` orientation so the result is aligned with
    # `img` (same extent/resolution/CRS) -- no transpose/flip is needed
    # downstream.
    nr = dim(img)[1]   # height -> tube circumference (sinusoid axis)
    nc = dim(img)[2]   # width  -> depth axis

    # Sinusoidal curvature offset, one value per circumference pixel (row).
    # Only a tube has one; a window and a tray are planar.
    if (geometry == "tube") {
      tube_thicc_tilted = round(tube_thicc * curve.factor, 3)
      circumference_px = tube_thicc * pi * dpi / 2.54
      if (circumference_px < 2)
        stop("Insufficient points for sine wave generation. Check tube_thicc and dpi values")
      # A scan covers at most one turn of the tube, so a taller image means the
      # diameter or the resolution is wrong.
      if (nr > round(circumference_px))
        stop("Image height (", nr, " px) exceeds one tube circumference (",
             round(circumference_px), " px at this tube_thicc and dpi). ",
             "Check tube_thicc and dpi.")

      # The term runs 0 (tube top, shallowest) to tube_thicc_tilted (tube
      # bottom); 1 - cos puts the trough at phase 0. center_offset places the
      # tube top on a row of the image, and the phase advances by one full turn
      # per circumference, so an image covering part of the tube carries the
      # matching stretch of the wave.
      # Up to 1 the value is a fraction of the image height, above 1 a row.
      top_row = if (center_offset <= 1) 1 + center_offset * (nr - 1) else center_offset
      if (top_row > nr)
        stop("center_offset = ", center_offset, " is beyond the last row (", nr, ")")
      phase   = 2 * pi * (seq_len(nr) - top_row) / circumference_px
      circ    = (1 - cos(phase)) * (tube_thicc_tilted / 2)
    } else {
      circ = rep(0, nr)
    }

    # Depth term per column (deeper to the right); curvature term per row.
    # M[row, col] = curvature(row) + depth(col) - start_soil, in img orientation.
    # Both terms are zero for "flat", which leaves every cell at -start_soil:
    # one depth for the whole image, so a profile over it has a single bin.
    depth_term = seq_len(nc) * px.to.cm.h * tilt.factor
    M = outer(circ, depth_term, FUN = "+") - start_soil   # dims [nr, nc]

    # Write onto a raster with the input's geometry, so extent/res/CRS match.
    masked.depthmap = terra::rast(if (is.null(template)) img[[1]] else template)
    terra::values(masked.depthmap) = as.vector(t(M))           # row-major for terra

    # Mask foreign objects. `mask` shares img's grid, so cell i of the mask and
    # cell i of the depth map are the same pixel.
    mv = as.vector(terra::values(mask))
    terra::values(masked.depthmap)[which(mv == 1)] <- NA

    return(masked.depthmap)

  }, error = function(e) {
    # Provide informative error message
    stop(paste("Error in create_depthmap:", e$message))
  })
}


#' Create a buffer halo) around non-zero pixels
#'
#' @param img SpatRaster/matrix/array - segmented image
#' @param width numeric - buffer width in pixels (default: 2)
#' @param halo.only logical - if TRUE, returns only the buffer zone (default: TRUE)
#' @param kernel character - shape of the thickening kernel: "circle" or "diamond"
#'
#' @return SpatRast - buffer zone around non-zero pixels
#' @export
#'
#'
#' @examples
#' data(seg_Oulanka2023_Session03_T067)
#' img <- terra::rast(seg_Oulanka2023_Session03_T067)
#' create_root_buffer(img, width = 2)
create_root_buffer = function(img, width=2, halo.only=TRUE, kernel="circle") {
  
  # Validation module
  tryCatch({
    # Check required input
    if (is.null(img)) {
      stop("Input image is required")
    }
    
    # Parameter validation
    if (!is.numeric(width) || width < 1) {
      stop("width must be a positive integer")
    }
    if (!is.logical(halo.only)) {
      stop("halo.only must be logical")
    }
    if (!kernel %in% c("circle", "diamond")) {
      stop("kernel must be either 'circle' or 'diamond'")
    }
    
    # Load and validate image
    im <- tryCatch({
      load_flexible_image(img, output_format="SpatRaster", scale = "binary")
    }, error = function(e) {
      stop("Failed to load image: ", e$message)
    })
    
    if (terra::ncell(im) == 0) {
      stop("Input image has no valid cells")
    }
    
    # Create kernel based on type
    k0 = if (kernel == "circle") {
      matrix(c(1,1,1,1,0,1,1,1,1), nrow=3, ncol=3)
    } else {
      matrix(c(0,1,0,1,0,1,0,1,0), nrow=3, ncol=3)
    }
    
    # Apply focal operation with error handling
    im2 = im
    for (itr in 1:width) {
      im2 <- tryCatch({
        terra::focal(im2, w=k0, fun="sum", na.policy="omit")
      }, error = function(e) {
        stop(sprintf("Focal operation failed at iteration %d: %s", itr, e$message))
      })
    }
    
    # Create output image
    out.im = (im2 >= 1) * 1
    
    if (halo.only) {
      out.im = out.im - im
      out.im = terra::subst(out.im, from=-1, to=0)
    }
    
    # Validate output
    if (terra::ncell(out.im) == 0) {
      warning("Output image has no valid cells")
    }
    
    return(out.im)
    
  }, error = function(e) {
    stop(paste("Error in Halo:", e$message))
  }, warning = function(w) {
    warning(paste("Warning in Halo:", w$message))
  })
}



#' Bin continuous depth values into discrete intervals
#'
#' @param depthmap SpatRaster/matrix/array - continuous depth values
#' @param nn numeric - bin width
#' @param round.option character - binning method: "rounding", "ceiling", or "floor"
#' @return SpatRaster - binned depth values
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' mask = img[[1]] - img[[2]]
#' mask[mask == 255] <- NA
#' img = img
#' depthmap = create_depthmap(img,mask,start.soil = 2.9 )
#' binned.map = binning(depthmap,nn = 5)
binning = function(depthmap, nn, round.option="rounding") {
  
  # Validation module
  tryCatch({
    # Check required inputs
    if (is.null(depthmap)) {
      stop("depthmap is required")
    }
    if (missing(nn)) {
      stop("bin width (nn) is required")
    }
    
    # Parameter validation
    if (!is.numeric(nn) || nn <= 0) {
      stop("bin width (nn) must be positive numeric")
    }
    if (!round.option %in% c("rounding", "ceiling", "floor")) {
      stop("round.option must be one of: 'rounding', 'ceiling', 'floor'")
    }
    
    # Load and validate image
    img <- tryCatch({
      load_flexible_image(depthmap, output_format="spatrast", scale = "none")
    }, error = function(e) {
      stop("Failed to load depthmap: ", e$message)
    })
    
    if (terra::ncell(img) == 0) {
      stop("Depthmap has no valid cells")
    }
    
    # Check for infinite or NA values
    if (any(is.infinite(terra::values(img)), na.rm=TRUE)) {
      warning("Infinite values detected in depthmap")
    }
    
    # Perform binning based on selected method
    im = tryCatch({
      if (round.option == "rounding") {
        nn * round(depthmap / nn, 0)
      } else if (round.option == "ceiling") {
        nn * ceiling(depthmap / nn)
      } else {
        nn * floor(depthmap / nn)
      }
    }, error = function(e) {
      stop("Binning operation failed: ", e$message)
    })
    
    # Validate output
    if (terra::ncell(im) == 0) {
      warning("Output has no valid cells")
    }
    
    return(im)
    
  }, error = function(e) {
    stop(paste("Error in binning:", e$message))
  }, warning = function(w) {
    warning(paste("Warning in binning:", w$message))
  })
}



#' Slice a scan along the rotation (circumferential) axis
#'
#' Splits a minirhizotron strip into \code{n} equal, contiguous, non-overlapping
#' slices along the rotation axis (image rows; depth runs along the columns).
#' Slice 1 is the top row band. Use to measure a trait per circumferential
#' position, e.g. as input to \code{rhythmicity()}.
#'
#' @param img A terra::SpatRaster.
#' @param n Number of slices around the circumference.
#' @return A list of \code{n} SpatRasters, in rotation order.
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' seg <- terra::rast(seg_Oulanka2023_Session01_T067)
#' slices <- slice_rotation(seg, 16)
slice_rotation <- function(img, n) {
  if (!inherits(img, "SpatRaster")) stop("img must be a terra::SpatRaster")
  if (!is.numeric(n) || length(n) != 1 || n < 1) stop("n must be a positive integer")
  nr <- terra::nrow(img)
  if (n > nr) stop("n (", n, ") cannot exceed the number of rotation pixels (", nr, ")")
  
  brk  <- round(seq(0, nr, length.out = n + 1))   # row breakpoints, exact coverage
  e    <- terra::ext(img)
  yres <- terra::res(img)[2]
  ymax <- e[4]
  
  lapply(seq_len(n), function(i) {
    r0 <- brk[i] + 1; r1 <- brk[i + 1]            # rows in this slice (1 = top)
    ytop <- ymax - (r0 - 1) * yres                # row -> y (row 1 sits at ymax)
    ybot <- ymax - r1 * yres
    terra::crop(img, terra::ext(e[1], e[2], ybot, ytop))
  })
}


#' Mask a scan to a depth and/or rotation zone (zone masking)
#'
#' \code{zoning()} performs \emph{zone masking}: it returns the input image with
#' every pixel \emph{outside} the requested zone set to \code{NA}, keeping the
#' original grid and extent. This lets you run any per-pixel trait function
#' (length, diameter, colour, landscape metrics) on a single depth slice without
#' splitting the raster into separate objects.
#'
#' The zone is defined by one or both of:
#' \itemize{
#'   \item \strong{depth} (\code{mode = "depth"}): keep only pixels whose binned
#'     depth (from \code{depth_map}) matches \code{depth}. A single value selects
#'     the closest available bin; a consecutive sequence selects an inclusive
#'     range; discrete values select each closest bin.
#'   \item \strong{rotation} (\code{mode = "rotation"}): crop to a contiguous band
#'     of circumferential slices. For simply splitting a tube into \code{n} equal
#'     slices, prefer \code{\link{slice_rotation}}; use the rotation mode here only
#'     when you need to combine a rotation band with depth masking in one call
#'     (\code{mode = "both"}).
#' }
#'
#' Note that depth masking sets out-of-zone pixels to \code{NA} (extent
#' unchanged), whereas rotation masking \emph{crops} the extent to the selected
#' band.
#'
#' @param img A \code{terra::SpatRaster}.
#' @param depth_map A \code{terra::SpatRaster} of binned depth values (see
#'   \code{\link{binning}}), aligned to \code{img}. Required for
#'   \code{mode = "depth"} or \code{"both"}.
#' @param depth Numeric. Target depth bin(s) to keep. Required for
#'   \code{mode = "depth"} or \code{"both"}.
#' @param select_layer Integer. Optionally return a single layer of the result.
#' @param crop_extent Numeric length-4 \code{c(xmin, xmax, ymin, ymax)}. Optional
#'   crop applied before masking.
#' @param mode One of \code{"rotation"}, \code{"depth"}, or \code{"both"}.
#' @param rotation_slices Numeric length-2 \code{c(from, to)} band of slices to
#'   keep (1 = top). Required for \code{mode = "rotation"} or \code{"both"}.
#' @param rotation_total_slices Numeric. Total number of slices the circumference
#'   is divided into. Required for \code{mode = "rotation"} or \code{"both"}.
#' @return A \code{terra::SpatRaster} masked (and/or cropped) to the requested
#'   zone.
#' @seealso \code{\link{binning}}, \code{\link{slice_rotation}}
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img <- terra::rast(seg_Oulanka2023_Session01_T067)
#' mask <- img[[1]] - img[[2]]; mask[mask == 255] <- NA
#' depth_map  <- create_depthmap(img, mask, start.soil = 2.9, dpi = 150 )
#' depth_bins <- binning(depth_map, nn = 5)
#' depth_bins <- terra::flip(terra::t(depth_bins))
#' # Keep only the root layer pixels that fall in the 10 cm depth bin
#' slice_10cm <- zoning(img[[2]], mode = "depth",
#'                      depth_map = depth_bins, depth = 10)
zoning <- function(
    img,
    depth_map = NULL,
    depth = NULL,
    select_layer = NULL,
    crop_extent = NULL,
    mode = c("rotation", "depth", "both"),
    rotation_slices = NULL,
    rotation_total_slices = NULL
) {
  mode <- match.arg(mode)

  if (!inherits(img, "SpatRaster")) {
    stop("img must be a terra::SpatRaster object.")
  }

  # Validate mode-specific requirements
  if (mode %in% c("depth", "both") && is.null(depth_map)) {
    stop("depth_map must be supplied for mode = 'depth' or 'both'.")
  }
  if (mode %in% c("depth", "both") && is.null(depth)) {
    stop("depth must be supplied for mode = 'depth' or 'both'.")
  }

  if (mode %in% c("rotation", "both")) {
    if (is.null(rotation_slices) || length(rotation_slices) != 2) {
      stop("rotation_slices must be numeric vector length 2 for mode 'rotation' or 'both'")
    }
    if (is.null(rotation_total_slices) || length(rotation_total_slices) != 1 || rotation_total_slices < 1) {
      stop("rotation_total_slices must be positive numeric scalar for mode 'rotation' or 'both'")
    }
    if (rotation_slices[1] >= rotation_slices[2]) {
      stop("rotation_slices must have increasing values: rotation_slices[1] < rotation_slices[2]")
    }
    if (any(rotation_slices > rotation_total_slices)) {
      stop("Values in rotation_slices cannot be greater than rotation_total_slices")
    }
  }

  # Validate spatial compatibility between img and depth_map
  if (!is.null(depth_map)) {
    if (!inherits(depth_map, "SpatRaster")) {
      stop("depth_map must be a terra::SpatRaster object.")
    }
    if (!terra::compareGeom(img, depth_map, stopOnError = FALSE)) {
      stop("img and depth_map must have compatible spatial properties (extent, resolution, CRS).")
    }
  }

  # Resolve requested depth value(s) to the closest available bin(s)
  depth_range_from_input <- function(depth, available) {
    available <- sort(unique(available[!is.na(available)]))
    if (length(available) == 0) {
      return(numeric(0))
    }
    if (length(depth) == 1) {
      matched <- available[which.min(abs(available - depth))]
    } else {
      sorted_depth <- sort(depth)
      diffs <- diff(sorted_depth)
      is_consecutive <- length(unique(round(diffs, 10))) == 1 && all(diffs > 0)
      if (is_consecutive && length(depth) > 2) {
        range_min <- min(depth)
        range_max <- max(depth)
        matched <- available[available >= range_min & available <= range_max]
      } else {
        matched <- sapply(depth, function(x) {
          available[which.min(abs(available - x))]
        })
        matched <- sort(unique(matched[!is.na(matched)]))
      }
    }
    return(matched)
  }

  original_depth_map <- depth_map

  # --- Optional crop first, for all modes ---
  if (!is.null(crop_extent)) {
    if (length(crop_extent) != 4) {
      stop("crop_extent must be c(xmin, xmax, ymin, ymax)")
    }
    if (crop_extent[1] >= crop_extent[2] || crop_extent[3] >= crop_extent[4]) {
      stop("Invalid crop_extent: xmin must be < xmax and ymin must be < ymax")
    }
    ext <- terra::ext(crop_extent[1], crop_extent[2], crop_extent[3], crop_extent[4])

    if (mode %in% c("depth", "both") && !is.null(depth_map)) {
      depth_vals_vec_original <- terra::values(original_depth_map, mat = FALSE)
      available_original <- sort(unique(depth_vals_vec_original[!is.na(depth_vals_vec_original)]))
      target_depths <- depth_range_from_input(depth, available_original)

      depth_map_cropped <- terra::crop(depth_map, ext)
      depth_vals_vec_cropped <- terra::values(depth_map_cropped, mat = FALSE)
      available_cropped <- sort(unique(depth_vals_vec_cropped[!is.na(depth_vals_vec_cropped)]))

      lost_depths <- setdiff(target_depths, available_cropped)
      if (length(lost_depths) > 0) {
        warning(sprintf("Cropping will remove pixels with desired depth values: %s. Consider adjusting crop_extent or depth parameters.",
                        paste(lost_depths, collapse = ", ")))
      }
      depth_map <- depth_map_cropped
    }
    img <- terra::crop(img, ext)
  }

  # --- Depth zone masking: NA outside the selected bin(s), extent unchanged ---
  if (mode %in% c("depth", "both")) {
    depth_vals_vec <- terra::values(depth_map, mat = FALSE)
    available <- sort(unique(depth_vals_vec[!is.na(depth_vals_vec)]))
    depth_vals <- depth_range_from_input(depth, available)

    if (length(depth_vals) == 0) {
      stop("No matching depths found in depth_map for the provided depth argument.")
    }

    if (length(depth_vals) == 1) {
      mask_condition <- depth_map == depth_vals[1]
    } else {
      mask_condition <- depth_map == depth_vals[1]
      for (i in 2:length(depth_vals)) {
        mask_condition <- mask_condition | (depth_map == depth_vals[i])
      }
    }
    img <- terra::ifel(mask_condition, img, NA)
  }

  # --- Rotation zone masking: crop to a contiguous band of slices ---
  if (mode %in% c("rotation", "both")) {
    img_dims <- dim(img)
    total_rows <- img_dims[1]

    slice_start <- floor((rotation_slices[1] - 1) * total_rows / rotation_total_slices) + 1
    slice_end <- ceiling(rotation_slices[2] * total_rows / rotation_total_slices)
    slice_start <- max(1, slice_start)
    slice_end <- min(total_rows, slice_end)

    x_res <- terra::res(img)[1]
    img_extent <- terra::ext(img)

    x_min_coord <- img_extent[1] + (slice_start - 1) * x_res
    x_max_coord <- img_extent[1] + slice_end * x_res

    rotation_crop_extent <- terra::ext(x_min_coord, x_max_coord,
                                       img_extent[3], img_extent[4])
    img <- terra::crop(img, rotation_crop_extent)
  }

  if (!is.null(select_layer)) {
    if (!is.numeric(select_layer) || length(select_layer) != 1 || select_layer < 1) {
      stop("select_layer must be a single positive integer.")
    }
    if (terra::nlyr(img) < select_layer) {
      stop("select_layer exceeds number of layers in img.")
    }
    img <- terra::subset(img, select_layer)
  }

  if (all(is.na(terra::values(img, mat = FALSE)))) {
    warning("Resulting image contains only NA values. Check your depth, rotation, or crop parameters.")
  }

  return(img)
}

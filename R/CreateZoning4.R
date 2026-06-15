
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
      load_flexible_image(img, output_format="SpatRaster", normalize=TRUE, binarize = TRUE)
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
      load_flexible_image(depthmap, output_format="spatrast", normalize=FALSE)
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

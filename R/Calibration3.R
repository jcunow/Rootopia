##############################
#' Estimates rotation from tape coverage
#'
#' This function analyzes image data to determine rotation based on tape coverage,
#' assuming more tape is present on the upper side of the tube.
#'
#' @param img Input image as raster, file name, or array
#' @param tape.brightness Brightness threshold for tape detection (0-1)
#' @param tape.quantile Quantile used to align brightness with tape (0-1)
#' @param extra.rows Additional rows to add for analysis
#' @param search.area Proportion of image to analyze (0-1)
#' @param nclasses Number of classes for pixel clustering
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `NULL`.
#' @return numeric Position of the center of extruding tape
#' @export
#'
#' @examples
#' img = seg_Oulanka2023_Session01_T067
#' r0 = estimate_rotation_center(img)
estimate_rotation_center = function(img, tape.brightness=0.66, extra.rows=100, search.area=0.45,
                     tape.quantile=0.98, nclasses=3, select.layer=NULL) {

  # Validation module
  tryCatch({
    # Check required input
    if (is.null(img)) {
      stop("Input image is required")
    }

    # Parameter validation
    if (!is.numeric(tape.brightness) || tape.brightness < 0 || tape.brightness > 1) {
      stop("tape.brightness must be numeric between 0 and 1")
    }
    if (!is.numeric(search.area) || search.area <= 0 || search.area > 1) {
      stop("search.area must be numeric between 0 and 1")
    }
    if (!is.numeric(tape.quantile) || tape.quantile <= 0 || tape.quantile > 1) {
      stop("tape.quantile must be numeric between 0 and 1")
    }
    if (!is.numeric(extra.rows) || extra.rows < 0) {
      stop("extra.rows must be a positive numeric value")
    }
    if (!is.numeric(nclasses) || nclasses < 2) {
      stop("nclasses must be numeric and at least 2")
    }
    if (!is.null(select.layer) && (!is.numeric(select.layer) || select.layer < 1)) {
      stop("select.layer must be NULL or a positive integer")
    }

    # Load and validate image
    im <- load_flexible_image(img, select.layer=select.layer,
                              output_format="array", normalize=TRUE)
    if (is.null(im)) {
      stop("Failed to load image")
    }
    if (length(dim(im)) != 3) {
      stop("Input image must be 3-dimensional array (RGB)")
    }

    # Main function logic with additional error handling
    red.line = array(dim=c(dim(im)[1], extra.rows, dim(im)[3]))
    red.line[,,1:dim(im)[3]] <- stats::quantile(im[,,1], tape.quantile, na.rm=TRUE)
    img1 = abind2(red.line, im, along=2)

    r.img1 = terra::rast(img1)
    r.img1 = terra::crop(r.img1, raster::extent(0, search.area*terra::ext(r.img1)[2],
                                                0, terra::ext(r.img1)[4]))

    # Clustering validation
    if (terra::ncell(r.img1) == 0) {
      stop("No valid pixels after cropping")
    }

    r1 = tryCatch({
      RStoolbox::unsuperClass(r.img1, nClasses=nclasses)
    }, error = function(e) {
      stop("Clustering failed: ", e$message)
    })

    clust.center = apply(r1$model$centers, 1, mean)
    threshold = tape.brightness * terra::global(r.img1, "max")[[1]]
    valid_clusters = clust.center[clust.center > threshold]

    if (length(valid_clusters) == 0) {
      warning("No clusters found above brightness threshold")
      return(NA)
    }

    clust = which(clust.center == max(valid_clusters))
    rr1 = r1$map == clust
    rr1 = rr1 * 1

    # Final calculations
    rsums = rowSums(terra::as.array(rr1), na.rm=TRUE)
    if (all(is.na(rsums))) {
      warning("No valid pixels found for rotation calculation")
      return(NA)
    }

    bin = dplyr::ntile(rsums, 2)
    zero.rotation.center = stats::median(which(bin == 2))

    return(zero.rotation.center)

  }, error = function(e) {
    stop(paste("Error in RotationE:", e$message))
  }, warning = function(w) {
    warning(paste("Warning in RotationE:", w$message))
  })
}

### Rotation Censoring

# Session = sampling campaign with scans taken in short succession
# the CI-600 doesn't have a full 360 degree rotation and thus sequential images have an region with no overlap
# these regions needs to be censored (next implementation)
# keep in mind that inner and outer tube diameter differ and a structural underestimation of rootlength needs a resize coefficient!!!

## input rotation matters for th conclusion!

#' Detect rotation shift between two images
#'
#' Calculates the rotation shift between two sequential images using either
#' cross-correlation or phase correlation methods.
#'
#' @param img1 Reference image (3-channel RGB)
#' @param img2 Subsequent image to compare (3-channel RGB)
#' @param cor.type Correlation type: "ccf" (cross) or "phase" (frequency domain)
#' @param fixed.depth.pixel Depth range to analyze c(start, end)
#' @param fixed.width Width of analysis region in pixels
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `NULL`.
#' @return Vector of shifts (x,y) in pixels
#' @export
#'
#' @examples
#' img1 = seg_Oulanka2023_Session01_T067
#' img2 = seg_Oulanka2023_Session03_T067
#' y.lag = estimate_rotation_shift(img1,img2,"phase")
estimate_rotation_shift = function(img1, img2, cor.type="phase",
                       fixed.depth.pixel=c(1000,4000),
                       fixed.width=NULL, select.layer=NULL) {

  # Validation module
  tryCatch({
    # Check required inputs
    if (is.null(img1) || is.null(img2)) {
      stop("Both input images are required")
    }

    # Parameter validation
    if (!cor.type %in% c("phase", "ccf")) {
      stop("cor.type must be either 'phase' or 'ccf'")
    }
    if (!is.numeric(fixed.depth.pixel) || length(fixed.depth.pixel) != 2) {
      stop("fixed.depth.pixel must be numeric vector of length 2")
    }
    if (fixed.depth.pixel[1] >= fixed.depth.pixel[2]) {
      stop("fixed.depth.pixel[1] must be less than fixed.depth.pixel[2]")
    }
    if (!is.null(fixed.width) && (!is.numeric(fixed.width) || fixed.width <= 0)) {
      stop("fixed.width must be NULL or positive numeric")
    }
    if (!is.null(select.layer) && (!is.numeric(select.layer) || select.layer < 1)) {
      stop("select.layer must be NULL or positive integer")
    }

    # Load and validate images
    im1 <- load_flexible_image(img1, select.layer=select.layer,
                               output_format="array", normalize=FALSE)
    im2 <- load_flexible_image(img2, select.layer=select.layer,
                               output_format="array", normalize=FALSE)

    if (is.null(im1) || is.null(im2)) {
      stop("Failed to load one or both images")
    }
    if (length(dim(im1)) != 3 || length(dim(im2)) != 3) {
      stop("Both inputs must be 3-dimensional arrays (RGB)")
    }

    # Set fixed width if NULL
    if (is.null(fixed.width)) {
      fixed.width = min(dim(img1)[1], dim(img2)[2])
    }

    # Convert to grayscale
    img11 <- im1[,,1]*0.21 + im1[,,2]*0.72 + im1[,,3]*0.07
    img22 <- im2[,,1]*0.21 + im2[,,2]*0.72 + im2[,,3]*0.07

    # Handle dimension differences
    dif.dim1 = nrow(img11) - nrow(img22)
    dif.dim2 = ncol(img11) - ncol(img22)
    if (dif.dim1 != 0 || dif.dim2 != 0) {
      warning(sprintf("Images differ in size by %d x %d pixels", dif.dim1, dif.dim2))

      # Fix rotation dimension
      rot.dim = round((nrow(img11)/2)-(fixed.width/2)) :
        round((nrow(img11)/2)+(fixed.width/2)-1)

      if (any(rot.dim < 1) || any(rot.dim > nrow(img11))) {
        stop("Invalid rotation dimensions after adjustment")
      }

      img11 = img11[rot.dim, ]
      img22 = img22[rot.dim, ]

      # Fix length dimension
      if (any(fixed.depth.pixel > ncol(img11)) ||
          any(fixed.depth.pixel > ncol(img22))) {
        stop("fixed.depth.pixel exceeds image dimensions")
      }

      img11 = img11[, fixed.depth.pixel]
      img22 = img22[, fixed.depth.pixel]
    }

    img11 = terra::as.array(terra::rast(img11))
    img22 = terra::as.array(terra::rast(img22))

    # Correlation analysis
    a = if (cor.type == "phase") {
      imagefx::pcorr3d(img11[,,1], img22[,,1])
    } else {
      imagefx::xcorr3d(img11[,,1], img22[,,1])
    }

    y.lag = c(a$max.shifts[1])
    x.lag = c(a$max.shifts[2])

    if (abs(x.lag) > 10) {
      warning(sprintf("Large depth shift detected: %d pixels", x.lag))
    }

    return(c(x.lag, y.lag))

  }, error = function(e) {
    stop(paste("Error in estimate_rotation_shift:", e$message))
  }, warning = function(w) {
    warning(paste("Warning in estimate_rotation_shift:", w$message))
  })
}


### Cuts images along rotational center
#' Censor image edges based on rotation
#'
#' Crops image edges to handle non-overlapping regions between sequential scans.
#'
#' @param img Input image to censor
#' @param center.offset Rotation shift in rows (from estimate_rotation_shift())
#' @param cut.buffer Proportion of image to cut when fixed_rotation=FALSE
#' @param fixed.rotation Use fixed output dimensions
#' @param fixed.width Output width when fixed_rotation=TRUE
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `NULL`.
#' @return Cropped raster image
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' censored.raster = rotation_censor(img,
#'                          center.offset = 120,
#'                          cut.buffer = 0.02,
#'                          fixed.rotation = FALSE)
#'                          
#' censored.raster = rotation_censor(img,
#'                          center.offset = 220,
#'                          cut.buffer = 0.02,
#'                          fixed.width = 1000,
#'                          fixed.rotation = TRUE)
rotation_censor = function(img, center.offset=0, cut.buffer=0.02,
                     fixed.rotation=TRUE, fixed.width=500, select.layer=NULL) {

  # Validation module
  tryCatch({
    # Check required input
    if (is.null(img)) {
      stop("Input image is required")
    }

    # Parameter validation
    if (!is.numeric(center.offset)) {
      stop("center.offset must be numeric")
    }
    if (!is.numeric(cut.buffer) || cut.buffer < 0 || cut.buffer > 1) {
      stop("cut.buffer must be numeric between 0 and 1")
    }
    if (!is.logical(fixed.rotation)) {
      stop("fixed.rotation must be logical")
    }
    if (!is.numeric(fixed.width) || fixed.width <= 0) {
      stop("fixed.width must be positive numeric")
    }
    if (!is.null(select.layer) && (!is.numeric(select.layer) || select.layer < 1)) {
      stop("select.layer must be NULL or positive integer")
    }

    # Load and validate image
    img.c <- load_flexible_image(img, select.layer=select.layer,
                                 output_format="spatrast", normalize=FALSE,binarize = FALSE)
    if (is.null(img.c)) {
      stop("Failed to load image")
    }

    # Calculate ratios
    offset.ratio = abs(center.offset) / dim(img.c)[1]
    cut.ratio = offset.ratio + cut.buffer

    # Validate ratios
    if (cut.ratio >= 1) {
      stop("Cut ratio too large, would remove entire image")
    }

    
    if (!fixed.rotation) {
      # Handle variable rotation cases
      ex = terra::ext(img.c)

      if (offset.ratio >= 0.5) {
        ex[4] = (1-cut.ratio) * ex[4]
      } else if (offset.ratio < 0.5) {
        ex[3] = cut.ratio * ex[3]
      } else {
        ex[3] = (cut.ratio/2) * ex[3]
        ex[4] = (1-cut.ratio/2) * ex[4]
      }

      img.cc = terra::crop(img.c, ex)
      
      if (terra::ncell(img.cc) == 0) {
        warning("Resulting image is empty after censoring")
        return(NULL)
      }
      return(img.cc)

      
    } else {
      # Handle fixed rotation
      center.row = center.offset
      ex = terra::ext(img.c)
      ex[2] = dim(img.c)[2]
      ex[4] = dim(img.c)[1]
      terra::ext(img.c) <- ex

      new_mid_point = center.row #+ (dim(img.c)[1]/2) 
      ex[3] = new_mid_point - fixed.width/2
      ex[4] = new_mid_point + fixed.width/2


      img.cc = terra::crop(img.c, ex)
      
      # Validate dimensions
      if (ex[3] < 0 || ex[4] > dim(img.c)[1]) {
        message("New image dimension: ",dim(img.cc)[1] ," is smaller than specified fixed.width: ", fixed.width, ". Too strong offset for this fixed.width. Consider adjusting the fixed.width.")
      }
      
      if (terra::ncell(img.cc) == 0) {
        warning("Resulting image is empty after censoring")
        return(NULL)
      }
      return(img.cc)
    }


  

  }, error = function(e) {
    stop(paste("Error in RotCensor:", e$message))
  })
}



#' Estimate soil surface position using tape markers
#'
#' Detects the soil surface by analyzing tape coverage patterns in the image.
#'
#' @param img Input image (raster, filename, or array)
#' @param search.area Proportion of image to analyze
#' @param tape.tresh Minimum tape coverage ratio
#' @param dpi Image resolution
#' @param nclasses Number of clustering classes
#' @param inverse Invert detection for dark markers
#' @param tape.overlap Safety margin for tape (cm)
#' @param tape.brightness Brightness threshold for tape
#' @param extra.rows Additional analysis rows
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `NULL`.
#' @param tape.quantile Brightness alignment quantile
#' @return data.frame with soil surface and tape end positions
#' @export
#'
#' @examples
#' img = rgb_Oulanka2023_Session03_T067
#' Soil0Estimates = estimate_soil_surface(img)
estimate_soil_surface = function(img, search.area=0.45, tape.tresh=0.33, dpi=150, nclasses=3,
                     inverse=FALSE, tape.overlap=0.5, tape.brightness=0.6,
                     extra.rows=100, tape.quantile=0.98, select.layer=NULL) {

  # Validation module
  tryCatch({
    # Check required input
    if (is.null(img)) {
      stop("Input image is required")
    }

    # Parameter validation
    if (!is.numeric(search.area) || search.area <= 0 || search.area > 1) {
      stop("search.area must be numeric between 0 and 1")
    }
    if (!is.numeric(tape.tresh) || tape.tresh <= 0 || tape.tresh > 1) {
      stop("tape.tresh must be numeric between 0 and 1")
    }
    if (!is.numeric(dpi) || dpi <= 0) {
      stop("dpi must be positive numeric")
    }
    if (!is.numeric(nclasses) || nclasses < 2) {
      stop("nclasses must be numeric and at least 2")
    }
    if (!is.logical(inverse)) {
      stop("inverse must be logical")
    }
    if (!is.numeric(tape.overlap) || tape.overlap < 0) {
      stop("tape.overlap must be non-negative numeric")
    }
    if (!is.numeric(tape.brightness) || tape.brightness <= 0) {
      stop("tape.brightness must be positive numeric")
    }
    if (!is.numeric(extra.rows) || extra.rows < 0) {
      stop("extra.rows must be non-negative numeric")
    }
    if (!is.numeric(tape.quantile) || tape.quantile <= 0 || tape.quantile > 1) {
      stop("tape.quantile must be numeric between 0 and 1")
    }
    if (!is.null(select.layer) && (!is.numeric(select.layer) || select.layer < 1)) {
      stop("select.layer must be NULL or positive integer")
    }

    # Load and validate image
    im <- load_flexible_image(img, select.layer=select.layer,
                              output_format="array", normalize=FALSE)
    if (is.null(im)) {
      stop("Failed to load image")
    }
    if (length(dim(im)) != 3) {
      stop("Input image must be 3-dimensional array (RGB)")
    }

    # Adjust parameters for inverse mode
    if (inverse) {
      tape.quantile = 1 - tape.quantile
      tape.brightness = 1 / tape.brightness
    }

    # Create red line array with validation
    red.line = array(dim=c(dim(im)[1], extra.rows, dim(im)[3]))
    quantile_value = tryCatch({
      stats::quantile(im[,,1], tape.quantile, na.rm=TRUE)
    }, error = function(e) {
      stop("Failed to calculate quantile: ", e$message)
    })

    if (is.na(quantile_value) || !is.finite(quantile_value)) {
      stop("Invalid quantile value calculated")
    }

    red.line[,,1:dim(im)[3]] <- quantile_value
    img1 = tryCatch({
      abind2(red.line, im, along=2)
    }, error = function(e) {
      stop("Failed to combine arrays: ", e$message)
    })

    # Convert to raster and crop
    r.img1 = terra::rast(img1)
    r.img1 = terra::crop(r.img1, c(0, search.area*terra::ext(r.img1)[2],
                                   0, terra::ext(r.img1)[4]))

    if (terra::ncell(r.img1) == 0) {
      stop("No valid pixels after cropping")
    }

    # Clustering
    r1 = tryCatch({
      RStoolbox::unsuperClass(r.img1, nClasses=nclasses)
    }, error = function(e) {
      stop("Clustering failed: ", e$message)
    })

    # Calculate cluster centers
    clust.center = apply(r1$model$centers, 1, mean)

    # Find appropriate cluster
    if (inverse) {
      global_min = terra::global(r.img1, "min")[[1]][1:3]
      threshold = tape.brightness * global_min
      valid_clusters = clust.center[clust.center < threshold]
      if (length(valid_clusters) == 0) {
        warning("No clusters found below brightness threshold")
        return(data.frame(soil0=NA, tape.end=NA))
      }
      clust = which(clust.center == min(valid_clusters))
    } else {
      global_max = terra::global(r.img1, "max")[[1]][1:3]
      threshold = tape.brightness * global_max
      valid_clusters = clust.center[clust.center > threshold]
      if (length(valid_clusters) == 0) {
        warning("No clusters found above brightness threshold")
        return(data.frame(soil0=NA, tape.end=NA))
      }
      clust = which(clust.center == max(valid_clusters))
    }

    # Create binary map
    rr1 = r1$map == clust
    rr1 = rr1 * 1

    # Iterate over depth with safety checks
    i = 1
    max_iter = dim(rr1)[2]  # Prevent infinite loop
    iter_count = 0

    while (i <= max_iter && iter_count < max_iter) {
      current_ratio = sum(rr1[,i,], na.rm=TRUE) / ncol(rr1)
      ratio_24 = if (i+24 <= max_iter) sum(rr1[,i+24,], na.rm=TRUE) / ncol(rr1) else 0
      ratio_48 = if (i+48 <= max_iter) sum(rr1[,i+48,], na.rm=TRUE) / ncol(rr1) else 0

      if (current_ratio <= tape.tresh &&
          ratio_24 <= tape.tresh &&
          ratio_48 <= tape.tresh) {
        break
      }

      i = i + 1
      iter_count = iter_count + 1

      if (i > dim(rr1)[2]) {
        warning("Loop exceeded image limits, resetting to 1")
        i = 1
      }
    }

    # Calculate final positions
    rw.ind = i
    rw.ind = rw.ind - extra.rows - round(tape.overlap * dpi/2.54)
    rw.ind = round(rw.ind)
    tape.end = rw.ind + round(tape.overlap * dpi/2.54)
    tape.end = round(tape.end)

    # Validate final results
    if (rw.ind < 0 || tape.end < 0) {
      warning("Negative position values calculated")
    }

    out = data.frame(soil0=rw.ind, tape.end=tape.end)
    return(out)

  }, error = function(e) {
    stop(paste("Error in SoilSurfE:", e$message))
  }, warning = function(w) {
    warning(paste("Warning in SoilSurfE:", w$message))
  })
}






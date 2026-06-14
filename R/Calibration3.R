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
  tryCatch({
    if (is.null(img)) stop("Input image is required")
    if (!is.numeric(tape.brightness) || tape.brightness < 0 || tape.brightness > 1)
      stop("tape.brightness must be numeric between 0 and 1")
    if (!is.numeric(search.area) || search.area <= 0 || search.area > 1)
      stop("search.area must be numeric between 0 and 1")
    if (!is.numeric(tape.quantile) || tape.quantile <= 0 || tape.quantile > 1)
      stop("tape.quantile must be numeric between 0 and 1")
    if (!is.numeric(extra.rows) || extra.rows < 0) stop("extra.rows must be a positive numeric value")
    if (!is.numeric(nclasses) || nclasses < 2) stop("nclasses must be numeric and at least 2")
    if (!is.null(select.layer) && (!is.numeric(select.layer) || select.layer < 1))
      stop("select.layer must be NULL or a positive integer")
    
    im <- load_flexible_image(img, select.layer=select.layer,
                              output_format="array", normalize=TRUE)
    if (is.null(im)) stop("Failed to load image")
    if (length(dim(im)) != 3) stop("Input image must be 3-dimensional array (RGB)")
    
    # bright reference band, then crop to the shallow-depth search area
    red.line = array(dim=c(dim(im)[1], extra.rows, dim(im)[3]))
    red.line[,,1:dim(im)[3]] <- stats::quantile(im[,,1], tape.quantile, na.rm=TRUE)
    img1   = abind2(red.line, im, along=2)
    r.img1 = terra::rast(img1)
    r.img1 = terra::crop(r.img1, terra::ext(0, search.area*terra::ext(r.img1)[2],
                                            0, terra::ext(r.img1)[4]))
    if (terra::ncell(r.img1) == 0) stop("No valid pixels after cropping")
    
    vals      = terra::as.array(r.img1)                 # [rows, cols, layers]
    maxval    = terra::global(r.img1, "max", na.rm=TRUE)[[1]]
    threshold = tape.brightness * maxval
    
    # how many distinct colours are actually present? (cap k-means accordingly)
    px = matrix(vals, ncol = dim(vals)[3])
    px = px[stats::complete.cases(px), , drop = FALSE]
    if (nrow(px) > 2e5) px = px[sample(nrow(px), 2e5), , drop = FALSE]
    n.distinct = nrow(unique(round(px, 4)))
    
    mask.mat = NULL
    
    # try clustering only if the data can support >= 2 clusters; never abort on failure
    if (n.distinct >= 2) {
      k = max(2L, min(as.integer(nclasses), n.distinct - 1L))
      r1 = tryCatch(RStoolbox::unsuperClass(r.img1, nClasses = k),
                    error = function(e) NULL)
      if (!is.null(r1)) {
        clust.center = apply(r1$model$centers, 1, mean)
        valid        = which(clust.center > threshold)
        if (length(valid)) {
          clust    = valid[which.max(clust.center[valid])]
          mask.mat = terra::as.array(r1$map == clust)[, , 1] * 1
        }
      }
    }
    
    # fallback: direct brightness threshold (robust for segmented / low-entropy images)
    if (is.null(mask.mat)) {
      bright   = apply(vals, c(1, 2), function(z) mean(z, na.rm = TRUE))
      mask.mat = (bright > threshold) * 1
    }
    
    if (sum(mask.mat, na.rm = TRUE) == 0) {
      warning("No tape pixels detected above threshold; returning NA")
      return(NA)
    }
    
    rsums = rowSums(mask.mat, na.rm = TRUE)
    if (all(rsums == 0)) { warning("No valid pixels for rotation calculation"); return(NA) }
    
    bin = dplyr::ntile(rsums, 2)
    return(round(stats::median(which(bin == 2)),0))
    
  }, error = function(e) stop(paste("Error in estimate_rotation_center:", e$message)))
}




### Rotation Censoring
#' Estimate rotational/depth shift between two root scans
#'
#' @param img1,img2 Image inputs (paths, arrays, or rasters).
#' @param cor.type "phase" (phase correlation) or "ccf" (normalized cross-corr).
#' @param fixed.depth.pixel Depth band along COLUMNS. Length-2 = range start:end;
#'   longer = explicit column indices; NULL = use full width.
#' @param fixed.width Optional: restrict the ROTATION axis (rows), centered.
#' @param select.layer Layer to use for multi-band inputs.
#' @param window Demean + Hann-window before FFT to suppress edge artefacts.
#' @param overlay If TRUE, also draw a before/after magenta-green overlay.
#' @param overlay.layer Layer to display in the overlay (root mask, default 2).
#' @return Named numeric vector: depth (column lag), rotation (row lag), peak.
#' @export
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' data(seg_Oulanka2023_Session03_T067)
#' img1 <- terra::rast(seg_Oulanka2023_Session01_T067)
#' img2 <- terra::rast(seg_Oulanka2023_Session03_T067)
#' estimate_rotation_shift(img1, img2, cor.type = "phase", select.layer = 2)
estimate_rotation_shift <- function(
    img1, img2,
    cor.type = "phase",
    fixed.depth.pixel = NULL,
    fixed.width = NULL,
    select.layer = NULL,
    window = TRUE,
    overlay = FALSE,
    overlay.layer = 2
) {
  tryCatch({
    
    if (is.null(img1) || is.null(img2)) stop("Both input images are required")
    if (!cor.type %in% c("phase", "ccf")) stop("cor.type must be 'phase' or 'ccf'")
    
    im1 <- load_flexible_image(img1, select.layer = select.layer,
                               output_format = "array", normalize = FALSE)
    im2 <- load_flexible_image(img2, select.layer = select.layer,
                               output_format = "array", normalize = FALSE)
    if (length(dim(im1)) != 3 || length(dim(im2)) != 3)
      stop("Inputs must be 3D RGB arrays")
    
    # luma projection (Rec. 709). rows = rotation axis, cols = depth axis.
    g1 <- im1[, , 1] * 0.21 + im1[, , 2] * 0.72 + im1[, , 3] * 0.07
    g2 <- im2[, , 1] * 0.21 + im2[, , 2] * 0.72 + im2[, , 3] * 0.07
    
    # common extent on BOTH axes (pcorr3d/xcorr3d need identical dims)
    rr <- seq_len(min(nrow(g1), nrow(g2)))
    cc <- seq_len(min(ncol(g1), ncol(g2)))
    if (nrow(g1) != nrow(g2) || ncol(g1) != ncol(g2))
      warning("Image size mismatch detected; cropping to common extent")
    g1 <- g1[rr, cc, drop = FALSE]; g2 <- g2[rr, cc, drop = FALSE]
    
    # depth BAND on columns. length-2 is a RANGE, not two literal indices.
    if (!is.null(fixed.depth.pixel)) {
      d <- if (length(fixed.depth.pixel) == 2) {
        if (fixed.depth.pixel[1] >= fixed.depth.pixel[2])
          stop("fixed.depth.pixel[1] must be < fixed.depth.pixel[2]")
        seq(fixed.depth.pixel[1], fixed.depth.pixel[2])
      } else fixed.depth.pixel
      d <- d[d >= 1 & d <= ncol(g1)]
      if (length(d) < 8) stop("depth band too narrow after clipping to image")
      g1 <- g1[, d, drop = FALSE]; g2 <- g2[, d, drop = FALSE]
    }
    
    # optional rotation-axis (rows) restriction, centered
    if (!is.null(fixed.width)) {
      w  <- min(fixed.width, nrow(g1))
      st <- floor((nrow(g1) - w) / 2) + 1
      rsub <- st:(st + w - 1)
      g1 <- g1[rsub, , drop = FALSE]; g2 <- g2[rsub, , drop = FALSE]
    }
    
    # demean + separable Hann window -> stabilises the correlation peak
    if (isTRUE(window)) {
      hann <- function(n) if (n == 1) 1 else 0.5 - 0.5 * cos(2 * pi * (0:(n - 1)) / (n - 1))
      win  <- outer(hann(nrow(g1)), hann(ncol(g1)))
      g1 <- (g1 - mean(g1)) * win
      g2 <- (g2 - mean(g2)) * win
    }
    
    a <- if (cor.type == "phase") imagefx::pcorr3d(g1, g2) else imagefx::xcorr3d(g1, g2)
    
    # pcorr3d/xcorr3d: max.shifts[1] = ROW lag (rotation), [2] = COLUMN lag (depth)
    out <- c(depth    = as.numeric(a$max.shifts[2]),
             rotation = as.numeric(a$max.shifts[1]),
             peak     = as.numeric(a$max.cor))
    
    if (isTRUE(overlay)) {
      f1 <- im1[, , overlay.layer]; f2 <- im2[, , overlay.layer]   # roots, 0 = bg
      if (max(c(f1, f2), na.rm = TRUE) > 1) { f1 <- f1 / 255; f2 <- f2 / 255 }
      f1[is.na(f1)] <- 0; f2[is.na(f2)] <- 0
      
      dxr <- as.integer(round(out["depth"]))      # column (depth) shift
      dyr <- as.integer(round(out["rotation"]))   # row (rotation) shift
      
      shift_mat <- function(m, dx, dy) {
        nr <- nrow(m); nc <- ncol(m); o <- matrix(0, nr, nc)
        rf <- max(1, 1 + dy):min(nr, nr + dy); rt <- max(1, 1 - dy):min(nr, nr - dy)
        cf <- max(1, 1 + dx):min(nc, nc + dx); ct <- max(1, 1 - dx):min(nc, nc - dx)
        o[rt, ct] <- m[rf, cf]; o
      }
      place <- function(m, nr, nc) {
        o <- matrix(0, nr, nc); o[seq_len(nrow(m)), seq_len(ncol(m))] <- m; o
      }
      comp <- function(x, y) {
        nr <- max(nrow(x), nrow(y)); nc <- max(ncol(x), ncol(y))
        z <- array(0, c(nr, nc, 3))
        z[, , 1] <- place(x, nr, nc); z[, , 2] <- place(y, nr, nc); z[, , 3] <- place(x, nr, nc)
        grDevices::as.raster(z, max = 1)
      }
      draw <- function(r, main) {
        d <- dim(r)
        graphics::plot(NA, xlim = c(1, d[2]), ylim = c(d[1], 1), asp = 1,
                       axes = FALSE, xlab = "", ylab = "", main = main)
        graphics::rasterImage(r, 1, 1, d[2], d[1])
      }
      op <- graphics::par(mfrow = c(2, 1)); on.exit(graphics::par(op), add = TRUE)
      draw(comp(f1, f2),                      "Raw (magenta = img1, green = img2)")
      draw(comp(f1, shift_mat(f2, dxr, dyr)), sprintf("Aligned (depth=%d, rotation=%d)", dxr, dyr))
    }
    
    
    if (isTRUE(overlay)) invisible(out) else out  
    
  }, error = function(e) stop("estimate_rotation_shift failed: ", e$message, call. = FALSE))
}



#' Censor image edges based on rotation
#'
#' Crops the rotation axis (rows) of a root scan. In fixed mode it returns a
#' fixed-width window centred on a given row (e.g. the rotation centre from
#' \code{estimate_rotation_center()}); in variable mode it trims a band sized by
#' a measured offset. Optionally previews what is kept versus cut.
#'
#' @param img Input image as raster, file name, or array.
#' @param center.offset Numeric. Meaning depends on \code{fixed.rotation}:
#'   when \code{TRUE}, the row to centre the kept window on (an absolute row,
#'   e.g. from \code{estimate_rotation_center()}); when \code{FALSE}, the
#'   rotation shift in rows to trim (e.g. from \code{estimate_rotation_shift()}).
#' @param cut.buffer Extra proportion of the rotation axis to trim (variable mode).
#' @param fixed.rotation Logical. If \code{TRUE}, return a fixed-width window.
#' @param fixed.width Output width in rows when \code{fixed.rotation = TRUE}.
#' @param select.layer Integer or \code{NULL}. Layer to use for multi-band inputs.
#' @param overlay Logical. If \code{TRUE}, plot the full image with the kept
#'   window (green outline) and discarded margins (red shading) before cropping.
#'   Default \code{FALSE}.
#' @param ... Passed to the underlying \code{terra} plotting call when
#'   \code{overlay = TRUE}.
#' @return A cropped \code{SpatRaster} (returned invisibly when
#'   \code{overlay = TRUE}), or \code{NULL} if the result is empty.
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img <- terra::rast(seg_Oulanka2023_Session01_T067)
#' r0  <- estimate_rotation_center(img)
#' rotation_censor(img, center.offset = r0, fixed.width = 800,
#'                 fixed.rotation = TRUE, overlay = TRUE)
rotation_censor <- function(img, center.offset = 0, cut.buffer = 0.02,
                            fixed.rotation = TRUE, fixed.width = 500,
                            select.layer = NULL, overlay = FALSE, ...) {
  tryCatch({
    
    if (is.null(img)) stop("Input image is required")
    if (!is.numeric(center.offset)) stop("center.offset must be numeric")
    if (!is.numeric(cut.buffer) || cut.buffer < 0 || cut.buffer > 1)
      stop("cut.buffer must be numeric between 0 and 1")
    if (!is.logical(fixed.rotation)) stop("fixed.rotation must be logical")
    if (!is.numeric(fixed.width) || fixed.width <= 0)
      stop("fixed.width must be positive numeric")
    if (!is.null(select.layer) && (!is.numeric(select.layer) || select.layer < 1))
      stop("select.layer must be NULL or positive integer")
    
    img.c <- load_flexible_image(img, select.layer = select.layer,
                                 output_format = "spatrast",
                                 normalize = FALSE, binarize = FALSE)
    if (is.null(img.c)) stop("Failed to load image")
    
    nr <- dim(img.c)[1]; nc <- dim(img.c)[2]    # rows = rotation axis
    offset    <- round(center.offset)
    buffer.px <- round(cut.buffer * nr)
    
    # --- determine the kept row window (lo:hi, 1 = top) ---
    if (!fixed.rotation) {
      cut.px <- abs(offset) + buffer.px
      if (cut.px >= nr) stop("Cut size too large, would remove entire image")
      if (offset > 0)      { lo <- cut.px + 1; hi <- nr }
      else if (offset < 0) { lo <- 1;          hi <- nr - cut.px }
      else { h <- floor(buffer.px / 2); lo <- 1 + h; hi <- nr - h }
    } else {
      mid      <- offset
      half     <- fixed.width / 2
      max.half <- min(mid - 1, nr - mid)
      if (half > max.half)
        message("fixed.width = ", fixed.width, " cannot be centred symmetrically on row ",
                round(mid), " (image is ", nr, " rows). Max symmetric width here is ",
                max(0, floor(2 * max.half)), " px; clamping to image bounds.")
      lo <- max(1,  round(mid - half))
      hi <- min(nr, round(mid + half) - 1)
    }
    
    # --- optional preview of kept (green) vs cut (red) on the full image ---
    if (isTRUE(overlay)) {
      disp <- img.c
      terra::ext(disp) <- c(0, nc, 0, nr)        # pixel-unit coords; row 1 = ymax
      y_keep_bottom <- nr - hi
      y_keep_top    <- nr - (lo - 1)
      if (terra::nlyr(disp) >= 3) {
        mx <- terra::global(disp, "max", na.rm = TRUE)[[1]]
        terra::plotRGB(disp, r = 1, g = 2, b = 3, scale = max(1, mx), ...)
      } else {
        terra::plot(disp, ...)
      }
      red <- grDevices::adjustcolor("red", 0.35)
      graphics::rect(0, y_keep_top,    nc, nr,            col = red, border = NA)      # cut above
      graphics::rect(0, 0,             nc, y_keep_bottom, col = red, border = NA)      # cut below
      graphics::rect(0, y_keep_bottom, nc, y_keep_top,    border = "green3", lwd = 2)  # kept
      graphics::legend("topright", bty = "n",
                       legend = c("kept (censored output)", "cut"),
                       col = c("green3", red), pch = c(0, 15))
    }
    
    # --- apply the crop (exact pixel rows) ---
    rows <- lo:hi
    if (length(rows) == 0) { warning("Resulting image is empty after censoring"); return(NULL) }
    arr <- terra::as.array(img.c)[rows, , , drop = FALSE]
    
    out <- terra::rast(arr)
    terra::ext(out) <- c(0, ncol(out), 0, nrow(out))
    
    if (isTRUE(overlay)) invisible(out) else out
    
  }, error = function(e) stop(paste("Error in rotation_censor:", e$message)))
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






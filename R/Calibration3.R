##############################
#' Estimates rotation from tape coverage
#'
#' This function analyzes image data to determine rotation based on tape coverage,
#' assuming more tape is present on the upper side of the tube.
#'
#' @param img Input image as raster, file name, or array
#' @param tape_brightness Brightness threshold for tape detection (0-1)
#' @param tape_quantile Quantile used to align brightness with tape (0-1)
#' @param extra_rows Additional rows to add for analysis
#' @param search_area Proportion of image to analyze (0-1)
#' @param nclasses Number of classes for pixel clustering
#' @param select_layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `NULL`.
#' @return numeric Position of the center of extruding tape
#' @export
#'
#' @examples
#' img = seg_Oulanka2023_Session01_T067
#' r0 = estimate_rotation_center(img)
estimate_rotation_center = function(img, tape_brightness=0.66, extra_rows=100, search_area=0.45,
                                    tape_quantile=0.98, nclasses=3, select_layer=NULL) {
  if (!requireNamespace("RStoolbox", quietly = TRUE))
    stop("Package 'RStoolbox' is required for estimate_rotation_center(). ",
         "Install it with: install.packages(\"RStoolbox\")")
  tryCatch({
    if (is.null(img)) stop("Input image is required")
    if (!is.numeric(tape_brightness) || tape_brightness < 0 || tape_brightness > 1)
      stop("tape_brightness must be numeric between 0 and 1")
    if (!is.numeric(search_area) || search_area <= 0 || search_area > 1)
      stop("search_area must be numeric between 0 and 1")
    if (!is.numeric(tape_quantile) || tape_quantile <= 0 || tape_quantile > 1)
      stop("tape_quantile must be numeric between 0 and 1")
    if (!is.numeric(extra_rows) || extra_rows < 0) stop("extra_rows must be a positive numeric value")
    if (!is.numeric(nclasses) || nclasses < 2) stop("nclasses must be numeric and at least 2")
    if (!is.null(select_layer) && (!is.numeric(select_layer) || select_layer < 1))
      stop("select_layer must be NULL or a positive integer")
    
    im <- load_flexible_image(img, select_layer=select_layer,
                              output_format="array", scale = "to_01")
    if (is.null(im)) stop("Failed to load image")
    if (length(dim(im)) != 3) stop("Input image must be 3-dimensional array (RGB)")
    
    # bright reference band, then crop to the shallow-depth search area
    red.line = array(dim=c(dim(im)[1], extra_rows, dim(im)[3]))
    red.line[,,1:dim(im)[3]] <- stats::quantile(im[,,1], tape_quantile, na.rm=TRUE)
    img1   = abind2(red.line, im, along=2)
    r.img1 = terra::rast(img1)
    r.img1 = terra::crop(r.img1, terra::ext(0, search_area*terra::ext(r.img1)[2],
                                            0, terra::ext(r.img1)[4]))
    if (terra::ncell(r.img1) == 0) stop("No valid pixels after cropping")
    
    vals      = terra::as.array(r.img1)                 # [rows, cols, layers]
    maxval    = terra::global(r.img1, "max", na.rm=TRUE)[[1]]
    threshold = tape_brightness * maxval
    
    # how many distinct colors are actually present? (cap k-means accordingly)
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
#' @param cor_type "phase" (phase correlation) or "ccf" (normalized cross-corr).
#' @param fixed_depth_pixel Depth band along COLUMNS. Length-2 = range start:end;
#'   longer = explicit column indices; NULL = use full width.
#' @param fixed_width Optional: restrict the ROTATION axis (rows), centered.
#' @param select_layer Layer to use for multi-band inputs.
#' @param window Demean + Hann-window before FFT to suppress edge artifacts.
#' @param overlay If TRUE, also draw a before/after magenta-green overlay.
#' @param overlay_layer Layer to display in the overlay (root mask, default 2).
#' @return Named numeric vector: depth (column lag), rotation (row lag), peak.
#' @export
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' data(seg_Oulanka2023_Session03_T067)
#' img1 <- terra::rast(seg_Oulanka2023_Session01_T067)
#' img2 <- terra::rast(seg_Oulanka2023_Session03_T067)
#' estimate_rotation_shift(img1, img2, cor_type = "phase", select_layer = 2)
estimate_rotation_shift <- function(
    img1, img2,
    cor_type = "phase",
    fixed_depth_pixel = NULL,
    fixed_width = NULL,
    select_layer = NULL,
    window = TRUE,
    overlay = FALSE,
    overlay_layer = 2
) {
  tryCatch({
    
    if (is.null(img1) || is.null(img2)) stop("Both input images are required")
    if (!cor_type %in% c("phase", "ccf")) stop("cor_type must be 'phase' or 'ccf'")
    
    im1 <- load_flexible_image(img1, select_layer = select_layer,
                               output_format = "array", scale = "none")
    im2 <- load_flexible_image(img2, select_layer = select_layer,
                               output_format = "array", scale = "none")
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
    if (!is.null(fixed_depth_pixel)) {
      d <- if (length(fixed_depth_pixel) == 2) {
        if (fixed_depth_pixel[1] >= fixed_depth_pixel[2])
          stop("fixed_depth_pixel[1] must be < fixed_depth_pixel[2]")
        seq(fixed_depth_pixel[1], fixed_depth_pixel[2])
      } else fixed_depth_pixel
      d <- d[d >= 1 & d <= ncol(g1)]
      if (length(d) < 8) stop("depth band too narrow after clipping to image")
      g1 <- g1[, d, drop = FALSE]; g2 <- g2[, d, drop = FALSE]
    }
    
    # optional rotation-axis (rows) restriction, centered
    if (!is.null(fixed_width)) {
      w  <- min(fixed_width, nrow(g1))
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
    
    a <- if (cor_type == "phase") imagefx::pcorr3d(g1, g2) else imagefx::xcorr3d(g1, g2)
    
    # pcorr3d/xcorr3d: max.shifts[1] = ROW lag (rotation), [2] = COLUMN lag (depth)
    out <- c(depth    = as.numeric(a$max.shifts[2]),
             rotation = as.numeric(a$max.shifts[1]),
             peak     = as.numeric(a$max.cor))
    
    if (isTRUE(overlay)) {
      f1 <- im1[, , overlay_layer]; f2 <- im2[, , overlay_layer]   # roots, 0 = bg
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
#' fixed-width window centered on a given row (e.g. the rotation center from
#' \code{estimate_rotation_center()}); in variable mode it trims a band sized by
#' a measured offset. Optionally previews what is kept versus cut.
#'
#' @param img Input image as raster, file name, or array.
#' @param center_offset Numeric or character. Where to center the kept window
#'   (in fixed mode), given in one of three forms:
#'   \itemize{
#'     \item \strong{Absolute row} - a number \code{> 1}: the exact row to
#'       center on (e.g. from \code{estimate_rotation_center()}). 
#'     \item \strong{Fraction} - a number in \code{[0, 1]}: a fraction of the
#'       image height. \code{0} = top, \code{0.25} = a quarter down,
#'       \code{0.5} = middle, \code{1} = bottom (so the center row is
#'       \code{center_offset * nrow}).
#'     \item \strong{Keyword} - \code{"top"} (= 0), \code{"middle"} /
#'       \code{"center"} / \code{"center"} (= 0.5), or \code{"bottom"} (= 1).
#'   }
#'   The default \code{0} centers on the top row.
#'   When \code{fixed_rotation = FALSE} the resolved value is used as the
#'   rotation shift in rows to trim (e.g. from
#'   \code{estimate_rotation_shift()}); pass an absolute number there.
#' @param cut_buffer Extra proportion of the rotation axis to trim (variable mode).
#' @param fixed_rotation Logical. If \code{TRUE}, return a fixed-width window.
#' @param fixed_width Output width in rows when \code{fixed_rotation = TRUE}.
#' @param select_layer Integer or \code{NULL}. Layer to use for multi-band inputs.
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
#' rotation_censor(img, center_offset = r0, fixed_width = 800,
#'                 fixed_rotation = TRUE, overlay = TRUE)
#'
#' # Same window centered on the middle of the image, via a fraction or keyword
#' rotation_censor(img, center_offset = 0.5,      fixed_width = 800)
#' rotation_censor(img, center_offset = "middle", fixed_width = 800)
#' # Top of the tube a quarter of the way down
#' rotation_censor(img, center_offset = 0.25, fixed_width = 800)
rotation_censor <- function(img, center_offset = 0, cut_buffer = 0.02,
                            fixed_rotation = TRUE, fixed_width = 500,
                            select_layer = NULL, overlay = FALSE, ...) {
  tryCatch({
    
    if (is.null(img)) stop("Input image is required")
    if (!(is.numeric(center_offset) || is.character(center_offset)) ||
        length(center_offset) != 1L)
      stop("center_offset must be a single number (absolute row or fraction in [0, 1]) or a keyword")
    if (!is.numeric(cut_buffer) || cut_buffer < 0 || cut_buffer > 1)
      stop("cut_buffer must be numeric between 0 and 1")
    if (!is.logical(fixed_rotation)) stop("fixed_rotation must be logical")
    if (!is.numeric(fixed_width) || fixed_width <= 0)
      stop("fixed_width must be positive numeric")
    if (!is.null(select_layer) && (!is.numeric(select_layer) || select_layer < 1))
      stop("select_layer must be NULL or positive integer")
    
    img.c <- load_flexible_image(img, select_layer = select_layer,
                                 output_format = "spatrast",
                                 scale = "none")
    if (is.null(img.c)) stop("Failed to load image")
    
    nr <- dim(img.c)[1]; nc <- dim(img.c)[2]    # rows = rotation axis

    # --- resolve center_offset (keyword / fraction / absolute row) to a row ---
    if (is.character(center_offset)) {
      key <- tolower(trimws(center_offset))
      center_offset <- switch(key,
        top    = 0,
        middle = 0.5, center = 0.5,
        bottom = 1,
        stop("center_offset keyword must be one of 'top', 'middle'/'center', ",
             "or 'bottom'; got '", center_offset, "'"))
    }
    # numbers in [0, 1] are a fraction of image height; > 1 is an absolute row
    if (center_offset >= 0 && center_offset <= 1)
      center_offset <- center_offset * nr

    offset    <- round(center_offset)
    buffer.px <- round(cut_buffer * nr)
    
    # --- determine the kept row window (lo:hi, 1 = top) ---
    if (!fixed_rotation) {
      cut.px <- abs(offset) + buffer.px
      if (cut.px >= nr) stop("Cut size too large, would remove entire image")
      if (offset > 0)      { lo <- cut.px + 1; hi <- nr }
      else if (offset < 0) { lo <- 1;          hi <- nr - cut.px }
      else { h <- floor(buffer.px / 2); lo <- 1 + h; hi <- nr - h }
    } else {
      mid      <- offset
      half     <- fixed_width / 2
      max.half <- min(mid - 1, nr - mid)
      if (half > max.half)
        message("fixed_width = ", fixed_width, " cannot be centered symmetrically on row ",
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
#' @param search_area Proportion of image to analyze
#' @param tape_thresh Minimum tape coverage ratio
#' @param dpi Image resolution
#' @param nclasses Number of clustering classes
#' @param inverse Invert detection for dark markers
#' @param tape_overlap Safety margin for tape (cm)
#' @param tape_brightness Brightness threshold for tape
#' @param extra_rows Additional analysis rows
#' @param select_layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `NULL`.
#' @param tape_quantile Brightness alignment quantile
#' @return data.frame with soil surface and tape end positions
#' @export
#'
#' @examples
#' img = rgb_Oulanka2023_Session03_T067
#' Soil0Estimates = estimate_soil_surface(img)
estimate_soil_surface = function(img, search_area=0.45, tape_thresh=0.33, dpi=150, nclasses=3,
                     inverse=FALSE, tape_overlap=0.5, tape_brightness=0.6,
                     extra_rows=100, tape_quantile=0.98, select_layer=NULL) {
  if (!requireNamespace("RStoolbox", quietly = TRUE))
    stop("Package 'RStoolbox' is required for estimate_soil_surface(). ",
         "Install it with: install.packages(\"RStoolbox\")")
  # Validation module
  tryCatch({
    # Check required input
    if (is.null(img)) {
      stop("Input image is required")
    }

    # Parameter validation
    if (!is.numeric(search_area) || search_area <= 0 || search_area > 1) {
      stop("search_area must be numeric between 0 and 1")
    }
    if (!is.numeric(tape_thresh) || tape_thresh <= 0 || tape_thresh > 1) {
      stop("tape_thresh must be numeric between 0 and 1")
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
    if (!is.numeric(tape_overlap) || tape_overlap < 0) {
      stop("tape_overlap must be non-negative numeric")
    }
    if (!is.numeric(tape_brightness) || tape_brightness <= 0) {
      stop("tape_brightness must be positive numeric")
    }
    if (!is.numeric(extra_rows) || extra_rows < 0) {
      stop("extra_rows must be non-negative numeric")
    }
    if (!is.numeric(tape_quantile) || tape_quantile <= 0 || tape_quantile > 1) {
      stop("tape_quantile must be numeric between 0 and 1")
    }
    if (!is.null(select_layer) && (!is.numeric(select_layer) || select_layer < 1)) {
      stop("select_layer must be NULL or positive integer")
    }

    # Load and validate image
    im <- load_flexible_image(img, select_layer=select_layer,
                              output_format="array", scale = "none")
    if (is.null(im)) {
      stop("Failed to load image")
    }
    if (length(dim(im)) != 3) {
      stop("Input image must be 3-dimensional array (RGB)")
    }

    # Adjust parameters for inverse mode
    if (inverse) {
      tape_quantile = 1 - tape_quantile
      tape_brightness = 1 / tape_brightness
    }

    # Create red line array with validation
    red.line = array(dim=c(dim(im)[1], extra_rows, dim(im)[3]))
    quantile_value = tryCatch({
      stats::quantile(im[,,1], tape_quantile, na.rm=TRUE)
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
    r.img1 = terra::crop(r.img1, c(0, search_area*terra::ext(r.img1)[2],
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
      threshold = tape_brightness * global_min
      valid_clusters = clust.center[clust.center < threshold]
      if (length(valid_clusters) == 0) {
        warning("No clusters found below brightness threshold")
        return(data.frame(soil0=NA, tape.end=NA))
      }
      clust = which(clust.center == min(valid_clusters))
    } else {
      global_max = terra::global(r.img1, "max")[[1]][1:3]
      threshold = tape_brightness * global_max
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

      if (current_ratio <= tape_thresh &&
          ratio_24 <= tape_thresh &&
          ratio_48 <= tape_thresh) {
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
    rw.ind = rw.ind - extra_rows - round(tape_overlap * dpi/2.54)
    rw.ind = round(rw.ind)
    tape.end = rw.ind + round(tape_overlap * dpi/2.54)
    tape.end = round(tape.end)

    # Validate final results
    if (rw.ind < 0 || tape.end < 0) {
      warning("Negative position values calculated")
    }

    out = data.frame(soil0=rw.ind, tape.end=tape.end)
    return(out)

  }, error = function(e) {
    stop(paste("Error in SoilSurfE:", e$message))
  })
}






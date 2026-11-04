########################

#' Orthogonal and diagonal connections per skeleton pixel
#'
#' Every connection between two 8-connected skeleton pixels is counted once,
#' at one of its two pixels: each pixel counts its connections to the left and
#' upper neighbor as orthogonal, and to the upper-left and lower-left neighbor
#' as diagonal. A pixel at a bend or a junction can therefore hold two
#' connections of one kind. Because each connection sits on exactly one pixel,
#' summing the maps over any set of pixels -- a depth bin, a diameter class --
#' gives that set's counts, and the counts of the parts add up to the whole.
#'
#' @param img Single-layer 0/1 SpatRaster.
#' @return Two-layer SpatRaster, \code{orth} and \code{diag}, holding 0, 1 or 2
#'   per pixel.
#' @keywords internal
#' @noRd
.skeleton_links <- function(img) {
  # Each kernel holds the pixel and one neighbor; a focal sum of 2 means both
  # are skeleton, i.e. one connection.
  left       <- matrix(c(0,1,0, 0,1,0, 0,0,0), 3, 3)
  up         <- matrix(c(0,0,0, 1,1,0, 0,0,0), 3, 3)
  up_left    <- matrix(c(1,0,0, 0,1,0, 0,0,0), 3, 3)
  down_left  <- matrix(c(0,0,1, 0,1,0, 0,0,0), 3, 3)
  pair <- function(w) terra::focal(img, w = w, fun = "sum", na.rm = TRUE) == 2

  orth <- pair(left) + pair(up)
  diag <- pair(up_left) + pair(down_left)
  out <- c(orth, diag)
  names(out) <- c("orth", "diag")
  out
}


#' Root length in pixels from diagonal and orthogonal connection counts
#'
#' Vectorized over \code{nd} and \code{no}, so it turns per-bin counts into
#' per-bin lengths in one call.
#'
#' @param nd,no Diagonal and orthogonal connection counts.
#' @param method One of \code{"kimura2"}, \code{"kimura1"},
#'   \code{"freeman_basic"}, \code{"freeman_corrected"}.
#' @return Length in pixels.
#' @keywords internal
#' @noRd
.length_estimate <- function(nd, no, method = "kimura2") {
  switch(method,
    freeman_basic     = sqrt(2) * nd + no,
    freeman_corrected = 0.948 * (sqrt(2) * nd + no),
    kimura1           = sqrt(nd^2 + (nd + no)^2),
    kimura2           = sqrt(nd^2 + (nd + no / 2)^2) + no / 2,
    stop("Unknown length method: ", method)
  )
}


#' Root length estimation from skeleton images
#'
#' Root length is estimated from a skeletonized binary image using
#' either Freeman chain-code based estimators or Kimura estimators.
#'
#' Let Nd be the number of diagonal pixel connections and No the number
#' of orthogonal pixel connections in the skeleton. Every connection between
#' two 8-connected skeleton pixels is counted once, so a pixel at a bend or a
#' junction contributes each of its connections.
#'
#' Freeman methods treat the skeleton as a discrete chain-code path:
#' - freeman_basic:
#'   L = sqrt(2) * Nd + No
#'
#' - freeman_corrected:
#'   L = 0.948 * (sqrt(2) * Nd + No)
#'
#' Kimura methods treat the skeleton as a discretized representation of
#' an underlying continuous curve and reduce orientation bias:
#'
#' - kimura1:
#'   L = sqrt(Nd^2 + (Nd + No)^2)
#'
#' - kimura2 (default):
#'   L = sqrt(Nd^2 + (Nd + No/2)^2) + No/2
#'
#' The Kimura2 estimator is generally preferred due to improved stability
#' across object orientations and curvature distributions.
#'
#' @param method Character. One of:
#'   "freeman_basic", "freeman_corrected", "kimura1", "kimura2"
#' @param img Skeletonized binary raster image. If \code{skeletonize = TRUE},
#'   a segmented (non-skeleton) mask can be supplied instead.
#' @param unit Character. Output unit: \code{"cm"}, \code{"inch"}, or \code{"px"}.
#' @param dpi Numeric. Scan resolution (dots per inch). Required when
#'   \code{unit} is \code{"cm"} or \code{"inch"}; not used for \code{"px"}.
#' @param select_layer Numeric. Which layer to measure if \code{img} has
#'   multiple layers. Without it a multi-layer image is measured on layer 1,
#'   with a warning: in the bundled RootDetector images layer 1 also carries
#'   the tape, and the skeleton is layer 2.
#' @param show_messages Logical. If \code{TRUE}, prints informational messages during processing.
#' @param skeletonize Logical. If \code{TRUE}, \code{img} is treated as a
#'   segmented mask and reduced to a skeleton internally via
#'   \code{skeletonize_image()} before computing length. Default
#'   \code{FALSE} (assumes \code{img} is already a skeleton).
#' @return Root length in pixels or converted units
#' @export
#' @examples
#' data(skl_Oulanka2023_Session01_T067)
#' skl <- terra::rast(skl_Oulanka2023_Session01_T067)
#' root_length(skl, unit = "cm", dpi = 150, select_layer = 2)
root_length <- function(img,
                        unit = "cm",
                        dpi = NULL,
                        select_layer = NULL,
                        method = c("kimura2",
                                   "kimura1",
                                   "freeman_basic",
                                   "freeman_corrected"),
                        show_messages = TRUE,
                        skeletonize = FALSE) {
  
  method <- match.arg(method)
  
  tryCatch({

    if (missing(img)) {
      stop("Image input is required")
    }
    
    if (!unit %in% c("px", "cm", "inch")) {
      stop("unit must be 'px', 'cm', or 'inch'")
    }
    
    if (unit %in% c("cm", "inch")) {
      if (is.null(dpi))
        stop("`dpi` is required for unit = \"", unit, "\": pass the scan ",
             "resolution of the image, e.g. dpi = 300, or use unit = \"px\"")
      if (!is.numeric(dpi) || length(dpi) != 1L || is.na(dpi) || dpi <= 0)
        stop("`dpi` must be a single positive number")
    }
    
    # -----------------------------
    # Load image
    # -----------------------------
    img <- load_flexible_image(
      img,
      select_layer = select_layer,
      output_format = "spatrast",
      scale = "binary"
    )
    
    if (is.null(img) || terra::nlyr(img) < 1) {
      stop("Invalid image after loading")
    }

    if (terra::nlyr(img) > 1) {
      warning(sprintf(paste("img has %d layers and select_layer is not set;",
                            "measuring layer 1. Layer 1 is often not the",
                            "skeleton: set select_layer to the layer that",
                            "holds it."), terra::nlyr(img)), call. = FALSE)
      img <- img[[1]]
    }
    
    # -----------------------------
    # Early exit: no roots present
    # -----------------------------
    rmax <- terra::global(img, "max", na.rm = TRUE)[1, 1]
    
    if (is.na(rmax) || rmax == 0) {
      return(0)
    }

    # -----------------------------
    # Skeletonize if a segmented (non-skeleton) mask was supplied
    # -----------------------------
    if (skeletonize) {
      img <- skeletonize_image(img, verbose = FALSE)
    }

    vals <- unique(terra::values(img))
    vals <- vals[!is.na(vals)]
    
    if (!all(vals %in% c(0, 1))) {
      warning("Non-binary image detected after preprocessing")
    }
    
    # -----------------------------
    # Connection counts and length estimate
    # -----------------------------
    links <- .skeleton_links(img)
    n_orth <- terra::global(links[["orth"]], "sum", na.rm = TRUE)[1, 1]
    n_diag <- terra::global(links[["diag"]], "sum", na.rm = TRUE)[1, 1]

    if (show_messages) {
      message("Diagonal: ", n_diag, " | Orthogonal: ", n_orth)
    }

    rootlength_px <- .length_estimate(n_diag, n_orth, method)
    
    # -----------------------------
    # Unit conversion
    # -----------------------------
    rootlength <- switch(
      unit,
      
      px = rootlength_px,
      cm = rootlength_px / dpi * 2.54,
      inch = rootlength_px / dpi
    )
    
    rootlength <- as.numeric(rootlength)
    
    if (is.na(rootlength) || rootlength < 0) {
      stop("Invalid root length computed")
    }
    
    return(rootlength)
    
  }, error = function(e) {
    stop("root_length failed: ", e$message)
  })
}


## RootScape Metrics

#' RootScapeMetric relies on Landscapemetrics to extract 'Root Scape' Features
#' akin to landscape analysis.
#'
#' @param img Segmented raster (values = 0, 1). Consider whether a skeletonized
#'   raster is more appropriate for your use case.
#' @param index_d Depth index for the output column "depth". Useful when called
#'   inside a loop over depth bins.
#' @param metrics Character vector of metrics to calculate; must be valid names
#'   from `landscapemetrics::list_lsm()`.
#' @param select_layer Integer. Specifies which layer to use if the input is a
#'   multi-band image. Default is `NULL` (single-layer expected).
#' @return A data frame in wide format: columns \code{object} and
#'   \code{depth}, then one column per metric named as in \code{metrics}.
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img <- terra::rast(seg_Oulanka2023_Session01_T067)
#' if (requireNamespace("landscapemetrics", quietly = TRUE)) {
#'   RootScapeObject <- root_scape_metrics(img, index_d = 80, select_layer = 2,
#'                                         metrics = c("lsm_c_ca"))
#' }
root_scape_metrics <- function(img, index_d = NA, select_layer = NULL,
                                metrics = c("lsm_c_ca", "lsm_l_ent", "lsm_c_pd",
                                            "lsm_c_np", "lsm_c_pland",
                                            "lsm_c_area_mn", "lsm_c_area_cv",
                                            "lsm_c_enn_mn", "lsm_c_enn_cv")) {
  if (!requireNamespace("landscapemetrics", quietly = TRUE))
    stop("Package 'landscapemetrics' is required for root_scape_metrics(). ",
         "Install it with: install.packages(\"landscapemetrics\")")
  tryCatch({
    if (missing(img)) stop("Image input is required")

    if (!is.character(metrics) || length(metrics) == 0) {
      stop("metrics must be a non-empty character vector")
    }

    available_metrics <- landscapemetrics::list_lsm()$function_name
    invalid_metrics   <- metrics[!metrics %in% available_metrics]
    if (length(invalid_metrics) > 0) {
      stop("Invalid metrics specified: ", paste(invalid_metrics, collapse = ", "))
    }

    .validate_select_layer(select_layer)

    img <- load_flexible_image(img, select_layer = select_layer,
                               output_format = "spatrast", scale = "none")

    if (is.null(img) || terra::nlyr(img) < 1) {
      stop("Invalid or empty image after loading")
    }

    rsm <- tryCatch(
      landscapemetrics::calculate_lsm(img, directions = 8,
                                      neighbourhood = 8, what = metrics),
      error = function(e) stop("Error calculating landscape metrics: ", e$message)
    )

    if (is.null(rsm) || nrow(rsm) == 0) stop("No valid metrics could be calculated")

    t.object     <- ifelse(rsm$class == 0, "deletable", "root")
    t.object     <- ifelse(is.na(rsm$class), "root", t.object)
    rsm$object   <- t.object
    rsm$depth    <- index_d

    rsm <- dplyr::distinct(rsm)
    rsm <- dplyr::filter(rsm, rsm$object != "deletable")
    rsm$id    <- NULL
    rsm$class <- NULL
    rsm$level <- NULL
    rsm$layer <- NULL
    
    # wide format
    rsm_wide <- tidyr::pivot_wider(
      rsm,
      names_from = "metric",
      values_from = "value",
      values_fill = NA_real_
    )

    return(rsm_wide)

  }, error = function(e) {
    stop("Error in root_scape_metrics: ", e$message)
  })
}


#' Count all pixels in a segmented image
#'
#' @param img A single-layer raster image (SpatRaster or compatible format).
#' @return A numeric value — the sum of all non-NA pixel values.
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img <- terra::rast(seg_Oulanka2023_Session01_T067)[[2]]
#' rootpixel <- count_pixels(img)
count_pixels <- function(img) {
  tryCatch({
    if (missing(img)) stop("Image input is required")

    img <- load_flexible_image(img, output_format = "spatrast", scale = "binary")

    if (is.null(img) || terra::nlyr(img) < 1) {
      stop("Invalid or empty image after loading")
    }

    srpx <- terra::global(img, "sum", na.rm = TRUE)

    if (is.null(srpx) || length(srpx) == 0) stop("Error calculating pixel sum")
    if (any(srpx[[1]] < 0)) {
      warning("Negative pixel sum calculated, which may indicate an issue with the input image")
    }

    return(srpx[[1]])

  }, error = function(e) {
    stop("Error in count_pixels: ", e$message)
  })
}



#' Calculate Image Coloration Metrics
#'
#' @param img Three-band raster (RGB) or path to image.
#' @param r Red channel luminosity weight. Default follows ITU-R BT.709.
#' @param g Green channel luminosity weight.
#' @param b Blue channel luminosity weight.
#' @return A data frame with columns: rcc, gcc, bcc, hue, saturation,
#'   luminosity, red, green, blue. \code{hue} and \code{saturation} come from
#'   the HSV transform of the mean color; \code{luminosity} is the
#'   \code{r}/\code{g}/\code{b}-weighted mean luminance (not HSV value).
#' @export
#'
#' @examples
#' data(rgb_Oulanka2023_Session03_T067)
#' img <- terra::rast(rgb_Oulanka2023_Session03_T067)
#' colorvector <- tube_coloration(img)
tube_coloration <- function(img, r = 0.2126, g = 0.7152, b = 0.0722) {
  # These are the exact Rec. 709 luma coefficients, not the rounded
  # 0.21/0.72/0.07 that rgb2gray() defaults to, so the two stay separate even
  # though the arithmetic is the same.
  tryCatch({
    if (missing(img)) stop("Image input is required")

    if (!is.numeric(r) || !is.numeric(g) || !is.numeric(b) ||
        r < 0 || g < 0 || b < 0 || abs((r + g + b) - 1) > 1e-9) {
      stop("RGB weights must be non-negative numbers that sum to 1")
    }

    img <- load_flexible_image(img, output_format = "spatrast",
                               scale = "none",
                               select_layer = NULL)

    if (is.null(img) || terra::nlyr(img) != 3) {
      stop("Input must be a valid three-band (RGB) image")
    }

    vr <- terra::values(img[[1]])
    vg <- terra::values(img[[2]])
    vb <- terra::values(img[[3]])

    if (all(is.na(vr)) || all(is.na(vg)) || all(is.na(vb))) {
      stop("One or more color channels contain no valid data")
    }

    mean.r <- mean(vr, na.rm = TRUE)
    mean.g <- mean(vg, na.rm = TRUE)
    mean.b <- mean(vb, na.rm = TRUE)

    if (any(c(mean.r, mean.g, mean.b) > 255) ||
        any(c(mean.r, mean.g, mean.b) < 0)) {
      warning("Color values outside expected range [0, 255]")
    }

    hsl <- tryCatch(
      grDevices::rgb2hsv(r = mean.r, g = mean.g, b = mean.b),
      error = function(e) stop("Error converting RGB to HSV: ", e$message)
    )

    intensity <- vr + vg + vb
    if (any(intensity == 0, na.rm = TRUE)) {
      message("Some pixels have zero intensity, which may affect color calculations")
    }

    # Rec.709 luminance from the caller's channel weights. HSV *value* (the
    # max channel) is not luminance and ignores r/g/b, so `luminosity` is the
    # weighted mean instead.
    lum.gray <- vr * r + vg * g + vb * b
    mean.lum <- round(mean(lum.gray, na.rm = TRUE), 4)

    rcc <- round(mean(vr / intensity, na.rm = TRUE), 4)
    gcc <- round(mean(vg / intensity, na.rm = TRUE), 4)
    bcc <- round(mean(vb / intensity, na.rm = TRUE), 4)

    colordf <- data.frame(
      rcc        = rcc,
      gcc        = gcc,
      bcc        = bcc,
      hue        = hsl[1],
      saturation = hsl[2],
      luminosity = mean.lum,
      red        = mean.r,
      green      = mean.g,
      blue       = mean.b
    )

    return(colordf)

  }, error = function(e) {
    stop("Error in tube_coloration: ", e$message)
  })
}


#' Texture calculation using Gray-Level Co-occurrence Matrix (GLCM)
#'
#' @param img.color Three-band (RGB) image: a `SpatRaster`, `RasterBrick`,
#'   `RasterStack`, or any other input \code{\link{load_flexible_image}}
#'   accepts, including a file path. Internally converted to a
#'   `raster::RasterBrick` as required by the `glcm` package.
#' @param grays Number of gray levels for GLCM quantization. Must be between
#'   2 and 255. Default is 7.
#' @param window Window size for GLCM calculation as a length-2 vector of odd
#'   positive integers, e.g. `c(9, 9)`. Default is `c(9, 9)`.
#' @param metrics Character vector of GLCM texture statistics to calculate.
#'   Valid options: "mean", "variance", "homogeneity", "contrast",
#'   "dissimilarity", "entropy", "second_moment", "correlation".
#' @return A RasterLayer (or RasterBrick for multiple metrics) with texture
#'   values.
#' @export
#'
#' @examples
#' data(rgb_Oulanka2023_Session03_T067)
#' img <- raster::brick(rgb_Oulanka2023_Session03_T067)
#' analyze_soil_texture(img, grays = 7, window = c(9, 9),
#'                      metrics = "second_moment")
analyze_soil_texture <- function(img.color, grays = 7, window = c(9, 9),
                                  metrics = c("variance", "second_moment")) {
  tryCatch({
    if (missing(img.color)) stop("Image input is required")

    if (!is.numeric(grays) || grays < 2 || grays > 255) {
      stop("grays must be a numeric value between 2 and 255")
    }

    if (!is.numeric(window) || length(window) != 2 ||
        any(window < 1) || any(window %% 2 == 0)) {
      stop("window must be a vector of two odd positive integers")
    }

    valid_metrics   <- c("mean", "variance", "homogeneity", "contrast",
                         "dissimilarity", "entropy", "second_moment", "correlation")
    invalid_metrics <- metrics[!metrics %in% valid_metrics]
    if (length(invalid_metrics) > 0) {
      stop("Invalid metrics specified: ", paste(invalid_metrics, collapse = ", "))
    }

    # glcm requires a raster::RasterBrick. Raster objects are converted
    # directly, which keeps their extent; anything else (file path, array,
    # cimg, ...) goes through load_flexible_image().
    if (inherits(img.color, c("SpatRaster", "RasterStack"))) {
      img.color <- raster::brick(img.color)
    } else if (!inherits(img.color, c("RasterLayer", "RasterBrick"))) {
      img.color <- load_flexible_image(img.color, output_format = "brick",
                                       scale = "none")
    }

    if (raster::nlayers(img.color) != 3) {
      stop("Input must be a valid three-band (RGB) image")
    }

    # Same weights as rgb2gray(), but inlined: this operates on a raster::brick
    # and divides by the value range, which rgb2gray()'s SpatRaster path does
    # not do. glcm::glcm() needs the [0, 1] input, hence the /mx.
    mx      <- max(raster::values(img.color), na.rm = TRUE)
    mx      <- if (mx > 1) 255 else 1
    img.gray <- (img.color[[1]] * 0.21 + img.color[[2]] * 0.72 +
                   img.color[[3]] * 0.07) / mx

    if (all(is.na(raster::values(img.gray)))) {
      stop("Grayscale conversion resulted in invalid data")
    }

    tx.im <- tryCatch(
      glcm::glcm(img.gray, n_grey = grays, window = window, statistics = metrics),
      error = function(e) stop("Error calculating texture metrics: ", e$message)
    )

    if (is.null(tx.im)) stop("Texture calculation produced no valid results")

    return(tx.im)

  }, error = function(e) {
    stop("Error in analyze_soil_texture: ", e$message)
  })
}




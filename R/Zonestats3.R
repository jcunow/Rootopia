########################

#' Root length estimation from skeleton images
#'
#' Root length is estimated from a skeletonized binary image using
#' either Freeman chain-code based estimators or Kimura estimators.
#'
#' Let Nd be the number of diagonal pixel connections and No the number
#' of orthogonal pixel connections in the skeleton.
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
#' @param skeletonize Logical. If \code{TRUE}, \code{img} is treated as a
#'   segmented mask and reduced to a skeleton internally via
#'   \code{skeletonize_image()} before computing length. Default
#'   \code{FALSE} (assumes \code{img} is already a skeleton).
#' @return Root length in pixels or converted units
root_length <- function(img,
                        unit = "cm",
                        dpi = 300,
                        select.layer = NULL,
                        method = c("kimura2",
                                   "kimura1",
                                   "freeman_basic",
                                   "freeman_corrected"),
                        show_messages = TRUE,
                        skeletonize = FALSE) {
  
  method <- match.arg(method)
  
  tryCatch({
    
    # -----------------------------
    # Input validation
    # -----------------------------
    if (missing(img)) {
      stop("Image input is required")
    }
    
    if (!unit %in% c("px", "cm", "inch")) {
      stop("unit must be 'px', 'cm', or 'inch'")
    }
    
    if (unit %in% c("cm", "inch")) {
      if (missing(dpi) || !is.numeric(dpi) || dpi <= 0) {
        stop("Valid positive dpi required for cm/inch conversion")
      }
    }
    
    # -----------------------------
    # Load image
    # -----------------------------
    img <- load_flexible_image(
      img,
      select.layer = select.layer,
      output_format = "spatrast",
      normalize = TRUE,
      binarize = TRUE
    )
    
    if (is.null(img) || terra::nlyr(img) < 1) {
      stop("Invalid image after loading")
    }

    # ensure single layer safely
    if (terra::nlyr(img) > 1) {
      img <- img[[1]]
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
    # Connectivity kernels
    # -----------------------------
    k0 <- matrix(c(0,1,0,
                   0,1,0,
                   0,0,0), 3, 3)
    
    k1 <- matrix(c(0,0,0,
                   1,1,0,
                   0,0,0), 3, 3)
    
    g0 <- matrix(c(1,0,0,
                   0,1,0,
                   0,0,0), 3, 3)
    
    g1 <- matrix(c(0,0,1,
                   0,1,0,
                   0,0,0), 3, 3)
    
    # -----------------------------
    # Focal operations
    # -----------------------------
    r0 <- terra::focal(img, w = k0, fun = "sum", na.rm = TRUE)
    r1 <- terra::focal(img, w = k1, fun = "sum", na.rm = TRUE)
    
    u0 <- terra::focal(img, w = g0, fun = "sum", na.rm = TRUE)
    u1 <- terra::focal(img, w = g1, fun = "sum", na.rm = TRUE)
    
    # -----------------------------
    # IMPORTANT: scalar extraction (fixes your error source)
    # -----------------------------
    orth.img <- (r0 == 2) | (r1 == 2)
    diag.img <- (u0 == 2) | (u1 == 2)
    
    kimura.sum.orth <- terra::global(orth.img, "sum", na.rm = TRUE)[1,1]
    kimura.sum.diag <- terra::global(diag.img, "sum", na.rm = TRUE)[1,1]
    
    if (show_messages) {
      message(
        "Diagonal: ", kimura.sum.diag,
        " | Orthogonal: ", kimura.sum.orth
      )
    }
    
    # -----------------------------
    # Length estimators
    # -----------------------------
    rootlength_px <- switch(
      method,
      
      freeman_basic =
        sqrt(2) * kimura.sum.diag + kimura.sum.orth,
      
      freeman_corrected =
        0.948 * (sqrt(2) * kimura.sum.diag + kimura.sum.orth),
      
      kimura1 =
        sqrt(kimura.sum.diag^2 +
               (kimura.sum.diag + kimura.sum.orth)^2),
      
      kimura2 =
        sqrt(kimura.sum.diag^2 +
               (kimura.sum.diag + kimura.sum.orth/2)^2) +
        kimura.sum.orth/2
    )
    
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
#' @param indexD Depth index for the output column "depth". Useful when called
#'   inside a loop over depth bins.
#' @param metrics Character vector of metrics to calculate; must be valid names
#'   from `landscapemetrics::list_lsm()`.
#' @param select.layer Integer. Specifies which layer to use if the input is a
#'   multi-band image. Default is `NULL` (single-layer expected).
#' @import dplyr
#' @return A data frame of metric values with columns: metric, value, object,
#'   depth.
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img <- terra::rast(seg_Oulanka2023_Session01_T067)
#' RootScapeObject <- root_scape_metrics(img, indexD = 80, select.layer = 2,
#'                                        metrics = c("lsm_c_ca"))
root_scape_metrics <- function(img, indexD = NA, select.layer = NULL,
                                metrics = c("lsm_c_ca", "lsm_l_ent", "lsm_c_pd",
                                            "lsm_c_np", "lsm_c_pland",
                                            "lsm_c_area_mn", "lsm_c_area_cv",
                                            "lsm_c_enn_mn", "lsm_c_enn_cv")) {
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

    if (!is.null(select.layer)) {
      if (!is.numeric(select.layer) || select.layer < 1) {
        stop("select.layer must be a positive integer")
      }
    }

    img <- load_flexible_image(img, select.layer = select.layer,
                               output_format = "spatrast", normalize = FALSE)

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
    rsm$depth    <- indexD

    rsm <- dplyr::distinct(rsm)
    rsm <- dplyr::filter(rsm, rsm$object != "deletable")
    rsm$id    <- NULL
    rsm$class <- NULL
    rsm$level <- NULL
    rsm$layer <- NULL

    return(rsm)

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

    img <- load_flexible_image(img, output_format = "spatrast", normalize = TRUE)

    if (is.null(img) || terra::nlyr(img) < 1) {
      stop("Invalid or empty image after loading")
    }

    srpx <- terra::global(img, "sum", na.rm = TRUE)

    if (is.null(srpx) || length(srpx) == 0) stop("Error calculating pixel sum")
    if (srpx[[1]] < 0) {
      warning("Negative pixel sum calculated, which may indicate an issue with the input image")
    }

    return(srpx[[1]])

  }, error = function(e) {
    stop("Error in count_pixels: ", e$message)
  })
}

#' Count pixels (deprecated alias for count_pixels)
#'
#' @description `px.sum()` is a deprecated alias for [count_pixels()]. Please
#'   update your code to use `count_pixels()` instead.
#' @param img A single-layer raster image.
#' @param layer Ignored. Provided for backward compatibility only; use
#'   [terra::subset()] to select a layer before calling this function.
#' @return A numeric value — the sum of all non-NA pixel values.
#' @export
px.sum <- function(img, layer = NULL) {
  .Deprecated("count_pixels",
              msg = "px.sum() is deprecated. Use count_pixels() instead. To select a layer, use terra::subset(img, layer) before calling count_pixels().")
  if (!is.null(layer)) {
    img <- terra::subset(img, layer)
  }
  count_pixels(img)
}


#' Calculate Image Coloration Metrics
#'
#' @param img Three-band raster (RGB) or path to image.
#' @param r Red channel luminosity weight. Default follows ITU-R BT.709.
#' @param g Green channel luminosity weight.
#' @param b Blue channel luminosity weight.
#' @return A data frame with columns: rcc, gcc, bcc, hue, saturation,
#'   luminosity, red, green, blue.
#' @export
#'
#' @examples
#' data(rgb_Oulanka2023_Session03_T067)
#' img <- terra::rast(rgb_Oulanka2023_Session03_T067)
#' colorvector <- tube_coloration(img)
tube_coloration <- function(img, r = 0.2126, g = 0.7152, b = 0.0722) {
  tryCatch({
    if (missing(img)) stop("Image input is required")

    if (!is.numeric(r) || !is.numeric(g) || !is.numeric(b) ||
        r < 0 || g < 0 || b < 0 || abs((r + g + b) - 1) > 1e-9) {
      stop("RGB weights must be non-negative numbers that sum to 1")
    }

    img <- load_flexible_image(img, output_format = "spatrast",
                               normalize = FALSE, binarize = FALSE,
                               select.layer = NULL)

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
      message("Some pixels have zero intensity, which may affect colour calculations")
    }

    lum.gray       <- vr * r + vg * g + vb * b
    mean.intensity <- round(mean(intensity,  na.rm = TRUE), 4)
    mean.lum       <- round(mean(lum.gray,   na.rm = TRUE), 4)

    rcc <- round(mean(vr / intensity, na.rm = TRUE), 4)
    gcc <- round(mean(vg / intensity, na.rm = TRUE), 4)
    bcc <- round(mean(vb / intensity, na.rm = TRUE), 4)

    colordf <- data.frame(
      rcc        = rcc,
      gcc        = gcc,
      bcc        = bcc,
      hue        = hsl[1],
      saturation = hsl[2],
      luminosity = hsl[3],
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
#' @param img.color Three-band raster (RGB) or path to image. Internally
#'   converted to a `raster::RasterBrick` as required by the `glcm` package.
#' @param grays Number of gray levels for GLCM quantization. Must be between
#'   2 and 255. Default is 7.
#' @param window Window size for GLCM calculation as a length-2 vector of odd
#'   positive integers, e.g. `c(9, 9)`. Default is `c(9, 9)`.
#' @param metrics Character vector of GLCM texture statistics to calculate.
#'   Valid options: "mean", "variance", "homogeneity", "contrast",
#'   "dissimilarity", "entropy", "second_moment", "correlation".
#' @return A RasterLayer (or RasterBrick for multiple metrics) with texture
#'   values.
#' @import raster
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

    # glcm requires a raster object; coerce from terra if needed
    if (inherits(img.color, "SpatRaster")) {
      img.color <- raster::brick(img.color)
    }

    if (!inherits(img.color, c("RasterLayer", "RasterBrick", "RasterStack"))) {
      stop("img.color must be a raster object (RasterLayer, RasterBrick, RasterStack, or SpatRaster)")
    }

    if (raster::nlayers(img.color) != 3) {
      stop("Input must be a valid three-band (RGB) image")
    }

    mx      <- max(raster::values(img.color), na.rm = TRUE)
    mx      <- if (mx > 1) 255 else 1
    img.gray <- (img.color[[1]] * 0.21 + img.color[[2]] * 0.72 +
                   img.color[[3]] * 0.07) / mx

    if (all(is.na(raster::values(img.gray)))) {
      stop("Grayscale conversion resulted in invalid data")
    }

    tx.im <- tryCatch(
      glcm::glcm(img.gray, n_grey = grays, statistics = metrics),
      error = function(e) stop("Error calculating texture metrics: ", e$message)
    )

    if (is.null(tx.im)) stop("Texture calculation produced no valid results")

    return(tx.im)

  }, error = function(e) {
    stop("Error in analyze_soil_texture: ", e$message)
  })
}


#' Approximate average root thickness (deprecated)
#'
#' @description `root_thickness()` is a naive estimator that back-calculates
#'   an average diameter from a total root length and a total root pixel
#'   count (`area / length`), assuming every pixel belongs to a single root of
#'   uniform width. It ignores branching, overlapping roots, and the actual
#'   local width distribution.
#'
#'   [root_diameter()] computes per-pixel diameters directly from the
#'   distance transform of the segmented mask, and returns `mean_diameter`,
#'   `median_diameter`, `quantiles`, `root_volume`, and `root_surface_area`
#'   (lateral surface area, assuming cylindrical root segments). Use those
#'   instead of `root_thickness()`.
#'
#' @param kimuralength Total root length in cm (e.g. from [root_length()]).
#' @param rootpx Total number of root pixels in the image section.
#' @param dpi Image resolution in dots per inch. Default is 300.
#' @return A numeric value in cm representing approximate average root diameter.
#' @export
#'
#' @examples
#' root.thicc <- root_thickness(kimuralength = 300, rootpx = 9500, dpi = 300)
root_thickness <- function(kimuralength, rootpx, dpi = 300) {
  .Deprecated("root_diameter",
              msg = paste("root_thickness() is a naive area/length estimator and is deprecated.",
                          "Use root_diameter() for per-pixel diameter, root_volume, and",
                          "root_surface_area computed from the distance transform."))
  tryCatch({
    if (missing(kimuralength) || missing(rootpx)) {
      stop("Both kimuralength and rootpx are required")
    }

    if (!is.numeric(kimuralength) || !is.numeric(rootpx) || !is.numeric(dpi)) {
      stop("All inputs must be numeric values")
    }

    if (kimuralength <= 0 || rootpx <= 0 || dpi <= 0) {
      stop("All inputs must be positive values")
    }

    px.length <- kimuralength * (dpi / 2.54)
    if (px.length <= 0) stop("Invalid pixel length calculation")

    thiccness <- rootpx / px.length
    if (thiccness <= 0) stop("Invalid thickness calculation")

    px.thicc <- thiccness / (dpi / 2.54)

    if (px.thicc <= 1e-9) {
      warning("Calculated thickness is very small, which may indicate measurement issues")
    }

    return(px.thicc)

  }, error = function(e) {
    stop("Error in root_thickness: ", e$message)
  })
}

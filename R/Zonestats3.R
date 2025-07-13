
########################
#' Calculate Root Length using Kimura's Method with optimizations
#'
#' @param img A skeletonized root image raster
#' @param unit Output unit ("px" or "cm")
#' @param dpi Image resolution (required when unit = "cm" or "inch")
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `2`.
#' @return Numeric value representing root length in specified unit
#' @export
#'
#' @examples
#' data(skl_Oulanka2023_Session01_T067)
#' img = terra::rast(skl_Oulanka2023_Session01_T067)
#' RL = root_length(img = img,unit = "cm", dpi = 300, select.layer = 2)
root_length = function(img, unit="cm", dpi=300, select.layer = 1) {
  # Input validation module
  tryCatch({
    # Check if img is provided
    if (missing(img)) {
      stop("Image input is required")
    }
    
    # Check if image is valid after loading
    img <- load_flexible_image(img, select.layer = select.layer,
                               output_format = "spatrast", normalize = TRUE, binarize = TRUE)
    

    # Validate unit parameter
    if (!unit %in% c("px", "cm", "inch")) {
      stop("Unit must be either 'px', 'inch' or 'cm'")
    }

    # Validate DPI when unit is cm
    if (unit == "cm" | unit =="inch") {
      if (missing(dpi) || !is.numeric(dpi) || dpi <= 0) {
        stop("Valid positive numeric dpi value is required when unit = 'cm'")
      }
    }



    if((is.null(select.layer) || select.layer<1) && terra::nlyr(img)>1){
      stop("Multiple layers present, select.layer should be a positive integer")
    }


    if (is.null(img) || terra::nlyr(img) < 1) {
      stop("Invalid or empty image after loading")
    }

    # Check if image contains any non-NA values
    if (all(is.na(terra::values(img)))) {
      stop("Image contains no valid data (all NA values)")
    }

    # Check if image contains only binary values after normalization
    vals <- unique(terra::values(img))
    vals <- vals[!is.na(vals)]
    if (!all(vals %in% c(0, 1))) {
      warning("Image may not be properly skeletonized (contains non-binary values)")
    }

    ## Original function logic starts here
    k0 = matrix(c(0,1,0,0,1,0,0,0,0), nrow = 3, ncol = 3)
    k1 = matrix(c(0,0,0,1,1,0,0,0,0), nrow = 3, ncol = 3)

    # Add error handling for focal operations
    r0 <- tryCatch({
      terra::focal(img, w = k0, fun = "sum", na.rm=TRUE)
    }, error = function(e) {
      stop("Error in focal operation with k0 matrix: ", e$message)
    })

    r1 <- tryCatch({
      terra::focal(img, w = k1, fun = "sum", na.rm=TRUE)
    }, error = function(e) {
      stop("Error in focal operation with k1 matrix: ", e$message)
    })

    orth.img = sum((r0 == 2) | (r1 == 2))

    g0 = matrix(c(1,0,0,0,1,0,0,0,0), nrow = 3, ncol = 3)
    g1 = matrix(c(0,0,1,0,1,0,0,0,0), nrow = 3, ncol = 3)

    u0 <- tryCatch({
      terra::focal(img, w = g0, fun = "sum", na.rm=TRUE)
    }, error = function(e) {
      stop("Error in focal operation with g0 matrix: ", e$message)
    })

    u1 <- tryCatch({
      terra::focal(img, w = g1, fun = "sum", na.rm=TRUE)
    }, error = function(e) {
      stop("Error in focal operation with g1 matrix: ", e$message)
    })

    diag.img = sum((u0 == 2) | (u1 == 2))

    # Add error handling for global operations
    kimura.sum.diag <- tryCatch({
      terra::global(diag.img, "sum", na.rm=TRUE)
    }, error = function(e) {
      stop("Error calculating diagonal sum: ", e$message)
    })

    kimura.sum.orth <- tryCatch({
      terra::global(orth.img, "sum", na.rm=TRUE)
    }, error = function(e) {
      stop("Error calculating orthogonal sum: ", e$message)
    })

    # Calculate root length 
    if(unit == "px") {
      rootlength = round(( kimura.sum.diag^2 + (kimura.sum.diag + kimura.sum.orth/2)^2 )^0.5 +
                           kimura.sum.orth/2)[[1]]
    } else if(unit == "cm") {
      rootlength = round((( kimura.sum.diag^2 + (kimura.sum.diag + kimura.sum.orth/2)^2 )^0.5 +
                            kimura.sum.orth/2) / dpi/2.54, 3)[[1]]
    }else if(unit == "inch"){
      rootlength = round((( kimura.sum.diag^2 + (kimura.sum.diag + kimura.sum.orth/2)^2 )^0.5 +
                            kimura.sum.orth/2) / dpi, 3)[[1]]
    }

    # Validate final result
    if (is.null(rootlength) || length(rootlength) == 0 || !is.numeric(rootlength)) {
      stop("Error calculating root length")
    }

    if (rootlength < 0) {
      warning("Calculated root length is negative, which may indicate an issue with the input image")
    }

    return(rootlength[[1]])

  }, error = function(e) {
    stop("Error in RootLength function: ", e$message)
  })
}

## RootScape Metrics

# input image should be segmented raster

#' RootScapeMetric relies on Landscapemetrics to extract 'Root Scape' Features akin to landscape analysis.
#'
#' @param img segmented raster  (values = 0,1). Consider whether skeletonized raster is appropriate.
#' @param indexD please specify depth. Will only affect the output column = "depth". Useful when used in a loop.
#' @param metrics which ,metrics should be calculated from the available ones in 'landscapemetrics::calculate_lsm()'.
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `2`.
#' @import dplyr
#' @return a bunch of metric values
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' RootScapeObject  = root_scape_metrics(img,indexD = 80, select.layer = 2,  metrics = c("lsm_c_ca"))
root_scape_metrics = function(img, indexD=NA, select.layer = NULL,
                            metrics = c("lsm_c_ca", "lsm_l_ent", "lsm_c_pd", "lsm_c_np", "lsm_c_pland",
                                        "lsm_c_area_mn", "lsm_c_area_cv", "lsm_c_enn_mn", "lsm_c_enn_cv")) {
  tryCatch({
    # Input validation
    if (missing(img)) {
      stop("Image input is required")
    }

    # Validate metrics parameter
    if (!is.character(metrics) || length(metrics) == 0) {
      stop("metrics must be a non-empty character vector")
    }

    # Check if specified metrics are available in landscapemetrics
    available_metrics <- landscapemetrics::list_lsm()$function_name
    invalid_metrics <- metrics[!metrics %in% available_metrics]
    if (length(invalid_metrics) > 0) {
      stop("Invalid metrics specified: ", paste(invalid_metrics, collapse = ", "))
    }

    # Validate select.layer if provided
    if (!is.null(select.layer)) {
      if (!is.numeric(select.layer) || select.layer < 1) {
        stop("select.layer must be a positive integer")
      }
    }

    # Load and validate image
    img <- load_flexible_image(img, select.layer = select.layer,
                               output_format = "spatrast", normalize = FALSE)

    if (is.null(img) || terra::nlyr(img) < 1) {
      stop("Invalid or empty image after loading")
    }

    # Calculate landscape metrics with error handling
    rsm <- tryCatch({
      landscapemetrics::calculate_lsm(img, directions = 8,
                                      neighbourhood = 8,
                                      what = metrics)
    }, error = function(e) {
      stop("Error calculating landscape metrics: ", e$message)
    })

    if (is.null(rsm) || nrow(rsm) == 0) {
      stop("No valid metrics could be calculated")
    }

    # Process results
    t.object <- ifelse(rsm$class == 0, "deletable", "root")
    t.object <- ifelse(is.na(rsm$class), "root", t.object)
    rsm$object <- t.object
    rsm$depth <- indexD

    # Clean and filter results
    rsm <- dplyr::distinct(rsm)
    rsm <- dplyr::filter(rsm, rsm$object != "deletable")
    rsm$id <- NULL
    rsm$class <- NULL
    rsm$level <- NULL
    rsm$layer <- NULL

    return(rsm)

  }, error = function(e) {
    stop("Error in RootScapeMetrics: ", e$message)
  })
}

#' counts all pixels in a segmented image
#'
#' @param img one layer image
#'
#' @return a numeric value
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)[[2]]
#' rootpixel  = count_pixels(img)
count_pixels = function(img) {
  tryCatch({
    # Input validation
    if (missing(img)) {
      stop("Image input is required")
    }

    # Load and validate image
    img <- load_flexible_image(img, output_format = "spatrast", normalize = TRUE)

    if (is.null(img) || terra::nlyr(img) < 1) {
      stop("Invalid or empty image after loading")
    }

    # Calculate pixel sum with validation
    srpx <- terra::global(img, "sum", na.rm = TRUE)

    if (is.null(srpx) || length(srpx) == 0) {
      stop("Error calculating pixel sum")
    }

    if (srpx[[1]] < 0) {
      warning("Negative pixel sum calculated, which may indicate an issue with the input image")
    }

    return(srpx[[1]])

  }, error = function(e) {
    stop("Error in count_pixels: ", e$message)
  })
}

#' Calculate Image Coloration Metrics
#'
#' @param img Three-band raster (RGB) or path to image
#' @param r Red channel weight
#' @param g Green channel weight
#' @param b Blue channel weight
#' @return Data frame of color metrics
#' @export
#'
#' @examples
#' data(rgb_Oulanka2023_Session03_T067)
#' img = terra::rast(rgb_Oulanka2023_Session03_T067)
#' colorvector = tube_coloration(img)
tube_coloration = function(img, r=0.2126, g=0.7152, b=0.0722) {
  tryCatch({
    # Input validation
    if (missing(img)) {
      stop("Image input is required")
    }

    # Validate weight parameters
    if (!is.numeric(r) || !is.numeric(g) || !is.numeric(b) ||
        r < 0 || g < 0 || b < 0 || (r + g + b) != 1) {
      stop("RGB weights must be non-negative numbers that sum to 1")
    }

    # Load and validate image
    img <- load_flexible_image(img, output_format = "spatrast", normalize = FALSE, binarize = FALSE, select.layer = NULL)

    if (is.null(img) || terra::nlyr(img) != 3) {
      stop("Input must be a valid three-band (RGB) image")
    }

    # Extract and validate RGB values
    vr <- terra::values(img[[1]])
    vg <- terra::values(img[[2]])
    vb <- terra::values(img[[3]])

    if (all(is.na(vr)) || all(is.na(vg)) || all(is.na(vb))) {
      stop("One or more color channels contain no valid data")
    }

    # Calculate means with validation
    mean.r <- mean(vr, na.rm=TRUE)
    mean.g <- mean(vg, na.rm=TRUE)
    mean.b <- mean(vb, na.rm=TRUE)

    # Validate color values
    if (any(c(mean.r, mean.g, mean.b) > 255) ||
        any(c(mean.r, mean.g, mean.b) < 0)) {
      warning("Color values outside expected range [0, 255]")
    }

    # Calculate HSL with error handling
    hsl <- tryCatch({
      grDevices::rgb2hsv(r = mean.r, g = mean.g, b = mean.b)
    }, error = function(e) {
      stop("Error converting RGB to HSV: ", e$message)
    })

    # Calculate intensity and color metrics
    intensity <- vr + vg + vb
    if (any(intensity == 0, na.rm = TRUE)) {
      message("Some pixels have zero intensity, which may affect color calculations")
    }

    lum.gray <- vr * r + vg * g + vb * b
    mean.intensity <- round(mean(intensity, na.rm=TRUE), 4)
    mean.lum <- round(mean(lum.gray, na.rm=TRUE), 4)

    # Calculate color channels with validation
    rcc <- round(mean(vr / intensity, na.rm = TRUE), 4)
    gcc <- round(mean(vg / intensity, na.rm = TRUE), 4)
    bcc <- round(mean(vb / intensity, na.rm = TRUE), 4)

    # Create and validate output dataframe
    colordf <- data.frame(
      rcc = rcc, gcc = gcc, bcc = bcc,
      hue = hsl[1], saturation = hsl[2], luminosity = hsl[3],
      red = mean.r, green = mean.g, blue = mean.b
    )

    return(colordf)

  }, error = function(e) {
    stop("Error in Tube.coloration: ", e$message)
  })
}

#' Texture calculation
#'
#' @param img.color Three-band raster or path to image
#' @param grays Number of gray levels
#' @param window Window size for GLCM
#' @param metrics Texture metrics to calculate
#' @return Raster with texture metrics
#' @export
#'
#' @import raster
#'
#' @examples
#' data(rgb_Oulanka2023_Session03_T067)
#' img = raster::brick(rgb_Oulanka2023_Session03_T067)
#' analyze_soil_texture(img, 7, c(9,9), metrics = "second_moment")
analyze_soil_texture = function(img.color, grays = 7, window = c(9,9),
                   metrics = c("variance","second_moment")) {
  tryCatch({
    # Input validation
    if (missing(img.color)) {
      stop("Image input is required")
    }

    # Validate grays parameter
    if (!is.numeric(grays) || grays < 2 || grays > 255) {
      stop("grays must be a numeric value between 2 and 255")
    }

    # Validate window parameter
    if (!is.numeric(window) || length(window) != 2 ||
        any(window < 1) || any(window %% 2 == 0)) {
      stop("window must be a vector of two odd positive integers")
    }

    # Validate metrics parameter
    valid_metrics <- c("mean", "variance", "homogeneity", "contrast",
                       "dissimilarity", "entropy", "second_moment",
                       "correlation")
    invalid_metrics <- metrics[!metrics %in% valid_metrics]
    if (length(invalid_metrics) > 0) {
      stop("Invalid metrics specified: ", paste(invalid_metrics, collapse = ", "))
    }

    # Load and validate image
    img <- load_flexible_image(img.color, output_format = "spatrast",
                               normalize = FALSE)

    if (is.null(img) || terra::nlyr(img) != 3) {
      stop("Input must be a valid three-band (RGB) image")
    }

    # Determine maximum value
    mx <- max(raster::values(img.color), na.rm=TRUE)
    mx <- if(mx > 1) 255 else 1

    # Calculate grayscale image
    img.gray <- (img.color[[1]]*0.21 + img.color[[2]]*0.72 +
                   img.color[[3]]*0.07) / mx

    # Validate grayscale conversion
    if (all(is.na(raster::values(img.gray)))) {
      stop("Grayscale conversion resulted in invalid data")
    }

    # Calculate texture with error handling
    tx.im <- tryCatch({
      glcm::glcm(img.gray,
                 n_grey = grays,
                 statistics = metrics)
    }, error = function(e) {
      stop("Error calculating texture metrics: ", e$message)
    })

    if (is.null(tx.im)) {
      stop("Texture calculation produced no valid results")
    }

    return(tx.im)

  }, error = function(e) {
    stop("Error in texture function: ", e$message)
  })
}


#' Approximate average Root Thickness
#'
#' @param kimuralength length of roots in image section, input unit is cm.
#' @param rootpx amount of rootpx in the image section
#' @param dpi image resolution
#'
#' @return a value in units cm
#' @export
#'
#' @examples root.ticc = root_thickness(kimuralength = 300,rootpx = 9500, dpi = 300)
root_thickness = function(kimuralength, rootpx, dpi=300) {
  tryCatch({
    # Input validation
    if (missing(kimuralength) || missing(rootpx)) {
      stop("Both kimuralength and rootpx are required")
    }

    # Validate numeric inputs
    if (!is.numeric(kimuralength) || !is.numeric(rootpx) || !is.numeric(dpi)) {
      stop("All inputs must be numeric values")
    }

    # Validate positive values
    if (kimuralength <= 0 || rootpx <= 0 || dpi <= 0) {
      stop("All inputs must be positive values")
    }

    # Calculate with validation
    px.length <- kimuralength * (dpi/2.54)
    if (px.length <= 0) {
      stop("Invalid pixel length calculation")
    }

    thiccness <- rootpx / px.length
    if (thiccness <= 0) {
      stop("Invalid thickness calculation")
    }

    px.thicc <- thiccness / (dpi/2.54)

    # Validate final result
    if (px.thicc <= 0.000000001) {
      warning("Calculated thickness is very small, which may indicate measurement issues")
    }

    return(px.thicc)

  }, error = function(e) {
    stop("Error in root.thickness: ", e$message)
  })
}

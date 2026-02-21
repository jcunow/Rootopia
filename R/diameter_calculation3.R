#' Estimate Root Diameters
#'
#' This function estimates root diameters and root volume from an input image using skeletonization and distance transform methods.
#' The input can be a file path, raster, image object, or array, which is converted to a binary image before processing.
#'
#' @param img A character string (file path), `SpatRaster`, `RasterBrick`, `RasterLayer`, `cimg`, `magick-image`, or array.
#'   The input image to process.
#' @param diagnostics Logical. If `TRUE`, enables diagnostic plots and logging. Default is `FALSE`.
#' @param skeleton_method Character. The method to use for skeletonization. Default is `"MAT"`. Will be skipped if skeleton `SpatRaster`is provided.
#' @param skeleton.img A character string (file path), `SpatRaster`, `RasterBrick`, `RasterLayer`, `cimg`, `magick-image`, or array. Uses this object instead of computing it from scratch.
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `2`.
#' @param unit output in pixel 'px', 'inch' or in 'cm'
#' @param dpi scan resolution. Only used if unit = 'cm' or 'inch'
#'
#' @details
#' The function works as follows:
#' - Converts the input image to a binary format (`cimg`).
#' - Applies a distance transform to compute the Euclidean distance for the foreground (root) pixels.
#' - Skeletonizes the binary image to identify root centerlines.
#' - Filters distance values to retain only those corresponding to the skeletonized regions.
#' - Computes diameter statistics, including quantiles, mean, and median diameters.
#'
#' The function supports various input formats and normalizes image values to the range [0, 1] if needed. It uses the `terra` package for raster operations and the `imager` package for image processing.
#'
#' @return A list containing:
#' \describe{
#'   \item{quantiles}{Numeric vector of diameter quantiles (10th to 100th percentile).}
#'   \item{mean_diameter}{Numeric. The mean root diameter.}
#'   \item{median_diameter}{Numeric. The median root diameter.}
#'   \item{diameters}{Numeric vector of all diameter values in the skeletonized regions.}
#'   \item{skeleton_rast}{`SpatRaster`. Binary raster mask of skeletonized regions.}
#'   \item{diameter_rast}{`SpatRaster`. Raster showing diameters in the skeletonized regions.}
#'   \item{distance_map_rast}{`SpatRaster`. Raster showing the distance transform values.}
#'   \item{root_volume}{Numeric. The sum of root volume - assuming cylindrical roots}
#' }
#'
#'@examples
#' # Example usage:
#' data(seg_Oulanka2023_Session01_T067)
#' result <- root_diameter(img = seg_Oulanka2023_Session01_T067,
#'   skeleton_method = "MAT", select.layer = 2, unit = "px",
#'   diagnostics = TRUE)
#'
#' # Access results:
#' print(result$mean_diameter)
#' terra::plot(result$skeleton_rast)
#'
#' @export
root_diameter <- function(img,  skeleton_method = "MAT", skeleton.img = NULL, select.layer = NULL, 
                          diagnostics = FALSE, unit = "cm", dpi = 300) {
  # Input validation and error handling module
  tryCatch({
    # Validate input parameters
    if (missing(img)) {
      stop("Input image is required")
    }


    # Validate select.layer
    if ( (!is.numeric(select.layer) || select.layer < 1) && !is.null(select.layer)) {
      stop("select.layer must be a positive integer or NULL")
    }

    # Validate diagnostics
    if (!is.logical(diagnostics)) {
      stop("diagnostics must be TRUE or FALSE")
    }
    
    # Validate unit parameter
    if (!unit %in% c("px", "cm","inch")) {
      stop("Unit must be either 'px', 'cm', or 'inch'")
    }
    
    # Validate DPI when unit is cm or inch
    if (unit == "cm" | unit == "inch") {
      if (is.null(dpi) || !is.numeric(dpi) || dpi <= 0) {
        stop("Valid positive numeric dpi value is required when unit = 'cm' or 'inch'.")
      }
    }

    # Load and validate image
    tryCatch({
      img <- load_flexible_image(img,
                                        select.layer = select.layer,
                                        output_format = "cimg",
                                        normalize = TRUE, binarize = TRUE)

      # Check if image is empty or all NA
      if (all(is.na(img)) || length(img) == 0) {
        stop("Input image is empty or contains only NA values")
      }

    }, error = function(e) {
      stop(sprintf("Failed to load image: %s", e$message))
    })

    # Distance transform function with error checking
    distance_transform <- function(img) {
      if (!inherits(img, "cimg")) {
        stop("Input to distance_transform must be a cimg object")
      }

      dt <- tryCatch({
        imager::as.cimg(imager::distance_transform(imager::as.cimg(img), value = 0, metric = 2L))
      }, error = function(e) {
        stop(sprintf("Distance transform failed: %s", e$message))
      })

      if (all(dt == 0)) {
        warning("Distance transform produced all zero values - check input image")
      }

      return(dt)
    }

    # Main processing with error checking
    distance_map <- distance_transform(img = img)
    # radius to diameter
    diameters <- distance_map * 2



    # Convert to SpatRast with validation
    Ds <- tryCatch({
      load_flexible_image(diameters, output_format = "SpatRaster", normalize = FALSE, binarize = FALSE, select.layer = NULL)


    }, error = function(e) {
      stop(sprintf("Failed to convert distance map to SpatRaster: %s", e$message))
    })

    IM <- tryCatch({
      load_flexible_image(img, output_format = "SpatRaster",normalize = FALSE, binarize = FALSE, select.layer = NULL)
    }, error = function(e) {
      stop(sprintf("Failed to convert image to SpatRaster: %s", e$message))
    })

    if(is.null(skeleton.img) ){
      # Skeletonization with validation
      IMS <- tryCatch({
        skeleton <- skeletonize_image(IM, methods = skeleton_method,select.layer = NULL)
        if (all(terra::values( skeleton) == 0)) {
          warning("Skeletonization produced empty result - check input image")
        }
        load_flexible_image(skeleton, output_format = "SpatRaster",  normalize = FALSE, binarize = FALSE, select.layer = NULL)
      }, error = function(e) {
        stop(sprintf("Skeletonization failed: %s", e$message))
      })      
    }else {
      IMS = load_flexible_image(skeleton.img, output_format = "SpatRaster",  normalize = FALSE, binarize = TRUE, select.layer = NULL)
    }

    
    # Filter root regions
    DsSKL <- Ds
    terra::ext(IMS) <- terra::ext(DsSKL) 
    DsSKL[[1:dim(Ds)[3]]][IMS == 0] <- NA  
    # DsSKL[IMS == 0] <- NA

     # remove non root distances
     DsSKL[DsSKL == 0] <- NA
     

    # Check if we have any valid measurements
    if (all(is.na(terra::values(DsSKL)))) {
      stop("No valid diameter measurements found after filtering")
    }

    # Create binary skeleton mask
    skl <- DsSKL
    skl[skl > 0] <- 1

    # outpur unit conversion
    if(unit == "cm") {
      DsSKL = DsSKL / (dpi / 2.54)
    }else if(unit == "inch"){
      DsSKL = DsSKL / (dpi)
    }else if(unit == "px"){
      DsSKL = DsSKL 
    }
    

    # Compute statistics with validation
    tryCatch({
      quantile_diameter <- stats::quantile(terra::values(DsSKL),
                                           probs = seq(0, 1, 0.1),
                                           na.rm = TRUE)
      mean_diameter <- terra::global(DsSKL,
                                     fun = function(x) base::mean(x, na.rm = TRUE))$global
      median_diameter <- terra::global(DsSKL,
                                       fun = function(x) stats::median(x, na.rm = TRUE))$global
      diameters <- terra::values(DsSKL, na.rm = TRUE)
      
      root.volume <- sum(diameters**2 * pi, na.rm = TRUE)
  
      # Validate statistics
      if (any(is.na(c(mean_diameter, median_diameter)))) {
        stop("Failed to compute diameter statistics")
      }

    }, error = function(e) {
      stop(sprintf("Statistical calculations failed: %s", e$message))
    })

    # Optional diagnostics output
    if (diagnostics) {
      message("Processing complete. Summary statistics:")
      message(sprintf("Mean diameter: %.2f", mean_diameter))
      message(sprintf("Median diameter: %.2f", median_diameter))
      message(sprintf("Number of valid measurements: %d", sum(!is.na(diameters))))
    }

    # Return results
    return(list(
      quantiles = quantile_diameter,
      mean_diameter = mean_diameter,
      median_diameter = median_diameter,
      diameters = diameters,
      skeleton_rast = skl,
      diameter_rast = DsSKL,
      distance_map_rast = Ds,
      root_volume = root.volume
    ))

  }, error = function(e) {
    # Main error handler
    stop(sprintf("Error in root.diameters: %s", e$message))
  }, warning = function(w) {
    # Warning handler
    warning(sprintf("Warning in root.diameters: %s", w$message))
  })
}





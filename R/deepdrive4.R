

#' Assess Root Growth Direction Relative to Depth Gradient
#'
#' This function calculates the proportion of root pixels that grow in the direction
#' of steepest local depth descent — a metric for evaluating how efficiently roots
#' explore the vertical soil profile.
#'
#' It uses a depth raster (`DepthMap`) and a root direction raster (`AngleMap`)
#' or estimates one from a binary root presence map (`RootMap`). Directions are
#' compared using a D8 neighborhood scheme.
#'
#' @param DepthMap A SpatRast object representing local depth (e.g., distance from surface or tube wall).
#' @param RootMap Optional. A binary SpatRast indicating root presence. Used to infer `AngleMap` if not provided.
#' @param AngleMap Optional. A SpatRast of root angles in D8 format (0, 45, ..., 315). If missing, inferred from `RootMap` and `DepthMap`.
#' @param select.layerRM Integer. Which layer to use from `RootMap` if it has multiple bands. Default is `2`.
#' @param select.layerDM Integer. Which layer to use from `DepthMap`. Default is `NULL`.
#' @param select.layerAM Integer. Which layer to use from `AngleMap`. Default is `NULL`.
#' @param return Character. `"value"` (default) returns a single numeric proportion. `"all"` returns a list with spatial outputs for visualization.
#'
#' @return
#' - If `return = "value"`: A numeric value between 0 and 1.
#'   - `1` means all root pixels align with the local steepest downward direction.
#'   - `0` means none of the root pixels align with local depth gradients.
#'
#' - If `return = "all"`: A named list with:
#'   - `deep_drive`: numeric proportion (same as above)
#'   - `angle_map`: the root direction map used (SpatRast)
#'   - `optimal_angle_map`: the steepest downward direction from `DepthMap` (SpatRast)
#'   - `aligned_roots`: binary SpatRast showing root pixels aligned with the steepest local depth descent
#'
#'
#' @examples
#' data(skl_Oulanka2023_Session01_T067)
#' im <- ceiling(terra::rast(skl_Oulanka2023_Session01_T067) / 255)
#' DepthMap <- terra::t(create_depthmap(im, center.offset = 0, tube.thicc = 3.5))
#'
#' # Just the deep drive score
#' deep_drive(DepthMap = DepthMap, RootMap = im, select.layerRM = 2)
#'
#' # Get spatial outputs too
#' res <- deep_drive(DepthMap = DepthMap, RootMap = im, select.layerRM = 2, return = "all")
#' plot(res$aligned_roots)
#'
#' @export
deep_drive <- function(DepthMap,
                       RootMap = NULL,
                       AngleMap = NULL,
                       select.layerRM = NULL,
                       select.layerDM = NULL,
                       select.layerAM = NULL,
                       return = c("value", "all")) {
  
  return <- match.arg(return)
  
  # Validation module
  tryCatch({
    # Check required inputs
    if (is.null(DepthMap)) {
      stop("DepthMap is required")
    }
    
    # Input type validation
    if (!is.null(select.layerRM) && !is.numeric(select.layerRM)) {
      stop("select.layerRM must be numeric")
    }
    if (!is.null(select.layerDM) && !is.numeric(select.layerDM)) {
      stop("select.layerDM must be numeric")
    }
    if (!is.null(select.layerAM) && !is.numeric(select.layerAM)) {
      stop("select.layerAM must be numeric")
    }
    
    # Layer selection validation
    if (!is.null(select.layerRM) && select.layerRM < 1) {
      stop("select.layerRM must be positive")
    }
    if (!is.null(select.layerDM) && select.layerDM < 1) {
      stop("select.layerDM must be positive")
    }
    if (!is.null(select.layerAM) && select.layerAM < 1) {
      stop("select.layerAM must be positive")
    }
    
    # Check if RootMap is provided when AngleMap is NULL
    if (is.null(AngleMap) && is.null(RootMap)) {
      stop("Either AngleMap or RootMap must be provided")
    }
    
    # Load and validate DepthMap
    DepthMap <- load_flexible_image(DepthMap, select.layer=select.layerDM,
                                    output_format="spatrast", normalize=FALSE, binarize = FALSE)
    if (is.null(DepthMap)) {
      stop("Failed to load DepthMap")
    }
    
    # Ensure DepthMap has valid values
    if (terra::global(DepthMap, "isNA", na.rm=TRUE)[1] == terra::ncell(DepthMap)) {
      stop("DepthMap contains only NA values")
    }
    
  }, error = function(e) {
    stop(paste("Validation error:", e$message))
  })
  
  # Main function logic with additional error handling
  tryCatch({
    if(is.null(AngleMap)){
      RootMap <- load_flexible_image(RootMap, select.layer=select.layerRM,
                                     output_format="spatrast", normalize=T, binarize = T)
      if (is.null(RootMap)) {
        stop("Failed to load RootMap")
      }
      
      # align extents
      terra::ext(RootMap) <- terra::ext(DepthMap)
      
      # ensure positive Depth increments
      DepthMap_for_angle = abs(DepthMap)
      
      # the D8 flowdir algorithm needs decreasing values
      dem = -DepthMap_for_angle
      dem[RootMap != 1] <- NA
      dem = terra::t(terra::flip(dem))
      AngleMap = terra::terrain(dem, v="flowdir")
      AngleMap = terra::subst(AngleMap,
                              from=c(0,1,2,4,8,16,32,64,128),
                              to=c(NA,90,135,180,225,270,315,0,45))
    } else {
      AngleMap <- load_flexible_image(AngleMap, select.layer=select.layerAM,
                                      output_format="spatrast", normalize=FALSE)
      if (is.null(AngleMap)) {
        stop("Failed to load AngleMap")
      }
    }
    
    # align orientation with AngleMap
    DepthMap = terra::t(terra::flip(DepthMap))
    
    # which pixel to go to reach the next deepest pixel in 8px neighbourhood
    mxslope <- terra::focal(DepthMap, w=c(3,3), fun=function(x)x[c(1:4,6:9)] - x[5])
    # correct for diagonal length 
    mxslope[[c(1,3,6,8)]] = mxslope[[c(1,3,6,8)]] / sqrt(2)
    # finding the steepest descend
    gg = terra::which.max(mxslope)
    # give meaningful labels
    gg = terra::subst(gg, from = c(1:8), to = c(315,0,45,270,90,225,180,135))
    
    # Simplified and more efficient calculation approach
    # Count aligned root pixels directly
    aligned_mask <- AngleMap == gg
    depthdrivepx <- sum(aligned_mask[], na.rm = TRUE)
    
    # Count total root pixels
    zonepx <- sum(!is.na(AngleMap[]), na.rm = TRUE)
    
    # Add validation for final calculation
    if (zonepx == 0) {
      warning("No valid pixels found for calculation")
      return(if (return == "all") list(deep_drive = NA, angle_map = AngleMap, optimal_angle_map = gg, aligned_roots = NULL, misaligned_roots = NULL) else NA)
    }
    
    deep.drive.v = depthdrivepx / zonepx
    
    # Prepare optional output layers
    if (return == "all") {
      aligned_roots <- aligned_mask
      
      return(list(
        deep_drive = deep.drive.v,
        angle_map = AngleMap,
        optimal_angle_map = gg,
        aligned_roots = terra::mask(aligned_roots, AngleMap),
      ))
    }
    
    return(deep.drive.v)
    
  }, error = function(e) {
    stop(paste("Error in main function:", e$message))
  }, warning = function(w) {
    warning(paste("Warning in main function:", w$message))
  })
}

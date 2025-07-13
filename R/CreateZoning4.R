
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
#' depthmap = create_depthmap(img,mask,start.soil = 2.9, select.layer = 2 )
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



#' Zone Image Data by Depth and/or Rotation Slices with Optional Spatial Cropping
#'
#' This function extracts and zones image data based on depth bins and/or rotation slices.
#' It supports three modes of operation:
#' \describe{
#'   \item{\code{"depth"}}{Zones image pixels according to binned depth values using a provided depth map and selected depth indices.}
#'   \item{\code{"rotation"}}{Zones the image by slicing it according to specified rotation slices along the x-axis (rows).}
#'   \item{\code{"both"}}{Applies depth zoning first, then applies rotation slicing on the depth-zoned image.}
#' }
#'
#' The function supports flexible cropping of the input image using a specified spatial extent,
#' which is applied before any zoning operations. When depth mode is used, the function warns
#' if cropping would remove pixels containing desired depth values.
#'
#' @param img A \code{terra::SpatRaster} object representing the input image. Can be multi-layer.
#' @param depth_map A \code{terra::SpatRaster} object with binned depth values corresponding spatially to \code{img}. Required if \code{mode} is \code{"depth"} or \code{"both"}. Must have compatible spatial properties with \code{img}.
#' @param depth Numeric vector or range specifying depth values to select from \code{depth_map}. Supports single values, sequences (e.g. \code{3:6}), or arbitrary numeric vectors. Required if \code{mode} is \code{"depth"} or \code{"both"}.
#' @param select_layer Integer scalar. Selects which layer to extract from \code{img} if multi-layer. Use \code{NULL} to keep all layers.
#' @param crop_extent Numeric vector of length 4 in format \code{c(xmin, xmax, ymin, ymax)} for spatial cropping applied before zoning operations. Use \code{NULL} for no cropping.
#' @param mode Character string specifying the zoning mode. One of \code{"depth"}, \code{"rotation"}, or \code{"both"}. Default is \code{"rotation"}.
#' @param rotation_slices Numeric vector of length 2 specifying start and end slice indices for rotation zoning (e.g., \code{c(2, 4)}). Values must be between 1 and \code{rotation_total_slices}. Required if \code{mode} is \code{"rotation"} or \code{"both"}.
#' @param rotation_total_slices Integer scalar specifying the total number of conceptual slices along the x-axis for rotation zoning. Required if \code{mode} is \code{"rotation"} or \code{"both"}.
#'
#' @return A \code{terra::SpatRaster} object containing the zoned and optionally cropped image data. Pixels not matching the zoning criteria are set to \code{NA}.
#'
#' @details
#' The function processes operations in the following order:
#' \enumerate{
#'   \item \strong{Spatial cropping}: If \code{crop_extent} is provided, both \code{img} and \code{depth_map} are cropped first.
#'   \item \strong{Depth zoning}: If mode includes "depth", pixels are masked based on depth values in \code{depth_map}.
#'   \item \strong{Rotation zoning}: If mode includes "rotation", the image is sliced along the x-axis (rows) according to \code{rotation_slices}.
#'   \item \strong{Layer selection}: If \code{select_layer} is specified, only that layer is retained.
#' }
#'
#' \strong{Depth Processing}: 
#' \itemize{
#'   \item If \code{depth} is a single value, selects the closest matching depth in \code{depth_map}
#'   \item If \code{depth} is a sequence (consecutive values), selects all depths within the range
#'   \item If \code{depth} is an arbitrary vector, selects the closest unique depths for each value
#' }
#'
#' \strong{Rotation Processing}:
#' \itemize{
#'   \item Divides the image rows into \code{rotation_total_slices} conceptual slices
#'   \item Extracts rows corresponding to slices \code{rotation_slices[1]} through \code{rotation_slices[2]}
#'   \item Useful for analyzing angular segments in circular or rotational patterns
#' }
#'
#' \strong{Performance Optimization}:
#' \itemize{
#'   \item Uses vectorized operations with terra::ifel() for fast depth filtering
#'   \item Direct logical operations for creating depth masks
#'   \item Eliminates intermediate raster creation for better performance
#' }
#'
#' \strong{Warnings}: The function issues warnings when:
#' \itemize{
#'   \item Cropping would remove pixels containing desired depth values
#'   \item The resulting image contains only NA values
#' }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' depth_map = terra::t(create_depthmap(img))
#' 
#' # Depth zoning example - select depths 3 through 6
#' zone_img <- zone_image(
#'   img = img,
#'   mode = "depth",
#'   depth_map = depth_map,
#'   depth = 3:6
#' )
#'
#' # Rotation zoning example - extract slices 2-4 out of 8 total
#' zone_img <- zone_image(
#'   img = img,
#'   mode = "rotation",
#'   rotation_slices = c(2, 4),
#'   rotation_total_slices = 8
#' )
#'
#' # Combined zoning with cropping
#' zone_img <- zone_image(
#'   img = img,
#'   mode = "both",
#'   depth_map = depth_map,
#'   depth = c(4.5, 5.2, 6.0),
#'   rotation_slices = c(1, 3),
#'   rotation_total_slices = 5,
#'   crop_extent = c(100, 500, 200, 400)
#' )
#'
#' # Select only the first layer and crop spatially
#' zone_img <- zone_image(
#'   img = img,
#'   mode = "depth",
#'   depth_map = depth_map,
#'   depth = 5,
#'   select_layer = 1,
#'   crop_extent = c(0, 1000, 0, 1000)
#' )
#'
#' # Mixed depth selection: range plus individual values
#' zone_img <- zone_image(
#'   img = img,
#'   mode = "depth",
#'   depth_map = depth_map,
#'   depth = c(3:8, 13, 15.5)  # Range 3-8 plus closest to 13 and 15.5
#' )
#' }
#'
#'
#' @export
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
    
    # Check spatial compatibility
    if (!terra::compareGeom(img, depth_map, stopOnError = FALSE)) {
      stop("img and depth_map must have compatible spatial properties (extent, resolution, CRS).")
    }
  }
  
  # Fast vectorized depth selection function
  depth_range_from_input <- function(depth, available) {
    available <- sort(unique(available[!is.na(available)]))
    if (length(available) == 0) {
      return(numeric(0))
    }
    
    # Handle different depth input patterns
    if (length(depth) == 1) {
      # Single value - find closest available depth
      matched <- available[which.min(abs(available - depth))]
    } else {
      # Multiple values - check if it's a consecutive sequence
      sorted_depth <- sort(depth)
      diffs <- diff(sorted_depth)
      is_consecutive <- length(unique(round(diffs, 10))) == 1 && all(diffs > 0)
      
      if (is_consecutive && length(depth) > 2) {
        # Treat as range: everything >= min and <= max
        range_min <- min(depth)
        range_max <- max(depth)
        matched <- available[available >= range_min & available <= range_max]
      } else {
        # Treat as discrete values - find closest available depth for each
        matched <- sapply(depth, function(x) {
          available[which.min(abs(available - x))]
        })
        matched <- sort(unique(matched[!is.na(matched)]))
      }
    }
    
    return(matched)
  }
  
  # Store original depth_map for potential operations
  original_depth_map <- depth_map
  
  # --- APPLY CROPPING FIRST FOR ALL MODES ---
  if (!is.null(crop_extent)) {
    if (length(crop_extent) != 4) {
      stop("crop_extent must be c(xmin, xmax, ymin, ymax)")
    }
    
    # Validate crop_extent values
    if (crop_extent[1] >= crop_extent[2] || crop_extent[3] >= crop_extent[4]) {
      stop("Invalid crop_extent: xmin must be < xmax and ymin must be < ymax")
    }
    
    # Build terra extent object from crop_extent vector
    ext <- terra::ext(crop_extent[1], crop_extent[2], crop_extent[3], crop_extent[4])
    
    # Check if cropping will affect depth analysis
    if (mode %in% c("depth", "both") && !is.null(depth_map)) {
      # Get available depths before cropping
      depth_vals_vec_original <- terra::values(original_depth_map, mat = FALSE)
      available_original <- sort(unique(depth_vals_vec_original[!is.na(depth_vals_vec_original)]))
      
      # Get target depths
      target_depths <- depth_range_from_input(depth, available_original)
      
      # Crop depth_map to see what depths remain
      depth_map_cropped <- terra::crop(depth_map, ext)
      depth_vals_vec_cropped <- terra::values(depth_map_cropped, mat = FALSE)
      available_cropped <- sort(unique(depth_vals_vec_cropped[!is.na(depth_vals_vec_cropped)]))
      
      # Check if any target depths are lost due to cropping
      lost_depths <- setdiff(target_depths, available_cropped)
      if (length(lost_depths) > 0) {
        warning(sprintf("Cropping will remove pixels with desired depth values: %s. Consider adjusting crop_extent or depth parameters.", 
                        paste(lost_depths, collapse = ", ")))
      }
      
      # Update depth_map to cropped version
      depth_map <- depth_map_cropped
    }
    
    # Crop img
    img <- terra::crop(img, ext)
  }
  
  # --- Mode 1: Depth zoning (FIXED VERSION) ---
  if (mode %in% c("depth", "both")) {
    # Get available depth values from (potentially cropped) depth_map
    depth_vals_vec <- terra::values(depth_map, mat = FALSE)
    available <- sort(unique(depth_vals_vec[!is.na(depth_vals_vec)]))
    depth_vals <- depth_range_from_input(depth, available)
    
    if (length(depth_vals) == 0) {
      stop("No matching depths found in depth_map for the provided depth argument.")
    }
    
    # FIXED: Use logical combination for multiple depths
    if (length(depth_vals) == 1) {
      mask_condition <- depth_map == depth_vals[1]
    } else {
      # For multiple depths, combine with OR operations
      mask_condition <- depth_map == depth_vals[1]
      for (i in 2:length(depth_vals)) {
        mask_condition <- mask_condition | (depth_map == depth_vals[i])
      }
    }
    img <- terra::ifel(mask_condition, img, NA)
  }
  
  # --- Mode 2: Rotation zoning ---
  if (mode %in% c("rotation", "both")) {
    # Get current image dimensions in pixels
    img_dims <- dim(img)  # c(nrow, ncol, nlyr)
    total_rows <- img_dims[1]
    
    # Calculate slice boundaries along x-axis (rows)
    slice_start <- floor((rotation_slices[1] - 1) * total_rows / rotation_total_slices) + 1
    slice_end <- ceiling(rotation_slices[2] * total_rows / rotation_total_slices)
    slice_start <- max(1, slice_start)
    slice_end <- min(total_rows, slice_end)
    
    # Create extent for x-axis slicing
    x_res <- terra::res(img)[1]
    img_extent <- terra::ext(img)
    
    x_min_coord <- img_extent[1] + (slice_start - 1) * x_res
    x_max_coord <- img_extent[1] + slice_end * x_res
    
    rotation_crop_extent <- terra::ext(x_min_coord, x_max_coord,
                                       img_extent[3], img_extent[4])
    
    # Apply rotation slice cropping
    img <- terra::crop(img, rotation_crop_extent)
  }
  
  # Select specific layer if requested
  if (!is.null(select_layer)) {
    if (!is.numeric(select_layer) || length(select_layer) != 1 || select_layer < 1) {
      stop("select_layer must be a single positive integer.")
    }
    if (terra::nlyr(img) < select_layer) {
      stop("select_layer exceeds number of layers in img.")
    }
    img <- terra::subset(img, select_layer)
  }
  
  # Check if result is empty (all NA values) and warn user
  if (all(is.na(terra::values(img, mat = FALSE)))) {
    warning("Resulting image contains only NA values. Check your depth, rotation, or crop parameters.")
  }
  
  return(img)
}
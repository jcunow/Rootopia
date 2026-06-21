####### create the depthmap

#' Create A Phase-Shifted, Tilt-Amplitude Sine Depth Map
#'
#' This function generates a depth map for minirhizotron images, accounting for tube geometry
#' and insertion angle.
#'
#' @param img Input image (accepts terra SpatRaster, matrix, array, or file path). For multi-band
#'           images, specify band_index parameter
#' @param mask Raster mask indicating foreign objects (1 = mask, 0 or NA = keep)
#' @param sinoid Logical; if TRUE, accounts for tube curvature in depth calculation
#' @param tube.thicc Numeric; diameter of minirhizotron tube in cm
#' @param tilt Numeric; minirhizotron tube insertion angle in degrees (typically 30-45)
#' @param dpi Numeric; image resolution in dots per inch
#' @param start.soil Numeric; soil surface boundary in cm (0 = surface)
#' @param center.offset Numeric; rotational center offset (0 = centered, 1 = edge)
#' @param progress Message; indicates how mny rows have been processed  
#'
#' @return terra raster object containing the depth map
#' @export
#'
#' @author Johannes Cunow, Robert Weigel
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' seg_Oulanka2023_Session01_T067 = terra::rast(seg_Oulanka2023_Session01_T067)
#' img = seg_Oulanka2023_Session01_T067
#' mask = seg_Oulanka2023_Session01_T067[[1]] - seg_Oulanka2023_Session01_T067[[2]]
#' mask[mask == 255] <- NA
#' map = create_depthmap(img,mask,start.soil = 0.1,
#'   sinoid = TRUE,
#'   tube.thicc = 7,
#'   tilt = 45,
#'   dpi = 300,
#'   center.offset = 0.1 )
create_depthmap = function(img, mask = NULL, sinoid = TRUE,
                           tube.thicc = 7, tilt = 45, dpi = 300,
                           start.soil = 0, center.offset = 0.5, progress = FALSE) {

  # Input validation module
  tryCatch({
    # Validate numeric parameters
    if (!is.numeric(tube.thicc) || tube.thicc <= 0)
      stop("tube.thicc must be a positive number")
    if (!is.numeric(tilt) || tilt <= 0 || tilt >= 90)
      stop("tilt must be between 0 and 90 degrees")
    if (!is.numeric(dpi) || dpi <= 0)
      stop("dpi must be a positive number")
    if (!is.numeric(start.soil))
      stop("start.soil must be a numeric value")
    if (!is.numeric(center.offset) || center.offset < 0 || center.offset > 1)
      stop("center.offset must be between 0 and 1")

    # Validate logical parameters
    if (!is.logical(sinoid))
      stop("sinoid must be TRUE or FALSE")

    # Image validation
    if (is.null(img))
      stop("Input image cannot be NULL")

    # all layers have the same dimensions, we select layer 1 to handle multi-layer and single-layer raster
    select.layer = 1
    
    # Try loading the image with error handling
    img <- tryCatch({
      load_flexible_image(img, select.layer = select.layer,
                          output_format = "spatrast", scale = "none")
    }, error = function(e) {
      stop(paste("Failed to load image:", e$message))
    })

    # Validate image dimensions
    if (any(dim(img)[1:2] <= 1))
      stop("Image dimensions must be greater than 1x1")

    # Mask validation and processing
    if (length(mask) == 0) {
      mask = img
      terra::values(mask) <- 0
    } else {
      # Check if mask dimensions match image
      if (!all(dim(mask)[1:2] == dim(img)[1:2]))
        stop("Mask dimensions must match image dimensions")
    }


    # Core computation with additional safety checks
    radiant = pi/180
    tilt.factor = sin((180-tilt) * radiant)
    tube.thicc.tilted = round(tube.thicc * tilt.factor, 3)
    px.to.cm.h = (2.54/dpi)

    target.col = dim(img)[1]
    target.row = dim(img)[2]

    if (sinoid) {
      df1 = seq(0*pi, 2*pi, length =  (round(tube.thicc*pi*dpi/2.54, 0)-1))
      if (length(df1) < 2)
        stop("Insufficient points for sine wave generation. Check tube.thicc and dpi values")

      df00 = (cos(df1+(pi*(center.offset))))*(tube.thicc.tilted/2) + (tube.thicc.tilted/2)

      # Safely handle vector subsetting
      if (target.col > length(df00))
        stop("Calculated sine wave shorter than image width. Wrong Tube Diameter")
      df11 = df00[1:target.col]
    } else {
      df11 = rep(0, target.col)
    }

    # Safe array allocation
    df = try(array(dim = c(target.row, target.col)))
    if (inherits(df, "try-error"))
      stop("Failed to allocate memory for depth map array")

    # Row processing with progress monitoring
    for (ii in 1:target.row) {
      df[ii,] = df11 + (ii*px.to.cm.h * tilt.factor)
      if (progress && ii %% 100 == 0)
        message(sprintf("Processing row %d of %d", ii, target.row))
    }

    df.depthmap = df - (start.soil)
    masked.depthmap = terra::rast(df.depthmap)

    # Final mask application with error checking
    terra::values(masked.depthmap)[terra::values(mask) == 1] <- NA

    return(masked.depthmap)

  }, error = function(e) {
    # Provide informative error message
    stop(paste("Error in create.depthmap:", e$message))
  })
}

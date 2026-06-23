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
#' @param tube_thicc Numeric; diameter of minirhizotron tube in cm
#' @param tilt Numeric; minirhizotron tube insertion angle in degrees (typically 30-45)
#' @param dpi Numeric; image resolution in dots per inch
#' @param start_soil Numeric; soil surface boundary in cm (0 = surface)
#' @param center_offset Numeric; rotational center offset (0 = centered, 1 = edge)
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
#' map = create_depthmap(img,mask,start_soil = 0.1,
#'   sinoid = TRUE,
#'   tube_thicc = 7,
#'   tilt = 45,
#'   dpi = 300,
#'   center_offset = 0.1 )
create_depthmap = function(img, mask = NULL, sinoid = TRUE,
                           tube_thicc = 7, tilt = 45, dpi = 300,
                           start_soil = 0, center_offset = 0.5, progress = FALSE) {

  # Input validation module
  tryCatch({
    # Validate numeric parameters
    if (!is.numeric(tube_thicc) || tube_thicc <= 0)
      stop("tube_thicc must be a positive number")
    if (!is.numeric(tilt) || tilt <= 0 || tilt >= 90)
      stop("tilt must be between 0 and 90 degrees")
    if (!is.numeric(dpi) || dpi <= 0)
      stop("dpi must be a positive number")
    if (!is.numeric(start_soil))
      stop("start_soil must be a numeric value")
    if (!is.numeric(center_offset) || center_offset < 0 || center_offset > 1)
      stop("center_offset must be between 0 and 1")

    # Validate logical parameters
    if (!is.logical(sinoid))
      stop("sinoid must be TRUE or FALSE")

    # Image validation
    if (is.null(img))
      stop("Input image cannot be NULL")

    # all layers have the same dimensions, we select layer 1 to handle multi-layer and single-layer raster
    select_layer = 1
    
    # Try loading the image with error handling
    img <- tryCatch({
      load_flexible_image(img, select_layer = select_layer,
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
    tube_thicc_tilted = round(tube_thicc * tilt.factor, 3)
    px.to.cm.h = (2.54/dpi)

    # Image geometry. Depth increases along the image WIDTH (columns); the
    # sinusoidal tube-curvature offset varies along the HEIGHT (rows). The map
    # is built directly in `img` orientation so the result is aligned with
    # `img` (same extent/resolution/CRS) -- no transpose/flip is needed
    # downstream.
    nr = dim(img)[1]   # height -> tube circumference (sinusoid axis)
    nc = dim(img)[2]   # width  -> depth axis

    # Sinusoidal curvature offset, one value per circumference pixel (row).
    if (sinoid) {
      df1 = seq(0*pi, 2*pi, length =  (round(tube_thicc*pi*dpi/2.54, 0)-1))
      if (length(df1) < 2)
        stop("Insufficient points for sine wave generation. Check tube_thicc and dpi values")

      df00 = (cos(df1+(pi*(center_offset))))*(tube_thicc_tilted/2) + (tube_thicc_tilted/2)

      # Safely handle vector subsetting
      if (nr > length(df00))
        stop("Calculated sine wave shorter than image height. Wrong Tube Diameter")
      circ = df00[1:nr]
    } else {
      circ = rep(0, nr)
    }

    # Depth term per column (deeper to the right); curvature term per row.
    # M[row, col] = curvature(row) + depth(col) - start_soil, in img orientation.
    depth_term = seq_len(nc) * px.to.cm.h * tilt.factor
    M = outer(rev(circ), depth_term, FUN = "+") - start_soil   # dims [nr, nc]

    # Write onto a raster that copies img's geometry so extent/res/CRS match.
    masked.depthmap = terra::rast(img[[1]])
    terra::values(masked.depthmap) = as.vector(t(M))           # row-major for terra

    # Mask foreign objects -- now correctly aligned cell-for-cell with img.
    mv = as.vector(terra::values(mask))
    terra::values(masked.depthmap)[which(mv == 1)] <- NA

    return(masked.depthmap)

  }, error = function(e) {
    # Provide informative error message
    stop(paste("Error in create.depthmap:", e$message))
  })
}

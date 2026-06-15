

#' Get list of supported conversion & input formats as a string
#' @keywords internal
supported_formats_string <- function(type = "Input") {
  if(type == "Input"){
    paste(
      "Supported formats:",
      "- File path (.jpg, .jpeg, .png, .tif, .tiff, .bmp)",
      "- SpatRaster",
      "- RasterLayer",
      "- RasterBrick",
      "- cimg",
      "- magick-image",
      "- matrix",
      "- array",
      sep = "\n"
    )
  }
  else if(type == "Output"){
    paste(
      "Supported formats:",
      "- SpatRaster",
      "- RasterLayer",
      "- RasterBrick",
      "- cimg",
      "- magick-image",
      "- matrix",
      "- array",
      sep = "\n"
    )
  }
}

#' Validate conversion parameters
#' @keywords internal
validate_conversion_params <- function(input, normalize, select.layer, binarize) {
  if (missing(input)) stop("Input is required")
  if (!is.logical(normalize)) stop("normalize must be TRUE or FALSE")
  if (!is.logical(binarize)) stop("binarize must be TRUE or FALSE")

  # Validate file path if input is a string
  if (is.character(input)) {
    if (length(input) != 1) stop("File path must be a single string")
    if (!file.exists(input)) stop(sprintf("File not found: %s", input))
    ext <- tolower(tools::file_ext(input))
    valid_extensions <- c("jpg", "jpeg", "png", "tif", "tiff", "bmp")
    if (!ext %in% valid_extensions) {
      stop(sprintf("Unsupported file extension: %s\nSupported extensions: %s", ext, paste(valid_extensions, collapse = ", ")))
    }
  }
}

#' Normalize or binarize the array or raster
#' @keywords internal
normalize_array <- function(arr, normalize, binarize) {
  if (normalize | binarize) {
    max_val <- max(as.vector(arr), na.rm = TRUE)
    if (is.na(max_val)) warning("Cannot normalize: all values are NA")
    if (max_val > 255) warning("Maximum value exceeds 255, normalization might be incorrect")
    if (binarize) arr <- ceiling(arr / max_val) else arr <- arr / 255
  }
  return(arr)
}

#' Load an image flexibly from file or convert from memory
#' @param input File path or image object
#' @param output_format Character, one out of "cimg", "spatrast", "matrix", "array", "brick", "raster", "spatrast", "magick-image". Other spellings are accepted.
#' @param select.layer Numeric, which layer to select if input has multiple layers
#' @param normalize Logical, whether to normalize values to 0-1 range if they're in 0-255
#' @param binarize Logical, whether the output is strictly 0 and 1. Overwrites normalize
#' @keywords internal
load_flexible_image <- function(input, output_format = "cimg", normalize = TRUE, select.layer = NULL, binarize = FALSE) {
  validate_conversion_params(input, normalize, select.layer, binarize)


  if (inherits(input, output_format) && is.null(select.layer)) {
    return(input)   # already in the requested form
  }

  arr <- NULL

  # Handle file input (image or raster)
  if (is.character(input) && file.exists(input)) {
    ext <- tolower(tools::file_ext(input))
    if (ext %in% c("tif", "tiff")) {
      arr <- as.array(tiff::readTIFF(input))  # 3D array
    } else {
      arr <- as.array(imager::load.image(input))  # 4D array
    }
  }
  else if (inherits(input,"matrix")) {
    arr <- array(input, dim = c(dim(input)))  # Convert matrix to 3D array
  }
  else if (inherits(input,"array")) {
    arr <- input
  }
  else if (inherits(input, c("RasterLayer", "RasterBrick", "SpatRaster"))) {
    arr <- terra::as.array(input)
  }
  else if (inherits(input, c("cimg"))) {
    arr <- as.array(input)[,,1,]
  }
  else if (inherits(input, c("magick-image"))) {
    arr <- imager::magick2cimg(input)
    arr <- as.array(arr)[,,1,]
  }
  else if(!inherits(input, c("RasterLayer", "RasterBrick", "SpatRaster","magick-image", "cimg","array","matrix","character"))) {
    stop("Unsupported input format.\n", supported_formats_string(type = "Input"))
  }


  dims <- dim(arr)

  # special case #1 for RasterLayer
  if (output_format == "RasterLayer" && is.null(select.layer)) {
    if(length(dim(arr)) > 2) {
      warning("select.layer must be specified if output_format is 'RasterLayer'. Default to layer 1." )
      select.layer = 1
    }}

   # special case #2 for RasterLayer to RasterBrick conversion
   if(inherits(input,"RasterLayer") && output_format == "RasterBrick" ) {
     # not possible, will force output format (raster::brick does it by default)
     output_format  <- "RasterLayer"
     warning("RasterLayer to RasterBrick not possible. Changed output_format to RasterLayer")
     select.layer = 1
   }

  # Handle layer selection for 3D or 4D arrays (select a specific layer/channel)
  if (!is.null(select.layer) && length(dims) > 2) {
    arr <- arr[,,select.layer]  # Select layer
  }

  # Normalize or binarize the data
  arr <- normalize_array(arr, normalize, binarize)



  ## Return based on desired output format
  # Overwrite 3D matrix choice
  if((output_format == "matrix" || 
      output_format == "Matrix" ||
      output_format == "MATRIX") && dims[3] > 1 ){
    warning("You cannot convert a 3D image to 2D matrix. An array is returned instead. Consider specifying 'select.layer' if you want to return a matrix.")
    output_format = "array"
  }
  
  # flexible output  conversion
  if (output_format == "SpatRaster" ||
      output_format == "spatraster" ||
      output_format == "spatrast" ||
      output_format == "SpatRast"||
      output_format == "SPATRASTER" ||
      output_format == "SPATRAST" ) {
    r <- terra::rast(arr)
    # Register RGB(A) channel metadata so terra::plotRGB() works on 3/4-band
    # outputs without a separate terra::RGB(r) <- ... call by the caller.
    if (terra::nlyr(r) %in% c(3L, 4L)) terra::RGB(r) <- seq_len(terra::nlyr(r))
    return(r)
  }
  else if (output_format == "RasterBrick" ||
           output_format == "rasterbrick" ||
           output_format == "raster" ||
           output_format == "brick" ||
           output_format == "BRICK" ||
           output_format == "RASTERBRICK") {
    return(raster::brick(terra::rast(arr)))
  }
  else if (output_format == "RasterLayer" ||
           output_format == "rasterlayer" ||
           output_format == "RASTERLAYER" ||
           output_format == "layer" ||
           output_format == "LAYER") {
    return(raster::raster(arr))
  }
  
    
  else if (output_format == "array" ||
           output_format == "Array" ||
           output_format == "ARRAY") {
    array(arr,dim = dims)
  }
  else if (output_format == "matrix" ||
           output_format == "Matrix" ||
           output_format == "MATRIX") {
    return(matrix(arr, dims))
  }
  else if (output_format == "cimg" ||
           output_format == "Cimg" ||
           output_format == "CIMG") {
    return(suppressWarnings(imager::as.cimg(arr)))
  }
  else if(output_format == "magick-image" ||
          output_format == "magick" ||
          output_format == "Magick" ||
          output_format == "Magick-Image" ||
          output_format == "MAGICK" ||
          output_format == "MAGICK-IMAGE" ){
    return(imager::cimg2magick(suppressWarnings(imager::as.cimg(arr))))
  }
  else {
    stop("Unsupported output format.\n", supported_formats_string(type = "Output"))
  }
}

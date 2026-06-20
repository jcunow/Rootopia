

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

# The set of accepted `scale` values, defined once so the signature default,
# match.arg(), and the deprecation shim all stay in sync.
.scale_choices <- c("to_01", "to_255", "binary", "none")

#' Resolve the single `scale` argument, mapping the deprecated logical flags
#'
#' `normalize`, `binarize`, and `denormalize` are kept only for backward
#' compatibility. Exactly one source of truth is allowed: either `scale` or the
#' old flags, never both.
#' @keywords internal
resolve_scale <- function(scale, normalize, binarize, denormalize, scale_given) {
  flags_given <- !is.null(normalize) || !is.null(binarize) || !is.null(denormalize)
  if (flags_given) {
    if (scale_given) {
      stop("Specify either `scale` or the deprecated normalize/binarize/denormalize flags, not both.")
    }
    warning("`normalize`, `binarize`, and `denormalize` are deprecated; use ",
            "`scale` (\"to_01\", \"to_255\", \"binary\", \"none\") instead.",
            call. = FALSE)
    normalize   <- isTRUE(normalize)
    binarize    <- isTRUE(binarize)
    denormalize <- isTRUE(denormalize)
    if (normalize && denormalize) {
      stop("normalize and denormalize are opposite transforms; set at most one to TRUE.")
    }
    return(if (binarize) "binary"
           else if (denormalize) "to_255"
           else if (normalize) "to_01"
           else "none")
  }
  match.arg(scale, .scale_choices)
}

#' Validate conversion parameters
#' @keywords internal
validate_conversion_params <- function(input, scale, select.layer) {
  if (missing(input)) stop("Input is required")
  if (!is.character(scale) || length(scale) != 1 || !scale %in% .scale_choices) {
    stop("scale must be one of: ", paste(.scale_choices, collapse = ", "))
  }

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

#' Rescale the array according to `scale`
#'
#' All conversions use fixed factors (255), never a per-image max, and are
#' guarded so a conversion is a no-op when the data is already in the target
#' range:
#' \itemize{
#'   \item `"none"`   leave values untouched
#'   \item `"to_01"`  0-255 -> 0-1 (divide by 255; skipped if already <= 1)
#'   \item `"to_255"` 0-1 -> 0-255 (multiply by 255; skipped if already > 1)
#'   \item `"binary"` strictly 0/1 via ceiling(arr / max)
#' }
#' @keywords internal
normalize_array <- function(arr, scale = "none") {
  if (identical(scale, "none")) return(arr)
  max_val <- max(as.vector(arr), na.rm = TRUE)
  if (is.na(max_val)) {
    warning("Cannot rescale: all values are NA")
    return(arr)
  }
  if (max_val > 255) warning("Maximum value exceeds 255, rescaling might be incorrect")
  switch(scale,
    binary = ceiling(arr / max_val),
    to_01  = if (max_val > 1)  arr / 255 else arr,
    to_255 = if (max_val <= 1) arr * 255 else arr,
    arr)
}

#' Load an image flexibly from file or convert from memory
#' @param input File path or image object
#' @param output_format Character, one out of "cimg", "spatrast", "matrix", "array", "brick", "raster", "spatrast", "magick-image". Other spellings are accepted.
#' @param scale Character, the value rescaling to apply. One of
#'   `"to_01"` (0-255 -> 0-1, the default), `"to_255"` (0-1 -> 0-255),
#'   `"binary"` (strictly 0/1), or `"none"` (leave values untouched). Each
#'   conversion is a no-op if the data is already in the target range.
#' @param select.layer Numeric, which layer to select if input has multiple layers
#' @param normalize,binarize,denormalize Deprecated logical flags kept for
#'   backward compatibility; use `scale` instead. `normalize = TRUE` maps to
#'   `scale = "to_01"`, `denormalize = TRUE` to `"to_255"`, and
#'   `binarize = TRUE` to `"binary"`. Supplying these together with `scale`, or
#'   setting both `normalize` and `denormalize`, is an error.
#' @export
load_flexible_image <- function(input, output_format = "cimg",
                                scale = c("to_01", "to_255", "binary", "none"),
                                select.layer = NULL,
                                normalize = NULL, binarize = NULL, denormalize = NULL) {
  scale <- resolve_scale(scale, normalize, binarize, denormalize,
                         scale_given = !missing(scale))
  validate_conversion_params(input, scale, select.layer)


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

  # Apply the requested value rescaling
  arr <- normalize_array(arr, scale)



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

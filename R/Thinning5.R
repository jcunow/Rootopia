##################################
# Skeletonization pipeline (clean)
# Zhang-Suen LUT thinning + endpoints + overlay
##################################



#' Zhang-Suen thinning using lookup table (LUT implementation)
#'
#' Performs iterative skeletonization of a binary raster using a
#' lookup-table encoding of 3x3 neighbourhood configurations.
#'
#' The algorithm operates as follows:
#'
#' 1. Raster values are extracted once into a numeric vector.
#' 2. Each iteration reconstructs a matrix representation of the image.
#' 3. A zero-padded border is added around the image.
#' 4. For each pixel equal to 1:
#'    - The 3x3 neighbourhood is extracted
#'    - A weighted sum (mask encoding) produces a unique neighbourhood code
#' 5. The code is used as an index into a 256-entry lookup table:
#'    - LUT value 1 or 3 -> pixel removed in first sub-step
#'    - LUT value 2 or 3 -> pixel removed in second sub-step
#' 6. Pixels are updated in-place in the vector representation
#' 7. Iteration stops when no pixels are removed or max_iter is reached
#'
#' Important properties of implementation:
#' - Raster is only written back once at the end
#' - All intermediate steps operate on in-memory matrix/vector data
#' - Neighborhood encoding is performed via explicit loops over image pixels
#'
#' @param img Binary SpatRaster (single layer only)
#' @param max_iter Maximum number of thinning iterations
#' @param verbose Logical. If TRUE, prints iteration count and pixel reduction
#'
#' @return Binary SpatRaster representing skeletonized image
#'
#' @keywords internal
lut_thin_fast <- function(img, max_iter = 200L, verbose = FALSE) {
  
  if (terra::nlyr(img) != 1L) stop("single-layer SpatRaster required")
  
  dims <- dim(img)
  nr <- dims[1]; nc <- dims[2]
  
  v <- as.integer(terra::values(img))
  n_start <- sum(v == 1L, na.rm = TRUE)
  
  lut <- as.integer(c(
    0,0,0,1,0,0,1,3,0,0,3,1,1,0,1,3,0,0,0,0,0,0,0,0,2,0,2,0,3,0,3,3,
    0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,0,0,0,3,0,2,2,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    2,0,0,0,0,0,0,0,2,0,0,0,2,0,0,0,3,0,0,0,0,0,0,0,3,0,0,0,3,0,2,0,
    0,1,3,1,0,0,1,3,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
    3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    2,3,1,3,0,0,1,3,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    2,3,0,1,0,0,0,1,0,0,0,0,0,0,0,0,3,3,0,1,0,0,0,0,2,2,0,0,2,0,0,0
  ))
  
  mask <- matrix(c(
    1,2,4,
    128,0,8,
    64,32,16
  ), nrow = 3, byrow = TRUE)
  
  iter <- 0L
  
  repeat {
    iter <- iter + 1L
    if (iter > max_iter) break
    
    changed <- FALSE
    
    m <- matrix(v, nrow = nr, ncol = nc, byrow = TRUE)
    
    mp <- matrix(0L, nr + 2, nc + 2)
    mp[2:(nr+1), 2:(nc+1)] <- m

    # vectorized neighborhood code (replaces the per-pixel double loop)
    nw<-mp[1:nr,1:nc];      nn<-mp[1:nr,2:(nc+1)];      ne<-mp[1:nr,3:(nc+2)]
    ww<-mp[2:(nr+1),1:nc];                              ee<-mp[2:(nr+1),3:(nc+2)]
    sw<-mp[3:(nr+2),1:nc];  ss<-mp[3:(nr+2),2:(nc+1)];  se<-mp[3:(nr+2),3:(nc+2)]
    Cmat <- 1L*nw + 2L*nn + 4L*ne + 128L*ww + 8L*ee + 64L*sw + 32L*ss + 16L*se
    Cmat[m != 1L] <- 0L
    code <- as.vector(t(Cmat))
    lut_code <- lut[code + 1L]
    
    rm1 <- lut_code %in% c(1L, 3L)
    rm2 <- lut_code %in% c(2L, 3L)
    
    if (any(rm1)) { v[rm1] <- 0L; changed <- TRUE }
    if (any(rm2)) { v[rm2] <- 0L; changed <- TRUE }
    
    if (!changed) break
  }
  
  if (verbose) {
    message(sprintf("LUT thinning: %d iters | %d -> %d px",
                    iter, n_start, sum(v == 1L)))
  }
  
  out <- img
  terra::values(out) <- v
  out
}





#' Detect endpoints and branching points in a skeleton image
#'
#' Computes local connectivity of each foreground pixel using an
#' 8-neighbourhood (Moore neighbourhood).
#'
#' The computation proceeds as follows:
#'
#' 1. Raster values are converted into a numeric matrix.
#' 2. A 1-pixel zero-padding border is added around the matrix.
#' 3. For each pixel in the original image:
#'    - A 3x3 window is extracted from the padded matrix
#'    - The center pixel is excluded
#'    - Neighbor count is computed as the sum of remaining 8 values
#'
#' Classification rules:
#' - Endpoint: pixel == 1 AND neighbor count == 1
#' - Branch point: pixel == 1 AND neighbor count > 2
#'
#' Outputs are converted back into SpatRaster objects.
#'
#' @param img Binary skeleton image. If \code{skeletonize = TRUE}, a segmented
#'   (non-skeleton) mask can be supplied instead.
#' @param select_layer Layer index for multi-layer rasters
#' @param skeletonize Logical. If \code{TRUE}, \code{img} is treated as a
#'   segmented mask and reduced to a skeleton internally via
#'   \code{skeletonize_image()} before detecting points. Default
#'   \code{FALSE} (assumes \code{img} is already a skeleton).
#'
#' @return List with:
#' \item{endpoints}{SpatRaster marking pixels with exactly one neighbor}
#' \item{branching_points}{SpatRaster marking pixels with more than two neighbors}
#'
#' @export
#' @examples
#' # Load example binary segmentation
#' data(seg_Oulanka2023_Session01_T067)
#'
#' # Ensure single-layer raster if needed
#' img <- seg_Oulanka2023_Session01_T067
#'
#' # Skeletonize
#' skel <- skeletonize_image(img, select_layer = 2, verbose = FALSE)
#' skel.points <- detect_skeleton_points(skel)
#' 
detect_skeleton_points <- function(img, select_layer = NULL, skeletonize = FALSE) {

  img <- load_flexible_image(img, output_format = "spatrast", scale = "binary",
                             select_layer = select_layer)
  if (terra::nlyr(img) > 1) img <- img[[1]]
  if (skeletonize) img <- skeletonize_image(img, verbose = FALSE)

  img <- matrix(terra::values(img), nrow = nrow(img), byrow = TRUE)

  nr <- nrow(img); nc <- ncol(img)
  padded <- matrix(0, nr + 2, nc + 2)
  padded[2:(nr+1), 2:(nc+1)] <- img
  # vectorized 8-neighbor count (replaces the per-pixel double loop)
  n <- padded[1:nr,1:nc]     + padded[1:nr,2:(nc+1)]     + padded[1:nr,3:(nc+2)] +
       padded[2:(nr+1),1:nc] +                            padded[2:(nr+1),3:(nc+2)] +
       padded[3:(nr+2),1:nc] + padded[3:(nr+2),2:(nc+1)] + padded[3:(nr+2),3:(nc+2)]
  
  list(
    endpoints = terra::rast((img == 1 & n == 1) * 1),
    branching_points = terra::rast((img == 1 & n > 2) * 1)
  )
}





#' Skeletonize binary image
#'
#' Applies iterative morphological thinning using a LUT-based
#' implementation of the Zhang-Suen skeletonization algorithm.
#'
#' Processing steps:
#'
#' 1. Input is converted to a single-layer binary SpatRaster using
#'    `load_flexible_image()`.
#' 2. Foreground pixel count is computed.
#' 3. Skeletonization is performed using `lut_thin_fast()`:
#'    - iterative removal of pixels based on 3x3 neighbourhood codes
#'    - lookup table determines pixel deletions in two sub-steps per iteration
#' 4. Output is the final thinned binary raster.
#' 5. Optionally, an overlay image is generated:
#'    - original image marked as base layer
#'    - skeleton pixels overlaid in a separate class
#'    - saved using base R PNG plotting
#'
#' @param img Input image (SpatRaster, matrix, or compatible format)
#' @param verbose Logical. Print summary statistics
#' @param select_layer Layer index for multi-layer rasters
#' @param overlay_png_path Optional file path for saving overlay visualization
#'
#' @return Binary SpatRaster representing skeletonized image
#' 
#' @examples
#' library(terra)
#'
#' # Load example binary segmentation
#' data(seg_Oulanka2023_Session01_T067)
#'
#' # Ensure single-layer raster if needed
#' img <- seg_Oulanka2023_Session01_T067
#'
#' # Skeletonize
#' skel <- skeletonize_image(img, select_layer = 2, verbose = FALSE)
#'
#' # Visual check
#' \dontrun{
#' skeletonize_image(img, select_layer = 2, overlay_png_path = "overlay.png")
#' }
#'
#' @export
skeletonize_image <- function(img,
                              verbose = TRUE,
                              select_layer = NULL,
                              overlay_png_path = NULL) {
  
  if (verbose) cat("Skeletonizing...\n")
  
  img_rast <- load_flexible_image(
    img,
    output_format = "spatrast",
    scale = "binary",
    select_layer = select_layer
  )[[1]]
  
  n_root <- sum(terra::values(img_rast) == 1, na.rm = TRUE)
  if (n_root == 0) stop("No foreground pixels")
  
  result <- lut_thin_fast(img_rast, verbose = verbose)
  
  n_skel <- sum(terra::values(result) == 1, na.rm = TRUE)
  
  if (verbose) {
    cat(sprintf("Root: %d -> Skeleton: %d\n", n_root, n_skel))
  }
  
  # -------------------------
  # overlay
  # -------------------------
  if (!is.null(overlay_png_path)) {
    
    overlay_png_path <- path.expand(overlay_png_path)
    
    # ----------------------------
    # extract matrices once
    # ----------------------------
    clean <- matrix(
      terra::values(img_rast),
      nrow = nrow(img_rast),
      byrow = TRUE
    )
    
    skel <- matrix(
      terra::values(result),
      nrow = nrow(result),
      byrow = TRUE
    )
    
    # ----------------------------
    # build categorical map
    # ----------------------------
    overlay <- matrix(0L, nrow(clean), ncol(clean))
    overlay[clean == 1] <- 1
    overlay[skel == 1]  <- 2
    
    # reverse rows for plotting orientation
    overlay <- overlay[nrow(overlay):1, , drop = FALSE]
    
    # ----------------------------
    # color palette (scientific / neutral)
    # ----------------------------
    pal <- c(
      "white",     # 0 background
      "gray65",    # 1 root structure
      "#0072B2" # 2 skeleton (Okabe-Ito blue)
    )
    
    overlay_col <- matrix(
      pal[overlay + 1],
      nrow = nrow(overlay),
      byrow = FALSE
    )
    
    # ----------------------------
    # write pixel-exact PNG
    # ----------------------------
    grDevices::png(
      filename = overlay_png_path,
      width = ncol(overlay_col),
      height = nrow(overlay_col),
      bg = "white",
      type = "cairo-png"
    )
    
    graphics::par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
    graphics::plot.new()
    
    graphics::rasterImage(
      grDevices::as.raster(overlay_col),
      0, 0, 1, 1
    )
    
    dev.off()
  }
  
  result
}

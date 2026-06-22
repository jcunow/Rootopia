# =============================================================================
# soil_rgb_classification.R
#
# Functions for classifying minirhizotron RGB images into soil material classes
# based on fixed LAB centroids derived from manual color calibration.
#
# Public functions:
#   classify_soil_rgb()       -- classify a SpatRaster, return map + metrics
#   build_soil_centroids()    -- derive centroids from user RGB picks
#   plot_soil_classification() -- visualize output of classify_soil_rgb()
#
# Internal helpers (not exported):
#   .rgb_to_lab()
#   .lab_to_rgb_hex()
#   .assign_nearest_centroid()
#   .default_soil_centroids()
# =============================================================================


# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' @keywords internal
.rgb_to_lab <- function(mat) {
  # mat: N x 3 matrix, RGB values 0-255
  # returns N x 3 matrix in CIE LAB (D65 illuminant)
  mat <- matrix(mat, ncol = 3) / 255   # guard against dimension drop
  lin <- ifelse(mat <= 0.04045, mat / 12.92, ((mat + 0.055) / 1.055)^2.4)
  lin <- matrix(lin, ncol = 3)
  M   <- matrix(c(0.4124564, 0.3575761, 0.1804375,
                  0.2126729, 0.7151522, 0.0721750,
                  0.0193339, 0.1191920, 0.9503041),
                nrow = 3, byrow = TRUE)
  xyz      <- matrix(lin %*% t(M), ncol = 3)  # explicit matrix coercion
  xyz[, 1] <- xyz[, 1] / 0.95047
  xyz[, 3] <- xyz[, 3] / 1.08883
  f <- function(t) ifelse(t > 0.008856, t^(1/3), 7.787 * t + 16/116)
  fx <- f(xyz[, 1]); fy <- f(xyz[, 2]); fz <- f(xyz[, 3])
  matrix(cbind(L = 116 * fy - 16, A = 500 * (fx - fy), B = 200 * (fy - fz)), ncol = 3)
}


#' @keywords internal
.lab_to_rgb_hex <- function(L, A, B) {
  # single LAB triplet -> "#RRGGBB" hex string
  fy    <- (L + 16) / 116
  fx    <- A / 500 + fy
  fz    <- fy - B / 200
  f_inv <- function(t) ifelse(t^3 > 0.008856, t^3, (t - 16/116) / 7.787)
  xyz   <- c(f_inv(fx) * 0.95047, f_inv(fy), f_inv(fz) * 1.08883)
  M_inv <- matrix(c( 3.2404542, -1.5371385, -0.4985314,
                     -0.9692660,  1.8760108,  0.0415560,
                     0.0556434, -0.2040259,  1.0572252),
                  nrow = 3, byrow = TRUE)
  lin  <- pmax(0, pmin(1, as.vector(M_inv %*% xyz)))
  srgb <- ifelse(lin <= 0.0031308, 12.92 * lin, 1.055 * lin^(1/2.4) - 0.055)
  srgb <- round(pmax(0, pmin(255, srgb * 255)))
  sprintf("#%02X%02X%02X", srgb[1], srgb[2], srgb[3])
}


#' @keywords internal
.assign_nearest_centroid <- function(lab_pixels, centroids_lab, max_dist_vec) {
  # Vectorised nearest-centroid assignment with per-class distance thresholds.
  # Uses: ||px - c||^2 = ||px||^2 - 2*(px.c) + ||c||^2
  # Returns integer vector: 1..K or 0 (unclassified)
  N <- nrow(lab_pixels)
  K <- nrow(centroids_lab)
  C <- as.matrix(centroids_lab)
  
  # N x K squared distance matrix -- keep explicit matrix() at every step
  # to prevent pmax/arithmetic from silently dropping dimensions
  px_sq  <- rowSums(lab_pixels^2)                              # length N
  cross  <- matrix(tcrossprod(lab_pixels, C), nrow = N)        # N x K
  c_sq   <- matrix(rep(rowSums(C^2), each = N), nrow = N)      # N x K
  sq_d   <- matrix(pmax(0, px_sq - 2 * cross + c_sq), nrow = N) # N x K
  distmat <- matrix(sqrt(sq_d), nrow = N)                      # N x K
  
  # Mask classes beyond their per-class threshold
  thresh  <- matrix(rep(max_dist_vec, each = N), nrow = N)     # N x K
  masked  <- distmat
  masked[distmat > thresh] <- Inf
  
  # Nearest unmasked class; 0 if all masked
  nearest <- max.col(-masked, ties.method = "first")
  nearest[matrix(rowSums(matrix(is.finite(masked), nrow = N)), ncol = 1) == 0L] <- 0L
  nearest
}


#' @keywords internal
.default_soil_centroids <- function() {
  # Default centroids calibrated from manual RGB picks on Oulanka 2023
  # minirhizotron images (Blended scans, Session 03).
  #
  # LAB space (D65). MAX_DIST is per-class Euclidean LAB threshold;
  # pixels beyond this from all centroids are labeled "unclassified".
  #
  # Close pairs to be aware of:
  #   dark_soil <-> red_soil  dist ~10 LAB units  (monitor unclassified %)
  #
  # These values are specific to this scanner / site / protocol.
  # Use build_soil_centroids() to derive centroids for your own data.
  data.frame(
    class    = c("dark_soil", "red_soil", "root",  "silver_tape", "coarse_debris"),
    L        = c( 8.2,         14.3,       35.2,    66.7,           23.1         ),
    A        = c(-2.3,          5.4,        2.1,    -3.0,           11.0         ),
    B        = c( 0.8,          5.9,        7.5,    -2.4,           13.1         ),
    MAX_DIST = c(14,            14,         26,      28,             11           ),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# build_soil_centroids()
# =============================================================================

#' Build soil class centroids from manual RGB color picks
#'
#' Converts raw RGB pixel picks (collected e.g. in QGIS, FIJI, or ImageJ) to
#' CIE LAB centroids and returns a centroid table ready to pass to
#' \code{\link{classify_soil_rgb}}. Also prints diagnostic summaries including
#' intra-class spread, per-class coverage, and inter-class distance warnings.
#'
#' @param picks A named list of RGB pick matrices. Each element corresponds to
#'   one class and must be a numeric matrix with 3 columns (R, G, B), values
#'   0-255, with one row per pick. Names become class names in the output.
#'   See \code{\link{classify_soil_rgb}}'s \strong{Building your own
#'   centroids (picks)} section for a worked example of constructing these
#'   matrices from an image (e.g. by cropping representative patches and
#'   calling \code{terra::values()}).
#' @param max_dist A named numeric vector of per-class LAB distance thresholds,
#'   matched by name to \code{picks}. Pixels further than this from a class
#'   centroid cannot be assigned to that class.
#' @param prior Optional \code{data.frame} of existing centroids (same format
#'   as the output of this function, or \code{.default_soil_centroids()}).
#'   When supplied, the new centroids derived from \code{picks} are blended
#'   with the prior centroids using \code{alpha}. Only classes present in both
#'   \code{picks} and \code{prior} are blended; new classes in \code{picks}
#'   that are absent from \code{prior} are added as-is.
#' @param alpha Numeric in [0, 1]. Blend weight for the prior centroids.
#'   \code{alpha = 0} (default) ignores the prior entirely -- centroids are
#'   derived purely from the new picks. \code{alpha = 1} returns the prior
#'   unchanged. \code{alpha = 0.5} weights old and new equally. Decrease
#'   \code{alpha} progressively as you collect more new picks to gradually
#'   shift calibration toward the new dataset.
#' @param verbose Logical. Print per-class summaries and inter-class distance
#'   matrix. Default \code{TRUE}.
#'
#' @return A \code{data.frame} with columns \code{class}, \code{L}, \code{A},
#'   \code{B}, \code{MAX_DIST}.
#'
#' @examples
#' \dontrun{
#' # Clean break -- new picks only
#' cents <- build_soil_centroids(new_picks, max_dist)
#'
#' # Blend: 30% old calibration, 70% new picks
#' cents <- build_soil_centroids(new_picks, max_dist,
#'                               prior = .default_soil_centroids(),
#'                               alpha = 0.3)
#'
#' # Iterative refinement across sessions:
#' # session 1
#' cents <- build_soil_centroids(picks_s1, max_dist)
#' # session 2 -- downweight session 1 to 20%
#' cents <- build_soil_centroids(picks_s2, max_dist, prior = cents, alpha = 0.2)
#' # session 3 -- downweight accumulated prior to 10%
#' cents <- build_soil_centroids(picks_s3, max_dist, prior = cents, alpha = 0.1)
#' }
#'
#' @export
build_soil_centroids <- function(picks, max_dist, prior = NULL, alpha = 0,
                                 verbose = TRUE) {
  
  if (!is.list(picks) || is.null(names(picks)))
    stop("picks must be a named list of RGB matrices")
  if (!is.numeric(max_dist) || is.null(names(max_dist)))
    stop("max_dist must be a named numeric vector")
  if (!all(names(picks) %in% names(max_dist)))
    stop("max_dist must have an entry for every class in picks: ",
         paste(setdiff(names(picks), names(max_dist)), collapse = ", "))
  if (!is.null(prior)) {
    if (!is.data.frame(prior) || !all(c("class","L","A","B") %in% names(prior)))
      stop("prior must be a data.frame with columns class, L, A, B")
    if (!is.numeric(alpha) || alpha < 0 || alpha > 1)
      stop("alpha must be a number between 0 and 1")
  }
  
  # Compute new centroids from picks
  results <- lapply(names(picks), function(cls) {
    mat <- picks[[cls]]
    if (!is.matrix(mat) || ncol(mat) != 3)
      stop("picks[['", cls, "']] must be a matrix with 3 columns (R, G, B)")
    
    lab  <- .rgb_to_lab(mat)
    lm   <- colMeans(lab)
    lsd  <- if (nrow(lab) > 1) apply(lab, 2, stats::sd) else c(0, 0, 0)
    to_c <- sqrt(rowSums(sweep(lab, 2, lm)^2))
    
    if (nrow(lab) > 1) {
      pairs      <- utils::combn(nrow(lab), 2)
      pair_dists <- apply(pairs, 2, function(ij)
        sqrt(sum((lab[ij[1], ] - lab[ij[2], ])^2)))
      spread_mean <- mean(pair_dists)
      spread_max  <- max(pair_dists)
    } else {
      spread_mean <- spread_max <- 0
    }
    
    md <- max_dist[cls]
    if (verbose) {
      cat(sprintf("\n%-16s  n=%d\n", cls, nrow(mat)))
      cat(sprintf("  LAB centroid : L=%.1f  A=%.1f  B=%.1f  hex=%s\n",
                  lm[1], lm[2], lm[3], .lab_to_rgb_hex(lm[1], lm[2], lm[3])))
      cat(sprintf("  LAB sd       : L=%.1f  A=%.1f  B=%.1f\n",
                  lsd[1], lsd[2], lsd[3]))
      cat(sprintf("  Intra spread : mean=%.1f  max=%.1f\n",
                  spread_mean, spread_max))
      cat(sprintf("  Dist to ctr  : mean=%.1f  max=%.1f  p90=%.1f\n",
                  mean(to_c), max(to_c), stats::quantile(to_c, 0.9)))
      cat(sprintf("  MAX_DIST     : %g  (covers %.0f%% of picks)\n",
                  md, 100 * mean(to_c <= md)))
    }
    list(centroid = lm)
  })
  
  lms <- do.call(rbind, lapply(results, `[[`, "centroid"))
  
  # Blend with prior if supplied
  if (!is.null(prior) && alpha > 0) {
    for (i in seq_along(names(picks))) {
      cls <- names(picks)[i]
      prior_row <- prior[prior$class == cls, , drop = FALSE]
      if (nrow(prior_row) == 1) {
        prior_lab <- as.numeric(prior_row[, c("L", "A", "B")])
        lms[i, ]  <- (1 - alpha) * lms[i, ] + alpha * prior_lab
        if (verbose)
          cat(sprintf("\n  [%s] blended with prior (alpha=%.2f)  ->  L=%.1f  A=%.1f  B=%.1f\n",
                      cls, alpha, lms[i, 1], lms[i, 2], lms[i, 3]))
      }
    }
  }
  
  # Inter-class distances + warnings
  if (verbose) {
    cat("\n\n=== Inter-class LAB distances ===\n")
    d <- round(as.matrix(stats::dist(lms)), 1)
    rownames(d) <- colnames(d) <- names(picks)
    print(d)
    close <- which(d < 15 & d > 0, arr.ind = TRUE)
    close <- close[close[, 1] < close[, 2], , drop = FALSE]
    if (nrow(close) > 0) {
      cat("\nWARNING -- close centroid pairs (< 15 LAB units):\n")
      for (i in seq_len(nrow(close))) {
        r <- close[i, 1]; cc <- close[i, 2]
        cat(sprintf("  %-16s <-> %-16s  %.1f\n",
                    names(picks)[r], names(picks)[cc], d[r, cc]))
      }
    }
  }
  
  out <- data.frame(
    class    = names(picks),
    L        = round(lms[, 1], 1),
    A        = round(lms[, 2], 1),
    B        = round(lms[, 3], 1),
    MAX_DIST = unname(max_dist[names(picks)]),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}


# =============================================================================
# classify_soil_rgb()
# =============================================================================

#' Classify soil material classes from a minirhizotron RGB raster
#'
#' Assigns each pixel of an RGB \code{SpatRaster} to a class
#' (e.g. dark soil, red soil, root, silver tape, coarse debris) by nearest-
#' centroid assignment in CIE LAB color space. Pixels beyond the per-class
#' distance threshold are labeled "unclassified". 
#'
#' @param img An RGB image with at least 3 layers/channels interpreted as
#'   R, G, B (in that order). Values may be 0-255 or 0-1 (auto-detected).
#'   Converted internally via \code{load_flexible_image()}, so any of its
#'   supported formats are accepted: file path (.jpg, .jpeg, .png, .tif,
#'   .tiff, .bmp), \code{SpatRaster}, \code{RasterLayer}/\code{RasterBrick},
#'   \code{cimg}, \code{magick-image}, \code{matrix}, or \code{array}.
#' @param centroids A \code{data.frame} with columns \code{class}, \code{L},
#'   \code{A}, \code{B}, \code{MAX_DIST} (see \strong{Centroid table format}
#'   below). Defaults to a set of centroids calibrated on Oulanka 2023
#'   minirhizotron scans -- see \strong{Building your own centroids} below
#'   and \code{\link{build_soil_centroids}} to derive centroids for your own
#'   data.
#' @param downsample_fact Integer spatial aggregation factor applied before
#'   classification for speed. \code{NULL} (default) uses full resolution.
#'   The output map is always disaggregated back to match the input resolution
#'   and extent exactly (nearest-neighbor, no interpolation).
#' @param compute_metrics Logical. Compute per-class pixel counts, area
#'   fractions, LAB statistics, and mean distance to centroid. Default
#'   \code{TRUE}. Set \code{FALSE} in tight batch loops where only the map is
#'   needed.
#' @param verbose Logical. Print progress messages and summary table. Default
#'   \code{TRUE}.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{map}}{A \code{SpatRaster} of integer class IDs with factor
#'     levels set to class names. Level 0 = "unclassified". Same CRS, extent,
#'     and resolution as \code{img}. Use directly with \code{terra::zonal()},
#'     \code{terra::mask()}, \code{terra::freq()}.}
#'   \item{\code{metrics}}{A \code{data.frame} with per-class pixel counts,
#'     area fractions (\%), LAB and RGB means and SDs, mean distance to
#'     centroid, and the actual mean color rendered as hex. \code{NULL} if
#'     \code{compute_metrics = FALSE}.}
#'   \item{\code{inter_dist}}{Numeric matrix of pairwise LAB distances between
#'     class centroids.}
#'   \item{\code{centroids}}{The centroid table used (useful when the default
#'     was applied).}
#' }
#'
#' @section Centroid table format:
#' The \code{centroids} table has one row per material class and these
#' columns:
#' \describe{
#'   \item{\code{class}}{Character. Name of the class (e.g.
#'     \code{"dark_soil"}, \code{"root"}). Becomes the factor level in
#'     \code{result$map} and the \code{class} column of \code{result$metrics}.}
#'   \item{\code{L}, \code{A}, \code{B}}{Numeric. The class centroid's
#'     coordinates in CIE LAB color space (D65 illuminant): \code{L} is
#'     lightness (0 = black, 100 = white), \code{A} is the green-red axis
#'     (negative = green, positive = red), and \code{B} is the blue-yellow
#'     axis (negative = blue, positive = yellow). These are the *mean* LAB
#'     values of the representative pixels for that class.}
#'   \item{\code{MAX_DIST}}{Numeric. The per-class assignment radius, in
#'     Euclidean LAB units. A pixel is assigned to the class whose centroid
#'     is nearest in LAB space, but only if that distance is \code{<=
#'     MAX_DIST}; otherwise the pixel is labeled \code{"unclassified"}.
#'     Larger values classify more pixels but risk merging visually distinct
#'     materials; smaller values leave more pixels unclassified. Typical
#'     values are roughly 10-30.}
#' }
#' Pixel RGB values are converted to LAB internally before distances are
#' computed, so \code{L}/\code{A}/\code{B} are never raw RGB numbers.
#'
#' @section Building your own centroids (picks):
#' The default \code{centroids} table was calibrated on one specific scanner
#' and site, so for other data you should derive your own via
#' \code{\link{build_soil_centroids}}.
#'
#' \code{build_soil_centroids()} takes \code{picks}: a named list with one
#' element per material class. Each element is a numeric matrix with exactly
#' 3 columns (R, G, B in 0-255), where every row is one color sample
#' believed to belong to that class. Classes may have different numbers of
#' rows. The function converts each matrix to LAB, averages it to a single
#' centroid, and returns a \code{data.frame} in the same format as
#' \code{.default_soil_centroids()} -- ready to pass straight back into
#' \code{classify_soil_rgb(centroids = ...)}.
#'
#' The simplest approach is to read representative RGB values off your scan
#' (e.g. using an image viewer's color picker) and enter them directly:
#' \preformatted{
#' picks <- list(
#'   dark_soil = matrix(c( 28,  22,  18,
#'                          32,  26,  21,
#'                          25,  19,  15), ncol = 3, byrow = TRUE),
#'   red_soil  = matrix(c( 80,  45,  35,
#'                          75,  42,  31), ncol = 3, byrow = TRUE),
#'   root      = matrix(c(180, 160, 130,
#'                        175, 155, 125,
#'                        185, 165, 135,
#'                        178, 158, 128), ncol = 3, byrow = TRUE),
#'   tape      = matrix(c(200, 205, 210,
#'                        195, 200, 205), ncol = 3, byrow = TRUE),
#'   debris    = matrix(c(100,  70,  45,
#'                         95,  65,  40), ncol = 3, byrow = TRUE)
#' )
#'
#' max_dist <- c(dark_soil = 14, red_soil = 14, root = 26,
#'               tape = 28, debris = 11)
#'
#' cents  <- build_soil_centroids(picks, max_dist)   # prints diagnostics
#' result <- classify_soil_rgb(img, centroids = cents)
#' }
#' Alternatively, picks can be extracted from known representative patches of
#' your image. \strong{Important:} \code{build_soil_centroids()} always treats
#' pick values as 0-255 (no auto-scaling). If your raster stores values in
#' 0-1, multiply by 255 before building picks so that the centroids and the
#' pixels seen inside \code{classify_soil_rgb()} are on the same scale.
#'
#' \code{build_soil_centroids()} prints diagnostics (intra-class spread,
#' inter-class distances, \code{MAX_DIST} coverage) to help you sanity-check
#' your choices. Provide a class for \emph{every} material present in your
#' scans: because classification is nearest-centroid, any material without its
#' own class is snapped into whichever defined class is closest.
#'
#' @seealso \code{\link{build_soil_centroids}}, \code{\link{plot_soil_classification}}
#'
#' @examples
#' \dontrun{
#' library(terra)
#' img    <- rast("scan.tiff")
#' result <- classify_soil_rgb(img)
#'
#' # Access outputs
#' terra::plot(result$map, maxcell = Inf)
#' result$metrics
#'
#' # Zonal stats with a depth-band raster
#' zones  <- rast("depth_zones.tif")
#' zstats <- terra::zonal(result$map, zones, fun = "freq")
#'
#' # Custom centroids -- see "Building your own centroids (picks)" above for
#' # how to construct `picks` and `max_dist`
#' cents  <- build_soil_centroids(picks, max_dist)
#' result <- classify_soil_rgb(img, centroids = cents)
#' }
#'
#' @importFrom terra nlyr values aggregate disagg resample
#' @export
classify_soil_rgb <- function(img,
                              centroids       = .default_soil_centroids(),
                              downsample_fact = NULL,
                              compute_metrics = TRUE,
                              verbose         = TRUE) {
  
  msg <- function(...) if (verbose) message(...)

  # --- Prepare raster ---------------------------------------------------------
  # Accepts file path, SpatRaster, RasterLayer/Brick, cimg, magick-image,
  # matrix, or array -- see load_flexible_image() for the full list.
  img <- load_flexible_image(img, output_format = "spatrast",
                             scale = "none")
  if (terra::nlyr(img) < 3)
    stop("img must have at least 3 layers (R, G, B)")

  required_cols <- c("class", "L", "A", "B", "MAX_DIST")
  if (!all(required_cols %in% names(centroids)))
    stop("centroids must have columns: ", paste(required_cols, collapse = ", "))

  # --- Prepare raster --------------------------------------------------------
  img_work <- img[[1:3]]
  names(img_work) <- c("R", "G", "B")
  msg(sprintf("Input: %d cols x %d rows = %d pixels",
              terra::ncol(img_work), terra::nrow(img_work), terra::ncell(img_work)))
  
  # --- Optional spatial downsampling ----------------------------------------
  img_orig <- NULL
  if (!is.null(downsample_fact) && downsample_fact > 1) {
    img_orig <- img_work
    img_work <- terra::aggregate(img_work, fact = downsample_fact,
                                 fun = "mean", na.rm = TRUE)
    msg(sprintf("Aggregated x%d -> %d cols x %d rows = %d pixels",
                downsample_fact,
                terra::ncol(img_work), terra::nrow(img_work), terra::ncell(img_work)))
  }
  
  # --- Extract pixel matrix --------------------------------------------------
  pixels     <- terra::values(img_work)
  valid_mask <- stats::complete.cases(pixels)
  pixels_rgb <- pixels[valid_mask, , drop = FALSE]
  msg(sprintf("Valid pixels: %d (%.1f%%)",
              sum(valid_mask), 100 * mean(valid_mask)))
  
  # Auto-scale 0-1 -> 0-255
  if (max(pixels_rgb, na.rm = TRUE) <= 1)
    pixels_rgb <- pixels_rgb * 255
  pixels_rgb <- matrix(pixels_rgb, ncol = 3)  # guard against dimension drop
  
  # --- Convert to LAB --------------------------------------------------------
  msg("Converting RGB -> LAB...")
  pixels_lab    <- .rgb_to_lab(pixels_rgb)
  centroids_lab <- as.matrix(centroids[, c("L", "A", "B")])
  
  # --- Assign pixels to classes ----------------------------------------------
  msg("Assigning pixels to classes...")
  labels_valid <- .assign_nearest_centroid(pixels_lab, centroids_lab,
                                           centroids$MAX_DIST)
  n_unc <- sum(labels_valid == 0L)
  msg(sprintf("  Unclassified: %d (%.2f%%)",
              n_unc, 100 * n_unc / length(labels_valid)))
  
  # --- Map labels back to raster ---------------------------------------------
  labels_full <- rep(NA_integer_, terra::ncell(img_work))
  labels_full[valid_mask] <- labels_valid
  
  class_map        <- img_work[[1]]
  terra::values(class_map) <- labels_full
  names(class_map) <- "class_id"
  
  # Disaggregate to original resolution (nearest-neighbor, no interpolation)
  if (!is.null(img_orig)) {
    msg("Disaggregating to original resolution...")
    class_map <- terra::disagg(class_map, fact = downsample_fact, method = "near")
    class_map <- terra::resample(class_map, img_orig, method = "near")
    names(class_map) <- "class_id"
  }
  
  # Set factor levels (makes raster self-documenting for zonal/freq)
  lvls <- data.frame(ID    = c(0L, seq_len(nrow(centroids))),
                     class = c("unclassified", centroids$class))
  levels(class_map) <- list(lvls)
  
  # --- Inter-class distance matrix -------------------------------------------
  inter_dist <- as.matrix(stats::dist(centroids_lab))
  rownames(inter_dist) <- colnames(inter_dist) <- centroids$class
  
  # Warn about close pairs
  warn_dist <- 20
  close <- which(inter_dist < warn_dist & inter_dist > 0, arr.ind = TRUE)
  close <- close[close[, 1] < close[, 2], , drop = FALSE]
  if (nrow(close) > 0) {
    msg(sprintf("\nWARNING: centroid pairs closer than %g LAB units:", warn_dist))
    for (i in seq_len(nrow(close))) {
      r <- close[i, 1]; cc <- close[i, 2]
      msg(sprintf("  %s  <->  %s   dist = %.1f",
                  centroids$class[r], centroids$class[cc], inter_dist[r, cc]))
    }
  }
  
  # --- Per-class metrics -----------------------------------------------------
  metrics <- NULL
  if (compute_metrics) {
    msg("Computing metrics...")
    all_ids   <- c(seq_len(nrow(centroids)), 0L)
    all_names <- c(centroids$class, "unclassified")
    n_total   <- length(labels_valid)
    
    metrics <- do.call(rbind, lapply(seq_along(all_ids), function(i) {
      id    <- all_ids[i]
      cname <- all_names[i]
      idx   <- which(labels_valid == id)
      n_px  <- length(idx)
      
      if (n_px == 0L) {
        return(data.frame(
          class = cname, n_pixels = 0L, pct_pixels = 0,
          L_mean = NA_real_, A_mean = NA_real_, B_mean = NA_real_,
          L_sd = NA_real_, A_sd = NA_real_, B_sd = NA_real_,
          R_mean = NA_real_, G_mean = NA_real_, B_mean_rgb = NA_real_,
          mean_dist_to_centroid = NA_real_, hex_actual = NA_character_,
          stringsAsFactors = FALSE
        ))
      }
      
      px_lab <- pixels_lab[idx, , drop = FALSE]
      px_rgb <- pixels_rgb[idx, , drop = FALSE]
      lm     <- colMeans(px_lab)
      ls     <- apply(px_lab, 2, stats::sd)
      rm     <- colMeans(px_rgb)
      md     <- if (id > 0L)
        mean(sqrt(rowSums(sweep(px_lab, 2, centroids_lab[id, ])^2)))
      else
        NA_real_
      
      data.frame(
        class                 = cname,
        n_pixels              = n_px,
        pct_pixels            = round(100 * n_px / n_total, 3),
        L_mean                = round(lm[1], 1),
        A_mean                = round(lm[2], 1),
        B_mean                = round(lm[3], 1),
        L_sd                  = round(ls[1], 2),
        A_sd                  = round(ls[2], 2),
        B_sd                  = round(ls[3], 2),
        R_mean                = round(rm[1], 1),
        G_mean                = round(rm[2], 1),
        B_mean_rgb            = round(rm[3], 1),
        mean_dist_to_centroid = round(md, 2),
        hex_actual            = .lab_to_rgb_hex(lm[1], lm[2], lm[3]),
        stringsAsFactors      = FALSE
      )
    }))
    rownames(metrics) <- NULL
    
    if (verbose) {
      message("\n===== CLASS METRICS =====")
      print(metrics[, c("class", "n_pixels", "pct_pixels",
                        "mean_dist_to_centroid", "hex_actual")])
    }
  }
  
  list(map        = class_map,
       metrics    = metrics,
       inter_dist = inter_dist,
       centroids  = centroids)
}


# =============================================================================
# plot_soil_classification()
# =============================================================================

#' Plot the output of classify_soil_rgb
#'
#' Produces a three-panel figure: the classified raster map, a legend showing
#' class names with area fractions and mean distances to centroids, and a
#' dendrogram of inter-class LAB distances.
#'
#' @param result A list returned by \code{\link{classify_soil_rgb}}.
#' @param color_mode Character. One of \code{"contrast"} (default),
#'   \code{"vibrant"}, or \code{"centroid"}.
#'   \describe{
#'     \item{\code{"contrast"}}{Muted, distinguishable colors defined in
#'       \code{class_colors}.}
#'     \item{\code{"vibrant"}}{Saturated, high-contrast colors defined in
#'       \code{vibrant_colors}.}
#'     \item{\code{"centroid"}}{Each class rendered as its actual LAB centroid
#'       color -- useful for sanity-checking centroid values.}
#'   }
#' @param class_colors Named character vector of hex colors for
#'   \code{"contrast"} mode. Names must match class names in
#'   \code{result$centroids}. Defaults to a built-in palette.
#' @param vibrant_colors Named character vector of hex colors for
#'   \code{"vibrant"} mode. Same naming requirement. Defaults to a built-in
#'   palette.
#' @param save_png File path for PNG output, or \code{NULL} (default) to plot
#'   to the active graphics device.
#' @param width,height PNG dimensions in pixels. Only used if
#'   \code{save_png} is not \code{NULL}.
#'
#' @return Invisibly \code{NULL}. Called for its side effect (plot).
#'
#' @seealso \code{\link{classify_soil_rgb}}
#'
#' @importFrom terra values
#' @importFrom stats dist hclust as.dendrogram
#' @importFrom graphics par layout plot legend axis abline mtext
#' @importFrom grDevices png dev.off
#' @export
plot_soil_classification <- function(result,
                                     color_mode    = "contrast",
                                     class_colors  = NULL,
                                     vibrant_colors = NULL,
                                     save_png      = NULL,
                                     width         = 1800,
                                     height        = 900) {
  
  class_map  <- result$map
  metrics    <- result$metrics
  centroids  <- result$centroids
  inter_dist <- result$inter_dist
  
  # Default palettes
  default_contrast <- c(
    dark_soil     = "#1A1A1A",
    red_soil      = "#8B3A2A",
    root          = "#D4B483",
    silver_tape   = "#A8B8BA",
    coarse_debris = "#C47A35",
    unclassified  = "#FF69B4"
  )
  default_vibrant <- c(
    dark_soil     = "#3D00FF",
    red_soil      = "#FF2200",
    root          = "#FFD700",
    silver_tape   = "#00E5FF",
    coarse_debris = "#FF6600",
    unclassified  = "#FF00AA"
  )
  if (is.null(class_colors))    class_colors   <- default_contrast
  if (is.null(vibrant_colors))  vibrant_colors <- default_vibrant
  
  # Build id -> hex lookup
  pal <- switch(
    color_mode,
    centroid = {
      cols <- mapply(.lab_to_rgb_hex, centroids$L, centroids$A, centroids$B)
      c(stats::setNames(unname(cols), as.character(seq_len(nrow(centroids)))),
        "0" = unname(class_colors["unclassified"]))
    },
    vibrant = {
      cols <- sapply(centroids$class, function(n)
        if (n %in% names(vibrant_colors)) vibrant_colors[[n]] else "#CCCCCC")
      c(stats::setNames(unname(cols), as.character(seq_len(nrow(centroids)))),
        "0" = unname(vibrant_colors["unclassified"]))
    },
    {  # "contrast" (default)
      cols <- sapply(centroids$class, function(n)
        if (n %in% names(class_colors)) class_colors[[n]] else "#CCCCCC")
      c(stats::setNames(unname(cols), as.character(seq_len(nrow(centroids)))),
        "0" = unname(class_colors["unclassified"]))
    }
  )
  
  # Dendrogram for legend ordering + panel 3
  hc         <- stats::hclust(stats::dist(as.matrix(inter_dist)), method = "average")
  dend_order <- hc$order
  
  # Legend entries: named classes in dendrogram order, unclassified last
  leg_ids    <- c(as.character(dend_order), "0")
  leg_names  <- c(centroids$class[dend_order], "unclassified")
  leg_pct    <- if (!is.null(metrics))
    metrics$pct_pixels[match(leg_names, metrics$class)]
  else
    rep(NA_real_, length(leg_names))
  leg_dist   <- if (!is.null(metrics))
    metrics$mean_dist_to_centroid[match(leg_names, metrics$class)]
  else
    rep(NA_real_, length(leg_names))
  leg_labels <- sprintf("%-14s %5.1f%%  d=%.1f",
                        leg_names,
                        ifelse(is.na(leg_pct),  0, leg_pct),
                        ifelse(is.na(leg_dist), 0, leg_dist))
  leg_fills  <- unname(pal[leg_ids])
  
  # Map colors for IDs present in raster
  map_ids  <- sort(unique(stats::na.omit(terra::values(class_map))))
  map_cols <- unname(pal[as.character(map_ids)])
  
  # Open device
  if (!is.null(save_png))
    grDevices::png(save_png, width = width, height = height, res = 150)
  
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  
  graphics::layout(matrix(1:3, nrow = 1), widths = c(3, 1.5, 1.8))
  
  # Panel 1: class map
  graphics::par(mar = c(2, 2, 3, 0.5))
  # plot() on a numeric raster with a pre-built color vector keyed to sorted unique IDs
  plot_ids <- sort(unique(stats::na.omit(as.vector(terra::values(class_map)))))
  plot_cols <- unname(pal[as.character(plot_ids)])
  terra::plot(class_map,
              col    = plot_cols,
              type   = "classes",
              legend = FALSE,
              axes   = FALSE,
              main   = color_mode,
              mar    = c(2, 2, 3, 0.5))  # terra::plot mar arg prevents par reset
  
  # Panel 2: legend on white
  graphics::par(mar = c(2, 2, 3, 0.5))
  graphics::plot.new()
  graphics::legend("center",
                   legend = leg_labels,
                   fill   = leg_fills,
                   border = "gray70",
                   bty    = "n",
                   cex    = 0.80,
                   title  = expression(bold("Class          %area  mean d")),
                   title.adj = 0)
  
  # Panel 3: centroid distance dendrogram
  graphics::par(mar = c(4, 0.5, 3, 7),
                mgp = c(0.5, 0.2, 0))   # pull x-axis title upward)
  dend <- stats::as.dendrogram(hc)
  graphics::plot(dend,
                 horiz    = TRUE,
                 axes     = TRUE,
                 leaflab  = "none",
                 main     = "",
                 yaxt     = "n",
                 xlab     = "Euclidean distance (LAB units)",
                 cex.axis = 0.8,
                 cex.lab  = 0.8)
  
  graphics::title("Centroid distances (LAB)",
                  line = 0.5,   # smaller = closer to plot
                  cex.main = 1.1)
  graphics::axis(1, lwd = 2.5, lwd.ticks = 1.5, cex.axis = 0.85, col = "black", pos = 0.7)
  for (i in seq_along(labels(dend))) {
    cname <- labels(dend)[i]
    cid   <- as.character(which(centroids$class == cname))
    graphics::axis(4, at = i, labels = cname,
                   col.axis = pal[cid], las = 1,
                   tick = FALSE, cex.axis = 0.82, font = 2,
                   xpd = TRUE)
  }
  
  
  if (!is.null(save_png)) grDevices::dev.off()
  invisible(NULL)
}
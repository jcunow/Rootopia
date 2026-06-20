#' Compute root traits over a depth profile from segmented (mini)rhizotron images
#'
#' @description
#' Processes a set of segmented rhizotron or minirhizotron images and returns a
#' tidy data frame of root traits summarised per depth interval.  Supports both
#' cylindrical tube geometry (minirhizotrons) and flat window geometry
#' (rhizotron panels).
#'
#' Each metric group is toggled independently.  If a block fails for one image
#' it is replaced with \code{NA} columns and a message is printed; processing
#' always continues to the next image.  If an entire image cannot be loaded it
#' is dropped and listed in a warning at the end.
#'
#' @section Image ordering:
#' Files are matched by position after \code{list.files()} sorts them
#' alphabetically.  The segmented, skeleton, and RGB directories must therefore
#' contain files whose alphabetical order corresponds to the same physical
#' sample.  Use \code{skl_file_index} and \code{rgb_file_index} to subset those
#' directories if necessary.
#'
#' @section Which paths are required for which metrics:
#' \describe{
#'   \item{\code{path.seg}}{Always required.}
#'   \item{\code{path.skl}}{Optional. Used by \code{calc_root_length},
#'     \code{calc_diameter_stats}, \code{calc_diameter_quantiles},
#'     \code{calc_root_angles}, and \code{calc_root_order_metrics} when
#'     supplied. If \code{path.skl} is \code{NULL} or a skeleton file is
#'     missing for an image, the skeleton is computed internally from the
#'     segmented image via \code{skeletonize_image()}.}
#'   \item{\code{path.rgb}}{Required when \code{calc_color_metrics} is
#'     \code{TRUE}.}
#' }
#'
#' @param path.seg Character. Path to directory of binary segmented images
#'   (foreground/root pixel = 1, background = 0).
#' @param path.skl Character or \code{NULL}. Path to directory of skeletonised
#'   images (one-pixel-wide centrelines of roots), used for length, diameter,
#'   angle, and branching-order metrics. If \code{NULL} or a file is missing,
#'   the skeleton is computed internally via \code{skeletonize_image()}.
#'   Default \code{NULL}.
#' @param path.rgb Character or \code{NULL}. Path to directory of blended
#'   RGB images, aligned to the segmented images.  Required for colour metrics.
#'   Default \code{NULL}.
#' @param seg_file_index Integer vector or \code{NULL}. Optional subset index
#'   applied to \code{list.files(path.seg)}, e.g. \code{37:72}.  Default
#'   \code{NULL} (use all files).
#' @param skl_file_index Integer vector or \code{NULL}. Optional subset index
#'   applied to \code{list.files(path.skl)}.  Default \code{NULL}.
#' @param rgb_file_index Integer vector or \code{NULL}. Optional subset index
#'   applied to \code{list.files(path.rgb)}.  Default \code{NULL}.
#'
#' @section Per-image metadata:
#' All metadata arguments below accept either a single value (recycled to all
#' images) or a vector of length equal to the number of images.
#'
#' @param insertion_angles Numeric. Insertion angle of the tube or window from
#'   vertical, in \strong{degrees}.  \code{0} = perfectly vertical,
#'   \code{30} = tilted 30 degrees from vertical (a common minirhizotron angle).
#'   Used by \code{create_depthmap()} to correct the depth scale.
#'   Default \code{0}.
#' @param soil_starts Numeric. Pixel row (in the original, un-rotated image)
#'   at which the soil surface begins.  Used to set the zero-depth reference.
#'   Default \code{0}.
#' @param tube_names Character or \code{NULL}. Sample/tube identifiers added as
#'   the \code{Tube} column.  If \code{NULL}, names are derived from characters
#'   3-5 from the right of the segmented file name, prefixed with \code{"T"}
#'   (e.g. \code{"T042"}).  Adjust if your naming convention differs.
#'   Default \code{NULL}.
#' @param session Character. Session or campaign label added as the
#'   \code{Session} column, e.g. \code{"2022_02"}.  Default \code{""}.
#'
# --- Scan and geometry settings ---
#' @param dpi Numeric. Scanner resolution in \strong{dots per inch}.  Used to
#'   convert pixel distances to physical units (cm).  Default \code{300}.
#' @param tube_diameter_cm Numeric. Inner diameter of the minirhizotron tube in
#'   \strong{centimetres}.  Passed to \code{create_depthmap()} as
#'   \code{tube.thicc}.  Its role in flat-window geometry (\code{flat_geometry
#'   = TRUE}) is uncertain; it likely has a default inside \code{create_depthmap()}
#'   and may be ignored -- leave at the default unless you know it matters for
#'   your setup.  Default \code{7}.
#' @param depth_interval_cm Numeric. Size of each depth bin in
#'   \strong{centimetres}.  Passed as \code{nn} to \code{binning()}.
#'   Default \code{5}.
#' @param flat_geometry Logical.  If \code{FALSE} (default), images are treated
#'   as cylindrical minirhizotron tubes and a sinusoidal depth correction is
#'   applied (\code{sinoid = TRUE} in \code{create_depthmap()}).  Set to
#'   \code{TRUE} for flat rhizotron windows (e.g. glass-fronted boxes) where no
#'   sinusoidal correction is needed (\code{sinoid = FALSE}).  Default
#'   \code{FALSE}.
#'
# --- Core metrics (on by default) ---
#' @param calc_root_pixels Logical. Count \code{rootpx} (foreground pixels) and
#'   \code{voidpx} (background pixels) per depth bin.  Required for density
#'   metrics.  Default \code{TRUE}.
#' @param calc_root_length Logical. Estimate root length (cm) per depth bin
#'   from the skeleton using D8 connectivity (orthogonal steps = 1 pixel,
#'   diagonal steps = \eqn{\sqrt{2}} pixels, isolated pixels = 1 pixel).
#'   Required for length density and angle metrics.  Default \code{TRUE}.
#' @param calc_diameter_stats Logical. Compute per-bin mean, maximum, and
#'   variance of root diameter (cm) using the distance-transform approach in
#'   \code{root_diameter()}.  Default \code{TRUE}.
#'
# --- Extended metrics (off by default) ---
#' @param calc_diameter_quantiles Logical. Compute the diameter distribution
#'   percentiles set by \code{diameter_quantiles} per bin, conditional means
#'   above each quantile, threshold-based root lengths (see
#'   \code{diameter_thresholds}), and modal diameter peaks via
#'   \code{modal_peaks()}.  Default \code{FALSE}.
#' @param calc_landscape_metrics Logical. Compute patch-level landscape metrics
#'   per depth bin via \code{root_scape_metrics()}: nearest-neighbour distance
#'   (\code{enn_mn}), joint entropy (\code{joinent}), relative mutual
#'   information (\code{relmutinf}), number of patches (\code{np}), and
#'   contagion (\code{contag}).  \strong{Slow}: one call per depth bin per
#'   image.  Default \code{FALSE}.
#' @param calc_color_metrics Logical. Compute mean chromatic coordinates (rcc,
#'   gcc, bcc), hue, saturation, luminosity, and raw RGB channel means
#'   separately for root pixels and background pixels via
#'   \code{tube_coloration()}.  Requires \code{path.rgb}.  Default \code{FALSE}.
#' @param calc_root_angles Logical. Compute \code{deep_drive} (fraction of
#'   skeleton pixels whose D8 flow direction matches the locally optimal
#'   downward direction) and \code{mean.steepness.angle} /
#'   \code{sd.steepness.angle} (degrees, 0 = horizontal, 90 = vertical).
#'   Uses \code{deep_drive()}.  Default \code{FALSE}.
#' @param calc_root_order_metrics Logical. Build a per-image branching-order
#'   graph via \code{branch_order_map()} and summarise it both per depth bin
#'   and per tube.  Adds \code{mean.branch_order}, \code{max.branch_order},
#'   \code{mean.root_order}, and \code{lateral_root_fraction} per depth bin,
#'   plus tube-level \code{main_root.*} / \code{lateral_roots.*} columns
#'   (length, diameter, branching frequency, etc., split by
#'   \code{order_metrics(..., focal = "thickest")}) and \code{n_root_orders}
#'   (the highest branch order found).  Requires a skeleton.  \strong{Slow}:
#'   builds one segment graph per image.  Default \code{FALSE}.
#'
# --- Derived metrics ---
#' @param calc_density_metrics Logical. Compute \code{rootpx.density} (percent
#'   root area cover) and \code{rootlength.density} (cm root length per cm^2
#'   imaged area) per bin.  Auto-enables \code{calc_root_pixels} and
#'   \code{calc_root_length}.  Default \code{TRUE}.
#' @param calc_distribution_indices Logical. Compute tube-level indices:
#'   \code{mrd} (mean rooting depth), \code{rpi} (root
#'   penetration index), and \code{total.length.density} (summed length
#'   density over all bins, in cm root per cm^2 per cm depth).  Auto-enables
#'   \code{calc_density_metrics}.  Default \code{TRUE}.
#' @param calc_advanced_metrics Logical. Compute per-bin derived metrics:
#'   \code{rootlength.fraction} (each bin's length density as a fraction of the
#'   tube total), \code{mean.var.diameter} (mean of within-bin diameter
#'   variance), and \code{rootsurface_rootvolume_ratio} (lateral surface area
#'   over cylinder volume, summed over skeleton pixels in the bin and expressed
#'   as cm^2 per cm^3).  Auto-enables \code{calc_distribution_indices} and
#'   \code{calc_diameter_stats}.  Default \code{TRUE}.
#'
# --- Diameter threshold settings (used when \code{calc_diameter_quantiles = TRUE}) ---
#' @param diameter_thresholds Numeric vector. Diameter cut-offs for computing
#'   \code{rootlength.above.*} and \code{avg.diameter.above.*} columns.
#'   Units are set by \code{diameter_threshold_unit}.  Default \code{c(0.2, 0.5, 1)}.
#' @param diameter_threshold_unit Character. Unit of \code{diameter_thresholds}:
#'   \code{"mm"} (default), \code{"cm"}, or \code{"px"}.
#' @param diameter_quantiles Numeric vector of probabilities (each strictly
#'   between 0 and 1) for the per-bin diameter percentiles computed when
#'   \code{calc_diameter_quantiles = TRUE}.  Default \code{c(0.90, 0.95, 0.99)}.
#'   Output columns are named from the probabilities: e.g. \code{0.90} gives
#'   \code{rootdiameter.90} (the 90th percentile) and
#'   \code{avg.diameter.top10pct} (mean diameter above it).
#'
# --- Output ---
#' @param output_path Character or \code{NULL}. If provided, the result is saved
#'   as an \code{.RData} file at this path.  The exported object is named
#'   \code{root.depth.metrics}.  Parent directories are created if they do not
#'   exist.  Default \code{NULL} (no file written).
#' @param verbose Logical. Print per-image progress lines showing image index,
#'   per-image time, cumulative elapsed time, estimated remaining time, and
#'   predicted clock-time of completion.  Default \code{TRUE}.
#'
#' @return A data frame with one row per tube x depth-bin combination.  Always
#'   present columns: \code{Tube}, \code{Session}, \code{Plot}, \code{depth}.
#'   All other columns depend on which metric groups are enabled; disabled or
#'   failed metrics appear as \code{NA} rather than being absent.  Returns
#'   \code{NULL} invisibly if every image failed.
#'
#' @details
#' \strong{Surface-to-volume ratio.}  For each skeleton pixel, the root segment
#' is modelled as a cylinder of length \eqn{l_i} (the D8 path length of that
#' pixel in cm) and radius \eqn{r_i} (half the local diameter in cm).  The
#' lateral surface area is \eqn{2 \pi r_i l_i} and the volume is
#' \eqn{\pi r_i^2 l_i}.  Their ratio simplifies to \eqn{2 / r_i}.
#' \code{rootsurface_rootvolume_ratio} is the length-weighted mean of
#' \eqn{2 / r_i} over all skeleton pixels in the depth bin, in units of
#' cm\eqn{^{-1}} (cm^2 surface per cm^3 volume).  Thicker roots have a smaller
#' ratio; fine roots have a larger ratio.
#'
#' \strong{Fault tolerance.}  Every metric block is wrapped in
#' \code{tryCatch}.  Failures produce a \code{[RootScanR] SKIPPED} message and
#' \code{NA} values; they never abort the run.
#'
#' \strong{Dependency resolution.}  Enabling a higher-level metric silently
#' enables its prerequisites and prints a message listing what was auto-enabled.
#'
#' @examples
#' \dontrun{
#' # Minimal -- fast default metrics only
#' result <- root_depth_metrics(
#'   path.seg         = "scans/segmented/2022_02/",
#'   path.skl         = "scans/skeleton/2022_02/",
#'   insertion_angles = tube_meta$angle,
#'   session          = "2022_02"
#' )
#'
#' # With diameter quantiles and root angle metrics
#' result <- root_depth_metrics(
#'   path.seg                = "scans/segmented/2022_02/",
#'   path.skl                = "scans/skeleton/2022_02/",
#'   path.rgb                = "scans/blended/2022_02/",
#'   rgb_file_index          = 37:72,
#'   insertion_angles        = tube_meta$angle,
#'   soil_starts             = tube_meta$soil_row,
#'   session                 = "2022_02",
#'   calc_diameter_quantiles = TRUE,
#'   calc_root_angles        = TRUE,
#'   calc_color_metrics      = TRUE,
#'   diameter_thresholds     = c(0.2, 0.5, 1),
#'   diameter_threshold_unit = "mm",
#'   output_path             = "output/root_metrics_2022_02.RData"
#' )
#'
#' # Flat rhizotron window (no sinusoidal tube correction)
#' result <- root_depth_metrics(
#'   path.seg      = "scans/segmented/rhizotron_A/",
#'   path.skl      = "scans/skeleton/rhizotron_A/",
#'   flat_geometry = TRUE
#' )
#' }
#'
#' @importFrom dplyr group_by filter summarise full_join across rename_with
#'   select mutate cur_data
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_sub
#' @importFrom terra rast zonal values ext focal terrain subst flip t trim
#'   resample mask compareGeom levels
#'
#' @export
root_depth_metrics <- function(
    
  # ---------- paths ----------------------------------------------------------
  path.seg,
  path.skl                = NULL,
  path.rgb                = NULL,
  seg_file_index          = NULL,
  skl_file_index          = NULL,
  rgb_file_index          = NULL,
  
  # ---------- per-image metadata ---------------------------------------------
  insertion_angles        = 0,
  soil_starts             = 0,
  tube_names              = NULL,
  session                 = "",
  
  # ---------- scan / geometry ------------------------------------------------
  dpi                     = 300,
  tube_diameter_cm        = 7,
  depth_interval_cm       = 5,
  flat_geometry           = FALSE,
  
  # ---------- core metrics (on by default) -----------------------------------
  calc_root_pixels        = TRUE,
  calc_root_length        = TRUE,
  calc_diameter_stats     = TRUE,
  
  # ---------- extended metrics (off by default) ------------------------------
  calc_diameter_quantiles = FALSE,
  calc_landscape_metrics  = FALSE,
  calc_color_metrics      = FALSE,
  calc_root_angles        = FALSE,
  calc_root_order_metrics = FALSE,
  
  # ---------- derived metrics ------------------------------------------------
  calc_density_metrics       = TRUE,
  calc_distribution_indices  = TRUE,
  calc_advanced_metrics      = TRUE,
  
  # ---------- diameter threshold settings ------------------------------------
  diameter_thresholds      = c(0.2, 0.5, 1),
  diameter_threshold_unit  = "mm",
  diameter_quantiles       = c(0.90, 0.95, 0.99),
  
  # ---------- output ---------------------------------------------------------
  output_path             = NULL,
  verbose                 = TRUE
  
) {
  
  # ===========================================================================
  # 0.  Internal helpers (not exported)
  # ===========================================================================
  
  .msg <- function(...) if (verbose) message(sprintf(...))
  
  # Run expr safely; on error emit message and return fallback
  .safe <- function(label, expr, fallback = NULL) {
    tryCatch(expr,
             error = function(e) {
               message(sprintf("[RootScanR] SKIPPED '%s': %s", label, conditionMessage(e)))
               fallback
             }
    )
  }
  
  # Recycle x to length n, or stop if length is wrong
  .recycle <- function(x, n, name) {
    if (length(x) == 1L) return(rep(x, n))
    if (length(x) == n)  return(x)
    stop(sprintf("'%s' must have length 1 or %d (one per image), got %d.",
                 name, n, length(x)), call. = FALSE)
  }
  
  # ===========================================================================
  # 1.  Resolve file lists
  # ===========================================================================
  if (!dir.exists(path.seg))
    stop("'path.seg' does not exist: ", path.seg, call. = FALSE)
  
  all_seg <- list.files(path.seg)
  im.ls.seg <- if (!is.null(seg_file_index)) all_seg[seg_file_index] else all_seg
  if (length(im.ls.seg) == 0)
    stop("No files found in 'path.seg': ", path.seg, call. = FALSE)
  n_images <- length(im.ls.seg)
  
  im.ls.skl <- NULL
  if (!is.null(path.skl) && dir.exists(path.skl)) {
    all_skl <- list.files(path.skl)
    im.ls.skl <- if (!is.null(skl_file_index)) all_skl[skl_file_index] else all_skl
  }
  
  im.ls.rgb <- NULL
  if (!is.null(path.rgb) && dir.exists(path.rgb)) {
    all_rgb <- list.files(path.rgb)
    im.ls.rgb <- if (!is.null(rgb_file_index)) all_rgb[rgb_file_index] else all_rgb
  }
  
  # Recycle per-image metadata
  insertion_angles <- .recycle(insertion_angles, n_images, "insertion_angles")
  soil_starts      <- .recycle(soil_starts,      n_images, "soil_starts")
  if (!is.null(tube_names))
    tube_names <- .recycle(tube_names, n_images, "tube_names")
  
  # ===========================================================================
  # 2.  Global dependency resolution
  # ===========================================================================
  if (calc_root_angles && !calc_root_length) {
    message("[RootScanR] Auto-enabling calc_root_length (required for calc_root_angles).")
    calc_root_length <- TRUE
  }
  if (calc_diameter_quantiles && !calc_root_length) {
    message("[RootScanR] Auto-enabling calc_root_length (required for calc_diameter_quantiles).")
    calc_root_length <- TRUE
  }
  if (calc_root_order_metrics && !calc_root_length) {
    message("[RootScanR] Auto-enabling calc_root_length (required for calc_root_order_metrics).")
    calc_root_length <- TRUE
  }
  if (calc_density_metrics) {
    if (!calc_root_pixels) {
      message("[RootScanR] Auto-enabling calc_root_pixels (required for calc_density_metrics).")
      calc_root_pixels <- TRUE
    }
    if (!calc_root_length) {
      message("[RootScanR] Auto-enabling calc_root_length (required for calc_density_metrics).")
      calc_root_length <- TRUE
    }
  }
  if (calc_distribution_indices && !calc_density_metrics) {
    message("[RootScanR] Auto-enabling calc_density_metrics (required for calc_distribution_indices).")
    calc_density_metrics <- TRUE
  }
  if (calc_advanced_metrics) {
    if (!calc_density_metrics) {
      message("[RootScanR] Auto-enabling calc_density_metrics (required for calc_advanced_metrics).")
      calc_density_metrics <- TRUE
    }
    if (!calc_distribution_indices) {
      message("[RootScanR] Auto-enabling calc_distribution_indices (required for calc_advanced_metrics).")
      calc_distribution_indices <- TRUE
    }
    if (!calc_diameter_stats) {
      message("[RootScanR] Auto-enabling calc_diameter_stats (required for rootsurface_rootvolume_ratio).")
      calc_diameter_stats <- TRUE
    }
  }
  
  # Warn early about missing paths
  needs_skl <- calc_root_length || calc_diameter_stats || calc_diameter_quantiles ||
    calc_root_angles || calc_root_order_metrics
  if (needs_skl && is.null(im.ls.skl)) {
    message(paste(
      "[RootScanR] Skeleton directory (path.skl) not found or not supplied.",
      "Skeletons will be computed internally per image via skeletonize_image()."
    ))
  }
  if (calc_color_metrics && is.null(im.ls.rgb)) {
    warning("RGB directory (path.rgb) not found or not supplied. Disabling calc_color_metrics.",
            call. = FALSE)
    calc_color_metrics <- FALSE
  }
  
  # Diameter threshold conversion to cm
  thr_cm <- switch(diameter_threshold_unit,
                   "cm" = diameter_thresholds,
                   "mm" = diameter_thresholds / 10,
                   "px" = diameter_thresholds / (dpi / 2.54),
                   stop("'diameter_threshold_unit' must be 'mm', 'cm', or 'px'.", call. = FALSE)
  )
  thr_names <- paste0(diameter_thresholds, diameter_threshold_unit)

  # Diameter quantile column names, derived from `diameter_quantiles` so the
  # output stays in step when the probabilities are customised. Defaults
  # c(0.90, 0.95, 0.99) reproduce the historical
  # rootdiameter.{90,95,99} / avg.diameter.top{10,5,1}pct columns.
  if (!is.numeric(diameter_quantiles) || length(diameter_quantiles) == 0 ||
      any(is.na(diameter_quantiles)) ||
      any(diameter_quantiles <= 0 | diameter_quantiles >= 1)) {
    stop("'diameter_quantiles' must be numeric probabilities strictly between 0 and 1.",
         call. = FALSE)
  }
  .fmt_pct   <- function(x) format(round(x, 6), trim = TRUE, scientific = FALSE)
  q_diam_names <- paste0("rootdiameter.", .fmt_pct(diameter_quantiles * 100))
  q_top_names  <- paste0("avg.diameter.top", .fmt_pct((1 - diameter_quantiles) * 100), "pct")

  # ===========================================================================
  # 3.  Main image loop
  # ===========================================================================
  root.list   <- vector("list", n_images)
  failed_imgs <- character(0)
  img_times   <- numeric(0)            # seconds per image, for rolling ETA
  t_start     <- proc.time()[["elapsed"]]
  
  for (l in seq_len(n_images)) {
    
    t_img <- proc.time()[["elapsed"]]
    
    seg_file <- im.ls.seg[l]
    tube     <- if (!is.null(tube_names)) tube_names[l] else
      paste0("T", stringr::str_sub(seg_file, start = -8, end = -6))
    angle    <- insertion_angles[l]
    soil0    <- soil_starts[l]
    
    # Local flag copies -- degraded per image without affecting other images
    do_pixels    <- calc_root_pixels
    do_length    <- calc_root_length
    do_diam_st   <- calc_diameter_stats
    do_diam_q    <- calc_diameter_quantiles
    do_landscape <- calc_landscape_metrics
    do_color     <- calc_color_metrics
    do_angles    <- calc_root_angles
    do_order     <- calc_root_order_metrics
    do_density   <- calc_density_metrics
    
    # -------------------------------------------------------------------------
    # 3a. Load images
    # -------------------------------------------------------------------------
    im <- .safe(sprintf("load segmented [%s]", seg_file), {
      img <- load_flexible_image(paste0(path.seg, seg_file),
                                 output_format = "spatrast",
                                 scale = "binary")
      if (dim(img)[3] > 3) img <- img[[1:3]]
      img
    })
    if (is.null(im)) {
      message(sprintf("[RootScanR] [%d/%d] %s: could not load segmented image -- skipping.",
                      l, n_images, tube))
      failed_imgs <- c(failed_imgs, seg_file)
      img_times   <- c(img_times, proc.time()[["elapsed"]] - t_img)
      next
    }
    
    im.skeleton <- NULL
    if (do_length || do_diam_st || do_diam_q || do_angles || do_order) {
      if (!is.null(im.ls.skl) && l <= length(im.ls.skl)) {
        im.skeleton <- .safe(sprintf("load skeleton [%s]", im.ls.skl[l]), {
          sk <- load_flexible_image(paste0(path.skl, im.ls.skl[l]),
                                    output_format = "spatrast",
                                    scale = "binary", select.layer = 2)
          if (dim(sk)[3] > 3) sk <- sk[[1:3]]
          sk
        })
      }
      if (is.null(im.skeleton)) {
        im.skeleton <- .safe(sprintf("skeletonize [%s]", seg_file), {
          skeletonize_image(im, verbose = FALSE)
        })
      }
      if (is.null(im.skeleton)) {
        message(sprintf("[RootScanR] %s: skeleton unavailable -- disabling length/diameter/angle/order metrics.", tube))
        do_length <- do_diam_st <- do_diam_q <- do_angles <- do_order <- FALSE
      }
    }
    
    im.rgb <- NULL
    if (do_color) {
      if (!is.null(im.ls.rgb) && l <= length(im.ls.rgb)) {
        im.rgb <- .safe(sprintf("load RGB [%s]", im.ls.rgb[l]),
                        terra::rast(paste0(path.rgb, im.ls.rgb[l])))
      }
      if (is.null(im.rgb)) {
        message(sprintf("[RootScanR] %s: RGB image unavailable -- disabling calc_color_metrics.", tube))
        do_color <- FALSE
      }
    }
    
    # -------------------------------------------------------------------------
    # 3b. Rotation censor (crop to tube interior)
    # -------------------------------------------------------------------------
    r0 <- round(dim(im)[1] / 2, 0)
    
    im <- .safe("rotation_censor (seg)",
                rotation_censor(im, center.offset = r0, fixed.rotation = TRUE, fixed.width = 1800),
                fallback = im)
    
    if (!is.null(im.skeleton))
      im.skeleton <- .safe("rotation_censor (skl)",
                           rotation_censor(im.skeleton, center.offset = r0, fixed.rotation = TRUE, fixed.width = 1800),
                           fallback = im.skeleton)
    
    if (!is.null(im.rgb))
      im.rgb <- .safe("rotation_censor (rgb)",
                      rotation_censor(im.rgb, center.offset = r0, fixed.rotation = TRUE, fixed.width = 1800),
                      fallback = im.rgb)
    
    # Keep only the segmentation layer; align extents
    im <- im[[2]]
    if (!is.null(im.skeleton)) terra::ext(im.skeleton) <- terra::ext(im)
    if (!is.null(im.rgb))      terra::ext(im.rgb)      <- terra::ext(im)
    
    # -------------------------------------------------------------------------
    # 3c. Depth map and binning
    # -------------------------------------------------------------------------
    DepthMap <- .safe("create_depthmap", {
      dm <- create_depthmap(
        img         = im,
        sinoid      = !flat_geometry,
        dpi         = dpi,
        start.soil  = soil0,
        center.offset = 0,
        tilt        = angle,
        tube.thicc  = tube_diameter_cm
      )
      dm <- terra::flip(terra::t(dm))
      terra::ext(dm) <- terra::ext(im)
      dm
    })
    if (is.null(DepthMap)) {
      message(sprintf("[RootScanR] [%d/%d] %s: create_depthmap failed -- skipping.", l, n_images, tube))
      failed_imgs <- c(failed_imgs, seg_file)
      img_times   <- c(img_times, proc.time()[["elapsed"]] - t_img)
      next
    }
    
    bm    <- binning(depthmap = DepthMap, nn = depth_interval_cm, round.option = "rounding")
    roots <- data.frame(depth = sort(unique(terra::values(bm))))
    
    # -------------------------------------------------------------------------
    # 3d. Basic pixel counts
    # -------------------------------------------------------------------------
    if (do_pixels) {
      pd <- .safe("pixel counts", {
        rp  <- terra::zonal(im, bm, "sum", na.rm = TRUE)
        vd  <- im; terra::values(vd) <- 1 - terra::values(im)
        vp  <- terra::zonal(vd, bm, "sum", na.rm = TRUE)
        out <- merge(rp, vp, by = names(rp)[1])
        colnames(out) <- c("depth", "rootpx", "voidpx")
        out
      })
      if (!is.null(pd)) {
        roots <- merge(roots, pd, by = "depth", all.x = TRUE)
      } else {
        roots$rootpx <- roots$voidpx <- NA_real_
        do_pixels <- FALSE
      }
    }
    
    # -------------------------------------------------------------------------
    # 3e. Root length from skeleton (D8 path lengths)
    # -------------------------------------------------------------------------
    root.length.map <- NULL
    angles_map      <- NULL
    
    if (do_length) {
      rl_res <- .safe("root length", {
        
        # Flat depth map (sinoid = FALSE) used to derive flow directions for
        # length calculation -- the sinusoidal correction is for the depth axis
        # only, not for path-length geometry.
        dm_flat <- create_depthmap(
          img           = im,
          sinoid        = FALSE,
          dpi           = dpi,
          start.soil    = soil0,
          center.offset = 0.5,
          tilt          = angle,
          tube.thicc    = tube_diameter_cm
        )
        dem_flat <- terra::flip(terra::t(dm_flat))
        terra::ext(dem_flat) <- terra::ext(im)
        dem_flat[im.skeleton != 1] <- NA
        
        ang <- terra::terrain(dem_flat, v = "flowdir")
        ang <- terra::subst(ang,
                            from = c(  0,  1,   2,   4,   8,  16,  32, 64, 128),
                            to   = c( NA, 90, 135, 180, 225, 270, 315,  0,  45))
        orth <- (ang ==   0 | ang ==  90 | ang == 180 | ang == 270)
        diag <- (ang ==  45 | ang == 135 | ang == 225 | ang == 315)
        orth <- terra::t(terra::flip(orth))
        diag <- terra::t(terra::flip(diag))
        
        w        <- matrix(1, 3, 3)
        nb_sum   <- terra::focal(im.skeleton, w, fun = sum, na.policy = "omit")
        isolated <- (nb_sum == 1)
        
        rlm <- orth * 1 + diag * sqrt(2) + terra::t(terra::flip(isolated))
        rlm <- terra::flip(terra::t(rlm))
        
        rl_z <- terra::zonal(rlm, bm, "sum", na.rm = TRUE)
        rl_z[[2]] <- rl_z[[2]] / (dpi / 2.54)   # pixels -> cm
        colnames(rl_z) <- c("depth", "rootlength")
        
        list(rootlength = rl_z, rlm = rlm, ang = ang)
      })
      
      if (!is.null(rl_res)) {
        roots           <- merge(roots, rl_res$rootlength, by = "depth", all.x = TRUE)
        root.length.map <- rl_res$rlm
        angles_map      <- rl_res$ang
      } else {
        roots$rootlength <- NA_real_
        do_length <- do_angles <- do_diam_q <- FALSE
      }
    }
    
    # -------------------------------------------------------------------------
    # 3f. Root diameter map (shared between stats and quantile blocks)
    # -------------------------------------------------------------------------
    rd.map <- NULL
    if (do_diam_st || do_diam_q) {
      rd.map <- .safe("root_diameter map", {
        dm <- root_diameter(im, skeleton.img = im.skeleton, unit = "cm")$diameter_rast
        terra::ext(dm) <- terra::ext(bm)
        dm
      })
      if (is.null(rd.map)) {
        message(sprintf("[RootScanR] %s: diameter map failed -- disabling diameter metrics.", tube))
        do_diam_st <- do_diam_q <- FALSE
      }
    }
    
    # Diameter summary statistics (mean, max, variance per bin)
    if (do_diam_st && !is.null(rd.map)) {
      diam <- .safe("diameter stats", {
        avg <- terra::zonal(rd.map, bm, "mean", na.rm = TRUE); colnames(avg) <- c("depth", "avg.diameter")
        mx  <- terra::zonal(rd.map, bm, "max",  na.rm = TRUE); colnames(mx)  <- c("depth", "max.diameter")
        vr  <- terra::zonal(rd.map, bm, "var",  na.rm = TRUE); colnames(vr)  <- c("depth", "var.diameter")
        Reduce(function(a, b) merge(a, b, by = "depth"), list(avg, mx, vr))
      })
      if (!is.null(diam)) {
        roots <- merge(roots, diam, by = "depth", all.x = TRUE)
      } else {
        roots[c("avg.diameter", "max.diameter", "var.diameter")] <- NA_real_
      }
    }
    
    # -------------------------------------------------------------------------
    # 3f2. Root branching-order metrics (per bin + per tube)
    # -------------------------------------------------------------------------
    if (do_order && !is.null(im.skeleton)) {
      ord_res <- .safe("root order metrics", {

        bo <- branch_order_map(skel = im.skeleton, mask = im, order = "branch_order",
                               unit = "cm", dpi = dpi, return_map = TRUE,
                               template = im.skeleton, verbose = FALSE)
        et <- bo$edges

        bo_map <- bo$class_map
        terra::ext(bo_map) <- terra::ext(bm)
        ro_map <- order_classification_map(et, im.skeleton, value = "root_order")
        terra::ext(ro_map) <- terra::ext(bm)

        bo_mean <- terra::zonal(bo_map, bm, "mean", na.rm = TRUE); colnames(bo_mean) <- c("depth", "mean.branch_order")
        bo_max  <- terra::zonal(bo_map, bm, "max",  na.rm = TRUE); colnames(bo_max)  <- c("depth", "max.branch_order")
        ro_mean <- terra::zonal(ro_map, bm, "mean", na.rm = TRUE); colnames(ro_mean) <- c("depth", "mean.root_order")

        lateral_px <- (bo_map > 1) * 1
        ordered_px <- (!is.na(bo_map)) * 1
        lat_z <- terra::zonal(lateral_px, bm, "sum", na.rm = TRUE); colnames(lat_z) <- c("depth", "lateral_px")
        tot_z <- terra::zonal(ordered_px, bm, "sum", na.rm = TRUE); colnames(tot_z) <- c("depth", "ordered_px")

        per_bin <- Reduce(function(a, b) merge(a, b, by = "depth", all = TRUE),
                          list(bo_mean, bo_max, ro_mean, lat_z, tot_z))
        per_bin$lateral_root_fraction <- ifelse(per_bin$ordered_px > 0,
                                                 per_bin$lateral_px / per_bin$ordered_px, NA_real_)
        per_bin$lateral_px <- per_bin$ordered_px <- NULL

        # Tube-level main-root vs lateral-root summary (thickest order = main root)
        om <- order_metrics(bo, focal = "thickest")
        metric_cols <- c("n_segments", "n_tips", "n_branch_points", "total_length",
                         "length_fraction", "mean_segment_length", "branching_frequency",
                         "mean_diameter", "median_diameter")
        grp_map <- c(focal = "main_root", rest = "lateral_roots")
        tube_row <- as.list(stats::setNames(rep(NA_real_, length(metric_cols) * length(grp_map)),
                                     as.vector(outer(grp_map, metric_cols, paste, sep = "."))))
        for (i in seq_len(nrow(om))) {
          g <- grp_map[[om$group[i]]]
          for (m in metric_cols) tube_row[[paste0(g, ".", m)]] <- om[[m]][i]
        }
        tube_row$n_root_orders <- suppressWarnings(max(et$branch_order, na.rm = TRUE))
        if (!is.finite(tube_row$n_root_orders)) tube_row$n_root_orders <- NA_real_

        list(per_bin = per_bin, tube = tube_row)
      })

      if (!is.null(ord_res)) {
        roots <- merge(roots, ord_res$per_bin, by = "depth", all.x = TRUE)
        for (nm in names(ord_res$tube)) roots[[nm]] <- ord_res$tube[[nm]]
      } else {
        roots[c("mean.branch_order", "max.branch_order", "mean.root_order",
               "lateral_root_fraction")] <- NA_real_
      }
    }

    roots$Tube <- tube

    # -------------------------------------------------------------------------
    # 3g. Density metrics (Level 1 derived)
    # -------------------------------------------------------------------------
    if (do_density) {
      .safe("density metrics", {
        if (do_pixels && all(c("rootpx", "voidpx") %in% names(roots))) {
          roots$rootpx.density <<- roots$rootpx / (roots$rootpx + roots$voidpx) * 100
        }
        if (do_length && do_pixels && all(c("rootlength", "rootpx", "voidpx") %in% names(roots))) {
          # rootlength.density: cm root length per cm^2 of imaged area per bin
          roots$rootlength.density <<-
            roots$rootlength / ((roots$rootpx + roots$voidpx) / (dpi / 2.54)^2)
        }
      })
    }
    
    # -------------------------------------------------------------------------
    # 3h. Per-depth-slice loop
    #     (landscape metrics, colour metrics, diameter quantiles, root angles)
    # -------------------------------------------------------------------------
    needs_slice <- do_landscape || do_color || do_diam_q || do_angles
    
    if (needs_slice) {
      
      depth.slices <- sort(unique(terra::values(bm)))
      
      lsm_names <- if (do_landscape)
        c("lsm_c_enn_mn", "lsm_l_joinent", "lsm_l_relmutinf", "lsm_l_np", "lsm_l_contag")
      else character(0)
      
      # Resample RGB to segmented extent once per image
      im.rgb.crop <- NULL
      if (do_color && !is.null(im.rgb))
        im.rgb.crop <- .safe("RGB resample", terra::resample(im.rgb, im, "bilinear"))
      
      # Pre-build deep_drive optimal-angle map once per image, reuse per slice
      bm_vals <- ang_vals <- gg_vals <- NULL
      
      if (do_angles && !is.null(angles_map)) {
        dd <- .safe("deep_drive (full image)", {
          adm <- terra::t(terra::flip(DepthMap))
          terra::ext(adm) <- terra::ext(angles_map)
          RootScanR::deep_drive(DepthMap = adm, AngleMap = angles_map, return = "all")
        })
        gg.full <- if (!is.null(dd)) dd$optimal_angle_map else NULL
        
        if (!is.null(gg.full)) {
          # Fix geometry if extents differ but dimensions match
          if (!terra::compareGeom(gg.full, bm, stopOnError = FALSE)) {
            if (all(dim(gg.full) == dim(bm))) {
              terra::ext(gg.full) <- terra::ext(bm)
            } else {
              message(sprintf("[RootScanR] %s: gg.full/bm dimension mismatch -- disabling calc_root_angles.", tube))
              gg.full   <- NULL
              do_angles <- FALSE
            }
          }
        } else {
          do_angles <- FALSE
        }
        
        if (do_angles) {
          bm_vals  <- terra::values(bm)
          ang_vals <- terra::values(angles_map)
          gg_vals  <- terra::values(gg.full)
        }
      }
      
      acc <- data.frame(depth = numeric(0))   # accumulator across slices
      
      for (d in depth.slices) {
        
        im.sl <- im; im.sl[bm != d] <- NA; im.sl <- terra::trim(im.sl)
        
        # --- Landscape metrics -----------------------------------------------
        base_row <- if (do_landscape && length(lsm_names) > 0) {
          .safe(sprintf("landscape (depth=%g)", d), {
            rs   <- root_scape_metrics(im.sl, metrics = lsm_names)
            rs$depth <- d
            wide <- tidyr::pivot_wider(rs, names_from = "metric", values_from = "value")
            wide$object <- NULL
            for (m in setdiff(c("depth", stringr::str_sub(lsm_names, start = 7)), names(wide)))
              wide[[m]] <- NA_real_
            wide
          }, fallback = data.frame(depth = d))
        } else {
          data.frame(depth = d)
        }
        
        # --- Colour metrics --------------------------------------------------
        rc <- pc <- data.frame()
        empty_col <- data.frame(rcc = NA, gcc = NA, bcc = NA, hue = NA,
                                saturation = NA, luminosity = NA,
                                red = NA, green = NA, blue = NA)
        if (do_color && !is.null(im.rgb.crop)) {
          cr <- .safe(sprintf("colour (depth=%g)", d), {
            sl  <- im.rgb.crop; sl[bm != d] <- NA; sl <- terra::trim(sl)
            ri  <- sl; ri[im.sl == 0] <- NA   # root pixels only
            pi_ <- sl; pi_[im.sl == 1] <- NA  # background pixels only
            rc_ <- tryCatch(RootScanR::tube_coloration(ri),  error = function(e) empty_col)
            pc_ <- tryCatch(RootScanR::tube_coloration(pi_), error = function(e) empty_col)
            list(
              rc = dplyr::rename_with(rc_, ~ paste0(.x, "_root")),
              pc = dplyr::rename_with(pc_, ~ paste0(.x, "_bg"))
            )
          })
          if (!is.null(cr)) { rc <- cr$rc; pc <- cr$pc }
        }
        
        # --- Diameter quantiles ----------------------------------------------
        q_row <- data.frame()
        if (do_diam_q && !is.null(rd.map) && !is.null(root.length.map)) {
          q_row <- .safe(sprintf("diameter quantiles (depth=%g)", d), {
            sl_rd  <- rd.map; sl_rd[bm != d] <- NA; sl_rd <- terra::trim(sl_rd)
            rd_v   <- terra::values(sl_rd, na.rm = FALSE)
            
            qv <- stats::quantile(rd_v, diameter_quantiles, na.rm = TRUE)

            res <- data.frame(depth = d)
            for (qi in seq_along(diameter_quantiles)) {
              res[[q_diam_names[qi]]] <- as.vector(qv[qi])
              res[[q_top_names[qi]]]  <- mean(rd_v[rd_v > qv[qi]], na.rm = TRUE)
            }
            
            sl_rl <- root.length.map; sl_rl[bm != d] <- NA
            for (ti in seq_along(thr_cm)) {
              ab  <- sl_rl; ab[rd.map < thr_cm[ti]] <- NA
              res[[paste0("rootlength.above.",    thr_names[ti])]] <-
                sum(terra::values(ab), na.rm = TRUE) / (dpi / 2.54)
              res[[paste0("avg.diameter.above.",  thr_names[ti])]] <-
                mean(rd_v[!is.na(rd_v) & rd_v > thr_cm[ti]], na.rm = TRUE)
            }
            
            rd_clean <- rd_v[!is.na(rd_v) & rd_v > 0]
            pk <- tryCatch({
              mp <- RootScanR::modal_peaks(rd_clean, display_type = "none",
                                           prominence_threshold = length(rd_clean) / sqrt(length(rd_clean)),
                                           mclust = FALSE)
              np <- length(mp$peak_x)
              data.frame(
                n.diameter.peaks = np,
                diameter.peak.1  = if (np >= 1) mp$peak_x[1] else NA_real_,
                diameter.peak.2  = if (np >= 2) mp$peak_x[2] else NA_real_,
                diameter.peak.3  = if (np >= 3) mp$peak_x[3] else NA_real_
              )
            }, error = function(e)
              data.frame(n.diameter.peaks = NA_real_, diameter.peak.1 = NA_real_,
                         diameter.peak.2  = NA_real_, diameter.peak.3  = NA_real_))
            
            cbind(res, pk)
          }, fallback = data.frame(depth = d))
        }
        
        # --- Root angles (deep_drive + steepness) ----------------------------
        ang_row <- data.frame()
        if (do_angles && !is.null(gg_vals)) {
          ang_row <- .safe(sprintf("root angles (depth=%g)", d), {
            idx   <- bm_vals == d
            valid <- idx & !is.na(ang_vals) & !is.na(gg_vals)
            if (sum(valid) == 0L) {
              data.frame(depth = d, deep_drive = NA_real_,
                         mean.steepness.angle = NA_real_, sd.steepness.angle = NA_real_)
            } else {
              as_  <- ang_vals[valid]
              gg_  <- gg_vals[valid]
              dv   <- mean(as_ == gg_)
              dev  <- (as_ %% 180) - (gg_ %% 180)
              st   <- abs(sin(dev * pi / 180)) * 90
              data.frame(depth = d, deep_drive = dv,
                         mean.steepness.angle = mean(st, na.rm = TRUE),
                         sd.steepness.angle   = stats::sd(st,   na.rm = TRUE))
            }
          }, fallback = data.frame(depth = d, deep_drive = NA_real_,
                                   mean.steepness.angle = NA_real_,
                                   sd.steepness.angle   = NA_real_))
        }
        
        # --- Merge all slice results -----------------------------------------
        row <- base_row
        for (df in Filter(
          function(x) is.data.frame(x) && nrow(x) > 0 && "depth" %in% names(x),
          list(rc, pc, q_row, ang_row)
        )) row <- merge(row, df, by = "depth", all.x = TRUE, sort = FALSE)
        
        # Safe column-aligning rbind
        if (nrow(acc) == 0L) {
          acc <- row
        } else {
          for (m in setdiff(names(acc), names(row))) row[[m]] <- NA_real_
          for (m in setdiff(names(row), names(acc))) acc[[m]] <- NA_real_
          acc <- rbind(acc, row[names(acc)])
        }
        
      } # end depth-slice loop
      
      roots <- merge(roots, acc, by = "depth", all.x = TRUE)
      
      if (do_density && do_landscape && do_pixels && "np" %in% names(roots))
        roots$np_density <- roots$np /
        ((roots$rootpx + roots$voidpx) / (dpi / 2.54)^2)
      
    } # end needs_slice block
    
    root.list[[l]] <- roots
    
    # -------------------------------------------------------------------------
    # 3i. Progress + rolling ETA
    # -------------------------------------------------------------------------
    img_times <- c(img_times, proc.time()[["elapsed"]] - t_img)
    if (verbose) {
      n_done   <- length(img_times)
      n_left   <- n_images - n_done
      cum_secs <- sum(img_times)
      avg_secs <- mean(img_times)
      eta_secs <- avg_secs * n_left
      finish   <- Sys.time() + eta_secs
      message(sprintf(
        "[RootScanR] [%d/%d] %s | img: %.0fs | elapsed: %.0fs | remaining: ~%.0fs | done ~%s",
        l, n_images, tube,
        img_times[n_done],
        cum_secs,
        eta_secs,
        format(finish, "%H:%M")
      ))
    }
    
  } # end image loop
  
  # ===========================================================================
  # 4.  Post-processing
  # ===========================================================================
  valid <- Filter(Negate(is.null), root.list)
  
  if (length(failed_imgs) > 0)
    warning(sprintf("[RootScanR] %d image(s) skipped entirely: %s",
                    length(failed_imgs),
                    paste(failed_imgs, collapse = ", ")),
            call. = FALSE)
  
  if (length(valid) == 0L) {
    warning("[RootScanR] No images were processed successfully.", call. = FALSE)
    return(invisible(NULL))
  }
  
  root.frame <- array2DF(array(valid))
  root.frame$Var1    <- NULL
  root.frame$Plot    <- suppressWarnings(
    as.numeric(stringr::str_sub(root.frame$Tube, start = 3)))
  root.frame$Session <- session
  
  # NaN -> NA across all numeric columns
  for (col in names(root.frame)[sapply(root.frame, is.numeric)])
    root.frame[[col]] <- ifelse(is.nan(root.frame[[col]]), NA_real_, root.frame[[col]])
  
  # ---------------------------------------------------------------------------
  # Level 2: tube-level distribution indices
  # ---------------------------------------------------------------------------
  r <- root.frame
  if (calc_distribution_indices && calc_density_metrics &&
      "rootlength.density" %in% names(r)) {
    coag1 <- .safe("distribution indices", {
      r |>
        dplyr::group_by(dplyr::.data$Tube) |>
        dplyr::filter(!is.na(dplyr::.data$depth) & !is.na(dplyr::.data$rootlength.density)) |>
        dplyr::summarise(
          mrd   = RootScanR::MRD(w = dplyr::.data$depth, roots = dplyr::.data$rootlength.density),
          rpi   = RootScanR::RPI(w = dplyr::.data$depth, roots = dplyr::.data$rootlength.density),
          # total.length.density: sum of (length density x bin size) across all bins
          # units: cm root per cm^2 (integrated over the full depth profile)
          total.length.density = sum(dplyr::.data$rootlength.density * depth_interval_cm, na.rm = TRUE),
          .groups = "drop"
        )
    })
    if (!is.null(coag1)) r <- dplyr::full_join(r, coag1, by = "Tube")
  }

  # ---------------------------------------------------------------------------
  # Level 3: per-bin advanced metrics
  # ---------------------------------------------------------------------------
  rr <- r
  if (calc_advanced_metrics && calc_distribution_indices &&
      "rootlength.density" %in% names(r)) {

    coag2 <- .safe("advanced metrics", {

      r |>
        dplyr::group_by(dplyr::.data$Tube, dplyr::.data$depth) |>
        dplyr::filter(!is.na(dplyr::.data$depth) & !is.na(dplyr::.data$rootlength.density)) |>
        dplyr::summarise(

          # Fraction of the tube's total length density contributed by this bin
          rootlength.fraction = if (calc_density_metrics && "total.length.density" %in% names(dplyr::pick(dplyr::everything())))
            dplyr::.data$rootlength.density / dplyr::.data$total.length.density else NA_real_,

          # Joint entropy per unit root pixel density
          # (meaningful only when landscape metrics are available)
          ent_per_rootpx = if (calc_landscape_metrics && "joinent" %in% names(dplyr::pick(dplyr::everything())))
            dplyr::.data$joinent / dplyr::.data$rootpx.density else NA_real_,

          # Number of root patches per unit root pixel density
          patch_density_norm = if ("np_density" %in% names(dplyr::pick(dplyr::everything())))
            dplyr::.data$np_density / dplyr::.data$rootpx.density else NA_real_,

          # Mean within-bin diameter variance (averaged over any sub-grouping)
          mean.var.diameter = if (calc_diameter_stats && "var.diameter" %in% names(dplyr::pick(dplyr::everything())))
            mean(dplyr::.data$var.diameter, na.rm = TRUE) else NA_real_,

          # Surface-to-volume ratio (cm^-1).
          # Each skeleton pixel i contributes a cylinder of length l_i (cm)
          # and radius r_i = avg.diameter / 2.
          # Lateral surface area = 2*pi*r_i*l_i
          # Volume               =   pi*r_i^2*l_i
          # Ratio per pixel      = 2 / r_i   (l_i cancels)
          # Here we use the bin-mean diameter as a proxy for r_i.
          # TODO: ideally this should be computed pixel-by-pixel using the full
          # diameter raster and root-length-map before zonal aggregation, so
          # that the length-weighted mean of 2/r_i is returned rather than
          # 2 / mean(r_i). Flag for future revision.
          rootsurface_rootvolume_ratio = if (
            calc_diameter_stats && calc_density_metrics &&
            all(c("avg.diameter", "rootlength.density") %in%
                names(dplyr::pick(dplyr::everything()))) &&
            !is.na(dplyr::.data$avg.diameter) && dplyr::.data$avg.diameter > 0
          ) 2 / (dplyr::.data$avg.diameter / 2) else NA_real_,

          .groups = "drop"
        )
    })
    if (!is.null(coag2)) rr <- dplyr::full_join(r, coag2, by = c("Tube", "depth"))
  }
  
  root.depth.metrics <- rr
  
  # ===========================================================================
  # 5.  Optional file output
  # ===========================================================================
  if (!is.null(output_path)) {
    .safe("save output file", {
      out_dir <- dirname(output_path)
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      save(root.depth.metrics, file = output_path)
      .msg("[RootScanR] Results saved to: %s", output_path)
    })
  }
  
  if (verbose) {
    total_min <- round(sum(img_times) / 60, 1)
    n_ok      <- n_images - length(failed_imgs)
    message(sprintf(
      "[RootScanR] Done. %d/%d images processed in %.1f min.",
      n_ok, n_images, total_min
    ))
  }
  
  invisible(root.depth.metrics)
}
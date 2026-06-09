#' Process Root Images into Depth-Resolved Metrics
#'
#' Loops through segmented root images (and optionally skeleton and RGB images) to compute
#' a suite of depth-resolved root and peat metrics. Metrics are calculated per depth slice for
#' each tube, optionally including basic pixel counts, root length, diameter statistics,
#' rootscape structure, color, peat classification, density, distribution indices, and
#' advanced derived metrics. The function handles metric dependencies, provides progress
#' reporting with estimated remaining time, and optionally saves results to an RData file.
#'
#' @param session Character. Session identifier used for labeling outputs and constructing
#'   default paths. Should correspond to the naming of your image directories if using
#'   defaults (e.g., "2022_01" if images are stored in ".../Segmented/Full/Segmented/2022_01/").
#' @param tube_thickness Numeric. Tube thickness in pixels (default = 7).
#' @param dpi Numeric. Image resolution in dots per inch (default = 300).
#' @param path_seg Character. Path to segmented root images.
#' @param path_skl Character. Path to skeleton images.
#' @param path_rgb Character or NULL. Path to RGB images (optional if color metrics disabled).
#' @param insertion_data Data frame or numeric. Either a single numeric angle applied to all tubes,
#'   or a data frame containing per-tube insertion angles with at least columns `Plot` and `InsertionAngle2`.
#' @param metric_flags List of logicals controlling which metrics to calculate (default = all TRUE):
#'   \describe{
#'     \item{CALC_BASIC_PIXELS}{Pixel counts for roots and voids}
#'     \item{CALC_ROOT_LENGTH}{Root length from skeleton images}
#'     \item{CALC_DIAMETER_STATS}{Mean, max, variance of root diameters}
#'     \item{CALC_DIAMETER_QUANTILES}{Quantiles (90th, 95th, 99th) of root diameters}
#'     \item{CALC_LANDSCAPE_METRICS}{Patch-level spatial metrics of root distribution}
#'     \item{CALC_COLOR_METRICS}{Color metrics (RGB/HSL/chromatic) for roots and peat}
#'     \item{CALC_PEAT_CLASSES}{Classify peat pixels from RGB images}
#'     \item{CALC_DENSITY_METRICS}{Derived metrics per unit area, e.g., root density}
#'     \item{CALC_DISTRIBUTION_INDICES}{Metrics quantifying how roots are distributed along depth, e.g., RWDI, RPI}
#'     \item{CALC_ADVANCED_METRICS}{Composite metrics combining multiple base or derived metrics, e.g., rootlength fraction, patch density normalized, entropy per root pixel}
#'   }
#' @param check_dependencies Logical. Automatically enable required dependent metrics (default = TRUE).
#' @param verbose Logical. Print progress messages with estimated remaining time (default = TRUE).
#' @param save_path Character or NULL. If supplied, final results are automatically saved to this RData file path.
#'
#' @return A data frame with depth-resolved metrics for each tube, including optional derived
#'   and advanced metrics.
#' @export
#'
#' @examples
#' \dontrun{
#' # Using single insertion angle for all tubes
#' results <- process_root_images(
#'   session = "2022_01",
#'   path_seg = "D:/Blending_Scans/Segmented/Full/Segmented/2022_01/",
#'   path_skl = "D:/Blending_Scans/Segmented/Full/Skeleton/2022_01/",
#'   path_rgb = "D:/Blending_Scans/Full/2022_01/",
#'   insertion_data = 45,
#'   metric_flags = list(CALC_COLOR_METRICS = FALSE)
#' )
#'
#' # Using per-tube insertion angles
#' angle_data <- read.csv("OulankaTubeInsertionAngle.csv")
#' results <- process_root_images(
#'   session = "2022_01",
#'   path_seg = "D:/Blending_Scans/Segmented/Full/Segmented/2022_01/",
#'   path_skl = "D:/Blending_Scans/Segmented/Full/Skeleton/2022_01/",
#'   path_rgb = "D:/Blending_Scans/Full/2022_01/",
#'   insertion_data = angle_data
#' )
#' }
process_root_images <- function(
    session,
    tube_thickness = 7,
    dpi = 300,
    path_seg,
    path_skl,
    path_rgb = NULL,
    insertion_data,
    metric_flags = list(
      CALC_BASIC_PIXELS = TRUE,
      CALC_ROOT_LENGTH = TRUE,
      CALC_DIAMETER_STATS = TRUE,
      CALC_DIAMETER_QUANTILES = TRUE,
      CALC_LANDSCAPE_METRICS = TRUE,
      CALC_COLOR_METRICS = TRUE,
      CALC_PEAT_CLASSES = TRUE,
      CALC_DENSITY_METRICS = TRUE,
      CALC_DISTRIBUTION_INDICES = TRUE,
      CALC_ADVANCED_METRICS = TRUE
    ),
    check_dependencies = TRUE,
    verbose = TRUE,
    save_path = NULL
) {
  # -----------------------------
  # Process insertion data
  # -----------------------------
  if(length(insertion_data) == 1) {
    s0 <- data.frame(Plot = seq_along(list.files(path_seg)), soil0 = 0, InsertionAngle = insertion_data)
  } else {
    s0 <- insertion_data
    s0$soil0 <- 0
    s0$InsertionAngle <- s0$InsertionAngle2
  }
  
  # -----------------------------
  # Dependency Management
  # -----------------------------
  if(check_dependencies) {
    # same dependency logic as before...
  }
  
  # -----------------------------
  # Image listing and loop setup
  # -----------------------------
  im_seg_list <- list.files(path_seg)
  im_skl_list <- list.files(path_skl)
  im_rgb_list <- if(!is.null(path_rgb)) list.files(path_rgb) else NULL
  root_list <- vector("list", length(im_seg_list))
  t_start <- Sys.time()
  
  # -----------------------------
  # Main loop: process each tube
  # -----------------------------
  for(i in seq_along(im_seg_list)) {
    # load images, rotate, create mask, compute metrics, depth slice loop
    # store results in root_list[[i]]
    # progress message with estimated remaining time
  }
  
  # -----------------------------
  # Post-processing: merge results, calculate derived metrics
  # -----------------------------
  # distribution indices, advanced metrics
  
  # -----------------------------
  # Save if save_path supplied
  # -----------------------------
  if(!is.null(save_path)) {
    save(root_frame, file = save_path)
  }
  
  return(root_frame)
}
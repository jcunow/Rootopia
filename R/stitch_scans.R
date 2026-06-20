## ---------------------------------------------------------------------------
## Sequential image stitching for (mini)rhizotron / flatbed scan sequences
##
## Ported from the Python "ImageStitcher" (jcunow/ImageStitcher). The pipeline
## groups frames that belong to the same tube, sorts them, and sequentially
## stitches each group into one long mosaic. Consecutive frames are aligned with
## an edge-based FFT phase correlation and composited with a linear feather
## blend across the overlap band.
##
## Geometry: img2 is placed at column (w1 - edge_width + dx) so that
##   overlap = edge_width - dx
## with dx the (negated) horizontal correlation peak and dy the vertical offset
## applied to img2; the canvas grows to absorb dy. These signs reproduce the
## validated reference stitch (the Python original's working configuration).
##
## Axis convention: frames are stitched along the image WIDTH (left -> right).
## For sequences acquired ALONG the tube set direction = "vertical" (frames are
## transposed internally, stitched the same way, and transposed back).
##
## The user-facing functions (list_tubes, list_scan_files, stitch_root_scans,
## stitch_image_sequence, stitch_image_pair) are exported and documented with
## examples; the alignment engine and helpers below them are @keywords internal.
## ---------------------------------------------------------------------------


# =============================================================================
# Public API
# =============================================================================

#' List the tubes (groups) found in a scan folder
#'
#' Summarises the unique groups (tubes) that \code{\link{stitch_root_scans}}
#' would build, so you can see the tube names and pick a range to stitch
#' (e.g. \code{tubes = 1:36}).
#'
#' @inheritParams list_scan_files
#' @return A data frame with columns \code{index} (1-based, the value to pass to
#'   \code{stitch_root_scans(tubes = ...)}), \code{tube} (group id) and
#'   \code{n_frames}.
#' @seealso \code{\link{stitch_root_scans}}, \code{\link{list_scan_files}}
#' @keywords internal
list_tubes <- function(input, pattern = NULL, group_regex = "T0\\d{2}") {
  files <- stitch_discover_files(input, pattern)
  g <- stitch_group_of(files, group_regex)
  g <- g[!is.na(g)]
  if (length(g) == 0)
    return(data.frame(index = integer(0), tube = character(0), n_frames = integer(0)))
  tubes <- sort(unique(g))
  counts <- table(g)
  data.frame(index = seq_along(tubes), tube = tubes,
             n_frames = as.integer(counts[tubes]), stringsAsFactors = FALSE)
}


#' List scan files with index and group id
#'
#' Discovers the image files that \code{\link{stitch_root_scans}} would process
#' and returns them as an indexed table, so you can see what is in a folder and
#' pick a range of files to stitch (e.g. \code{select = 1:36}). For tube-level
#' selection see \code{\link{list_tubes}}. Uses the same discovery and sort order
#' as \code{stitch_root_scans}, so the indices line up.
#'
#' @param input Either a directory (searched recursively) or a character vector
#'   of image file paths.
#' @param pattern Optional substring used to keep only matching file names
#'   (e.g. \code{".tiff"}). \code{NULL} keeps all files.
#' @param group_regex Regular expression identifying the group id within each
#'   path. Default \code{"T0\\d{2}"} matches tube labels such as \code{T067}.
#' @return A data frame with columns \code{index} (1-based position), \code{file}
#'   (full path) and \code{group} (matched id, or \code{NA}).
#' @seealso \code{\link{list_tubes}}, \code{\link{stitch_root_scans}}
#' @keywords internal
list_scan_files <- function(input, pattern = NULL, group_regex = "T0\\d{2}") {
  files <- stitch_discover_files(input, pattern)
  data.frame(
    index = seq_along(files),
    file  = files,
    group = stitch_group_of(files, group_regex),
    stringsAsFactors = FALSE
  )
}


#' Batch-stitch grouped scan sequences (tubes) into mosaics
#'
#' High-level driver ported from the Python \code{ImageStitcher} \code{main} /
#' \code{process_subset}. Files are discovered, optionally subset, grouped by an
#' id pattern (one group per tube), sorted within each group, and stitched into
#' one mosaic per group. Single-frame groups are passed through unchanged.
#' Mosaics are returned and, optionally, written to disk as PNGs.
#'
#' Call \code{\link{list_tubes}} first to see the tube names, then pass a range
#' to \code{tubes} (e.g. \code{tubes = 1:36}) to stitch just those tubes.
#'
#' @param input Either a directory (searched recursively) or a character vector
#'   of image file paths.
#' @param pattern Optional substring used to keep only matching file names
#'   (e.g. \code{".tiff"}). \code{NULL} keeps all files.
#' @param group_regex Regular expression identifying the group id within each
#'   path. Default \code{"T0\\d{2}"} matches tube labels such as \code{T067}.
#'   Use \code{NULL} to stitch every file into a single mosaic.
#' @param select Optional integer vector of indices into the (sorted) \emph{file}
#'   list, e.g. \code{1:36}. See \code{\link{list_scan_files}}. \code{NULL} uses
#'   all files. Applied before grouping.
#' @param tubes Optional \emph{tube} selection: integer indices into the sorted
#'   tube list (e.g. \code{1:36}, see \code{\link{list_tubes}}), a character
#'   vector of tube names (e.g. \code{c("T037", "T040")}), or the string
#'   \code{"ask"} to print the tubes and choose a range interactively in one
#'   call (interactive sessions only). \code{NULL} keeps all tubes.
#' @param out_dir Optional directory to write one mosaic per tube to (named
#'   \code{<out_prefix><tube>.<ext>}). Created if needed. \code{NULL} (default)
#'   returns mosaics only.
#' @param out_prefix Filename prefix for written mosaics.
#' @param out_format Output image format when \code{out_dir} is set: \code{"png"}
#'   (default) or \code{"tiff"}. Requires the corresponding package
#'   (\pkg{png} or \pkg{tiff}).
#' @inheritParams stitch_image_pair
#' @param report Logical. If \code{TRUE}, return a list with both the mosaics
#'   and a per-step performance table instead of just the mosaics (see Value).
#' @param verbose Logical; print per-tube progress and a mean/min alignment peak.
#' @return If \code{report = FALSE} (default), invisibly a named list of mosaics
#'   (one numeric \code{(H, W, C)} array per tube). If \code{report = TRUE}, a
#'   list \code{list(mosaics, report)} where \code{report} is a data frame with
#'   columns \code{tube}, \code{step}, \code{dx}, \code{dy}, \code{peak}
#'   (confidence; higher is better) and \code{overlap} (\code{= edge_width - dx}).
#' @seealso \code{\link{list_tubes}}, \code{\link{list_scan_files}},
#'   \code{\link{stitch_image_sequence}}
#' @export
#' @examples
#' \dontrun{
#' # 1) See the tubes (names + frame counts)
#' list_tubes("path/to/scans", pattern = ".tiff")
#'
#' # 2) Stitch the first 36 tubes, with a performance report and a preprocess
#' res <- stitch_root_scans("path/to/scans", pattern = ".tiff",
#'                          tubes = 1:36, preprocess = "grad", report = TRUE)
#' res$report
#' aggregate(peak ~ tube, res$report, mean)
#'
#' # 3) A named subset, written straight to PNG
#' stitch_root_scans("path/to/scans", pattern = ".tiff",
#'                   tubes = c("T037", "T040"), out_dir = "path/to/output")
#' }
stitch_root_scans <- function(input, pattern = NULL, group_regex = "T0\\d{2}",
                              select = NULL, tubes = NULL, out_dir = NULL,
                              out_prefix = "", out_format = "png", method = "phase",
                              edge_width = 250, vertical_region = 1000,
                              vertical_offset = 300, direction = "horizontal",
                              preprocess = "none", blend = "linear", blend_width = NULL,
                              report = FALSE, verbose = TRUE) {
  tryCatch({
    out_format <- match.arg(out_format, c("png", "tiff"))
    out_ext <- if (out_format == "tiff") "tif" else "png"
    files <- stitch_discover_files(input, pattern)
    if (length(files) == 0) stop("No input files found")

    if (!is.null(select)) {
      if (!is.numeric(select)) stop("'select' must be numeric indices, e.g. 1:36")
      select <- as.integer(select)
      if (any(is.na(select)) || any(select < 1L) || any(select > length(files)))
        stop("'select' out of range: there are ", length(files),
             " files. Use list_scan_files() to see valid indices.")
      files <- files[select]
    }

    groups <- stitch_group_of(files, group_regex)
    keep <- !is.na(groups)
    if (any(!keep) && verbose)
      message("Skipping ", sum(!keep), " file(s) without a group id matching ",
              "'", group_regex, "'")
    files <- files[keep]; groups <- groups[keep]
    if (length(files) == 0) stop("No files matched group_regex '", group_regex, "'")

    unique_groups <- sort(unique(groups))

    if (is.character(tubes) && length(tubes) == 1L && tubes == "ask") {
      if (!interactive())
        stop("tubes = 'ask' requires an interactive session; ",
             "pass tube indices (e.g. 1:36) or names instead.")
      tubes <- stitch_prompt_tubes(unique_groups, groups)
    }

    if (!is.null(tubes)) {
      if (is.numeric(tubes)) {
        ti <- as.integer(tubes)
        if (any(is.na(ti)) || any(ti < 1L) || any(ti > length(unique_groups)))
          stop("'tubes' index out of range: ", length(unique_groups),
               " tube(s) available. Use list_tubes() to see them.")
        unique_groups <- unique_groups[ti]
      } else if (is.character(tubes)) {
        miss <- setdiff(tubes, unique_groups)
        if (length(miss)) stop("Unknown tube(s): ", paste(miss, collapse = ", "))
        unique_groups <- unique_groups[unique_groups %in% tubes]
      } else stop("'tubes' must be numeric indices (e.g. 1:36) or tube names")
    }

    if (!is.null(out_dir) && !dir.exists(out_dir))
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    if (verbose) message("Stitching ", length(unique_groups), " tube(s): ",
                         paste(unique_groups, collapse = ", "))

    results <- vector("list", length(unique_groups))
    names(results) <- unique_groups
    report_rows <- list()

    for (gi in seq_along(unique_groups)) {
      g <- unique_groups[gi]
      members <- sort(files[groups == g])
      if (verbose) message("[", gi, "/", length(unique_groups), "] ", g,
                           " - ", length(members), " frame(s)")

      if (length(members) == 1) {
        results[[gi]] <- stitch_to_hwc(members[[1]])
      } else {
        sr <- stitch_image_sequence(members, method = method,
                                    edge_width = edge_width, vertical_region = vertical_region,
                                    vertical_offset = vertical_offset, direction = direction,
                                    preprocess = preprocess, blend = blend,
                                    blend_width = blend_width, return_offsets = TRUE)
        results[[gi]] <- sr$mosaic
        if (nrow(sr$offsets)) {
          report_rows[[length(report_rows) + 1L]] <-
            data.frame(tube = g, sr$offsets, stringsAsFactors = FALSE)
          if (verbose)
            message("      alignment peak: mean=", round(mean(sr$offsets$peak), 4),
                    " min=", round(min(sr$offsets$peak), 4))
        }
      }

      if (!is.null(out_dir))
        stitch_write_image(results[[gi]],
                           file.path(out_dir, paste0(out_prefix, g, ".", out_ext)))
    }

    if (report) {
      rep_df <- if (length(report_rows)) {
        do.call(rbind, report_rows)[, c("tube", "step", "dx", "dy", "peak", "overlap")]
      } else {
        data.frame(tube = character(0), step = integer(0), dx = numeric(0),
                   dy = numeric(0), peak = numeric(0), overlap = numeric(0))
      }
      return(invisible(list(mosaics = results, report = rep_df)))
    }
    invisible(results)
  }, error = function(e) stop("stitch_root_scans failed: ", e$message, call. = FALSE))
}


#' Sequentially stitch a sequence of scans into one mosaic
#'
#' Frames are stitched in the given order: the first pair is composited, then
#' each subsequent frame is stitched onto the growing mosaic, in memory (no
#' temporary PNG round-trips).
#'
#' @param images A character vector of image file paths, or a list of image
#'   objects accepted by \code{\link{load_flexible_image}}. Stitched in the
#'   order given, so sort beforehand if needed.
#' @inheritParams stitch_image_pair
#' @param return_offsets Logical. If \code{TRUE}, return a list with the mosaic
#'   and a per-step performance data frame (\code{step}, \code{dx}, \code{dy},
#'   \code{peak}, \code{overlap}) instead of just the mosaic.
#' @return A numeric \code{(height, width, channel)} array (0-255); or, when
#'   \code{return_offsets = TRUE}, a list \code{list(mosaic, offsets)}.
#' @seealso \code{\link{stitch_image_pair}}, \code{\link{stitch_root_scans}}
#' @keywords internal
#' @examples
#' \dontrun{
#' set.seed(1)
#' img <- array(runif(70 * 200 * 3) * 255, dim = c(70, 200, 3))
#' frames <- list(img[, 1:90, , drop = FALSE],
#'                img[, 71:160, , drop = FALSE],
#'                img[, 141:200, , drop = FALSE])
#' res <- stitch_image_sequence(frames, edge_width = 25, vertical_region = 70,
#'                              vertical_offset = 0, return_offsets = TRUE)
#' dim(res$mosaic)
#' res$offsets
#' }
stitch_image_sequence <- function(images, method = "phase",
                                  edge_width = 250, vertical_region = 1000,
                                  vertical_offset = 300, direction = "horizontal",
                                  preprocess = "none", blend = "linear", blend_width = NULL,
                                  return_offsets = FALSE) {
  tryCatch({
    if (is.null(images)) stop("'images' is required")
    if (is.character(images)) images <- as.list(images)
    n <- length(images)
    if (n == 0) stop("'images' is empty")
    empty_offsets <- data.frame(step = integer(0), dx = numeric(0), dy = numeric(0),
                                peak = numeric(0), overlap = numeric(0))
    if (n == 1) {
      acc <- stitch_to_hwc(images[[1]])
      if (return_offsets) return(list(mosaic = acc, offsets = empty_offsets))
      return(acc)
    }

    acc <- stitch_to_hwc(images[[1]])
    offs <- vector("list", n - 1)
    for (i in 2:n) {
      nxt <- stitch_to_hwc(images[[i]])
      res <- stitch_compose_pair(acc, nxt, method = method,
                                 edge_width = edge_width, vertical_region = vertical_region,
                                 vertical_offset = vertical_offset, direction = direction,
                                 preprocess = preprocess, blend = blend, blend_width = blend_width)
      acc <- res$mosaic
      offs[[i - 1]] <- data.frame(step = i - 1L, dx = res$dx, dy = res$dy,
                                  peak = res$peak, overlap = res$overlap)
    }
    if (return_offsets) return(list(mosaic = acc, offsets = do.call(rbind, offs)))
    acc
  }, error = function(e) stop("stitch_image_sequence failed: ", e$message, call. = FALSE))
}


#' Stitch two scans into one mosaic with linear feather blending
#'
#' Core compositing step ported from the Python \code{ImageStitcher}
#' \code{stitch_two_images}. \code{img2} is placed to the right of \code{img1}
#' at the offset estimated by an edge-based FFT phase correlation, so that
#' \code{overlap = edge_width - dx}; the overlap band is blended with a
#' horizontal alpha ramp (1 -> 0) and the vertical offset is applied by growing
#' the canvas.
#'
#' @param img1,img2 Image inputs (file path, \code{SpatRaster}, \code{Raster*},
#'   \code{cimg}, \code{magick-image}, matrix, or array). See
#'   \code{\link{load_flexible_image}} for supported formats.
#' @param method Alignment method. Currently only \code{"phase"} (FFT phase
#'   correlation) is supported. \code{"feature"} (SIFT/ORB keypoint matching in
#'   the Python original) requires OpenCV and has no CRAN-compatible R backend;
#'   it raises an informative error.
#' @param edge_width Width in pixels of the edge band used for alignment.
#'   Clamped to the image width. Set it close to the true overlap - roughly
#'   1-2x (so the overlap is between about half and all of \code{edge_width}) -
#'   for the most reliable peak.
#' @param vertical_region Height in pixels of the vertical band used for
#'   alignment.
#' @param vertical_offset Starting row (from the top) of the vertical band.
#' @param direction \code{"horizontal"} (default, frames side by side) or
#'   \code{"vertical"} (frames stacked top to bottom, e.g. minirhizotron strips
#'   acquired down the tube). Vertical transposes the inputs, runs the same
#'   horizontal core, and transposes the result back.
#' @param preprocess Preprocessing applied to the edge bands before correlation:
#'   one of \code{"none"} (default), \code{"center"} (subtract mean),
#'   \code{"norm"} (divide by SD), \code{"center_norm"}, \code{"hann"}
#'   (demean + Hann window), \code{"grad"} (gradient magnitude) or
#'   \code{"grad_norm"}. \code{"hann"}, \code{"grad"} and \code{"center"} help on
#'   scans with uneven lighting or sparse texture; \code{"norm"} alone barely
#'   changes phase correlation (it is already scale-invariant). On well-textured
#'   scans \code{"none"} is usually fine.
#' @param blend How the overlap band is combined: \code{"linear"} (default,
#'   alpha ramp 1 -> 0, good for colour scans), \code{"overlay"} (img2 hides
#'   img1), \code{"overlay_first"} (img1 hides img2), \code{"max"} (lighten /
#'   union - recommended for segmented/binary masks, where averaging would make
#'   fractional values and ghost thin roots) or \code{"min"} (darken).
#' @param blend_width Optional width (px) of the linear ramp, centred in the
#'   overlap (hard img1 to its left, hard img2 to its right); \code{NULL}
#'   (default) ramps across the whole overlap. A smaller value reduces root
#'   ghosting on colour scans. Ignored unless \code{blend = "linear"}.
#' @return A numeric \code{(height, width, channel)} array in 0-255. The
#'   estimated \code{c(dx, dy, peak)} is attached as \code{attr(., "offset")}
#'   (\code{dx} the horizontal shift, with \code{overlap = edge_width - dx}).
#'   Convert for plotting with e.g. \code{terra::rast()} /
#'   \code{terra::plotRGB()}.
#' @seealso \code{\link{stitch_image_sequence}}, \code{\link{stitch_root_scans}}
#' @keywords internal
#' @examples
#' \dontrun{
#' set.seed(1)
#' img   <- array(runif(80 * 160 * 3) * 255, dim = c(80, 160, 3))
#' left  <- img[, 1:100, , drop = FALSE]
#' right <- img[, 81:160, , drop = FALSE]          # 20 px overlap
#' mosaic <- stitch_image_pair(left, right, edge_width = 30,
#'                             vertical_region = 80, vertical_offset = 0)
#' dim(mosaic)
#' attr(mosaic, "offset")
#' }
stitch_image_pair <- function(img1, img2, method = "phase",
                              edge_width = 250, vertical_region = 1000,
                              vertical_offset = 300, direction = "horizontal",
                              preprocess = "none", blend = "linear", blend_width = NULL) {
  tryCatch({
    a1 <- stitch_to_hwc(img1)
    a2 <- stitch_to_hwc(img2)
    res <- stitch_compose_pair(a1, a2, method = method,
                               edge_width = edge_width, vertical_region = vertical_region,
                               vertical_offset = vertical_offset, direction = direction,
                               preprocess = preprocess, blend = blend, blend_width = blend_width)
    out <- res$mosaic
    attr(out, "offset") <- c(dx = res$dx, dy = res$dy, peak = res$peak)
    out
  }, error = function(e) stop("stitch_image_pair failed: ", e$message, call. = FALSE))
}


# =============================================================================
# Internal: alignment engine
# =============================================================================

#' Edge-based phase-correlation alignment of two scans
#'
#' Estimates the relative pixel shift between the right edge of \code{img1} and
#' the left edge of \code{img2} using FFT phase correlation - the alignment step
#' behind \code{\link{stitch_image_pair}}.
#'
#' Only a vertical band of the edges is used (rows
#' \code{vertical_offset} .. \code{vertical_offset + vertical_region}, clamped to
#' the image) and only the outermost \code{edge_width} columns, optionally
#' \code{preprocess}ed. The cross-power spectrum is
#' \eqn{F_1 \bar{F_2} / (|F_1||F_2| + \epsilon)} and the shift is read from the
#' correlation peak (with wrap-around past the half dimension). The returned
#' \code{dx} is the placement shift: \code{overlap = edge_width - dx}. \code{dy}
#' is the vertical shift to apply to \code{img2}.
#'
#' @param img1,img2 Image inputs accepted by \code{\link{load_flexible_image}}.
#' @param edge_width Width in pixels of the edge band used for alignment
#'   (clamped to the image width).
#' @param vertical_region Height in pixels of the vertical band used for
#'   alignment.
#' @param vertical_offset Starting row (from the top) of the vertical band.
#' @param preprocess Preprocessing of the edge bands: one of \code{"none"},
#'   \code{"center"}, \code{"norm"}, \code{"center_norm"}, \code{"hann"},
#'   \code{"grad"} or \code{"grad_norm"}.
#' @return Named numeric vector \code{c(dx, dy, peak)}: \code{dx} horizontal
#'   placement shift (\code{overlap = edge_width - dx}), \code{dy} vertical shift
#'   (pixels) and \code{peak} the normalised correlation peak height.
#' @seealso \code{\link{stitch_image_pair}}, \code{\link{estimate_rotation_shift}}
#' @keywords internal
align_phase_correlation <- function(img1, img2, edge_width = 250,
                                    vertical_region = 1000, vertical_offset = 300,
                                    preprocess = "none") {
  tryCatch({
    if (is.null(img1) || is.null(img2)) stop("Both input images are required")
    preprocess <- match.arg(preprocess,
                            c("none", "center", "norm", "center_norm",
                              "hann", "grad", "grad_norm"))
    a1 <- stitch_to_hwc(img1); a2 <- stitch_to_hwc(img2)
    g1 <- stitch_luma(a1);     g2 <- stitch_luma(a2)

    h1 <- nrow(g1); w1 <- ncol(g1)
    h2 <- nrow(g2); w2 <- ncol(g2)

    # vertical band (0-based, end-exclusive), clamped to each image
    ys1 <- min(vertical_offset, max(0, h1 - vertical_region))
    ys2 <- min(vertical_offset, max(0, h2 - vertical_region))
    ye1 <- ys1 + min(vertical_region, h1 - ys1)
    ye2 <- ys2 + min(vertical_region, h2 - ys2)
    ew  <- min(edge_width, w1, w2)
    if (ew < 2) stop("edge_width too small after clamping to image width")

    e1 <- g1[(ys1 + 1):ye1, (w1 - ew + 1):w1, drop = FALSE]   # right edge of img1
    e2 <- g2[(ys2 + 1):ye2, 1:ew, drop = FALSE]               # left  edge of img2

    # phase correlation needs identical dims; crop to the common band height
    nr <- min(nrow(e1), nrow(e2))
    e1 <- stitch_preprocess(e1[1:nr, , drop = FALSE], preprocess)
    e2 <- stitch_preprocess(e2[1:nr, , drop = FALSE], preprocess)
    nc <- ncol(e1)

    f1  <- stats::fft(e1)
    f2  <- stats::fft(e2)
    cps <- (f1 * Conj(f2)) / (Mod(f1) * Mod(f2) + 1e-10)
    corr <- Mod(stats::fft(cps, inverse = TRUE))

    lin    <- which.max(corr) - 1L          # column-major linear index, 0-based
    dy_raw <- lin %% nr                      # row    (Y) peak position
    dx_raw <- lin %/% nr                     # column (X) peak position
    if (dy_raw > nr %/% 2) dy_raw <- dy_raw - nr
    if (dx_raw > nc %/% 2) dx_raw <- dx_raw - nc

    # Placement-ready shifts: img2 at (w1 - edge_width + dx, dy). The horizontal
    # peak is negated, the vertical is not - this reproduces the validated stitch.
    c(dx = -dx_raw, dy = dy_raw, peak = max(corr) / sum(corr))
  }, error = function(e) stop("align_phase_correlation failed: ", e$message, call. = FALSE))
}


#' Compose two already-loaded (H, W, C) arrays into a mosaic
#'
#' The compositing core shared by \code{\link{stitch_image_pair}} and
#' \code{\link{stitch_image_sequence}}.
#'
#' @param a1,a2 Already-loaded \code{(H, W, C)} numeric arrays (0-255).
#' @param method,edge_width,vertical_region,vertical_offset,direction,preprocess,blend,blend_width
#'   As described for \code{\link{stitch_image_pair}}.
#' @return A list \code{list(mosaic, dx, dy, peak, overlap)}.
#' @keywords internal
stitch_compose_pair <- function(a1, a2, method, edge_width,
                                vertical_region, vertical_offset,
                                direction = "horizontal", preprocess = "none",
                                blend = "linear", blend_width = NULL) {
  direction <- match.arg(direction, c("horizontal", "vertical"))
  blend <- match.arg(blend, c("linear", "overlay", "overlay_first", "max", "min"))
  if (direction == "vertical") {           # transpose, stitch horizontally, transpose back
    a1 <- aperm(a1, c(2, 1, 3)); a2 <- aperm(a2, c(2, 1, 3))
  }

  ch <- stitch_match_channels(a1, a2); a1 <- ch[[1]]; a2 <- ch[[2]]
  ct <- dim(a1)[3]

  if (identical(method, "feature")) {
    stop("method = 'feature' (SIFT/ORB keypoint matching) needs OpenCV and has ",
         "no CRAN-compatible R backend. Use method = 'phase'. If you require ",
         "feature matching, run the original Python ImageStitcher, or bind ",
         "OpenCV via the (non-CRAN) 'Rvision' package.", call. = FALSE)
  }
  if (!identical(method, "phase")) stop("method must be 'phase' or 'feature'")

  al <- align_phase_correlation(a1, a2, edge_width = edge_width,
                                vertical_region = vertical_region,
                                vertical_offset = vertical_offset, preprocess = preprocess)
  dx <- as.numeric(al["dx"]); dy <- as.numeric(al["dy"]); peak <- as.numeric(al["peak"])

  h1 <- dim(a1)[1]; w1 <- dim(a1)[2]
  h2 <- dim(a2)[1]; w2 <- dim(a2)[2]
  ew <- min(edge_width, w1, w2)

  # placement: img2 at column (w1 - ew + dx) so that overlap = ew - dx
  img2_x <- as.integer(round(w1 - ew + dx))
  if (img2_x < 0)
    stop("Estimated horizontal placement is negative (dx = ", round(dx),
         "); check edge_width / frame order.")
  actual_overlap <- max(0L, w1 - img2_x)
  img2_y <- as.integer(round(dy))

  if (img2_y >= 0) {
    img1_ys <- 0L; img2_ys <- img2_y
    out_h <- max(h1, h2 + img2_y)
  } else {
    img1_ys <- -img2_y; img2_ys <- 0L
    out_h <- max(h1 - img2_y, h2)
  }
  out_w <- img2_x + w2

  if (actual_overlap <= 0)
    warning("Frames do not overlap (dx = ", round(dx),
            "); mosaic may contain gaps.", call. = FALSE)

  stitched <- array(0, dim = c(out_h, out_w, ct))
  stitched[(img1_ys + 1):(img1_ys + h1), 1:w1, ] <- a1   # place img1

  # combine the overlap band according to `blend`
  ov_start <- img2_x; ov_end <- w1; ov_w <- ov_end - ov_start
  by_start <- max(img1_ys, img2_ys)
  by_end   <- min(img1_ys + h1, img2_ys + h2)
  if (by_end > by_start && ov_w > 0) {
    bh <- by_end - by_start
    r1 <- by_start - img1_ys
    r2 <- by_start - img2_ys
    reg1 <- a1[(r1 + 1):(r1 + bh), (ov_start + 1):ov_end, , drop = FALSE]
    reg2 <- a2[(r2 + 1):(r2 + bh), 1:ov_w, , drop = FALSE]
    mh <- min(dim(reg1)[1], dim(reg2)[1])
    mw <- min(dim(reg1)[2], dim(reg2)[2])
    if (mh > 0 && mw > 0) {
      reg1 <- reg1[1:mh, 1:mw, , drop = FALSE]
      reg2 <- reg2[1:mh, 1:mw, , drop = FALSE]
      stitched[(by_start + 1):(by_start + mh), (ov_start + 1):(ov_start + mw), ] <-
        stitch_blend_overlap(reg1, reg2, blend, blend_width)
    }
  }

  # copy the non-overlapping remainder of img2
  if (w1 < out_w) {
    src_x <- actual_overlap
    cw <- w2 - actual_overlap
    if (cw > 0) {
      stitched[(img2_ys + 1):(img2_ys + h2), (w1 + 1):(w1 + cw), ] <-
        a2[, (src_x + 1):(src_x + cw), , drop = FALSE]
    }
  }

  if (direction == "vertical") stitched <- aperm(stitched, c(2, 1, 3))
  list(mosaic = stitched, dx = dx, dy = dy, peak = peak, overlap = actual_overlap)
}


# =============================================================================
# Internal: array / IO helpers
# =============================================================================

#' Preprocess an edge band before correlation
#'
#' @param m A numeric matrix (edge band).
#' @param method One of \code{"none"}, \code{"center"} (subtract mean),
#'   \code{"norm"} (divide by SD), \code{"center_norm"}, \code{"hann"}
#'   (demean + separable Hann window), \code{"grad"} (central-difference
#'   gradient magnitude) or \code{"grad_norm"}.
#' @return The preprocessed matrix.
#' @keywords internal
stitch_preprocess <- function(m, method = "none") {
  if (identical(method, "none")) return(m)
  nr <- nrow(m); nc <- ncol(m)
  grad_mag <- function(x) {
    gx <- matrix(0, nr, nc); gy <- matrix(0, nr, nc)
    if (nc >= 3) gx[, 2:(nc - 1)] <- (x[, 3:nc] - x[, 1:(nc - 2)]) / 2
    if (nr >= 3) gy[2:(nr - 1), ] <- (x[3:nr, ] - x[1:(nr - 2), ]) / 2
    sqrt(gx^2 + gy^2)
  }
  if (method %in% c("grad", "grad_norm")) m <- grad_mag(m)
  if (method %in% c("center", "center_norm")) m <- m - mean(m)
  if (method == "hann") {
    hn <- function(n) if (n <= 1) rep(1, max(n, 1L)) else 0.5 - 0.5 * cos(2 * pi * (0:(n - 1)) / (n - 1))
    m <- (m - mean(m)) * outer(hn(nr), hn(nc))
  }
  if (method %in% c("norm", "center_norm", "grad_norm"))
    m <- m / (stats::sd(as.vector(m)) + 1e-10)
  m
}

#' Combine the overlap band of two aligned regions
#'
#' @param reg1,reg2 Numeric \code{(h, w, c)} arrays - the aligned overlap of
#'   img1 (left) and img2 (right).
#' @param blend One of \code{"linear"} (alpha ramp 1 -> 0), \code{"overlay"}
#'   (img2 on top), \code{"overlay_first"} (img1 on top), \code{"max"}
#'   (lighten / union) or \code{"min"} (darken).
#' @param blend_width Optional ramp width (px) for \code{"linear"}, centred in
#'   the overlap; \code{NULL} ramps across the whole overlap.
#' @return The combined \code{(h, w, c)} array.
#' @keywords internal
stitch_blend_overlap <- function(reg1, reg2, blend = "linear", blend_width = NULL) {
  d <- dim(reg1); mh <- d[1]; mw <- d[2]; ct <- d[3]
  if (blend == "overlay")       return(reg2)        # img2 hides img1
  if (blend == "overlay_first") return(reg1)        # img1 hides img2
  if (blend == "max") { o <- pmax(reg1, reg2); dim(o) <- d; return(o) }
  if (blend == "min") { o <- pmin(reg1, reg2); dim(o) <- d; return(o) }
  # linear feather: alpha 1 (left, img1) -> 0 (right, img2)
  alpha <- seq.int(1, 0, length.out = mw)
  if (!is.null(blend_width)) {                       # ramp only over a centred band
    bw <- max(1L, min(as.integer(blend_width), mw))
    left <- (mw - bw) %/% 2L
    alpha <- c(rep(1, left), seq.int(1, 0, length.out = bw), rep(0, mw - left - bw))
  }
  aw <- array(rep(matrix(alpha, mh, mw, byrow = TRUE), ct), dim = d)
  reg1 * aw + reg2 * (1 - aw)
}

#' Interactively prompt for a tube selection
#'
#' @param unique_groups Sorted character vector of tube ids.
#' @param groups Character vector of every kept file's tube id (for counts).
#' @return Integer indices into \code{unique_groups}, or \code{NULL} for all.
#' @keywords internal
stitch_prompt_tubes <- function(unique_groups, groups) {
  counts <- as.integer(table(groups)[unique_groups])
  tab <- data.frame(index = seq_along(unique_groups), tube = unique_groups,
                    n_frames = counts, stringsAsFactors = FALSE)
  message("Tubes found:")
  print(tab, row.names = FALSE)
  ans <- trimws(readline("Select tubes (e.g. 1:36 or 1,4,9; blank = all): "))
  if (!nzchar(ans)) return(NULL)
  tryCatch(as.integer(eval(parse(text = paste0("c(", ans, ")")), envir = baseenv())),
           error = function(e) stop("Could not parse tube selection: ", ans, call. = FALSE))
}

#' Discover and sort the image files to stitch
#'
#' @param input A directory path or a character vector of file paths.
#' @param pattern Optional substring filter on file names.
#' @return A sorted character vector of file paths.
#' @keywords internal
stitch_discover_files <- function(input, pattern = NULL) {
  if (is.character(input) && length(input) == 1 && dir.exists(input)) {
    files <- list.files(input, recursive = TRUE, full.names = TRUE)
  } else if (is.character(input)) {
    files <- input
  } else {
    stop("'input' must be a directory path or a character vector of file paths")
  }
  if (!is.null(pattern)) files <- files[grepl(pattern, files, fixed = TRUE)]
  sort(files)
}

#' Load any supported image to a (height, width, channel) array in 0-255
#'
#' Normalises orientation so that dimension 1 is rows (Y), dimension 2 is
#' columns (X) and dimension 3 is channels - matching the (H, W, C) layout used
#' by the Python reference (OpenCV) so the stitching math ports one-to-one.
#'
#' @param input An image in any format accepted by
#'   \code{\link{load_flexible_image}}.
#' @return A numeric \code{(H, W, C)} array scaled to 0-255.
#' @keywords internal
stitch_to_hwc <- function(input) {
  a <- load_flexible_image(input, output_format = "array", scale = "none")
  d <- dim(a)
  if (length(d) == 4) {
    # imager path returns (width, height, depth = 1, channels)
    a <- a[, , 1, , drop = TRUE]
    if (length(dim(a)) == 2) a <- array(a, c(dim(a), 1))
    a <- aperm(a, c(2, 1, 3))            # (W, H, C) -> (H, W, C)
  } else if (length(d) == 2) {
    a <- array(a, c(d, 1))               # grayscale (H, W) -> (H, W, 1)
  } else if (length(d) != 3) {
    stop("Unsupported image dimensions for stitching: ", paste(d, collapse = "x"))
  }
  # everything else (tiff, SpatRaster, plain arrays) is assumed (H, W, C)
  mx <- suppressWarnings(max(a, na.rm = TRUE))
  if (is.finite(mx) && mx <= 1) a <- a * 255    # bring 0-1 sources to 0-255
  storage.mode(a) <- "double"
  a
}

#' Rec. 601 luma projection of an (H, W, C) array
#'
#' OpenCV's COLOR_BGR2GRAY uses the BT.601 weights (0.299, 0.587, 0.114); the
#' reference is ported with the same weights so the correlation peak matches.
#'
#' @param a A numeric \code{(H, W, C)} array.
#' @return A numeric \code{(H, W)} luma matrix.
#' @keywords internal
stitch_luma <- function(a) {
  cc <- dim(a)[3]
  if (is.na(cc) || cc < 3) return(a[, , 1])
  0.299 * a[, , 1] + 0.587 * a[, , 2] + 0.114 * a[, , 3]
}

#' Force two (H, W, C) arrays to share the same channel count
#'
#' @param a1,a2 Numeric \code{(H, W, C)} arrays.
#' @return A list of the two arrays expanded to a common channel count.
#' @keywords internal
stitch_match_channels <- function(a1, a2) {
  c1 <- dim(a1)[3]; c2 <- dim(a2)[3]
  ct <- max(c1, c2)
  grow <- function(a, ct) {
    c <- dim(a)[3]
    if (c == ct) return(a)
    if (c == 1) return(array(rep(as.vector(a), ct), c(dim(a)[1], dim(a)[2], ct)))
    stop("Incompatible channel counts (", c1, " vs ", c2, ")")
  }
  list(grow(a1, ct), grow(a2, ct))
}

#' Extract a group id (e.g. tube label) from each path
#'
#' @param x Character vector of paths.
#' @param group_regex Regular expression for the group id, or \code{NULL} to
#'   place everything in a single group.
#' @return Character vector of group ids (\code{NA} where no match).
#' @keywords internal
stitch_group_of <- function(x, group_regex) {
  if (is.null(group_regex)) return(rep("mosaic", length(x)))
  m <- regexpr(group_regex, x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  out[m != -1L] <- regmatches(x, m)
  out
}

#' Write an (H, W, C) 0-255 array to PNG or TIFF (chosen by file extension)
#'
#' @param a A numeric \code{(H, W, C)} array (0-255).
#' @param path Output path; the extension (\code{.png}, \code{.tif}/\code{.tiff})
#'   selects the format.
#' @return The path, invisibly.
#' @keywords internal
stitch_write_image <- function(a, path) {
  img <- a / 255
  img[img < 0] <- 0; img[img > 1] <- 1
  ext <- tolower(tools::file_ext(path))
  if (ext == "png") {
    if (!requireNamespace("png", quietly = TRUE))
      stop("Package 'png' is required to write PNG output. ",
           "Install it with install.packages(\"png\").")
    png::writePNG(img, target = path)
  } else if (ext %in% c("tif", "tiff")) {
    if (!requireNamespace("tiff", quietly = TRUE))
      stop("Package 'tiff' is required to write TIFF output. ",
           "Install it with install.packages(\"tiff\").")
    tiff::writeTIFF(img, where = path)
  } else {
    stop("Unsupported output extension '", ext, "'; use 'png' or 'tiff'.")
  }
  invisible(path)
}

#' Compute root traits over a depth profile from root scans
#'
#' @description
#' Processes a directory of root images and returns a tidy data frame of root
#' traits summarized per depth interval.  Handles flatbed scans and flat
#' rhizotron windows by default, and cylindrical minirhizotron tubes when tube
#' parameters are supplied (see \strong{Geometry}).  Input may be already
#' segmented or a raw grayscale/RGB scan, which is binarized on the way in (see
#' \strong{Binarization}).
#'
#' \code{batch_root_traits()} is an alias for the same function.
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
#'   \item{\code{path_seg}}{Always required.}
#'   \item{\code{path_skl}}{Optional. Used by \code{calc_root_length},
#'     \code{calc_diameter_stats}, \code{calc_diameter_quantiles},
#'     \code{calc_root_angles}, and \code{calc_root_order_metrics} when
#'     supplied. If \code{path_skl} is \code{NULL} or a skeleton file is
#'     missing for an image, the skeleton is computed internally from the
#'     segmented image via \code{skeletonize_image()}.}
#'   \item{\code{path_rgb}}{Required when \code{calc_color_metrics} is
#'     \code{TRUE}.}
#' }
#'
#' @param path_seg Character. Path to directory of binary segmented images
#'   (foreground/root pixel = 1, background = 0).
#' @param path_skl Character or \code{NULL}. Path to directory of skeletonised
#'   images (one-pixel-wide centerlines of roots), used for length, diameter,
#'   angle, and branching-order metrics. If \code{NULL} or a file is missing,
#'   the skeleton is computed internally via \code{skeletonize_image()}.
#'   Default \code{NULL}.
#' @param path_rgb Character or \code{NULL}. Path to directory of blended
#'   RGB images, aligned to the segmented images.  Required for color metrics.
#'   Default \code{NULL}.
#' @param seg_file_index Integer vector or \code{NULL}. Optional subset index
#'   applied to \code{list.files(path_seg)}, e.g. \code{37:72}.  Default
#'   \code{NULL} (use all files).
#' @param skl_file_index Integer vector or \code{NULL}. Optional subset index
#'   applied to \code{list.files(path_skl)}.  Default \code{NULL}.
#' @param rgb_file_index Integer vector or \code{NULL}. Optional subset index
#'   applied to \code{list.files(path_rgb)}.  Default \code{NULL}.
#'
#' @section Per-image metadata:
#' All metadata arguments below accept either a single value (recycled to all
#' images) or a vector of length equal to the number of images.
#'
#' @param insertion_angles Numeric or \code{NULL}. Insertion angle of the tube,
#'   in \strong{degrees measured from horizontal} -- \code{90} is a vertical
#'   tube, \code{45} a tube pushed in at 45 degrees.  This is the convention
#'   \code{create_depthmap()} uses for \code{tilt}, and it is the opposite of
#'   the "degrees from vertical" wording used by some minirhizotron software,
#'   so convert before passing it in.  Values must be strictly between 0 and 90.
#'   Supplying this argument switches the run to minirhizotron geometry (see
#'   \strong{Geometry}).  Default \code{NULL} (flatbed).
#' @param soil_starts Numeric. Pixel row (in the original, un-rotated image)
#'   at which the soil surface begins.  Used to set the zero-depth reference.
#'   Default \code{0}.
#' @param tube_names Character or \code{NULL}. Sample/tube identifiers added as
#'   the \code{Tube} column.  If \code{NULL}, names are derived from characters
#'   3-5 from the right of the segmented file name, prefixed with \code{"T"}
#'   (e.g. \code{"T042"}).  Adjust if your naming convention differs.
#'
#'   Names must be \strong{unique}, one per image, and the run stops if they are
#'   not.  They key the tube-level joins, where two images sharing a name are
#'   read as one tube and their rows are multiplied out.  The derived default
#'   collides whenever the files share a suffix -- a directory of
#'   \code{"CLAS1A10-15fine.tif"} names every image \code{"Tfin"} -- so flatbed
#'   scans generally need this argument.
#'   Default \code{NULL}.
#' @param session Character. Session or campaign label added as the
#'   \code{Session} column, e.g. \code{"2022_02"}.  Default \code{""}.
#'
# --- Scan and geometry settings ---
#' @inheritParams create_depthmap
#' @param tube_diameter_cm Numeric or \code{NULL}. Inner diameter of the
#'   minirhizotron tube in \strong{centimetres}.  Passed to
#'   \code{create_depthmap()} as \code{tube_thicc}, where it sets the amplitude
#'   and wavelength of the sinusoidal curvature correction.  Supplying this
#'   argument switches the run to minirhizotron geometry (see
#'   \strong{Geometry}).  Default \code{NULL} (flatbed).
#' @param tube_center_offset Numeric in \code{[0, 1]}. Phase of the sinusoidal
#'   curvature correction: where the top of the tube falls across the image
#'   height, as a fraction.  \code{0} (default) puts it at the first row.  Set it
#'   if your scanner's rotational reference differs, or the depth assigned to a
#'   root will be offset by up to half a tube diameter.  Minirhizotron only.
#' @param depth_interval_cm Numeric or \code{NULL}. Size of each depth bin in
#'   \strong{centimetres}.  Passed as \code{nn} to \code{binning()}.
#'   \code{NULL} switches on whole-image mode, where the scan is treated as a
#'   single bin and summarized in one row (see \strong{Whole-image mode}).
#'   Default \code{5}.
#' @param bin_round Character. How \code{binning()} assigns a depth to a bin:
#'   \code{"rounding"} (default, \code{nn * round(depth/nn)}), \code{"floor"},
#'   or \code{"ceiling"}.  \strong{Note what "rounding" does to the top bin}: with
#'   \code{depth_interval_cm = 5} it spans 0-2.5 cm while every other bin spans
#'   5 cm, because the label is the bin's center rather than its top edge.  Per-bin
#'   densities are unaffected (they divide by each bin's own measured area), but
#'   \code{mrd} and \code{total.length.density} multiply by
#'   \code{depth_interval_cm} as though every bin were full width, so they are
#'   biased by the half-width top bin.  \code{"floor"} gives the soil-science
#'   convention -- 0-5, 5-10, labeled by the shallower edge -- and is the better
#'   choice for a new analysis; the default is kept for continuity with existing
#'   ones.
#' @param rotation_fixed_width Numeric. Width in \strong{rows} that each image is
#'   cropped to along the rotation axis, centered on the middle row, before any
#'   trait is measured (see \code{rotation_censor()}).  This trims the tube
#'   edges, where the curvature of the tube distorts what the scanner sees.  An
#'   image with fewer rows than this cannot be cropped symmetrically, so
#'   \code{rotation_censor()} clamps to the image bounds and says so -- the
#'   image is then used at full width.
#'   Only applied under minirhizotron geometry; a flatbed scan is measured at
#'   full width.  Default \code{1800}.
#'
#' @section Geometry:
#' There is no geometry switch.  The geometry follows from whether you supply
#' tube parameters:
#' \describe{
#'   \item{Flatbed (default)}{Neither \code{insertion_angles} nor
#'     \code{tube_diameter_cm} is supplied.  The scan is treated as a flat
#'     surface imaged head-on: no sinusoidal curvature correction
#'     (\code{sinoid = FALSE}), no foreshortening of the depth axis (one pixel
#'     along the image width is one pixel of depth), and no
#'     \code{rotation_censor()} crop, since there is no tube interior to crop
#'     to.}
#'   \item{Minirhizotron}{Either argument is supplied.  The sinusoidal
#'     curvature correction is switched on, the depth axis is foreshortened by
#'     \code{sin(insertion_angles)}, and each image is cropped to
#'     \code{rotation_fixed_width} rows about its center.  An argument you
#'     leave out falls back to a neutral default: \code{tube_diameter_cm = 7}
#'     and \code{insertion_angles = 90} (vertical tube), both announced in the
#'     run log.}
#' }
#' In both cases depth runs along the image \strong{width} (left to right),
#' which is the orientation minirhizotron scanners produce.  A flatbed scan with
#' the soil surface at the top must be rotated 90 degrees before it is passed
#' in, or the depth profile will be built across the wrong axis.
#'
#' @section Whole-image mode:
#' A flatbed scan of washed roots in a tray has no depth axis -- where a root
#' lies on the tray says nothing about where it grew.  Setting
#' \code{depth_interval_cm = NULL} treats each image as one bin, so the result
#' carries one row per image, with \code{depth = 0}, and the density metrics
#' become whole-scan quantities: \code{rootpx.density} is the percentage of the
#' scan covered by root and \code{rootlength.density} is cm of root per
#' cm\eqn{^2} of scanned area -- how densely the tray was packed.  Both already
#' divide by the bin's own measured pixel area (\code{rootpx + voidpx}), so
#' there is nothing to bin by and no depth map is built, unless
#' \code{calc_root_angles = TRUE} asks for one.  Nothing is masked or trimmed
#' per slice either.  On a large scan both are worth skipping: a depth map is a
#' full-size double raster.  (\code{calc_root_length} still builds a flat
#' surface of its own to read flow directions from; that one belongs to the
#' length estimator, not to the binning, and is built either way.)
#'
#' \code{calc_distribution_indices} and \code{calc_advanced_metrics} are
#' switched off in this mode: mean rooting depth, and each bin's share of the
#' profile, mean nothing when there is only one bin.
#'
#' @section Which parameters to set:
#' Most runs set five or six of these arguments, and which five depends on the
#' geometry: a third of them describe a minirhizotron tube and do nothing at all
#' to a flatbed scan, with no warning when you tune one that is inert.  The
#' \emph{Batch Processing} vignette has the full tier list for both geometries;
#' in short:
#' \describe{
#'   \item{Flatbed}{Set \code{dpi}, \code{tube_names} and
#'     \code{depth_interval_cm = NULL}, plus \code{binarize_threshold} and
#'     \code{dark_roots} for unsegmented scans.  \code{insertion_angles},
#'     \code{tube_diameter_cm}, \code{tube_center_offset},
#'     \code{rotation_fixed_width}, \code{soil_starts} and \code{bin_round}
#'     do nothing here.}
#'   \item{Minirhizotron}{Set \code{dpi}, \code{insertion_angles} (degrees from
#'     \emph{horizontal}), \code{tube_diameter_cm}, \code{soil_starts},
#'     \code{tube_names}, \code{depth_interval_cm} and \code{bin_round}, and
#'     look at \code{rotation_fixed_width} -- how much of the tube's curved edge
#'     to trim is a property of your scanner that no default can guess.}
#' }
#' The per-image arguments -- \code{insertion_angles}, \code{soil_starts},
#' \code{binarize_threshold}, \code{dark_roots}, \code{tube_names} -- each take
#' either one value for the whole run or one per image, in the order of
#' \code{list.files(path_seg)}.
#'
#' @section Binarization:
#' Flatbed scans are usually delivered as grayscale or RGB with the full 0-255
#' range, not as a segmented mask.  Such an image is reduced to a single
#' grayscale layer (\code{rgb2gray()} for RGB input) and cut at
#' \code{binarize_threshold}, in the same way RhizoVision Explorer does it.
#' Already-segmented input passes through untouched: with the default
#' \code{binarize = "auto"} the cut is only applied when the image actually
#' has more than two gray levels.
#'
#' @param binarize Either \code{"auto"} (default), \code{TRUE}, or
#'   \code{FALSE}.  \code{"auto"} thresholds an image only if it has more
#'   than two distinct gray levels, so binary masks are left alone and raw scans
#'   are binarized.  \code{TRUE} always thresholds, \code{FALSE} never does
#'   (any non-zero pixel is then taken as root).  Note that \code{TRUE} on an
#'   image already coded 0/1 will \emph{invert} it, because 0 counts as dark;
#'   this is what \code{"auto"} exists to prevent.
#' @param binarize_threshold Numeric. Gray level at which a scan is cut into
#'   root and background, on the \strong{0-255} scale (the RhizoVision Explorer
#'   convention).  Images that load on a 0-1 scale get the same cut-off
#'   rescaled, so the number means the same thing either way.  Default
#'   \code{200}.
#' @param dark_roots Logical.  \code{TRUE} (default) means roots are
#'   \emph{darker} than the background, as on a flatbed scan of washed roots on
#'   a white tray: pixels at or below \code{binarize_threshold} become root.
#'   \code{FALSE} inverts this for bright-roots-on-dark-background images.  If
#'   more than half the pixels come out as root, the polarity is probably wrong
#'   and a warning says so.
#' @param clean_max_hole_size Numeric. Fill enclosed background holes of up to
#'   this many pixels before anything is measured; \code{0} (default) fills
#'   none, \code{Inf} fills every hole.  Pinholes inside a painted root read as
#'   background, which eats into the distance transform and so into diameter.
#' @param clean_max_artifact_size Numeric. Drop disconnected foreground blobs of
#'   up to this many pixels; \code{0} (default) drops none, \code{Inf} drops
#'   everything not touching the largest structure.  Specks count as root area,
#'   as isolated skeleton pixels, and as root tips.  Both use
#'   \code{\link{clean_image}} and need the \pkg{imager} package.
#' @param seg_layer Integer or \code{NULL}. Which layer of the segmented image
#'   to measure.  \code{NULL} (default) picks automatically: a single-layer
#'   image is used as is, and a 3- or 4-layer image is converted to grayscale
#'   with \code{rgb2gray()}.  Set this if your files carry the segmentation in
#'   one specific band.
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
#'   \code{root_diameter()}, plus \code{root.surface.area} (lateral surface
#'   area, cm\eqn{^2}) and \code{root.volume} (cm\eqn{^3}) per bin, and
#'   \code{rootsurface_rootvolume_ratio} (cm\eqn{^{-1}}).  \code{root.surface.area}
#'   is the curved wall wrapping each cylindrical root (its soil-contact area),
#'   \emph{not} the flat root area visible in the image.
#'   \code{rootsurface_rootvolume_ratio} is the length-weighted mean of the
#'   \emph{local} ratio \eqn{2 / r_i} over skeleton pixels, so it is dominated by
#'   fine roots and is deliberately not equal to
#'   \code{root.surface.area / root.volume} (a bulk ratio dominated by thick
#'   roots); the two answer different questions.  Default \code{TRUE}.
#'
# --- Extended metrics (off by default) ---
#' @param calc_diameter_quantiles Logical. Compute the diameter distribution
#'   percentiles set by \code{diameter_quantiles} per bin, conditional means
#'   above each quantile, and threshold-based root lengths (see
#'   \code{diameter_thresholds}).  Default \code{FALSE}.
#' @param calc_modal_peaks Logical. Compute modal diameter peaks per bin via
#'   \code{modal_peaks()} (\code{n.diameter.peaks}, \code{diameter.peak.1/2/3}).
#'   Auto-enables \code{calc_diameter_quantiles}, since it reuses the same
#'   per-bin diameter raster.  \strong{Slow}: one call per depth bin per
#'   image.  Default \code{FALSE}.
#' @param calc_landscape_metrics Logical. Compute patch-level landscape metrics
#'   per depth bin via \code{root_scape_metrics()}: nearest-neighbor distance
#'   (\code{enn_mn}), joint entropy (\code{joinent}), relative mutual
#'   information (\code{relmutinf}), number of patches (\code{np}), and
#'   contagion (\code{contag}).  \strong{Slow}: one call per depth bin per
#'   image.  Default \code{FALSE}.
#' @param calc_color_metrics Logical. Compute mean chromatic coordinates (rcc,
#'   gcc, bcc), hue, saturation, luminosity, and raw RGB channel means
#'   separately for root pixels and background pixels via
#'   \code{tube_coloration()}.  Requires \code{path_rgb}.  Default \code{FALSE}.
#' @param calc_root_angles Logical. Compute \code{deep_drive} (fraction of
#'   skeleton pixels whose D8 flow direction matches the locally optimal
#'   downward direction) and \code{mean.steepness.angle} /
#'   \code{sd.steepness.angle} (degrees, 0 = horizontal, 90 = vertical).
#'   Uses \code{deep_drive()}.  Default \code{FALSE}.
#' @param calc_root_order_metrics Logical. Build a per-image branching-order
#'   graph via \code{branch_order_map()} and summarize it both per depth bin
#'   and per tube.  Adds \code{mean.<scheme>} and \code{max.<scheme>} for the
#'   scheme named by \code{order_scheme} (e.g. \code{mean.strahler_order}),
#'   plus \code{mean.root_order} and \code{lateral_root_fraction} per depth
#'   bin, tube-level \code{main_root.*} / \code{lateral_roots.*} columns
#'   (length, diameter, branching frequency, etc.), and \code{n_root_orders}
#'   (the highest order found).  Requires a skeleton.  \strong{Slow}:
#'   builds one segment graph per image.  Default \code{FALSE}.
#'
# --- Branching-order settings (used when calc_root_order_metrics = TRUE) ---
#' @param order_scheme Character. Which ordering labels the per-bin order
#'   columns and the main-root split.  \code{"strahler_order"} (default) is the
#'   convention of the fine-root literature (Pregitzer et al. 2002; Fitter):
#'   every distal unbranched root is order 1, and the order rises only where two
#'   roots of equal order meet.  \code{"tip_order"} is the same leaf-peeling
#'   with the order raised at \emph{every} junction, so an axis carrying ten
#'   laterals reaches order 11.  \code{"branch_order"} counts the other way --
#'   the thickest root of each component is 1 and its laterals 2 -- and
#'   \code{"root_order"} gives each continuous root the maximum
#'   \code{tip_order} along it.  All four are computed regardless; this only
#'   picks which one is reported.  See \code{\link{branch_order_map}}.
#' @param diam_weight Numeric >= 0. At a junction, which two arms are read as
#'   the same root continuing: \code{straightness + diam_weight *
#'   diameter_similarity}.  \code{0} uses the angle alone, larger values let
#'   thickness decide.  Affects \code{root_order}, \code{branch_order} and the
#'   segment grouping, not \code{strahler_order} or \code{tip_order}.
#'   Default \code{0.5}.
#' @param prune_spur_length_cm Numeric. Remove terminal skeleton branches
#'   ("spurs") shorter than this, in \strong{centimetres}, before any trait is
#'   measured; \code{0} (default) prunes nothing.  Thinning leaves short stubs
#'   where roots are wide or ragged, and each one is counted as a root tip and
#'   adds a little length, so the pruning applies to the skeleton every metric
#'   is measured from -- length and diameter as well as the order graph.  Start
#'   around two to three times the width of your thickest root and check one
#'   image before trusting a batch.  Uses \code{\link{prune_skeleton}}.
#' @param prune_spur_iter Integer. Pruning passes, so that a spur exposed by
#'   removing another is caught too.  More passes eat further into real roots.
#'   Default \code{1}.
#'
# --- Derived metrics ---
#' @param calc_density_metrics Logical. Compute \code{rootpx.density} (percent
#'   root area cover) and \code{rootlength.density} (cm root length per cm^2
#'   imaged area) per bin.  Auto-enables \code{calc_root_pixels} and
#'   \code{calc_root_length}.  Default \code{TRUE}.
#' @param calc_distribution_indices Logical. Compute tube-level indices:
#'   \code{mrd} (mean rooting depth) and \code{total.length.density} (summed
#'   length density over all bins, in cm root per cm^2 per cm depth).
#'   Auto-enables \code{calc_density_metrics}.  Default \code{TRUE}.
#' @param calc_advanced_metrics Logical. Compute per-bin derived metrics:
#'   \code{rootlength.fraction} (each bin's length density as a fraction of the
#'   tube total) and \code{mean.var.diameter} (mean of within-bin diameter
#'   variance).  Auto-enables \code{calc_distribution_indices} and
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
#' \strong{Surface, volume, and their ratio.}  For each skeleton pixel, the root
#' segment is modeled as a cylinder of length \eqn{l_i} (one pixel edge in cm)
#' and radius \eqn{r_i} (half the local diameter in cm).  Its \emph{lateral}
#' surface area is \eqn{2 \pi r_i l_i} -- the curved wall wrapping the root, i.e.
#' the area in contact with the soil, not the flat root area seen in the image --
#' and its volume is \eqn{\pi r_i^2 l_i}.  \code{root.surface.area} (cm\eqn{^2})
#' and \code{root.volume} (cm\eqn{^3}) are these quantities summed over all
#' skeleton pixels in the depth bin.
#'
#' The per-pixel surface-to-volume ratio simplifies to \eqn{2 / r_i}, and
#' \code{rootsurface_rootvolume_ratio} (cm\eqn{^{-1}}) is the length-weighted
#' mean of \eqn{2 / r_i} over the bin's skeleton pixels.  Because it averages the
#' \emph{local} ratio, it is dominated by fine roots (small \eqn{r_i} give large
#' \eqn{2 / r_i}).  This is deliberately not the same as the bulk ratio
#' \code{root.surface.area / root.volume}, which is dominated by thick roots
#' (they hold most of the volume); the two summarize different things and will
#' not match unless every root in the bin has the same diameter.
#'
#' \strong{Fault tolerance.}  Every metric block is wrapped in
#' \code{tryCatch}.  Failures produce a \code{[Rootopia] SKIPPED} message and
#' \code{NA} values; they never abort the run.
#'
#' \strong{Dependency resolution.}  Enabling a higher-level metric silently
#' enables its prerequisites and prints a message listing what was auto-enabled.
#'
#' @examples
#' \dontrun{
#' # Flatbed scans, not yet binarized -- the default path.
#' # No tube arguments, so this is flat geometry: no curvature correction,
#' # no foreshortening, no tube crop. Roots are dark on a bright tray, so
#' # everything at or below gray level 200 is taken as root.
#' result <- batch_root_traits(
#'   path_seg           = "scans/flatbed/tray_01/",
#'   dpi                = 600,
#'   binarize_threshold = 200,
#'   session            = "2024_spring"
#' )
#'
#' # Minirhizotron -- supplying an insertion angle switches the geometry
#' result <- root_depth_metrics(
#'   path_seg         = "scans/segmented/2022_02/",
#'   path_skl         = "scans/skeleton/2022_02/",
#'   insertion_angles = tube_meta$angle,   # degrees from horizontal
#'   tube_diameter_cm = 7,
#'   session          = "2022_02"
#' )
#'
#' # With diameter quantiles and root angle metrics
#' result <- root_depth_metrics(
#'   path_seg                = "scans/segmented/2022_02/",
#'   path_skl                = "scans/skeleton/2022_02/",
#'   path_rgb                = "scans/blended/2022_02/",
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
#' # Flat rhizotron window, already segmented. binarize = "auto" sees a
#' # two-level image and leaves it alone.
#' result <- root_depth_metrics(
#'   path_seg = "scans/segmented/rhizotron_A/",
#'   path_skl = "scans/skeleton/rhizotron_A/"
#' )
#' }
#'
#' @importFrom dplyr group_by filter summarize full_join rename_with
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_sub
#' @importFrom terra rast zonal values ext focal terrain subst flip t trim
#'   resample mask compareGeom levels
#'
#' @export
root_depth_metrics <- function(
    
  # ---------- paths ----------------------------------------------------------
  path_seg,
  path_skl                = NULL,
  path_rgb                = NULL,
  seg_file_index          = NULL,
  skl_file_index          = NULL,
  rgb_file_index          = NULL,
  
  # ---------- per-image metadata ---------------------------------------------
  insertion_angles        = NULL,
  soil_starts             = 0,
  tube_names              = NULL,
  session                 = "",
  
  # ---------- scan / geometry ------------------------------------------------
  dpi                     = 300,
  tube_diameter_cm        = NULL,
  tube_center_offset      = 0,
  depth_interval_cm       = 5,
  bin_round               = c("rounding", "floor", "ceiling"),
  rotation_fixed_width    = 1800,
  
  # ---------- binarization (raw flatbed scans) -------------------------------
  binarize                = "auto",
  binarize_threshold      = 200,
  dark_roots              = TRUE,
  seg_layer               = NULL,
  clean_max_hole_size     = 0,
  clean_max_artifact_size = 0,
  
  # ---------- core metrics (on by default) -----------------------------------
  calc_root_pixels        = TRUE,
  calc_root_length        = TRUE,
  calc_diameter_stats     = TRUE,
  
  # ---------- extended metrics (off by default) ------------------------------
  calc_diameter_quantiles = FALSE,
  calc_modal_peaks        = FALSE,
  calc_landscape_metrics  = FALSE,
  calc_color_metrics      = FALSE,
  calc_root_angles        = FALSE,
  calc_root_order_metrics = FALSE,

  # ---------- branching-order settings ---------------------------------------
  order_scheme            = c("strahler_order", "branch_order",
                              "root_order", "tip_order"),
  diam_weight             = 0.5,
  prune_spur_length_cm    = 0,
  prune_spur_iter         = 1L,
  
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
               message(sprintf("[Rootopia] SKIPPED '%s': %s", label, conditionMessage(e)))
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
  
  # Reduce whatever came off disk to a single-layer 0/1 root mask.
  #
  # A flatbed scan arrives as grayscale or RGB spanning the full 0-255 range, so
  # it has to be cut at a gray level before anything downstream can treat a
  # pixel as root. An already-segmented image arrives as 0/1 or 0/255 and must
  # pass through untouched -- that is what binarize = "auto" is for: it only
  # cuts when the image really has more than two gray levels.
  #
  # Everything except the color metrics works on this single layer, which is
  # why the reduction happens here at load time rather than per metric.
  .as_root_mask <- function(img, label, thr255, dark) {
    
    nl <- terra::nlyr(img)
    
    gray <- if (!is.null(seg_layer)) {
      if (seg_layer > nl)
        stop(sprintf("seg_layer = %d but '%s' has only %d layer(s).",
                     seg_layer, label, nl), call. = FALSE)
      img[[seg_layer]]
    } else if (nl == 1L) {
      img
    } else if (nl >= 3L) {
      rgb2gray(img[[1:3]])            # drops any alpha band along the way
    } else {
      img[[1]]                        # 2 layers: nothing sensible to weight
    }
    
    vals   <- terra::values(gray)
    n_lev  <- length(unique(vals[!is.na(vals)]))
    do_bin <- if (identical(binarize, "auto")) n_lev > 2L else isTRUE(binarize)
    
    if (!do_bin) return((gray > 0) * 1)
    
    # binarize_threshold is a gray level on the 0-255 scale. Images that loaded
    # on a 0-1 scale get the same cut-off rescaled, so the number the user typed
    # means the same thing either way.
    mx  <- max(vals, na.rm = TRUE)
    thr <- if (mx <= 1) thr255 / 255 else thr255
    
    mask <- if (dark) (gray <= thr) * 1 else (gray >= thr) * 1
    
    # A mask that is mostly foreground nearly always means the polarity is
    # inverted, and it is far cheaper to say so here than to wonder about the
    # numbers three hours later.
    frac <- mean(terra::values(mask), na.rm = TRUE)
    if (!is.na(frac) && frac > 0.5)
      warning(sprintf(
        paste0("'%s': %.0f%% of pixels classified as root at binarize_threshold = %s. ",
               "If roots are brighter than the background here, set dark_roots = FALSE."),
        label, frac * 100, format(thr255)), call. = FALSE)
    
    mask
  }
  
  # ===========================================================================
  # 0b. Geometry: flatbed unless tube parameters were supplied
  # ===========================================================================
  # There is no toggle. Asking for a tube diameter or an insertion angle is
  # what makes this a minirhizotron run; asking for neither makes it a flatbed
  # run. The two paths differ in three places: the sinusoidal curvature
  # correction, the foreshortening of the depth axis, and the rotation_censor()
  # crop to the tube interior.
  tube_geometry <- !is.null(insertion_angles) || !is.null(tube_diameter_cm)
  
  if (tube_geometry) {
    if (is.null(tube_diameter_cm)) {
      tube_diameter_cm <- 7
      .msg("[Rootopia] Minirhizotron geometry: no tube_diameter_cm given, using %g cm.",
           tube_diameter_cm)
    }
    if (is.null(insertion_angles)) {
      insertion_angles <- 90
      .msg("[Rootopia] Minirhizotron geometry: no insertion_angles given, assuming a vertical tube (90 degrees from horizontal).")
    }
    if (!is.numeric(insertion_angles) || any(is.na(insertion_angles)) ||
        any(insertion_angles <= 0 | insertion_angles > 90))
      stop(paste0("'insertion_angles' must be numeric and in (0, 90], measured from ",
                  "horizontal (90 = vertical tube). Note this is the opposite of the ",
                  "'degrees from vertical' convention some software uses."),
           call. = FALSE)
    .msg("[Rootopia] Minirhizotron geometry: %g cm tube, insertion angle(s) %s degrees from horizontal.",
         tube_diameter_cm, paste(unique(insertion_angles), collapse = ", "))
  } else {
    # Flat surface imaged head-on: one pixel along the width is one pixel of
    # depth, so the depth axis must not be foreshortened. create_depthmap()
    # scales it by sin(tilt), which makes 90 the neutral value. tube_thicc is
    # unused when sinoid = FALSE but must still be a positive number.
    insertion_angles <- 90
    tube_diameter_cm <- 7
    .msg(paste0("[Rootopia] Flatbed geometry: neither insertion_angles nor tube_diameter_cm given, ",
                "so no tube curvature correction and no tube-interior crop. ",
                "Supply either one for minirhizotron scans."))
  }
  
  # ===========================================================================
  # 0c. Depth binning: a profile, unless depth_interval_cm is NULL
  # ===========================================================================
  # NULL says this image has no depth axis, as for a tray of washed roots. Every
  # pixel then falls in the same bin and the densities become whole-scan
  # numbers, because they already divide by the bin's own pixel area.
  order_scheme <- match.arg(order_scheme)
  bin_round    <- match.arg(bin_round)
  if (!is.numeric(diam_weight) || length(diam_weight) != 1L ||
      is.na(diam_weight) || diam_weight < 0)
    stop("'diam_weight' must be a single number >= 0.", call. = FALSE)
  if (!is.numeric(prune_spur_length_cm) || length(prune_spur_length_cm) != 1L ||
      is.na(prune_spur_length_cm) || prune_spur_length_cm < 0)
    stop("'prune_spur_length_cm' must be a single number >= 0 (0 = no pruning).",
         call. = FALSE)

  whole_image <- is.null(depth_interval_cm)
  if (whole_image) {
    .msg(paste0("[Rootopia] Whole-image mode (depth_interval_cm = NULL): one row per image, ",
                "densities per cm^2 of scanned area."))
  } else if (!is.numeric(depth_interval_cm) || length(depth_interval_cm) != 1L ||
             is.na(depth_interval_cm) || depth_interval_cm <= 0) {
    stop("'depth_interval_cm' must be a single positive number, or NULL for whole-image mode.",
         call. = FALSE)
  }

  # ===========================================================================
  # 1.  Resolve file lists
  # ===========================================================================
  if (!dir.exists(path_seg))
    stop("'path_seg' does not exist: ", path_seg, call. = FALSE)
  
  all_seg <- list.files(path_seg)
  im.ls.seg <- if (!is.null(seg_file_index)) all_seg[seg_file_index] else all_seg
  if (length(im.ls.seg) == 0)
    stop("No files found in 'path_seg': ", path_seg, call. = FALSE)
  n_images <- length(im.ls.seg)
  
  im.ls.skl <- NULL
  if (!is.null(path_skl) && dir.exists(path_skl)) {
    all_skl <- list.files(path_skl)
    im.ls.skl <- if (!is.null(skl_file_index)) all_skl[skl_file_index] else all_skl
  }
  
  im.ls.rgb <- NULL
  if (!is.null(path_rgb) && dir.exists(path_rgb)) {
    all_rgb <- list.files(path_rgb)
    im.ls.rgb <- if (!is.null(rgb_file_index)) all_rgb[rgb_file_index] else all_rgb
  }
  
  # Recycle per-image metadata
  insertion_angles   <- .recycle(insertion_angles,   n_images, "insertion_angles")
  soil_starts        <- .recycle(soil_starts,        n_images, "soil_starts")
  # Exposure drifts across a scanning session, so these are per image too.
  binarize_threshold <- .recycle(binarize_threshold, n_images, "binarize_threshold")
  dark_roots         <- .recycle(dark_roots,         n_images, "dark_roots")
  if (!is.null(tube_names))
    tube_names <- .recycle(tube_names, n_images, "tube_names")
  else
    tube_names <- paste0("T", stringr::str_sub(im.ls.seg, start = -8, end = -6))

  # Tube names label the rows and key the tube-level joins at the end, where two
  # images sharing a name are one tube to dplyr: the join multiplies their rows
  # out, and the run reports a row count nobody asked for. The default names are
  # three characters cut from the file name, so a set of files ending in the same
  # suffix collapses onto one name. Cheaper to say so here than to explain the
  # row count later.
  dup <- unique(tube_names[duplicated(tube_names)])
  if (length(dup) > 0L) {
    shown <- paste(dup[seq_len(min(5L, length(dup)))], collapse = ", ")
    if (length(dup) > 5L) shown <- paste0(shown, ", ...")
    stop(sprintf(paste0(
      "Tube names must be unique, one per image, but %d name(s) are used more than once: %s.\n",
      "Without 'tube_names' a name is cut from characters 3-5 from the right of each file ",
      "name, so files sharing a suffix (\"fine.tif\") all land on the same name. Pass ",
      "'tube_names' explicitly, one name per file, in the order of list.files(path_seg)."),
      length(dup), shown), call. = FALSE)
  }

  # ===========================================================================
  # 2.  Global dependency resolution
  # ===========================================================================
  if (calc_root_angles && !calc_root_length) {
    message("[Rootopia] Auto-enabling calc_root_length (required for calc_root_angles).")
    calc_root_length <- TRUE
  }
  if (calc_diameter_quantiles && !calc_root_length) {
    message("[Rootopia] Auto-enabling calc_root_length (required for calc_diameter_quantiles).")
    calc_root_length <- TRUE
  }
  if (calc_modal_peaks && !calc_diameter_quantiles) {
    message("[Rootopia] Auto-enabling calc_diameter_quantiles (required for calc_modal_peaks).")
    calc_diameter_quantiles <- TRUE
  }
  if (calc_root_order_metrics && !calc_root_length) {
    message("[Rootopia] Auto-enabling calc_root_length (required for calc_root_order_metrics).")
    calc_root_length <- TRUE
  }
  if (calc_density_metrics) {
    if (!calc_root_pixels) {
      message("[Rootopia] Auto-enabling calc_root_pixels (required for calc_density_metrics).")
      calc_root_pixels <- TRUE
    }
    if (!calc_root_length) {
      message("[Rootopia] Auto-enabling calc_root_length (required for calc_density_metrics).")
      calc_root_length <- TRUE
    }
  }
  if (calc_distribution_indices && !calc_density_metrics) {
    message("[Rootopia] Auto-enabling calc_density_metrics (required for calc_distribution_indices).")
    calc_density_metrics <- TRUE
  }
  if (calc_advanced_metrics) {
    if (!calc_density_metrics) {
      message("[Rootopia] Auto-enabling calc_density_metrics (required for calc_advanced_metrics).")
      calc_density_metrics <- TRUE
    }
    if (!calc_distribution_indices) {
      message("[Rootopia] Auto-enabling calc_distribution_indices (required for calc_advanced_metrics).")
      calc_distribution_indices <- TRUE
    }
    if (!calc_diameter_stats) {
      message("[Rootopia] Auto-enabling calc_diameter_stats (required for mean.var.diameter).")
      calc_diameter_stats <- TRUE
    }
  }
  # Both of these describe how roots are spread over a profile, so a run with a
  # single bin has nothing for them to describe.
  if (whole_image && (calc_distribution_indices || calc_advanced_metrics)) {
    message(paste("[Rootopia] Whole-image mode: disabling calc_distribution_indices and",
                  "calc_advanced_metrics (they need more than one depth bin)."))
    calc_distribution_indices <- FALSE
    calc_advanced_metrics     <- FALSE
  }
  
  # Warn early about missing paths
  needs_skl <- calc_root_length || calc_diameter_stats || calc_diameter_quantiles ||
    calc_root_angles || calc_root_order_metrics
  if (needs_skl && is.null(im.ls.skl)) {
    message(paste(
      "[Rootopia] Skeleton directory (path_skl) not found or not supplied.",
      "Skeletons will be computed internally per image via skeletonize_image()."
    ))
  }
  if (calc_color_metrics && is.null(im.ls.rgb)) {
    warning("RGB directory (path_rgb) not found or not supplied. Disabling calc_color_metrics.",
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
  # output stays in step when the probabilities are customised. The default
  # c(0.90, 0.95, 0.99) yields the rootdiameter.{90,95,99} and
  # avg.diameter.top{10,5,1}pct columns.
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
    tube     <- tube_names[l]
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
    # scale = "none" on purpose: "binary" is ceiling(x / max), which turns every
    # non-zero gray value into root and would silently destroy an unbinarized
    # flatbed scan. .as_root_mask() does the cutting, with a real threshold.
    im <- .safe(sprintf("load segmented [%s]", seg_file), {
      img <- load_flexible_image(file.path(path_seg, seg_file),
                                 output_format = "spatrast",
                                 scale = "none")
      .as_root_mask(img, seg_file, binarize_threshold[l], dark_roots[l])
    })
    if (is.null(im)) {
      message(sprintf("[Rootopia] [%d/%d] %s: could not load segmented image -- skipping.",
                      l, n_images, tube))
      failed_imgs <- c(failed_imgs, seg_file)
      img_times   <- c(img_times, proc.time()[["elapsed"]] - t_img)
      next
    }
    
    # -------------------------------------------------------------------------
    # 3a2. Segmentation clean-up (optional)
    # -------------------------------------------------------------------------
    # Specks and pinholes in a painted segmentation become root pixels, root
    # tips and a little length. Cleaning before the skeleton is traced keeps
    # every metric measuring the same image.
    if (clean_max_hole_size > 0 || clean_max_artifact_size > 0) {
      im <- .safe("clean_image", {
        clean_image(im,
                    max_hole_size     = clean_max_hole_size,
                    max_artifact_size = clean_max_artifact_size,
                    output_format     = "spatrast")
      }, fallback = im)
    }

    im.skeleton <- NULL
    if (do_length || do_diam_st || do_diam_q || do_angles || do_order) {
      if (!is.null(im.ls.skl) && l <= length(im.ls.skl)) {
        im.skeleton <- .safe(sprintf("load skeleton [%s]", im.ls.skl[l]), {
          sk <- load_flexible_image(file.path(path_skl, im.ls.skl[l]),
                                    output_format = "spatrast",
                                    scale = "binary")
          # A skeleton file is binary by nature, so "binary" scaling is safe
          # here. Layer 2 only exists for multi-band files; single-band
          # skeletons are the norm for flatbed output.
          if (terra::nlyr(sk) > 1L) sk <- sk[[2]]
          sk
        })
      }
      if (is.null(im.skeleton)) {
        im.skeleton <- .safe(sprintf("skeletonize [%s]", seg_file), {
          skeletonize_image(im, verbose = FALSE)
        })
      }
      if (is.null(im.skeleton)) {
        message(sprintf("[Rootopia] %s: skeleton unavailable -- disabling length/diameter/angle/order metrics.", tube))
        do_length <- do_diam_st <- do_diam_q <- do_angles <- do_order <- FALSE
      }
    }
    
    im.rgb <- NULL
    if (do_color) {
      if (!is.null(im.ls.rgb) && l <= length(im.ls.rgb)) {
        im.rgb <- .safe(sprintf("load RGB [%s]", im.ls.rgb[l]),
                        terra::rast(file.path(path_rgb, im.ls.rgb[l])))
      }
      if (is.null(im.rgb)) {
        message(sprintf("[Rootopia] %s: RGB image unavailable -- disabling calc_color_metrics.", tube))
        do_color <- FALSE
      }
    }
    
    # -------------------------------------------------------------------------
    # 3b. Rotation censor (crop to tube interior) -- minirhizotron only
    # -------------------------------------------------------------------------
    # A flatbed scan has no tube interior, so cropping it to a fixed number of
    # rows about the center would just throw away real data.
    if (tube_geometry) {
      r0 <- round(dim(im)[1] / 2, 0)
      
      im <- .safe("rotation_censor (seg)",
                  rotation_censor(im, center_offset = r0, fixed_rotation = TRUE,
                                  fixed_width = rotation_fixed_width),
                  fallback = im)
      
      if (!is.null(im.skeleton))
        im.skeleton <- .safe("rotation_censor (skl)",
                             rotation_censor(im.skeleton, center_offset = r0, fixed_rotation = TRUE,
                                             fixed_width = rotation_fixed_width),
                             fallback = im.skeleton)
      
      if (!is.null(im.rgb))
        im.rgb <- .safe("rotation_censor (rgb)",
                        rotation_censor(im.rgb, center_offset = r0, fixed_rotation = TRUE,
                                        fixed_width = rotation_fixed_width),
                        fallback = im.rgb)
      
      # rotation_censor() hands back every layer it was given; the traits below
      # need the single mask layer again.
      if (terra::nlyr(im) > 1L) im <- im[[min(2L, terra::nlyr(im))]]
      if (!is.null(im.skeleton) && terra::nlyr(im.skeleton) > 1L)
        im.skeleton <- im.skeleton[[min(2L, terra::nlyr(im.skeleton))]]
    }
    
    # -------------------------------------------------------------------------
    # 3b2. Spur pruning (optional)
    # -------------------------------------------------------------------------
    # Thinning leaves short terminal spurs where roots are wide or ragged, and
    # they are counted as root tips and as a little extra length. Pruning here,
    # before anything measures the skeleton, keeps length, diameter and order
    # metrics telling the same story. Pruning after the crop, so that ends the
    # crop itself created are treated like any other terminal.
    if (prune_spur_length_cm > 0 && !is.null(im.skeleton)) {
      im.skeleton <- .safe("prune spurs", {
        prune_skeleton(im.skeleton, mask = im,
                       min_length = prune_spur_length_cm * dpi / 2.54,
                       iter = prune_spur_iter, output = "skeleton",
                       verbose = FALSE)
      }, fallback = im.skeleton)
    }

    # Align extents
    if (!is.null(im.skeleton)) terra::ext(im.skeleton) <- terra::ext(im)
    if (!is.null(im.rgb))      terra::ext(im.rgb)      <- terra::ext(im)
    
    # -------------------------------------------------------------------------
    # 3c. Depth map and binning
    # -------------------------------------------------------------------------
    # Whole-image mode has nothing to bin by, so this map is only worth building
    # when deep_drive() will read it; skipping it saves a full-size double
    # raster per image. (The length block below builds a flat surface of its own
    # for flow directions -- that one is part of the estimator, not the binning.)
    need_depthmap <- !whole_image || do_angles
    
    DepthMap <- if (need_depthmap) .safe("create_depthmap", {
      dm <- create_depthmap(
        img         = im,
        sinoid      = tube_geometry,
        dpi         = dpi,
        start_soil  = soil0,
        center_offset = tube_center_offset,
        tilt        = angle,
        tube_thicc  = tube_diameter_cm
      )
      # create_depthmap() returns a map on `im`'s grid, so only the extent has
      # to be carried over -- no transpose.
      terra::ext(dm) <- terra::ext(im)
      dm
    })
    if (need_depthmap && is.null(DepthMap)) {
      if (whole_image) {
        # The bins do not depend on it here, so only the angles are lost.
        message(sprintf("[Rootopia] %s: create_depthmap failed -- disabling calc_root_angles.", tube))
        do_angles <- FALSE
      } else {
        message(sprintf("[Rootopia] [%d/%d] %s: create_depthmap failed -- skipping.", l, n_images, tube))
        failed_imgs <- c(failed_imgs, seg_file)
        img_times   <- c(img_times, proc.time()[["elapsed"]] - t_img)
        next
      }
    }
    
    bm <- if (whole_image) {
      # One bin, numbered 0, holding every pixel of the scan.
      b <- terra::rast(im)
      terra::values(b) <- 0
      b
    } else {
      binning(depthmap = DepthMap, nn = depth_interval_cm, round_option = bin_round)
    }
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
        #
        # This builds a per-pixel length map rather than calling root_length()
        # per depth bin, and that is deliberate: root_length() uses the Kimura2
        # estimator, which is not additive -- the sum over bins does not equal
        # the whole. A per-pixel Freeman length can be zonal-summed, which is
        # what a depth profile needs, at the cost of Kimura's correction.
        dm_flat <- create_depthmap(
          img           = im,
          sinoid        = FALSE,
          dpi           = dpi,
          start_soil    = soil0,
          center_offset = 0.5,
          tilt          = angle,
          tube_thicc    = tube_diameter_cm
        )
        # create_depthmap() is already aligned with `im` (no transpose).
        dem_flat <- dm_flat
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
        # dpi is not optional here: root_diameter() converts pixels to cm with
        # 2.54 / dpi, and its own default is 300, so leaving it out silently
        # scales every diameter as though the scan were 300 dpi.
        dm <- root_diameter(im, skeleton_img = im.skeleton,
                            unit = "cm", dpi = dpi)$diameter_rast
        terra::ext(dm) <- terra::ext(bm)
        dm
      })
      if (is.null(rd.map)) {
        message(sprintf("[Rootopia] %s: diameter map failed -- disabling diameter metrics.", tube))
        do_diam_st <- do_diam_q <- FALSE
      }
    }
    
    # Diameter summary statistics (mean, max, variance per bin)
    if (do_diam_st && !is.null(rd.map)) {
      diam <- .safe("diameter stats", {
        avg <- terra::zonal(rd.map, bm, "mean", na.rm = TRUE); colnames(avg) <- c("depth", "avg.diameter")
        mx  <- terra::zonal(rd.map, bm, "max",  na.rm = TRUE); colnames(mx)  <- c("depth", "max.diameter")
        vr  <- terra::zonal(rd.map, bm, "var",  na.rm = TRUE); colnames(vr)  <- c("depth", "var.diameter")

        # Lateral root surface area and volume per bin. Each skeleton pixel is a
        # cylinder of diameter d (cm, from rd.map) and length one pixel edge
        # (px_cm = 2.54 / dpi). The *lateral* surface is the curved wall that
        # wraps around the cylindrical root (pi * d * length) -- the soil-contact
        # area of the whole root, not the flat root area visible in the image.
        # Volume is pi * (d/2)^2 * length. Summing per pixel and scaling by the
        # pixel length gives cm^2 (surface) and cm^3 (volume) per depth bin.
        px_cm <- 2.54 / dpi
        sa <- terra::zonal(pi * rd.map, bm, "sum", na.rm = TRUE)
        sa[, 2] <- sa[, 2] * px_cm;        colnames(sa) <- c("depth", "root.surface.area")
        vo <- terra::zonal(pi * (rd.map / 2)^2, bm, "sum", na.rm = TRUE)
        vo[, 2] <- vo[, 2] * px_cm;        colnames(vo) <- c("depth", "root.volume")

        # Surface-to-volume ratio per bin (cm^-1): the length-weighted mean of
        # the *local* per-pixel ratio 2 / r = 4 / diameter. Because skeleton
        # pixels share equal unit length, the plain mean of 4 / diameter over the
        # bin is that length-weighted mean. It is computed on the diameter raster
        # before aggregation, so it is unbiased (4 / mean(diameter) would not be,
        # since mean(1 / d) != 1 / mean(d)).
        # NOTE: this is the mean of local ratios, so it is dominated by fine
        # roots and is deliberately NOT equal to root.surface.area /
        # root.volume (the bulk ratio, which is dominated by thick roots). The
        # two answer different questions; do not expect them to match.
        sv  <- terra::zonal(4 / rd.map, bm, "mean", na.rm = TRUE)
        colnames(sv) <- c("depth", "rootsurface_rootvolume_ratio")
        Reduce(function(a, b) merge(a, b, by = "depth"), list(avg, mx, vr, sa, vo, sv))
      })
      if (!is.null(diam)) {
        roots <- merge(roots, diam, by = "depth", all.x = TRUE)
      } else {
        roots[c("avg.diameter", "max.diameter", "var.diameter",
                "root.surface.area", "root.volume",
                "rootsurface_rootvolume_ratio")] <- NA_real_
      }
    }
    
    # -------------------------------------------------------------------------
    # 3f2. Root branching-order metrics (per bin + per tube)
    # -------------------------------------------------------------------------
    if (do_order && !is.null(im.skeleton)) {
      ord_res <- .safe("root order metrics", {

        bo <- branch_order_map(skel = im.skeleton, mask = im, order = order_scheme,
                               unit = "cm", dpi = dpi, return_map = TRUE,
                               template = im.skeleton, diam_weight = diam_weight,
                               verbose = FALSE)
        et <- bo$edges

        bo_map <- bo$class_map
        terra::ext(bo_map) <- terra::ext(bm)
        ro_map <- order_classification_map(et, im.skeleton, value = "root_order")
        terra::ext(ro_map) <- terra::ext(bm)

        # Named for the scheme in use, so a column can never be read under the
        # wrong definition.
        bo_mean <- terra::zonal(bo_map, bm, "mean", na.rm = TRUE); colnames(bo_mean) <- c("depth", paste0("mean.", order_scheme))
        bo_max  <- terra::zonal(bo_map, bm, "max",  na.rm = TRUE); colnames(bo_max)  <- c("depth", paste0("max.",  order_scheme))
        ro_mean <- terra::zonal(ro_map, bm, "mean", na.rm = TRUE); colnames(ro_mean) <- c("depth", "mean.root_order")

        # Tube-level main-root vs lateral-root summary. The main root is the
        # order class with the largest length-weighted diameter, which under
        # branch_order is order 1 and under Strahler the highest order -- so the
        # split is taken from the diameters rather than from the order numbers,
        # and means the same thing whichever scheme is in use.
        om    <- order_metrics(bo, focal = "thickest")
        focal <- attr(om, "focal_orders")

        main_px <- bo_map                      # one vectorised pass, not per cell
        terra::values(main_px) <- as.numeric(terra::values(bo_map) %in% focal)
        lateral_px <- (!is.na(bo_map)) * 1 - main_px
        ordered_px <- (!is.na(bo_map)) * 1
        lat_z <- terra::zonal(lateral_px, bm, "sum", na.rm = TRUE); colnames(lat_z) <- c("depth", "lateral_px")
        tot_z <- terra::zonal(ordered_px, bm, "sum", na.rm = TRUE); colnames(tot_z) <- c("depth", "ordered_px")

        per_bin <- Reduce(function(a, b) merge(a, b, by = "depth", all = TRUE),
                          list(bo_mean, bo_max, ro_mean, lat_z, tot_z))
        per_bin$lateral_root_fraction <- ifelse(per_bin$ordered_px > 0,
                                                 per_bin$lateral_px / per_bin$ordered_px, NA_real_)
        per_bin$lateral_px <- per_bin$ordered_px <- NULL

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
        tube_row$n_root_orders <- suppressWarnings(max(et[[order_scheme]], na.rm = TRUE))
        if (!is.finite(tube_row$n_root_orders)) tube_row$n_root_orders <- NA_real_

        list(per_bin = per_bin, tube = tube_row)
      })

      if (!is.null(ord_res)) {
        roots <- merge(roots, ord_res$per_bin, by = "depth", all.x = TRUE)
        for (nm in names(ord_res$tube)) roots[[nm]] <- ord_res$tube[[nm]]
      } else {
        roots[c(paste0("mean.", order_scheme), paste0("max.", order_scheme),
                "mean.root_order", "lateral_root_fraction")] <- NA_real_
      }
    }

    roots$Tube <- tube

    # -------------------------------------------------------------------------
    # 3g. Density metrics (Level 1 derived)
    # -------------------------------------------------------------------------
    if (do_density) {
      .safe("density metrics", {
        if (do_pixels && all(c("rootpx", "voidpx") %in% names(roots))) {
          roots$rootpx.density <- roots$rootpx / (roots$rootpx + roots$voidpx) * 100
        }
        if (do_length && do_pixels && all(c("rootlength", "rootpx", "voidpx") %in% names(roots))) {
          # rootlength.density: cm root length per cm^2 of imaged area per bin
          roots$rootlength.density <-
            roots$rootlength / ((roots$rootpx + roots$voidpx) / (dpi / 2.54)^2)
        }
      })
    }
    
    # -------------------------------------------------------------------------
    # 3h. Per-depth-slice loop
    #     (landscape metrics, color metrics, diameter quantiles, root angles)
    # -------------------------------------------------------------------------
    needs_slice <- do_landscape || do_color || do_diam_q || do_angles
    
    if (needs_slice) {
      
      depth.slices <- sort(unique(terra::values(bm)))
      # With a single bin every pixel is already in the slice, so the masking
      # and trimming below would do nothing but copy full-size rasters.
      one_bin      <- length(depth.slices) == 1L
      
      lsm_names <- if (do_landscape)
        c("lsm_c_enn_mn", "lsm_l_joinent", "lsm_l_relmutinf", "lsm_l_np", "lsm_l_contag")
      else character(0)
      
      # Resample RGB to segmented extent once per image
      im.rgb.crop <- NULL
      if (do_color && !is.null(im.rgb))
        im.rgb.crop <- .safe("RGB resample", terra::resample(im.rgb, im, "bilinear"))
      
      # Pre-build deep_drive optimal-angle map once per image, reuse per slice
      bm_vals <- ang_vals <- gg_vals <- NULL
      
      if (do_angles && !is.null(angles_map) && !is.null(DepthMap)) {
        dd <- .safe("deep_drive (full image)", {
          adm <- terra::t(terra::flip(DepthMap))
          terra::ext(adm) <- terra::ext(angles_map)
          Rootopia::deep_drive(DepthMap = adm, AngleMap = angles_map, return = "all")
        })
        gg.full <- if (!is.null(dd)) dd$optimal_angle_map else NULL
        
        if (!is.null(gg.full)) {
          # Fix geometry if extents differ but dimensions match
          if (!terra::compareGeom(gg.full, bm, stopOnError = FALSE)) {
            if (all(dim(gg.full) == dim(bm))) {
              terra::ext(gg.full) <- terra::ext(bm)
            } else {
              message(sprintf("[Rootopia] %s: gg.full/bm dimension mismatch -- disabling calc_root_angles.", tube))
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
        
        im.sl <- im
        if (!one_bin) { im.sl[bm != d] <- NA; im.sl <- terra::trim(im.sl) }
        
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
        
        # --- Color metrics --------------------------------------------------
        rc <- pc <- data.frame()
        empty_col <- data.frame(rcc = NA, gcc = NA, bcc = NA, hue = NA,
                                saturation = NA, luminosity = NA,
                                red = NA, green = NA, blue = NA)
        if (do_color && !is.null(im.rgb.crop)) {
          cr <- .safe(sprintf("color (depth=%g)", d), {
            sl  <- im.rgb.crop
            if (!one_bin) { sl[bm != d] <- NA; sl <- terra::trim(sl) }
            ri  <- sl; ri[im.sl == 0] <- NA   # root pixels only
            pi_ <- sl; pi_[im.sl == 1] <- NA  # background pixels only
            rc_ <- tryCatch(Rootopia::tube_coloration(ri),  error = function(e) empty_col)
            pc_ <- tryCatch(Rootopia::tube_coloration(pi_), error = function(e) empty_col)
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
            sl_rd  <- rd.map
            if (!one_bin) { sl_rd[bm != d] <- NA; sl_rd <- terra::trim(sl_rd) }
            rd_v   <- terra::values(sl_rd, na.rm = FALSE)
            
            qv <- stats::quantile(rd_v, diameter_quantiles, na.rm = TRUE)

            res <- data.frame(depth = d)
            for (qi in seq_along(diameter_quantiles)) {
              res[[q_diam_names[qi]]] <- as.vector(qv[qi])
              res[[q_top_names[qi]]]  <- mean(rd_v[rd_v >= qv[qi]], na.rm = TRUE)
            }
            
            sl_rl <- root.length.map
            if (!one_bin) sl_rl[bm != d] <- NA
            for (ti in seq_along(thr_cm)) {
              ab  <- sl_rl; ab[rd.map < thr_cm[ti]] <- NA
              res[[paste0("rootlength.above.",    thr_names[ti])]] <-
                sum(terra::values(ab), na.rm = TRUE) / (dpi / 2.54)
              res[[paste0("avg.diameter.above.",  thr_names[ti])]] <-
                mean(rd_v[!is.na(rd_v) & rd_v >= thr_cm[ti]], na.rm = TRUE)
            }
            
            pk <- if (calc_modal_peaks) {
              rd_clean <- rd_v[!is.na(rd_v) & rd_v > 0]
              tryCatch({
                mp <- Rootopia::modal_peaks(rd_clean, display_type = "none",
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
            } else {
              data.frame(n.diameter.peaks = NA_real_, diameter.peak.1 = NA_real_,
                         diameter.peak.2  = NA_real_, diameter.peak.3  = NA_real_)
            }

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
        "[Rootopia] [%d/%d] %s | img: %.0fs | elapsed: %.0fs | remaining: ~%.0fs | done ~%s",
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
    warning(sprintf("[Rootopia] %d image(s) skipped entirely: %s",
                    length(failed_imgs),
                    paste(failed_imgs, collapse = ", ")),
            call. = FALSE)
  
  if (length(valid) == 0L) {
    warning("[Rootopia] No images were processed successfully.", call. = FALSE)
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
        dplyr::group_by(Tube) |>
        dplyr::filter(!is.na(depth) & !is.na(rootlength.density)) |>
        dplyr::summarize(
          mrd   = MRD(w = depth, roots = rootlength.density),
          # total.length.density: sum of (length density x bin size) across all bins
          # units: cm root per cm^2 (integrated over the full depth profile)
          total.length.density = sum(rootlength.density * depth_interval_cm, na.rm = TRUE),
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

      # Which optional columns exist is a property of the whole frame, not of
      # each group, so it is settled once here rather than re-tested per group.
      has <- function(col) col %in% names(r)

      r |>
        dplyr::group_by(Tube, depth) |>
        dplyr::filter(!is.na(depth) & !is.na(rootlength.density)) |>
        dplyr::summarize(

          # Fraction of the tube's total length density contributed by this bin
          rootlength.fraction = if (calc_density_metrics && has("total.length.density"))
            rootlength.density / total.length.density else NA_real_,

          # Joint entropy per unit root pixel density
          # (meaningful only when landscape metrics are available)
          ent_per_rootpx = if (calc_landscape_metrics && has("joinent"))
            joinent / rootpx.density else NA_real_,

          # Number of root patches per unit root pixel density
          patch_density_norm = if (has("np_density"))
            np_density / rootpx.density else NA_real_,

          # Mean within-bin diameter variance (averaged over any sub-grouping)
          mean.var.diameter = if (calc_diameter_stats && has("var.diameter"))
            mean(var.diameter, na.rm = TRUE) else NA_real_,

          # NB: rootsurface_rootvolume_ratio is computed per bin in section 3f
          # (directly on the diameter raster) so it is unbiased; see there.

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
      .msg("[Rootopia] Results saved to: %s", output_path)
    })
  }
  
  if (verbose) {
    total_min <- round(sum(img_times) / 60, 1)
    n_ok      <- n_images - length(failed_imgs)
    message(sprintf(
      "[Rootopia] Done. %d/%d images processed in %.1f min.",
      n_ok, n_images, total_min
    ))
  }
  
  invisible(root.depth.metrics)
}

# Column names referenced inside dplyr verbs above. Declaring them keeps
# R CMD check from reading non-standard evaluation as undefined globals.
utils::globalVariables(c(
  "Tube", "depth", "rootlength.density", "total.length.density",
  "joinent", "rootpx.density", "np_density", "var.diameter"
))


#' @rdname root_depth_metrics
#' @export
batch_root_traits <- root_depth_metrics

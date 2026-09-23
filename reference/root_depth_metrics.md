# Compute root traits over a depth profile from root scans

Processes a directory of root images and returns a tidy data frame of
root traits summarized per depth interval. Handles flatbed scans and
flat rhizotron windows by default, and cylindrical minirhizotron tubes
when tube parameters are supplied (see **Geometry**). Input may be
already segmented or a raw grayscale/RGB scan, which is binarized on the
way in (see **Binarization**).

`batch_root_traits()` is an alias for the same function.

Each metric group is toggled independently. If a block fails for one
image it is replaced with `NA` columns and a message is printed;
processing always continues to the next image. If an entire image cannot
be loaded it is dropped and listed in a warning at the end.

## Usage

``` r
root_depth_metrics(
  path_seg,
  path_skl = NULL,
  path_rgb = NULL,
  seg_file_index = NULL,
  skl_file_index = NULL,
  rgb_file_index = NULL,
  insertion_angles = NULL,
  soil_starts = 0,
  tube_names = NULL,
  session = "",
  dpi = 300,
  tube_diameter_cm = NULL,
  tube_center_offset = 0,
  depth_interval_cm = 5,
  bin_round = c("rounding", "floor", "ceiling"),
  rotation_fixed_width = 1800,
  binarize = "auto",
  binarize_threshold = 200,
  dark_roots = TRUE,
  seg_layer = NULL,
  clean_max_hole_size = 0,
  clean_max_artifact_size = 0,
  calc_root_pixels = TRUE,
  calc_root_length = TRUE,
  calc_diameter_stats = TRUE,
  calc_diameter_quantiles = FALSE,
  calc_modal_peaks = FALSE,
  calc_landscape_metrics = FALSE,
  calc_color_metrics = FALSE,
  calc_root_angles = FALSE,
  calc_root_order_metrics = FALSE,
  order_scheme = c("strahler_order", "branch_order", "root_order", "tip_order"),
  diam_weight = 0.5,
  prune_spur_length_cm = 0,
  prune_spur_iter = 1L,
  calc_density_metrics = TRUE,
  calc_distribution_indices = TRUE,
  calc_advanced_metrics = TRUE,
  diameter_thresholds = c(0.2, 0.5, 1),
  diameter_threshold_unit = "mm",
  diameter_quantiles = c(0.9, 0.95, 0.99),
  output_path = NULL,
  verbose = TRUE
)

batch_root_traits(
  path_seg,
  path_skl = NULL,
  path_rgb = NULL,
  seg_file_index = NULL,
  skl_file_index = NULL,
  rgb_file_index = NULL,
  insertion_angles = NULL,
  soil_starts = 0,
  tube_names = NULL,
  session = "",
  dpi = 300,
  tube_diameter_cm = NULL,
  tube_center_offset = 0,
  depth_interval_cm = 5,
  bin_round = c("rounding", "floor", "ceiling"),
  rotation_fixed_width = 1800,
  binarize = "auto",
  binarize_threshold = 200,
  dark_roots = TRUE,
  seg_layer = NULL,
  clean_max_hole_size = 0,
  clean_max_artifact_size = 0,
  calc_root_pixels = TRUE,
  calc_root_length = TRUE,
  calc_diameter_stats = TRUE,
  calc_diameter_quantiles = FALSE,
  calc_modal_peaks = FALSE,
  calc_landscape_metrics = FALSE,
  calc_color_metrics = FALSE,
  calc_root_angles = FALSE,
  calc_root_order_metrics = FALSE,
  order_scheme = c("strahler_order", "branch_order", "root_order", "tip_order"),
  diam_weight = 0.5,
  prune_spur_length_cm = 0,
  prune_spur_iter = 1L,
  calc_density_metrics = TRUE,
  calc_distribution_indices = TRUE,
  calc_advanced_metrics = TRUE,
  diameter_thresholds = c(0.2, 0.5, 1),
  diameter_threshold_unit = "mm",
  diameter_quantiles = c(0.9, 0.95, 0.99),
  output_path = NULL,
  verbose = TRUE
)
```

## Arguments

- path_seg:

  Character. Path to directory of binary segmented images
  (foreground/root pixel = 1, background = 0).

- path_skl:

  Character or `NULL`. Path to directory of skeletonised images
  (one-pixel-wide centerlines of roots), used for length, diameter,
  angle, and branching-order metrics. If `NULL` or a file is missing,
  the skeleton is computed internally via
  [`skeletonize_image()`](https://jcunow.github.io/Rootopia/reference/skeletonize_image.md).
  Default `NULL`.

- path_rgb:

  Character or `NULL`. Path to directory of blended RGB images, aligned
  to the segmented images. Required for color metrics. Default `NULL`.

- seg_file_index:

  Integer vector or `NULL`. Optional subset index applied to
  `list.files(path_seg)`, e.g. `37:72`. Default `NULL` (use all files).

- skl_file_index:

  Integer vector or `NULL`. Optional subset index applied to
  `list.files(path_skl)`. Default `NULL`.

- rgb_file_index:

  Integer vector or `NULL`. Optional subset index applied to
  `list.files(path_rgb)`. Default `NULL`.

- insertion_angles:

  Numeric or `NULL`. Insertion angle of the tube, in **degrees measured
  from horizontal** – `90` is a vertical tube, `45` a tube pushed in at
  45 degrees. This is the convention
  [`create_depthmap()`](https://jcunow.github.io/Rootopia/reference/create_depthmap.md)
  uses for `tilt`, and it is the opposite of the "degrees from vertical"
  wording used by some minirhizotron software, so convert before passing
  it in. Values must be strictly between 0 and 90. Supplying this
  argument switches the run to minirhizotron geometry (see
  **Geometry**). Default `NULL` (flatbed).

- soil_starts:

  Numeric. Pixel row (in the original, un-rotated image) at which the
  soil surface begins. Used to set the zero-depth reference. Default
  `0`.

- tube_names:

  Character or `NULL`. Sample/tube identifiers added as the `Tube`
  column. If `NULL`, names are derived from characters 3-5 from the
  right of the segmented file name, prefixed with `"T"` (e.g. `"T042"`).
  Adjust if your naming convention differs.

  Names must be **unique**, one per image, and the run stops if they are
  not. They key the tube-level joins, where two images sharing a name
  are read as one tube and their rows are multiplied out. The derived
  default collides whenever the files share a suffix – a directory of
  `"CLAS1A10-15fine.tif"` names every image `"Tfin"` – so flatbed scans
  generally need this argument. Default `NULL`.

- session:

  Character. Session or campaign label added as the `Session` column,
  e.g. `"2022_02"`. Default `""`.

- dpi:

  Numeric; image resolution in dots per inch

- tube_diameter_cm:

  Numeric or `NULL`. Inner diameter of the minirhizotron tube in
  **centimetres**. Passed to
  [`create_depthmap()`](https://jcunow.github.io/Rootopia/reference/create_depthmap.md)
  as `tube_thicc`, where it sets the amplitude and wavelength of the
  sinusoidal curvature correction. Supplying this argument switches the
  run to minirhizotron geometry (see **Geometry**). Default `NULL`
  (flatbed).

- tube_center_offset:

  Numeric in `[0, 1]`. Phase of the sinusoidal curvature correction:
  where the top of the tube falls across the image height, as a
  fraction. `0` (default) puts it at the first row. Set it if your
  scanner's rotational reference differs, or the depth assigned to a
  root will be offset by up to half a tube diameter. Minirhizotron only.

- depth_interval_cm:

  Numeric or `NULL`. Size of each depth bin in **centimetres**. Passed
  as `nn` to
  [`binning()`](https://jcunow.github.io/Rootopia/reference/binning.md).
  `NULL` switches on whole-image mode, where the scan is treated as a
  single bin and summarized in one row (see **Whole-image mode**).
  Default `5`.

- bin_round:

  Character. How
  [`binning()`](https://jcunow.github.io/Rootopia/reference/binning.md)
  assigns a depth to a bin: `"rounding"` (default,
  `nn * round(depth/nn)`), `"floor"`, or `"ceiling"`. **Note what
  "rounding" does to the top bin**: with `depth_interval_cm = 5` it
  spans 0-2.5 cm while every other bin spans 5 cm, because the label is
  the bin's center rather than its top edge. Per-bin densities are
  unaffected (they divide by each bin's own measured area), but `mrd`
  and `total.length.density` multiply by `depth_interval_cm` as though
  every bin were full width, so they are biased by the half-width top
  bin. `"floor"` gives the soil-science convention – 0-5, 5-10, labeled
  by the shallower edge – and is the better choice for a new analysis;
  the default is kept for continuity with existing ones.

- rotation_fixed_width:

  Numeric. Width in **rows** that each image is cropped to along the
  rotation axis, centered on the middle row, before any trait is
  measured (see
  [`rotation_censor()`](https://jcunow.github.io/Rootopia/reference/rotation_censor.md)).
  This trims the tube edges, where the curvature of the tube distorts
  what the scanner sees. An image with fewer rows than this cannot be
  cropped symmetrically, so
  [`rotation_censor()`](https://jcunow.github.io/Rootopia/reference/rotation_censor.md)
  clamps to the image bounds and says so – the image is then used at
  full width. Only applied under minirhizotron geometry; a flatbed scan
  is measured at full width. Default `1800`.

- binarize:

  Either `"auto"` (default), `TRUE`, or `FALSE`. `"auto"` thresholds an
  image only if at least one of its layers has more than two distinct
  values, so binary masks (single- or multi-layer) are left alone and
  raw scans are binarized. `TRUE` always thresholds, `FALSE` never does
  (a pixel is then root when it is non-zero in every layer). Note that
  `TRUE` on an image already coded 0/1 will *invert* it, because 0
  counts as dark; this is what `"auto"` exists to prevent.

- binarize_threshold:

  Numeric. Gray level at which a scan is cut into root and background,
  on the **0-255** scale (the RhizoVision Explorer convention). Images
  that load on a 0-1 scale get the same cut-off rescaled, so the number
  means the same thing either way. Default `200`.

- dark_roots:

  Logical. `TRUE` (default) means roots are *darker* than the
  background, as on a flatbed scan of washed roots on a white tray:
  pixels at or below `binarize_threshold` become root. `FALSE` inverts
  this for bright-roots-on-dark-background images. If more than half the
  pixels come out as root, the polarity is probably wrong and a warning
  says so.

- seg_layer:

  Integer or `NULL`. Which layer of the segmented image to measure.
  `NULL` (default) uses every layer: a single-layer image as is, the
  first three layers of a 3- or 4-layer image (any alpha band is
  dropped), and the first layer of a 2-layer image. Thresholded images
  are converted to grayscale with
  [`rgb2gray()`](https://jcunow.github.io/Rootopia/reference/rgb2gray.md);
  binary ones take root as non-zero in every layer. Set this if your
  files carry the segmentation in one specific band.

- clean_max_hole_size:

  Numeric. Fill enclosed background holes of up to this many pixels
  before anything is measured; `0` (default) fills none, `Inf` fills
  every hole. Pinholes inside a painted root read as background, which
  eats into the distance transform and so into diameter.

- clean_max_artifact_size:

  Numeric. Drop disconnected foreground blobs of up to this many pixels;
  `0` (default) drops none, `Inf` drops everything not touching the
  largest structure. Specks count as root area, as isolated skeleton
  pixels, and as root tips. Both use
  [`clean_image`](https://jcunow.github.io/Rootopia/reference/clean_image.md)
  and need the imager package.

- calc_root_pixels:

  Logical. Count `rootpx` (foreground pixels) and `voidpx` (background
  pixels) per depth bin. Required for density metrics. Default `TRUE`.

- calc_root_length:

  Logical. Estimate root length (cm) per depth bin from the skeleton
  using D8 connectivity (orthogonal steps = 1 pixel, diagonal steps =
  \\\sqrt{2}\\ pixels, isolated pixels = 1 pixel). Required for length
  density and angle metrics. Default `TRUE`.

- calc_diameter_stats:

  Logical. Compute per-bin mean, maximum, and variance of root
  diameter (cm) using the distance-transform approach in
  [`root_diameter()`](https://jcunow.github.io/Rootopia/reference/root_diameter.md),
  plus `root.surface.area` (lateral surface area, cm\\^2\\) and
  `root.volume` (cm\\^3\\) per bin, and `rootsurface_rootvolume_ratio`
  (cm\\^{-1}\\). `root.surface.area` is the curved wall wrapping each
  cylindrical root (its soil-contact area), *not* the flat root area
  visible in the image. `rootsurface_rootvolume_ratio` is the
  length-weighted mean of the *local* ratio \\2 / r_i\\ over skeleton
  pixels, so it is dominated by fine roots and is deliberately not equal
  to `root.surface.area / root.volume` (a bulk ratio dominated by thick
  roots); the two answer different questions. Default `TRUE`.

- calc_diameter_quantiles:

  Logical. Compute the diameter distribution percentiles set by
  `diameter_quantiles` per bin, conditional means above each quantile,
  and threshold-based root lengths (see `diameter_thresholds`). Default
  `FALSE`.

- calc_modal_peaks:

  Logical. Compute modal diameter peaks per bin via
  [`modal_peaks()`](https://jcunow.github.io/Rootopia/reference/modal_peaks.md)
  (`n.diameter.peaks`, `diameter.peak.1/2/3`). Auto-enables
  `calc_diameter_quantiles`, since it reuses the same per-bin diameter
  raster. **Slow**: one call per depth bin per image. Default `FALSE`.

- calc_landscape_metrics:

  Logical. Compute patch-level landscape metrics per depth bin via
  [`root_scape_metrics()`](https://jcunow.github.io/Rootopia/reference/root_scape_metrics.md):
  nearest-neighbor distance (`enn_mn`), joint entropy (`joinent`),
  relative mutual information (`relmutinf`), number of patches (`np`),
  and contagion (`contag`). **Slow**: one call per depth bin per image.
  Default `FALSE`.

- calc_color_metrics:

  Logical. Compute mean chromatic coordinates (rcc, gcc, bcc), hue,
  saturation, luminosity, and raw RGB channel means separately for root
  pixels and background pixels via
  [`tube_coloration()`](https://jcunow.github.io/Rootopia/reference/tube_coloration.md).
  Requires `path_rgb`. Default `FALSE`.

- calc_root_angles:

  Logical. Compute `deep_drive` (fraction of skeleton pixels whose D8
  flow direction matches the locally optimal downward direction) and
  `mean.steepness.angle` / `sd.steepness.angle` (degrees, 0 =
  horizontal, 90 = vertical). Uses
  [`deep_drive()`](https://jcunow.github.io/Rootopia/reference/deep_drive.md).
  Default `FALSE`.

- calc_root_order_metrics:

  Logical. Build a per-image branching-order graph via
  [`branch_order_map()`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md)
  and summarize it both per depth bin and per tube. Adds `mean.<scheme>`
  and `max.<scheme>` for the scheme named by `order_scheme` (e.g.
  `mean.strahler_order`), plus `mean.root_order` and
  `lateral_root_fraction` per depth bin, tube-level `main_root.*` /
  `lateral_roots.*` columns (length, diameter, branching frequency,
  etc.), and `n_root_orders` (the highest order found). Requires a
  skeleton. **Slow**: builds one segment graph per image. Default
  `FALSE`.

- order_scheme:

  Character. Which ordering labels the per-bin order columns and the
  main-root split. `"strahler_order"` (default) is the convention of the
  fine-root literature (Pregitzer et al. 2002; Fitter): every distal
  unbranched root is order 1, and the order rises only where two roots
  of equal order meet. `"tip_order"` is the same leaf-peeling with the
  order raised at *every* junction, so an axis carrying ten laterals
  reaches order 11. `"branch_order"` counts the other way – the thickest
  root of each component is 1 and its laterals 2 – and `"root_order"`
  gives each continuous root the maximum `tip_order` along it. All four
  are computed regardless; this only picks which one is reported. See
  [`branch_order_map`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md).

- diam_weight:

  Numeric \>= 0. At a junction, which two arms are read as the same root
  continuing: `straightness + diam_weight * diameter_similarity`. `0`
  uses the angle alone, larger values let thickness decide. Affects
  `root_order`, `branch_order` and the segment grouping, not
  `strahler_order` or `tip_order`. Default `0.5`.

- prune_spur_length_cm:

  Numeric. Remove terminal skeleton branches ("spurs") shorter than
  this, in **centimetres**, before any trait is measured; `0` (default)
  prunes nothing. Thinning leaves short stubs where roots are wide or
  ragged, and each one is counted as a root tip and adds a little
  length, so the pruning applies to the skeleton every metric is
  measured from – length and diameter as well as the order graph. Start
  around two to three times the width of your thickest root and check
  one image before trusting a batch. Uses
  [`prune_skeleton`](https://jcunow.github.io/Rootopia/reference/prune_skeleton.md).

- prune_spur_iter:

  Integer. Pruning passes, so that a spur exposed by removing another is
  caught too. More passes eat further into real roots. Default `1`.

- calc_density_metrics:

  Logical. Compute `rootpx.density` (percent root area cover) and
  `rootlength.density` (cm root length per cm^2 imaged area) per bin.
  Auto-enables `calc_root_pixels` and `calc_root_length`. Default
  `TRUE`.

- calc_distribution_indices:

  Logical. Compute tube-level indices: `mrd` (mean rooting depth) and
  `total.length.density` (summed length density over all bins, in cm
  root per cm^2 per cm depth). Auto-enables `calc_density_metrics`.
  Default `TRUE`.

- calc_advanced_metrics:

  Logical. Compute per-bin derived metrics: `rootlength.fraction` (each
  bin's length density as a fraction of the tube total) and
  `mean.var.diameter` (mean of within-bin diameter variance).
  Auto-enables `calc_distribution_indices` and `calc_diameter_stats`.
  Default `TRUE`.

- diameter_thresholds:

  Numeric vector. Diameter cut-offs for computing `rootlength.above.*`
  and `avg.diameter.above.*` columns. Units are set by
  `diameter_threshold_unit`. Default `c(0.2, 0.5, 1)`.

- diameter_threshold_unit:

  Character. Unit of `diameter_thresholds`: `"mm"` (default), `"cm"`, or
  `"px"`.

- diameter_quantiles:

  Numeric vector of probabilities (each strictly between 0 and 1) for
  the per-bin diameter percentiles computed when
  `calc_diameter_quantiles = TRUE`. Default `c(0.90, 0.95, 0.99)`.
  Output columns are named from the probabilities: e.g. `0.90` gives
  `rootdiameter.90` (the 90th percentile) and `avg.diameter.top10pct`
  (mean diameter above it).

- output_path:

  Character or `NULL`. If provided, the result is saved as an `.RData`
  file at this path. The exported object is named `root.depth.metrics`.
  Parent directories are created if they do not exist. Default `NULL`
  (no file written).

- verbose:

  Logical. Print per-image progress lines showing image index, per-image
  time, cumulative elapsed time, estimated remaining time, and predicted
  clock-time of completion. Default `TRUE`.

## Value

A data frame with one row per tube x depth-bin combination. Always
present columns: `Tube`, `Session`, `Plot`, `depth`. All other columns
depend on which metric groups are enabled; disabled or failed metrics
appear as `NA` rather than being absent. Returns `NULL` invisibly if
every image failed.

## Details

**Surface, volume, and their ratio.** For each skeleton pixel, the root
segment is modeled as a cylinder of length \\l_i\\ (one pixel edge in
cm) and radius \\r_i\\ (half the local diameter in cm). Its *lateral*
surface area is \\2 \pi r_i l_i\\ – the curved wall wrapping the root,
i.e. the area in contact with the soil, not the flat root area seen in
the image – and its volume is \\\pi r_i^2 l_i\\. `root.surface.area`
(cm\\^2\\) and `root.volume` (cm\\^3\\) are these quantities summed over
all skeleton pixels in the depth bin.

The per-pixel surface-to-volume ratio simplifies to \\2 / r_i\\, and
`rootsurface_rootvolume_ratio` (cm\\^{-1}\\) is the length-weighted mean
of \\2 / r_i\\ over the bin's skeleton pixels. Because it averages the
*local* ratio, it is dominated by fine roots (small \\r_i\\ give large
\\2 / r_i\\). This is deliberately not the same as the bulk ratio
`root.surface.area / root.volume`, which is dominated by thick roots
(they hold most of the volume); the two summarize different things and
will not match unless every root in the bin has the same diameter.

**Fault tolerance.** Every metric block is wrapped in `tryCatch`.
Failures produce a `[Rootopia] SKIPPED` message and `NA` values; they
never abort the run.

**Dependency resolution.** Enabling a higher-level metric silently
enables its prerequisites and prints a message listing what was
auto-enabled.

## Image ordering

Files are matched by position after
[`list.files()`](https://rdrr.io/r/base/list.files.html) sorts them
alphabetically. The segmented, skeleton, and RGB directories must
therefore contain files whose alphabetical order corresponds to the same
physical sample. Use `skl_file_index` and `rgb_file_index` to subset
those directories if necessary.

## Which paths are required for which metrics

- `path_seg`:

  Always required.

- `path_skl`:

  Optional. Used by `calc_root_length`, `calc_diameter_stats`,
  `calc_diameter_quantiles`, `calc_root_angles`, and
  `calc_root_order_metrics` when supplied. If `path_skl` is `NULL` or a
  skeleton file is missing for an image, the skeleton is computed
  internally from the segmented image via
  [`skeletonize_image()`](https://jcunow.github.io/Rootopia/reference/skeletonize_image.md).

- `path_rgb`:

  Required when `calc_color_metrics` is `TRUE`.

## Per-image metadata

All metadata arguments below accept either a single value (recycled to
all images) or a vector of length equal to the number of images.

## Geometry

There is no geometry switch. The geometry follows from whether you
supply tube parameters:

- Flatbed (default):

  Neither `insertion_angles` nor `tube_diameter_cm` is supplied. The
  scan is treated as a flat surface imaged head-on: no sinusoidal
  curvature correction (`sinoid = FALSE`), no foreshortening of the
  depth axis (one pixel along the image width is one pixel of depth),
  and no
  [`rotation_censor()`](https://jcunow.github.io/Rootopia/reference/rotation_censor.md)
  crop, since there is no tube interior to crop to.

- Minirhizotron:

  Either argument is supplied. The sinusoidal curvature correction is
  switched on, the depth axis is foreshortened by
  `sin(insertion_angles)`, and each image is cropped to
  `rotation_fixed_width` rows about its center. An argument you leave
  out falls back to a neutral default: `tube_diameter_cm = 7` and
  `insertion_angles = 90` (vertical tube), both announced in the run
  log.

In both cases depth runs along the image **width** (left to right),
which is the orientation minirhizotron scanners produce. A flatbed scan
with the soil surface at the top must be rotated 90 degrees before it is
passed in, or the depth profile will be built across the wrong axis.

## Whole-image mode

A flatbed scan of washed roots in a tray has no depth axis – where a
root lies on the tray says nothing about where it grew. Setting
`depth_interval_cm = NULL` treats each image as one bin, so the result
carries one row per image, with `depth = 0`, and the density metrics
become whole-scan quantities: `rootpx.density` is the percentage of the
scan covered by root and `rootlength.density` is cm of root per cm\\^2\\
of scanned area – how densely the tray was packed. Both already divide
by the bin's own measured pixel area (`rootpx + voidpx`), so there is
nothing to bin by and no depth map is built, unless
`calc_root_angles = TRUE` asks for one. Nothing is masked or trimmed per
slice either. On a large scan both are worth skipping: a depth map is a
full-size double raster. (`calc_root_length` still builds a flat surface
of its own to read flow directions from; that one belongs to the length
estimator, not to the binning, and is built either way.)

`calc_distribution_indices` and `calc_advanced_metrics` are switched off
in this mode: mean rooting depth, and each bin's share of the profile,
mean nothing when there is only one bin.

## Which parameters to set

Most runs set five or six of these arguments, and which five depends on
the geometry: a third of them describe a minirhizotron tube and do
nothing at all to a flatbed scan, with no warning when you tune one that
is inert. The *Batch Processing* vignette has the full tier list for
both geometries; in short:

- Flatbed:

  Set `dpi`, `tube_names` and `depth_interval_cm = NULL`, plus
  `binarize_threshold` and `dark_roots` for unsegmented scans.
  `insertion_angles`, `tube_diameter_cm`, `tube_center_offset`,
  `rotation_fixed_width`, `soil_starts` and `bin_round` do nothing here.

- Minirhizotron:

  Set `dpi`, `insertion_angles` (degrees from *horizontal*),
  `tube_diameter_cm`, `soil_starts`, `tube_names`, `depth_interval_cm`
  and `bin_round`, and look at `rotation_fixed_width` – how much of the
  tube's curved edge to trim is a property of your scanner that no
  default can guess.

The per-image arguments – `insertion_angles`, `soil_starts`,
`binarize_threshold`, `dark_roots`, `tube_names` – each take either one
value for the whole run or one per image, in the order of
`list.files(path_seg)`.

## Binarization

Flatbed scans are usually delivered as grayscale or RGB with the full
0-255 range, not as a segmented mask. Such an image is reduced to a
single grayscale layer
([`rgb2gray()`](https://jcunow.github.io/Rootopia/reference/rgb2gray.md)
for RGB input) and cut at `binarize_threshold`, in the same way
RhizoVision Explorer does it. Already-segmented input passes through
untouched: with the default `binarize = "auto"` the cut is only applied
when at least one layer of the image has more than two distinct values.
The check runs on each layer before any conversion to gray, so a
multi-class mask that is binary in every layer – e.g. RootDetector
output with white root, red second class and black background – is not
thresholded. In an image that is not thresholded, a pixel is root when
it is non-zero in every layer, i.e. white in an RGB mask. Use
`seg_layer` when root is coded differently.

## Examples

``` r
if (FALSE) { # \dontrun{
# Flatbed scans, not yet binarized -- the default path.
# No tube arguments, so this is flat geometry: no curvature correction,
# no foreshortening, no tube crop. Roots are dark on a bright tray, so
# everything at or below gray level 200 is taken as root.
result <- batch_root_traits(
  path_seg           = "scans/flatbed/tray_01/",
  dpi                = 600,
  binarize_threshold = 200,
  session            = "2024_spring"
)

# Minirhizotron -- supplying an insertion angle switches the geometry
result <- root_depth_metrics(
  path_seg         = "scans/segmented/2022_02/",
  path_skl         = "scans/skeleton/2022_02/",
  insertion_angles = tube_meta$angle,   # degrees from horizontal
  tube_diameter_cm = 7,
  session          = "2022_02"
)

# With diameter quantiles and root angle metrics
result <- root_depth_metrics(
  path_seg                = "scans/segmented/2022_02/",
  path_skl                = "scans/skeleton/2022_02/",
  path_rgb                = "scans/blended/2022_02/",
  rgb_file_index          = 37:72,
  insertion_angles        = tube_meta$angle,
  soil_starts             = tube_meta$soil_row,
  session                 = "2022_02",
  calc_diameter_quantiles = TRUE,
  calc_root_angles        = TRUE,
  calc_color_metrics      = TRUE,
  diameter_thresholds     = c(0.2, 0.5, 1),
  diameter_threshold_unit = "mm",
  output_path             = "output/root_metrics_2022_02.RData"
)

# Flat rhizotron window, already segmented. binarize = "auto" sees a
# two-level image and leaves it alone.
result <- root_depth_metrics(
  path_seg = "scans/segmented/rhizotron_A/",
  path_skl = "scans/skeleton/rhizotron_A/"
)
} # }
```

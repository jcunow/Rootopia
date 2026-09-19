# Compute root traits over a depth profile from root scans

Processes a directory of root images and returns a tidy data frame of
root traits summarized per depth interval. Handles flatbed scans and
flat rhizotron windows by default, and cylindrical minirhizotron tubes
when tube parameters are supplied (see **Geometry**). Input may be
already segmented or a raw greyscale/RGB scan, which is binarized on the
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
  depth_interval_cm = 5,
  rotation_fixed_width = 1800,
  binarize = "auto",
  binarize_threshold = 200,
  dark_roots = TRUE,
  seg_layer = NULL,
  calc_root_pixels = TRUE,
  calc_root_length = TRUE,
  calc_diameter_stats = TRUE,
  calc_diameter_quantiles = FALSE,
  calc_modal_peaks = FALSE,
  calc_landscape_metrics = FALSE,
  calc_color_metrics = FALSE,
  calc_root_angles = FALSE,
  calc_root_order_metrics = FALSE,
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
  depth_interval_cm = 5,
  rotation_fixed_width = 1800,
  binarize = "auto",
  binarize_threshold = 200,
  dark_roots = TRUE,
  seg_layer = NULL,
  calc_root_pixels = TRUE,
  calc_root_length = TRUE,
  calc_diameter_stats = TRUE,
  calc_diameter_quantiles = FALSE,
  calc_modal_peaks = FALSE,
  calc_landscape_metrics = FALSE,
  calc_color_metrics = FALSE,
  calc_root_angles = FALSE,
  calc_root_order_metrics = FALSE,
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
  Adjust if your naming convention differs. Default `NULL`.

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

- depth_interval_cm:

  Numeric or `NULL`. Size of each depth bin in **centimetres**. Passed
  as `nn` to
  [`binning()`](https://jcunow.github.io/Rootopia/reference/binning.md).
  `NULL` switches on whole-image mode, where the scan is treated as a
  single bin and summarised in one row (see **Whole-image mode**).
  Default `5`.

- rotation_fixed_width:

  Numeric. Width in **rows** that each image is cropped to along the
  rotation axis, centred on the middle row, before any trait is measured
  (see
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
  image only if it has more than two distinct grey levels, so binary
  masks are left alone and raw scans are binarized. `TRUE` always
  thresholds, `FALSE` never does (any non-zero pixel is then taken as
  root). Note that `TRUE` on an image already coded 0/1 will *invert*
  it, because 0 counts as dark; this is what `"auto"` exists to prevent.

- binarize_threshold:

  Numeric. Grey level at which a scan is cut into root and background,
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
  `NULL` (default) picks automatically: a single-layer image is used as
  is, and a 3- or 4-layer image is converted to greyscale with
  [`rgb2gray()`](https://jcunow.github.io/Rootopia/reference/rgb2gray.md).
  Set this if your files carry the segmentation in one specific band.

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
  and summarize it both per depth bin and per tube. Adds
  `mean.branch_order`, `max.branch_order`, `mean.root_order`, and
  `lateral_root_fraction` per depth bin, plus tube-level `main_root.*` /
  `lateral_roots.*` columns (length, diameter, branching frequency,
  etc., split by `order_metrics(..., focal = "thickest")`) and
  `n_root_orders` (the highest branch order found). Requires a skeleton.
  **Slow**: builds one segment graph per image. Default `FALSE`.

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
segment is modelled as a cylinder of length \\l_i\\ (one pixel edge in
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
(they hold most of the volume); the two summarise different things and
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
  `rotation_fixed_width` rows about its centre. An argument you leave
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

## Binarization

Flatbed scans are usually delivered as greyscale or RGB with the full
0-255 range, not as a segmented mask. Such an image is reduced to a
single greyscale layer
([`rgb2gray()`](https://jcunow.github.io/Rootopia/reference/rgb2gray.md)
for RGB input) and cut at `binarize_threshold`, in the same way
RhizoVision Explorer does it. Already-segmented input passes through
untouched: with the default `binarize = "auto"` the cut is only applied
when the image actually has more than two grey levels.

## Examples

``` r
if (FALSE) { # \dontrun{
# Flatbed scans, not yet binarized -- the default path.
# No tube arguments, so this is flat geometry: no curvature correction,
# no foreshortening, no tube crop. Roots are dark on a bright tray, so
# everything at or below grey level 200 is taken as root.
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

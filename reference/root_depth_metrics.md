# Compute root traits over a depth profile from segmented (mini)rhizotron images

Processes a set of segmented rhizotron or minirhizotron images and
returns a tidy data frame of root traits summarized per depth interval.
Supports both cylindrical tube geometry (minirhizotrons) and flat window
geometry (rhizotron panels).

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
  insertion_angles = 0,
  soil_starts = 0,
  tube_names = NULL,
  session = "",
  dpi = 300,
  tube_diameter_cm = 7,
  depth_interval_cm = 5,
  flat_geometry = FALSE,
  calc_root_pixels = TRUE,
  calc_root_length = TRUE,
  calc_diameter_stats = TRUE,
  calc_diameter_quantiles = FALSE,
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

  Numeric. Insertion angle of the tube or window from vertical, in
  **degrees**. `0` = perfectly vertical, `30` = tilted 30 degrees from
  vertical (a common minirhizotron angle). Used by
  [`create_depthmap()`](https://jcunow.github.io/Rootopia/reference/create_depthmap.md)
  to correct the depth scale. Default `0`.

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

  Numeric. Scanner resolution in **dots per inch**. Used to convert
  pixel distances to physical units (cm). Default `300`.

- tube_diameter_cm:

  Numeric. Inner diameter of the minirhizotron tube in **centimetres**.
  Passed to
  [`create_depthmap()`](https://jcunow.github.io/Rootopia/reference/create_depthmap.md)
  as `tube_thicc`. Its role in flat-window geometry
  (`flat_geometry = TRUE`) is uncertain; it likely has a default inside
  [`create_depthmap()`](https://jcunow.github.io/Rootopia/reference/create_depthmap.md)
  and may be ignored – leave at the default unless you know it matters
  for your setup. Default `7`.

- depth_interval_cm:

  Numeric. Size of each depth bin in **centimetres**. Passed as `nn` to
  [`binning()`](https://jcunow.github.io/Rootopia/reference/binning.md).
  Default `5`.

- flat_geometry:

  Logical. If `FALSE` (default), images are treated as cylindrical
  minirhizotron tubes and a sinusoidal depth correction is applied
  (`sinoid = TRUE` in
  [`create_depthmap()`](https://jcunow.github.io/Rootopia/reference/create_depthmap.md)).
  Set to `TRUE` for flat rhizotron windows (e.g. glass-fronted boxes)
  where no sinusoidal correction is needed (`sinoid = FALSE`). Default
  `FALSE`.

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
  [`root_diameter()`](https://jcunow.github.io/Rootopia/reference/root_diameter.md).
  Default `TRUE`.

- calc_diameter_quantiles:

  Logical. Compute the diameter distribution percentiles set by
  `diameter_quantiles` per bin, conditional means above each quantile,
  threshold-based root lengths (see `diameter_thresholds`), and modal
  diameter peaks via
  [`modal_peaks()`](https://jcunow.github.io/Rootopia/reference/modal_peaks.md).
  Default `FALSE`.

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

  Logical. Compute tube-level indices: `mrd` (mean rooting depth), `rpi`
  (root penetration index), and `total.length.density` (summed length
  density over all bins, in cm root per cm^2 per cm depth). Auto-enables
  `calc_density_metrics`. Default `TRUE`.

- calc_advanced_metrics:

  Logical. Compute per-bin derived metrics: `rootlength.fraction` (each
  bin's length density as a fraction of the tube total),
  `mean.var.diameter` (mean of within-bin diameter variance), and
  `rootsurface_rootvolume_ratio` (lateral surface area over cylinder
  volume, summed over skeleton pixels in the bin and expressed as cm^2
  per cm^3). Auto-enables `calc_distribution_indices` and
  `calc_diameter_stats`. Default `TRUE`.

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

**Surface-to-volume ratio.** For each skeleton pixel, the root segment
is modelled as a cylinder of length \\l_i\\ (the D8 path length of that
pixel in cm) and radius \\r_i\\ (half the local diameter in cm). The
lateral surface area is \\2 \pi r_i l_i\\ and the volume is \\\pi r_i^2
l_i\\. Their ratio simplifies to \\2 / r_i\\.
`rootsurface_rootvolume_ratio` is the length-weighted mean of \\2 /
r_i\\ over all skeleton pixels in the depth bin, in units of cm\\^{-1}\\
(cm^2 surface per cm^3 volume). Thicker roots have a smaller ratio; fine
roots have a larger ratio.

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

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal -- fast default metrics only
result <- root_depth_metrics(
  path_seg         = "scans/segmented/2022_02/",
  path_skl         = "scans/skeleton/2022_02/",
  insertion_angles = tube_meta$angle,
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

# Flat rhizotron window (no sinusoidal tube correction)
result <- root_depth_metrics(
  path_seg      = "scans/segmented/rhizotron_A/",
  path_skl      = "scans/skeleton/rhizotron_A/",
  flat_geometry = TRUE
)
} # }
```

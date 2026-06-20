# Batch Processing with root_depth_metrics()

## Batch Processing with `root_depth_metrics()`

### Who this vignette is for

If you have a folder of segmented minirhizotron images and want
depth-resolved root traits for all of them in one call — without writing
loops yourself —
[`root_depth_metrics()`](https://jcunow.github.io/Rootopia/reference/root_depth_metrics.md)
is the right starting point.

It handles file loading, depth-map construction, per-bin zonal
statistics, derived metrics, fault tolerance, and optional file output.
You toggle which metric groups to compute; everything else is taken care
of.

For a line-by-line explanation of the underlying steps, see the
[Minirhizotron
Scans](https://jcunow.github.io/Rootopia/articles/MinirhizotronScans_vignettes.md)
vignette. If multiple images come from the same tube, consider stitching
them beforehand [Image
Stitching](https://jcunow.github.io/Rootopia/articles/Stitching_vignette.md)

------------------------------------------------------------------------

### Minimum requirements

| What you need | Where it comes from |
|----|----|
| Segmented images (binary, root = 1) | [RootDetector](https://github.com/ExPlEcoGreifswald/RootDetector) or [RootPainter](https://github.com/Abe404/root_painter) |
| Skeletonised images (one-pixel centrelines) — *optional* | RootDetector (channel 2), or omit `path.skl` and [`root_depth_metrics()`](https://jcunow.github.io/Rootopia/reference/root_depth_metrics.md) will compute skeletons internally via [`skeletonize_image()`](https://jcunow.github.io/Rootopia/reference/skeletonize_image.md) — needed for length, diameter, angle, and branching-order metrics |
| RGB images aligned to segmented images | Original scan — needed for colour metrics only |
| Tube insertion angle per image | Field metadata |
| Soil-surface pixel row per image — *optional*, defaults to `0` | Field metadata or [`estimate_soil_surface()`](https://jcunow.github.io/Rootopia/reference/estimate_soil_surface.md); only needed if you want depths reported relative to the true soil surface rather than the top of the image |

Images in the three directories must be in the **same alphabetical
order** as each other; each position corresponds to the same physical
tube/window.

------------------------------------------------------------------------

### Installation

``` r

# install.packages("remotes")
remotes::install_github("jcunow/Rootopia")

library(Rootopia)
library(tidyverse)   # for downstream plotting
```

------------------------------------------------------------------------

### Quick start — default metrics only

This is the fastest way to get going. Only `path.seg` is required;
everything else has sensible defaults. `path.skl` is optional — if
omitted, skeletons needed for length, diameter, and angle metrics are
computed internally via
[`skeletonize_image()`](https://jcunow.github.io/Rootopia/reference/skeletonize_image.md).

``` r

result <- root_depth_metrics(
  path.seg = "scans/segmented/2022_02/",
  path.skl = "scans/skeleton/2022_02/",   # optional — computed internally if omitted
  session  = "2022_02"
)

head(result)
```

The result is a data frame with one row per tube × depth-bin
combination. Always-present columns: `Tube`, `Session`, `Plot`, `depth`.
Default metric columns: `rootpx`, `voidpx`, `rootlength`,
`avg.diameter`, `max.diameter`, `var.diameter`, `rootpx.density`,
`rootlength.density`, `mrd`, `rpi`, `total.length.density`,
`rootlength.fraction`, `mean.var.diameter`,
`rootsurface_rootvolume_ratio`.

------------------------------------------------------------------------

### Using per-image metadata

When tubes differ in insertion angle or soil-surface position, supply a
vector of values (one per image). Both accept a single recycled value
when all images share the same setting.

``` r

# Load your field metadata table
tube_meta <- read.csv("metadata/tube_metadata_2022_02.csv")
# Expected columns: filename, angle, soil_row, tube_id

result <- root_depth_metrics(
  path.seg         = "scans/segmented/2022_02/",
  path.skl         = "scans/skeleton/2022_02/",
  insertion_angles = tube_meta$angle,       # degrees from vertical
  soil_starts      = tube_meta$soil_row,    # pixel row of soil surface
  tube_names       = tube_meta$tube_id,     # e.g. "T042"
  session          = "2022_02",
  dpi              = 300
)
```

If your tube names follow a different convention from the default
(characters 3–5 from the right of the filename), always supply
`tube_names` explicitly.

------------------------------------------------------------------------

### Toggling metric groups

Each metric group is a logical argument. Disabled groups produce `NA`
columns rather than absent columns, so your downstream code can always
reference the same column names regardless of what was computed.

``` r

result <- root_depth_metrics(
  path.seg  = "scans/segmented/2022_02/",
  path.skl  = "scans/skeleton/2022_02/",
  path.rgb  = "scans/blended/2022_02/",
  session   = "2022_02",

  # --- core (on by default, shown explicitly for clarity) ---
  calc_root_pixels    = TRUE,
  calc_root_length    = TRUE,
  calc_diameter_stats = TRUE,

  # --- derived (on by default) ---
  calc_density_metrics      = TRUE,
  calc_distribution_indices = TRUE,
  calc_advanced_metrics     = TRUE,

  # --- extended (off by default — enable as needed) ---
  calc_diameter_quantiles = TRUE,   # 90th/95th/99th percentile + modal peaks
  calc_color_metrics      = TRUE,   # RGB colour per root vs background pixels
  calc_root_angles        = TRUE,   # deep_drive + steepness angle distribution
  calc_root_order_metrics = TRUE,   # branching order (main vs lateral roots)
  calc_landscape_metrics  = FALSE   # slow — enable only when needed
)
```

#### What each group adds

| Argument | New columns | Notes |
|----|----|----|
| `calc_root_pixels` | `rootpx`, `voidpx` | Needed for density metrics |
| `calc_root_length` | `rootlength` | Needed for density + angle metrics |
| `calc_diameter_stats` | `avg.diameter`, `max.diameter`, `var.diameter` | Needed for surface:volume ratio |
| `calc_diameter_quantiles` | `rootdiameter.90/95/99`, `avg.diameter.top*pct`, `rootlength.above.*`, `n.diameter.peaks`, `diameter.peak.*` | See `diameter_quantiles`, `diameter_thresholds` |
| `calc_landscape_metrics` | `enn_mn`, `joinent`, `relmutinf`, `np`, `contag`, `np_density` | **Slow** — one call per depth bin |
| `calc_color_metrics` | `rcc_root`, `gcc_root`, …, `rcc_bg`, `gcc_bg`, … | Requires `path.rgb` |
| `calc_root_angles` | `deep_drive`, `mean.steepness.angle`, `sd.steepness.angle` | Requires `calc_root_length` |
| `calc_root_order_metrics` | Per bin: `mean.branch_order`, `max.branch_order`, `mean.root_order`, `lateral_root_fraction`. Per tube: `main_root.*` / `lateral_roots.*` (length, diameter, branching frequency, …), `n_root_orders` | **Slow** — builds one segment graph per image. Requires `calc_root_length` |
| `calc_density_metrics` | `rootpx.density`, `rootlength.density` | Auto-enables pixels + length |
| `calc_distribution_indices` | `mrd`, `rpi`, `total.length.density` | Tube-level summary |
| `calc_advanced_metrics` | `rootlength.fraction`, `mean.var.diameter`, `rootsurface_rootvolume_ratio` | Auto-enables distribution indices + diameter stats |

------------------------------------------------------------------------

### Diameter thresholds

When `calc_diameter_quantiles = TRUE`, you can specify custom diameter
thresholds that separate fine roots from coarse roots:

``` r

result <- root_depth_metrics(
  path.seg                = "scans/segmented/2022_02/",
  path.skl                = "scans/skeleton/2022_02/",
  calc_diameter_quantiles = TRUE,
  diameter_thresholds     = c(0.2, 0.5, 1.0),   # fine / medium / coarse
  diameter_threshold_unit = "mm"                 # "mm", "cm", or "px"
)
# Resulting columns include: rootlength.above.0.2mm, rootlength.above.0.5mm, ...
```

------------------------------------------------------------------------

### Root branching order (main vs. lateral roots)

`calc_root_order_metrics = TRUE` builds a segment graph from the
skeleton of each image (via
[`branch_order_map()`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md))
and classifies every root segment by its **branch order** — the
thickest, most central root in each connected component is order 1 (the
“main” axis), its laterals are order 2, their laterals order 3, and so
on. See the [Minirhizotron
Scans](https://jcunow.github.io/Rootopia/articles/MinirhizotronScans_vignettes.html#6-root-branching-order)
vignette for a worked example of the underlying pipeline.

``` r

result <- root_depth_metrics(
  path.seg                = "scans/segmented/2022_02/",
  path.skl                = "scans/skeleton/2022_02/",
  session                 = "2022_02",
  calc_root_order_metrics = TRUE
)
```

This adds two kinds of columns:

**Per depth bin** — how branching order varies with depth:

| Column | Meaning |
|----|----|
| `mean.branch_order` | Mean branch order of skeleton pixels in the bin |
| `max.branch_order` | Highest branch order reached in the bin |
| `mean.root_order` | Mean root order (max tip-order along each root) in the bin |
| `lateral_root_fraction` | Fraction of skeleton pixels belonging to branch order \> 1 (laterals) |

**Per tube** (repeated on every row of that tube) — main root
vs. lateral roots, split by `order_metrics(..., focal = "thickest")`:

| Column prefix | Meaning |
|----|----|
| `main_root.*` | Totals/means for the thickest (order-1) root(s): `total_length`, `mean_diameter`, `median_diameter`, `n_segments`, `n_tips`, `n_branch_points`, `length_fraction`, `mean_segment_length`, `branching_frequency` |
| `lateral_roots.*` | The same set of statistics for everything else (order \>= 2) |
| `n_root_orders` | Highest branch order found in the image |

``` r

library(ggplot2)

ggplot(result, aes(x = depth, y = lateral_root_fraction, colour = Tube)) +
  geom_line() +
  geom_point() +
  coord_flip() +
  scale_x_reverse() +
  theme_minimal() +
  labs(
    title = "Lateral Root Fraction by Depth",
    x     = "Soil Depth (cm)",
    y     = "Fraction of skeleton pixels with branch order > 1"
  )
```

> **Slow.** This builds one segment graph per image, which is
> considerably more expensive than the other metric groups. Disable it
> (the default) for quick exploratory runs and enable it only when you
> need architectural detail.

------------------------------------------------------------------------

### Flat rhizotron windows

For glass-fronted rhizotron boxes (flat geometry, no tube curvature),
set `flat_geometry = TRUE`. This disables the sinusoidal depth
correction that is otherwise applied to cylindrical minirhizotron tubes.

``` r

result <- root_depth_metrics(
  path.seg      = "scans/segmented/rhizotron_A/",
  path.skl      = "scans/skeleton/rhizotron_A/",
  flat_geometry = TRUE,
  insertion_angles = 0   # vertical windows
)
```

------------------------------------------------------------------------

### Subsetting files

If your skeleton and RGB directories contain more files than the
segmented directory (e.g. multiple sessions in one folder), use
`*_file_index` to select the matching subset:

``` r

# Use only files 37–72 from the RGB directory
result <- root_depth_metrics(
  path.seg       = "scans/segmented/2022_02/",
  path.skl       = "scans/skeleton/2022_02/",
  path.rgb       = "scans/blended/all_sessions/",
  rgb_file_index = 37:72,
  session        = "2022_02"
)
```

------------------------------------------------------------------------

### Saving output

Pass `output_path` to save the result as an `.RData` file (the object is
named `root.depth.metrics` inside the file). The directory is created if
it does not exist.

``` r

result <- root_depth_metrics(
  path.seg    = "scans/segmented/2022_02/",
  path.skl    = "scans/skeleton/2022_02/",
  session     = "2022_02",
  output_path = "output/root_metrics_2022_02.RData"
)

# Reload later:
# load("output/root_metrics_2022_02.RData")
# head(root.depth.metrics)
```

------------------------------------------------------------------------

### Progress and fault tolerance

By default (`verbose = TRUE`) the function prints a progress line after
each image showing per-image time, cumulative elapsed time, estimated
remaining time, and predicted clock-time of completion:

    [Rootopia] [3/48] T042 | img: 12s | elapsed: 38s | remaining: ~546s | done ~14:22

If a metric block fails for a single image (e.g. a corrupted skeleton
file), it is skipped with a `[Rootopia] SKIPPED` message and the
corresponding columns are filled with `NA`. Processing always continues.
Images that cannot be loaded at all are skipped entirely and listed in a
warning at the end.

To suppress all progress output:

``` r

result <- root_depth_metrics(
  path.seg = "scans/segmented/2022_02/",
  path.skl = "scans/skeleton/2022_02/",
  verbose  = FALSE
)
```

------------------------------------------------------------------------

### Example: visualising depth profiles

``` r

library(ggplot2)

# Root length density by depth, one line per tube
ggplot(result, aes(x = depth, y = rootlength.density, colour = Tube)) +
  geom_line() +
  geom_point() +
  coord_flip() +
  scale_x_reverse() +
  theme_minimal() +
  labs(
    title = "Root Length Density by Depth",
    x     = "Soil Depth (cm)",
    y     = "Root Length Density (cm root / cm² area)"
  )

# Average root diameter by depth
ggplot(result, aes(x = depth, y = avg.diameter)) +
  geom_smooth(method = "loess", se = TRUE) +
  coord_flip() +
  scale_x_reverse() +
  theme_minimal() +
  labs(
    title = "Average Root Diameter with Depth",
    x     = "Soil Depth (cm)",
    y     = "Average Diameter (cm)"
  )

# mrd distribution across tubes
ggplot(result |> dplyr::distinct(Tube, mrd), aes(x = mrd)) +
  geom_histogram(bins = 15, fill = "steelblue4", colour = "white") +
  theme_minimal() +
  labs(
    title = "Mean Root Depth Distribution (MRD)",
    x     = "MRD",
    y     = "Count"
  )
```

------------------------------------------------------------------------

### Full example with all options

``` r

tube_meta <- read.csv("metadata/tube_metadata_2022_02.csv")

result <- root_depth_metrics(

  # Paths
  path.seg       = "scans/segmented/2022_02/",
  path.skl       = "scans/skeleton/2022_02/",
  path.rgb       = "scans/blended/2022_02/",

  # Per-image metadata
  insertion_angles = tube_meta$angle,
  soil_starts      = tube_meta$soil_row,
  tube_names       = tube_meta$tube_id,
  session          = "2022_02",

  # Geometry
  dpi               = 300,
  tube_diameter_cm  = 7,
  depth_interval_cm = 5,
  flat_geometry     = FALSE,

  # Core metrics
  calc_root_pixels    = TRUE,
  calc_root_length    = TRUE,
  calc_diameter_stats = TRUE,

  # Extended metrics
  calc_diameter_quantiles = TRUE,
  calc_color_metrics      = TRUE,
  calc_root_angles        = TRUE,
  calc_root_order_metrics = TRUE,    # main vs. lateral root architecture
  calc_landscape_metrics  = FALSE,   # enable if you need patch metrics

  # Derived metrics
  calc_density_metrics      = TRUE,
  calc_distribution_indices = TRUE,
  calc_advanced_metrics     = TRUE,

  # Diameter thresholds (for quantile block)
  diameter_thresholds     = c(0.2, 0.5, 1.0),
  diameter_threshold_unit = "mm",

  # Output
  output_path = "output/root_metrics_2022_02.RData",
  verbose     = TRUE
)
```

------------------------------------------------------------------------

### What to read next

- [Minirhizotron
  Scans](https://jcunow.github.io/Rootopia/articles/MinirhizotronScans_vignettes.md)
  — step-by-step explanation of each underlying function
- [Flatbed
  Scans](https://jcunow.github.io/Rootopia/articles/FlatBedScans_vignettes.md)
  — trait extraction from flatbed scanner images (no depth dimension)
- [Rotation
  Bias](https://jcunow.github.io/Rootopia/articles/Rotation_Bias_vignettes.md)
  — correcting for tube rotation artefacts before analysis
- [Function
  reference](https://jcunow.github.io/Rootopia/reference/index.md) —
  full documentation for every exported function

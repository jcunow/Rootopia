# Classify peat material classes from a minirhizotron RGB raster

Assigns each pixel of an RGB `SpatRaster` to a peat material class (e.g.
dark peat, red peat, root, silver tape, coarse debris) by nearest-
centroid assignment in CIE LAB colour space. Pixels beyond the per-class
distance threshold are labelled "unclassified".

## Usage

``` r
classify_peat_rgb(
  img,
  centroids = .default_peat_centroids(),
  downsample_fact = NULL,
  compute_metrics = TRUE,
  verbose = TRUE
)
```

## Arguments

- img:

  A `SpatRaster` with at least 3 layers interpreted as R, G, B (in that
  order). Values may be 0–255 or 0–1 (auto-detected).

- centroids:

  A `data.frame` with columns `class`, `L`, `A`, `B`, `MAX_DIST`.
  Defaults to a set of centroids calibrated on Oulanka 2023
  minirhizotron scans — see
  [`build_peat_centroids`](https://jcunow.github.io/RootScanR/reference/build_peat_centroids.md)
  to derive centroids for your own data.

- downsample_fact:

  Integer spatial aggregation factor applied before classification for
  speed. `NULL` (default) uses full resolution. The output map is always
  disaggregated back to match the input resolution and extent exactly
  (nearest-neighbour, no interpolation).

- compute_metrics:

  Logical. Compute per-class pixel counts, area fractions, LAB
  statistics, and mean distance to centroid. Default `TRUE`. Set `FALSE`
  in tight batch loops where only the map is needed.

- verbose:

  Logical. Print progress messages and summary table. Default `TRUE`.

## Value

A named list with elements:

- `map`:

  A `SpatRaster` of integer class IDs with factor levels set to class
  names. Level 0 = "unclassified". Same CRS, extent, and resolution as
  `img`. Use directly with
  [`terra::zonal()`](https://rspatial.github.io/terra/reference/zonal.html),
  [`terra::mask()`](https://rspatial.github.io/terra/reference/mask.html),
  [`terra::freq()`](https://rspatial.github.io/terra/reference/freq.html).

- `metrics`:

  A `data.frame` with per-class pixel counts, area fractions (%), LAB
  and RGB means and SDs, mean distance to centroid, and the actual mean
  colour rendered as hex. `NULL` if `compute_metrics = FALSE`.

- `inter_dist`:

  Numeric matrix of pairwise LAB distances between class centroids.

- `centroids`:

  The centroid table used (useful when the default was applied).

## See also

[`build_peat_centroids`](https://jcunow.github.io/RootScanR/reference/build_peat_centroids.md),
[`plot_peat_classification`](https://jcunow.github.io/RootScanR/reference/plot_peat_classification.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(terra)
img    <- rast("scan.tiff")
result <- classify_peat_rgb(img)

# Access outputs
terra::plot(result$map)
result$metrics

# Zonal stats with a depth-band raster
zones  <- rast("depth_zones.tif")
zstats <- terra::zonal(result$map, zones, fun = "freq")

# Custom centroids
cents  <- build_peat_centroids(my_picks, my_max_dist)
result <- classify_peat_rgb(img, centroids = cents)
} # }
```

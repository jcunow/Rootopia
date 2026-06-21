# Classify soil material classes from a minirhizotron RGB raster

Assigns each pixel of an RGB `SpatRaster` to a class (e.g. dark soil,
red soil, root, silver tape, coarse debris) by nearest- centroid
assignment in CIE LAB color space. Pixels beyond the per-class distance
threshold are labeled "unclassified".

## Usage

``` r
classify_soil_rgb(
  img,
  centroids = .default_soil_centroids(),
  downsample_fact = NULL,
  compute_metrics = TRUE,
  verbose = TRUE
)
```

## Arguments

- img:

  An RGB image with at least 3 layers/channels interpreted as R, G, B
  (in that order). Values may be 0-255 or 0-1 (auto-detected). Converted
  internally via
  [`load_flexible_image()`](https://jcunow.github.io/Rootopia/reference/load_flexible_image.md),
  so any of its supported formats are accepted: file path (.jpg, .jpeg,
  .png, .tif, .tiff, .bmp), `SpatRaster`, `RasterLayer`/`RasterBrick`,
  `cimg`, `magick-image`, `matrix`, or `array`.

- centroids:

  A `data.frame` with columns `class`, `L`, `A`, `B`, `MAX_DIST` (see
  **Centroid table format** below). Defaults to a set of centroids
  calibrated on Oulanka 2023 minirhizotron scans – see **Building your
  own centroids** below and
  [`build_soil_centroids`](https://jcunow.github.io/Rootopia/reference/build_soil_centroids.md)
  to derive centroids for your own data.

- downsample_fact:

  Integer spatial aggregation factor applied before classification for
  speed. `NULL` (default) uses full resolution. The output map is always
  disaggregated back to match the input resolution and extent exactly
  (nearest-neighbor, no interpolation).

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
  color rendered as hex. `NULL` if `compute_metrics = FALSE`.

- `inter_dist`:

  Numeric matrix of pairwise LAB distances between class centroids.

- `centroids`:

  The centroid table used (useful when the default was applied).

## Centroid table format

The `centroids` table has one row per material class and these columns:

- `class`:

  Character. Name of the class (e.g. `"dark_soil"`, `"root"`). Becomes
  the factor level in `result$map` and the `class` column of
  `result$metrics`.

- `L`, `A`, `B`:

  Numeric. The class centroid's coordinates in CIE LAB color space (D65
  illuminant): `L` is lightness (0 = black, 100 = white), `A` is the
  green-red axis (negative = green, positive = red), and `B` is the
  blue-yellow axis (negative = blue, positive = yellow). These are the
  \*mean\* LAB values of the representative pixels for that class.

- `MAX_DIST`:

  Numeric. The per-class assignment radius, in Euclidean LAB units. A
  pixel is assigned to the class whose centroid is nearest in LAB space,
  but only if that distance is `<= MAX_DIST`; otherwise the pixel is
  labeled `"unclassified"`. Larger values classify more pixels but risk
  merging visually distinct materials; smaller values leave more pixels
  unclassified. Typical values are roughly 10-30.

Pixel RGB values are converted to LAB internally before distances are
computed, so `L`/`A`/`B` are never raw RGB numbers.

## Building your own centroids (picks)

The default `centroids` table was calibrated on one specific scanner and
site, so for other data you should derive your own via
[`build_soil_centroids`](https://jcunow.github.io/Rootopia/reference/build_soil_centroids.md).

[`build_soil_centroids()`](https://jcunow.github.io/Rootopia/reference/build_soil_centroids.md)
takes `picks`: a named list with one element per material class. Each
element is a numeric matrix with exactly 3 columns (R, G, B in 0-255),
where every row is one color sample believed to belong to that class.
Classes may have different numbers of rows. The function converts each
matrix to LAB, averages it to a single centroid, and returns a
`data.frame` in the same format as `.default_soil_centroids()` – ready
to pass straight back into `classify_soil_rgb(centroids = ...)`.

The simplest approach is to read representative RGB values off your scan
(e.g. using an image viewer's color picker) and enter them directly:


    picks <- list(
      dark_soil = matrix(c( 28,  22,  18,
                             32,  26,  21,
                             25,  19,  15), ncol = 3, byrow = TRUE),
      red_soil  = matrix(c( 80,  45,  35,
                             75,  42,  31), ncol = 3, byrow = TRUE),
      root      = matrix(c(180, 160, 130,
                           175, 155, 125,
                           185, 165, 135,
                           178, 158, 128), ncol = 3, byrow = TRUE),
      tape      = matrix(c(200, 205, 210,
                           195, 200, 205), ncol = 3, byrow = TRUE),
      debris    = matrix(c(100,  70,  45,
                            95,  65,  40), ncol = 3, byrow = TRUE)
    )

    max_dist <- c(dark_soil = 14, red_soil = 14, root = 26,
                  tape = 28, debris = 11)

    cents  <- build_soil_centroids(picks, max_dist)   # prints diagnostics
    result <- classify_soil_rgb(img, centroids = cents)

Alternatively, picks can be extracted from known representative patches
of your image. **Important:**
[`build_soil_centroids()`](https://jcunow.github.io/Rootopia/reference/build_soil_centroids.md)
always treats pick values as 0-255 (no auto-scaling). If your raster
stores values in 0-1, multiply by 255 before building picks so that the
centroids and the pixels seen inside `classify_soil_rgb()` are on the
same scale.

[`build_soil_centroids()`](https://jcunow.github.io/Rootopia/reference/build_soil_centroids.md)
prints diagnostics (intra-class spread, inter-class distances,
`MAX_DIST` coverage) to help you sanity-check your choices. Provide a
class for *every* material present in your scans: because classification
is nearest-centroid, any material without its own class is snapped into
whichever defined class is closest.

## See also

[`build_soil_centroids`](https://jcunow.github.io/Rootopia/reference/build_soil_centroids.md),
[`plot_soil_classification`](https://jcunow.github.io/Rootopia/reference/plot_soil_classification.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(terra)
img    <- rast("scan.tiff")
result <- classify_soil_rgb(img)

# Access outputs
terra::plot(result$map)
result$metrics

# Zonal stats with a depth-band raster
zones  <- rast("depth_zones.tif")
zstats <- terra::zonal(result$map, zones, fun = "freq")

# Custom centroids -- see "Building your own centroids (picks)" above for
# how to construct `picks` and `max_dist`
cents  <- build_soil_centroids(picks, max_dist)
result <- classify_soil_rgb(img, centroids = cents)
} # }
```

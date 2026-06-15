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

  An RGB image with at least 3 layers/channels interpreted as R, G, B
  (in that order). Values may be 0-255 or 0-1 (auto-detected). Converted
  internally via
  [`load_flexible_image()`](https://jcunow.github.io/RootScanR/reference/load_flexible_image.md),
  so any of its supported formats are accepted: file path (.jpg, .jpeg,
  .png, .tif, .tiff, .bmp), `SpatRaster`, `RasterLayer`/`RasterBrick`,
  `cimg`, `magick-image`, `matrix`, or `array`.

- centroids:

  A `data.frame` with columns `class`, `L`, `A`, `B`, `MAX_DIST` (see
  **Centroid table format** below). Defaults to a set of centroids
  calibrated on Oulanka 2023 minirhizotron scans – see **Building your
  own centroids** below and
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

## Centroid table format

The `centroids` table has one row per material class and these columns:

- `class`:

  Character. Name of the class (e.g. `"dark_peat"`, `"root"`). Becomes
  the factor level in `result$map` and the `class` column of
  `result$metrics`.

- `L`, `A`, `B`:

  Numeric. The class centroid's coordinates in CIE LAB colour space (D65
  illuminant): `L` is lightness (0 = black, 100 = white), `A` is the
  green-red axis (negative = green, positive = red), and `B` is the
  blue-yellow axis (negative = blue, positive = yellow). These are the
  \*mean\* LAB values of the representative pixels for that class.

- `MAX_DIST`:

  Numeric. The per-class assignment radius, in Euclidean LAB units. A
  pixel is assigned to the class whose centroid is nearest in LAB space,
  but only if that distance is `<= MAX_DIST`; otherwise the pixel is
  labelled `"unclassified"`. Larger values classify more pixels but risk
  merging visually distinct materials; smaller values leave more pixels
  unclassified. Typical values are roughly 10-30.

Pixel RGB values are converted to LAB internally before distances are
computed, so `L`/`A`/`B` are never raw RGB numbers.

## Building your own centroids (picks)

The default `centroids` table was calibrated on one specific scanner and
site, so for other datasets you should derive your own using
[`build_peat_centroids`](https://jcunow.github.io/RootScanR/reference/build_peat_centroids.md).

That function requires `picks`: a named list where each element is an
`N x 3` numeric matrix of RGB values (0–255), with one row per manually
selected pixel that is known to belong to a given material class. The
more representative pixels supplied for each class, the more robust the
resulting centroid will be.


    picks <- list(
      dark_peat = rbind(
        c(32, 28, 25),
        c(35, 30, 27),
        c(29, 25, 22),
        c(31, 27, 24)
      ),
      root = rbind(
        c(215, 180, 130),
        c(220, 185, 135),
        c(210, 175, 128),
        c(218, 182, 132)
      )
    )

    max_dist <- c(
      dark_peat = 14,
      root      = 26
    )

    cents <- build_peat_centroids(picks, max_dist)

    result <- classify_peat_rgb(
      img,
      centroids = cents
    )

Each row of a pick matrix represents a single RGB observation:


    #      R    G    B
    c(215, 180, 130)

Typically, picks are obtained by interactively sampling representative
pixels from an image and grouping them by material class. The resulting
centroids are calculated in CIE LAB colour space, not RGB space.

[`build_peat_centroids()`](https://jcunow.github.io/RootScanR/reference/build_peat_centroids.md)
converts the supplied RGB picks to LAB, computes a centroid for each
class, and reports intra-class spread and inter-class distances to help
assess whether the chosen `MAX_DIST` values are appropriate.

## See also

[`build_peat_centroids`](https://jcunow.github.io/RootScanR/reference/build_peat_centroids.md),
[`plot_peat_classification`](https://jcunow.github.io/RootScanR/reference/plot_peat_classification.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(terra)
img    <- seg_Oulanka2023_Session01_T067
result <- classify_peat_rgb(img)

# Access outputs
terra::plot(result$map)
result$metrics


# Custom centroids -- see "Building your own centroids (picks)" above for
# how to construct `picks` and `max_dist`
cents  <- build_peat_centroids(picks, max_dist)
result <- classify_peat_rgb(img, centroids = cents)
} # }
```

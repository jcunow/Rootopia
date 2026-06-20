# Build soil class centroids from manual RGB colour picks

Converts raw RGB pixel picks (collected e.g. in QGIS, FIJI, or ImageJ)
to CIE LAB centroids and returns a centroid table ready to pass to
[`classify_soil_rgb`](https://jcunow.github.io/RootScanR/reference/classify_soil_rgb.md).
Also prints diagnostic summaries including intra-class spread, per-class
coverage, and inter-class distance warnings.

## Usage

``` r
build_soil_centroids(picks, max_dist, prior = NULL, alpha = 0, verbose = TRUE)
```

## Arguments

- picks:

  A named list of RGB pick matrices. Each element corresponds to one
  class and must be a numeric matrix with 3 columns (R, G, B), values
  0-255, with one row per pick. Names become class names in the output.
  See
  [`classify_soil_rgb`](https://jcunow.github.io/RootScanR/reference/classify_soil_rgb.md)'s
  **Building your own centroids (picks)** section for a worked example
  of constructing these matrices from an image (e.g. by cropping
  representative patches and calling
  [`terra::values()`](https://rspatial.github.io/terra/reference/values.html)).

- max_dist:

  A named numeric vector of per-class LAB distance thresholds, matched
  by name to `picks`. Pixels further than this from a class centroid
  cannot be assigned to that class.

- prior:

  Optional `data.frame` of existing centroids (same format as the output
  of this function, or `.default_soil_centroids()`). When supplied, the
  new centroids derived from `picks` are blended with the prior
  centroids using `alpha`. Only classes present in both `picks` and
  `prior` are blended; new classes in `picks` that are absent from
  `prior` are added as-is.

- alpha:

  Numeric in \[0, 1\]. Blend weight for the prior centroids. `alpha = 0`
  (default) ignores the prior entirely – centroids are derived purely
  from the new picks. `alpha = 1` returns the prior unchanged.
  `alpha = 0.5` weights old and new equally. Decrease `alpha`
  progressively as you collect more new picks to gradually shift
  calibration toward the new dataset.

- verbose:

  Logical. Print per-class summaries and inter-class distance matrix.
  Default `TRUE`.

## Value

A `data.frame` with columns `class`, `L`, `A`, `B`, `MAX_DIST`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Clean break -- new picks only
cents <- build_soil_centroids(new_picks, max_dist)

# Blend: 30% old calibration, 70% new picks
cents <- build_soil_centroids(new_picks, max_dist,
                              prior = .default_soil_centroids(),
                              alpha = 0.3)

# Iterative refinement across sessions:
# session 1
cents <- build_soil_centroids(picks_s1, max_dist)
# session 2 -- downweight session 1 to 20%
cents <- build_soil_centroids(picks_s2, max_dist, prior = cents, alpha = 0.2)
# session 3 -- downweight accumulated prior to 10%
cents <- build_soil_centroids(picks_s3, max_dist, prior = cents, alpha = 0.1)
} # }
```

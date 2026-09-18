# Rasterise a per-segment value onto the image grid

Writes any per-segment column (default the order class) back onto the
full image grid, aligned to `template`, for masking and zonal
statistics. Background pixels are `NA`.

## Usage

``` r
order_classification_map(et, template, value = "branch_order")
```

## Arguments

- et:

  An `edges` table carrying `attr(., "segments")` and
  `attr(., "crop_offset")` (run with `keep_segments = TRUE`).

- template:

  `SpatRaster` (or matrix) defining the output grid/extent.

- value:

  Column of `et` to rasterise (e.g. `"branch_order"`,
  `"mean_diameter"`).

## Value

A `SpatRaster` (or matrix) of `value` per root pixel, `NA` elsewhere,
aligned to `template`.

## Details

Only pixels that belong to a segment are painted. Junction contraction
dissolves the interior of each branch-point cluster, so the map carries
one unpainted pixel per branch point and a pixel count taken off the map
runs that much below the skeleton. Read lengths from `et$length`, which
accounts for the contracted pixels; the map is for masking and zonal
statistics, not for measuring.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- branch_order_map(skel, mask, order = "branch_order", unit = "px")
cmap <- order_classification_map(res$edges, template = skel,
                                 value = "branch_order")
} # }
```

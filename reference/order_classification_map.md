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

## Examples

``` r
if (FALSE) { # \dontrun{
res <- branch_order_map(skel, mask, order = "branch_order", unit = "px")
cmap <- order_classification_map(res$edges, template = skel,
                                 value = "branch_order")
} # }
```

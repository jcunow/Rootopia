# Plot the output of classify_soil_rgb

Produces a three-panel figure: the classified raster map, a legend
showing class names with area fractions and mean distances to centroids,
and a dendrogram of inter-class LAB distances.

## Usage

``` r
plot_soil_classification(
  result,
  color_mode = "contrast",
  class_colors = NULL,
  vibrant_colors = NULL,
  save_png = NULL,
  width = 1800,
  height = 900
)
```

## Arguments

- result:

  A list returned by
  [`classify_soil_rgb`](https://jcunow.github.io/RootScanR/reference/classify_soil_rgb.md).

- color_mode:

  Character. One of `"contrast"` (default), `"vibrant"`, or
  `"centroid"`.

  `"contrast"`

  :   Muted, distinguishable colours defined in `class_colors`.

  `"vibrant"`

  :   Saturated, high-contrast colours defined in `vibrant_colors`.

  `"centroid"`

  :   Each class rendered as its actual LAB centroid colour – useful for
      sanity-checking centroid values.

- class_colors:

  Named character vector of hex colours for `"contrast"` mode. Names
  must match class names in `result$centroids`. Defaults to a built-in
  palette.

- vibrant_colors:

  Named character vector of hex colours for `"vibrant"` mode. Same
  naming requirement. Defaults to a built-in palette.

- save_png:

  File path for PNG output, or `NULL` (default) to plot to the active
  graphics device.

- width, height:

  PNG dimensions in pixels. Only used if `save_png` is not `NULL`.

## Value

Invisibly `NULL`. Called for its side effect (plot).

## See also

[`classify_soil_rgb`](https://jcunow.github.io/RootScanR/reference/classify_soil_rgb.md)

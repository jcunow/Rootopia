# Paint a per-segment value onto the root system (PNG)

Colors each segment by any value you supply, over the mask. Segments
with `NA` are drawn in red so failures stay visible.

## Usage

``` r
render_order_overlay(
  segs,
  values,
  dims,
  file,
  mask = NULL,
  max_side = 2000,
  continuous = NA,
  value_label = NULL,
  n_bins = 64L,
  tip_order = NULL
)
```

## Arguments

- segs:

  Segment list (kept via `keep_segments = TRUE`).

- values:

  Numeric vector, one value per segment: an order column, or any
  per-segment quantity such as those added by
  [`root_transport`](https://jcunow.github.io/Rootopia/reference/root_transport.md).

- dims:

  Integer length-2 vector, the (cropped) raster dimensions.

- file:

  Output PNG path.

- mask:

  Optional 0/1 mask drawn as a gray background.

- max_side:

  Long-side cap (px) for the output image.

- continuous:

  `NA` (default) picks the legend from the data: continuous when any
  value has a fractional part or there are more than 24 distinct values,
  discrete otherwise. `TRUE` or `FALSE` forces it.

- value_label:

  Optional title drawn above the colour bar, e.g. `"path resistance"`.
  Continuous legend only.

- n_bins:

  Number of colour steps in the continuous ramp.

- tip_order:

  Deprecated name for `values`.

## Value

Invisibly, `file`.

## Details

Two legends are available and the right one is chosen from the data. A
branch order is a handful of small integers and reads as named classes,
so it gets one labelled swatch per order. A transport quantity – path
resistance, path length from the base, the root length a segment serves
– is continuous and gets a colour ramp with a scale bar instead. Pass
`continuous` to override the choice.

## See also

[`root_transport`](https://jcunow.github.io/Rootopia/reference/root_transport.md)
for quantities worth painting.

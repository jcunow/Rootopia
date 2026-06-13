# Native-resolution validation of a sub-window

Renders a small window at `scale`x magnification (no downsampling) so
the order-coloured graph can be checked against the skeleton. Skeleton
pixels are grey; traced/ordered pixels are coloured on top, so bare grey
reveals skeleton the graph missed. Coordinates are in the original image
frame.

## Usage

``` r
plot_order_window(
  et,
  skel,
  r_range,
  c_range,
  scale = 3,
  file = "window.png",
  order_col = "root_order"
)
```

## Arguments

- et:

  An `edges` table carrying `attr(., "segments")` and
  `attr(., "crop_offset")` (run with `keep_segments = TRUE`).

- skel:

  The original skeleton (`SpatRaster` or matrix) for the grey
  background.

- r_range, c_range:

  Integer length-2 row/column ranges (original coordinates).

- scale:

  Magnification factor (device px per image px).

- file:

  Output PNG path.

- order_col:

  Which order column to colour by (default `"root_order"`).

## Value

Invisibly, `file`.
